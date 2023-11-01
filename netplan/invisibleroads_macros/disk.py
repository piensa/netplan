import codecs
import fnmatch
import os
import re
import tarfile
import tempfile
import zipfile
from contextlib import contextmanager
from os import chdir, close, getcwd, listdir, makedirs, remove, walk
from os.path import (
    abspath, basename, dirname, exists, expanduser, isdir, islink, join,
    realpath, relpath, sep)
from shutil import copy2, copyfileobj, move, rmtree
from tempfile import _RandomNameSequence, mkdtemp, mkstemp

from .exceptions import BadArchive, BadFormat, BadPath
from .security import make_random_string, ALPHABET
from .text import unicode_safely


ARCHIVE_EXTENSIONS = '.tar.gz', '.tar.xz', '.zip'
COMMAND_LINE_HOME = '%UserProfile%' if os.name == 'nt' else '~'
HOME_FOLDER = expanduser('~')
TEMPORARY_FOLDER = expanduser('~/.tmp')
_MINIMUM_UNIQUE_LENGTH = 10


class TemporaryStorage(object):

    def __init__(self, parent_folder=None, suffix='', prefix=''):
        if parent_folder is None:
            parent_folder = make_folder(TEMPORARY_FOLDER)
        self.folder = make_unique_folder(parent_folder, suffix, prefix)

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback):
        remove_safely(self.folder)


class _CustomRandomNameSequence(_RandomNameSequence):

    def next(self):
        return self.__next__()

    def __next__(self):
        choose = self.rng.choice
        return ''.join(choose(ALPHABET) for x in range(_MINIMUM_UNIQUE_LENGTH))


def make_folder(folder):
    'Make sure a folder exists without raising an exception'
    try:
        makedirs(folder)
    except OSError:
        pass
    return folder


def make_unique_folder(
        parent_folder=None, suffix='', prefix='',
        length=_MINIMUM_UNIQUE_LENGTH):
    if parent_folder:
        make_folder(parent_folder)
    suffix = _prepare_suffix(suffix, length)
    return mkdtemp(suffix, prefix, parent_folder)


def make_unique_path(
        parent_folder=None, suffix='', prefix='',
        length=_MINIMUM_UNIQUE_LENGTH):
    if parent_folder:
        make_folder(parent_folder)
    suffix = _prepare_suffix(suffix, length)
    descriptor, path = mkstemp(suffix, prefix, parent_folder)
    close(descriptor)
    return path


def clean_folder(folder):
    'Remove folder contents but keep folder'
    for x_name in listdir(make_folder(folder)):
        x_path = join(folder, x_name)
        remove_safely(x_path)
    return folder


def copy_folder(target_folder, source_folder):
    'Copy contents without removing target_folder'
    from os import readlink, symlink
    if exists(target_folder):
        clean_folder(target_folder)
    else:
        make_folder(target_folder)
    for old_name in listdir(source_folder):
        old_path = join(source_folder, old_name)
        new_path = join(target_folder, old_name)
        if islink(old_path):
            symlink(readlink(old_path), new_path)
        elif isdir(old_path):
            copy_folder(new_path, old_path)
        else:
            copy2(old_path, target_folder)
    return target_folder


def move_folder(target_folder, source_folder):
    'Move contents without removing target_folder'
    if not exists(target_folder):
        move(source_folder, target_folder)
        return target_folder
    clean_folder(target_folder)
    for x_name in listdir(source_folder):
        x_path = join(source_folder, x_name)
        move(x_path, target_folder)
    remove_safely(source_folder)
    return target_folder


def remove_safely(folder_or_path):
    'Make sure a path or folder does not exist without raising an exception'
    try:
        rmtree(folder_or_path)
    except OSError:
        try:
            remove(folder_or_path)
        except OSError:
            pass
    return folder_or_path


def find_path(folder, file_name):
    'Locate file in folder'
    for root_folder, folder_names, file_names in walk(folder):
        if file_name in file_names:
            file_path = join(root_folder, file_name)
            break
    else:
        raise IOError('cannot find {0} in {1}'.format(file_name, folder))
    return file_path


def find_paths(folder, include_expression='*', exclude_expression=''):
    'Locate files in folder matching expression'
    return [
        unicode_safely(join(root_folder, file_name))
        for root_folder, folder_names, file_names in walk(folder)
        for file_name in fnmatch.filter(file_names, include_expression)
        if not fnmatch.fnmatch(file_name, exclude_expression)]


def get_relative_path(
        absolute_or_folder_relative_path, folder, external_folders=None,
        resolve_links=True):
    if not absolute_or_folder_relative_path:
        return absolute_or_folder_relative_path
    expanded_folder = expanduser(folder)
    absolute_folder = abspath(expanded_folder)
    absolute_path = get_absolute_path(
        absolute_or_folder_relative_path, folder, external_folders,
        resolve_links)
    return relpath(absolute_path, absolute_folder)


def get_absolute_path(
        absolute_or_folder_relative_path, folder, external_folders=None,
        resolve_links=True):
    if not absolute_or_folder_relative_path:
        return absolute_or_folder_relative_path
    expanded_path = expanduser(absolute_or_folder_relative_path)
    expanded_folder = expanduser(folder)
    absolute_path = abspath(join(expanded_folder, expanded_path))
    if external_folders == '*':
        return absolute_path
    absolute_folder = abspath(expanded_folder)
    get_path = realpath if resolve_links else lambda x: x
    real_path = get_path(absolute_path)
    real_folder = get_path(absolute_folder)
    for external_folder in external_folders or []:
        external_folder = get_path(expanduser(external_folder))
        if real_path.startswith(external_folder):
            break
    else:
        if relpath(real_path, real_folder).startswith('..'):
            raise BadPath('%s is not in %s' % (real_path, real_folder))
    return absolute_path


def has_archive_extension(path):
    for extension in ARCHIVE_EXTENSIONS:
        if path.endswith(extension):
            return True
    return False


def compress(
        source_folder, target_path=None, external_folders=None, excludes=None):
    """Compress folder.
    Specify archive extension (.tar.gz .tar.xz .zip) in target_path."""
    if not target_path:
        target_path = source_folder + '.tar.gz'
    if target_path.endswith('.tar.gz') or target_path.endswith('.tar.xz'):
        compress_tar(source_folder, target_path, external_folders, excludes)
    elif target_path.endswith('.zip'):
        compress_zip(source_folder, target_path, external_folders, excludes)
    else:
        raise BadFormat('compression format not supported (%s)' % target_path)
    return target_path


def compress_tar(
        source_folder, target_path=None, external_folders=None, excludes=None):
    'Compress folder as tar archive'
    if not target_path:
        target_path = source_folder + '.tar.gz'
    compression_format = 'xz' if target_path.endswith('.xz') else 'gz'
    folder = realpath(source_folder)
    with tarfile.open(
        target_path, 'w:' + compression_format, dereference=True,
    ) as target_file:
        _process_folder(folder, excludes, external_folders, target_file.add)
    return target_path


def compress_zip(
        source_folder, target_path=None, external_folders=None, excludes=None):
    'Compress folder as zip archive'
    if not target_path:
        target_path = source_folder + '.zip'
    source_folder = realpath(source_folder)
    with zipfile.ZipFile(
        target_path, 'w', zipfile.ZIP_DEFLATED, allowZip64=True,
    ) as target_file:
        _process_folder(
            source_folder, excludes, external_folders, target_file.write)
    return target_path


def uncompress(source_path, target_folder=None):
    if not exists(source_path):
        raise IOError('file not found (%s)' % source_path)
    if source_path.endswith('.tar.gz') or source_path.endswith('.tar.xz'):
        compression_format = 'xz' if source_path.endswith('.xz') else 'gz'
        try:
            source_file = tarfile.open(source_path, 'r:' + compression_format)
        except tarfile.ReadError:
            raise BadArchive('archive unreadable (%s)' % source_path)
        extension_expression = r'\.tar\.%s$' % compression_format
    elif source_path.endswith('.zip'):
        try:
            source_file = zipfile.ZipFile(source_path, 'r')
        except zipfile.BadZipfile:
            raise BadArchive('archive unreadable (%s)' % source_path)
        extension_expression = r'\.zip$'
    else:
        raise BadFormat('compression format not supported (%s)' % source_path)
    default_target_folder = re.sub(extension_expression, '', source_path)
    target_folder = target_folder or default_target_folder
    source_file.extractall(target_folder)
    source_file.close()
    return target_folder


def are_same_path(path1, path2):
    return realpath(expand_path(path1)) == realpath(expand_path(path2))


def is_x_parent_of_y(folder, path):
    return realpath(path).startswith(realpath(folder))


def has_name_match(path, expressions):
    name = basename(str(path))
    for expression in expressions or []:
        if fnmatch.fnmatch(name, expression):
            return True
    return False


@contextmanager
def cd(target_folder):
    source_folder = getcwd()
    try:
        chdir(target_folder)
        yield
    finally:
        chdir(source_folder)


def make_enumerated_folder_for(script_path, first_index=1):
    script_name = get_file_stem(script_path)
    if 'run' == script_name:
        script_name = get_file_stem(dirname(abspath(script_path)))
    return make_enumerated_folder(join(sep, 'tmp', script_name), first_index)


def make_enumerated_folder(base_folder, first_index=1):
    'Make a unique enumerated folder in base_folder'

    def suggest_folder(index):
        return join(base_folder, str(index))

    target_index = first_index
    target_folder = suggest_folder(target_index)
    while True:
        try:
            makedirs(target_folder)
            break
        except OSError:
            target_index += 1
            target_folder = suggest_folder(target_index)
    return target_folder


def change_owner_and_group_recursively(target_folder, target_username):
    'Change uid and gid of folder and its contents, treating links as files'
    from os import lchown     # Undefined in Windows
    from pwd import getpwnam  # Undefined in Windows
    pw_record = getpwnam(target_username)
    target_uid = pw_record.pw_uid
    target_gid = pw_record.pw_gid
    for root_folder, folders, names in walk(target_folder):
        for folder in folders:
            lchown(join(root_folder, folder), target_uid, target_gid)
        for name in names:
            lchown(join(root_folder, name), target_uid, target_gid)
    lchown(target_folder, target_uid, target_gid)


def replace_file_extension(path, extension):
    parent_folder = dirname(path)
    file_basename, file_extension = basename(path).split('.', 1)
    return join(parent_folder, file_basename + extension)


def strip_file_extension(path):
    parent_folder = dirname(path)
    file_basename = basename(path).split('.', 1)[0]
    return join(parent_folder, file_basename)


def get_file_stem(path):
    'Return file name without extension (x/y/z/file.txt.zip -> file)'
    return basename(path).split('.', 1)[0]


def get_file_extension(path, max_length=16):
    # Extract extension
    try:
        file_extension = basename(path).split('.', 1)[1]
    except IndexError:
        return ''
    # Sanitize characters
    file_extension = ''.join(x for x in file_extension if x.isalnum() or x in [
        '.', '-', '_',
    ]).rstrip()
    # Limit length
    return '.' + file_extension[-max_length:]


def load_text(source_path):
    return codecs.open(source_path, 'r', encoding='utf-8').read()


def copy_text(target_path, source_text):
    prepare_path(target_path)
    codecs.open(target_path, 'w', encoding='utf-8').write(source_text)
    return target_path


def copy_file(target_path, source_file):
    prepare_path(target_path)
    copyfileobj(source_file, open(target_path, 'wb'))
    return target_path


def copy_path(target_path, source_path):
    prepare_path(target_path)
    copy2(source_path, target_path)
    return target_path


def link_safely(target_folder_or_path, source_folder_or_path):
    try:
        return make_hard_link(target_folder_or_path, source_folder_or_path)
    except (OSError, ValueError):
        return make_soft_link(target_folder_or_path, source_folder_or_path)


def make_hard_link(target_path, source_path):
    target_path = expanduser(target_path)
    source_path = expanduser(source_path)
    if not exists(source_path):
        raise IOError('file not found (%s)' % source_path)
    if isdir(source_path):
        raise ValueError('link not valid (%s)' % source_path)
    try:
        f = os.link
    except AttributeError:
        # Copy because the function is not available in Windows
        return copy_path(target_path, source_path)
    prepare_path(target_path)
    f(source_path, target_path)
    return target_path


def make_soft_link(target_path, source_path):
    target_path = expanduser(target_path)
    source_path = expanduser(source_path)
    if not exists(source_path):
        raise IOError('file not found (%s)' % source_path)
    if realpath(target_path) == source_path:
        return target_path
    if is_x_parent_of_y(target_path, source_path):
        raise ValueError(
            'link not valid (target_path="%s", source_path="%s")' % (
                target_path, source_path))
    try:
        f = os.symlink
    except AttributeError:
        # Copy because the function is not available in Windows
        return copy_path(target_path, source_path)
    prepare_path(target_path)
    f(source_path, target_path)
    return target_path


def move_path(target_path, source_path):
    prepare_path(target_path)
    move(source_path, target_path)
    return target_path


def expand_path(path):
    return abspath(expanduser(path))


def prepare_path(path):
    make_folder(dirname(remove_safely(path)))
    return path


def _prepare_suffix(suffix, length):
    if length < _MINIMUM_UNIQUE_LENGTH:
        raise ValueError(
            'length must be greater than %s' % _MINIMUM_UNIQUE_LENGTH)
    return make_random_string(length - _MINIMUM_UNIQUE_LENGTH) + suffix


def _process_folder(source_folder, excludes, external_folders, write_path):
    for root_folder, folders, names in walk(source_folder, followlinks=True):
        for source_name in folders + names:
            if has_name_match(source_name, excludes):
                continue
            source_path = join(root_folder, source_name)
            try:
                target_path = get_relative_path(
                    source_path, source_folder, external_folders)
            except BadPath:
                continue
            write_path(realpath(source_path), target_path)


tempfile._name_sequence = _CustomRandomNameSequence()
