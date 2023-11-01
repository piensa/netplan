import attr
import codecs
import configparser
import functools
import re
import shlex
from argparse import ArgumentError, ArgumentParser
from collections import OrderedDict
from importlib import import_module
from os.path import dirname, expanduser, realpath
from six import string_types

from .calculator import get_int
from .disk import expand_path, get_absolute_path, get_relative_path
from .exceptions import BadPath
from .iterable import merge_dictionaries
from .log import format_path, format_summary, get_log


L = get_log(__name__)
SECTION_TEMPLATE = '[%s]\n%s\n'


class Settings(dict):

    def set(self, settings, prefix, key, default=None, parse=None):
        value = set_default(settings, prefix + key, default, parse)
        self[key] = value
        return value


class StoicArgumentParser(ArgumentParser):

    def add_argument(self, *args, **kw):
        try:
            return super(StoicArgumentParser, self).add_argument(*args, **kw)
        except ArgumentError:
            pass


class TerseArgumentParser(StoicArgumentParser):

    def add(self, argument_name, *args, **kw):
        d = {}
        if argument_name.endswith('_path'):
            d['metavar'] = 'PATH'
        elif '_as_percent_of_' in argument_name:
            d['metavar'] = 'FLOAT'
            d['type'] = float
        elif argument_name.endswith('_year'):
            d['metavar'] = 'YEAR'
            d['type'] = int
        elif argument_name.endswith('_in_years'):
            d['metavar'] = 'INTEGER'
            d['type'] = int
        d.update(kw)
        self.add_argument('--' + argument_name, *args, **d)


class RawCaseSensitiveConfigParser(configparser.RawConfigParser):
    optionxform = str


def set_default(settings, key, default=None, parse=None):
    value = settings.get(key, default)
    if key not in settings:
        L.warning('using default %s = %s' % (key, value))
    elif value in ('', None):
        L.warning('missing %s' % key)
    elif parse:
        value = parse(value)
    settings[key] = value
    return value


def save_settings(
        configuration_path, settings_by_section_name,
        suffix_format_packs=None):
    configuration_file = codecs.open(configuration_path, 'w', encoding='utf-8')
    configuration_file.write(format_settings(
        settings_by_section_name, suffix_format_packs) + '\n')
    return configuration_path


def load_relative_settings(
        configuration_path, section_name, external_folders=None):
    settings = OrderedDict()
    configuration_folder = dirname(expand_path(configuration_path))
    for k, v in load_settings(configuration_path, section_name).items():
        if k.endswith('_path') or k.endswith('_folder'):
            v = get_absolute_path(v, configuration_folder, external_folders)
        settings[k] = v
    return settings


def load_settings(configuration_path, section_name):
    configuration = RawCaseSensitiveConfigParser()
    configuration.read(configuration_path, 'utf-8')
    try:
        items = configuration.items(section_name)
    except configparser.NoSectionError:
        items = []
    return OrderedDict(items)


def format_settings(settings_by_section_name, suffix_format_packs=None):
    configuration_parts = []
    for section_name, settings in settings_by_section_name.items():
        configuration_parts.append(SECTION_TEMPLATE % (
            section_name, format_summary(settings, suffix_format_packs)))
    return '\n'.join(configuration_parts).strip()


def gather_settings(settings, prefix, parse_setting=None):
    d = {}
    prefix_pattern = re.compile('^' + prefix.replace('.', r'\.'))
    if not parse_setting:
        parse_setting = parse_raw_setting
    for k, v in settings.items():
        if not k.startswith(prefix):
            continue
        d = merge_dictionaries(d, parse_setting(prefix_pattern.sub('', k), v))
    return d


def parse_raw_setting(k, v):
    return {k: v}


def resolve_attribute(attribute_spec):
    # Modified from pkg_resources.EntryPoint.resolve()
    if not attribute_spec or not hasattr(attribute_spec, 'split'):
        return attribute_spec
    module_url, attributes_string = attribute_spec.split(':')
    module = import_module(module_url)
    try:
        attribute = functools.reduce(
            getattr, attributes_string.split('.'), module)
    except AttributeError:
        raise ImportError('could not resolve attribute (%s)' % module_url)
    return attribute


def split_arguments(command_string):
    try:
        xs = shlex.split(command_string)
    except UnicodeEncodeError:
        xs = shlex.split(command_string.encode('utf-8'))
    return [x.strip() for x in xs]


def make_absolute_paths(d, folder, external_folders=False):
    d = OrderedDict(d)
    for k, v in d.items():
        if isinstance(v, dict):
            v = make_absolute_paths(v, folder, external_folders)
        elif is_path_key(k) and v:
            try:
                v = get_absolute_path(v, folder, external_folders)
            except BadPath:
                _log_bad_path(v)
                v = ''
        d[k] = v
    return d


def make_relative_paths(d, folder, external_folders=None):
    d = OrderedDict(d)
    for k, v in d.items():
        if isinstance(v, dict):
            v = make_relative_paths(v, folder, external_folders)
        elif is_path_key(k) and v:
            try:
                v = get_relative_path(v, folder, external_folders)
            except BadPath:
                _log_bad_path(v)
                v = ''
        d[k] = v
    return d


def is_path_key(k):
    return k.endswith('_path') or k.endswith('_folder')


def encode_object(o):
    """Use with json.dumps to serialize classes to JSON.

    x = json.dumps(d, default=encode_object)
    """
    if hasattr(o, '__attrs_attrs__'):
        d = attr.asdict(o)
        d['__class__'] = o.__class__.__name__
        return d
    raise TypeError(repr(o) + ' is not JSON serializable')


def define_decode_object(class_by_name):
    """Use with json.loads to deserialize classes from JSON.

    decode_object = define_decode_object(globals())
    d = json.loads(x, object_hook=decode_object)
    """
    def decode_object(d):
        if '__class__' in d:
            class_name = d.pop('__class__')
            Class = class_by_name[class_name]
            return Class(**d)
        return d
    return decode_object


def define_gather_numbers(expression):

    def gather_numbers(settings):
        numbers = []
        number_pattern = re.compile(expression)
        for k, v in settings.items():
            match = number_pattern.match(k)
            if not match:
                continue
            numbers.append(int(match.group(1)))
        return sorted(numbers)

    return gather_numbers


def parse_integers(x):
    return [int(y) for y in parse_list(x)]


def parse_list(x):
    if isinstance(x, string_types):
        x = x.split()
    return x


def parse_minute_count(x):
    return get_int(parse_second_count(x) / 60.)


def parse_second_count(x):
    if isinstance(x, int):
        return x
    if not isinstance(x, string_types):
        raise ValueError
    try:
        x_count, x_unit = re.match('(\d+)([hms])', x.strip()).groups()
    except AttributeError:
        raise ValueError
    x_count = int(x_count)
    if x_unit == 'h':
        x_count *= 60 * 60
    elif x_unit == 'm':
        x_count *= 60
    return x_count


def _log_bad_path(path):
    real_path = realpath(expanduser(path))
    L.warning('bad path ignored (%s -> %s)' % (
        format_path(path),
        format_path(real_path)))
