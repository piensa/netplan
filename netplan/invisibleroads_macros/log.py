from __future__ import print_function

try:
    from logging import NullHandler
except ImportError:  # Python 2.6
    from logging import Handler

    class NullHandler(Handler):

        def emit(self, record):
            pass
finally:
    from logging import getLogger

import re
import traceback
from collections import OrderedDict
from os.path import expanduser
from six import string_types
from sys import stderr

from .disk import COMMAND_LINE_HOME


INDENT = ' ' * 2


class LogDictionary(dict):

    def __setitem__(self, k, v):
        super(LogDictionary, self).__setitem__(k, v)
        print('%s = %s' % (k, v))


def get_log(name):
    log = getLogger(name)
    log.addHandler(NullHandler())
    return log


def filter_nested_dictionary(value_by_key, f):
    d = OrderedDict()
    for k, v in value_by_key.items():
        if f(k):
            continue
        if isinstance(v, dict):
            v = filter_nested_dictionary(v, f)
        d[k] = v
    return d


def print_error(x, *args):
    'Print to standard error stream'
    if not isinstance(x, string_types):
        x = str(x)
    print(x % args, file=stderr)


def stylize_dictionary(value_by_key, suffix_format_packs):
    d = {}
    for key, value in value_by_key.items():
        for suffix, format_value in suffix_format_packs:
            if key.endswith(suffix):
                value = format_value(value)
                break
        d[key] = value
    return d


def format_summary(value_by_key, suffix_format_packs=None):
    suffix_format_packs = list(suffix_format_packs or [])
    suffix_format_packs.extend([
        ('_folder', format_path),
        ('_path', format_path),
    ])
    return format_nested_dictionary(OrderedDict(
        value_by_key), suffix_format_packs)


def format_nested_dictionary(
        value_by_key, suffix_format_packs=None, prefix=''):
    parts = []
    for key, value in value_by_key.items():
        left_hand_side = prefix + str(key)
        if isinstance(value, dict):
            parts.append(format_nested_dictionary(
                value, suffix_format_packs, left_hand_side + '.'))
            continue
        for suffix, format_value in suffix_format_packs or []:
            if key.endswith(suffix):
                parts.append(format_assignment(
                    left_hand_side, format_value(value)))
                break
        else:
            if not isinstance(value, string_types):
                value = str(value)
            if '\n' in value:
                value = format_indented_block(value)
            parts.append(format_assignment(
                left_hand_side, value))
    return '\n'.join(parts)


def format_assignment(left_hand_side, right_hand_side):
    left_hand_side = left_hand_side.strip()
    if right_hand_side.startswith('\n'):
        operator = ' ='
    else:
        operator = ' = '
        right_hand_side = right_hand_side.strip()
    return left_hand_side + operator + right_hand_side


def format_path(x):
    return re.sub(r'^' + expanduser('~'), COMMAND_LINE_HOME, x)


def format_indented_block(x, indent=INDENT):
    return '\n' + '\n'.join(indent + line for line in x.splitlines())


def format_decimal(x, fractional_digit_count=2):
    template = '{:,.%sf}' % fractional_digit_count
    return template.format(x)


def format_delta(x):
    template = '{:+,}'
    return template.format(x)


def format_number(x):
    return format_decimal(x, fractional_digit_count=0)


def get_nested_dictionary(nested_lists):
    d = OrderedDict()
    for k, v in nested_lists:
        if k.endswith('_'):
            try:
                v = get_nested_dictionary(v)
            except TypeError:
                pass
            else:
                k = k[:-1]
        d[k] = v
    return d


def get_nested_lists(nested_dictionary):
    xs = []
    for k, v in nested_dictionary.items():
        if isinstance(v, dict):
            v = get_nested_lists(v)
            k = k + '_'
        xs.append((k, v))
    return xs


def parse_nested_dictionary_from(raw_dictionary, max_depth=float('inf')):
    value_by_key = OrderedDict()
    for key, value in OrderedDict(raw_dictionary).items():
        key_parts = key.split('.')
        d = value_by_key
        depth = 0
        while key_parts:
            key_part = key_parts.pop(0)
            if len(key_parts) and depth < max_depth:
                if key_part not in d:
                    d[key_part] = OrderedDict()
                d = d[key_part]
                depth += 1
            else:
                d['.'.join([key_part] + key_parts)] = value
                break
    return value_by_key


def parse_nested_dictionary(
        text, is_key=lambda x: True, indent=INDENT, max_depth=float('inf')):
    raw_dictionary = parse_raw_dictionary(text, is_key, indent)
    return parse_nested_dictionary_from(raw_dictionary, max_depth)


def parse_raw_dictionary(text, is_key=lambda x: True, indent=INDENT):
    raw_dictionary, key = OrderedDict(), None
    for line in text.splitlines():
        if line.startswith(indent):
            if key is not None:
                value = line[2:].rstrip()
                raw_dictionary[key].append(value)
            continue
        try:
            key, value = line.split(' = ', 1)
            if not is_key(key):
                key = None
        except ValueError:
            key = None
        if not key:
            continue
        key = key.strip()
        value = value.strip()
        raw_dictionary[key] = [value]
    for k, v in raw_dictionary.items():
        raw_dictionary[k] = '\n'.join(v).strip()
    return raw_dictionary


def log_traceback(log, d=None):
    log.error(traceback.format_exc() + format_summary(d or {}))
