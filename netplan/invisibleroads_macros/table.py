import re
from functools import partial

from .text import compact_whitespace


UPPER_LOWER_PATTERN = re.compile(r'(.)([A-Z][a-z]+)')
LOWER_UPPER_PATTERN = re.compile(r'([a-z0-9])([A-Z])')
LETTER_DIGIT_PATTERN = re.compile(r'([a-z])([0-9])')
DIGIT_LETTER_PATTERN = re.compile(r'([0-9])([a-z])')


def duplicate_selected_column_names(selected_column_names, column_names):
    suffix = '*'
    while set(column_names).intersection(
            x + suffix for x in selected_column_names):
        suffix += '*'
    return list(column_names) + [x + suffix for x in selected_column_names]


def load_csv_safely(path, **kw):
    try:
        from pandas import read_csv
    except ImportError:
        import pip
        pip.main('install pandas'.split())
        from pandas import read_csv
    if 'skipinitialspace' not in kw:
        kw['skipinitialspace'] = True
    f = partial(read_csv, **kw)
    try:
        return f(path, encoding='utf-8')
    except UnicodeDecodeError:
        pass
    try:
        return f(path, encoding='latin-1')
    except UnicodeDecodeError:
        pass
    return f(open(path, errors='replace'))


def normalize_key(
        x, word_separator=' ', separate_camel_case=False,
        separate_letter_digit=False):
    """
    Normalize key using a variation of the method described in
    http://stackoverflow.com/a/1176023/192092

    ONETwo   one two
    OneTwo   one two
    one-two  one two
    one_two  one two
    one2     one 2
    1two     1 two
    """
    if separate_camel_case:
        x = UPPER_LOWER_PATTERN.sub(r'\1 \2', x)
        x = LOWER_UPPER_PATTERN.sub(r'\1 \2', x)
    x = x.lower()
    if separate_letter_digit:
        x = LETTER_DIGIT_PATTERN.sub(r'\1 \2', x)
        x = DIGIT_LETTER_PATTERN.sub(r'\1 \2', x)
    word_separators = ['-', '_', ' ']
    if word_separator not in word_separators:
        word_separators.append(word_separator)
    word_separator_expression = '[' + ''.join(word_separators) + ']'
    word_separator_pattern = re.compile(word_separator_expression)
    x = word_separator_pattern.sub(' ', x)
    x = compact_whitespace(x)
    return x.replace(' ', word_separator)
