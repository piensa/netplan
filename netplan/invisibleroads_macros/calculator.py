from __future__ import absolute_import
from math import ceil, floor


def define_normalize(xs, ys):
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)

    x_width = (x_max - x_min) or 1
    y_width = (y_max - y_min) or 1

    def normalize(x):
        unit_x = (x - x_min) / float(x_width)
        return (unit_x * y_width) + y_min

    return normalize


def divide_safely(numerator, denominator, default):
    if not denominator:
        if isinstance(default, Exception):
            raise default
        return default
    return numerator / float(denominator)


def get_int(x):
    return int(ceil(x)) if x > 0 else int(floor(x))


def get_percent_change(new_value, old_value):
    if not old_value and not new_value:
        return 0
    if not old_value:
        sign = 1 if new_value > 0 else -1
        return sign * float('inf')
    change_difference = new_value - old_value
    change_ratio = change_difference / float(old_value)
    return 100 * change_ratio


def round_number(x):
    return int(round(abs(x)))
