from .security import RANDOM


def make_random_color():
    return '#%06x' % RANDOM.randint(0, 0xFFFFFF)
