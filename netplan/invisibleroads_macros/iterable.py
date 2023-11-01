from collections import OrderedDict, defaultdict
from collections.abc import Callable, MutableSet
from copy import deepcopy

from .security import RANDOM


class OrderedDefaultDict(OrderedDict):
    # http://stackoverflow.com/a/6190500/192092

    def __init__(self, default_factory=None, *a, **kw):
        if (default_factory is not None and
                not isinstance(default_factory, Callable)):
            raise TypeError('first argument must be callable')
        OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
        return type(self), args, None, None, self.items()

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return type(self)(self.default_factory, self)

    def __deepcopy__(self, memo):
        return type(self)(self.default_factory, deepcopy(self.items()))

    def __repr__(self):
        return 'OrderedDefaultDict(%s, %s)' % (
            self.default_factory, OrderedDict.__repr__(self))


class OrderedSet(MutableSet):
    'Written by Raymond Hettinger'
    # http://code.activestate.com/recipes/576694/

    def __init__(self, iterable=None):
        self.end = end = []
        end += [None, end, end]         # sentinel node for doubly linked list
        self.map = {}                   # key --> [key, prev, next]
        if iterable is not None:
            self |= iterable

    def __len__(self):
        return len(self.map)

    def __contains__(self, key):
        return key in self.map

    def add(self, key):
        if key not in self.map:
            end = self.end
            curr = end[1]
            curr[2] = end[1] = self.map[key] = [key, curr, end]

    def discard(self, key):
        if key in self.map:
            key, prev, next = self.map.pop(key)
            prev[2] = next
            next[1] = prev

    def __iter__(self):
        end = self.end
        curr = end[2]
        while curr is not end:
            yield curr[0]
            curr = curr[2]

    def __reversed__(self):
        end = self.end
        curr = end[1]
        while curr is not end:
            yield curr[0]
            curr = curr[1]

    def pop(self, last=True):
        if not self:
            raise KeyError('set is empty')
        key = self.end[1][0] if last else self.end[2][0]
        self.discard(key)
        return key

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self))

    def __eq__(self, other):
        if isinstance(other, OrderedSet):
            return len(self) == len(other) and list(self) == list(other)
        return set(self) == set(other)


def drop_duplicates(xs):
    # https://stackoverflow.com/a/39835527/192092
    return list(OrderedDict.fromkeys(xs))


def flatten_dictionaries(dictionary_of_dictionaries):
    'Combined nested dictionary into a simple dictionary'
    d = OrderedDict()
    for outer_k, outer_v in dictionary_of_dictionaries.items():
        if isinstance(outer_v, dict):
            d[outer_k] = outer_v
            continue
        for inner_k, inner_v in flatten_dictionaries(outer_v).items():
            d['%s.%s' % (outer_k, inner_k)] = inner_v
    return d


def flatten_lists(list_of_lists):
    # http://stackoverflow.com/a/952952/192092
    return [item for sublist in list_of_lists for item in sublist]


def get_lists_from_tuples(xs):
    'Convert tuples to lists'
    # http://stackoverflow.com/a/1014669
    if isinstance(xs, (list, tuple)):
        return list(map(get_lists_from_tuples, xs))
    return xs


def get_tuples_from_lists(xs):
    'Convert tuples to lists'
    # http://stackoverflow.com/a/1014669
    if isinstance(xs, (list, tuple)):
        return tuple(map(get_tuples_from_lists, xs))
    return xs


def make_tree():
    return defaultdict(make_tree)


def merge_dictionaries(*dictionaries):
    'Combine multiple nested dictionaries and overwrite duplicate keys'
    d = OrderedDict()
    for outer_d in dictionaries:
        for outer_k, outer_v in outer_d.items():
            d1 = d.get(outer_k, {})
            d2 = outer_v
            if isinstance(d1, dict) and isinstance(d2, dict):
                d[outer_k] = merge_dictionaries(d1, d2)
            else:
                d[outer_k] = d2
    return d


def sort_dictionary(value_by_key, sorted_keys):
    'Sort dictionary by keys'
    d = OrderedDict()
    for key in sorted_keys:
        try:
            d[key] = value_by_key[key]
        except KeyError:
            pass
    return d


def shuffled(items):
    items = list(items)
    RANDOM.shuffle(items)
    return items
