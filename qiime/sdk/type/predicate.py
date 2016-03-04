# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
from qiime.sdk.type import Type


class PredicateMeta(type):
    def __getitem__(cls, *restrictions):
        return cls(restrictions)


class PredicateConstructor(metaclass=PredicateMeta):
    def __init__(self, predicate):
        self._predicate = predicate

    def __invert__(self):
        return self(lambda x: return not self._predicate(x))



class Domain(Predicate):
    def __init__(self, args):
        arg_types = set(map(type, args))


        def predicate(value):
            pass

        super().__init__(predicate)

class Property(Predicate):
    pass
