# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
from qiime.sdk.type import Type

class PredicateMeta(type):
    def __getitem__(cls, restrictions):
        return cls(restrictions)

class Predicate(metaclass=PredicateMeta):
    def __init__(self, restrictions):
        pass

    def __invert__(self):
        return self

class Domain(Predicate):
    pass

class Property(Predicate):
    pass
