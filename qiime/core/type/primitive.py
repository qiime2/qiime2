# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from qiime.core.type import BaseType


class Map(BaseType, fields=('Key', 'Value'), variant_of=BaseType.Primitive):
    class Key:
        pass

    class Value:
        pass


class List(BaseType, fields='Contents', variant_of=BaseType.Primitive):
    class Contents:
        pass


class Bag(BaseType, fields='Contents', variant_of=BaseType.Primitive):
    class Contents:
        pass


class Set(BaseType, fields='Contents', variant_of=BaseType.Primitive):
    class Contents:
        pass


_variants = (Map.Key, Map.Value, List.Contents, Bag.Contents, Set.Contents,
             BaseType.Primitive)


class Str(BaseType, variant_of=_variants):

    def from_string(self, string):
        # TODO: is this cast necessary?
        return str(string)

    def to_string(self, data):
        # TODO: is this cast necessary?
        return str(data)


class Int(BaseType, variant_of=_variants):

    def from_string(self, string):
        return int(string)

    def to_string(self, data):
        return str(data)


class Float(BaseType, variant_of=_variants):

    def from_string(self, string):
        return float(string)

    def to_string(self, data):
        return str(data)


class Column(BaseType, variant_of=BaseType.Primitive):
    pass
