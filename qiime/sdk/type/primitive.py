# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from qiime.sdk.type import Type


class Map(Type, fields=('Key', 'Value'), variant_of=Type.Primitive):
    class Key:
        pass

    class Value:
        pass


class List(Type, fields='Contents', variant_of=Type.Primitive):
    class Contents:
        pass


class Bag(Type, fields='Contents', variant_of=Type.Primitive):
    class Contents:
        pass


class Set(Type, fields='Contents', variant_of=Type.Primitive):
    class Contents:
        pass


_variants = (Map.Key, Map.Value, List.Contents, Bag.Contents, Set.Contents,
             Type.Primitive)


class Str(Type, variant_of=_variants):
    pass


class Int(Type, variant_of=_variants):

    def from_string(self, string):
        return int(string)

    def to_string(self, data):
        return str(data)


class Float(Type, variant_of=_variants):
    pass

class Column(Type, variant_of=Type.Primitive):
    pass
