# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json


from qiime2.core.type.template import TypeTemplate, instantiate

def is_collection_type(x):
    raise TypeError

class _CollectionBase(TypeTemplate):
    def __init__(self, fields=()):
        self.fields = fields

    def __eq__(self, other):
        return type(self) is type(other) and self.fields == other.fields

    def get_name(self):
        return self.__class__.__name__

    def get_kind(self):
        if self.fields:
            return self.fields[0].template.get_kind()
        return "collection"

    def specialize(self, fields):
        return self.__class__(fields).template

    def is_variant(self):
        return False

    def validate_predicate(self, predicate, expr):
        raise TypeError("Predicates cannot be applied to %r" % expr)

    def validate_fields(self, fields, expr):
        try:
            field, = fields
        except ValueError:
            raise TypeError

        if isinstance(field, self.__class__):
            raise TypeError

    def validate_union(self, other):
        if type(other) is not type(self):
            raise TypeError

        if other.kind != self.kind:
            raise TypeError

    def is_element(self, value, expr):
        if isinstance(value, self.view) and len(value) > 0:
            return all(v in expr.fields[0] for v in value)
        return False

    def validate_intersection(self, other):
        pass


class _1DCollectionBase(_CollectionBase):
    def get_field_names(self):
        return ['type']



@instantiate
class Set(_1DCollectionBase):
    view = set

@instantiate
class List(_1DCollectionBase):
    view = list

@instantiate
class Tuple(_CollectionBase):
    view = tuple

    def get_field_names(self):
        return ['*types']

    def validate_fields(self, fields, expr):
        # Tuples may contain anything, and as many fields as desired
        pass


