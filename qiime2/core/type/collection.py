# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json

from qiime2.core.type.template import TypeTemplate


class _CollectionBase(TypeTemplate):
    public_proxy = 'encode', 'decode'

    def __init__(self):
        # For semantic types
        self.variant_of = frozenset()

    def __eq__(self, other):
        return type(self) is type(other)

    def get_name(self):
        return self.__class__.__name__[1:]  # drop `_`

    def get_kind_expr(self, self_expr):
        if self_expr.fields:
            return self_expr.fields[0].kind
        return ""

    def get_kind(self):
        raise NotImplementedError

    def is_variant(self, self_expr, varfield):
        return False

    def validate_predicate(self, predicate):
        raise TypeError("Predicates cannot be applied to %r" % self.get_name())

    def is_element_expr(self, self_expr, value):
        contained_expr = self_expr.fields[0]
        if isinstance(value, self._view) and len(value) > 0:
            return all(v in contained_expr for v in value)
        return False

    def is_element(self, value):
        raise NotImplementedError

    def get_union_membership_expr(self, self_expr):
        return self.get_name() + '-' + self.get_kind_expr(self_expr)

    # For primitive types
    def encode(self, value):
        return json.dumps(list(value))

    def decode(self, string):
        return self._view(json.loads(string))


class _1DCollectionBase(_CollectionBase):
    def validate_field(self, name, field):
        if isinstance(field, _1DCollectionBase):
            raise TypeError("Cannot nest collection types.")
        if field.get_name() in {'MetadataColumn', 'Metadata'}:
            raise TypeError("Cannot use %r with metadata." % self.get_name())

    def get_field_names(self):
        return ['type']


class _Set(_1DCollectionBase):
    _view = set


class _List(_1DCollectionBase):
    _view = list


class _Tuple(_CollectionBase):
    _view = tuple

    def get_kind_expr(self, self_expr):
        return ""

    def get_field_names(self):
        return ['*types']

    def validate_field_count(self, count):
        if not count:
            raise TypeError("Tuple type must contain at least one element.")

    def validate_field(self, name, field):
        # Tuples may contain anything, and as many fields as desired
        pass


Set = _Set()
List = _List()
Tuple = _Tuple()
