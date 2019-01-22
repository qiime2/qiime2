# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json

from . import grammar
from .primitive import _PrimitiveBase, is_primitive_type, Metadata
from .semantic import _SemanticMixin, is_semantic_type


def is_collection_type(type_):
    return isinstance(type_, (_CollectionBase, _CollectionExpression))


class _CollectionBase(grammar.CompositeType):
    def __init__(self, name, view):
        # Only 1d collections are supported for now
        self._view = view
        super().__init__(name, field_names=['elements'])

    def _apply_fields_(self, fields):
        elements, = fields
        if is_collection_type(elements):
            # This must be the first branch, as is_primitive/semantic is true
            # for collections types as well.
            raise TypeError("Cannot nest collection types.")
        elif is_semantic_type(elements):
            return _CollectionSemantic(self.name, self._view, fields=fields)
        elif is_primitive_type(elements):
            # TODO consider making an `is_metadata_type` helper if this check
            # is repeated in other parts of the codebase
            if elements is Metadata or elements.name == 'MetadataColumn':
                raise TypeError("Cannot use collections on metadata.")
            return _CollectionPrimitive(self.name, self._view, fields=fields)
        else:
            raise NotImplementedError


class _CollectionExpression(grammar.TypeExpression):
    def __init__(self, name, view, fields=(), predicate=None):
        self._view = view
        super().__init__(name, fields, predicate)

    def _is_element_(self, value):
        contained_type, = self.fields
        if not isinstance(value, self._view):
            return False
        if not value:
            # Collections of size 1 are not allowed as it overlaps with
            # optional values in an awkward way that is difficult to explain
            # to end-users and make interfaces more difficult.
            return False
        for el in value:
            if el not in contained_type:
                return False

        return True

    def to_ast(self):
        ast = super().to_ast()
        ast['type'] = 'collection'
        return ast

    def _validate_union_(self, other, handshake=False):
        # This is actually handled really well by the grammar, but interfaces
        # would need to observe unions which may have multiple collections
        # or mixes of collections and singletons, which is more complicated
        # and we don't need it.
        # If this is supported in the future, we should check that the type
        # of self and other are the same so that we don't mix primitive and
        # semantic types.
        raise TypeError("Unions of collection types are not supported.")

    def _validate_predicate_(self, predicate):
        # Something like `MinLen(...)` might be handy...
        raise TypeError("There are no predicates for collection types.")

    def _apply_fields_(self, fields):
        # Just pass the `view` along.
        return self.__class__(self.name, self._view, fields=fields,
                              predicate=self.predicate)

    def is_concrete(self):
        # Prevent use as an output or artifact type
        return False


class _CollectionPrimitive(_CollectionExpression, _PrimitiveBase):
    def encode(self, value):
        return json.dumps(list(value))

    def decode(self, string):
        return self._view(json.loads(string))


class _CollectionSemantic(_CollectionExpression, _SemanticMixin):
    def is_variant(self, varfield):
        # Collections shouldn't be part of a composite type
        return False


List = _CollectionBase('List', list)
Set = _CollectionBase('Set', set)
