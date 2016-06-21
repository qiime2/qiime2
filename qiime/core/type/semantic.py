# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import types
import collections
import itertools

from . import grammar
from qiime.core.type.util import overrides


_RESERVED_NAMES = {
    # Predicates:
    'range', 'choice', 'properties', 'arguments',
    # Primitives:
    'integer', 'int', 'string', 'str', 'metadata', 'metadatacolumn', 'column',
    'float', 'double', 'number', 'set', 'list', 'bag', 'map', 'dict',
    # Type System:
    'semantictype', 'propertymap', 'propertiesmap', 'typemap', 'typevariable',
    'predicate'
    # TODO: expand this set.
}


def SemanticType(name, field_names=None, field_members=None, variant_of=None):
    """
    """
    # TODO: IMPORTANT! Refactor this monstrous atrocity of a factory.

    if type(name) is not str:
        raise TypeError("Names of semantic types must be strings, not %r."
                        % name)
    if name.lower() in _RESERVED_NAMES:
        raise ValueError("%r is a reserved name." % name)

    if variant_of is None:
        variant_of = ()
    elif isinstance(variant_of, VariantField):
        variant_of = (variant_of,)
    else:
        variant_of = tuple(variant_of)
        for variant in variant_of:
            if not isinstance(variant, VariantField):
                raise ValueError("Element %r of %r is not a variant field"
                                 " (ExampleType.field['name'])."
                                 % (variant, variant_of))

    if field_names is not None:

        if type(field_names) is str:
            field_names = (field_names,)
        else:
            field_names = tuple(field_names)
            for field_name in field_names:
                if type(field_name) is not str:
                    raise ValueError("Field name %r from %r is not a string."
                                     % (field_name, field_names))
            if len(set(field_names)) != len(field_names):
                raise ValueError("Duplicate field names in %r." % field_names)

        if field_members is None:
            field_members = {}
        if not isinstance(field_members, collections.Mapping):
            raise ValueError("")
        else:
            field_members = field_members.copy()
        for key, value in field_members.items():
            if key not in field_names:
                raise ValueError("Field members %r:%r are not a part of"
                                 " `field_names` (%r)."
                                 % (key, value, field_names))
            else:
                if isinstance(value, (_IncompleteSemanticType, _SemanticType)):
                    # TODO: this is disgusting
                    if (isinstance(value, _SemanticType) and
                            not value.is_concrete()):
                        raise ValueError("")
                    field_members[key] = (value,)
                else:
                    value = tuple(value)
                    for v in value:
                        if (not isinstance(v, _SemanticType) or
                                not v.is_concrete()):
                            raise ValueError("")
                    field_members[key] = value
        for key in field_names:
            if key not in field_members:
                field_members[key] = ()

        return _IncompleteSemanticType(name, field_names, field_members,
                                       variant_of)
    return _SemanticType(name, variant_of)


def is_semantic_type(type_):
    return isinstance(type_, (_SemanticMixin, _IncompleteSemanticType))


class VariantField:
    def __init__(self, type_name, field_name, field_members):
        self.type_name = type_name
        self.field_name = field_name
        # TODO verify that field_members does not have unions
        # or incomplete types in it
        self.field_members = field_members

    def is_member(self, semantic_type):
        sans_predicate = semantic_type._apply_predicate_(None)
        # TODO: Make this not awful.
        for field in self.field_members:
            if isinstance(field, _IncompleteSemanticType):
                try:
                    field = field[sans_predicate.fields]
                except TypeError:
                    pass
            if field == sans_predicate:
                return True
        return False

    def __repr__(self):
        return "%s.field[%r]" % (self.type_name, self.field_name)


class _IncompleteSemanticType(grammar.CompositeType):
    def __init__(self, name, field_names, field_members, variant_of):
        self.field_names = field_names
        self.field = types.MappingProxyType({
            f: VariantField(name, f, field_members[f])
            for f in self.field_names})
        self.variant_of = variant_of

        super().__init__(name, field_names)

    @overrides(grammar.CompositeType)
    def _validate_field_(self, name, value):
        super()._validate_field_(name, value)

        varfield = self.field[name]
        if not value.is_variant(varfield):
            raise TypeError("%r is not a variant of %r." % (value, varfield))

    @overrides(grammar.CompositeType)
    def _apply_fields_(self, fields):
        return _SemanticType(self.name, self.variant_of, fields=fields)

    def is_concrete(self):
        return False


class _SemanticMixin:
    @property
    def _concrete(self):
        return False

    def __le__(self, other):
        if type(self) is not type(other):
            return NotImplemented
        if self.name != other.name:
            return False
        for f1, f2 in itertools.zip_longest(self.fields, other.fields):
            if not (f1 <= f2):
                return False
        if other.predicate is not None:
            if not (self.predicate <= other.predicate):
                return False
        return True

    def __ge__(self, other):
        if type(self) is not type(other):
            return NotImplemented
        if self.name != other.name:
            return False
        for f1, f2 in itertools.zip_longest(self.fields, other.fields):
            if not (f1 >= f2):
                return False
        if other.predicate is not None:
            if not (self.predicate >= other.predicate):
                return False
        return True

    def is_variant(self, varfield):
        return varfield in self.variant_of or varfield.is_member(self)

    def is_concrete(self):
        if not self._concrete:
            return False
        if self.predicate is not None and not self.predicate.is_concrete():
            return False
        # all is vacuously true when self.fields is empty
        if not all(f.is_concrete() for f in self.fields):
            return False

        return True


class _SemanticType(grammar.TypeExpression, _SemanticMixin):
    @property
    @overrides(_SemanticMixin)
    def _concrete(self):
        return True

    def __init__(self, name, variant_of, **kwargs):
        self.variant_of = variant_of

        super().__init__(name, **kwargs)

    def _validate_intersection_(self, other, handshake=False):
        raise TypeError("Cannot intersect %r and %r. (Only semantic type"
                        " variables can be intersected.)" % (self, other))

    def _build_union_(self, members):
        return _SemanticUnionType(members)

    def _validate_predicate_(self, predicate):
        raise NotImplementedError("TODO")

    def _apply_predicate_(self, predicate):
        return self.__class__(self.name, self.variant_of,
                              fields=self.fields, predicate=predicate)


class _SemanticUnionType(grammar.UnionTypeExpression, _SemanticMixin):
    @overrides(_SemanticMixin)
    def is_variant(self, varfield):
        return all(m.is_variant(varfield) for m in self.members)

    def __le__(self, other):
        if type(self) is type(other):
            return self.members <= other.members
        return all(m <= other for m in self.members)

    def __ge__(self, other):
        if type(self) is type(other):
            return self.members >= other.members
        return any(m >= other for m in self.members)


class _SemanticIntersectionType(grammar.IntersectionTypeExpression,
                                _SemanticMixin):
    def __init__(self):
        # IMPORTANT! Remove this __init__ method.
        raise NotImplementedError("TODO")

    @overrides(_SemanticMixin)
    def is_variant(self, varfield):
        return all(m.is_variant(varfield) for m in self.members)


class SemanticMap(grammar.MappingTypeExpression, _SemanticMixin):
    def __init__(self):
        # IMPORTANT! Remove this __init__ method.
        raise NotImplementedError("TODO")

    @overrides(_SemanticMixin)
    def is_variant(self, varfield):
        # A TypeMap may be on either side of the signature, we only know
        # which side once it is placed in a signature. Otherwise, either
        # construction is completely valid as long as all members agree.
        return (all(m.is_variant(varfield) for m in self.mapping) or
                all(m.is_variant(varfield) for m in self.mapping.values()))

    def _validate_member_(self, member):
        if not member.is_concrete():
            raise ValueError("")

    def _build_intersection_(self, members):
        return _SemanticIntersectionType(members)


class Properties(grammar.Predicate):
    def __init__(self, include=(), exclude=()):
        raise NotImplementedError("TODO")
        self.include = tuple(sorted(include))
        self.exclude = tuple(sorted(exclude))

        super().__init__(include, exclude)

    def __hash__(self):
        pass

    def __eq__(self):
        pass

    def __repr__(self):
        args = []
        if self.include:
            args.append("include=%r" % self.include)
        if self.exclude:
            args.append("exclude=%r" % self.exclude)

        return "%s(%s)" % (self.__class.__name__, ', '.join(args))

    def is_concrete(self):
        return True


class PropertyMap(Properties):
    def __init__(self, name, mapping):
        raise NotImplementedError("TODO")
        self.name = name
        self.mapping = mapping

    def __repr__(self):
        return self.name

    @overrides(Properties)
    def is_concrete(self):
        pass
