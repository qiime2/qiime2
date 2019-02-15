# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import types
import collections
import itertools

from . import grammar
from qiime2.core.util import overrides

_RESERVED_NAMES = {
    # Predicates:
    'range', 'choice', 'properties', 'arguments',
    # Primitives:
    'integer', 'int', 'string', 'str', 'metadata', 'metadatacolumn',
    'categoricalmetadatacolumn', 'numericmetadatacolumn', 'column',
    'categoricalcolumn', 'numericcolumn', 'metacol', 'categoricalmetacol',
    'numericmetacol', 'metadatacategory', 'float', 'double', 'number', 'set',
    'list', 'bag', 'multiset', 'map', 'dict', 'nominal', 'ordinal',
    'categorical', 'numeric', 'interval', 'ratio', 'continuous', 'discrete',
    # Type System:
    'semantictype', 'propertymap', 'propertiesmap', 'typemap', 'typevariable',
    'predicate'
}


def _validate_name(name):
    if type(name) is not str:
        raise TypeError("Names of semantic types must be strings, not %r."
                        % name)
    if name.lower() in _RESERVED_NAMES:
        raise ValueError("%r is a reserved name." % name)


def SemanticType(name, field_names=None, field_members=None, variant_of=None):
    """Create a new semantic type.

    Parameters
    ----------
    name : str
        The name of the semantic type, this should match the variable which it
        is assigned to.
    field_names : str, iterable of str, optional
        Name(s) of the fields where member types and can be placed. This makes
        the type a composite type, meaning that fields must be provided to be
        a realized semantic type. These names will define an ad-hoc variant
        types accessible as `name`.field[`field_names` member].
    field_members : mapping, optional
        A mapping of strings in `field_names` to one or more semantic types
        which are known to be members of the field (the variant type).
    variant_of : VariantField, iterable of VariantField, optional
        Define the semantic type to be a member of one or more variant types
        allowing it to be placed in the respective fields defined by those
        variant types.

    Returns
    -------
    A Semantic Type
        There are several (private) types which may be returned, but anything
        returned by this factory will cause `is_semantic_type` to return True.

    """
    _validate_name(name)
    variant_of = _munge_variant_of(variant_of)

    if field_names is not None or field_members is not None:
        field_names = _munge_field_names(field_names)
        field_members = _munge_field_members(field_names, field_members)

        return _IncompleteSemanticType(name, field_names, field_members,
                                       variant_of)
    return _SemanticType(name, variant_of)


def _munge_variant_of(variant_of):
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
    return variant_of


def _munge_field_names(field_names):
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

    return field_names


def _munge_field_members(field_names, field_members):
    fixed = {k: () for k in field_names}

    if field_members is None:
        return fixed
    if not isinstance(field_members, collections.Mapping):
        raise ValueError("")

    fixed.update(field_members)

    for key, value in field_members.items():
        if key not in field_names:
            raise ValueError("Field member key: %r is not in `field_names`"
                             " (%r)." % (key, field_names))
        if is_semantic_type(value):
            fixed[key] = (value,)
        else:
            value = tuple(value)
            for v in value:
                if not is_semantic_type(v):
                    raise ValueError("Field member: %r (of field %r) is not a"
                                     " semantic type." % (v, key))
            fixed[key] = value
    return fixed


def is_semantic_type(type_):
    return isinstance(type_, (_SemanticMixin, _IncompleteSemanticType))


class VariantField:
    def __init__(self, type_name, field_name, field_members):
        self.type_name = type_name
        self.field_name = field_name
        self.field_members = field_members

    def is_member(self, semantic_type):
        for field_member in self.field_members:
            if isinstance(field_member, _IncompleteSemanticType):
                # Pseudo-subtyping like Foo[X] <= Foo[Any].
                # (_IncompleteSemanticType will never have __le__ because you
                #  are probably doing something wrong with it (this totally
                #                                              doesn't count!))
                if semantic_type.name == field_member.name:
                    return True
                # ... it doesn't count because this is a way of restricting our
                # ontology and isn't really crucial. Where it matters would be
                # in function application where the semantics must be defined
                # precisely and Foo[Any] is anything but precise.
            else:
                if semantic_type <= field_member:
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


class _SemanticMixin:
    def is_variant(self, varfield):
        return varfield in self.variant_of or varfield.is_member(self)


class _SemanticType(grammar.TypeExpression, _SemanticMixin):
    def __init__(self, name, variant_of, **kwargs):
        self.variant_of = variant_of

        super().__init__(name, **kwargs)

    def _apply_fields_(self, fields):
        return self.__class__(self.name, self.variant_of, fields=fields,
                              predicate=self.predicate)

    def _validate_intersection_(self, other, handshake=False):
        pass

    def _build_union_(self, *members):
        return _SemanticUnionType(*members)

    def _validate_predicate_(self, predicate):
        if not isinstance(predicate, Properties):
            raise TypeError()

    def _apply_predicate_(self, predicate):
        return self.__class__(self.name, self.variant_of,
                              fields=self.fields, predicate=predicate)

    def _is_element_(self, value):
        import qiime2.sdk
        if not isinstance(value, qiime2.sdk.Artifact):
            return False
        return value.type <= self

    def to_ast(self):
        ast = super().to_ast()
        ast['type'] = 'semantic-type'
        ast['concrete'] = self.is_concrete()
        return ast


class _SemanticUnionType(grammar.UnionTypeExpression, _SemanticMixin):
    @overrides(_SemanticMixin)
    def is_variant(self, varfield):
        return all(m.is_variant(varfield) for m in self)


class Properties(grammar.Predicate):
    def __init__(self, include=(), exclude=()):
        if type(include) is str:
            include = (include,)
        if type(exclude) is str:
            exclude = (exclude,)

        self.include = list(include)
        self.exclude = list(exclude)
        for prop in itertools.chain(self.include, self.exclude):
            if type(prop) is not str:
                raise TypeError("%r in %r is not a string." % (prop, self))

        super().__init__(include, exclude)

    def __hash__(self):
        return hash(tuple(self.include)) ^ hash(tuple(self.exclude))

    def __eq__(self, other):
        return (type(self) is type(other) and
                self.include == other.include and
                self.exclude == other.exclude)

    def __repr__(self):
        args = []
        if self.include:
            args.append("%r" % self.include)
        if self.exclude:
            args.append("exclude=%r" % self.exclude)

        return "%s(%s)" % (self.__class__.__name__, ', '.join(args))

    def _is_subtype_(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return (set(other.include) <= set(self.include) and
                set(other.exclude) <= set(self.exclude))

    def to_ast(self):
        ast = super().to_ast()
        ast['include'] = self.include
        ast['exclude'] = self.exclude
        return ast
