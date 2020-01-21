# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import types
import collections
import itertools

from qiime2.core.type.grammar import IncompleteExp, UnionExp, IntersectionExp
from qiime2.core.type.template import TypeTemplate, PredicateTemplate
from qiime2.core.type.util import is_semantic_type, is_qiime_type

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
    'tuple', 'row', 'record',
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
        The name of the semantic type: this should match the variable to which
        the semantic type is assigned.
    field_names : str, iterable of str, optional
        Name(s) of the fields where member types can be placed. This makes
        the type a composite type, meaning that fields must be provided to
        produce realized semantic types. These names will define ad-hoc
        variant types accessible as `name`.field[`field_names` member].
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
    field_names = _munge_field_names(field_names)
    field_members = _munge_field_members(field_names, field_members)

    return SemanticTemplate(name, field_names, field_members, variant_of)


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
    if field_names is None:
        return ()

    if type(field_names) is str:
        return (field_names,)

    field_names = tuple(field_names)
    for field_name in field_names:
        if type(field_name) is not str:
            raise ValueError("Field name %r from %r is not a string."
                             % (field_name, field_names))
    if len(set(field_names)) != len(field_names):
        raise ValueError("Duplicate field names in %r." % field_names)

    return field_names


def _munge_field_members(field_names, field_members):
    if field_names is None:
        return {}

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
        if is_qiime_type(value) and is_semantic_type(value):
            fixed[key] = (value,)
        else:
            value = tuple(value)
            for v in value:
                if not is_semantic_type(v):
                    raise ValueError("Field member: %r (of field %r) is not a"
                                     " semantic type." % (v, key))
            fixed[key] = value
    return fixed


class VariantField:
    def __init__(self, type_name, field_name, field_members):
        self.type_name = type_name
        self.field_name = field_name
        self.field_members = field_members

    def is_member(self, semantic_type):
        for field_member in self.field_members:
            if isinstance(field_member, IncompleteExp):
                # Pseudo-subtyping like Foo[X] <= Foo[Any].
                # (IncompleteExp will never have __le__ because you
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


class SemanticTemplate(TypeTemplate):
    public_proxy = 'field',

    def __init__(self, name, field_names, field_members, variant_of):
        self.name = name
        self.field_names = field_names
        self.__field = {f: VariantField(name, f, field_members[f])
                        for f in self.field_names}
        self.variant_of = variant_of

    @property
    def field(self):
        return types.MappingProxyType(self.__field)

    def __eq__(self, other):
        return (type(self) is type(other)
                and self.name == other.name
                and self.fields == other.fields
                and self.variant_of == other.variant_of)

    def __hash__(self):
        return (hash(type(self)) ^ hash(self.name)
                ^ hash(self.fields) ^ hash(self.variant_of))

    def get_kind(self):
        return 'semantic-type'

    def get_name(self):
        return self.name

    def get_field_names(self):
        return self.field_names

    def is_element_expr(self, self_expr, value):
        import qiime2.sdk
        if not isinstance(value, qiime2.sdk.Artifact):
            return False
        return value.type <= self_expr

    def is_element(self, value):
        raise NotImplementedError

    def validate_field(self, name, field):
        raise NotImplementedError

    def validate_fields_expr(self, self_expr, fields_expr):
        self.validate_field_count(len(fields_expr))
        for expr, varf in zip(fields_expr,
                              [self.field[n] for n in self.field_names]):
            if (expr.template is not None
                    and hasattr(expr.template, 'is_variant')):
                check = expr.template.is_variant
            else:
                check = self.is_variant
            if not check(expr, varf):
                raise TypeError("%r is not a variant of %r" % (expr, varf))

    @classmethod
    def is_variant(cls, expr, varf):
        if isinstance(expr, UnionExp):
            return all(cls.is_variant(e, varf) for e in expr.members)
        if isinstance(expr, IntersectionExp):
            return any(cls.is_variant(e, varf) for e in expr.members)
        return varf.is_member(expr) or varf in expr.template.variant_of

    def validate_predicate(self, predicate):
        if not isinstance(predicate, Properties):
            raise TypeError()

    def update_ast(self, ast):
        ast['builtin'] = False


class Properties(PredicateTemplate):
    def __init__(self, *include, exclude=()):
        if len(include) == 1 and isinstance(include[0],
                                            (list, tuple, set, frozenset)):
            include = tuple(include[0])

        if type(exclude) is str:
            exclude = (exclude,)

        self.include = tuple(include)
        self.exclude = tuple(exclude)
        for prop in itertools.chain(self.include, self.exclude):
            if type(prop) is not str:
                raise TypeError("%r in %r is not a string." % (prop, self))

    def __hash__(self):
        return hash(frozenset(self.include)) ^ hash(frozenset(self.exclude))

    def __eq__(self, other):
        return (type(self) is type(other) and
                set(self.include) == set(other.include) and
                set(self.exclude) == set(other.exclude))

    def __repr__(self):
        args = []
        if self.include:
            args.append(', '.join(repr(s) for s in self.include))
        if self.exclude:
            args.append("exclude=%r" % list(self.exclude))

        return "%s(%s)" % (self.__class__.__name__, ', '.join(args))

    def is_symbol_subtype(self, other):
        if type(self) is not type(other):
            return False

        return (set(other.include) <= set(self.include) and
                set(other.exclude) <= set(self.exclude))

    def is_symbol_supertype(self, other):
        if type(self) is not type(other):
            return False

        return (set(other.include) >= set(self.include) and
                set(other.exclude) >= set(self.exclude))

    def collapse_intersection(self, other):
        if type(self) is not type(other):
            return None

        new_include_set = set(self.include) | set(other.include)
        new_exclude_set = set(self.exclude) | set(other.exclude)

        new_include = []
        new_exclude = []
        for inc in itertools.chain(self.include, other.include):
            if inc in new_include_set:
                new_include.append(inc)
                new_include_set.remove(inc)
        for exc in itertools.chain(self.exclude, other.exclude):
            if exc in new_exclude_set:
                new_exclude.append(exc)
                new_exclude_set.remove(exc)

        return self.__class__(*new_include, exclude=new_exclude).template

    def get_kind(self):
        return 'semantic-type'

    def get_name(self):
        return self.__class__.__name__

    def is_element(self, expr):
        return True  # attached TypeExp checks this

    def get_union_membership_expr(self, self_expr):
        return 'predicate-' + self.get_name()

    def update_ast(self, ast):
        ast['include'] = list(self.include)
        ast['exclude'] = list(self.exclude)
