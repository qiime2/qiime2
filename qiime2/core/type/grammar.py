# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import itertools

from qiime2.core.util import tuplize
from ..util import ImmutableBase


def common_subtypes(*types):
    smallest = []
    greatest = []
    for t in types:
        for i, p in enumerate(greatest):
            if t <= p:
                smallest[i] = t
                break
            elif t >= p:
                greatest[i] = t
                break
        else:  # evaled if the above loop did not break
            smallest.append(t)
            greatest.append(t)

    return smallest, greatest


class _TypeBase(ImmutableBase):
    """Provides reflexive methods."""

    def __ne__(self, other):
        return not self == other

    def __rmod__(self, predicate):
        raise TypeError("Predicate must be applied to the right-hand side of"
                        " a type expression.")

    def __ror__(self, other):
        return self | other  # union should be associative

    def __rand__(self, other):
        return self & other  # intersection should be associative


class CompositeType(_TypeBase):
    def __init__(self, name, field_names):
        # These classes aren't user-facing, but some light validation avoids
        # accidental issues. However, we don't want to waste a lot of time with
        # validation here, validation should happen elsewhere.
        if not len(field_names):
            raise ValueError("`field_names` cannot be an empty array-like.")

        self.name = name
        self.field_names = field_names

        self._freeze_()

    def __contains__(self, value):
        raise TypeError("Cannot check membership of %r, %r is missing"
                        " arguments for its fields." % (value, self))

    def __mod__(self, predicate):
        raise TypeError("Cannot apply predicate %r, %r is missing arguments"
                        " for its fields." % (predicate, self))

    def __or__(self, other):
        raise TypeError("Cannot union with %r, %r is missing arguments"
                        " for its fields." % (other, self))

    def __and__(self, other):
        raise TypeError("Cannot intersect with %r, %r is missing arguments"
                        " for its fields." % (other, self))

    def __getitem__(self, fields):
        fields = tuplize(fields)
        if len(fields) != len(self.field_names):
            raise TypeError("%r takes %d field(s), %d provided."
                            % (self, len(self.field_names), len(fields)))
        for args in zip(self.field_names, fields):
            self._validate_field_(*args)

        return self._apply_fields_(fields=fields)

    def __repr__(self):
        return "%s[%s]" % (self.name,
                           ', '.join('{%s}' % f for f in self.field_names))

    def _validate_field_(self, name, value):
        """Called when a field is provided to a `CompositeType`.

        This method is designed to be overridden to influence the behavior of
        the grammar. It is recommended to call super as the default
        implementation includes useful type checks and errors.

        Parameters
        ----------
        name : str
            The name of the field being set
        value : TypeExpression
            The value of the field being set

        Raises
        ------
        TypeError
            Raised when the field is rejected. By default this is when a field
            is not provided a `TypeExpression`.

        """
        if not isinstance(value, TypeExpression):
            if isinstance(value, self.__class__):
                raise TypeError("Incomplete type %r provided as a field to %r"
                                % (value, self))
            raise TypeError("%r cannot be used as a field to %r (not a type)."
                            % (value, self))

    def _apply_fields_(self, fields):
        """Called when a `CompositeType` is promoted to a `TypeExpression`.

        This method is designed to be overridden to influence the behaviour of
        the grammar. An overriding method should ensure that `self.name` is
        propogated and that the provided `fields` are passed.

        Parameters
        ----------
        fields : Tuple[TypeExpression, ...]
            The fields which should be provided to the `TypeExpression`

        Returns
        -------
        TypeExpression
            Typically this will return a subclass of `TypeExpression`.
        """
        return TypeExpression(self.name, fields=fields)

    def iter_symbols(self):
        yield self.name

    def is_concrete(self):
        return False


class _Subtypeable(_TypeBase):
    def __iter__(self):
        yield self

    def __le__(self, other):
        r = self._is_subtype_(other)
        if r is NotImplemented:
            r = other._is_supertype_(self)
            if r is NotImplemented:
                return False
        return r

    def __ge__(self, other):
        r = self._is_supertype_(other)
        if r is NotImplemented:
            r = other._is_subtype_(self)
            if r is NotImplemented:
                return False
        return r

    def _aug_is_subtype(self, other):
        r = self._is_subtype_(other)
        if r is NotImplemented:
            r = other._is_supertype_(self)
            if r is NotImplemented:
                return False
        return r

    def _is_subtype_(self, other):
        raise NotImplementedError('_is_subtype_')

    def _is_supertype_(self, other):
        raise NotImplementedError('_is_supertype_')

    def __or__(self, other):
        self._aug_validate_union(other)

        if self >= other:
            return self
        if self <= other:
            return other

        _, greatest_common_subtypes = common_subtypes(*self, *other)
        return self._build_union_(*greatest_common_subtypes)

    def _aug_validate_union(self, other):
        self._validate_union_(other)
        if isinstance(other, _Subtypeable):
            other._validate_union_(self)

    def _validate_union_(self, other):
        raise NotImplementedError('_validate_union_')

    def _build_union_(self, *members):
        raise NotImplementedError('_build_union_')

    def __and__(self, other):
        self._aug_validate_intersection(other)

        if self >= other:
            return other
        if self <= other:
            return self

        t = self.get_bottom()
        if isinstance(other, type(t)) or isinstance(self, type(t)):
            for s, o in itertools.product(self, other):
                t |= s & o
            return t

        try:
            IntersectionType = type(self.get_top())
        except NotImplementedError:
            return self._build_intersection_(self, other)

        this = self
        # "unpack" intersections as their "iter" yields only themselves, not
        # their components
        if isinstance(self, IntersectionType):
            self = self.members
        if isinstance(other, IntersectionType):
            other = other.members
        least_common_subtypes, _ = common_subtypes(*self, *other)
        return this._build_intersection_(*least_common_subtypes)

    def _aug_validate_intersection(self, other):
        self._validate_intersection_(other)
        if isinstance(other, _Subtypeable):
            other._validate_intersection_(self)

    def _validate_intersection_(self, other):
        raise NotImplementedError('_validate_intersection_')

    def _build_intersection_(self, *members):
        raise NotImplementedError('_build_intersection_')

    def get_bottom(self):
        return self._build_union_()

    def is_bottom(self):
        return self.equals(self.get_bottom())

    def get_top(self):
        return self._build_intersection_()

    def is_top(self):
        return self.equals(self.get_top())

    def equals(self, other):
        return self <= other <= self


class TypeExpression(_Subtypeable):
    def __init__(self, name, fields=(), predicate=None):
        self.name = name
        self.fields = fields

        if predicate is not None and predicate.is_top():
            predicate = None
        self.predicate = predicate

        self._freeze_()

    @property
    def evaled_predicate(self):
        if self.predicate is None:
            return IntersectionPredicate()
        return self.predicate

    def __hash__(self):
        return (hash(self.__class__.__name__) ^
                hash(self.name) ^
                hash(self.predicate) ^
                hash(self.fields))

    def __eq__(self, other):
        # Deep equality, but not semantic equality.
        if type(self) is not type(other):
            return NotImplemented
        return (self.name == other.name and
                self.predicate == other.predicate and
                self.fields == other.fields)

    def __repr__(self):
        result = self.name
        if self.fields:
            result += '[%s]' % ', '.join(repr(f) for f in self.fields)
        if self.predicate:
            result += ' %% %r' % self.predicate
        return result

    def __getitem__(self, fields):
        raise TypeError("%r has no empty fields (not subscriptable)." % self)

    def _apply_fields_(self, fields):
        return self.__class__(self.name, fields=fields,
                              predicate=self.predicate)

    def __contains__(self, value):
        return self._is_element_(value) and value in self.evaled_predicate

    def _is_element_(self, value):
        return False

    def __mod__(self, predicate):
        if self.predicate:
            raise TypeError("%r already has a predicate." % self)
        if predicate is None:
            return self

        self._validate_predicate_(predicate)
        return self._apply_predicate_(predicate=predicate)

    def _validate_predicate_(self, predicate):
        if not isinstance(predicate, Predicate):
            raise TypeError("%r is not a predicate." % predicate)

    def _apply_predicate_(self, predicate):
        if isinstance(predicate, type(predicate.get_bottom())):
            return self._build_union_(*(
                self.__class__(self.name, fields=self.fields, predicate=p)
                for p in predicate.members))
        return self.__class__(self.name, fields=self.fields,
                              predicate=predicate)

    def _validate_union_(self, other):
        if not isinstance(other, TypeExpression):
            if isinstance(other, CompositeType):
                raise TypeError("Cannot union an incomplete type %r with %r."
                                % (other, self))
            else:
                raise TypeError("%r is not a type expression." % other)

    def _build_union_(self, *members):
        return UnionTypeExpression(*members)

    def _validate_intersection_(self, other):
        if not isinstance(other, TypeExpression):
            if isinstance(other, CompositeType):
                raise TypeError("Cannot intersect an incomplete type %r with"
                                " %r." % (other, self))
            else:
                raise TypeError("%r is not a type expression." % other)

    def _build_intersection_(self, *members):
        try:
            self, other = members
        except ValueError:
            raise NotImplementedError

        if self.name != other.name:
            return self.get_bottom()

        new_fields = tuple(s & o for s, o in itertools.zip_longest(
            self.fields, other.fields, fillvalue=self.get_bottom()))
        if any(f.is_bottom() for f in new_fields):
            return self.get_bottom()

        new_predicate = self.evaled_predicate & other.evaled_predicate
        if new_predicate.is_bottom():
            return self.get_bottom()

        return self._apply_fields_(new_fields)._apply_predicate_(new_predicate)

    def _is_subtype_(self, other):
        if not hasattr(other, 'name') or self.name != other.name:
            return NotImplemented

        for f1, f2 in itertools.zip_longest(
                self.fields, other.fields, fillvalue=self.get_bottom()):
            if not (f1 <= f2):
                return False

        if not (self.evaled_predicate <= other.evaled_predicate):
            return False

        return True

    def _is_supertype_(self, other):
        return NotImplemented

    def iter_atomics(self):
        seen = set()
        for elem in self:
            for fields in itertools.product(*elem.fields):
                for predicate in elem.evaled_predicate:
                    base = elem._apply_fields_(fields)
                    refined = base._apply_predicate_(predicate)
                    if refined not in seen:
                        seen.add(refined)
                        yield refined

    def is_concrete(self):
        """Deprecated, use is_atomic"""
        return self.is_atomic()

    def is_atomic(self):
        return len(list(self.iter_atomics())) == 1

    def iter_symbols(self):
        yield self.name
        for field in self.fields:
            yield from field.iter_symbols()

    def to_ast(self):
        return {
            "type": 'expression',
            "name": self.name,
            "predicate": self.predicate.to_ast() if self.predicate else {},
            "fields": [field.to_ast() for field in self.fields]
        }


class UnionTypeExpression(TypeExpression):
    def __init__(self, *members):
        # Constructing a Union directly will allow redundant/subtype members
        # this is useful for structure types
        self.members = members
        super().__init__(name='')

    def __hash__(self):
        return super().__hash__() ^ hash(self.members)

    def __eq__(self, other):
        return (super().__eq__(other)
                and set(self.members) == set(other.members))

    def __iter__(self):
        yield from self.members

    def __repr__(self):
        return " | ".join(repr(m) for m in self.members)

    def _build_union_(self, *members):
        return self.__class__(*members)

    def _is_element_(self, value):
        return any(value in s for s in self.members)

    def _is_subtype_(self, other):
        if isinstance(other, self.__class__):
            return all(any(s <= o for o in other.members)
                       for s in self.members)
        return all(s <= other for s in self.members)

    def _is_supertype_(self, other):
        if isinstance(other, self.__class__):
            return all(any(s >= o for s in self.members)
                       for o in other.members)
        return any(s >= other for s in self.members)

    def _validate_predicate_(self, predicate):
        raise TypeError("Cannot apply predicates to union types.")

    def to_ast(self):
        return {
            'type': 'union',
            'members': [m.to_ast() for m in self.members]
        }


class Predicate(_Subtypeable):
    def __init__(self, *args, **kwargs):
        truthy = any(map(bool, args)) or any(map(bool, kwargs.values()))
        if not truthy:
            raise TypeError("Predicate %r has no arguments or cannot "
                            "constrain the type." % self.__class__.__name__)

        self._freeze_()

    def __hash__(self):
        # This trivially satisfies the property:
        # x == x => hash(x) == hash(x)
        # Subclasses ought to override this with something less... collision-y.
        return 0

    def __eq__(self, other):
        return NotImplemented

    def __contains__(self, value):
        return self._is_element_(value)

    def _is_element_(self, value):
        raise NotImplementedError('_is_element_')

    def _is_subtype_(self, other):
        return NotImplemented

    def _is_supertype_(self, other):
        return NotImplemented

    def _validate_union_(self, other):
        if not isinstance(other, Predicate):
            raise TypeError("Cannot union %r and %r" % (self, other))

    def _build_union_(self, *members):
        return UnionPredicate(*members)

    def _validate_intersection_(self, other):
        if not isinstance(other, Predicate):
            raise TypeError("Cannot intersect %r and %r" % (self, other))

    def _build_intersection_(self, *members):
        return IntersectionPredicate(*members)

    def to_ast(self):
        return {
            'type': 'predicate',
            'name': self.__class__.__name__
        }


class _IdentityPredicateBase(Predicate):
    _operator = '?'

    def __init__(self, *members):
        self.members = members
        self._freeze_()

    def __hash__(self):
        return super().__hash__() ^ hash(self.members)

    def __eq__(self, other):
        return super().__eq__(other) and self.members == other.members

    def __repr__(self):
        if not self.members:
            return ""
        return "(%s)" % (" %s " % self._operator).join(
            repr(m) for m in self.members)


class UnionPredicate(_IdentityPredicateBase):
    _operator = '|'

    def __iter__(self):
        yield from self.members

    def _build_union_(self, *members):
        return self.__class__(*members)

    def _is_element_(self, value):
        return any(value in s for s in self.members)

    def _is_subtype_(self, other):
        if isinstance(other, self.__class__):
            return all(any(s <= o for o in other.members)
                       for s in self.members)
        return all(s <= other for s in self.members)

    def _is_supertype_(self, other):
        if isinstance(other, self.__class__):
            return all(any(s >= o for s in self.members)
                       for o in other.members)
        return any(s >= other for s in self.members)

    def to_ast(self):
        return {
            'type': 'union',
            'members': [m.to_ast() for m in self.members]
        }


class IntersectionPredicate(_IdentityPredicateBase):
    _operator = '&'

    def _build_intersection_(self, *members):
        return self.__class__(*members)

    def _is_element_(self, value):
        return all(value in s for s in self.members)

    def _is_subtype_(self, other):
        if isinstance(other, self.__class__):
            return all(any(s <= o for s in self.members)
                       for o in other.members)
        if isinstance(other, type(other.get_bottom())):
            return other >= self
        return any(m <= other for m in self.members)

    def _is_supertype_(self, other):
        if isinstance(other, self.__class__):
            return all(any(s >= o for o in other.members)
                       for s in self.members)
        if isinstance(other, type(other.get_bottom())):
            return other <= self
        return all(m >= other for m in self.members)

    def to_ast(self):
        return {
            'type': 'intersection',
            'members': [m.to_ast() for m in self.members]
        }


class TestPredicate(Predicate):
    def __init__(self, name):
        self.name = name

        self._freeze_()

    def __repr__(self):
        return self.name

    def _is_subtype_(self, other):
        if isinstance(other, self.__class__):
            return self.name == other.name
        return NotImplemented

    def _is_supertype_(self, other):
        if isinstance(other, self.__class__):
            return self.name == other.name
        return NotImplemented
