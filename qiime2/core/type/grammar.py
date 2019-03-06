# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import itertools
from abc import ABCMeta, abstractmethod

from qiime2.core.util import tuplize
from ..util import ImmutableBase


def least_common_subtypes(*types):
    least = []
    for t in types:
        for idx, partition_rep in enumerate(least):
            if t <= partition_rep:
                least[idx] = t
                break
        else:  # evaled if the above loop did not break
            least.append(t)

    return least


def greatest_common_subtypes(*types):
    greatest = []
    for t in types:
        for idx, partition_rep in enumerate(greatest):
            if t >= partition_rep:
                greatest[idx] = t
                break
        else:  # evaled if the above loop did not break
            greatest.append(t)

    return greatest


class _ExpBase(metaclass=ABCMeta):
    def __init__(self, template):
        self.template = template

    def __getattr__(self, name):
        if name in self.template.public_proxy:
            return getattr(self.template, name)

        raise AttributeError("%r object has no attribute %r"
                             % (type(self), name))

    @property
    def name(self):
        return self.template.get_name()

    @property
    def kind(self):
        return self.template.get_kind()

    @abstractmethod
    def __eq__(self, other):
        raise NotImplementedError

    def __ne__(self, other):
        return not self == other

    @abstractmethod
    def __le__(self, other):
        raise NotImplementedError

    @abstractmethod
    def __ge__(self, other):
        raise NotImplementedError

    @abstractmethod
    def __or__(self, other):
        raise NotImplementedError

    def __ror__(self, other):
        return self | other

    @abstractmethod
    def __and__(self, other):
        raise NotImplementedError

    def __rand__(self, other):
        return self & other

    @abstractmethod
    def equals(self, other):
        raise NotImplementedError

    def is_concrete(self):
        return False


class IncompleteExp(_ExpBase):
    def __eq__(self, other):
        return type(self) is type(other) and self.name == other.name

    def __hash__(self):
        return hash(type(self) ^ hash(self.name))

    def __repr__(self):
        fields = ', '.join('{%s}' % f for f in self.template.get_field_names())
        return self.name + ('[%s]' % fields)

    def __le__(self, other):
        raise TypeError("Cannot compare subtype, %r is missing arguments"
                        " for its fields." % (self,))

    def __ge__(self, other):
        raise TypeError("Cannot compare supertype, %r is missing arguments"
                        " for its fields." % (self,))

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
        self.template.validate_fields(tuple(f.template for f in fields),
                                      fields)
        template = self.template.specialize(fields)
        return TypeExp(template, fields=fields)

    def equals(self, other):
        return self == other


class _AlgebraicExpBase(_ExpBase):
    def __le__(self, other):
        first = self._is_subtype_(other)
        if first is not NotImplemented:
            return first

        second = other._is_supertype_(self)
        if second is not NotImplemented:
            return second

        return False

    def __ge__(self, other):
        first = self._is_supertype_(other)
        if first is not NotImplemented:
            return first

        second = other._is_subtype_(self)
        if second is not NotImplemented:
            return second

        return False

    def __or__(self, other):
        self.template.validate_union(other.template)
        other.template.validate_union(self.template)

        if self >= other:
            return self
        if self <= other:
            return other

        members = greatest_common_subtypes(*self.unpack_union(),
                                           *other.unpack_union())
        return UnionExp(members, template=self.template)

    def __and__(self, other):
        self.template.validate_intersection(other.template)
        other.template.validate_intersection(self.template)

        # inverse of __or__
        if self >= other:
            return other
        if self <= other:
            return self

        # Distribute over union
        t = UnionExp(template=self.template)
        if isinstance(self, UnionExp) or isinstance(other, UnionExp):
            for s, o in itertools.product(self.unpack_union(),
                                          other.unpack_union()):
                t |= s & o
            return t

        # Give the expression a chance to collapse, as many intersections
        # are contradictions
        collapse = self._collapse_intersection_(other)
        if collapse is not NotImplemented:
            return collapse

        # Back to the regularly scheduled inverse of __or__
        members = least_common_subtypes(*self.unpack_intersection(),
                                        *other.unpack_intersection())
        return IntersectionExp(members, template=self.template)

    def _collapse_intersection_(self, other):
        return NotImplemented

    def equals(self, other):
        return self <= other <= self

    def is_concrete(self):
        return len(list(self.unpack_union())) == 1

    # These methods are to be overridden by UnionExp
    def is_bottom(self):
        return False

    def unpack_union(self):
        yield self

    # These methods are to be overridden by IntersectionExp
    def is_top(self):
        return False

    def unpack_intersection(self):
        yield self


class TypeExp(_AlgebraicExpBase):
    def __init__(self, template, fields=(), predicate=None):
        super().__init__(template)
        self.fields = tuple(fields)
        self.predicate = predicate


    @property
    def full_predicate(self):
        if self.predicate is None:
            return IntersectionExp()
        return self.predicate

    def __eq__(self, other):
        return (type(self) is type(other)
                and self.kind == other.kind
                and self.name == other.name
                and self.fields == other.fields
                and self.predicate == other.predicate)

    def __hash__(self):
        return (hash(type(self))
                ^ hash(self.kind) ^ hash(self.name)
                ^ hash(self.fields) ^ hash(self.predicate))

    def __repr__(self):
        result = self.name
        if self.fields:
            result += '[%s]' % ', '.join(repr(f) for f in self.fields)
        if self.predicate:
            predicate = repr(self.predicate)
            if isinstance(predicate, (IntersectionExp, UnionExp)):
                predicate = '(%s)' % predicate
            result += ' % ' + predicate
        return result

    def __contains__(self, value):
        return (self.template.is_element(value, self)
                and value in self.full_predicate)

    def _is_subtype_(self, other):
        if not self.template.is_subtype(other.template):
            return NotImplemented

        for f1, f2 in itertools.zip_longest(self.fields, other.fields,
                                            fillvalue=IntersectionExp()):
            if not (f1 <= f2):
                return False

        if not (self.full_predicate <= other.full_predicate):
            return False

        return True

    def _is_supertype_(self, other):
        return NotImplemented

    def __mod__(self, predicate):
        if self.predicate:
            raise TypeError("%r already has a predicate, will not add %r"
                            % (self, predicate))
        if predicate is None:
            return self

        if isinstance(predicate, UnionExp):
            return UnionExp((self.duplicate(predicate=p)
                             for p in predicate.members),
                            template=self.template)

        return self.duplicate(predicate=predicate)

    def __rmod__(self, other):
        raise TypeError("Predicate (%r) must be applied to the right-hand side"
                        " of a type expression." % (other,))

    def duplicate(self, fields=(), predicate=None):
        if fields == ():
            fields = self.fields
        else:
            self.template.validate_fields((f.template for f in fields), self)

        if predicate is None:
            predicate = self.predicate
        else:
            self.template.validate_predicate(predicate.template, self)

        return self.__class__(self.template, fields=fields,
                              predicate=predicate)

    def _collapse_intersection_(self, other):
        if self.name != other.name:
            return UnionExp(template=self.template)

        new_fields = tuple(s & o for s, o in itertools.zip_longest(
            self.fields, other.fields,
            fillvalue=UnionExp(template=self.template)))
        if any(f.is_bottom() for f in new_fields):
            return UnionExp(template=self.template)

        new_predicate = self.full_predicate & other.full_predicate
        if new_predicate.is_bottom():
            return UnionExp(template=self.template)

        return self.duplicate(fields=new_fields, predicate=new_predicate)

    def is_concrete(self):
        if not super().is_concrete():
            return False

        if not self.full_predicate.is_concrete():
            return False

        if any(not f.is_concrete() for f in self.fields):
            return False

        return True


class PredicateExp(_AlgebraicExpBase):
    def __eq__(self, other):
        return self.template == other.template

    def __hash__(self):
        return hash(self.template)

    def __contains__(self, value):
        return self.template.is_element(value, self)

    def __repr__(self):
        return repr(self.template)

    def _is_subtype_(self, other):
        if self.template.is_subtype(other.template):
            return True

        return NotImplemented

    def _is_supertype_(self, other):
        if self.template.is_supertype(other.template):
            return True

        return NotImplemented

    def _collapse_intersection_(self, other):
        first = self.template.collapse_intersection(other.template)
        if first is None:
            return UnionExp()
        elif first is not NotImplemented:
            return self.__class__(first)

        second = other.template.collapse_intersection(self.template)
        if second is None:
            return UnionExp()
        elif second is not NotImplemented:
            return self.__class__(second)

        return NotImplemented


class _IdentityExpBase(_AlgebraicExpBase):
    _operator = ' ? '

    @property
    def kind(self):
        if self.template is not None:
            return super().kind

        return "<identity>"

    @property
    def name(self):
        return ""

    def __init__(self, members=(), template=None):
        super().__init__(template)
        self.members = tuple(members)

    def __eq__(self, other):
        return (type(self) is type(other)
                and set(self.members) == set(other.members))

    def __hash__(self):
        return hash(type(self)) ^ hash(frozenset(self.members))

    def __repr__(self):
        if not self.members:
            return self.__class__.__name__ + "()"

        return self._operator.join(repr(m) for m in self.members)

    def _is_subtype_(self, other):
        return self <= other

    def _is_supertype_(self, other):
        return self >= other


class UnionExp(_IdentityExpBase):
    _operator = ' | '

    def __contains__(self, value):
        return any(value in s for s in self.members)

    def __le__(self, other):
        if isinstance(other, UnionExp):
            return all(any(s <= o for o in other.members)
                       for s in self.members)
        return all(s <= other for s in self.members)

    def __ge__(self, other):
        if isinstance(other, self.__class__):
            return all(any(s >= o for s in self.members)
                       for o in other.members)
        return any(s >= other for s in self.members)

    def is_bottom(self):
        return not self.members

    def unpack_union(self):
        yield from self.members


class IntersectionExp(_IdentityExpBase):
    _operator = ' & '

    def __contains__(self, value):
        return all(value in s for s in self.members)

    def __le__(self, other):
        if isinstance(other, UnionExp):
            return other >= self

        if isinstance(other, IntersectionExp):
            return all(any(s <= o for s in self.members)
                       for o in other.members)
        return any(m <= other for m in self.members)

    def __ge__(self, other):
        if isinstance(other, UnionExp):
            return other <= self

        if isinstance(other, IntersectionExp):
            return all(any(s >= o for o in other.members)
                       for s in self.members)
        return all(m >= other for m in self.members)

    def is_top(self):
        return not self.members

    def unpack_intersection(self):
        yield from self.members
