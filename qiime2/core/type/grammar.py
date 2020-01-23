# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import itertools
from abc import ABCMeta, abstractmethod

from qiime2.core.util import tuplize, ImmutableBase


def maximal_antichain(*types):
    maximal_elements = {}  # easy to delete, retains order
    for t in types:
        placed = False
        for e in list(maximal_elements):
            if e <= t:
                # Delete first! Duplicate keys would disappear otherwise
                del maximal_elements[e]
                maximal_elements[t] = None
                placed = True
        if not placed:
            maximal_elements[t] = None
    return tuple(maximal_elements)


def minimal_antichain(*types):
    minimal_elements = {}  # easy to delete, retains order
    for t in types:
        placed = False
        for e in list(minimal_elements):
            if t <= e:
                # Delete first! Duplicate keys would disappear otherwise
                del minimal_elements[e]
                minimal_elements[t] = None
                placed = True
        if not placed:
            minimal_elements[t] = None
    return tuple(minimal_elements)


class _ExpBase(ImmutableBase, metaclass=ABCMeta):
    def __init__(self, template):
        # Super basic smoke-test
        assert template is None or template.is_template
        self.template = template

    def __getattr__(self, name):
        if ('template' in self.__dict__
                and self.template is not None
                and name in self.template.public_proxy):
            return getattr(self.template, name)

        raise AttributeError("%r object has no attribute %r"
                             % (type(self), name))

    # Prevent infinite recursion when pickling due to __getattr__
    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, state):
        self.__dict__ = state

    @property
    def name(self):
        return self.template.get_name_expr(self)

    @property
    def kind(self):
        return self.template.get_kind_expr(self)

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

    def iter_symbols(self):
        yield self.name


class IncompleteExp(_ExpBase):
    def __init__(self, template):
        super().__init__(template)

        if (self.template is None
                or not list(self.template.get_field_names_expr(self))):
            raise ValueError("Template %r has no fields, should not be used"
                             " with a IncompleteExp." % (template,))

    def __eq__(self, other):
        if type(self) is not type(other):
            return NotImplemented

        return (self.name == other.name
                and tuple(self.template.get_field_names_expr(self))
                == tuple(other.template.get_field_names_expr(self)))

    def __hash__(self):
        return (hash(type(self))
                ^ hash(self.name)
                ^ hash(tuple(self.template.get_field_names_expr(self))))

    def __repr__(self):
        fields = ', '.join(
            '{%s}' % f for f in self.template.get_field_names_expr(self))
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
        for field in fields:
            if not isinstance(field, _AlgebraicExpBase):
                raise TypeError("Field %r is not complete type expression."
                                % (field,))
        self.template.validate_fields_expr(self, fields)
        return TypeExp(self.template, fields=fields)

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
        if not ((self.is_bottom() or other.is_bottom())
                or (self.get_union_membership() == other.get_union_membership()
                    and self.get_union_membership() is not None)):
            raise TypeError("Cannot union %r and %r" % (self, other))

        if self >= other:
            return self
        if self <= other:
            return other

        union = UnionExp((*self.unpack_union(), *other.unpack_union()))
        return union.normalize()

    def __and__(self, other):
        if (not self.can_intersect() or not other.can_intersect()
                or (self.kind != other.kind
                    and not (self.is_top() or other.is_top()))):
            raise TypeError("Cannot intersect %r and %r" % (self, other))

        # inverse of __or__
        if self >= other:
            return other
        if self <= other:
            return self

        # Distribute over union
        if isinstance(self, UnionExp) or isinstance(other, UnionExp):
            m = []
            for s, o in itertools.product(self.unpack_union(),
                                          other.unpack_union()):
                m.append(s & o)
            return UnionExp(m).normalize()

        elements = list(itertools.chain(self.unpack_intersection(),
                                        other.unpack_intersection()))

        if len(elements) > 1:
            # Give the expression a chance to collapse, as many intersections
            # are contradictions
            collapse = elements[0]._collapse_intersection_(elements[1])
            if collapse is not NotImplemented:
                for e in elements[2:]:
                    collapse = collapse._collapse_intersection_(e)

                return collapse

        # Back to the regularly scheduled inverse of __or__
        members = minimal_antichain(*self.unpack_intersection(),
                                    *other.unpack_intersection())
        return IntersectionExp(members)

    def _collapse_intersection_(self, other):
        return NotImplemented

    def equals(self, other):
        return self <= other <= self

    def is_concrete(self):
        return len(list(self.unpack_union())) == 1

    def get_union_membership(self):
        if self.template is not None:
            return self.template.get_union_membership_expr(self)

        return True

    def can_intersect(self):
        return True

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
        if predicate is not None and predicate.is_top():
            predicate = None

        self.fields = tuple(fields)
        self.predicate = predicate

        super()._freeze_()

    @property
    def full_predicate(self):
        if self.predicate is None:
            return IntersectionExp()
        return self.predicate

    def __eq__(self, other):
        if type(self) is not type(other):
            return NotImplemented

        return (self.kind == other.kind
                and self.name == other.name
                and self.fields == other.fields
                and self.full_predicate == other.full_predicate)

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
            if self.predicate.template is None:  # is _IdentityExpBase
                predicate = '(%s)' % predicate
            result += ' % ' + predicate
        return result

    def __getitem__(self, fields):
        raise TypeError("Cannot apply fields (%r) to %r,"
                        " fields already present." % (fields, self))

    def __contains__(self, value):
        return (self.template.is_element_expr(self, value)
                and value in self.full_predicate)

    def __iter__(self):
        yield from {self.duplicate(fields=fields)
                    for fields in itertools.product(*self.fields)}

    def iter_symbols(self):
        yield self.name
        for field in self.fields:
            yield from field.iter_symbols()

    def _is_subtype_(self, other):
        if other.template is None:
            return NotImplemented

        if not self.template.is_symbol_subtype_expr(self, other):
            return False
        for f1, f2 in itertools.zip_longest(self.fields, other.fields,
                                            # more fields = more specific
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
        if predicate is None or predicate.is_top():
            return self

        return self.duplicate(predicate=predicate)

    def __rmod__(self, other):
        raise TypeError("Predicate (%r) must be applied to the right-hand side"
                        " of a type expression." % (other,))

    def duplicate(self, fields=(), predicate=None):
        if fields == ():
            fields = self.fields
        else:
            self.template.validate_fields_expr(self, fields)

        if predicate is None:
            predicate = self.predicate
        elif predicate.is_top():
            predicate = None
        elif predicate.template is not None:
            self.template.validate_predicate_expr(self, predicate)

        return self.__class__(self.template, fields=fields,
                              predicate=predicate)

    def _collapse_intersection_(self, other):
        if self.name != other.name:
            return UnionExp()

        new_fields = tuple(
            s & o for s, o in itertools.zip_longest(self.fields, other.fields,
                                                    # same as a type mismatch
                                                    fillvalue=UnionExp()))
        if any(f.is_bottom() for f in new_fields):
            return UnionExp()

        new_predicate = self.full_predicate & other.full_predicate
        if new_predicate.is_bottom():
            return UnionExp()

        return self.duplicate(fields=new_fields, predicate=new_predicate)

    def is_concrete(self):
        return self._bool_attr_method('is_concrete')

    def _bool_attr_method(self, method_name):
        def method(s): return getattr(s, method_name)()

        if any(not method(f) for f in self.fields):
            return False
        if not method(self.full_predicate):
            return False

        return True

    def to_ast(self):
        ast = {
            "type": "expression",
            "builtin": True,
            "name": self.name,
            "predicate": self.predicate.to_ast() if self.predicate else None,
            "fields": [field.to_ast() for field in self.fields]
        }
        self.template.update_ast_expr(self, ast)
        return ast


class PredicateExp(_AlgebraicExpBase):
    def __init__(self, template):
        super().__init__(template)
        super()._freeze_()

    def __eq__(self, other):
        return self.template == other.template

    def __hash__(self):
        return hash(self.template)

    def __contains__(self, value):
        return self.template.is_element_expr(self, value)

    def __repr__(self):
        return repr(self.template)

    def _is_subtype_(self, other):
        if (other.template is not None
                and self.template.is_symbol_subtype_expr(self, other)):
            return True

        return NotImplemented

    def _is_supertype_(self, other):
        if (other.template is not None
                and self.template.is_symbol_supertype_expr(self, other)):
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

    def to_ast(self):
        ast = {
            "type": "predicate",
            "name": self.name,
        }
        self.template.update_ast_expr(self, ast)
        return ast


class _IdentityExpBase(_AlgebraicExpBase):
    """
    Base class for IntersectionExp and UnionExp.

    If there are no members, then they are Top or Bottom types respectively
    and represent identity values (like 1 for mul and 0 for add) for the
    type algebra.

    There is no template object for these expressions. That property will
    always be `None`.
    """
    _operator = ' ? '

    def __init__(self, members=()):
        super().__init__(template=None)
        self.members = tuple(members)

        super()._freeze_()

    @property
    def kind(self):
        if not self.members:
            return "identity"
        return self.members[0].kind

    @property
    def name(self):
        return ""

    def __eq__(self, other):
        return (type(self) is type(other)
                and set(self.members) == set(other.members))

    def __hash__(self):
        return hash(type(self)) ^ hash(frozenset(self.members))

    def __repr__(self):
        if not self.members:
            return self.__class__.__name__ + "()"
        return self._operator.join(repr(m) for m in self.members)

    def __iter__(self):
        for m in self.unpack_union():
            yield from m

    def iter_symbols(self):
        for m in self.unpack_union():
            yield from m.iter_symbols()

    def get_union_membership(self):
        if self.members:
            return self.members[0].get_union_membership()


class UnionExp(_IdentityExpBase):
    _operator = ' | '  # used by _IdentityExpBase.__repr__

    def __contains__(self, value):
        return any(value in s for s in self.members)

    def _is_subtype_(self, other):
        if (isinstance(other, self.__class__)
                and type(other) is not self.__class__):  # other is subclass
            return NotImplemented

        # if other isn't a union, becomes all(s <= other for s in self.members)
        return all(any(s <= o for o in other.unpack_union())
                   for s in self.unpack_union())

    def _is_supertype_(self, other):
        return all(any(s >= o for s in self.unpack_union())
                   for o in other.unpack_union())

    def is_bottom(self):
        return not self.members

    def unpack_union(self):
        yield from self.members

    def to_ast(self):
        return {
            "type": "union",
            "members": [m.to_ast() for m in self.members]
        }

    def normalize(self):
        elements = self.members

        groups = {}
        for e in elements:
            if type(e) is TypeExp:
                candidate = e.duplicate(predicate=IntersectionExp())
                if candidate in groups:
                    groups[candidate].append(e)
                else:
                    groups[candidate] = [e]
            else:
                # groups should be empty already, but don't even attempt
                # collapsing if its a union of type expressions and "other"
                groups = {}
                break

        if groups:
            elements = []
            for canidate, group in groups.items():
                if len(group) == 1:
                    elements.append(group[0])
                else:
                    predicate = UnionExp([t.full_predicate for t in group])
                    predicate = predicate.normalize()
                    elements.append(candidate.duplicate(predicate=predicate))

        if len(elements) < 20:
            members = maximal_antichain(*elements)
        else:
            members = elements

        if len(members) == 1:
            return members[0]
        return UnionExp(members)


class IntersectionExp(_IdentityExpBase):
    _operator = ' & '  # used by _IdentityExpBase.__repr__

    def __contains__(self, value):
        return all(value in s for s in self.members)

    def _is_subtype_(self, other):
        if isinstance(other, UnionExp):
            # Union will treat `self` as an atomic type, comparing
            # its elements against `self`. This in turn will recurse back to
            # `self` allowing it to check if it is a subtype of the union
            # elements. That check will ultimately compare the elements of
            # `self` against a single element of the union.
            return NotImplemented

        return all(any(s <= o for s in self.unpack_intersection())
                   for o in other.unpack_intersection())

    def _is_supertype_(self, other):
        if isinstance(other, UnionExp):
            return NotImplemented

        return all(any(s >= o for o in other.unpack_intersection())
                   for s in self.unpack_intersection())

    def is_top(self):
        return not self.members

    def unpack_intersection(self):
        yield from self.members

    def to_ast(self):
        return {
            "type": "intersection",
            "members": [m.to_ast() for m in self.members]
        }
