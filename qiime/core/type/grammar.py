# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import types
import itertools

from qiime.core.util import tuplize


def _immutable_error(*args):
    raise TypeError("Types are immutable.")


class _ImmutableBase:
    """Disables instance mutation (but not static/class mutation) using
    `_freeze_` method. This class also implements shared methods."""

    def _freeze_(self):
        """Disables __setattr__ when called. It is idempotent."""
        self._frozen = True  # The particular value doesn't matter

    __delattr__ = __setitem__ = __delitem__ = _immutable_error

    def __setattr__(self, *args):
        # This doesn't stop silly things like
        # object.__setattr__(obj, ...), but that's a pretty rude thing
        # to do anyways. We are just trying to avoid accidental mutation.
        if hasattr(self, '_frozen'):
            _immutable_error()
        super().__setattr__(*args)

    def __ne__(self, other):
        return not self == other

    def __rmod__(self, predicate):
        raise TypeError("Predicate must be applied to the right-hand side of"
                        " a type expression.")

    def __ror__(self, other):
        return self | other  # union should be associative

    def __rand__(self, other):
        return self & other  # intersection should be associative


class CompositeType(_ImmutableBase):
    def __init__(self, name, field_names):
        # These classes aren't user-facing, but some light validation avoids
        # accidental issues. However, we don't want to waste a lot of time with
        # validation here, validation should happen elsewhere.
        if not len(field_names):
            raise ValueError("`field_names` cannot be an empty array-like.")

        self.name = name
        self.field_names = field_names

        self._freeze_()

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


class TypeExpression(_ImmutableBase):
    def __init__(self, name, fields=(), predicate=None):
        self.name = name
        self.predicate = predicate
        self.fields = fields

        self._freeze_()

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

    def equals(self, other):
        # Different from __eq__ which has to match hashing but can't
        # consider semantic equality
        return self <= other <= self

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
        return (self._is_element_(value) and
                ((not self.predicate) or value in self.predicate))

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
        return self.__class__(self.name, fields=self.fields,
                              predicate=predicate)

    def __or__(self, other):
        self._validate_union_(other, handshake=False)
        if self == other:
            return self
        return self._build_union_((self, other))

    def _validate_union_(self, other, handshake=False):
        if not isinstance(other, TypeExpression):
            if isinstance(other, CompositeType):
                raise TypeError("Cannot union an incomplete type %r with %r."
                                % (other, self))
            else:
                raise TypeError("%r is not a type expression." % other)

        if not handshake:
            other._validate_union_(self, handshake=True)

    def _build_union_(self, members):
        return UnionTypeExpression(members)

    def __and__(self, other):
        self._validate_intersection_(other, handshake=False)
        if self == other:
            return other
        return self._build_intersection_((self, other))

    def _validate_intersection_(self, other, handshake=False):
        if not isinstance(other, TypeExpression):
            if isinstance(other, CompositeType):
                raise TypeError("Cannot intersect an incomplete type %r with"
                                " %r." % (other, self))
            else:
                raise TypeError("%r is not a type expression." % other)

        if not handshake:
            other._validate_intersection_(self, handshake=True)

    def _build_intersection_(self, members):
        return IntersectionTypeExpression(members)

    def __le__(self, other):
        return all(any(s._aug_is_subtype(o) for o in other) for s in self)

    def __ge__(self, other):
        return all(any(o._aug_is_subtype(s) for s in self) for o in other)

    def _aug_is_subtype(self, other):
        r = self._is_subtype_(other)
        if r is NotImplemented:
            return other._is_supertype_(self)
        return r

    def _is_subtype_(self, other):
        if self.name != other.name:
            return False
        for f1, f2 in itertools.zip_longest(self.fields, other.fields):
            if not (f1 <= f2):
                return False
        if other.predicate and not self.predicate <= other.predicate:
            return False
        return True

    def _is_supertype_(self, other):
        # Invoked only when `other`'s `_is_subtype_` returned `NotImplemented`
        # that really shouldn't be needed most of the time.
        raise NotImplementedError

    def __iter__(self):
        yield from set(self._apply_fields_(fields=fields)
                       for fields in itertools.product(*self.fields))

    def is_concrete(self):
        return len(list(self)) == 1

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


class _SetOperationBase(TypeExpression):
    _operator = '?'  # Used for repr only - ? chosen as it is not a Python op.

    def __init__(self, members):
        m = []
        for member in members:
            # We can flatten the object a little, which will avoid excessive
            # recursion (it would look like a cons-list otherwise)
            if type(member) is type(self):
                m.extend(member.members)
            else:
                m.append(member)

        self.members = frozenset(m)

        super().__init__('')  # Unions/intersections do not have a name

    def __hash__(self):
        return super().__hash__() ^ hash(self.members)

    def __eq__(self, other):
        super_eq = super().__eq__(other)
        if super_eq is NotImplemented:
            return NotImplemented
        return super_eq and self.members == other.members

    def __repr__(self):
        return (" %s " % self._operator) \
            .join(sorted([repr(m) for m in self.members]))

    def _validate_predicate_(self, predicate):
        raise TypeError("Cannot apply predicates to union/intersection types.")

    def to_ast(self):
        return {
            'members': [m.to_ast() for m in self.members]
        }

    def __iter__(self):
        yield from set(itertools.chain.from_iterable(self.members))


class UnionTypeExpression(_SetOperationBase):
    _operator = '|'

    def _validate_intersection_(self, other, handshake=False):
        raise TypeError("Cannot intersect %r with %r." % (self, other))

    def _build_union_(self, members):
        return self.__class__(members)

    def to_ast(self):
        r = super().to_ast()
        r['type'] = 'union'
        return r


class Predicate(_ImmutableBase):
    def __init__(self, *args, **kwargs):
        self._truthy = any(map(bool, args)) or any(map(bool, kwargs.values()))

        self._freeze_()

    def __hash__(self):
        # This trivially satisfies the property:
        # x == x => hash(x) == hash(x)
        # Subclasses ought to override this with something less... collision-y.
        return 0

    def __eq__(self, other):
        raise NotImplementedError

    def __contains__(self, value):
        return self._is_element_(value)

    def _is_element_(self, value):
        return True

    def __bool__(self):
        return self._truthy

    def __le__(self, other):
        if other is None:
            other = self.__class__()
        return self._is_subtype_(other)

    def __ge__(self, other):
        if other is None:
            other = self.__class__()
        return other._is_subtype_(self)

    def _aug_is_subtype(self, other):
        r = self._is_subtype_(other)
        if r is NotImplemented:
            return other._is_supertype_(self)
        return r

    def _is_subtype_(self, other):
        raise NotImplementedError

    def _is_supertype_(self, other):
        raise NotImplementedError

    def to_ast(self):
        return {
            'type': 'predicate',
            'name': self.__class__.__name__
        }


# TODO: finish these classes:
class IntersectionTypeExpression(_SetOperationBase):
    _operator = '&'

    def _validate_union_(self, other, handshake=False):
        raise TypeError("Cannot union %r with %r." % (self, other))

    def _build_intersection_(self, members):
        return self.__class__(members)

    def to_ast(self):
        r = super().to_ast()
        r['type'] = 'intersection'
        return r


class MappingTypeExpression(TypeExpression):
    def __init__(self, name, mapping):
        if type(mapping) is not dict:  # we really only want dict literals
            raise ValueError()

        if type(name) is not str:
            raise ValueError()

        for key in mapping:
            self._validate_member_(key)
        for value in mapping.values():
            self._validate_member_(value)
        # Read only proxy of mapping, mutation to `mapping` will be reflected
        # but there isn't much we can do about that. Good use of this object
        # would involve a dict literal anyway.
        self.mapping = types.MappingProxyType(mapping)

        super().__init__(name)

    def __hash__(self):
        return super().__hash__() ^ hash(frozenset(self.mapping.items()))

    def __eq__(self, other):
        super_eq = super().__eq__(other)
        if super_eq is NotImplemented:
            return NotImplemented
        return super_eq and (set(self.mapping.items()) ==
                             set(other.mapping.items()))

    def _validate_predicate_(self, predicate):
        raise TypeError("Cannot apply predicates to type variables.")

    def _validate_intersection_(self, other, handshake=False):
        if type(self) != type(other):
            raise TypeError()
        if set(self.mapping) != set(other.mapping):
            raise TypeError()

        super()._validate_intersection_(other, handshake=handshake)

    def _validate_union_(self, other, handshake=False):
        # This has a reasonable definition (ensure disjoint sets on left-hand)
        # the opposite of intersection, but there isn't really a good use-case
        # for it at this time.
        raise TypeError("Cannot union type variables.")

    def to_ast(self):
        return {
            "type": "map",
            "mapping": [list(item) for item in self.mapping.items()]
        }
