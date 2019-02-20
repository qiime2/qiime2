# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numbers
import functools
import itertools

from qiime2.core.type.grammar import (
    TypeExpression, CompositeType, Predicate, UnionTypeExpression,
    IntersectionPredicate, UnionPredicate)
import qiime2.metadata as metadata
import qiime2.core.util as util


def is_primitive_type(type_):
    return isinstance(type_, _PrimitiveBase)


class _PrimitiveBase(TypeExpression):
    def _validate_union_(self, other):
        pass

    def _validate_intersection_(self, other):
        pass


class _Primitive(_PrimitiveBase):

    def _build_union_(self, *members):
        return UnionPrimitiveType(*members)

    def _validate_predicate_(self, predicate):
        super()._validate_predicate_(predicate)
        # _valid_predicates will be on the class obj of the primitives to
        # make defining them easy
        valid = tuple(self._valid_predicates) + (
            UnionPrimitivePredicate, IntersectionPrimitivePredicate)
        if not isinstance(predicate, valid):
            raise TypeError("Cannot supply %r as a predicate to %r type,"
                            " permitted predicates: %r"
                            % (predicate, self,
                               {p.__name__ for p in self._valid_predicates}))

        for bound in predicate.iter_boundaries():
            if not self._is_element_(bound):
                raise TypeError("%r has the wrong types for %r."
                                % (predicate, self))

    def to_ast(self):
        ast = super().to_ast()
        ast['type'] = 'primitive'
        return ast


class UnionPrimitiveType(_Primitive, UnionTypeExpression):
    def _get_order(self):
        order = []
        for t in (Bool, Int, Float, Str):
            for m in self.members:
                if m.name == t.name:
                    order.append(m)

        return order

    def encode(self, value):
        for m in self._get_order():
            if value in m:
                return m.encode(value)

    def decode(self, string):
        for m in self._get_order():
            try:
                value = m.decode(string)
            except Exception:
                pass
            else:
                break
        return value


class _PrimitivePredicate(Predicate):
    def _build_union_(self, *members):
        return UnionPrimitivePredicate(*members)

    def _build_intersection_(self, *members):
        return IntersectionPrimitivePredicate(*members)


class UnionPrimitivePredicate(_PrimitivePredicate, UnionPredicate):
    def iter_boundaries(self):
        for m in self.members:
            yield from m.iter_boundaries()


class IntersectionPrimitivePredicate(_PrimitivePredicate,
                                     IntersectionPredicate):
    def iter_boundaries(self):
        for m in self.members:
            yield from m.iter_boundaries()


_RANGE_DEFAULT_START = float('-inf')
_RANGE_DEFAULT_END = float('inf')
_RANGE_DEFAULT_INCLUSIVE_START = True
_RANGE_DEFAULT_INCLUSIVE_END = False


class Range(_PrimitivePredicate):
    def __init__(self, *args, inclusive_start=_RANGE_DEFAULT_INCLUSIVE_START,
                 inclusive_end=_RANGE_DEFAULT_INCLUSIVE_END):
        if len(args) == 2:
            self.start, self.end = args
        elif len(args) == 1:
            self.start = _RANGE_DEFAULT_START
            self.end, = args
        elif len(args) == 0:
            self.start = _RANGE_DEFAULT_START
            self.end = _RANGE_DEFAULT_END
        else:
            raise ValueError("")
        self.inclusive_start = inclusive_start
        self.inclusive_end = inclusive_end

        if self.start is None:
            self.start = _RANGE_DEFAULT_START
        if self.end is None:
            self.end = _RANGE_DEFAULT_END

        if self.end < self.start:
            raise ValueError("End of range precedes start.")

        falsey = (self.start == _RANGE_DEFAULT_START
                  and self.end == _RANGE_DEFAULT_END)

        super().__init__(not falsey)

    def __hash__(self):
        return (hash(type(self)) ^
                hash(self.start) ^
                hash(self.end) ^
                hash(self.inclusive_start) ^
                hash(self.inclusive_end))

    def __eq__(self, other):
        return (type(self) is type(other) and
                self.start == other.start and
                self.end == other.end and
                self.inclusive_start == other.inclusive_start and
                self.inclusive_end == other.inclusive_end)

    def __repr__(self):
        args = []
        args.append(repr(self.start))
        args.append(repr(self.end))
        if self.inclusive_start is not _RANGE_DEFAULT_INCLUSIVE_START:
            args.append('inclusive_start=%r' % self.inclusive_start)
        if self.inclusive_end is not _RANGE_DEFAULT_INCLUSIVE_END:
            args.append('inclusive_end=%r' % self.inclusive_end)

        return "%s(%s)" % (self.__class__.__name__, ', '.join(args))

    def _is_element_(self, value):
        if self.start is not None:
            if self.inclusive_start:
                if value < self.start:
                    return False
            elif value <= self.start:
                return False

        if self.end is not None:
            if self.inclusive_end:
                if value > self.end:
                    return False
            elif value >= self.end:
                return False

        return True

    def _is_subtype_(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented

        if other.start > self.start:
            return False
        elif (other.start == self.start
                and (not other.inclusive_start)
                and self.inclusive_start):
            return False

        if other.end < self.end:
            return False
        elif (other.end == self.end
                and (not other.inclusive_end)
                and self.inclusive_end):
            return False

        return True

    def _is_supertype_(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented

        if other.start < self.start:
            return False
        elif (other.start == self.start
                and (not self.inclusive_start)
                and other.inclusive_start):
            return False

        if other.end > self.end:
            return False
        elif (other.end == self.end
                and (not self.inclusive_end)
                and other.inclusive_end):
            return False

        return True

    def _build_intersection_(self, *members):
        if not members:
            return super()._build_intersection_()

        for m in members:
            if not isinstance(m, self.__class__):
                return self.get_bottom()

        r = functools.reduce(self.intersect_two_ranges, members)
        if r is None:
            return self.get_bottom()
        return r

    @classmethod
    def intersect_two_ranges(cls, a, b):
        if a is None or b is None:
            return None

        if a.start < b.start:
            new_start = b.start
            new_inclusive_start = b.inclusive_start
        elif b.start < a.start:
            new_start = a.start
            new_inclusive_start = a.inclusive_start
        else:
            new_start = a.start
            new_inclusive_start = a.inclusive_start and b.inclusive_start

        if a.end > b.end:
            new_end = b.end
            new_inclusive_end = b.inclusive_end
        elif b.end > a.end:
            new_end = a.end
            new_inclusive_end = a.inclusive_end
        else:
            new_end = a.end
            new_inclusive_end = a.inclusive_end and b.inclusive_end

        if new_end < new_start:
            return None
        if (new_start == new_end
                and not (new_inclusive_start and new_inclusive_end)):
            return None

        return cls(new_start, new_end, inclusive_start=new_inclusive_start,
                   inclusive_end=new_inclusive_end)

    def iter_boundaries(self):
        if self.start != float('-inf'):
            yield self.start
        if self.end != float('inf'):
            yield self.end

    def to_ast(self):
        ast = super().to_ast()
        ast['start'] = self.start
        ast['end'] = self.end
        ast['inclusive-start'] = self.inclusive_start
        ast['inclusive-end'] = self.inclusive_end
        return ast


def Start(start, inclusive=_RANGE_DEFAULT_INCLUSIVE_START):
    return Range(start, _RANGE_DEFAULT_END, inclusive_start=inclusive)


def End(end, inclusive=_RANGE_DEFAULT_INCLUSIVE_END):
    return Range(_RANGE_DEFAULT_START, end, inclusive_end=inclusive)


class Choices(_PrimitivePredicate):
    def __init__(self, *choices):
        if not choices:
            raise ValueError("'Choices' cannot be instantiated with an empty"
                             " set.")

        # Backwards compatibility with old Choices({1, 2, 3}) syntax
        if len(choices) == 1:
            if isinstance(choices[0], (list, set, tuple, frozenset)):
                choices = choices[0]

        self.choices = tuple(choices)
        if len(choices) != len(set(choices)):
            raise ValueError("Duplicates found in choices: %r"
                             % util.find_duplicates(choices))

        super().__init__(choices)

    def __hash__(self):
        return hash(type(self)) ^ hash(frozenset(self.choices))

    def __eq__(self, other):
        return (type(self) == type(other)
                and set(self.choices) == set(other.choices))

    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__,
                           repr(list(self.choices))[1:-1])

    def _is_element_(self, value):
        return value in self.choices

    def _is_subtype_(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented

        return set(self.choices) <= set(other.choices)

    def _is_supertype_(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented

        return set(self.choices) >= set(other.choices)

    def _build_intersection_(self, *members):
        if not members:
            return super()._build_intersection_()

        for m in members:
            if not isinstance(m, self.__class__):
                return self.get_bottom()

        members_iter = iter(members)

        new_choices_set = set(next(members_iter).choices)
        for m in members_iter:
            new_choices_set &= set(m.choices)

        if not new_choices_set:
            return self.get_bottom()

        # order by appearance:
        new_choices = []
        for c in itertools.chain.from_iterable(m.choices for m in members):
            if c in new_choices_set:
                new_choices.append(c)
                new_choices_set.remove(c)

        return self.__class__(new_choices)

    def iter_boundaries(self):
        yield from self.choices

    def to_ast(self):
        ast = super().to_ast()
        ast['choices'] = list(self.choices)
        return ast


class _Int(_Primitive):
    _valid_predicates = {Range}

    def _is_element_(self, value):
        # Works with numpy just fine.
        return (isinstance(value, numbers.Integral)
                and value is not True
                and value is not False)

    def _is_subtype_(self, other):
        if (isinstance(other, type(Float)) and
                self.evaled_predicate <= other.evaled_predicate):
            return True
        return super()._is_subtype_(other)

    def decode(self, string):
        return int(string)

    def encode(self, value):
        return str(value)


class _Str(_Primitive):
    _valid_predicates = {Choices}
    decode = encode = lambda self, arg: arg

    def _is_element_(self, value):
        # No reason for excluding bytes other than extreme prejudice.
        return isinstance(value, str)


class _Float(_Primitive):
    _valid_predicates = {Range}

    def _is_element_(self, value):
        # Works with numpy just fine.
        return isinstance(value, numbers.Real)

    def decode(self, string):
        return float(string)

    def encode(self, value):
        return str(value)


class _Bool(_Primitive):
    _valid_predicates = {Choices}

    def _is_element_(self, value):
        return isinstance(value, bool)

    def decode(self, string):
        if string not in ('false', 'true'):
            raise TypeError("%s is neither 'true' or 'false'" % string)

        return string == 'true'

    def encode(self, value):
        if value:
            return 'true'
        else:
            return 'false'


class _Metadata(_Primitive):
    _valid_predicates = set()

    def _is_element_(self, value):
        return isinstance(value, metadata.Metadata)

    def decode(self, metadata):
        # This interface should have already retrieved this object.
        if not self._is_element_(metadata):
            raise TypeError("`Metadata` must be provided by the interface"
                            " directly.")
        return metadata

    def encode(self, value):
        # TODO: Should this be the provenance representation? Does that affect
        # decode?
        return value


class _MetadataColumn(CompositeType):
    def _validate_field_(self, name, value):
        if not isinstance(value, (_MetadataColumnType,
                                  _MetadataColumnTypeUnion)):
            raise TypeError("Unsupported type in field: %r" % value)

    def _apply_fields_(self, fields):
        return _MetadataColumnExpression(self.name, fields=fields)


class _MetadataColumnExpression(_Primitive):
    _valid_predicates = set()

    def _is_element_(self, value):
        return value in self.fields[0]

    def decode(self, metadata_column):
        # This interface should have already retrieved this object.
        if metadata_column not in self:
            raise TypeError("`MetadataColumn` must be provided by the "
                            "interface directly.")
        return metadata_column

    def encode(self, value):
        # TODO: Should this be the provenance representation? Does that affect
        # decode?
        return value


class _MetadataColumnType(_Primitive):
    _valid_predicates = set()

    def __init__(self, name, view, fields=(), predicate=None):
        self._view = view
        super().__init__(name, fields, predicate)

    def _is_element_(self, value):
        return isinstance(value, self._view)

    def _apply_fields_(self, fields):
        self.__class__(self.name, self._view, fields=fields,
                       predicate=self.predicate)

    def _validate_union_(self, other, handshake=False):
        if not isinstance(other, self.__class__):
            raise TypeError("Unsupported union: %r" % other)

    def _build_union_(self, *members):
        return _MetadataColumnTypeUnion(*members)


class _MetadataColumnTypeUnion(UnionTypeExpression):
    pass


Bool = _Bool('Bool')
Int = _Int('Int')
Float = _Float('Float')
Str = _Str('Str')
Metadata = _Metadata('Metadata')
MetadataColumn = _MetadataColumn('MetadataColumn', field_names=['type'])
Numeric = _MetadataColumnType('Numeric', metadata.NumericMetadataColumn)
Categorical = _MetadataColumnType('Categorical',
                                  metadata.CategoricalMetadataColumn)
