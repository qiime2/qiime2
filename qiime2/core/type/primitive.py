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

from qiime2.core.type.template import (
    TypeTemplate, PredicateTemplate, instantiate)
import qiime2.metadata as metadata
import qiime2.core.util as util


def is_primitive_type(expr):
    return expr.kind == 'primitive'


_RANGE_DEFAULT_START = float('-inf')
_RANGE_DEFAULT_END = float('inf')
_RANGE_DEFAULT_INCLUSIVE_START = True
_RANGE_DEFAULT_INCLUSIVE_END = False


class _PrimitivePredicateBase(PredicateTemplate):
    def get_kind(self):
        return 'primitive'

    def get_name(self):
        return self.__class__.__name__

    def validate_union(self, other):
        pass

    def validate_intersection(self, other):
        pass


class Range(_PrimitivePredicateBase):
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

        return "Range(%s)" % (', '.join(args),)

    def is_element(self, value, expr):
        if self.inclusive_start:
            if value < self.start:
                return False
        elif value <= self.start:
            return False

        if self.inclusive_end:
            if value > self.end:
                return False
        elif value >= self.end:
            return False

        return True

    def is_subtype(self, other):
        if type(self) is not type(other):
            return False

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

    def is_supertype(self, other):
        if type(self) is not type(other):
            return False

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

    def collapse_intersection(self, other):
        if type(self) is not type(other):
            return None

        if self.start < other.start:
            new_start = other.start
            new_inclusive_start = other.inclusive_start
        elif other.start < self.start:
            new_start = self.start
            new_inclusive_start = self.inclusive_start
        else:
            new_start = self.start
            new_inclusive_start = (
                self.inclusive_start and other.inclusive_start)

        if self.end > other.end:
            new_end = other.end
            new_inclusive_end = other.inclusive_end
        elif other.end > self.end:
            new_end = self.end
            new_inclusive_end = self.inclusive_end
        else:
            new_end = self.end
            new_inclusive_end = self.inclusive_end and other.inclusive_end

        if new_end < new_start:
            return None
        if (new_start == new_end
                and not (new_inclusive_start and new_inclusive_end)):
            return None

        return self.__class__(new_start, new_end,
                              inclusive_start=new_inclusive_start,
                              inclusive_end=new_inclusive_end)


    def iter_boundaries(self):
        if self.start != float('-inf'):
            yield self.start
        if self.end != float('inf'):
            yield self.end


def Start(start, inclusive=_RANGE_DEFAULT_INCLUSIVE_START):
    return Range(start, _RANGE_DEFAULT_END, inclusive_start=inclusive)


def End(end, inclusive=_RANGE_DEFAULT_INCLUSIVE_END):
    return Range(_RANGE_DEFAULT_START, end, inclusive_end=inclusive)


class Choices(_PrimitivePredicateBase):
    def __init__(self, *choices):
        if not choices:
            raise ValueError("'Choices' cannot be instantiated with an empty"
                             " set.")

        # Backwards compatibility with old Choices({1, 2, 3}) syntax
        if len(choices) == 1:
            if not isinstance(choices[0], (bool, str)):
                choices = choices[0]

        self.choices = choices = tuple(choices)
        if len(choices) != len(set(choices)):
            raise ValueError("Duplicates found in choices: %r"
                             % util.find_duplicates(choices))


    def __hash__(self):
        return hash(type(self)) ^ hash(frozenset(self.choices))

    def __eq__(self, other):
        return (type(self) == type(other)
                and set(self.choices) == set(other.choices))

    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__,
                           repr(list(self.choices))[1:-1])

    def is_element(self, value, expr):
        return value in self.choices

    def is_subtype(self, other):
        if type(self) is not type(other):
            return False

        return set(self.choices) <= set(other.choices)

    def is_supertype(self, other):
        if type(self) is not type(other):
            return False

        return set(self.choices) >= set(other.choices)

    def collapse_intersection(self, other):
        if type(self) is not type(other):
            return None

        new_choices_set = set(self.choices) & set(other.choices)
        if not new_choices_set:
            return None

        # order by appearance:
        new_choices = []
        for c in itertools.chain(self.choices, other.choices):
            if c in new_choices_set:
                new_choices.append(c)
                new_choices_set.remove(c)

        return self.__class__(new_choices)

    def iter_boundaries(self):
        yield from self.choices


class _PrimitiveTemplateBase(TypeTemplate):
    def __eq__(self, other):
        return type(self) is type(other)

    def __hash__(self):
        return hash(type(self))

    def get_name(self):
        return self.__class__.__name__

    def get_kind(self):
        return 'primitive'

    def get_field_names(self):
        return []

    def validate_fields(self, fields, expr):
        raise TypeError

    def validate_predicate(self, predicate, expr):
        return

        if type(predicate) not in self._valid_predicates:
            raise TypeError

        for bound in predicate.iter_boundaries():
            if not self.is_element(bound, expr):
                raise TypeError

    def validate_union(self, other):
        pass

    def validate_intersection(self, other):
        pass


@instantiate
class Int(_PrimitiveTemplateBase):
    _valid_predicates = {Range}

    def is_element(self, value, expr):
        return (value is not True and value is not False
                and isinstance(value, numbers.Integral))

    def is_subtype(self, other):
        if other.get_name() == 'Float':
            return True
        return super().is_subtype(other)


@instantiate
class Str(_PrimitiveTemplateBase):
    _valid_predicates = {Choices}

    def is_element(self, value, expr):
        return isinstance(value, str)


@instantiate
class Float(_PrimitiveTemplateBase):
    _valid_predicates = {Range}

    def is_element(self, value, expr):
        # Works with numpy just fine.
        return isinstance(value, numbers.Real)


@instantiate
class Bool(_PrimitiveTemplateBase):
    _valid_predicates = {Choices}

    def is_element(self, value, expr):
        return value is True or value is False


@instantiate
class Metadata(_PrimitiveTemplateBase):
    _valid_predicates = set()

    def is_element(self, value, expr):
        return isinstance(value, metadata.Metadata)


@instantiate
class MetadataColumn(_PrimitiveTemplateBase):
    _valid_predicates = set()

    def is_element(self, value, expr):
        return value in expr.fields[0]

    def get_field_names(self):
        return ["type"]

    def validate_fields(self, fields, expr):
        try:
            field, = fields
        except ValueError:
            raise TypeError("%r given too many fields: %r" % (self, fields))

        if field.get_name() not in ["Numeric", "Categorical"]:
            raise TypeError("Unsupported type in field: %r" % (field,))


@instantiate
class Categorical(_PrimitiveTemplateBase):
    _valid_predicates = set()

    def is_element(self, value, expr):
        return isinstance(value, metadata.CategoricalMetadataColumn)


@instantiate
class Numeric(_PrimitiveTemplateBase):
    _valid_predicates = set()

    def is_element(self, value, expr):
        return isinstance(value, metadata.NumericMetadataColumn)
