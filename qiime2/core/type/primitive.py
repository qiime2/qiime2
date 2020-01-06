# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numbers
import itertools

from qiime2.core.type.template import TypeTemplate, PredicateTemplate
import qiime2.metadata as metadata
import qiime2.core.util as util


_RANGE_DEFAULT_START = float('-inf')
_RANGE_DEFAULT_END = float('inf')
_RANGE_DEFAULT_INCLUSIVE_START = True
_RANGE_DEFAULT_INCLUSIVE_END = False


class _PrimitivePredicateBase(PredicateTemplate):
    def get_kind(self):
        return 'primitive'

    def get_name(self):
        return self.__class__.__name__


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
            raise ValueError("Too many arguments passed, expected 0, 1, or 2.")
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
        start = self.start
        if start == float('-inf'):
            start = None
        end = self.end
        if end == float('inf'):
            end = None

        args.append(repr(start))
        args.append(repr(end))
        if self.inclusive_start is not _RANGE_DEFAULT_INCLUSIVE_START:
            args.append('inclusive_start=%r' % self.inclusive_start)
        if self.inclusive_end is not _RANGE_DEFAULT_INCLUSIVE_END:
            args.append('inclusive_end=%r' % self.inclusive_end)

        return "Range(%s)" % (', '.join(args),)

    def is_element(self, value):
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

    def is_symbol_subtype(self, other):
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

    def is_symbol_supertype(self, other):
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
                              inclusive_end=new_inclusive_end).template

    def iter_boundaries(self):
        if self.start != float('-inf'):
            yield self.start
        if self.end != float('inf'):
            yield self.end

    def update_ast(self, ast):
        start = self.start
        if start == float('-inf'):
            start = None

        end = self.end
        if end == float('inf'):
            end = None

        ast['range'] = [start, end]
        ast['inclusive'] = [self.inclusive_start, self.inclusive_end]


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

    def is_element(self, value):
        return value in self.choices

    def is_symbol_subtype(self, other):
        if type(self) is not type(other):
            return False

        return set(self.choices) <= set(other.choices)

    def is_symbol_supertype(self, other):
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

        return self.__class__(new_choices).template

    def iter_boundaries(self):
        yield from self.choices

    def update_ast(self, ast):
        ast['choices'] = list(self.choices)

    def unpack_union(self):
        for c in self.choices:
            yield self.__class__(c)


class _PrimitiveTemplateBase(TypeTemplate):
    public_proxy = 'encode', 'decode'

    def __eq__(self, other):
        return type(self) is type(other)

    def __hash__(self):
        return hash(type(self))

    def get_name(self):
        return self.__class__.__name__[1:]  # drop `_`

    def get_kind(self):
        return 'primitive'

    def get_field_names(self):
        return []

    def validate_field(self, name, field):
        raise NotImplementedError

    def validate_predicate_expr(self, self_expr, predicate_expr):
        predicate = predicate_expr.template
        if type(predicate) not in self._valid_predicates:
            raise TypeError(str(predicate_expr))

        for bound in predicate.iter_boundaries():
            if not self.is_element_expr(self_expr, bound):
                raise TypeError(bound)

    def validate_predicate(self, predicate):
        raise NotImplementedError


class _Int(_PrimitiveTemplateBase):
    _valid_predicates = {Range}

    def is_element(self, value):
        return (value is not True and value is not False
                and isinstance(value, numbers.Integral))

    def is_symbol_subtype(self, other):
        if other.get_name() == 'Float':
            return True
        return super().is_symbol_subtype(other)

    def decode(self, string):
        return int(string)

    def encode(self, value):
        return str(value)


class _Str(_PrimitiveTemplateBase):
    _valid_predicates = {Choices}

    def is_element(self, value):
        return isinstance(value, str)

    def decode(self, string):
        return str(string)

    def encode(self, value):
        return str(value)


class _Float(_PrimitiveTemplateBase):
    _valid_predicates = {Range}

    def is_symbol_supertype(self, other):
        if other.get_name() == 'Int':
            return True
        return super().is_symbol_supertype(other)

    def is_element(self, value):
        # Works with numpy just fine.
        return (value is not True and value is not False
                and isinstance(value, numbers.Real))

    def decode(self, string):
        return float(string)

    def encode(self, value):
        return str(value)


class _Bool(_PrimitiveTemplateBase):
    _valid_predicates = {Choices}

    def is_element(self, value):
        return value is True or value is False

    def validate_predicate(self, predicate):
        if type(predicate) is Choices:
            if set(predicate.iter_boundaries()) == {True, False}:
                raise TypeError("Choices should be ommitted when "
                                "Choices(True, False).")

    def decode(self, string):
        if string not in ('false', 'true'):
            raise TypeError("%s is neither 'true' or 'false'" % string)

        return string == 'true'

    def encode(self, value):
        if value:
            return 'true'
        else:
            return 'false'


class _Metadata(_PrimitiveTemplateBase):
    _valid_predicates = set()

    def is_element(self, value):
        return isinstance(value, metadata.Metadata)

    def decode(self, metadata):
        # This interface should have already retrieved this object.
        if not self.is_element(metadata):
            raise TypeError("`Metadata` must be provided by the interface"
                            " directly.")
        return metadata

    def encode(self, value):
        # TODO: Should this be the provenance representation? Does that affect
        # decode?

        return value


class _MetadataColumn(_PrimitiveTemplateBase):
    _valid_predicates = set()

    def is_element_expr(self, self_expr, value):
        return value in self_expr.fields[0]

    def is_element(self, value):
        raise NotImplementedError

    def get_field_names(self):
        return ["type"]

    def validate_field(self, name, field):
        if field.get_name() not in ("Numeric", "Categorical"):
            raise TypeError("Unsupported type in field: %r"
                            % (field.get_name(),))

    def decode(self, value):
        # This interface should have already retrieved this object.
        if not isinstance(value, metadata.MetadataColumn):
            raise TypeError("`Metadata` must be provided by the interface"
                            " directly.")
        return value

    def encode(self, value):
        # TODO: Should this be the provenance representation? Does that affect
        # decode?

        return value


class _Categorical(_PrimitiveTemplateBase):
    _valid_predicates = set()

    def get_union_membership_expr(self, self_expr):
        return 'metadata-column'

    def is_element(self, value):
        return isinstance(value, metadata.CategoricalMetadataColumn)


class _Numeric(_PrimitiveTemplateBase):
    _valid_predicates = set()

    def get_union_membership_expr(self, self_expr):
        return 'metadata-column'

    def is_element(self, value):
        return isinstance(value, metadata.NumericMetadataColumn)


Int = _Int()
Float = _Float()
Bool = _Bool()
Str = _Str()
Metadata = _Metadata()
MetadataColumn = _MetadataColumn()
Categorical = _Categorical()
Numeric = _Numeric()


def infer_primitive_type(value):
    for t in (Int, Float):
        if value in t:
            return t % Range(value, value, inclusive_end=True)
    for t in (Bool, Str):
        if value in t:
            return t % Choices(value)
    for t in (Metadata, MetadataColumn[Categorical], MetadataColumn[Numeric]):
        if value in t:
            return t
    raise ValueError("Unknown primitive type: %r" % (value,))
