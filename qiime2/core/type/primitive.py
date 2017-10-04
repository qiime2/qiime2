# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json
import numbers
import re
import collections.abc

from qiime2.core.type.grammar import TypeExpression, CompositeType, Predicate
import qiime2.metadata as metadata


def is_primitive_type(type_):
    return isinstance(type_, _PrimitiveBase)


class _PrimitiveBase(TypeExpression):
    def _validate_union_(self, other, handshake=False):
        # It is possible we may want this someday: `Int | Str`, but order of
        # encode/decode dispatch wouldn't be straight-forward.
        # Order of the declaration could indicate "MRO", but then equality
        # must consider Int | Str != Str | Int, which feels weird.
        raise TypeError("Cannot union primitive types.")

    def _validate_intersection_(self, other, handshake=False):
        # This literally makes no sense for primitives. Even less sense than
        # C's Union type (which is actually an intersection type...)
        raise TypeError("Cannot intersect primitive types.")


class _Primitive(_PrimitiveBase):
    def _validate_predicate_(self, predicate):
        super()._validate_predicate_(predicate)
        # _valid_predicates will be on the class obj of the primitives to
        # make defining them easy
        if not isinstance(predicate, tuple(self._valid_predicates)):
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


class _Collection(CompositeType):
    def _validate_field_(self, name, value):
        if not isinstance(value, _Primitive):
            if isinstance(value, _CollectionPrimitive):
                raise TypeError("Cannot nest collection types.")
            else:
                raise TypeError("Collection type (%r) must be provided"
                                " primitives as arguments to its fields,"
                                " not %r" % (self, value))

        super()._validate_field_(name, value)

    def _apply_fields_(self, fields):
        return _CollectionPrimitive(self._is_element_.__func__,
                                    self.encode.__func__, self.decode.__func__,
                                    self.name, fields=fields)


class _CollectionPrimitive(_PrimitiveBase):
    def __init__(self, is_element, encode, decode, *args, **kwargs):
        # TODO: This is a nasty hack
        self._encode = encode
        self._decode = decode
        self._is_element = is_element

        super().__init__(*args, **kwargs)

    def encode(self, value):
        return self._encode(self, value)

    def decode(self, string):
        return self._decode(self, string)

    def _is_element_(self, value):
        return self._is_element(self, value)

    def _validate_predicate_(self, predicate):
        raise TypeError("Predicates cannot be applied directly to collection"
                        " types.")

    def _apply_fields_(self, fields):
        return _CollectionPrimitive(self._is_element, self._encode,
                                    self._decode, self.name, fields=fields)

    def to_ast(self):
        ast = super().to_ast()
        ast['type'] = 'collection'
        return ast


_RANGE_DEFAULT_START = None
_RANGE_DEFAULT_END = None
_RANGE_DEFAULT_INCLUSIVE_START = True
_RANGE_DEFAULT_INCLUSIVE_END = False


class Range(Predicate):
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

        super().__init__(args)

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

    def iter_boundaries(self):
        if self.start is not None:
            yield self.start
        if self.end is not None:
            yield self.end

    def to_ast(self):
        ast = super().to_ast()
        ast['start'] = self.start
        ast['end'] = self.end
        ast['inclusive-start'] = self.inclusive_start
        ast['inclusive-end'] = self.inclusive_end
        return ast


class Choices(Predicate):
    def __init__(self, choices):
        self.choices = set(choices)

        super().__init__(choices)

    def __hash__(self):
        return hash(type(self)) ^ hash(frozenset(self.choices))

    def __eq__(self, other):
        return type(self) == type(other) and self.choices == other.choices

    def __repr__(self):
        return "%s({%s})" % (self.__class__.__name__,
                             repr(sorted(self.choices))[1:-1])

    def _is_element_(self, value):
        return value in self.choices

    def iter_boundaries(self):
        yield from self.choices

    def to_ast(self):
        ast = super().to_ast()
        ast['choices'] = list(self.choices)
        return ast


class Arguments(Predicate):
    def __init__(self, parameter):
        self.parameter = parameter

        super().__init__(parameter)

    def __hash__(self):
        return hash(type(self)) ^ hash(self.parameter)

    def __eq__(self, other):
        return type(self) == type(other) and self.parameter == other.parameter

    def __repr__(self):
        return "%s(%r)" % (self.__class__.__name__, self.parameter)

    def _is_element_(self, value):
        raise NotImplementedError("Membership cannot be determined by this"
                                  " predicate directly.")

    def iter_boundaries(self):
        yield from []

    def to_ast(self):
        ast = super().to_ast()
        ast['parameter'] = self.parameter
        return ast


class _Dict(_Collection):
    def _is_element_(self, value):
        if not isinstance(value, collections.abc.Mapping):
            return False
        key_type, value_type = self.fields
        for k, v in value.items():
            if k not in key_type or v not in value_type:
                return False
        return True

    def decode(self, string):
        return json.loads(string)

    def encode(self, value):
        return json.dumps(value)


Dict = _Dict('Dict', field_names=['keys', 'values'])


class _List(_Collection):
    def _is_element_(self, value):
        if not isinstance(value, collections.abc.Sequence):
            return False
        element_type, = self.fields
        for v in value:
            if v not in element_type:
                return False
        return True

    def decode(self, string):
        return json.loads(string)

    def encode(self, value):
        return json.dumps(value)


List = _List('List', field_names=['elements'])


class _Set(_Collection):
    def _is_element_(self, value):
        if not isinstance(value, collections.abc.Set):
            return False
        element_type, = self.fields
        for v in value:
            if v not in element_type:
                return False
        return True

    def decode(self, string):
        return set(json.loads(string))

    def encode(self, value):
        return json.dumps(list(value))


Set = _Set('Set', field_names=['elements'])


class _Int(_Primitive):
    _valid_predicates = {Range, Arguments}

    def _is_element_(self, value):
        # Works with numpy just fine.
        return isinstance(value, numbers.Integral)

    def _is_subtype_(self, other):
        if (isinstance(other, type(Float)) and
                self.predicate <= other.predicate):
            return True
        return super()._is_subtype_(other)

    def decode(self, string):
        return int(string)

    def encode(self, value):
        return str(value)


Int = _Int('Int')


class _Str(_Primitive):
    _valid_predicates = {Choices, Arguments}
    decode = encode = lambda self, arg: arg

    def _is_element_(self, value):
        # No reason for excluding bytes other than extreme prejudice.
        return isinstance(value, str)


Str = _Str('Str')


class _Float(_Primitive):
    _valid_predicates = {Range, Arguments}

    def _is_element_(self, value):
        # Works with numpy just fine.
        return isinstance(value, numbers.Real)

    def decode(self, string):
        return float(string)

    def encode(self, value):
        return str(value)


Float = _Float('Float')


class _Bool(_Primitive):
    _valid_predicates = {Arguments, }

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


Bool = _Bool('Bool')


class _Color(type(Str)):
    def _is_element_(self, value):
        # Regex from: http://stackoverflow.com/a/1636354/579416
        return bool(re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', value))


Color = Colour = _Color('Color')


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


Metadata = _Metadata('Metadata')


class _MetadataCategory(_Primitive):
    _valid_predicates = set()

    def _is_element_(self, value):
        return isinstance(value, metadata.MetadataCategory)

    def decode(self, metadata_category):
        # This interface should have already retrieved this object.
        if not self._is_element_(metadata_category):
            raise TypeError("`MetadataCategory` must be provided by the"
                            " interface directly.")
        return metadata_category

    def encode(self, value):
        # TODO: Should this be the provenance representation? Does that affect
        # decode?
        return value


MetadataCategory = _MetadataCategory('MetadataCategory')
