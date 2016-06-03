# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from qiime.core.type.grammar import TypeExpression, CompositeType, Predicate
import json


def is_primitive_type(type_):
    return isinstance(type_, _PrimitiveBase)


def _symbolify(*args, **kwargs):
    """Used to bootstrap the primitive type classes.
    They don't have a factory like the semantic types do.

    """
    def decorator(cls):
        return cls(*args, **kwargs)
    return decorator


class _PrimitiveBase(TypeExpression):
    def __contains__(self, value):
        # TODO: Make this work correctly. It is stubbed as True as current
        # usage is through encode/decode, we will need this to work for
        # the Artifact API and predicates.
        return True

    def _validate_union_(self, other, handshake=False):
        # It is possible we may want this someday: `Int | Str`, but order of
        # encode/decode dispatch wouldn't be straight-forward.
        # Order of the declaration could indicate "MRO", but then equality
        # must consider Int | Str != Str | Int, which feels weird.
        raise TypeError("Cannot union primitive types.")

    def _validate_intersection_(self, other, handshake=False):
        # This literally makes no sense for primitives.
        raise TypeError("Cannot intersect primitive types.")


class _Primitive(_PrimitiveBase):
    def __init__(self, **kwargs):
        # Use class name to avoid duplicating the name twice
        super().__init__(self.__class__.__name__, **kwargs)

    def _validate_predicate_(self, predicate):
        super()._validate_predicate_(predicate)
        # _valid_predicates will be on the class obj of the primitives to
        # make defining them easy
        if not isinstance(predicate, tuple(self._valid_predicates)):
            raise TypeError("Cannot supply %r as a predicate to %r type,"
                            " permitted predicates: %r"
                            % (predicate, self, self._valid_predicates))

        # for bound in predicate.iter_boundaries():
        #     self._validate_predicate_bound_(bound)

    def _apply_predicate_(self, predicate):
        return self.__class__(fields=self.fields, predicate=predicate)


class _Collection(CompositeType):
    def __init__(self, field_names):
        super().__init__(self.__class__.__name__, field_names)

    def _validate_field_(self, name, value):
        if not isinstance(value, _Primitive):
            if isinstance(value, _CollectionPrimitive):
                raise TypeError("Cannot nest collection types.")
            else:
                raise TypeError("Collection type (%r) must be provided"
                                " primitives as arguments to its fields,"
                                " not %r" % (self, value))

        super()._validate_field_(name, value)

    def _build_expression_(self, fields):
        return _CollectionPrimitive(self.encode, self.decode, self.name,
                                    fields=fields)


class _CollectionPrimitive(_PrimitiveBase):
    def __init__(self, encode, decode, *args, **kwargs):
        # TODO: This is a nasty hack
        self._encode = encode
        self._decode = decode

    def encode(self, value):
        return self._encode(value)

    def decode(self, string):
        return self._decode(string)

    def _validate_predicate_(self, predicate):
        raise TypeError("Predicates cannot be applied directly to collection"
                        " types.")


_RANGE_DEFAULT_START = 0
_RANGE_DEFAULT_END = None
_RANGE_DEFAULT_INCLUSIVE_START = True
_RANGE_DEFAULT_INCLUSIVE_END = False


class Range(Predicate):
    def __init__(self, *args, inclusive_start=_RANGE_DEFAULT_INCLUSIVE_START,
                 inclusive_end=_RANGE_DEFAULT_INCLUSIVE_END):
        # TODO: Make this not silly.
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

    def __contains__(self, value):
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

    def __repr__(self):
        args = []
        if self.start is not _RANGE_DEFAULT_START:
            args.append(repr(self.start))
        if self.end is not _RANGE_DEFAULT_END:
            args.append(repr(self.end))
        if self.inclusive_start is not _RANGE_DEFAULT_INCLUSIVE_START:
            args.append('inclusive_start=%r' % self.inclusive_start)
        if self.inclusive_end is not _RANGE_DEFAULT_INCLUSIVE_END:
            args.append('inclusive_end=%r' % self.inclusive_end)

        return "%s(%s)" % (self.__class__.__name__, ', '.join(args))


class Choices(Predicate):
    def __init__(self, choices):
        self.choices = set(choices)

        super().__init__(choices)

    def __repr__(self):
        return "%s({%s})" % (self.__class__.__name__,
                             repr(sorted(self.choices))[1:-1])


class Arguments(Predicate):
    def __init__(self, parameter):
        self.parameter = parameter

        super().__init__(parameter)

    def __repr__(self):
        return "%s(%r)" % (self.__class__.__name__, self.parameter)


@_symbolify(field_names=['keys', 'values'])
class Dict(_Collection):
    def decode(self, string):
        return json.loads(string)

    def encode(self, value):
        return json.dumps(value)


@_symbolify(field_names=['elements'])
class List(_Collection):
    def decode(self, string):
        return json.loads(string)

    def encode(self, value):
        return json.dumps(value)


@_symbolify(field_names=['elements'])
class Set(_Collection):
    def decode(self, string):
        return set(json.loads(string))

    def encode(self, value):
        return json.dumps(list(value))


@_symbolify()
class Int(_Primitive):
    _valid_predicates = {Range, Arguments}

    def decode(self, string):
        return int(string)

    def encode(self, value):
        return str(value)


@_symbolify()
class Str(_Primitive):
    _valid_predicates = {Choices, Arguments}
    decode = encode = lambda self, x: x


@_symbolify()
class Float(_Primitive):
    _valid_predicates = {Range, Arguments}

    def decode(self, string):
        return float(string)

    def encode(self, value):
        return str(value)


@_symbolify()
class Color(type(Str)):
    pass

# @_symbolify
# class Column(_Primitive):
#     pass
