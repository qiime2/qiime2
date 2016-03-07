# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from collections import namedtuple
import operator
import inspect

from qiime.util import tuplize


class VariantType:
    def __init__(self, interface):
        self.interface = interface

    def is_member(self, other):
        if other is Any or other is Nil:
            return True
        return self in other().__variant__()

    def validate(self, cls):
        """Make sure a class meets the variant's interface"""
        # TODO: implement for realz
        return True


class FieldSignature(namedtuple('FieldSigTuple', ['variant_types'])):
    # def __getitem__(self, i):
    #     return self.variant_types[i]

    def __len__(self):
        return len(self.variant_types)

    def __iter__(self):
        return iter(self.variant_types)

    def validate(self, fields):
        if len(self.variant_types) != len(fields):
            raise TypeError()
        for field, variant in zip(fields, self.variant_types):
            if not variant.is_member(field):
                raise TypeError()

    def get_default_fields(self):
        return tuple([Any] * len(self))


class TypeMeta(type, object):
    def __new__(meta, *args, **kwargs):
        return super().__new__(meta, *args)

    def __init__(cls, name, bases, dct,
                 variant_of=(), fields=(),
                 generic=False, inheritance='invariant'):
        super().__init__(name, bases, dct)
        cls.__instance = None
        # TODO: check that __lt__ and __gt__ have been overriden if True
        cls._is_generic = generic

        # Please, let us never need to indicate contra/co/bi-variant.
        if inheritance != 'invariant':
            raise TypeError('No. Inheritance will remain invariant.')

        cls._set_variants(variant_of)
        cls._set_fields(fields)

    def _set_variants(cls, variant_types):
        variant_types = tuplize(variant_types)
        for variant in variant_types:
            variant.validate(cls)
        cls.variants = variant_types

    def _set_fields(cls, fields):
        wrapped = []
        for field in tuplize(fields):
            wrapped_field = VariantType(getattr(cls, field))
            setattr(cls, field, wrapped_field)
            wrapped.append(wrapped_field)
        cls.field_signature = FieldSignature(tuple(wrapped))

    @property
    def _instance(cls):
        if cls.__instance is None:
            try:
                cls.__instance = cls()
            except TypeError:
                # We are an ephemeral type like Not and don't exist without
                # an instance of another type. None of the class-level
                # operations will occur anyways.
                pass
        return cls.__instance

    def __getattr__(cls, attr):
        return getattr(cls._instance, attr)

    def __getitem__(cls, fields):
        return cls._instance.__getitem__(fields)

    def __mod__(cls, predicates):
        return cls._instance.__mod__(predicates)

    def __rmod__(cls, predicates):
        return cls._instance.__rmod__(predicates)

    def __repr__(cls):
        return cls._instance.__repr__()

    def __iter__(cls):
        return cls._instance.__iter__()

    def __invert__(cls):
        return cls._instance.__invert__()

    def __lt__(cls, other):
        return cls._instance.__lt__(other)

    def __gt__(cls, other):
        return cls._instance.__gt__(other)


class Type(metaclass=TypeMeta, fields=('Artifact', 'Primitive')):
    def __init__(self, fields=None, predicates=()):
        super().__init__()
        if fields is None:
            self.fields = self.field_signature.get_default_fields()
        else:
            self.field_signature.validate(fields)
            self.fields = fields

        self.predicates = predicates

    def __getitem__(self, fields):
        if self.fields != self.field_signature.get_default_fields():
            raise TypeError()
        return self(fields=tuplize(fields))

    def __mod__(self, predicates):
        if self.predicates:
            raise TypeError()
        return self(predicates=tuplize(predicates))

    def __rmod__(self, predicates):
        raise TypeError("No. Do it right.")

    def __iter__(self):
        yield self
        for f in self.fields:
            yield from f

    def __repr__(self):
        r = c = ""

        if self.predicates:
            predicates = (self.predicates,)
            if len(self.predicates) == 1:
                predicates = self.predicates
            r = " %% %r" % predicates

        if self.fields:
            right = -1
            if len(self.fields) == 1:
                right = -2
            c = "[%s]" % repr(self.fields)[1:right]

        return self.__class__.__name__ + c + r

    def __call__(self, fields=None, predicates=None):
        if fields is None:
            fields = self.fields
        if predicates is None:
            predicates = self.predicates
        return self.__class__(fields=fields, predicates=predicates)

    def __lt__(self, other):
        return self._cmp(other, operator.lt)

    def __gt__(self, other):
        return self._cmp(other, operator.gt)

    def _cmp(self, other, op):
        if other._is_generic:
            # The other type is more qualified to answer this question
            return NotImplemented
        elif type(self) is type(other()):  # noqa
            # Note: This assumes type-inheritance is invariant.
            return all([op(f1, f2) for f1, f2 in
                       zip(self.fields, other.fields)])
        return False

    def __eq__(self, other):
        return other > self > other

    def __ne__(self, other):
        # TODO: Should this be `(not self < other) and (not self > other)`?
        return not (self == other)

    def __invert__(self):
        return _NotType(self)

    def __variant__(self):
        return self.__class__.variants

    def is_concrete(self):
        for type_ in self:
            if type_._is_generic:
                return False
        return True

    def get_imports(self):
        imps = []
        for type_ in self:
            cls = type(type_())
            imp = inspect.getmodule(cls)
            imps.append((cls.__name__, imp.__name__))
        return set(imps)

    class Artifact:
        def save(self, data, data_writer):
            pass

        def load(data_reader):
            pass

    class Primitive:
        def from_string(str):
            pass

        def to_string(model):
            pass


class _SimpleNotUnion(Type, generic=True):
    def __init__(self, type_):
        self._instance = type_

    def __call__(self):
        return self

    def __invert__(self):
        return self._instance

    def __lt__(self, other):
        return type(other()) is self.__class__ \
            and other._instance == self._instance

    def __gt__(self, other):
        return type(other()) is self.__class__ \
            and other._instance == self._instance

    def __mod__(self, other):
        raise TypeError()

    def __getitem__(self, other):
        raise TypeError()

    def __repr__(self):
        return repr(self._instance)

    def __variant__(self):
        return self._instance.__variant__()


class _NotType(Type, generic=True):
    def __init__(self, type_):
        super().__init__()
        self._instance = type_()
        self._negate = type_(fields=tuple(~f for f in type_.fields))
        if not self._negate.fields:
            self._negate = _SimpleNotUnion(self._instance)

    def __call__(self):
        return self

    def __invert__(self):
        return self._instance

    def __lt__(self, other):
        # T > ~X[Y] <==> T in X'\X & T > X[~Y]
        if other is Any or other is Nil:
            return NotImplemented
        if type(other()) is self.__class__:
            return self._instance > other._instance
        return self._is_relative_complement(other) and self._negate < other

    def __gt__(self, other):
        # T < ~X[Y] <==> T in X'\X | T < X[~Y]
        if other is Any or other is Nil:
            return NotImplemented
        if type(other()) is self.__class__:
            return self._instance < other._instance
        return self._is_relative_complement(other) or self._negate > other

    def _is_relative_complement(self, t):
        return sibling_variants(self, t) and not t < type(self._instance)

    def __mod__(self, other):
        raise TypeError()

    def __getitem__(self, other):
        raise TypeError()

    def __repr__(self):
        if self._instance.predicates:
            return "~(%r)" % self._instance
        return "~%r" % self._instance

    def __variant__(self):
        return self._instance.__variant__()


class Any:
    def __new__(cls):
        if not hasattr(cls, '_instance'):
            cls._instance = super().__new__(cls)
        cls._is_generic = True
        return cls._instance

    def __call__(self):
        return self

    def __lt__(self, other):
        return other is self.__class__._instance

    def __gt__(self, other):
        return True

    def __repr__(self):
        return "Any"

    def __invert__(self):
        return Nil


class Nil:
    def __new__(cls):
        if not hasattr(cls, '_instance'):
            cls._instance = super().__new__(cls)
        cls._is_generic = True
        return cls._instance

    def __call__(self):
        return self

    def __lt__(self, other):
        return True

    def __gt__(self, other):
        return other is self.__class__._instance

    def __repr__(self):
        return "Nil"

    def __invert__(self):
        return Any


Any = Any()
Nil = Nil()


def sibling_variants(a, b):
    if a is Any or b is Any:
        return True
    elif a is Nil or b is Nil:
        return False
    a_var = a().__variant__()
    b_var = b().__variant__()
    # They share a variant-type
    if a_var == b_var or bool(set(a_var) & set(b_var)):
        # Do more here?
        return True
    return False
