# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import inspect

class _AutoInstanceMeta(type):
    def __new__(metaclass, name, bases, dict):
        """Instantiate the object, not the classobject!

        Typically `__new__` of a metaclass would return a class instance, but
        we go one step further and instantiate the class instance.

        __init__ is called in both cases to make sure everything is square.

        """
        cls = type.__new__(metaclass, name, bases, dict)
        type.__init__(cls, name, bases, dict)
        # cls is now initialized
        obj = cls.__new__(cls, name, bases, dict)
        cls.__init__(obj)
        # obj is now initialized: `obj` will be bound to the name of the class
        # when it is defined instead of `cls` as is typical
        return obj

    def __call__(self, name, bases, dict):
        """Override the () notation on the class.

        This effectively allows us to use an instance as a type. In other
        words the class can act as its own metaclass as its call() is defered
        to the metaclass __new__ (which creates a new instance).

        (In the expression: `class X(Y): pass` the metaclass of Y is used to
         construct X, i.e. type(Y) is used. Therefore, if Y was an
         instance, then type(Y) is its class. So by overriding __call__ of the
         metaclass, we can change the behavior of type(Y)(name, bases, dict),
         which is the signature provided to metaclasses by class definitions)

        """
        # We have to get the types of our bases as they are likely to be
        # instances.
        bases = tuple([b if inspect.isclass(b) else type(b) for b in bases])
        return _AutoInstanceMeta.__new__(type(self), name, bases, dict)

class EvilSymbol(metaclass=_AutoInstanceMeta):
    """Doom any decendants to bind their class name to an instance of class."""
    def __init__(self):
        """Needed for metaclass to work.

        `object` isn't a very real object, so it whines about not requiring
        arguments. This overrides object.__init__ so that self is expected
        instead of implied."""
        if hasattr(self, '_initialized'):
            raise Exception("Already initialized, if subclassing, use: "
                            "`type(%s).__init__(self)` instead."
                            % self.__class__.__name__)
        self._initialized = True

    def _constructor_(self, *args, **kwargs):
        """Since we don't have the call syntax anymore, recreate it."""
        cls = self.__class__
        obj = cls.__new__(cls)
        obj.__init__(*args, **kwargs)
        return obj

class BaseType(EvilSymbol):
    def __init__(self, subtypes=None, restrictions=None):
        # `super` doesn't work, use oldschool style, but remember EvilSymbol is
        # an instance, not a type.
        type(EvilSymbol).__init__(self)

        self._subtypes = subtypes
        self._restrictions = restrictions

    @property
    def subtypes(self):
        return self._subtypes

    @property
    def restrictions(self):
        return self._restrictions

    def __getitem__(self, subtypes):
        if self.subtypes is not None:
            raise Exception()
        return self._copy(subtypes=subtypes)

    def __mod__(self, restrictions):
        if self.restrictions is not None:
            raise Exception()
        return self._copy(restrictions=restrictions)

    def _copy(self, subtypes=None, restrictions=None):
        if subtypes is None:
            subtypes = self.subtypes
        if restrictions is None:
            restrictions = self.restrictions
        return self._constructor_(subtypes=subtypes, restrictions=restrictions)

class ArtifactType(BaseType):
    pass

class PrimitiveType(BaseType):
    pass

class Artifact:
    @classmethod
    def load(cls, qtf):
        reuturn cls()

    def __init__(self, value=None, directory_cls=Directory, **kwargs):
        self.__value = None
        self.provenance = Provenance()
        self.type = Type()
        self.uuid = UUID()
        self._directory = directory_cls()
        self._value = value

    @property
    def _value(self):
        return self.__value

    @_value.setter
    def _value(self, value):
        if value is not None:
            self.type.validate(value)
        self.__value = value

    @property
    def value(self):
        if self._value is None:
            self._value = self.type.deserialize(self._directory)
        return self._value

    def save(qtf):
        pass
