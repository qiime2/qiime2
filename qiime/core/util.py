# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import contextlib
import warnings

import decorator


def tuplize(x):
    if type(x) is not tuple:
        return (x,)
    return x


def overrides(cls):
    def decorator(func):
        if not hasattr(cls, func.__name__):
            raise AssertionError("%r does not override %r"
                                 % (func, cls.__name__))
        return func
    return decorator


@contextlib.contextmanager
def warning():
    def _warnformat(msg, category, filename, lineno, file=None, line=None):
        return '%s:%s: %s: %s\n' % (filename, lineno, category.__name__, msg)

    default_warn_format = warnings.formatwarning
    try:
        warnings.formatwarning = _warnformat
        yield warnings.warn
    finally:
        warnings.formatwarning = default_warn_format


# Descriptor protocol for creating an attribute that is bound to an
# (arbitrarily nested) attribute accessible to the instance at runtime.
class LateBindingAttribute:
    def __init__(self, attribute):
        self._attribute = attribute

    def __get__(self, obj, cls=None):
        attrs = self._attribute.split('.')
        curr_attr = obj
        for attr in attrs:
            curr_attr = getattr(curr_attr, attr)
        return staticmethod(curr_attr).__get__(obj, cls)


# Removes the first parameter from a callable's signature.
class DropFirstParameter(decorator.FunctionMaker):
    @classmethod
    def from_function(cls, function):
        return cls.create(function, "return None", {})

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.signature = self._remove_first_arg(self.signature)
        self.shortsignature = self._remove_first_arg(self.shortsignature)

    def _remove_first_arg(self, string):
        return ",".join(string.split(',')[1:])[1:]
