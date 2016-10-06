# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import contextlib
import warnings
import hashlib

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


def human_readable(relative_delta):
    attrs = ['years', 'months', 'days', 'hours', 'minutes', 'seconds', 'microseconds']
    results = []
    for attr in attrs:
        value = getattr(relative_delta, attr)
        if value != 0:
            if value == 1:
                attr = attr[:-1]
            results.append("%d %s" % (value, attr))
    text = results[-1]
    if results[:-1]:
        text = ', and '.join([', '.join(results[:-1]), text])
    return text


def md5sum(filepath):
    md5 = hashlib.md5()
    with open(str(filepath), mode='rb') as fh:
        for chunk in iter(lambda: fh.read(8192), b""):
            md5.update(chunk)
    return md5.hexdigest()


def md5sum_directory(directory):
    sums = {}
    for root, dirs, files in os.walk(str(directory), topdown=True):
        dirs[:] = [d for d in dirs if not d[0] == '.']
        for file in files:
            if file[0] == '.':
                continue

            path = os.path.join(root, file)
            sums[os.path.relpath(path, start=directory)] = md5sum(path)
    return sums


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
