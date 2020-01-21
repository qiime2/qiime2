# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import contextlib
import warnings
import hashlib
import os
import io
import collections

import decorator


def get_view_name(view):
    from .format import FormatBase
    if not isinstance(view, type):
        view = view.__class__

    if issubclass(view, FormatBase):
        # Not qualname because we don't have a notion of "nested" formats
        return view.__name__

    return ':'.join([view.__module__, view.__qualname__])


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


def superscript(number):
    table = {
        '0': chr(8304), '1': chr(185), '2': chr(178), '3': chr(179),
        **{str(i): chr(x) for i, x in enumerate(range(8308, 8314), 4)},
        'a': chr(7491), 'e': chr(7497), 'f': chr(7584), 'i': chr(8305),
        'n': chr(8319), '-': chr(8315), '.': chr(39), ',': chr(39)
    }
    return ''.join([table[d] for d in str(number)])


def find_duplicates(iterable):
    """Find duplicate values in an iterable.

    Parameters
    ----------
    iterable : iterable
        Iterable to search for duplicates.

    Returns
    -------
    set
        Values that are duplicated in `iterable`.

    Notes
    -----
    Values in `iterable` must be hashable.

    """
    # Modified from https://stackoverflow.com/a/9835819/3776794 to return
    # duplicates instead of remove duplicates from an iterable.
    seen = set()
    duplicates = set()
    for value in iterable:
        if value in seen:
            duplicates.add(value)
        else:
            seen.add(value)
    return duplicates


# Concept from: http://stackoverflow.com/a/11157649/579416
def duration_time(relative_delta):
    attrs = ['years', 'months', 'days', 'hours', 'minutes', 'seconds',
             'microseconds']
    results = []
    for attr in attrs:
        value = getattr(relative_delta, attr)
        if value != 0:
            if value == 1:
                # Remove plural 's'
                attr = attr[:-1]
            results.append("%d %s" % (value, attr))
    if results:
        text = results[-1]
        if results[:-1]:
            text = ', and '.join([', '.join(results[:-1]), text])
        return text
    else:
        # Great Scott! No time has passed!
        return '0 %s' % attrs[-1]


def md5sum(filepath):
    md5 = hashlib.md5()
    with open(str(filepath), mode='rb') as fh:
        for chunk in iter(lambda: fh.read(io.DEFAULT_BUFFER_SIZE), b""):
            md5.update(chunk)
    return md5.hexdigest()


def md5sum_directory(directory):
    directory = str(directory)
    sums = collections.OrderedDict()
    for root, dirs, files in os.walk(directory, topdown=True):
        dirs[:] = sorted([d for d in dirs if not d[0] == '.'])
        for file in sorted(files):
            if file[0] == '.':
                continue

            path = os.path.join(root, file)
            sums[os.path.relpath(path, start=directory)] = md5sum(path)
    return sums


def to_checksum_format(filepath, checksum):
    # see https://www.gnu.org
    # /software/coreutils/manual/html_node/md5sum-invocation.html
    if '\\' in filepath or '\n' in filepath:
        filepath = filepath.replace('\\', '\\\\').replace('\n', '\\n')
        checksum = '\\' + checksum

    return '%s  %s' % (checksum, filepath)


def from_checksum_format(line):
    line = line.rstrip('\n')
    parts = line.split('  ', 1)
    if len(parts) < 2:
        parts = line.split(' *', 1)

    checksum, filepath = parts

    if checksum[0] == '\\':
        chars = ''
        escape = False
        # Gross, but regular `.replace` will overlap with itself and
        # negative lookbehind in regex is *probably* harder than scanning
        for char in filepath:
            # 1) Escape next character
            if not escape and char == '\\':
                escape = True
                continue

            # 2) Handle escape sequence
            if escape:
                try:
                    chars += {'\\': '\\', 'n': '\n'}[char]
                except KeyError:
                    chars += '\\' + char  # Wasn't an escape after all
                escape = False
                continue

            # 3) Nothing interesting
            chars += char

        checksum = checksum[1:]
        filepath = chars

    return filepath, checksum


@contextlib.contextmanager
def warning():
    def _warnformat(msg, category, filename, lineno, file=None, line=None):
        return '%s:%s: %s: %s\n' % (filename, lineno, category.__name__, msg)

    default_warn_format = warnings.formatwarning
    try:
        warnings.formatwarning = _warnformat
        warnings.filterwarnings('always')
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


def _immutable_error(obj, *args):
    raise TypeError('%s is immutable.' % obj.__class__.__name__)


class ImmutableBase:
    def _freeze_(self):
        """Disables __setattr__ when called. It is idempotent."""
        self._frozen = True  # The particular value doesn't matter

    __delattr__ = __setitem__ = __delitem__ = _immutable_error

    def __setattr__(self, *args):
        # This doesn't stop silly things like
        # object.__setattr__(obj, ...), but that's a pretty rude thing
        # to do anyways. We are just trying to avoid accidental mutation.
        if hasattr(self, '_frozen'):
            _immutable_error(self)
        super().__setattr__(*args)
