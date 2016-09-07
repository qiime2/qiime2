# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import re
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


def parse_type(string, expect=None):
    """Convert a string into a TypeExpression

    Parameters
    ----------
    string : str
        The string type expression to convert into a TypeExpression
    expect : {'semantic', 'primitive', 'visualization'}, optional
        Will raise a TypeError if the resulting TypeExpression is not a member
        of `expect`.

    Returns
    -------
    TypeExpression

    Raises
    ------
    ValueError
        Raised when `expect` has an invalid value
    TypeError
        Raised when the expression contains invalid characters
    TypeError
        Raised when the expression does not result in a member of `expect`
    UnknownTypeError
        Raised when unkown types are present in the expression.

    """
    # Avoid circular imports
    import qiime.sdk
    import qiime.core.type as qtype

    if expect is not None and expect not in {'semantic', 'primitive',
                                             'visualization'}:
        raise ValueError("`expect` got %r, must be 'semantic', 'primitive',"
                         " 'visualization', or None." % (expect,))

    if '\n' in string or '\r' in string or ';' in string:
        raise TypeError("Found multiple statements in type expression %r. Will"
                        " not evaluate to avoid arbitrary code execution."
                        % string)

    pm = qiime.sdk.PluginManager()
    locals_ = {n: getattr(qtype, n) for n in qtype.__all__ if '_' not in n}
    locals_.update({k: v.semantic_type for k, v in pm.semantic_types.items()})

    try:
        type_expr = eval(string, {'__builtins__': {}}, locals_)
        if expect is None:
            pass
        elif expect == 'semantic' and qtype.is_semantic_type(type_expr):
            pass
        elif expect == 'primitive' and qtype.is_primitive_type(type_expr):
            pass
        elif expect == 'visualization' and type_expr == qtype.Visualization:
            pass
        else:
            raise TypeError("Type expression %r is not a %s type."
                            % (type_expr, expect))
        return type_expr
    except NameError as e:
        # http://stackoverflow.com/a/2270822/579416
        name, = re.findall("name '(\w+)' is not defined", str(e))
        raise UnknownTypeError("Name %r is not a defined QIIME type, a plugin"
                               " may be needed to define it." % name)


# Makes it possible to programmatically detect when a type doesn't exist.
class UnknownTypeError(TypeError):
    pass


def parse_format(format_str):
    # Avoid circular imports
    import qiime.sdk
    from qiime.plugin.model.base import FormatBase

    pm = qiime.sdk.PluginManager()
    for type_format_record in pm.type_formats:
        if type_format_record.format.__name__ == format_str:
            return type_format_record.format

    for input in pm.transformers:
        if issubclass(input, FormatBase) and input.__name__ == format_str:
            return input
        for output in pm.transformers[input]:
            if (issubclass(output, FormatBase) and
                    output.__name__ == format_str):
                return output

    raise TypeError("No format: %s" % format_str)


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
