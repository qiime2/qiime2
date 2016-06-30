# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import re


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
    locals_.update({k: v[1] for k, v in pm.semantic_types.items()})

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
