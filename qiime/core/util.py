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


def parse_type(string):
    # Avoid circular imports
    import qiime.sdk
    from .type import Visualization, Properties

    if len(string.split('\n')) != 1:
        raise TypeError("Found multiple lines in type expression %r. Will not "
                        "evaluate to avoid arbitrary code execution."
                        % string)

    if ';' in string:
        raise TypeError("Invalid type expression %r. Will not evaluate to"
                        " avoid arbitrary code execution." % string)

    pm = qiime.sdk.PluginManager()
    locals_ = {k: v[1] for k, v in pm.semantic_types.items()}
    locals_[Visualization.name] = Visualization
    locals_["Properties"] = Properties

    try:
        return eval(string, {'__builtins__': {}}, locals_)
    except NameError as e:
        # http://stackoverflow.com/a/2270822/579416
        name, = re.findall("name '(\w+)' is not defined", str(e))
        raise UnknownTypeError("Name %r is not a defined QIIME type, a plugin"
                               " may be needed to define it." % name)


# Makes it possible to programmatically detect when a type doesn't exist.
class UnknownTypeError(TypeError):
    pass
