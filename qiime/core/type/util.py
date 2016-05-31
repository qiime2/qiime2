# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


# TODO make this more robust (hand me a list and see what happens...)
def tuplize(x):
    if type(x) is not tuple:
        return (x,)
    return x


# TODO: Move this, but for now it provides some reasonable documentation to
# the type system
def overrides(cls):
    def decorator(func):
        if not hasattr(cls, func.__name__):
            raise AssertionError("%r does not override %r"
                                 % (func, cls.__name__))
        return func
    return decorator
