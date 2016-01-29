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
