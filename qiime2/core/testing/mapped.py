# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os

from .plugin import dummy_plugin, C1, C2, C3, Foo, Bar, Baz, EchoFormat
from qiime2.plugin import (
    TypeMap, TypeMatch, Properties, Visualization, Bool, Choices)


def constrained_input_visualization(output_dir: str, a: EchoFormat,
                                    b: EchoFormat):
    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write("<p>%s</p>" % a.path.read_text())
        fh.write("<p>%s</p>" % b.path.read_text())


T, U, V = TypeMap({
    (Foo, Foo): Visualization,
    (Bar, Bar): Visualization,
    (Baz, Baz): Visualization,
    (C1[Foo], C1[Foo]): Visualization,
    (C1[Bar], C1[Bar]): Visualization,
    (C1[Baz], C1[Baz]): Visualization
})
dummy_plugin.visualizers.register_function(
    function=constrained_input_visualization,
    inputs={
        'a': T,
        'b': U
    },
    parameters={},
    name="Constrained Input Visualization",
    description="Ensure Foo/Bar/Baz match"
)
del T, U, V


def combinatorically_mapped_method(a: EchoFormat, b: EchoFormat
                                   ) -> (EchoFormat, EchoFormat):
    return a, b


T, R = TypeMap({
    Foo: Bar,
    Bar: Baz,
    Baz: Foo
})
X, Y = TypeMap({
    C3[Foo | Bar | Baz, Foo | Bar | Baz, Foo]: Foo,
    C3[Foo | Bar | Baz, Foo | Bar | Baz, Bar]: Bar,
    C3[Foo | Bar | Baz, Foo | Bar | Baz, Baz]: Baz
})
dummy_plugin.methods.register_function(
    function=combinatorically_mapped_method,
    inputs={
        'a': C1[T],
        'b': X
    },
    parameters={},
    outputs=[
        ('x', C2[R, R]),
        ('y', Y)
    ],
    name="Combinatorically Mapped Method",
    description="Test that multiple typemaps can be used"
)
del T, R, X, Y


def double_bound_variable_method(a: EchoFormat, b: EchoFormat,
                                 extra: EchoFormat) -> EchoFormat:
    return extra


T, R = TypeMap({
    Foo: Bar,
    Bar: Baz,
    Baz: Foo
})
dummy_plugin.methods.register_function(
    function=double_bound_variable_method,
    inputs={
        'a': T,
        'b': T,
        'extra': Foo
    },
    parameters={},
    outputs=[
        ('x', R)
    ],
    name="Double Bound Variable Method",
    description="Test reuse of variables"
)
del T, R


def bool_flag_swaps_output_method(a: EchoFormat, b: bool) -> EchoFormat:
    return a


P, R = TypeMap({
    Choices(True): C1[Foo],
    Choices(False): Foo
})
dummy_plugin.methods.register_function(
    function=bool_flag_swaps_output_method,
    inputs={
        'a': Bar
    },
    parameters={
        'b': Bool % P
    },
    outputs=[
        ('x', R)
    ],
    name='Bool Flag Swaps Output Method',
    description='Test if a parameter can change output'
)
del P, R


def predicates_preserved_method(a: EchoFormat) -> EchoFormat:
    return a


P = TypeMatch([Properties('A'), Properties('B'), Properties('C'),
               Properties('X', 'Y')])
dummy_plugin.methods.register_function(
    function=predicates_preserved_method,
    inputs={
        'a': Foo % P
    },
    parameters={},
    outputs=[
        ('x', Foo % P)
    ],
    name='Predicates Preserved Method',
    description='Test that predicates are preserved'
)
del P
