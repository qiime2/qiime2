# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from qiime2.core.type.parse import ast_to_type, string_to_ast
from qiime2.core.testing.type import Foo, Bar, C1, C2
from qiime2.plugin import (Int, Float, Str, Bool, Range, Choices, TypeMap,
                           Properties, List, Set, Visualization, Metadata,
                           MetadataColumn, Categorical, Numeric)


class TestParsing(unittest.TestCase):
    def assert_roundtrip(self, type):
        ast = string_to_ast(repr(type))
        type1 = ast_to_type(ast)
        type2 = ast_to_type(type1.to_ast())

        self.assertEqual(type, type1)
        self.assertEqual(ast, type1.to_ast())
        self.assertEqual(type1, type2)

    def test_simple_semantic_type(self):
        self.assert_roundtrip(Foo)
        self.assert_roundtrip(Bar)
        self.assert_roundtrip(C1[Foo])

    def test_union_semantic_type(self):
        self.assert_roundtrip(Foo | Bar)
        self.assert_roundtrip(C1[Foo | Bar])

    def test_complicated_semantic_type(self):
        self.assert_roundtrip(C2[C1[Foo % Properties(["A", "B"]) | Bar],
                                 Foo % Properties("A")
                                 ] % Properties(exclude=["B", "C"]))

    def test_collection_semantic_type(self):
        self.assert_roundtrip(List[Foo | Bar])
        self.assert_roundtrip(Set[Bar])

    def test_visualization(self):
        self.assert_roundtrip(Visualization)

    def test_primitive_simple(self):
        self.assert_roundtrip(Int)
        self.assert_roundtrip(Float)
        self.assert_roundtrip(Str)
        self.assert_roundtrip(Bool)

    def test_primitive_predicate(self):
        self.assert_roundtrip(Int % Range(0, 10))
        self.assert_roundtrip(
            Int % (Range(0, 10) | Range(50, 100, inclusive_end=True)))
        self.assert_roundtrip(Float % Range(None, 10))
        self.assert_roundtrip(Float % Range(0, None))
        self.assert_roundtrip(Str % Choices("A"))
        self.assert_roundtrip(Str % Choices(["A"]))
        self.assert_roundtrip(Str % Choices("A", "B"))
        self.assert_roundtrip(Str % Choices(["A", "B"]))
        self.assert_roundtrip(Bool % Choices(True))
        self.assert_roundtrip(Bool % Choices(False))

    def test_collection_primitive(self):
        self.assert_roundtrip(Set[Str % Choices('A', 'B', 'C')])
        self.assert_roundtrip(List[Int % Range(1, 3, inclusive_end=True)
                                   | Str % Choices('A', 'B', 'C')])

    def test_metadata_primitive(self):
        self.assert_roundtrip(Metadata)
        self.assert_roundtrip(MetadataColumn[Numeric])
        self.assert_roundtrip(MetadataColumn[Categorical])
        self.assert_roundtrip(MetadataColumn[Numeric | Categorical])

    def test_typevars(self):
        T, U, V, W, X = TypeMap({
            (Foo, Bar, Str % Choices('A', 'B')): (C1[Foo], C1[Bar]),
            (Foo | Bar, Foo, Str): (C1[Bar], C1[Foo])
        })

        scope = {}
        T1 = ast_to_type(T.to_ast(), scope=scope)
        U1 = ast_to_type(U.to_ast(), scope=scope)
        V1 = ast_to_type(V.to_ast(), scope=scope)
        W1 = ast_to_type(W.to_ast(), scope=scope)
        X1 = ast_to_type(X.to_ast(), scope=scope)

        self.assertEqual(len(scope), 1)
        self.assertEqual(scope[id(T.mapping)], [T1, U1, V1, W1, X1])

        self.assertEqual(T1.mapping.lifted, T.mapping.lifted)

        self.assertIs(T1.mapping, U1.mapping)
        self.assertIs(U1.mapping, V1.mapping)
        self.assertIs(V1.mapping, W1.mapping)
        self.assertIs(W1.mapping, X1.mapping)

    def test_syntax_error(self):
        with self.assertRaisesRegex(ValueError, "could not be parsed"):
            string_to_ast('$')

    def test_bad_juju(self):
        with self.assertRaisesRegex(ValueError, "one type expression"):
            string_to_ast('import os; os.rmdir("something-important")')

    def test_more_bad(self):
        with self.assertRaisesRegex(ValueError, "Unknown expression"):
            string_to_ast('lambda x: x')

    def test_weird(self):
        with self.assertRaisesRegex(ValueError, "Unknown literal"):
            string_to_ast('FeatureTable(Foo + Bar)')


if __name__ == '__main__':
    unittest.main()
