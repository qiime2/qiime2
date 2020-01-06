# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from qiime2.core.type import (
    is_collection_type, is_primitive_type, is_semantic_type, Set, List,
    SemanticType, Int, Metadata, MetadataColumn, Categorical, Numeric, Range)


class TestIsTypes(unittest.TestCase):
    def test_list_semantic_type(self):
        Foo = SemanticType('Foo')

        self.assertTrue(is_collection_type(List[Foo]))
        self.assertTrue(is_semantic_type(List[Foo]))
        self.assertFalse(is_primitive_type(List[Foo]))

    def test_set_semantic_type(self):
        Foo = SemanticType('Foo')

        self.assertTrue(is_collection_type(Set[Foo]))
        self.assertTrue(is_semantic_type(Set[Foo]))
        self.assertFalse(is_primitive_type(Set[Foo]))

    def test_list_primitive_type(self):
        self.assertTrue(is_collection_type(List[Int % Range(5)]))
        self.assertTrue(is_primitive_type(List[Int % Range(5)]))
        self.assertFalse(is_semantic_type(List[Int % Range(5)]))

    def test_set_primitive_type(self):
        self.assertTrue(is_collection_type(Set[Int % Range(5)]))
        self.assertTrue(is_primitive_type(Set[Int % Range(5)]))
        self.assertFalse(is_semantic_type(Set[Int % Range(5)]))


class TestCollectionBase(unittest.TestCase):
    def test_no_list_metadata(self):
        with self.assertRaisesRegex(TypeError, 'metadata'):
            List[Metadata]

    def test_no_set_metadata(self):
        with self.assertRaisesRegex(TypeError, 'metadata'):
            List[Metadata]

    def test_no_list_metadata_column(self):
        with self.assertRaisesRegex(TypeError, 'metadata'):
            List[MetadataColumn[Categorical]]

        with self.assertRaisesRegex(TypeError, 'metadata'):
            List[MetadataColumn[Numeric]]

    def test_no_set_metadata_column(self):
        with self.assertRaisesRegex(TypeError, 'metadata'):
            Set[MetadataColumn[Categorical]]

        with self.assertRaisesRegex(TypeError, 'metadata'):
            Set[MetadataColumn[Numeric]]

    def test_no_nesting_list_list(self):
        with self.assertRaisesRegex(TypeError, 'nest'):
            List[List[Int]]

    def test_no_nesting_set_set(self):
        with self.assertRaisesRegex(TypeError, 'nest'):
            Set[Set[Int]]

    def test_no_nesting_mixed(self):
        with self.assertRaisesRegex(TypeError, 'nest'):
            List[Set[Int]]


class TestCollectionExpression(unittest.TestCase):
    def test_bad_union(self):
        with self.assertRaisesRegex(TypeError, 'not union'):
            List[Int] | Set[Int]

    def test_union_inside_collection(self):
        Foo = SemanticType('Foo')
        Bar = SemanticType('Bar')

        self.assertTrue(List[Foo] <= List[Foo | Bar])

    def test_no_predicate(self):
        with self.assertRaisesRegex(TypeError, 'cannot be applied'):
            List[Int] % Range(5)

    def is_concrete(self):
        Foo = SemanticType('Foo')

        self.assertFalse(List[Foo].is_concrete())
        self.assertFalse(Set[Int].is_concrete())

    def test_to_ast_semantic(self):
        Foo = SemanticType('Foo')

        ast = List[Foo].to_ast()
        self.assertEqual(ast['fields'][0], Foo.to_ast())

    def test_to_ast_primitive(self):
        ast = List[Int % Range(5)].to_ast()
        self.assertEqual(ast['fields'][0], (Int % Range(5)).to_ast())

    def test_contains_list_primitive(self):
        self.assertTrue([1, 2, 3] in List[Int])
        self.assertTrue([-1, 2, 3] in List[Int])

        self.assertFalse([-1, 2, 3] in List[Int % Range(0, 5)])
        self.assertFalse([1, 1.1, 1.11] in List[Int])
        self.assertFalse({1, 2, 3} in List[Int])
        self.assertFalse(object() in List[Int])

    def test_contains_set_primitive(self):
        self.assertTrue({1, 2, 3} in Set[Int])
        self.assertTrue({-1, 2, 3} in Set[Int])

        self.assertFalse({-1, 2, 3} in Set[Int % Range(0, 5)])
        self.assertFalse({1, 1.1, 1.11} in Set[Int])
        self.assertFalse([1, 2, 3] in Set[Int])
        self.assertFalse(object() in Set[Int])

    def test_variant_of_field_members(self):
        Bar = SemanticType('Bar')
        Foo = SemanticType('Foo', field_names='foo',
                           field_members={'foo': List[Bar]})
        with self.assertRaisesRegex(TypeError, 'is not a variant'):
            Foo[List[Bar]]

    def test_variant_of_alt(self):
        Foo = SemanticType('Foo', field_names='foo')
        Bar = SemanticType('Bar', variant_of=Foo.field['foo'])

        with self.assertRaisesRegex(TypeError, 'is not a variant'):
            Foo[Set[Bar]]

    def test_encode_decode_set(self):
        value = List[Int].decode("[1, 2, 3]")
        self.assertEqual(value, [1, 2, 3])

        json = List[Int].encode(value)
        self.assertEqual(json, "[1, 2, 3]")

    def test_encode_decode_list(self):
        value = Set[Int].decode("[1, 2, 3]")
        self.assertEqual(value, {1, 2, 3})

        json = Set[Int].encode(value)
        self.assertEqual(json, "[1, 2, 3]")


if __name__ == '__main__':
    unittest.main()
