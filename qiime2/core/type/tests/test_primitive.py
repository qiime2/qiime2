# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import pandas as pd

import qiime2.metadata as metadata
import qiime2.core.type.primitive as primitive
import qiime2.core.type.grammar as grammar


class TestIntersectTwoRanges(unittest.TestCase):

    def assertIntersectEqual(self, a, b, exp):
        r1 = a & b
        r2 = b & a

        self.assertEqual(r1, r2)
        self.assertEqual(r1, exp)

    def test_overlap_simple(self):
        a = primitive.Range(0, 10)
        b = primitive.Range(3, 7)

        self.assertIntersectEqual(a, b, b)

    def test_overlap_inclusive_point(self):
        a = primitive.Range(0, 5, inclusive_end=True)
        b = primitive.Range(5, 10)
        exp = primitive.Range(5, 5, inclusive_start=True, inclusive_end=True)

        self.assertIntersectEqual(a, b, exp)

    def test_disjoint_far(self):
        a = primitive.Range(-10, -5)
        b = primitive.Range(5, 10)

        self.assertIntersectEqual(a, b, grammar.UnionExp())

    def test_disjoint_exclusive_point(self):
        a = primitive.Range(0, 5, inclusive_end=False)
        b = primitive.Range(5, 9, inclusive_start=False)

        self.assertIntersectEqual(a, b, grammar.UnionExp())


class TestChoices(unittest.TestCase):
    def test_list_constructor(self):
        choices = primitive.Choices(['a', 'b', 'c'])

        self.assertEqual(choices.template.choices, ('a', 'b', 'c'))
        self.assertIn('a', choices)
        self.assertNotIn('x', choices)

    def test_set_constructor(self):
        choices = primitive.Choices({'a', 'b', 'c'})

        self.assertEqual(choices.template.choices, ('a', 'b', 'c'))
        self.assertIn('a', choices)
        self.assertNotIn('x', choices)

    def test_varargs_constructor(self):
        choices = primitive.Choices('a', 'b', 'c')

        self.assertEqual(choices.template.choices, ('a', 'b', 'c'))
        self.assertIn('a', choices)
        self.assertNotIn('x', choices)

    def test_union(self):
        a = primitive.Choices('a', 'b', 'c')
        b = primitive.Choices('x', 'y', 'z')

        r = a | b

        self.assertIn('a', r)
        self.assertIn('x', r)
        self.assertNotIn('foo', r)

    def test_intersection(self):
        a = primitive.Choices('a', 'b', 'c')
        b = primitive.Choices('a', 'c', 'z')

        r = a & b

        self.assertIn('a', r)
        self.assertIn('c', r)
        self.assertNotIn('b', r)
        self.assertNotIn('z', r)


class TestMetadataColumn(unittest.TestCase):

    def test_decode_categorical_value(self):
        value = pd.Series({'a': 'a', 'b': 'b', 'c': 'c'}, name='foo')
        value.index.name = 'id'
        cat_md = metadata.CategoricalMetadataColumn(value)

        res = primitive.MetadataColumn[primitive.Categorical].decode(cat_md)
        self.assertIs(res, cat_md)

    def test_decode_numeric_value(self):
        value = pd.Series({'a': 1, 'b': 2, 'c': 3}, name='foo')
        value.index.name = 'id'
        num_md = metadata.NumericMetadataColumn(value)

        res = primitive.MetadataColumn[primitive.Categorical].decode(num_md)
        self.assertIs(res, num_md)

    def test_decode_other(self):
        with self.assertRaisesRegex(TypeError, 'provided.*directly'):
            primitive.MetadataColumn[primitive.Categorical].decode(
                "<metadata>")


if __name__ == '__main__':
    unittest.main()
