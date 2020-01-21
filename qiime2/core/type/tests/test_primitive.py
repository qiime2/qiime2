# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
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
