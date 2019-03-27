# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

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


if __name__ == '__main__':
    unittest.main()
