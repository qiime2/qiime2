# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
from qiime.sdk.type import Type, Any, Nil

class IdentityTests:

    def test_equality(self):
        # X == Y iff Y < X < Y
        self.assertTrue(self.Foo < self.Foo)
        self.assertTrue(self.Foo > self.Foo)
        self.assertTrue(self.Foo == self.Foo)
        self.assertFalse(self.Foo != self.Foo)

        self.assertTrue(self.Bar < self.Bar)
        self.assertTrue(self.Bar > self.Bar)
        self.assertTrue(self.Bar == self.Bar)
        self.assertFalse(self.Bar != self.Bar)

    def test_subtype_identities(self):
        # Consider Foo[Bar], we can represent different types as sets
        # representing their respective domains:

        #  TT FF BB FF TT
        # (  (  (  )  )  )
        # Where T is any other member of F's sum-type, F is Foo, and B is Bar.

        # There are 8 ways to fill this in:
        # 0: (  (  (  )  )  )  Nil
        # 1: (  (  (%%)  )  )  Foo[ Bar]
        # 2: (  (%%(  )%%)  )  Foo[~Bar]
        # 3: (  (%%(%%)%%)  )  Foo
        # 4: (%%(  (  )  )%%) ~Foo
        # 5: (%%(  (%%)  )%%) ~Foo[~Bar]  this is an arguably useless selection
        # 6: (%%(%%(  )%%)%%) ~Foo[ Bar]
        # 7: (%%(%%(%%)%%)%%)  Any

        # From that, we can derive the following non-trivial identities:

        # 1 < 3; 1 < 5; 1 == 1
        self.assertTrue(self.Foo[self.Bar] < self.Foo)
        self.assertTrue(self.Foo[self.Bar] < ~self.Foo[~self.Bar])
        self.assertEqual(self.Foo[self.Bar], self.Foo[self.Bar])

        # 2 < 3; 2 < 6; 2 == 2
        self.assertTrue(self.Foo[~self.Bar] < self.Foo)
        self.assertTrue(self.Foo[~self.Bar] < ~self.Foo[self.Bar])
        self.assertEqual(self.Foo[~self.Bar], self.Foo[~self.Bar])

        # 3 > 2; 3 > 1; 3 == 3
        self.assertTrue(self.Foo > self.Foo[~self.Bar])
        self.assertTrue(self.Foo > self.Foo[self.Bar])
        self.assertEqual(self.Foo, self.Foo)

        # 4 < 5; 4 < 6; 4 == 4
        self.assertTrue(~self.Foo < ~self.Foo[~self.Bar])
        self.assertTrue(~self.Foo < ~self.Foo[self.Bar])
        self.assertEqual(~self.Foo, ~self.Foo)

        # 5 > 4; 5 > 1; 5 == 5
        self.assertTrue(~self.Foo[~self.Bar] > ~self.Foo)
        self.assertTrue(~self.Foo[~self.Bar] > self.Foo[self.Bar])
        self.assertEqual(~self.Foo[~self.Bar], ~self.Foo[~self.Bar])

        # 6 > 4; 6 > 2; 6 == 6
        self.assertTrue(~self.Foo[self.Bar] > ~self.Foo)
        self.assertTrue(~self.Foo[self.Bar] > self.Foo[~self.Bar])
        self.assertEqual(~self.Foo[self.Bar], ~self.Foo[self.Bar])

        # Because Foo is itself a product-type, the situation of Foo[Bar, Baz]
        # would be a cartesian product of the above cases. This would be 64
        # different ways to fill in the sets, which is more than anyone has
        # time for. Instead, we can infer that Baz is independent from Bar,
        # therefore as long every field is considered before declaring the
        # subtype-edness of Foo, the answer will still be correct. Therefore
        # the above tests are sufficient.

class TestAlgebraicDataTypes(unittest.TestCase, IdentityTests):
    def setUp(self):

        class Foo(Type, fields='X'):
            class X: pass

        class Bar(Type, variant_of=Foo.X):
            pass

        self.Foo = Foo
        self.Bar = Bar

class TestRefinementTypes:
    pass

class TestNegatedRefinementTypes:
    pass

if __name__ == '__main__':
    unittest.main()
