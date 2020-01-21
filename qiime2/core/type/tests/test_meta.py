# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import pickle

import qiime2.core.type.collection as col
import qiime2.core.type.meta as meta
from qiime2.core.type.tests.test_grammar import MockTemplate, MockPredicate


class TestSelect(unittest.TestCase):
    def test_select_simple(self):
        Foo = MockTemplate('Foo')
        Bar = MockTemplate('Bar')

        X, Y = meta.TypeMap({
            Foo: Bar
        })

        sel, = meta.select_variables(X)

        self.assertIs(sel(X), X)
        self.assertIs(sel(Foo), Foo)
        self.assertIs(sel(X, swap=Foo), Foo)

    def test_select_inside_field(self):
        C2 = MockTemplate('C2', fields=('a', 'b'))
        Foo = MockTemplate('Foo')
        Bar = MockTemplate('Bar')

        X, Y = meta.TypeMap({
            Foo: Bar
        })

        sel, = meta.select_variables(C2[X, Foo])

        self.assertIs(sel(C2[X, Bar]), X)
        self.assertIs(sel(C2[Bar, Foo]), Bar)
        self.assertEqual(sel(C2[X, Foo], swap=Foo), C2[Foo, Foo])

    def test_select_predicate(self):
        Foo = MockTemplate('Foo')
        Bar = MockTemplate('Bar')
        P = MockPredicate('P')
        Q = MockPredicate('Q')

        X, Y = meta.TypeMap({
            P & Q: Foo,
            P: Bar,
            Q: Foo
        })

        sel, = meta.select_variables(Foo % X)

        self.assertIs(sel(Foo % X), X)
        self.assertIs(sel(Foo % P), P)
        self.assertEqual(sel(Foo % X, swap=Q), Foo % Q)

    def test_multiple_select(self):
        C2 = MockTemplate('C2', fields=('a', 'b'))
        Foo = MockTemplate('Foo')
        Bar = MockTemplate('Bar')
        P = MockPredicate('P')
        Q = MockPredicate('Q')
        R = MockPredicate('R')

        X1, Y1 = meta.TypeMap({
            Foo: Bar
        })

        X2, Y2 = meta.TypeMap({
            P & Q: Foo,
            P: Bar,
            Q: Foo
        })

        expr = C2[X1, Foo % X2] % X2
        pred_sel, field_sel, field_pred_sel = meta.select_variables(expr)

        self.assertIs(pred_sel(expr), X2)
        self.assertIs(pred_sel(C2[Bar, Foo % Q] % P), P)
        self.assertEqual(pred_sel(expr, swap=R), C2[X1, Foo % X2] % R)

        self.assertIs(field_sel(expr), X1)
        self.assertIs(field_sel(C2[Bar, Foo]), Bar)
        self.assertEqual(field_sel(expr, swap=Foo), C2[Foo, Foo % X2] % X2)

        self.assertIs(field_pred_sel(expr), X2)
        self.assertIs(field_pred_sel(C2[Bar, Foo % Q] % P), Q)
        self.assertEqual(field_pred_sel(expr, swap=R), C2[X1, Foo % R] % X2)


class TestTypeMap(unittest.TestCase):
    def test_missing_branch_requested(self):
        P = MockPredicate('P')
        Q = MockPredicate('Q')

        with self.assertRaisesRegex(ValueError, 'Ambiguous'):
            meta.TypeMap({P: P, Q: Q})

    def test_mismatched_pieces(self):
        P = MockPredicate('P')
        Bar = MockTemplate('Bar')

        with self.assertRaisesRegex(ValueError, 'in the same'):
            meta.TypeMap({P: P, Bar: Bar})

    def test_iter_sorted(self):
        P = MockPredicate('P', alphabetize=True)
        Q = MockPredicate('Q', alphabetize=True)
        Other = MockPredicate('Other')

        X, Y = meta.TypeMap({
            P & Other: Other, P: P, Q & Other: Other, Q: Q, Other: Other
        })
        mapping = X.mapping

        self.assertEqual(
            list(mapping.lifted),
            [col.Tuple[P & Other], col.Tuple[P], col.Tuple[Q & Other],
             col.Tuple[Q], col.Tuple[Other]])

    def test_variables(self):
        P = MockPredicate('P')
        Q = MockPredicate('Q')
        R = MockPredicate('R')
        S = MockPredicate('S')

        X, Y = meta.TypeMap({
            P & Q: R & S,
            P: R,
            Q: S,
        })

        self.assertEqual(X.members, (P & Q, P, Q))
        self.assertEqual(Y.members, (R & S, R, S))
        self.assertEqual(X.index, 0)
        self.assertEqual(Y.index, 1)
        self.assertTrue(X.input)
        self.assertTrue(Y.output)
        self.assertFalse(X.output)
        self.assertFalse(Y.input)

        # subtyping
        self.assertFalse(S <= X)
        self.assertFalse(R <= X)
        self.assertLessEqual(P, X)
        self.assertLessEqual(Q, X)
        self.assertLessEqual(P & Q, X)
        self.assertLessEqual(P & S | P & R, X)

    def test_pickle(self):
        Foo = MockTemplate('Foo')
        Bar = MockTemplate('Bar')

        X, Y = meta.TypeMap({
            Foo: Bar
        })

        X1, Y1 = pickle.loads(pickle.dumps((X, Y)))  # Pickled together

        self.assertIs(X1.mapping, Y1.mapping)

        self.assertEqual(X1.index, X.index)
        self.assertEqual(Y1.index, Y.index)

    def test_subtype(self):
        P = MockPredicate('P')
        Q = MockPredicate('Q')
        R = MockPredicate('R')
        S = MockPredicate('S')

        T, U, Y = meta.TypeMap({
            (P & Q, P & Q): R & S,
            (P & Q, Q): R & S,
            (P, P): R,
            (Q, Q): S
        })

        self.assertLessEqual(P, T)
        self.assertLessEqual(Q, T)

        self.assertFalse(P | Q <= T)

        self.assertLessEqual(T, U)


class TestTypeMatch(unittest.TestCase):
    def test_missing_branch_provided(self):
        P = MockPredicate('P')
        Q = MockPredicate('Q')

        T = meta.TypeMatch([P, Q])

        self.assertEqual(T.members, (P & Q, P, Q))

    def test_variable(self):
        P = MockPredicate('P')
        Q = MockPredicate('Q')
        R = MockPredicate('R')

        # This strange list is for branch coverage mostly
        T = meta.TypeMatch([P & Q, P, Q, R])

        self.assertTrue(T.input)
        self.assertTrue(T.output)
        self.assertEqual(T.index, 1)  # it really is supposed to be 1


class TestMatch(unittest.TestCase):
    def test_single_variable(self):
        Foo = MockTemplate('Foo')
        Bar = MockTemplate('Bar')
        P = MockPredicate('P')
        X, Y = meta.TypeMap({
            Foo % P: Foo,
            Bar % P: Foo % P,
            Foo: Bar,
            Bar: Bar % P
        })
        input_signature = dict(input1=X)
        output_signature = dict(output1=Y)
        foop = dict(input1=Foo % P)
        barp = dict(input1=Bar % P)
        foo = dict(input1=Foo)
        bar = dict(input1=Bar)

        self.assertEqual(meta.match(foop, input_signature, output_signature),
                         dict(output1=Foo))
        self.assertEqual(meta.match(barp, input_signature, output_signature),
                         dict(output1=Foo % P))
        self.assertEqual(meta.match(foo, input_signature, output_signature),
                         dict(output1=Bar))
        self.assertEqual(meta.match(bar, input_signature, output_signature),
                         dict(output1=Bar % P))

    def test_nested_match(self):
        C2 = MockTemplate('C2', fields=('a', 'b'))
        Foo = MockTemplate('Foo')
        Bar = MockTemplate('Bar')
        P = MockPredicate('P')
        X, Y = meta.TypeMap({
            Foo % P: Foo,
            Bar % P: Foo % P,
            Foo: Bar,
            Bar: Bar % P
        })
        input_signature = dict(input1=C2[X, Bar])
        output_signature = dict(output1=C2[Bar, Y])
        foop = dict(input1=C2[Foo % P, Bar])
        barp = dict(input1=C2[Bar % P, Bar])
        foo = dict(input1=C2[Foo, Foo])
        bar = dict(input1=C2[Bar, Foo])

        self.assertEqual(meta.match(foop, input_signature, output_signature),
                         dict(output1=C2[Bar, Foo]))
        self.assertEqual(meta.match(barp, input_signature, output_signature),
                         dict(output1=C2[Bar, Foo % P]))
        self.assertEqual(meta.match(foo, input_signature, output_signature),
                         dict(output1=C2[Bar, Bar]))
        self.assertEqual(meta.match(bar, input_signature, output_signature),
                         dict(output1=C2[Bar, Bar % P]))

    def test_multiple_variables(self):
        C2 = MockTemplate('C2', fields=('a', 'b'))
        Foo = MockTemplate('Foo')
        Bar = MockTemplate('Bar')
        P = MockPredicate('P')
        A, B, C, Y, Z = meta.TypeMap({
            (Foo % P, Bar, Bar): (Foo, Foo),
            (Foo, Bar % P, Foo): (Bar, Foo),
            (Foo, Foo, Bar): (Foo, Bar),
            (Bar, Bar % P, Foo): (Bar, Bar)
        })

        input_signature = dict(input1=C2[A, B], input2=C)
        output_signature = dict(output1=C2[Y, Z])

        fbb = dict(input1=C2[Foo % P, Bar], input2=Bar)
        fbf = dict(input1=C2[Foo, Bar % P], input2=Foo)
        ffb = dict(input1=C2[Foo, Foo], input2=Bar % P)  # subtype on in2!
        bbf = dict(input1=C2[Bar % P, Bar % P], input2=Foo)  # subtype on in1

        self.assertEqual(meta.match(fbb, input_signature, output_signature),
                         dict(output1=C2[Foo, Foo]))
        self.assertEqual(meta.match(fbf, input_signature, output_signature),
                         dict(output1=C2[Bar, Foo]))
        self.assertEqual(meta.match(ffb, input_signature, output_signature),
                         dict(output1=C2[Foo, Bar]))
        self.assertEqual(meta.match(bbf, input_signature, output_signature),
                         dict(output1=C2[Bar, Bar]))

    def test_multiple_mappings(self):
        C2 = MockTemplate('C2', fields=('a', 'b'))
        Foo = MockTemplate('Foo')
        Bar = MockTemplate('Bar')
        P = MockPredicate('P')

        X, Y = meta.TypeMap({
            Foo % P: Foo,
            Bar % P: Foo % P,
            Foo: Bar,
            Bar: Bar % P
        })

        T, R = meta.TypeMap({
            Bar % P: Foo,
            Foo % P: Foo % P,
            Bar: Bar,
            Foo: Bar % P
        })

        input_signature = dict(input1=C2[X, T])
        output_signature = dict(output1=C2[R, Y])
        foop = dict(input1=C2[Foo % P, Bar])
        barp = dict(input1=C2[Bar % P, Bar % P])
        foo = dict(input1=C2[Foo, Foo])
        bar = dict(input1=C2[Bar, Foo])

        self.assertEqual(meta.match(foop, input_signature, output_signature),
                         dict(output1=C2[Bar, Foo]))
        self.assertEqual(meta.match(barp, input_signature, output_signature),
                         dict(output1=C2[Foo, Foo % P]))
        self.assertEqual(meta.match(foo, input_signature, output_signature),
                         dict(output1=C2[Bar % P, Bar]))
        self.assertEqual(meta.match(bar, input_signature, output_signature),
                         dict(output1=C2[Bar % P, Bar % P]))

    def test_no_solution(self):
        C2 = MockTemplate('C2', fields=('a', 'b'))
        Foo = MockTemplate('Foo')
        Bar = MockTemplate('Bar')
        P = MockPredicate('P')
        A, B, C, Y, Z = meta.TypeMap({
            (Foo % P, Bar, Bar): (Foo, Foo),
            (Foo, Bar % P, Foo): (Bar, Foo),
            (Foo, Foo, Bar): (Foo, Bar),
            (Bar, Bar % P, Foo): (Bar, Bar)
        })

        input_signature = dict(input1=C2[A, B], input2=C)
        output_signature = dict(output1=C2[Y, Z])

        with self.assertRaisesRegex(ValueError, 'No solution'):
            meta.match(dict(input1=C2[Foo, Foo], input2=Foo),
                       input_signature, output_signature)

    def test_inconsistent_binding(self):
        C2 = MockTemplate('C2', fields=('a', 'b'))
        Foo = MockTemplate('Foo')
        Bar = MockTemplate('Bar')
        P = MockPredicate('P')
        A, B, C, Y, Z = meta.TypeMap({
            (Foo % P, Bar, Bar): (Foo, Foo),
            (Foo, Bar % P, Foo): (Bar, Foo),
            (Foo, Foo, Bar): (Foo, Bar),
            (Bar, Bar % P, Foo): (Bar, Bar)
        })

        input_signature = dict(input1=C2[A, B], input2=C2[C, C])
        output_signature = dict(output1=C2[Y, Z])

        with self.assertRaisesRegex(ValueError, 'to match'):
            meta.match(dict(input1=C2[Foo, Bar % P], input2=C2[Foo, Bar]),
                       input_signature, output_signature)

    def test_consistent_subtype_binding(self):
        C2 = MockTemplate('C2', fields=('a', 'b'))
        Foo = MockTemplate('Foo')
        Bar = MockTemplate('Bar')
        P = MockPredicate('P')
        A, B, C, Y, Z = meta.TypeMap({
            (Foo % P, Bar, Bar): (Foo, Foo),
            (Foo, Bar % P, Foo): (Bar, Foo),
            (Foo, Foo, Bar): (Foo, Bar),
            (Bar, Bar % P, Foo): (Bar, Bar)
        })

        input_signature = dict(input1=C2[A, B], input2=C2[C, C])
        output_signature = dict(output1=C2[Y, Z])

        cons = dict(input1=C2[Foo, Bar % P], input2=C2[Foo, Foo % P])

        self.assertEqual(meta.match(cons, input_signature, output_signature),
                         dict(output1=C2[Bar, Foo]))

    def test_missing_variables(self):
        C2 = MockTemplate('C2', fields=('a', 'b'))
        Foo = MockTemplate('Foo')
        Bar = MockTemplate('Bar')
        P = MockPredicate('P')
        A, B, C, Y, Z = meta.TypeMap({
            (Foo % P, Bar, Bar): (Foo, Foo),
            (Foo, Bar % P, Foo): (Bar, Foo),
            (Foo, Foo, Bar): (Foo, Bar),
            (Bar, Bar % P, Foo): (Bar, Bar)
        })

        input_signature = dict(input1=C2[A, B], input2=Foo)
        output_signature = dict(output1=C2[Y, Z])

        with self.assertRaisesRegex(ValueError, 'Missing'):
            meta.match(dict(input1=C2[Foo, Foo], input2=Foo),
                       input_signature, output_signature)

    def test_no_variables(self):
        Foo = MockTemplate('Foo')
        Bar = MockTemplate('Bar')
        P = MockPredicate('P')

        input_signature = dict(input1=Foo, input2=Bar)
        output_signature = dict(output1=Bar % P, output2=Foo % P)

        given = dict(input1=Foo % P, input2=Bar)

        self.assertEqual(meta.match(given, input_signature, output_signature),
                         output_signature)

    def test_type_match(self):
        Foo = MockTemplate('Foo')
        Bar = MockTemplate('Bar')
        Baz = MockTemplate('Baz')
        P = MockPredicate('P')

        T = meta.TypeMatch([Baz, Foo, Bar])

        input_signature = dict(input1=T)
        output_signature = dict(output1=T)

        foop = dict(input1=Foo % P)
        barp = dict(input1=Bar % P)
        foo = dict(input1=Foo)
        bar = dict(input1=Bar)

        self.assertEqual(meta.match(foop, input_signature, output_signature),
                         dict(output1=Foo))
        self.assertEqual(meta.match(barp, input_signature, output_signature),
                         dict(output1=Bar))
        self.assertEqual(meta.match(foo, input_signature, output_signature),
                         dict(output1=Foo))
        self.assertEqual(meta.match(bar, input_signature, output_signature),
                         dict(output1=Bar))

    def test_type_match_auto_intersect(self):
        C1 = MockTemplate('C1', fields=('a',))
        Foo = MockTemplate('Foo')
        P = MockPredicate('P')
        Q = MockPredicate('Q')
        R = MockPredicate('R')
        S = MockPredicate('S')

        T = meta.TypeMatch([P, Q, R, S])

        input_signature = dict(input1=C1[Foo] % T)
        output_signature = dict(output1=Foo % T)

        pqrs = dict(input1=C1[Foo] % (P & Q & R & S))
        p = dict(input1=C1[Foo] % P)
        pr = dict(input1=C1[Foo] % (P & R))
        qs = dict(input1=C1[Foo] % (Q & S))

        self.assertEqual(meta.match(pqrs, input_signature, output_signature),
                         dict(output1=Foo % (P & Q & R & S)))
        self.assertEqual(meta.match(p, input_signature, output_signature),
                         dict(output1=Foo % P))
        self.assertEqual(meta.match(pr, input_signature, output_signature),
                         dict(output1=Foo % (P & R)))
        self.assertEqual(meta.match(qs, input_signature, output_signature),
                         dict(output1=Foo % (Q & S)))


if __name__ == '__main__':
    unittest.main()
