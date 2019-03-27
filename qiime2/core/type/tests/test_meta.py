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
            P: Bar,
            Q: Foo,
            P & Q: Foo
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
            P: Bar,
            Q: Foo,
            P & Q: Foo
        })

        expr = C2[X1, Foo % X2] % X2
        pred_sel, field_sel, field_pred_sel =  meta.select_variables(expr)

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
            P: P, Q: Q,
            Other: Other, P & Other: Other, Q & Other: Other
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
            P: R,
            Q: S,
            P & Q: R & S
        })

        self.assertEqual(X.members, (P & Q, P, Q))
        self.assertEqual(Y.members, (R & S, R, S))
        self.assertEqual(X.index, 0)
        self.assertEqual(Y.index, 1)
        self.assertTrue(X.input)
        self.assertTrue(Y.output)
        self.assertFalse(X.output)
        self.assertFalse(Y.input)

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
            (P, P): R,
            (Q, Q): S,
            (P & Q, Q): R & S,
            (P & Q, P & Q): R & S
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
    pass

if __name__ == '__main__':
    unittest.main()
