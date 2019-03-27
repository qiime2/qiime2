import unittest

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








if __name__ == '__main__':
    unittest.main()
