# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pickle
import unittest
import collections

import qiime2.core.type.grammar as grammar
import qiime2.core.type.template as template


class _MockBase:
    public_proxy = 'example',

    def __init__(self, name, fields=()):
        self.test_data = {}
        self.name = name
        self.fields = fields

    def track_call(func):
        def wrapped(self, *args, **kwargs):
            self.test_data[func.__name__] = True
            return func(self, *args, **kwargs)
        return wrapped

    @track_call
    def __eq__(self, other):
        return id(self) == id(other)

    @track_call
    def __hash__(self):
        return hash(id(self))

    @track_call
    def get_field_names(self):
        return self.fields

    @track_call
    def get_name(self):
        return self.name

    @track_call
    def get_kind(self):
        return "tester"

    @track_call
    def validate_union(self, other):
        pass

    @track_call
    def validate_intersection(self, other):
        pass

    @track_call
    def is_element(self, value):
        return bool(value)

    @track_call
    def collapse_intersection(self, other):
        return super().collapse_intersection(other)

    @track_call
    def is_symbol_subtype(self, other):
        return self.name == other.name

    @track_call
    def is_symbol_supertype(self, other):
        return self.name == other.name

    @track_call
    def update_ast(self, ast):
        ast['extra_junk'] = self.name

    def validate_field(self, name, field):
        self.test_data['validate_field'] = name
        if field.name == 'InvalidMember':
            raise TypeError('InvalidMember cannot be used')

    @track_call
    def validate_predicate(self, predicate):
        pass

    @track_call
    def example(self):
        return ...


class MockTemplate(_MockBase, template.TypeTemplate):
    pass


class MockPredicate(_MockBase, template.PredicateTemplate):
    def __init__(self, name, alphabetize=False):
        self.alphabetize = alphabetize
        super().__init__(name)

    def __repr__(self):
        return self.name

    def is_symbol_subtype(self, other):
        if not self.alphabetize:
            return super().is_symbol_subtype(other)

        return self.name <= other.name

    def is_symbol_supertype(self, other):
        if not self.alphabetize:
            return super().is_symbol_supertype(other)

        return self.name >= other.name


class TestIncompleteExp(unittest.TestCase):
    def IncompleteExp(self, name, fields):
        expr = MockTemplate(name, fields)
        self.assertIsInstance(expr, grammar.IncompleteExp)
        return expr

    def test_construction_sanity(self):
        expr = MockTemplate('foo')  # TypeExpr
        with self.assertRaisesRegex(ValueError, "no fields"):
            # template has no fields, so putting it in an IncompleteExp
            # doesn't make sense
            expr = grammar.IncompleteExp(expr.template)

    def test_mod(self):
        with self.assertRaisesRegex(TypeError, 'predicate'):
            self.IncompleteExp('foo', ('a',)) % ...

    def test_or(self):
        with self.assertRaisesRegex(TypeError, 'union'):
            self.IncompleteExp('foo', ('a',)) | ...

    def test_and(self):
        with self.assertRaisesRegex(TypeError, 'intersect'):
            self.IncompleteExp('foo', ('a',)) & ...

    def test_repr(self):
        self.assertEqual(repr(self.IncompleteExp('Example', ('foo',))),
                         'Example[{foo}]')

        self.assertEqual(repr(self.IncompleteExp('Example', ('f', 'b'))),
                         'Example[{f}, {b}]')

    def test_le(self):
        expr_a = self.IncompleteExp('Foo', ('a',))
        expr_b = self.IncompleteExp('Bar', ('b',))
        with self.assertRaisesRegex(TypeError, 'missing arguments'):
            expr_a <= expr_b

    def test_ge(self):
        expr_a = self.IncompleteExp('Foo', ('a',))
        expr_b = self.IncompleteExp('Bar', ('b',))
        with self.assertRaisesRegex(TypeError, 'missing arguments'):
            expr_a >= expr_b

    def test_in(self):
        expr_a = self.IncompleteExp('Foo', ('a',))
        with self.assertRaisesRegex(TypeError, 'missing arguments'):
            ... in expr_a

    def test_field_w_typeexp(self):
        expr_a = self.IncompleteExp('Foo', ('baz',))
        expr_inner = MockTemplate('Bar')

        result = expr_a[expr_inner]

        self.assertEqual(repr(result), 'Foo[Bar]')
        self.assertIsInstance(result, grammar.TypeExp)
        self.assertEqual(expr_a.template.test_data['validate_field'], 'baz')

    def test_field_w_incompleteexp(self):
        expr_a = self.IncompleteExp('Foo', ('a',))
        expr_b = self.IncompleteExp('Bar', ('b',))
        with self.assertRaisesRegex(TypeError, 'complete type expression'):
            expr_a[expr_b]

    def test_field_w_nonsense(self):
        expr_a = self.IncompleteExp('Foo', ('a',))
        with self.assertRaisesRegex(TypeError, 'complete type expression'):
            expr_a[...]

    def test_field_wrong_length(self):
        X = MockTemplate('X')
        C = self.IncompleteExp('C', ['foo', 'bar'])
        with self.assertRaisesRegex(TypeError, '1'):
            C[X]

        C = self.IncompleteExp('C', ['foo'])
        with self.assertRaisesRegex(TypeError, '2'):
            C[X, X]

    def test_field_nested_expression(self):
        X = MockTemplate('X')
        C = self.IncompleteExp('C', ['foo', 'bar'])
        self.assertEqual(repr(C[X, C[C[X, X], X]]), 'C[X, C[C[X, X], X]]')

    def test_field_invalid_member(self):
        C = self.IncompleteExp('C', ['foo'])
        InvalidMember = MockTemplate('InvalidMember')
        with self.assertRaisesRegex(TypeError, 'InvalidMember'):
            C[InvalidMember]

    def test_field_union(self):
        X = MockTemplate('X')
        Y = MockTemplate('Y')
        Z = MockTemplate('Z')
        C = self.IncompleteExp('C', ['foo'])

        result = C[X | Y | Z]

        self.assertEqual(repr(result), "C[X | Y | Z]")

    def test_field_invalid_union(self):
        X = MockTemplate('X')
        InvalidMember = MockTemplate('InvalidMember')
        Z = MockTemplate('Z')
        C = self.IncompleteExp('C', ['foo'])

        with self.assertRaisesRegex(TypeError, 'InvalidMember'):
            C[X | InvalidMember | Z]

    def test_field_insane(self):
        X = MockTemplate('X')
        Y = MockTemplate('Y')
        Z = MockTemplate('Z')
        InvalidIntersection = grammar.IntersectionExp(
            members=(MockTemplate('InvalidMember'), Y))
        C = self.IncompleteExp('C', ['foo'])

        with self.assertRaisesRegex(TypeError, 'InvalidMember'):
            C[X | InvalidIntersection | Z]

    def test_iter_symbols(self):
        expr = self.IncompleteExp('Example', ('foo',))

        self.assertEqual(list(expr.iter_symbols()), ['Example'])

    def test_is_concrete(self):
        expr = self.IncompleteExp('Example', ('foo',))

        self.assertFalse(expr.is_concrete())

    def test_pickle(self):
        expr = self.IncompleteExp('Example', ('foo',))

        clone = pickle.loads(pickle.dumps(expr))

        self.assertEqual(expr, clone)

    def test_proxy(self):
        expr = self.IncompleteExp('Example', ('foo',))

        self.assertIs(expr.example(), ...)
        self.assertTrue(expr.template.test_data['example'])

    def test_eq_nonsense(self):
        expr_a = self.IncompleteExp('Example', ('foo',))

        self.assertEqual(expr_a.__eq__(...), NotImplemented)

    def test_hash_eq_equals(self):
        expr_a = self.IncompleteExp('Example', ('foo',))
        expr_b = self.IncompleteExp('Example', ('foo',))

        self.assertEqual(hash(expr_a), hash(expr_b))
        self.assertEqual(expr_a, expr_b)
        self.assertTrue(expr_a.equals(expr_b))

    def test_not_hash_eq_equals_field_mismatch(self):
        expr_a = self.IncompleteExp('Example', ('foo',))
        expr_b = self.IncompleteExp('Example', ('something_else',))

        self.assertNotEqual(hash(expr_a), hash(expr_b))
        self.assertNotEqual(expr_a, expr_b)
        self.assertFalse(expr_a.equals(expr_b))


class TestTypeExp(unittest.TestCase):
    def test_hashable(self):
        X = MockTemplate('X')
        Y = MockTemplate('Y', fields=('a',))
        Z = MockTemplate('Z')
        P = MockPredicate('P')

        self.assertIsInstance(X, collections.Hashable)
        # There really shouldn't be a collision between these:
        self.assertNotEqual(hash(X), hash(Z % P))

        self.assertEqual(Y[X], Y[X])
        self.assertEqual(hash(Y[X]), hash(Y[X]))

    def test_eq_nonsense(self):
        X = MockTemplate('X')
        self.assertIs(X.__eq__(42), NotImplemented)
        self.assertFalse(X == 42)

    def test_eq_different_instances(self):
        X = MockTemplate('X')
        X_ = MockTemplate('X')
        self.assertIsNot(X, X_)
        self.assertEqual(X, X_)

    def test_field(self):
        X = MockTemplate('X')
        with self.assertRaisesRegex(TypeError, 'fields'):
            X['scikit-bio/assets/.no.gif']

        Y = MockTemplate('Y', fields=('foo',))[X]
        with self.assertRaisesRegex(TypeError, 'fields'):
            Y[';-)']

    def test_repr(self):
        Face0 = MockTemplate('(o_-)')
        Face1 = MockTemplate('-_-')
        Exclaim0 = MockTemplate('!')
        Exclaim1 = MockTemplate('!', fields=('a',))
        Exclaim2 = MockTemplate('!', fields=('a', 'b'))
        Face2 = MockPredicate('(o_o)')
        Face3 = grammar.IntersectionExp(
            (MockPredicate('='), MockPredicate('=')))  # repr -> "= & ="
        Face4 = grammar.UnionExp(
            (MockPredicate('<'), MockPredicate('<')))

        self.assertEqual(repr(Exclaim0), '!')
        self.assertEqual(repr(Exclaim1[Face1]), '![-_-]')
        self.assertEqual(repr(Exclaim2[Face1, Exclaim0]), '![-_-, !]')
        self.assertEqual(repr(Exclaim2[Face1, Exclaim0] % Face2),
                         '![-_-, !] % (o_o)')
        self.assertEqual(repr(Face0 % Face2), '(o_-) % (o_o)')
        self.assertEqual(repr(Face0 % Face3), '(o_-) % (= & =)')
        self.assertEqual(repr(Exclaim2[Face1, Exclaim0] % Face3),
                         '![-_-, !] % (= & =)')
        self.assertEqual(repr(Exclaim2[Face1, Exclaim0] % Face4),
                         '![-_-, !] % (< | <)')

    def test_full_predicate(self):
        expr = MockTemplate('Foo')
        predicate = MockPredicate('Bar')

        self.assertIs((expr % predicate).full_predicate, predicate)
        self.assertTrue(expr.full_predicate.is_top())

    def test_in(self):
        expr = MockTemplate('Foo')

        self.assertIn(True, expr)
        self.assertTrue(expr.template.test_data['is_element'])

        expr = MockTemplate('Bar') % MockPredicate('Baz')

        self.assertIn(True, expr)
        self.assertTrue(expr.template.test_data['is_element'])
        self.assertTrue(expr.predicate.template.test_data['is_element'])

    def test_not_in(self):
        expr = MockTemplate('Foo')

        self.assertNotIn(False, expr)

        expr = MockTemplate('Bar') % MockPredicate('Baz')

        self.assertNotIn(False, expr)

    def test_mod(self):
        Bar = MockTemplate('Bar')
        Baz = MockPredicate('Baz')

        noop = Bar % grammar.IntersectionExp()
        self.assertIs(Bar, noop)

        with self.assertRaisesRegex(TypeError, 'predicate'):
            (Bar % Baz) % Baz

        with self.assertRaisesRegex(TypeError, 'right-hand'):
            Baz % Bar

    def test_iter(self):
        Foo = MockTemplate('Foo')
        Bar = MockTemplate('Bar')
        C2 = MockTemplate('C2', fields=('a', 'b'))
        P = MockPredicate('P')
        Q = MockPredicate('Q')

        self.assertEqual(
            {
                 Foo,
                 C2[Foo, Foo],
                 C2[Foo, C2[Foo % (P & Q), Bar]],
                 C2[Foo, C2[Foo % (P & Q), Foo]],
                 C2[Foo, C2[Bar, Bar]],
                 C2[Foo, C2[Bar, Foo]]
            },
            set(
                Foo | C2[Foo, Foo | C2[Foo % (P & Q) | Bar, Bar | Foo]]
            )
        )

    def test_iter_symbols(self):
        Foo = MockTemplate('Foo')
        Bar = MockTemplate('Bar')
        C2 = MockTemplate('C2', fields=('a', 'b'))
        P = MockPredicate('P')
        Q = MockPredicate('Q')

        self.assertEqual(
            {'Foo', 'C2', 'Bar'},
            set((Foo | C2[Foo, Foo | C2[Foo % (P & Q) | Bar, Bar | Foo]]
                 ).iter_symbols()))

    def test_is_concrete(self):
        Foo = MockTemplate('Foo')
        Bar = MockTemplate('Bar')
        C2 = MockTemplate('C2', fields=('a', 'b'))
        P = MockPredicate('P')
        Q = MockPredicate('Q')

        self.assertTrue(Foo.is_concrete())
        self.assertTrue(C2[Foo, Bar].is_concrete())
        self.assertTrue((Foo % P).is_concrete())
        self.assertTrue((C2[Foo % P, Bar] % Q).is_concrete())

    def test_not_concrete(self):
        Foo = MockTemplate('Foo')
        Bar = MockTemplate('Bar')
        C2 = MockTemplate('C2', fields=('a', 'b'))
        P = MockPredicate('P')
        Q = MockPredicate('Q')
        AnnoyingToMake = grammar.TypeExp(
            Foo.template, predicate=grammar.UnionExp((P, Q)))

        self.assertFalse((Foo | Bar).is_concrete())
        self.assertFalse(C2[Foo | Bar, Bar].is_concrete())
        self.assertFalse(C2[Foo, Bar | Foo].is_concrete())
        self.assertFalse(AnnoyingToMake.is_concrete())

class TestIntersection(unittest.TestCase):
    def test_basic(self):
        P = MockPredicate('P')
        Q = MockPredicate('Q')

        result = P & Q
        self.assertEqual(repr(result), "P & Q")

    def test_subtype(self):
        P = MockPredicate('P', alphabetize=True)
        Q = MockPredicate('Q', alphabetize=True)

        self.assertIs(Q & P, P)
        self.assertIs(P & Q, P)


#    def test_validate_union_w_nonsense(self):
#        X = grammar.TypeExp('X')
#        with self.assertRaisesRegex(TypeError, 'expression'):
#            X._validate_union_(42)
#
#    def test_validate_union_w_composite_type(self):
#        X = grammar.TypeExp('X')
#        with self.assertRaisesRegex(TypeError, 'incomplete'):
#            X._validate_union_(grammar.CompositeType('A', field_names=('X',)))
#
#    def test_validate_union_w_valid(self):
#        X = grammar.TypeExp('X')
#        Y = grammar.TypeExp('Y')
#        X._validate_union_(Y)
#
#    def test_validate_union_implements_handshake(self):
#        local = {}
#        X = grammar.TypeExp('X')
#
#        class Example(grammar.TypeExp):
#            def _validate_union_(self, other):
#                local['other'] = other
#
#        X | Example('Example')
#        self.assertIs(local['other'], X)
#
#    def test_build_union(self):
#        X = grammar.TypeExp('X')
#        Y = grammar.TypeExp('Y')
#        union = X._build_union_(X, Y)
#        self.assertIsInstance(union, grammar.UnionExp)
#        self.assertEqual(union.members, (X, Y))
#
#    def test_validate_intersection_w_nonsense(self):
#        X = grammar.TypeExp('X')
#        with self.assertRaisesRegex(TypeError, 'expression'):
#            X._validate_intersection_(42)
#
#    def test_validate_intersection_w_composite_type(self):
#        X = grammar.TypeExp('X')
#        with self.assertRaisesRegex(TypeError, 'incomplete'):
#            X._validate_intersection_(
#                grammar.CompositeType('A', field_names=('X',)))
#
#    def test_validate_intersection_w_valid(self):
#        X = grammar.TypeExp('X')
#        Y = grammar.TypeExp('Y')
#        X._validate_intersection_(Y)
#
#    def test_validate_intersection_implements_handshake(self):
#        local = {}
#        X = grammar.TypeExp('X')
#
#        class Example(grammar.TypeExp):
#            def _validate_intersection_(self, other):
#                local['other'] = other
#
#        X & Example('Example')
#        self.assertIs(local['other'], X)
#
#    def test_build_intersection(self):
#        X = grammar.TypeExp('X')
#        Y = grammar.TypeExp('Y')
#        intersection = X._build_intersection_(X, Y)
#        # Should produce the bottom type:
#        self.assertEqual(intersection, grammar.UnionExp())
#
#    def test_validate_predicate_w_nonsense(self):
#        X = grammar.TypeExp('X')
#        with self.assertRaisesRegex(TypeError, 'predicate'):
#            X._validate_predicate_(42)
#
#    def test_validate_predicate_w_valid(self):
#        predicate = grammar.Predicate(True)
#        X = grammar.TypeExp('X')
#        X._validate_predicate_(predicate)
#        # Test passed.
#
#    def test_apply_predicate(self):
#        predicate = grammar.Predicate(True)
#        Y = grammar.TypeExp('Y')
#        X = grammar.TypeExp('X', fields=(Y,))
#
#        result = X._apply_predicate_(predicate)
#        self.assertIsInstance(result, grammar.TypeExp)
#        self.assertEqual(result.fields, (Y,))
#
#    def test_is_subtype_wrong_name(self):
#        Y = grammar.TypeExp('Y')
#        X = grammar.TypeExp('X')
#
#        self.assertIs(Y._is_subtype_(X), NotImplemented)
#        self.assertIs(X._is_subtype_(Y), NotImplemented)
#        self.assertFalse(Y >= X)
#        self.assertFalse(Y >= X)
#        self.assertFalse(X <= Y)
#        self.assertFalse(X <= Y)
#
#    def test_is_subtype_diff_fields(self):
#        F1 = grammar.TypeExp('F1')
#        F2 = grammar.TypeExp('F2')
#        X = grammar.TypeExp('X', fields=(F1,))
#        X_ = grammar.TypeExp('X', fields=(F2,))
#
#        self.assertFalse(X_._is_subtype_(X))
#        self.assertFalse(X._is_subtype_(X_))
#
#    def test_is_subtype_diff_predicates(self):
#        class Pred(grammar.Predicate):
#            def __init__(self, value):
#                self.value = value
#                super().__init__(value)
#
#            def _is_subtype_(self, other):
#                if isinstance(other, self.__class__):
#                    return self.value <= other.value
#                return NotImplemented
#
#        P1 = Pred(1)
#        P2 = Pred(2)
#        X = grammar.TypeExp('X', predicate=P1)
#        X_ = grammar.TypeExp('X', predicate=P2)
#
#        self.assertFalse(X_._is_subtype_(X))
#        self.assertTrue(X._is_subtype_(X_))
#
#    def test_is_subtype_matches(self):
#        X = grammar.TypeExp('X')
#        X_ = grammar.TypeExp('X')
#
#        self.assertTrue(X._is_subtype_(X))
#        self.assertTrue(X_._is_subtype_(X))
#        self.assertTrue(X._is_subtype_(X_))
#        self.assertTrue(X_._is_subtype_(X_))
#
#    def test_is_subtype_matches_w_fields(self):
#        F1 = grammar.TypeExp('F1')
#        F2 = grammar.TypeExp('F2')
#        X = grammar.TypeExp('X', fields=(F1,))
#        X_ = grammar.TypeExp('X', fields=(F2,))
#
#        self.assertFalse(X_._is_subtype_(X))
#        self.assertFalse(X._is_subtype_(X_))
#
#    def test_is_subtype_matches_w_predicate(self):
#        class Pred(grammar.Predicate):
#            def __init__(self, value=0):
#                self.value = value
#                super().__init__(value)
#
#            def _is_subtype_(self, other):
#                if isinstance(other, self.__class__):
#                    return self.value <= other.value
#                return NotImplemented
#
#        P1 = Pred(1)
#        P1_ = Pred(1)
#        X = grammar.TypeExp('X', predicate=P1)
#        X_ = grammar.TypeExp('X', predicate=P1_)
#
#        self.assertTrue(X._is_subtype_(X))
#        self.assertTrue(X_._is_subtype_(X))
#        self.assertTrue(X._is_subtype_(X_))
#        self.assertTrue(X_._is_subtype_(X_))
#
#
#class TestTypeExpressionMod(unittest.TestCase):
#    def setUp(self):
#        self.local = {}
#
#    def test_mod_w_existing_predicate(self):
#        X = grammar.TypeExp('X', predicate=grammar.Predicate('Truthy'))
#        with self.assertRaisesRegex(TypeError, 'predicate'):
#            X % grammar.Predicate('Other')
#
#    def test_mod_w_none_predicate(self):
#        X = grammar.TypeExp('X', predicate=None)
#        predicate = grammar.Predicate("Truthy")
#        self.assertIs((X % predicate).predicate, predicate)
#
#    def test_mod_w_none(self):
#        X = grammar.TypeExp('X')
#        self.assertEqual(X % None, X)
#
#    def test_validate_predicate_called(self):
#        class Example(grammar.TypeExp):
#            def _validate_predicate_(s, predicate):
#                self.local['predicate'] = predicate
#
#        example = Example('Example')
#        p = grammar.Predicate(...)
#        example % p
#        self.assertIs(self.local['predicate'], p)
#
#    def test_apply_predicate_called(self):
#        class Example(grammar.TypeExp):
#            def _validate_predicate_(s, predicate):
#                pass  # Let anything through
#
#            def _apply_predicate_(s, predicate):
#                self.local['predicate'] = predicate
#                return ...
#
#        example = Example('Example')
#        new_type_expr = example % 'Foo'
#        self.assertEqual(self.local['predicate'], 'Foo')
#        self.assertIs(new_type_expr, ...)
#
#
#class TestTypeExpressionOr(unittest.TestCase):
#    def setUp(self):
#        self.local = {}
#
#    def test_identity(self):
#        X = grammar.TypeExp('X')
#        X_ = grammar.TypeExp('X')
#        self.assertIs(X | X_, X)
#
#    def test_several(self):
#        X = grammar.TypeExp('X')
#        Y = grammar.TypeExp('Y')
#        Z = grammar.TypeExp('Z')
#
#        self.assertIsInstance(X | Y | Z, grammar.UnionExp)
#        self.assertEqual(X | Y | Z | X | Z, Y | Z | X)
#
#    def test_validate_union_called(self):
#        class Example(grammar.TypeExp):
#            def _validate_union_(s, other):
#                if 'other' in self.local:
#                    self.local['self'] = other
#                else:
#                    self.local['other'] = other
#
#        foo = Example('foo')
#        bar = Example('bar')
#        foo | bar
#        self.assertEqual(self.local['self'], foo)
#        self.assertEqual(self.local['other'], bar)
#
#    def test_build_union_called(self):
#        class Example(grammar.TypeExp):
#            def _validate_union_(s, other):
#                pass  # Let anything through
#
#            def _build_union_(s, *members):
#                self.local['members'] = members
#                return ...
#
#        foo = Example('foo')
#        bar = Example('bar')
#        new_type_expr = foo | bar
#
#        self.assertEqual(self.local['members'], (foo, bar))
#        self.assertIs(new_type_expr, ...)
#
#
#class TestTypeExpressionAnd(unittest.TestCase):
#    def setUp(self):
#        self.local = {}
#
#    def test_identity(self):
#        X = grammar.TypeExp('X')
#        X_ = grammar.TypeExp('X')
#        self.assertIs(X & X_, X_)
#
#    def test_several(self):
#        X = grammar.TypeExp('X')
#        Y = grammar.TypeExp('Y')
#        Z = grammar.TypeExp('Z')
#
#        self.assertEqual(X & Y & Z & X & Z, grammar.IntersectionExp())
#
#    def test_validate_intersection_called(self):
#        class Example(grammar.TypeExp):
#            def _validate_intersection_(s, other):
#                if 'other' in self.local:
#                    self.local['self'] = other
#                else:
#                    self.local['other'] = other
#
#        foo = Example('foo')
#        bar = Example('bar')
#        foo & bar
#        self.assertEqual(self.local['self'], foo)
#        self.assertEqual(self.local['other'], bar)
#
#    def test_build_intersection_called(self):
#        class Example(grammar.TypeExp):
#            def _validate_intersection_(s, other):
#                pass  # Let anything through
#
#            def _build_intersection_(s, *members):
#                self.local['members'] = members
#                return ...
#
#        foo = Example('foo')
#        bar = Example('bar')
#        new_type_expr = foo & bar
#
#        self.assertEqual(self.local['members'], (foo, bar))
#        self.assertIs(new_type_expr, ...)
#
#
#class TestTypeExpressionLE(unittest.TestCase):
#    def setUp(self):
#        self.local = {}
#
#    def test_is_subtype_called(self):
#        class Example(grammar.TypeExp):
#            def _is_subtype_(s, other):
#                self.local['other'] = other
#                return self.local['return']
#
#        example = Example('Example')
#        other = Example('Other')
#
#        self.local['return'] = True
#        result = example <= other
#        self.assertEqual(self.local['other'], other)
#        self.assertTrue(result)
#
#        self.local['return'] = False
#        result = example <= other
#        self.assertEqual(self.local['other'], other)
#        self.assertFalse(result)
#
#
#class TestTypeExpressionGE(unittest.TestCase):
#    def setUp(self):
#        self.local = {}
#
#    def test_is_subtype_called(self):
#        class Example(grammar.TypeExp):
#            def _is_subtype_(s, other):
#                self.local['other'] = other
#                return self.local['return']
#
#        example = Example('Example')
#        other = Example('Other')
#
#        self.local['return'] = True
#        result = example >= other
#        self.assertEqual(self.local['other'], example)
#        self.assertTrue(result)
#
#        self.local['return'] = False
#        result = example >= other
#        self.assertEqual(self.local['other'], example)
#        self.assertFalse(result)
#


if __name__ == '__main__':
    unittest.main()
