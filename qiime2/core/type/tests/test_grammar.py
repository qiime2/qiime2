# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import collections

import qiime2.core.type.grammar as grammar


class TestTypeBase(unittest.TestCase):
    def setUp(self):
        class Example(grammar._TypeBase):
            __getitem__ = __or__ = __and__ = lambda s, x: x

            def __eq__(self, other):
                return False

        self.Example = Example

    def test_ne(self):
        example = self.Example()
        self.assertNotEqual(example, 42)
        self.assertNotEqual(42, example)

    def test_rmod(self):
        example = self.Example()
        with self.assertRaisesRegex(TypeError, 'right-hand'):
            42 % example

    def test_rand(self):
        self.assertEqual(42 & self.Example(), 42)

    def test_ror(self):
        self.assertEqual(42 | self.Example(), 42)

    def test_delattr(self):
        example = self.Example()
        with self.assertRaisesRegex(TypeError, 'immutable'):
            del example.foo

    def test_setitem(self):
        example = self.Example()
        with self.assertRaisesRegex(TypeError, 'immutable'):
            example['foo'] = 1

    def test_delitem(self):
        example = self.Example()
        with self.assertRaisesRegex(TypeError, 'immutable'):
            del example['foo']

    def test_getitem(self):
        example = self.Example()
        self.assertEqual(example[1], 1)

    def test_freeze(self):
        example = self.Example()
        example.foo = 1
        self.assertEqual(example.foo, 1)

        example._freeze_()
        self.assertEqual(example.foo, 1)
        with self.assertRaisesRegex(TypeError, 'immutable'):
            example.foo = 1

        with self.assertRaisesRegex(TypeError, 'immutable'):
            example.bar = 1

    # These tests are not concerned with rewriting properties on the class,
    # that behaviour is left unspecified to match Python.


class TestCompositeType(unittest.TestCase):
    def test_immutable(self):
        # If this test fails, then the hiearchy has been rearranged and the
        # properties tested for `_TypeBase` should be tested for
        # this class.
        #     - Your Friendly Dead Man's Switch
        self.assertIsInstance(grammar.CompositeType('Example', ('foo',)),
                              grammar._TypeBase)

    def test_field_sanity(self):
        with self.assertRaisesRegex(ValueError, 'empty'):
            grammar.CompositeType('Example', ())

    def test_mod(self):
        with self.assertRaisesRegex(TypeError, 'predicate'):
            grammar.CompositeType('Example', ('foo',)) % None

    def test_or(self):
        with self.assertRaisesRegex(TypeError, 'union'):
            grammar.CompositeType('Example', ('foo',)) | None

    def test_and(self):
        with self.assertRaisesRegex(TypeError, 'intersect'):
            grammar.CompositeType('Example', ('foo',)) & None

    def test_repr(self):
        self.assertEqual(repr(grammar.CompositeType('Example', ('foo',))),
                         'Example[{foo}]')

        self.assertEqual(repr(grammar.CompositeType('Example', ('f', 'b'))),
                         'Example[{f}, {b}]')

    def test_validate_field_w_typeexp(self):
        Example = grammar.CompositeType('Example', ('foo',))
        # Check that no error is raised:
        Example._validate_field_('foo', grammar.TypeExpression('X'))
        # Test passed if we reach this line.

    def test_validate_field_w_comptype(self):
        Example = grammar.CompositeType('Example', ('foo',))
        with self.assertRaisesRegex(TypeError, 'Incomplete'):
            Example._validate_field_('foo', Example)

    def test_validate_field_w_nonsense(self):
        Example = grammar.CompositeType('Example', ('foo',))
        with self.assertRaisesRegex(TypeError, 'Ellipsis'):
            Example._validate_field_('foo', Ellipsis)

    def test_apply_fields(self):
        X = grammar.TypeExpression('X')
        Example = grammar.CompositeType('Example', ('foo',))

        result = Example._apply_fields_((X,))

        self.assertEqual(result.fields, (X,))
        self.assertEqual(result.name, 'Example')
        self.assertIsInstance(result, grammar.TypeExpression)

    def test_iter_symbols(self):
        Example = grammar.CompositeType('Example', ('foo',))

        self.assertEqual(list(Example.iter_symbols()), ['Example'])


class TestCompositeTypeGetItem(unittest.TestCase):
    def setUp(self):
        self.local = {}

    def test_wrong_length(self):
        X = grammar.TypeExpression('X')
        composite_type = grammar.CompositeType('C', ['foo', 'bar'])
        with self.assertRaisesRegex(TypeError, '1'):
            composite_type[X]

        composite_type = grammar.CompositeType('C', ['foo'])
        with self.assertRaisesRegex(TypeError, '2'):
            composite_type[X, X]

    def test_nested_expression(self):
        X = grammar.TypeExpression('X')
        C = grammar.CompositeType('C', ['foo', 'bar'])
        self.assertEqual(repr(C[X, C[C[X, X], X]]), 'C[X, C[C[X, X], X]]')

    def test_validate_field_called(self):
        class MyCompositeType(grammar.CompositeType):
            def _validate_field_(s, name, value):
                self.local['name'] = name
                self.local['value'] = value

        my_type = MyCompositeType('MyType', ['foo'])
        my_type[...]
        self.assertEqual(self.local['name'], 'foo')
        self.assertEqual(self.local['value'], ...)

    def test_apply_fields_called(self):
        class MyCompositeType(grammar.CompositeType):
            def _validate_field_(*args):
                pass  # Let anything through

            def _apply_fields_(s, fields):
                self.local['fields'] = fields
                return ...

        my_type = MyCompositeType('MyType', ['foo'])
        type_exp = my_type['!']  # '!' is not a `TypeExpression`
        self.assertEqual(self.local['fields'], ('!',))
        self.assertEqual(type_exp, ...)


class TestTypeExpression(unittest.TestCase):
    def test_immutable(self):
        # If this test fails, then the hiearchy has been rearranged and the
        # properties tested for `_TypeBase` should be tested for
        # this class.
        #     - Your Friendly Dead Man's Switch
        self.assertIsInstance(grammar.TypeExpression('X'),
                              grammar._TypeBase)

    def test_hashable(self):
        a = grammar.TypeExpression('X')
        b = grammar.TypeExpression('Y', fields=(a,))
        c = grammar.TypeExpression('Y', fields=(a,))
        d = grammar.TypeExpression('Z', predicate=grammar.Predicate())

        self.assertIsInstance(a, collections.Hashable)
        # There really shouldn't be a collision between these:
        self.assertNotEqual(hash(a), hash(d))

        self.assertEqual(b, c)
        self.assertEqual(hash(b), hash(c))

    # TODO: Test dictionaries work well

    def test_eq_nonsense(self):
        X = grammar.TypeExpression('X')
        self.assertIs(X.__eq__(42), NotImplemented)
        self.assertFalse(X == 42)

    def test_eq_different_instances(self):
        X = grammar.TypeExpression('X')
        X_ = grammar.TypeExpression('X')
        self.assertIsNot(X, X_)
        self.assertEqual(X, X_)

    # TODO: Add more equality tests

    def test_mod(self):
        X = grammar.TypeExpression('X')
        with self.assertRaisesRegex(TypeError, 'fields'):
            X['scikit-bio/assets/.no.gif']

        Y = grammar.TypeExpression('Y', fields=(X,))
        with self.assertRaisesRegex(TypeError, 'fields'):
            Y[';-)']

    def test_repr(self):
        # Subclass elements to demonstrate dispatch occurs correctly.
        class Face1(grammar.TypeExpression):
            def __repr__(self):
                return "-_-"

        class Exclaim(grammar.TypeExpression):
            def __repr__(self):
                return '!'

        class Face2(grammar.Predicate):
            def __repr__(self):
                return '(o_o)'

        self.assertEqual(
            repr(grammar.TypeExpression('!')),
            '!')
        self.assertEqual(
            repr(grammar.TypeExpression('!', fields=(Face1(''),))),
            '![-_-]')
        self.assertEqual(
            repr(grammar.TypeExpression('!',
                                        fields=(Face1(''), Exclaim('!')))),
            '![-_-, !]')
        self.assertEqual(
            repr(grammar.TypeExpression('!',
                                        fields=(Face1(''), Exclaim('!')),
                                        predicate=Face2(True))),
            '![-_-, !] % (o_o)')

        self.assertEqual(
            repr(grammar.TypeExpression('(o_-)',
                                        predicate=Face2(True))),
            '(o_-) % (o_o)')

    def test_validate_union_w_nonsense(self):
        X = grammar.TypeExpression('X')
        with self.assertRaisesRegex(TypeError, 'expression'):
            X._validate_union_(42)

    def test_validate_union_w_composite_type(self):
        X = grammar.TypeExpression('X')
        with self.assertRaisesRegex(TypeError, 'incomplete'):
            X._validate_union_(grammar.CompositeType('A', field_names=('X',)))

    def test_validate_union_w_valid(self):
        X = grammar.TypeExpression('X')
        Y = grammar.TypeExpression('Y')
        X._validate_union_(Y)

    def test_validate_union_implements_handshake(self):
        local = {}
        X = grammar.TypeExpression('X')

        class Example(grammar.TypeExpression):
            def _validate_union_(self, other, handshake=False):
                local['other'] = other
                local['handshake'] = handshake

        X._validate_union_(Example('Example'))
        self.assertIs(local['other'], X)
        self.assertTrue(local['handshake'])

    def test_build_union(self):
        X = grammar.TypeExpression('X')
        Y = grammar.TypeExpression('Y')
        union = X._build_union_((X, Y))
        self.assertIsInstance(union, grammar.UnionTypeExpression)
        self.assertEqual(union.members, frozenset({X, Y}))

    def test_validate_intersection_w_nonsense(self):
        X = grammar.TypeExpression('X')
        with self.assertRaisesRegex(TypeError, 'expression'):
            X._validate_intersection_(42)

    def test_validate_intersection_w_composite_type(self):
        X = grammar.TypeExpression('X')
        with self.assertRaisesRegex(TypeError, 'incomplete'):
            X._validate_intersection_(
                grammar.CompositeType('A', field_names=('X',)))

    def test_validate_intersection_w_valid(self):
        X = grammar.TypeExpression('X')
        Y = grammar.TypeExpression('Y')
        X._validate_intersection_(Y)

    def test_validate_intersection_implements_handshake(self):
        local = {}
        X = grammar.TypeExpression('X')

        class Example(grammar.TypeExpression):
            def _validate_intersection_(self, other, handshake=False):
                local['other'] = other
                local['handshake'] = handshake

        X._validate_intersection_(Example('Example'))
        self.assertIs(local['other'], X)
        self.assertTrue(local['handshake'])

    def test_build_intersection(self):
        X = grammar.TypeExpression('X')
        Y = grammar.TypeExpression('Y')
        intersection = X._build_intersection_((X, Y))
        self.assertIsInstance(intersection, grammar.IntersectionTypeExpression)
        self.assertEqual(intersection.members, frozenset({X, Y}))

    def test_validate_predicate_w_nonsense(self):
        X = grammar.TypeExpression('X')
        with self.assertRaisesRegex(TypeError, 'predicate'):
            X._validate_predicate_(42)

    def test_validate_predicate_w_valid(self):
        predicate = grammar.Predicate(True)
        X = grammar.TypeExpression('X')
        X._validate_predicate_(predicate)
        # Test passed.

    def test_apply_predicate(self):
        predicate = grammar.Predicate(True)
        Y = grammar.TypeExpression('Y')
        X = grammar.TypeExpression('X', fields=(Y,))

        result = X._apply_predicate_(predicate)
        self.assertIsInstance(result, grammar.TypeExpression)
        self.assertEqual(result.fields, (Y,))

    def test_is_subtype_wrong_name(self):
        Y = grammar.TypeExpression('Y')
        X = grammar.TypeExpression('X')

        self.assertFalse(Y._is_subtype_(X))
        self.assertFalse(X._is_subtype_(Y))

    def test_is_subtype_diff_fields(self):
        F1 = grammar.TypeExpression('F1')
        F2 = grammar.TypeExpression('F2')
        X = grammar.TypeExpression('X', fields=(F1,))
        X_ = grammar.TypeExpression('X', fields=(F2,))

        self.assertFalse(X_._is_subtype_(X))
        self.assertFalse(X._is_subtype_(X_))

    def test_is_subtype_diff_predicates(self):
        class Pred(grammar.Predicate):
            def __init__(self, value):
                self.value = value
                super().__init__(value)

            def _is_subtype_(self, other):
                return self.value <= other.value

        P1 = Pred(1)
        P2 = Pred(2)
        X = grammar.TypeExpression('X', predicate=P1)
        X_ = grammar.TypeExpression('X', predicate=P2)

        self.assertFalse(X_._is_subtype_(X))
        self.assertTrue(X._is_subtype_(X_))

    def test_is_subtype_matches(self):
        X = grammar.TypeExpression('X')
        X_ = grammar.TypeExpression('X')

        self.assertTrue(X._is_subtype_(X))
        self.assertTrue(X_._is_subtype_(X))
        self.assertTrue(X._is_subtype_(X_))
        self.assertTrue(X_._is_subtype_(X_))

    def test_is_subtype_matches_w_fields(self):
        F1 = grammar.TypeExpression('F1')
        F2 = grammar.TypeExpression('F2')
        X = grammar.TypeExpression('X', fields=(F1,))
        X_ = grammar.TypeExpression('X', fields=(F2,))

        self.assertFalse(X_._is_subtype_(X))
        self.assertFalse(X._is_subtype_(X_))

    def test_is_subtype_matches_w_predicate(self):
        class Pred(grammar.Predicate):
            def __init__(self, value=0):
                self.value = value
                super().__init__(value)

            def _is_subtype_(self, other):
                return self.value <= other.value

        P1 = Pred(1)
        P1_ = Pred(1)
        X = grammar.TypeExpression('X', predicate=P1)
        X_ = grammar.TypeExpression('X', predicate=P1_)

        self.assertTrue(X._is_subtype_(X))
        self.assertTrue(X_._is_subtype_(X))
        self.assertTrue(X._is_subtype_(X_))
        self.assertTrue(X_._is_subtype_(X_))


class TestTypeExpressionMod(unittest.TestCase):
    def setUp(self):
        self.local = {}

    def test_mod_w_existing_predicate(self):
        X = grammar.TypeExpression('X', predicate=grammar.Predicate('Truthy'))
        with self.assertRaisesRegex(TypeError, 'predicate'):
            X % grammar.Predicate('Other')

    def test_mod_w_falsy_predicate(self):
        X = grammar.TypeExpression('X', predicate=grammar.Predicate())
        predicate = grammar.Predicate("Truthy")
        self.assertIs((X % predicate).predicate, predicate)

    def test_mod_w_none(self):
        X = grammar.TypeExpression('X')
        self.assertEqual(X % None, X)

    def test_validate_predicate_called(self):
        class Example(grammar.TypeExpression):
            def _validate_predicate_(s, predicate):
                self.local['predicate'] = predicate

        example = Example('Example')
        example % 42
        self.assertEqual(self.local['predicate'], 42)

    def test_apply_predicate_called(self):
        class Example(grammar.TypeExpression):
            def _validate_predicate_(s, predicate):
                pass  # Let anything through

            def _apply_predicate_(s, predicate):
                self.local['predicate'] = predicate
                return ...

        example = Example('Example')
        new_type_expr = example % 'Foo'
        self.assertEqual(self.local['predicate'], 'Foo')
        self.assertIs(new_type_expr, ...)


class TestTypeExpressionOr(unittest.TestCase):
    def setUp(self):
        self.local = {}

    def test_identity(self):
        X = grammar.TypeExpression('X')
        X_ = grammar.TypeExpression('X')
        self.assertIs(X | X_, X)

    def test_several(self):
        X = grammar.TypeExpression('X')
        Y = grammar.TypeExpression('Y')
        Z = grammar.TypeExpression('Z')

        self.assertIsInstance(X | Y | Z, grammar.UnionTypeExpression)
        self.assertEqual(X | Y | Z | X | Z, Y | Z | X)

    def test_validate_union_called(self):
        class Example(grammar.TypeExpression):
            def _validate_union_(s, other, handshake):
                self.local['other'] = other
                self.local['handshake'] = handshake

        example = Example('Example')
        example | 42
        self.assertEqual(self.local['other'], 42)
        self.assertFalse(self.local['handshake'])

    def test_build_union_called(self):
        class Example(grammar.TypeExpression):
            def _validate_union_(s, other, handshake):
                pass  # Let anything through

            def _build_union_(s, members):
                self.local['members'] = members
                return ...

        example = Example('Example')
        new_type_expr = example | 42
        self.assertEqual(self.local['members'], (example, 42))
        self.assertIs(new_type_expr, ...)


class TestTypeExpressionAnd(unittest.TestCase):
    def setUp(self):
        self.local = {}

    def test_identity(self):
        X = grammar.TypeExpression('X')
        X_ = grammar.TypeExpression('X')
        self.assertIs(X & X_, X_)

    def test_several(self):
        X = grammar.TypeExpression('X')
        Y = grammar.TypeExpression('Y')
        Z = grammar.TypeExpression('Z')

        self.assertIsInstance(X & Y & Z, grammar.IntersectionTypeExpression)
        self.assertEqual(X & Y & Z & X & Z, Y & Z & X)

    def test_validate_intersection_called(self):
        class Example(grammar.TypeExpression):
            def _validate_intersection_(s, other, handshake):
                self.local['other'] = other
                self.local['handshake'] = handshake

        example = Example('Example')
        example & 42
        self.assertEqual(self.local['other'], 42)
        self.assertFalse(self.local['handshake'])

    def test_build_intersection_called(self):
        class Example(grammar.TypeExpression):
            def _validate_intersection_(s, other, handshake):
                pass  # Let anything through

            def _build_intersection_(s, members):
                self.local['members'] = members
                return ...

        example = Example('Example')
        new_type_expr = example & 42
        self.assertEqual(self.local['members'], (example, 42))
        self.assertIs(new_type_expr, ...)


class TestTypeExpressionLE(unittest.TestCase):
    def setUp(self):
        self.local = {}

    def test_is_subtype_called(self):
        class Example(grammar.TypeExpression):
            def _is_subtype_(s, other):
                self.local['other'] = other
                return self.local['return']

        example = Example('Example')
        other = Example('Other')

        self.local['return'] = True
        result = example <= other
        self.assertEqual(self.local['other'], other)
        self.assertTrue(result)

        self.local['return'] = False
        result = example <= other
        self.assertEqual(self.local['other'], other)
        self.assertFalse(result)


class TestTypeExpressionGE(unittest.TestCase):
    def setUp(self):
        self.local = {}

    def test_is_subtype_called(self):
        class Example(grammar.TypeExpression):
            def _is_subtype_(s, other):
                self.local['other'] = other
                return self.local['return']

        example = Example('Example')
        other = Example('Other')

        self.local['return'] = True
        result = example >= other
        self.assertEqual(self.local['other'], example)
        self.assertTrue(result)

        self.local['return'] = False
        result = example >= other
        self.assertEqual(self.local['other'], example)
        self.assertFalse(result)


# TODO: test the following:
# - _SetOperationBase
# - UnionTypeExpression
# - IntersectionTypeExpression
# - MappingTypeExpression
# - Predicate


if __name__ == '__main__':
    unittest.main()
