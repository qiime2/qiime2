# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import qiime2.core.type.semantic as semantic
import qiime2.core.type.grammar as grammar
import qiime2.core.type.primitive as primitive
import qiime2.core.type.visualization as visualization


class TestIsSemanticType(unittest.TestCase):
    def test_primitives_not_semantic(self):
        looped = False
        for element in dir(primitive):
            looped = True
            element = getattr(primitive, element)
            self.assertFalse(semantic.is_semantic_type(element))
        self.assertTrue(looped)

    def test_visualization_not_semantic(self):
        self.assertFalse(
            semantic.is_semantic_type(visualization.Visualization))

    def test_type_expr_not_semantic(self):
        TypeExpr = grammar.TypeExpression('TypeExpr')
        self.assertFalse(semantic.is_semantic_type(TypeExpr))

    def test_simple_semantic_type(self):
        A = semantic.SemanticType('A')
        X = semantic.SemanticType('X')
        Foo = semantic.SemanticType('Foo', field_names=['a', 'b'])

        self.assertTrue(semantic.is_semantic_type(A))
        self.assertTrue(semantic.is_semantic_type(X))
        self.assertTrue(semantic.is_semantic_type(Foo))

    def test_composite_semantic_type(self):
        Foo = semantic.SemanticType('Foo', field_names=['a', 'b'])
        A = semantic.SemanticType('A', variant_of=Foo.field['a'])
        B = semantic.SemanticType('B', variant_of=Foo.field['b'])

        self.assertTrue(semantic.is_semantic_type(A))
        self.assertTrue(semantic.is_semantic_type(B))
        self.assertTrue(semantic.is_semantic_type(Foo[A, B]))

#
# class TestSemanticTypeFactory(unittest.TestCase):
#     def test_name(self):
#         pass
#
#     def test_name_bad_input(self):
#         pass
#
#     def test_field_names_singlestr(self):
#         pass
#
#     def test_field_names_multistr(self):
#         pass
#
#     def test_field_names_bad_input(self):
#         pass
#
#     def test_field_members_no_field_names(self):
#         pass
#
#     def test_field_members_mismatch_field_names(self):
#         pass
#
#     def test_field_members_single_type(self):
#         pass
#
#     def test_field_members_mixed_type(self):
#         pass
#
#     def test_field_members_bad_input(self):
#         pass
#
#     def test_variant_of_single(self):
#         pass
#
#     def test_variant_of_multi(self):
#         pass
#
#     def test_variant_of_bad_input(self):
#         pass
#
#
# class TestVariantField(unittest.TestCase):
#     pass


class TestIncompleteSemanticType(unittest.TestCase):
    def test_composite_type(self):
        # If this test fails, then the hiearchy has been rearranged and the
        # properties tested for `CompositeType` should be tested for
        # this class.
        #     - Your Friendly Dead Man's Switch
        self.assertIsInstance(semantic.SemanticType('X', field_names='foo'),
                              grammar.CompositeType)

        self.assertIsInstance(
            semantic._IncompleteSemanticType('X', field_names=('foo',),
                                             field_members={'foo': ()},
                                             variant_of=()),
            grammar.CompositeType)

    # def test_validate_field(self):
    #     pass
    #
    # def test_apply_fields(self):
    #     pass


class TestSemanticType(unittest.TestCase):
    def test_composite_type(self):
        # If this test fails, then the hiearchy has been rearranged and the
        # properties tested for `TypeExpression` should be tested for
        # this class.
        #     - Your Friendly Dead Man's Switch
        self.assertIsInstance(semantic.SemanticType('X'),
                              grammar.TypeExpression)

        self.assertIsInstance(semantic._SemanticType('X', variant_of=()),
                              grammar.TypeExpression)


class TestSemanticUnionType(unittest.TestCase):
    def test_union_type_expr(self):
        # If this test fails, then the hiearchy has been rearranged and the
        # properties tested for `UnionTypeExpression` should be tested for
        # this class.
        #     - Your Friendly Dead Man's Switch
        X = semantic.SemanticType('X')
        Y = semantic.SemanticType('Y')
        self.assertIsInstance(X | Y, grammar.UnionTypeExpression)

        self.assertIsInstance(semantic._SemanticUnionType([X, Y]),
                              grammar.UnionTypeExpression)


class TestProperties(unittest.TestCase):
    def test_predicate(self):
        # If this test fails, then the hiearchy has been rearranged and the
        # properties tested for `Predicate` should be tested for
        # this class.
        #     - Your Friendly Dead Man's Switch
        self.assertIsInstance(semantic.Properties('X'), grammar.Predicate)

#
# class TestComplicatedExpression(unittest.TestCase):
#     def test_1(self):
#         pass
#
#     def test_2(self):
#         pass
#


if __name__ == '__main__':
    unittest.main()
