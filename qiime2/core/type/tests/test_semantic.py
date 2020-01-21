# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
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
            if isinstance(element, grammar._ExpBase):
                self.assertFalse(semantic.is_semantic_type(element))
        self.assertTrue(looped)

    def test_visualization_not_semantic(self):
        self.assertFalse(
            semantic.is_semantic_type(visualization.Visualization))

    def test_type_expr_not_semantic(self):
        TypeExpr = grammar.TypeExp(None)
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


if __name__ == '__main__':
    unittest.main()
