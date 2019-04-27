# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from qiime2.core.type import (
    parse_primitive, Int, Float, Bool, Str, List, Set)


class TestParsePrimitiveNonCollectionsSimple(unittest.TestCase):
    def test_int_type_int_value(self):
        obs = parse_primitive(Int, '42')
        self.assertEqual(obs, 42)
        self.assertIsInstance(obs, int)

    def test_float_type_int_value(self):
        obs = parse_primitive(Float, '42')
        self.assertEqual(obs, 42.0)
        self.assertIsInstance(obs, float)

    def test_bool_type_int_value(self):
        with self.assertRaisesRegex(ValueError, 'walk the plank'):
            parse_primitive(Bool, '42')

    def test_str_type_int_value(self):
        obs = parse_primitive(Str, '42')
        self.assertEqual(obs, '42')
        self.assertIsInstance(obs, str)

    def test_int_type_float_value(self):
        with self.assertRaisesRegex(ValueError, 'walk the plank'):
            parse_primitive(Int, '42.0')

    def test_float_type_float_value(self):
        obs = parse_primitive(Float, '42.0')
        self.assertEqual(obs, 42.0)
        self.assertIsInstance(obs, float)

    def test_bool_type_float_value(self):
        with self.assertRaisesRegex(ValueError, 'walk the plank'):
            parse_primitive(Bool, '42.0')

    def test_str_type_float_value(self):
        obs = parse_primitive(Str, '42.0')
        self.assertEqual(obs, '42.0')
        self.assertIsInstance(obs, str)

    def test_int_type_bool_value(self):
        with self.assertRaisesRegex(ValueError, 'walk the plank'):
            parse_primitive(Int, 'True')

    def test_float_type_bool_value(self):
        with self.assertRaisesRegex(ValueError, 'walk the plank'):
            parse_primitive(Float, 'True')

    def test_bool_type_bool_value(self):
        obs = parse_primitive(Bool, 'True')
        self.assertEqual(obs, True)
        self.assertIsInstance(obs, bool)

    def test_str_type_bool_value(self):
        obs = parse_primitive(Str, 'True')
        self.assertEqual(obs, 'True')
        self.assertIsInstance(obs, str)

    def test_int_type_str_value(self):
        with self.assertRaisesRegex(ValueError, 'walk the plank'):
            parse_primitive(Int, 'peanut')

    def test_float_type_str_value(self):
        with self.assertRaisesRegex(ValueError, 'walk the plank'):
            parse_primitive(Float, 'peanut')

    def test_bool_type_str_value(self):
        with self.assertRaisesRegex(ValueError, 'walk the plank'):
            parse_primitive(Bool, 'peanut')

    def test_str_type_str_value(self):
        obs = parse_primitive(Str, 'peanut')
        self.assertEqual(obs, 'peanut')
        self.assertIsInstance(obs, str)


class TestParsePrimitiveNonCollectionsSimpleUnions(unittest.TestCase):
    def setUp(self):
        super().setUp()

        self.exprs = [
            Int | Bool,
            Int | Str,
            Float | Bool,
            Float | Str,
            Bool | Str,
        ]

    def test_int_union_float_expr_int_value(self):
        # Int | Float == Float
        obs = parse_primitive(Int | Float, '42')
        self.assertEqual(obs, 42.0)
        self.assertIsInstance(obs, float)

    def test_int_union_float_expr_float_value(self):
        # Int | Float == Float
        obs = parse_primitive(Int | Float, '42.0')
        self.assertEqual(obs, 42.0)
        self.assertIsInstance(obs, float)

    def test_int_union_float_expr_bool_value(self):
        with self.assertRaisesRegex(ValueError, 'walk the plank'):
            # Int | Float == Float
            parse_primitive(Int | Float, 'True')

    def test_int_union_float_expr_str_value(self):
        with self.assertRaisesRegex(ValueError, 'walk the plank'):
            # Int | Float == Float
            parse_primitive(Int | Float, 'peanut')

    def test_simple_unions_with_int_value(self):
        for expr in self.exprs:
            with self.subTest(expr=expr):
                obs = parse_primitive(expr, '42')
                self.assertEqual(obs, 42)
                self.assertIsInstance(obs, int)

    def test_simple_unions_with_float_value(self):
        for expr in self.exprs:
            with self.subTest(expr=expr):
                obs = parse_primitive(expr, '42.1')
                self.assertEqual(obs, 42.1)
                self.assertIsInstance(obs, float)

    def test_simple_unions_with_bool_value(self):
        for expr in self.exprs:
            with self.subTest(expr=expr):
                obs = parse_primitive(expr, 'True')
                self.assertEqual(obs, True)
                self.assertIsInstance(obs, bool)

    def test_simple_unions_with_str_value(self):
        for expr in self.exprs:
            with self.subTest(expr=expr):
                obs = parse_primitive(expr, 'peanut')
                self.assertEqual(obs, 'peanut')
                self.assertIsInstance(obs, str)


class TestParsePrimitiveCollections(unittest.TestCase):
    def test_simple_int(self):
        pass

    def test_simple_float(self):
        pass

    def test_simple_bool(self):
        pass

    def test_simple_str(self):
        pass

    def test_monomorphic_int_or_float(self):
        pass

    def test_monomorphic_int_or_bool(self):
        pass

    def test_monomorphic_int_or_str(self):
        pass

    def test_monomorphic_float_or_bool(self):
        pass

    def test_monomorphic_float_or_str(self):
        pass

    def test_monomorphic_bool_or_str(self):
        pass

    def test_composite_int_or_float(self):
        pass

    def test_composite_int_or_bool(self):
        pass

    def test_composite_int_or_str(self):
        pass

    def test_composite_float_or_bool(self):
        pass

    def test_composite_float_or_str(self):
        pass

    def test_composite_bool_or_str(self):
        pass

    def test_complex_int_float_bool(self):
        pass

    def test_complex_int_float_str(self):
        pass

    def test_complex_float_bool_str(self):
        pass


if __name__ == '__main__':
    unittest.main()
