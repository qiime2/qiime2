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
        obs = parse_primitive(Bool, '42')
        self.assertEqual(obs, '42')
        self.assertIsInstance(obs, str)

    def test_str_type_int_value(self):
        obs = parse_primitive(Str, '42')
        self.assertEqual(obs, '42')
        self.assertIsInstance(obs, str)

    def test_int_type_float_value(self):
        obs = parse_primitive(Int, '42.0')
        self.assertEqual(obs, '42.0')
        self.assertIsInstance(obs, str)

    def test_float_type_float_value(self):
        obs = parse_primitive(Float, '42.0')
        self.assertEqual(obs, 42.0)
        self.assertIsInstance(obs, float)

    def test_bool_type_float_value(self):
        obs = parse_primitive(Bool, '42.0')
        self.assertEqual(obs, '42.0')
        self.assertIsInstance(obs, str)

    def test_str_type_float_value(self):
        obs = parse_primitive(Str, '42.0')
        self.assertEqual(obs, '42.0')
        self.assertIsInstance(obs, str)

    def test_int_type_bool_value(self):
        obs = parse_primitive(Int, 'True')
        self.assertEqual(obs, 'True')
        self.assertIsInstance(obs, str)

    def test_float_type_bool_value(self):
        obs = parse_primitive(Float, 'True')
        self.assertEqual(obs, 'True')
        self.assertIsInstance(obs, str)

    def test_bool_type_bool_value(self):
        obs = parse_primitive(Bool, 'True')
        self.assertEqual(obs, True)
        self.assertIsInstance(obs, bool)

    def test_str_type_bool_value(self):
        obs = parse_primitive(Str, 'True')
        self.assertEqual(obs, 'True')
        self.assertIsInstance(obs, str)

    def test_int_type_str_value(self):
        obs = parse_primitive(Int, 'peanut')
        self.assertEqual(obs, 'peanut')
        self.assertIsInstance(obs, str)

    def test_float_type_str_value(self):
        obs = parse_primitive(Float, 'peanut')
        self.assertEqual(obs, 'peanut')
        self.assertIsInstance(obs, str)

    def test_bool_type_str_value(self):
        obs = parse_primitive(Bool, 'peanut')
        self.assertEqual(obs, 'peanut')
        self.assertIsInstance(obs, str)

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
        # Int | Float == Float
        obs = parse_primitive(Int | Float, 'True')
        self.assertEqual(obs, 'True')
        self.assertIsInstance(obs, str)

    def test_int_union_float_expr_str_value(self):
        # Int | Float == Float
        obs = parse_primitive(Int | Float, 'peanut')
        self.assertEqual(obs, 'peanut')
        self.assertIsInstance(obs, str)

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


class TestParsePrimitiveCollectionsSimple(unittest.TestCase):
    def test_list_of_int(self):
        obs = parse_primitive(List[Int], ('1', '2', '3'))
        self.assertEqual(obs, [1, 2, 3])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], int)

    def test_set_of_int(self):
        obs = parse_primitive(Set[Int], ('1', '2', '3'))
        self.assertEqual(obs, {1, 2, 3})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), int)

    def test_list_of_float(self):
        obs = parse_primitive(List[Float], ('1.0', '2.0', '3.0'))
        self.assertEqual(obs, [1.0, 2.0, 3.0])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], float)

    def test_set_of_float(self):
        obs = parse_primitive(Set[Float], ('1.0', '2.0', '3.0'))
        self.assertEqual(obs, {1.0, 2.0, 3.0})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), float)

    def test_list_of_bool(self):
        obs = parse_primitive(List[Bool], ('True', 'False', 'True'))
        self.assertEqual(obs, [True, False, True])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], bool)

    def test_set_of_bool(self):
        obs = parse_primitive(Set[Bool], ('True', 'False'))
        self.assertEqual(obs, {True, False})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), bool)

    def test_list_of_str(self):
        obs = parse_primitive(List[Str], ('peanut', 'the', 'dog'))
        self.assertEqual(obs, ['peanut', 'the', 'dog'])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], str)

    def test_set_of_str(self):
        obs = parse_primitive(Set[Str], ('peanut', 'the', 'dog'))
        self.assertEqual(obs, {'peanut', 'the', 'dog'})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), str)

    # The next two tests _aren't_ monomorphic, because unions of Int and Float
    # always yield a Float (List[Int] | List[Float] == List[Float]).
    def test_list_int_or_float_with_int_value(self):
        obs = parse_primitive(List[Int] | List[Float], ('1', '2', '3'))
        self.assertEqual(obs, [1.0, 2.0, 3.0])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], float)

    def test_list_int_or_float_with_float_value(self):
        obs = parse_primitive(List[Int] | List[Float], ('1.1', '2.2', '3.3'))
        self.assertEqual(obs, [1.1, 2.2, 3.3])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], float)


class TestParsePrimitiveCollectionsMonomorphic(unittest.TestCase):
    def test_list_int_or_bool_with_int_value(self):
        obs = parse_primitive(List[Int] | List[Bool], ('1', '2', '3'))
        self.assertEqual(obs, [1, 2, 3])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], int)

    def test_list_int_or_bool_with_bool_value(self):
        obs = parse_primitive(List[Int] | List[Bool],
                              ('True', 'False', 'True'))
        self.assertEqual(obs, [True, False, True])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], bool)

    def test_set_int_or_bool_with_int_value(self):
        obs = parse_primitive(Set[Int] | Set[Bool], ('1', '2', '3'))
        self.assertEqual(obs, {1, 2, 3})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), int)

    def test_set_int_or_bool_with_bool_value(self):
        obs = parse_primitive(Set[Int] | Set[Bool], ('True', 'False'))
        self.assertEqual(obs, {True, False})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), bool)

    def test_list_int_or_str_with_int_value(self):
        obs = parse_primitive(List[Int] | List[Str], ('1', '2', '3'))
        self.assertEqual(obs, [1, 2, 3])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], int)

    def test_list_int_or_str_with_str_value(self):
        obs = parse_primitive(List[Int] | List[Str], ('peanut', 'the', 'dog'))
        self.assertEqual(obs, ['peanut', 'the', 'dog'])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], str)

    def test_set_int_or_str_with_int_value(self):
        obs = parse_primitive(Set[Int] | Set[Str], ('1', '2', '3'))
        self.assertEqual(obs, {1, 2, 3})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), int)

    def test_set_int_or_str_with_str_value(self):
        obs = parse_primitive(Set[Int] | Set[Str], ('peanut', 'the', 'dog'))
        self.assertEqual(obs, {'peanut', 'the', 'dog'})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), str)

    def test_list_float_or_bool_with_float_value(self):
        obs = parse_primitive(List[Float] | List[Bool], ('1.1', '2.2', '3.3'))
        self.assertEqual(obs, [1.1, 2.2, 3.3])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], float)

    def test_list_float_or_bool_with_bool_value(self):
        obs = parse_primitive(List[Float] | List[Bool],
                              ('True', 'False', 'True'))
        self.assertEqual(obs, [True, False, True])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], bool)

    def test_set_float_or_bool_with_float_value(self):
        obs = parse_primitive(Set[Float] | Set[Bool], ('1.1', '2.2', '3.3'))
        self.assertEqual(obs, {1.1, 2.2, 3.3})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), float)

    def test_set_float_or_bool_with_bool_value(self):
        obs = parse_primitive(Set[Float] | Set[Bool],
                              ('True', 'False', 'True'))
        self.assertEqual(obs, {True, False})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), bool)

    def test_list_float_or_str_with_float_value(self):
        obs = parse_primitive(List[Float] | List[Str], ('1.1', '2.2', '3.3'))
        self.assertEqual(obs, [1.1, 2.2, 3.3])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], float)

    def test_list_float_or_str_with_str_value(self):
        obs = parse_primitive(List[Float] | List[Str],
                              ('peanut', 'the', 'dog'))
        self.assertEqual(obs, ['peanut', 'the', 'dog'])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], str)

    def test_set_float_or_str_with_float_value(self):
        obs = parse_primitive(Set[Float] | Set[Str], ('1.1', '2.2', '3.3'))
        self.assertEqual(obs, {1.1, 2.2, 3.3})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), float)

    def test_set_float_or_str_with_str_value(self):
        obs = parse_primitive(Set[Float] | Set[Str], ('peanut', 'the', 'dog'))
        self.assertEqual(obs, {'peanut', 'the', 'dog'})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), str)

    def test_list_bool_or_str_with_bool_value(self):
        obs = parse_primitive(List[Bool] | List[Str],
                              ('True', 'False', 'True'))
        self.assertEqual(obs, [True, False, True])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], bool)

    def test_list_bool_or_str_with_str_value(self):
        obs = parse_primitive(List[Bool] | List[Str], ('peanut', 'the', 'dog'))
        self.assertEqual(obs, ['peanut', 'the', 'dog'])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], str)

    def test_set_bool_or_str_with_bool_value(self):
        obs = parse_primitive(Set[Bool] | Set[Str],
                              ('True', 'False', 'True'))
        self.assertEqual(obs, {True, False})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), bool)

    def test_set_bool_or_str_with_str_value(self):
        obs = parse_primitive(Set[Bool] | Set[Str], ('peanut', 'the', 'dog'))
        self.assertEqual(obs, {'peanut', 'the', 'dog'})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), str)


class TestParsePrimitiveCollectionsXYZ(unittest.TestCase):
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
