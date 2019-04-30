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
        with self.assertRaisesRegex(ValueError, 'Could not coerce'):
            parse_primitive(Bool, '42')

    def test_str_type_int_value(self):
        obs = parse_primitive(Str, '42')
        self.assertEqual(obs, '42')
        self.assertIsInstance(obs, str)

    def test_int_type_float_value(self):
        with self.assertRaisesRegex(ValueError, 'Could not coerce'):
            parse_primitive(Int, '42.0')

    def test_float_type_float_value(self):
        obs = parse_primitive(Float, '42.0')
        self.assertEqual(obs, 42.0)
        self.assertIsInstance(obs, float)

    def test_bool_type_float_value(self):
        with self.assertRaisesRegex(ValueError, 'Could not coerce'):
            parse_primitive(Bool, '42.0')

    def test_str_type_float_value(self):
        obs = parse_primitive(Str, '42.0')
        self.assertEqual(obs, '42.0')
        self.assertIsInstance(obs, str)

    def test_int_type_bool_value(self):
        with self.assertRaisesRegex(ValueError, 'Could not coerce'):
            parse_primitive(Int, 'True')

    def test_float_type_bool_value(self):
        with self.assertRaisesRegex(ValueError, 'Could not coerce'):
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
        with self.assertRaisesRegex(ValueError, 'Could not coerce'):
            parse_primitive(Int, 'peanut')

    def test_float_type_str_value(self):
        with self.assertRaisesRegex(ValueError, 'Could not coerce'):
            parse_primitive(Float, 'peanut')

    def test_bool_type_str_value(self):
        with self.assertRaisesRegex(ValueError, 'Could not coerce'):
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
        # Int | Float == Float
        with self.assertRaisesRegex(ValueError, 'Could not coerce'):
            parse_primitive(Int | Float, 'True')

    def test_int_union_float_expr_str_value(self):
        # Int | Float == Float
        with self.assertRaisesRegex(ValueError, 'Could not coerce'):
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

    # The next tests _aren't_ monomorphic, because unions of Int and Float
    # always yield a Float (List[Int] | List[Float] == List[Float]).
    def test_list_int_or_float_with_int_value(self):
        obs = parse_primitive(List[Int] | List[Float], ('1', '2', '3'))
        self.assertEqual(obs, [1.0, 2.0, 3.0])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], float)

    def test_set_int_or_float_with_int_value(self):
        obs = parse_primitive(Set[Int] | Set[Float], ('1', '2', '3'))
        self.assertEqual(obs, {1.0, 2.0, 3.0})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), float)

    def test_list_int_or_float_with_float_value(self):
        obs = parse_primitive(List[Int] | List[Float], ('1.1', '2.2', '3.3'))
        self.assertEqual(obs, [1.1, 2.2, 3.3])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], float)

    def test_set_int_or_float_with_float_value(self):
        obs = parse_primitive(Set[Int] | Set[Float], ('1.1', '2.2', '3.3'))
        self.assertEqual(obs, {1.1, 2.2, 3.3})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), float)

    def test_list_int_or_float_int_value(self):
        obs = parse_primitive(List[Int | Float], ('1', '2', '3'))
        self.assertEqual(obs, [1.0, 2.0, 3.0])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], float)

    def test_set_int_or_float_int_value(self):
        obs = parse_primitive(Set[Int | Float], ('1', '2', '3'))
        self.assertEqual(obs, {1.0, 2.0, 3.0})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), float)


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

    def test_list_int_or_bool_with_mixed_value_variant_a(self):
        pass

    def test_list_int_or_bool_with_mixed_value_variant_b(self):
        pass

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

    def test_list_int_or_str_with_mixed_value_variant_a(self):
        obs = parse_primitive(List[Int] | List[Str], ('1', 'the', 'dog'))
        self.assertEqual(obs, ['1', 'the', 'dog'])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], str)
        self.assertIsInstance(obs[1], str)

    def test_list_int_or_str_with_mixed_value_variant_b(self):
        obs = parse_primitive(List[Int] | List[Str], ('peanut', 'the', '1'))
        self.assertEqual(obs, ['peanut', 'the', '1'])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], str)
        self.assertIsInstance(obs[2], str)

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

    def test_list_float_or_bool_with_mixed_value_variant_a(self):
        pass

    def test_list_float_or_bool_with_mixed_value_variant_b(self):
        pass

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

    def test_list_float_or_str_with_mixed_value_variant_a(self):
        pass

    def test_list_float_or_str_with_mixed_value_variant_b(self):
        pass

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

    def test_list_bool_or_str_with_mixed_value_variant_a(self):
        pass

    def test_list_bool_or_str_with_mixed_value_variant_b(self):
        pass

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

    # TODO: test this branch
    # def test_list_bool_or_str_with_mixed_value(self):
    #     with self.assertRaisesRegex(ValueError, 'not all matching'):
    #         parse_primitive(List[Bool] | List[Str], ('peanut', 'the', 'True'))


class TestParsePrimitiveCollectionsComposite(unittest.TestCase):
    def test_list_int_or_bool_with_int_value(self):
        obs = parse_primitive(List[Int | Bool], ('1', '2', '3'))
        self.assertEqual(obs, [1, 2, 3])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], int)

    def test_list_int_or_bool_with_float_value(self):
        with self.assertRaisesRegex(ValueError, 'Could not coerce'):
            parse_primitive(List[Int | Bool], ('1.1', '2.2', '3.3'))

    def test_list_int_or_bool_with_bool_value(self):
        obs = parse_primitive(List[Int | Bool], ('True', 'False', 'True'))
        self.assertEqual(obs, [True, False, True])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], bool)

    def test_list_int_or_bool_with_str_value(self):
        with self.assertRaisesRegex(ValueError, 'Could not coerce'):
            parse_primitive(List[Int | Bool], ('peanut', 'the', 'dog'))

    def test_list_int_or_bool_with_mixed_value(self):
        obs = parse_primitive(List[Int | Bool], ('1', 'False', 2, 'True'))
        self.assertEqual(obs, [1, False, 2, True])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], int)
        self.assertIsInstance(obs[1], bool)

    def test_list_int_or_bool_with_mixed_value_variant_a(self):
        pass

    def test_list_int_or_bool_with_mixed_value_variant_b(self):
        pass

    def test_list_int_or_bool_with_bad_mix_value(self):
        with self.assertRaisesRegex(ValueError, 'Could not coerce'):
            parse_primitive(List[Int | Bool], ('1', 'True', 'dog'))

    def test_set_int_or_bool_with_int_value(self):
        obs = parse_primitive(Set[Int | Bool], ('1', '2', '3'))
        self.assertEqual(obs, {1, 2, 3})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), int)

    def test_set_int_or_bool_with_bool_value(self):
        obs = parse_primitive(Set[Int | Bool], ('True', 'False', 'True'))
        self.assertEqual(obs, {True, False})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), bool)

    def test_set_int_or_bool_with_mixed_value(self):
        obs = parse_primitive(Set[Int | Bool], ('1', 'False', '2', 'True'))
        self.assertEqual(obs, {1, False, 2, True})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), bool)
        self.assertIsInstance(obs.pop(), int)

    def test_list_int_or_str_with_int_value(self):
        obs = parse_primitive(List[Int | Str], ('1', '2', '3'))
        self.assertEqual(obs, [1, 2, 3])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], int)

    def test_list_int_or_str_with_str_value(self):
        obs = parse_primitive(List[Int | Str], ('peanut', 'the', 'dog'))
        self.assertEqual(obs, ['peanut', 'the', 'dog'])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], str)

    def test_list_int_or_str_with_mixed_value_variant_a(self):
        obs = parse_primitive(List[Int | Str], ('1', 'the', 'dog'))
        self.assertEqual(obs, [1, 'the', 'dog'])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], int)
        self.assertIsInstance(obs[1], str)

    def test_list_int_or_str_with_mixed_value_variant_b(self):
        obs = parse_primitive(List[Int | Str], ('peanut', 'the', '1'))
        self.assertEqual(obs, ['peanut', 'the', 1])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], str)
        self.assertIsInstance(obs[2], int)

    def test_set_int_or_str_with_int_value(self):
        obs = parse_primitive(Set[Int | Str], ('1', '2', '3'))
        self.assertEqual(obs, {1, 2, 3})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), int)

    def test_set_int_or_str_with_str_value(self):
        obs = parse_primitive(Set[Int | Str], ('peanut', 'the', 'dog'))
        self.assertEqual(obs, {'peanut', 'the', 'dog'})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), str)

    def test_set_int_or_str_with_mixed_value(self):
        obs = parse_primitive(Set[Int | Str], ('1', 'the', '2', 'dog'))
        self.assertEqual(obs, {1, 'the', 2, 'dog'})
        self.assertIsInstance(obs, set)

    def test_list_float_or_bool_with_float_value(self):
        obs = parse_primitive(List[Float | Bool], ('1.1', '2.2', '3.3'))
        self.assertEqual(obs, [1.1, 2.2, 3.3])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], float)

    def test_list_float_or_bool_with_bool_value(self):
        obs = parse_primitive(List[Float | Bool], ('True', 'False', 'True'))
        self.assertEqual(obs, [True, False, True])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], bool)

    def test_list_float_or_bool_with_mixed_value_variant_a(self):
        obs = parse_primitive(List[Float | Bool], ('True', '2.2', '3.3'))
        self.assertEqual(obs, [True, 2.2, 3.3])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], bool)
        self.assertIsInstance(obs[1], float)

    def test_list_float_or_bool_with_mixed_value_variant_b(self):
        obs = parse_primitive(List[Float | Bool], ('1.1', '2.2', 'False'))
        self.assertEqual(obs, [1.1, 2.2, False])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], float)
        self.assertIsInstance(obs[-1], bool)

    def test_list_float_or_bool_with_bad_mix_value(self):
        with self.assertRaisesRegex(ValueError, 'Could not coerce'):
            parse_primitive(List[Float | Bool], ('1.1', '2.2', 'peanut'))

    def test_set_float_or_bool_with_float_value(self):
        obs = parse_primitive(Set[Float | Bool], ('1.1', '2.2', '3.3'))
        self.assertEqual(obs, {1.1, 2.2, 3.3})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), float)

    def test_set_float_or_bool_with_bool_value(self):
        obs = parse_primitive(Set[Float | Bool], ('True', 'False', 'True'))
        self.assertEqual(obs, {True, False})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), bool)

    def test_list_float_or_str_with_float_value(self):
        obs = parse_primitive(List[Float | Str], ('1.1', '2.2', '3.3'))
        self.assertEqual(obs, [1.1, 2.2, 3.3])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], float)

    def test_list_float_or_str_with_str_value(self):
        obs = parse_primitive(List[Float | Str], ('peanut', 'the', 'dog'))
        self.assertEqual(obs, ['peanut', 'the', 'dog'])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], str)

    def test_list_float_or_str_with_mixed_value_variant_a(self):
        obs = parse_primitive(List[Float | Str], ('peanut', '2.2', '3.3'))
        self.assertEqual(obs, ['peanut', 2.2, 3.3])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], str)
        self.assertIsInstance(obs[1], float)

    def test_list_float_or_str_with_mixed_value_variant_b(self):
        obs = parse_primitive(List[Float | Str], ('1.1', '2.2', 'dog'))
        self.assertEqual(obs, [1.1, 2.2, 'dog'])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], float)
        self.assertIsInstance(obs[-1], str)

    def test_set_float_or_str_with_float_value(self):
        obs = parse_primitive(Set[Float | Str], ('1.1', '2.2', '3.3'))
        self.assertEqual(obs, {1.1, 2.2, 3.3})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), float)

    def test_set_float_or_str_with_str_value(self):
        obs = parse_primitive(Set[Float | Str], ('peanut', 'the', 'dog'))
        self.assertEqual(obs, {'peanut', 'the', 'dog'})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), str)

    def test_list_bool_or_str_with_bool_value(self):
        obs = parse_primitive(List[Bool | Str], ('True', 'False', 'True'))
        self.assertEqual(obs, [True, False, True])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], bool)

    def test_list_bool_or_str_with_str_value(self):
        obs = parse_primitive(List[Bool | Str], ('peanut', 'the', 'dog'))
        self.assertEqual(obs, ['peanut', 'the', 'dog'])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], str)

    def test_list_bool_or_str_with_mixed_value_variant_a(self):
        obs = parse_primitive(List[Bool | Str], ('True', 'the', 'dog'))
        self.assertEqual(obs, [True, 'the', 'dog'])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], bool)
        self.assertIsInstance(obs[-1], str)

    def test_list_bool_or_str_with_mixed_value_variant_b(self):
        obs = parse_primitive(List[Bool | Str], ('peanut', 'the', 'True'))
        self.assertEqual(obs, ['peanut', 'the', True])
        self.assertIsInstance(obs, list)
        self.assertIsInstance(obs[0], str)
        self.assertIsInstance(obs[-1], bool)

    def test_set_bool_or_str_with_bool_value(self):
        obs = parse_primitive(Set[Bool | Str], ('True', 'False', 'True'))
        self.assertEqual(obs, {True, False})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), bool)

    def test_set_bool_or_str_with_str_value(self):
        obs = parse_primitive(Set[Bool | Str], ('peanut', 'the', 'dog'))
        self.assertEqual(obs, {'peanut', 'the', 'dog'})
        self.assertIsInstance(obs, set)
        self.assertIsInstance(obs.pop(), str)


class TestParsePrimitiveCollectionsComplex(unittest.TestCase):
    def test_list_int_bool_or_list_float_with_bool_int_value(self):
        obs = parse_primitive(List[Int | Bool] | List[Float],
                              ('1', '2', 'True', 'False'))
        self.assertEqual(obs, [1, 2, True, False])

    def test_list_int_bool_or_list_float_with_float_value(self):
        obs = parse_primitive(List[Int | Bool] | List[Float],
                              ('1.1', '2.2', '3.3', '4.4'))
        self.assertEqual(obs, [1.1, 2.2, 3.3, 4.4])

    def test_list_int_bool_or_list_float_with_bad_value(self):
        with self.assertRaisesRegex(ValueError, 'Could not coerce'):
            parse_primitive(List[Int | Bool] | List[Float],
                            ('1', '2.2', 'True', 'False'))

    def test_list_int_str_or_list_float_with_str_int_value(self):
        obs = parse_primitive(List[Int | Str] | List[Float],
                              ('1', '2', 'peanut', 'the'))
        self.assertEqual(obs, [1, 2, 'peanut', 'the'])

    def test_list_int_str_or_list_float_with_float_value(self):
        obs = parse_primitive(List[Int | Str] | List[Float],
                              ('1.1', '2.2', '3.3', '4.4'))
        self.assertEqual(obs, [1.1, 2.2, 3.3, 4.4])

    def test_list_int_str_or_list_float_with_mixed_value(self):
        obs = parse_primitive(List[Int | Str] | List[Float],
                              ('1.1', '2', 'True', 'peanut'))
        self.assertEqual(obs, ['1.1', 2, 'True', 'peanut'])

    def test_list_float_bool_or_list_str_with_float_bool_value(self):
        obs = parse_primitive(List[Float | Bool] | List[Int],
                              ('1', '2', 'True', 'False'))
        self.assertEqual(obs, [1, 2, True, False])

    def test_list_float_bool_or_list_str_with_int_value(self):
        obs = parse_primitive(List[Float | Bool] | List[Int],
                              ('1', '2', '3', '4'))
        self.assertEqual(obs, [1, 2, 3, 4])

    def test_list_float_bool_or_list_str_with_bad_value(self):
        with self.assertRaisesRegex(ValueError, 'Could not coerce'):
            parse_primitive(List[Float | Bool] | List[Int],
                            ('1', '2.2', 'True', 'peanut'))

    def test_set_int_bool_or_list_float_with_bool_int_value(self):
        obs = parse_primitive(Set[Int | Bool] | Set[Float],
                              ('1', '2', 'True', 'False'))
        self.assertEqual(obs, {1, 2, True, False})


class TestParsePrimitiveCollectionsNonStringInputs(unittest.TestCase):
    pass


if __name__ == '__main__':
    unittest.main()
