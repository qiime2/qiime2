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
    def test_int_value(self):
        value = '42'

        obs = parse_primitive(Int, value)
        self.assertEqual(42, obs)
        self.assertEqual(int, type(obs))

        obs = parse_primitive(Float, value)
        self.assertEqual(42.0, obs)
        self.assertEqual(float, type(obs))

        obs = parse_primitive(Bool, value)
        self.assertEqual(True, obs)
        self.assertEqual(bool, type(obs))

        obs = parse_primitive(Str, value)
        self.assertEqual('42', obs)
        self.assertEqual(str, type(obs))

    def test_float_value(self):
        value = '42.0'

        with self.assertRaisesRegex(ValueError, 'walk the plank'):
            parse_primitive(Int, value)

        obs = parse_primitive(Float, value)
        self.assertEqual(42.0, obs)
        self.assertEqual(float, type(obs))

        obs = parse_primitive(Bool, value)
        self.assertEqual(True, obs)
        self.assertEqual(bool, type(obs))

        obs = parse_primitive(Str, value)
        self.assertEqual('42.0', obs)
        self.assertEqual(str, type(obs))

    def test_bool_value(self):
        value = 'True'

        with self.assertRaisesRegex(ValueError, 'walk the plank'):
            parse_primitive(Int, value)

        with self.assertRaisesRegex(ValueError, 'walk the plank'):
            parse_primitive(Float, value)

        obs = parse_primitive(Bool, value)
        self.assertEqual(True, obs)
        self.assertEqual(bool, type(obs))

        obs = parse_primitive(Str, value)
        self.assertEqual('True', obs)
        self.assertEqual(str, type(obs))

    def test_str_value(self):
        value = 'peanut'

        with self.assertRaisesRegex(ValueError, 'walk the plank'):
            parse_primitive(Int, value)

        with self.assertRaisesRegex(ValueError, 'walk the plank'):
            parse_primitive(Float, value)

        obs = parse_primitive(Bool, value)
        self.assertEqual(True, obs)
        self.assertEqual(bool, type(obs))

        obs = parse_primitive(Str, value)
        self.assertEqual('peanut', obs)
        self.assertEqual(str, type(obs))


class TestParsePrimitiveNonCollectionsSimpleUnions(unittest.TestCase):
    def test_int_or_float_value(self):
        int_value = '42'
        float_value = '42.0'

        # Int | Float == Float
        obs = parse_primitive(Int | Float, int_value)
        self.assertEqual(42.0, obs)
        self.assertEqual(float, type(obs))

        obs = parse_primitive(Int | Float, float_value)
        self.assertEqual(42.0, obs)
        self.assertEqual(float, type(obs))

    def test_int_or_bool_value(self):
        int_value = '42'
        float_value = '42.0'
        bool_value = 'True'

        obs = parse_primitive(Int | Bool, int_value)
        self.assertEqual(42, obs)
        self.assertEqual(int, type(obs))

        with self.assertRaisesRegex(ValueError, 'foo'):
            parse_primitive(Int | Bool, float_value)

        obs = parse_primitive(Int | Bool, bool_value)
        self.assertEqual(True, obs)
        self.assertEqual(bool, type(obs))

    def test_int_or_str_value(self):
        pass

    def test_float_or_bool_value(self):
        pass

    def test_float_or_str_value(self):
        pass

    def test_bool_or_str_value(self):
        pass


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
