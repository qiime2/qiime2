# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest


class TestParseParameterNonCollections(unittest.TestCase):
    def test_int(self):
        pass

    def test_float(self):
        pass

    def test_bool(self):
        pass

    def test_str(self):
        pass

    def test_int_or_float(self):
        pass

    def test_int_or_bool(self):
        pass

    def test_int_or_str(self):
        pass

    def test_float_or_bool(self):
        pass

    def test_float_or_str(self):
        pass

    def test_bool_or_str(self):
        pass


class TestParseParameterCollections(unittest.TestCase):
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
