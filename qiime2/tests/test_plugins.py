# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest
import pkg_resources

from qiime2.plugins import assert_no_nans_in_tables


class TestTableNaNs(unittest.TestCase):
    def get_data_path(self, filename):
        return pkg_resources.resource_filename('qiime2.sdk.tests',
                                               'data/%s' % filename)

    def test_table_does_not_have_nans(self):
        noNaN = self.get_data_path('no_nan.html')

        with open(noNaN) as fh:
            assert_no_nans_in_tables(fh)

    def test_table_has_nans(self):
        hasNaN = self.get_data_path('has_nan.html')

        with open(hasNaN) as fh:
            with self.assertRaises(AssertionError):
                assert_no_nans_in_tables(fh)


if __name__ == '__main__':
    unittest.main()
