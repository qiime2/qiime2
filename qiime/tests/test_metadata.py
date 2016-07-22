# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources
import unittest

import pandas as pd
import pandas.util.testing as pdt

import qiime


class TestMetadata(unittest.TestCase):
    def test_load_does_not_cast_index_or_column_types(self):
        fp = pkg_resources.resource_filename(
            'qiime.tests', 'data/metadata-no-type-cast.tsv')

        metadata = qiime.Metadata.load(fp)
        df = metadata.to_dataframe()

        exp_index = pd.Index(['0.000001', '0.004000', '0.000000'],
                             dtype=object, name='my-index')
        exp_df = pd.DataFrame({'col1': ['2', '1', '3'],
                               'col2': ['b', 'b', 'c'],
                               'col3': ['2.5', '4.2', '-9.999']},
                              index=exp_index, dtype=object)
        pdt.assert_frame_equal(
            df, exp_df, check_dtype=True, check_index_type=True,
            check_column_type=True, check_frame_type=True, check_names=True,
            check_exact=True)


if __name__ == '__main__':
    unittest.main()
