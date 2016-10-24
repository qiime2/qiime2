# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import pandas as pd
import pandas.util.testing as pdt

import qiime
from qiime.plugins import dummy_plugin


class TestProvenanceIntegration(unittest.TestCase):
    def test_chain_with_metadata(self):
        df = pd.DataFrame({'a': ['1', '2', '3']}, index=['0', '1', '2'])

        a = qiime.Artifact.import_data('IntSequence1', [1, 2, 3])
        m = qiime.Metadata(df)
        mc = qiime.MetadataCategory(df['a'])

        b = dummy_plugin.actions.identity_with_metadata(a, m).out
        c = dummy_plugin.actions.identity_with_metadata_category(b, mc).out

        p_dir = c._archiver.provenance_dir

        new_m = qiime.Metadata.load(
            str(p_dir / 'artifacts' / str(b.uuid) / 'action' / 'metadata.tsv'))

        pdt.assert_frame_equal(m.to_dataframe(), new_m.to_dataframe(),
                               check_names=False)

        with (p_dir / 'action' / 'metadata.tsv').open() as fh:
            self.assertEqual(fh.read(), '0\t1\n1\t2\n2\t3\n')


if __name__ == '__main__':
    unittest.main()
