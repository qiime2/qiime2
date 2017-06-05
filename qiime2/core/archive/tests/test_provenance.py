# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import pandas as pd
import pandas.util.testing as pdt

import qiime2
from qiime2.plugins import dummy_plugin


class TestProvenanceIntegration(unittest.TestCase):
    def test_chain_with_metadata(self):
        df = pd.DataFrame({'a': ['1', '2', '3']}, index=['0', '1', '2'])

        a = qiime2.Artifact.import_data('IntSequence1', [1, 2, 3])
        m = qiime2.Metadata(df)
        mc = qiime2.MetadataCategory(df['a'])

        b = dummy_plugin.actions.identity_with_metadata(a, m).out
        c = dummy_plugin.actions.identity_with_metadata_category(b, mc).out

        p_dir = c._archiver.provenance_dir

        new_m = qiime2.Metadata.load(
            str(p_dir / 'artifacts' / str(b.uuid) / 'action' / 'metadata.tsv'))

        pdt.assert_frame_equal(m.to_dataframe(), new_m.to_dataframe(),
                               check_names=False)

        with (p_dir / 'action' / 'metadata.tsv').open() as fh:
            self.assertEqual(fh.read(), '0\t1\n1\t2\n2\t3\n')

    def test_chain_with_artifact_metadata(self):
        metadata_artifact_1 = qiime2.Artifact.import_data(
            'Mapping', {'a': 'foo', 'b': 'bar'})
        metadata_artifact_2 = qiime2.Artifact.import_data(
            'Mapping', {'c': 'baz'})
        m = qiime2.Metadata.from_artifact(metadata_artifact_1)
        mc = qiime2.MetadataCategory.from_artifact(metadata_artifact_2, 'c')

        a = qiime2.Artifact.import_data('IntSequence1', [1, 2, 3])

        b = dummy_plugin.actions.identity_with_metadata(a, m).out
        c = dummy_plugin.actions.identity_with_metadata_category(b, mc).out

        p_dir = c._archiver.provenance_dir

        m_yaml_value = "%s:metadata.tsv" % metadata_artifact_1.uuid
        mc_yaml_value = "%s:metadata.tsv" % metadata_artifact_2.uuid

        # Check action files for uuid-metadata values
        with (p_dir / 'action' / 'action.yaml').open() as fh:
            self.assertIn(mc_yaml_value, fh.read())
        with (p_dir / 'artifacts' / str(b.uuid) / 'action' /
              'action.yaml').open() as fh:
            self.assertIn(m_yaml_value, fh.read())

        # Check that metadata is written out fully
        new_m = qiime2.Metadata.load(
            str(p_dir / 'artifacts' / str(b.uuid) / 'action' / 'metadata.tsv'))

        pdt.assert_frame_equal(m.to_dataframe(), new_m.to_dataframe(),
                               check_names=False)

        # Check that provenance of originating metadata artifact exists
        self.assertTrue((p_dir / 'artifacts' / str(metadata_artifact_1.uuid) /
                         'action' / 'action.yaml').exists())
        self.assertTrue((p_dir / 'artifacts' / str(metadata_artifact_2.uuid) /
                         'action' / 'action.yaml').exists())

    def test_chain_with_merged_artifact_metadata(self):
        md_artifact1 = qiime2.Artifact.import_data(
            'Mapping', {'a': 'foo', 'b': 'bar'})
        md_artifact2 = qiime2.Artifact.import_data(
            'Mapping', {'c': 'baz'})
        md1 = qiime2.Metadata.from_artifact(md_artifact1)
        md2 = qiime2.Metadata.from_artifact(md_artifact2)
        merged_md = md1.merge(md2)
        merged_mdc = merged_md.get_category('c')

        a = qiime2.Artifact.import_data('IntSequence1', [1, 2, 3])

        b = dummy_plugin.actions.identity_with_metadata(a, merged_md).out
        c = dummy_plugin.actions.identity_with_metadata_category(
            b, merged_mdc).out

        p_dir = c._archiver.provenance_dir

        yaml_value = "%s,%s:metadata.tsv" % (md_artifact1.uuid,
                                             md_artifact2.uuid)

        # Check action files for uuid-metadata values
        with (p_dir / 'action' / 'action.yaml').open() as fh:
            self.assertIn(yaml_value, fh.read())
        with (p_dir / 'artifacts' / str(b.uuid) / 'action' /
              'action.yaml').open() as fh:
            self.assertIn(yaml_value, fh.read())

        # Check that metadata is written out fully
        with (p_dir / 'action' / 'metadata.tsv').open() as fh:
            self.assertEqual(fh.read(), '0\tbaz\n')

        new_merged_md = qiime2.Metadata.load(
            str(p_dir / 'artifacts' / str(b.uuid) / 'action' / 'metadata.tsv'))
        pdt.assert_frame_equal(new_merged_md.to_dataframe(),
                               merged_md.to_dataframe(), check_names=False)

        # Check that provenance of originating metadata artifacts exists
        self.assertTrue((p_dir / 'artifacts' / str(md_artifact1.uuid) /
                         'action' / 'action.yaml').exists())
        self.assertTrue((p_dir / 'artifacts' / str(md_artifact2.uuid) /
                         'action' / 'action.yaml').exists())


if __name__ == '__main__':
    unittest.main()
