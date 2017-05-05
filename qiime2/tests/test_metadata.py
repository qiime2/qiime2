# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources
import sqlite3
import unittest

import pandas as pd
import pandas.util.testing as pdt

import qiime2
from qiime2.core.testing.util import get_dummy_plugin


class TestMetadataLoad(unittest.TestCase):
    def test_does_not_cast_index_or_column_types(self):
        fp = pkg_resources.resource_filename(
            'qiime2.tests', 'data/metadata-no-type-cast.tsv')

        metadata = qiime2.Metadata.load(fp)
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


class TestMetadataFromArtifact(unittest.TestCase):
    def setUp(self):
        get_dummy_plugin()

    def test_from_artifact(self):
        A = qiime2.Artifact.import_data('Mapping', {'a': '1', 'b': '3'})
        md = qiime2.Metadata.from_artifact(A)
        pdt.assert_frame_equal(md.to_dataframe(),
                               pd.DataFrame({'a': '1', 'b': '3'}, index=[0]))

    def test_from_bad_artifact(self):
        A = qiime2.Artifact.import_data('IntSequence1', [1, 2, 3, 4])
        with self.assertRaisesRegex(ValueError, 'metadata'):
            qiime2.Metadata.from_artifact(A)

    def test_artifact(self):
        A = qiime2.Artifact.import_data('Mapping', {'a': ['1', '2'],
                                                    'b': ['2', '3']})
        md = qiime2.Metadata.from_artifact(A)
        art = md.artifact
        self.assertIsInstance(art, qiime2.Artifact)
        self.assertEquals(art, A)


class TestIDs(unittest.TestCase):
    def test_default(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)

        actual = metadata.ids()
        expected = {'S1', 'S2', 'S3'}
        self.assertEqual(actual, expected)

    def test_incomplete_where(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=['S1', 'S2', 'S3'])
        metadata = qiime2.Metadata(df)

        where = "Subject='subject-1' AND SampleType="
        with self.assertRaises(ValueError):
            metadata.ids(where)

        where = "Subject="
        with self.assertRaises(ValueError):
            metadata.ids(where)

    def test_invalid_where(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=['S1', 'S2', 'S3'])
        metadata = qiime2.Metadata(df)

        where = "not-a-column-name='subject-1'"
        with self.assertRaises(ValueError):
            metadata.ids(where)

    def test_empty_result(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)

        where = "Subject='subject-3'"
        actual = metadata.ids(where)
        expected = set()
        self.assertEqual(actual, expected)

    def test_simple_expression(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)

        where = "Subject='subject-1'"
        actual = metadata.ids(where)
        expected = {'S1', 'S2'}
        self.assertEqual(actual, expected)

        where = "Subject='subject-2'"
        actual = metadata.ids(where)
        expected = {'S3'}
        self.assertEqual(actual, expected)

        where = "Subject='subject-3'"
        actual = metadata.ids(where)
        expected = set()
        self.assertEqual(actual, expected)

        where = "SampleType='gut'"
        actual = metadata.ids(where)
        expected = {'S1', 'S3'}
        self.assertEqual(actual, expected)

        where = "SampleType='tongue'"
        actual = metadata.ids(where)
        expected = {'S2'}
        self.assertEqual(actual, expected)

    def test_more_complex_expressions(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)

        where = "Subject='subject-1' OR Subject='subject-2'"
        actual = metadata.ids(where)
        expected = {'S1', 'S2', 'S3'}
        self.assertEqual(actual, expected)

        where = "Subject='subject-1' AND Subject='subject-2'"
        actual = metadata.ids(where)
        expected = set()
        self.assertEqual(actual, expected)

        where = "Subject='subject-1' AND SampleType='gut'"
        actual = metadata.ids(where)
        expected = {'S1'}
        self.assertEqual(actual, expected)

    def test_index_without_name(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=['S1', 'S2', 'S3'])
        metadata = qiime2.Metadata(df)

        actual = metadata.ids(where="SampleType='gut'")
        expected = {'S1', 'S3'}
        self.assertEqual(actual, expected)

    def test_index_with_column_name_clash(self):
        df = pd.DataFrame(
            {'Subject': ['subject-1', 'subject-1', 'subject-2'],
             'SampleType': ['gut', 'tongue', 'gut']},
            index=pd.Index(['S1', 'S2', 'S3'], name='SampleType'))
        metadata = qiime2.Metadata(df)

        with self.assertRaises(sqlite3.OperationalError):
            metadata.ids(where="Subject='subject-1'")

    def test_duplicate_columns(self):
        df = pd.DataFrame([['subject-1', 'gut'],
                           ['subject-1', 'tongue'],
                           ['subject-2', 'gut']],
                          index=['S1', 'S2', 'S3'], columns=['foo', 'foo'])
        metadata = qiime2.Metadata(df)

        with self.assertRaises(sqlite3.OperationalError):
            metadata.ids(where="foo='subject-2'")

    def test_query_by_index(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)

        actual = metadata.ids(where="id='S2' OR id='S1'")
        expected = {'S1', 'S2'}
        self.assertEqual(actual, expected)


if __name__ == '__main__':
    unittest.main()
