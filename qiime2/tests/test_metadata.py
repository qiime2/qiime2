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
import tempfile

import pandas as pd
import pandas.util.testing as pdt

import qiime2
from qiime2.core.testing.util import get_dummy_plugin, ReallyEqualMixin


class TestMetadata(unittest.TestCase):
    def setUp(self):
        self.illegal_chars = ['/', '\0', '\\', '*', '<', '>', '?', '|', '$']

    def test_valid_metadata(self):
        exp_index = pd.Index(['a', 'b', 'c'], dtype=object)
        exp_df = pd.DataFrame({'col1': ['2', '1', '3']},
                              index=exp_index, dtype=object)

        metadata = qiime2.Metadata(exp_df)
        df = metadata.to_dataframe()

        pdt.assert_frame_equal(
            df, exp_df, check_dtype=True, check_index_type=True,
            check_column_type=True, check_frame_type=True, check_names=True,
            check_exact=True)

    def test_valid_metadata_str(self):
        exp_index = pd.Index(['a', 'b', 'c'], dtype=str)
        exp_df = pd.DataFrame({'col1': ['2', '1', '3']},
                              index=exp_index, dtype=str)

        metadata = qiime2.Metadata(exp_df)
        df = metadata.to_dataframe()

        pdt.assert_frame_equal(
            df, exp_df, check_dtype=True, check_index_type=True,
            check_column_type=True, check_frame_type=True, check_names=True,
            check_exact=True)

    def test_valid_metadata_no_columns(self):
        exp_index = pd.Index(['a', 'b', 'c'], dtype=object)
        exp_df = pd.DataFrame({}, index=exp_index, dtype=object)

        metadata = qiime2.Metadata(exp_df)
        obs_df = metadata.to_dataframe()

        self.assertFalse(obs_df.index.empty)
        self.assertTrue(obs_df.columns.empty)
        pdt.assert_frame_equal(
            obs_df, exp_df, check_dtype=True, check_index_type=True,
            check_column_type=True, check_frame_type=True, check_names=True,
            check_exact=True)

    def test_artifacts(self):
        index = pd.Index(['a', 'b', 'c'], dtype=object)
        df = pd.DataFrame({'col1': ['2', '1', '3']}, index=index, dtype=object)

        metadata = qiime2.Metadata(df)

        self.assertEqual(metadata.artifacts, [])

    def test_empty_metadata(self):
        # No index, no columns.
        df = pd.DataFrame([], index=[])

        with self.assertRaisesRegex(ValueError, 'Metadata is empty'):
            qiime2.Metadata(df)

        # No index, has columns.
        df = pd.DataFrame([], index=[], columns=['a', 'b'])

        with self.assertRaisesRegex(ValueError, 'Metadata is empty'):
            qiime2.Metadata(df)

    def test_invalid_metadata_characters_in_category(self):
        for val in self.illegal_chars:
            index = pd.Index(['a', 'b', 'c'], dtype=object)
            df = pd.DataFrame({'col1%s' % val: ['2', '1', '3']},
                              index=index, dtype=object)

            with self.assertRaisesRegex(ValueError,
                                        'Invalid characters.*category'):
                qiime2.Metadata(df)

    def test_invalid_metadata_characters_in_index(self):
        for val in self.illegal_chars:
            index = pd.Index(['a', 'b%s' % val, 'c'], dtype=object)
            df = pd.DataFrame({'col1': ['2', '1', '3']},
                              index=index, dtype=object)

            with self.assertRaisesRegex(ValueError,
                                        'Invalid character.*index'):
                qiime2.Metadata(df)

    def test_invalid_columns_dtype(self):
        with self.assertRaisesRegex(ValueError, 'Non-string.*category label'):
            qiime2.Metadata(pd.DataFrame(['a', 'b', 'c']))

    def test_invalid_index_dtype(self):
        with self.assertRaisesRegex(ValueError, 'Non-string.*index values'):
            qiime2.Metadata(pd.DataFrame({'foo': ['a', 'b', 'c']}))

    def test_duplicate_categories(self):
        index = pd.Index(['a', 'b'], dtype=object)
        df = pd.DataFrame({'foo': [1, 2], 'bar': [3, 4]}, index=index)
        df.columns = ['foo', 'foo']

        with self.assertRaisesRegex(ValueError, 'Duplicate.*category'):
            qiime2.Metadata(df)

    def test_duplicate_indices(self):
        index = pd.Index(['b', 'b', 'b'], dtype=object)
        df = pd.DataFrame({'foo': [1, 2, 3]}, index=index)

        with self.assertRaisesRegex(ValueError, 'Duplicate.*index values'):
            qiime2.Metadata(df)


class TestMetadataLoad(unittest.TestCase):
    def test_no_columns(self):
        fp = pkg_resources.resource_filename(
            'qiime2.tests', 'data/metadata-no-columns.tsv')

        metadata = qiime2.Metadata.load(fp)
        obs_df = metadata.to_dataframe()

        exp_index = pd.Index(['a', 'b', 'id'], name='my-index', dtype=object)
        exp_df = pd.DataFrame({}, index=exp_index, dtype=object)

        self.assertFalse(obs_df.index.empty)
        self.assertTrue(obs_df.columns.empty)
        pdt.assert_frame_equal(
            obs_df, exp_df, check_dtype=True, check_index_type=True,
            check_column_type=True, check_frame_type=True, check_names=True,
            check_exact=True)

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

    def test_artifacts(self):
        fp = pkg_resources.resource_filename(
            'qiime2.tests', 'data/metadata.tsv')

        metadata = qiime2.Metadata.load(fp)

        self.assertEqual(metadata.artifacts, [])

    def test_invalid_metadata_characters_in_category(self):
        fp = pkg_resources.resource_filename(
            'qiime2.tests', 'data/metadata-illegal-categories-characters.tsv')

        with self.assertRaisesRegex(ValueError,
                                    'Invalid characters.*category'):
            qiime2.Metadata.load(fp)

    def test_invalid_metadata_characters_in_index(self):
        fp = pkg_resources.resource_filename(
            'qiime2.tests', 'data/metadata-illegal-index-characters.tsv')

        with self.assertRaisesRegex(ValueError,
                                    'Invalid characters.*index'):
            qiime2.Metadata.load(fp)

    def test_non_tsv_metadata_file(self):
        with tempfile.TemporaryFile() as bad_file:
            bad_file.write(b'\x07\x08\x07')
            with self.assertRaisesRegex(ValueError, 'No columns to parse'):
                qiime2.Metadata.load(bad_file)


class TestMetadataFromArtifact(unittest.TestCase):
    def setUp(self):
        get_dummy_plugin()

    def test_from_artifact(self):
        A = qiime2.Artifact.import_data('Mapping', {'a': '1', 'b': '3'})
        md = qiime2.Metadata.from_artifact(A)
        pdt.assert_frame_equal(md.to_dataframe(),
                               pd.DataFrame({'a': '1', 'b': '3'}, index=['0']))

    def test_from_bad_artifact(self):
        A = qiime2.Artifact.import_data('IntSequence1', [1, 2, 3, 4])
        with self.assertRaisesRegex(ValueError, 'Artifact has no metadata'):
            qiime2.Metadata.from_artifact(A)

    def test_invalid_metadata_characters_in_category(self):
        A = qiime2.Artifact.import_data('Mapping', {'a': '1', '>b': '3'})
        with self.assertRaisesRegex(ValueError, 'Invalid characters'):
            qiime2.Metadata.from_artifact(A)

    def test_artifacts(self):
        A = qiime2.Artifact.import_data('Mapping', {'a': ['1', '2'],
                                                    'b': ['2', '3']})
        md = qiime2.Metadata.from_artifact(A)
        obs = md.artifacts
        self.assertEqual(obs, [A])


class TestGetCategory(unittest.TestCase):
    def setUp(self):
        get_dummy_plugin()

    def test_artifacts_are_propagated(self):
        A = qiime2.Artifact.import_data('Mapping', {'a': '1', 'b': '3'})
        md = qiime2.Metadata.from_artifact(A)

        obs = md.get_category('b')

        self.assertEqual(obs.artifacts, [A])
        pdt.assert_series_equal(obs.to_series(),
                                pd.Series(['3'], index=['0'], name='b'))


class TestMerge(unittest.TestCase):
    def setUp(self):
        get_dummy_plugin()

    def test_merging_one(self):
        md = qiime2.Metadata(pd.DataFrame(
            {'a': [1, 2, 3], 'b': [4, 5, 6]}, index=['id1', 'id2', 'id3']))

        obs = md.merge()

        self.assertIsNot(obs, md)
        self.assertEqual(obs, md)

    def test_merging_two(self):
        md1 = qiime2.Metadata(pd.DataFrame(
            {'a': [1, 2, 3], 'b': [4, 5, 6]}, index=['id1', 'id2', 'id3']))
        md2 = qiime2.Metadata(pd.DataFrame(
            {'c': [7, 8, 9], 'd': [10, 11, 12]}, index=['id1', 'id2', 'id3']))

        obs = md1.merge(md2)

        exp = qiime2.Metadata(pd.DataFrame(
            {'a': [1, 2, 3], 'b': [4, 5, 6],
             'c': [7, 8, 9], 'd': [10, 11, 12]}, index=['id1', 'id2', 'id3']))
        self.assertEqual(obs, exp)

    def test_merging_three(self):
        md1 = qiime2.Metadata(pd.DataFrame(
            {'a': [1, 2, 3], 'b': [4, 5, 6]}, index=['id1', 'id2', 'id3']))
        md2 = qiime2.Metadata(pd.DataFrame(
            {'c': [7, 8, 9], 'd': [10, 11, 12]}, index=['id1', 'id2', 'id3']))
        md3 = qiime2.Metadata(pd.DataFrame(
            {'e': [13, 14, 15], 'f': [16, 17, 18]},
            index=['id1', 'id2', 'id3']))

        obs = md1.merge(md2, md3)

        exp = qiime2.Metadata(pd.DataFrame(
            {'a': [1, 2, 3], 'b': [4, 5, 6],
             'c': [7, 8, 9], 'd': [10, 11, 12],
             'e': [13, 14, 15], 'f': [16, 17, 18]},
            index=['id1', 'id2', 'id3']))
        self.assertEqual(obs, exp)

    def test_merging_unaligned_indices(self):
        md1 = qiime2.Metadata(pd.DataFrame(
            {'a': [1, 2, 3], 'b': [4, 5, 6]}, index=['id1', 'id2', 'id3']))
        md2 = qiime2.Metadata(pd.DataFrame(
            {'c': [9, 8, 7], 'd': [12, 11, 10]}, index=['id3', 'id2', 'id1']))
        md3 = qiime2.Metadata(pd.DataFrame(
            {'e': [13, 15, 14], 'f': [16, 18, 17]},
            index=['id1', 'id3', 'id2']))

        obs = md1.merge(md2, md3)

        exp = qiime2.Metadata(pd.DataFrame(
            {'a': [1, 2, 3], 'b': [4, 5, 6],
             'c': [7, 8, 9], 'd': [10, 11, 12],
             'e': [13, 14, 15], 'f': [16, 17, 18]},
            index=['id1', 'id2', 'id3']))
        self.assertEqual(obs, exp)

    def test_inner_join(self):
        md1 = qiime2.Metadata(pd.DataFrame(
            {'a': [1, 2, 3], 'b': [4, 5, 6]}, index=['id1', 'id2', 'id3']))
        md2 = qiime2.Metadata(pd.DataFrame(
            {'c': [7, 8, 9], 'd': [10, 11, 12]}, index=['id2', 'X', 'Y']))
        md3 = qiime2.Metadata(pd.DataFrame(
            {'e': [13, 14, 15], 'f': [16, 17, 18]}, index=['X', 'id3', 'id2']))

        # Single shared ID.
        obs = md1.merge(md2, md3)

        exp = qiime2.Metadata(pd.DataFrame(
            {'a': [2], 'b': [5], 'c': [7], 'd': [10], 'e': [15], 'f': [18]},
            index=['id2']))
        self.assertEqual(obs, exp)

        # Multiple shared IDs.
        obs = md1.merge(md3)

        exp = qiime2.Metadata(pd.DataFrame(
            {'a': [2, 3], 'b': [5, 6], 'e': [15, 14], 'f': [18, 17]},
            index=['id2', 'id3']))
        self.assertEqual(obs, exp)

    def test_index_and_column_merge_order(self):
        md1 = qiime2.Metadata(pd.DataFrame(
            [[1], [2], [3], [4]],
            index=['id1', 'id2', 'id3', 'id4'], columns=['a']))
        md2 = qiime2.Metadata(pd.DataFrame(
            [[5], [6], [7]], index=['id4', 'id3', 'id1'], columns=['b']))
        md3 = qiime2.Metadata(pd.DataFrame(
            [[8], [9], [10]], index=['id1', 'id4', 'id3'], columns=['c']))

        obs = md1.merge(md2, md3)

        exp = qiime2.Metadata(pd.DataFrame(
            [[1, 7, 8], [3, 6, 10], [4, 5, 9]],
            index=['id1', 'id3', 'id4'], columns=['a', 'b', 'c']))
        self.assertEqual(obs, exp)

        # Merging in different order produces different index/column order.
        obs = md2.merge(md1, md3)

        exp = qiime2.Metadata(pd.DataFrame(
            [[5, 4, 9], [6, 3, 10], [7, 1, 8]],
            index=['id4', 'id3', 'id1'], columns=['b', 'a', 'c']))
        self.assertEqual(obs, exp)

    def test_no_columns(self):
        md1 = qiime2.Metadata(pd.DataFrame({}, index=['id1', 'id2', 'id3']))
        md2 = qiime2.Metadata(pd.DataFrame({}, index=['id2', 'X', 'id1']))
        md3 = qiime2.Metadata(pd.DataFrame({}, index=['id1', 'id3', 'id2']))

        obs = md1.merge(md2, md3)

        exp = qiime2.Metadata(pd.DataFrame({}, index=['id1', 'id2']))
        self.assertEqual(obs, exp)

    def test_index_and_column_names(self):
        md1 = qiime2.Metadata(pd.DataFrame(
            {'a': [1, 2]},
            index=pd.Index(['id1', 'id2'], name='foo'),
            columns=pd.Index(['a'], name='abc')))
        md2 = qiime2.Metadata(pd.DataFrame(
            {'b': [3, 4]},
            index=pd.Index(['id1', 'id2'], name='bar'),
            columns=pd.Index(['b'], name='def')))

        obs = md1.merge(md2)

        exp = qiime2.Metadata(pd.DataFrame(
            {'a': [1, 2], 'b': [3, 4]}, index=['id1', 'id2']))
        self.assertEqual(obs, exp)
        self.assertIsNone(obs._dataframe.index.name)
        self.assertIsNone(obs._dataframe.columns.name)

    def test_no_artifacts(self):
        md1 = qiime2.Metadata(pd.DataFrame(
            {'a': [1, 2]}, index=['id1', 'id2']))
        md2 = qiime2.Metadata(pd.DataFrame(
            {'b': [3, 4]}, index=['id1', 'id2']))

        metadata = md1.merge(md2)

        self.assertEqual(metadata.artifacts, [])

    def test_with_artifacts(self):
        artifact1 = qiime2.Artifact.import_data('Mapping',
                                                {'a': '1', 'b': '2'})
        artifact2 = qiime2.Artifact.import_data('Mapping', {'d': '4'})

        md_from_artifact1 = qiime2.Metadata.from_artifact(artifact1)
        md_from_artifact2 = qiime2.Metadata.from_artifact(artifact2)
        md_no_artifact = qiime2.Metadata(pd.DataFrame(
            {'c': ['3', '42']}, index=['0', '1']))

        # Merge three metadata objects -- the first has an artifact, the second
        # does not, and the third has an artifact.
        obs = md_from_artifact1.merge(md_no_artifact, md_from_artifact2)

        exp = pd.DataFrame(
            {'a': '1', 'b': '2', 'c': '3', 'd': '4'}, index=['0'])
        pdt.assert_frame_equal(obs.to_dataframe(), exp)
        self.assertEqual(obs.artifacts, [artifact1, artifact2])

    def test_disjoint_indices(self):
        md1 = qiime2.Metadata(pd.DataFrame(
            {'a': [1, 2, 3], 'b': [4, 5, 6]}, index=['id1', 'id2', 'id3']))
        md2 = qiime2.Metadata(pd.DataFrame(
            {'c': [7, 8, 9], 'd': [10, 11, 12]}, index=['X', 'Y', 'Z']))

        with self.assertRaisesRegex(ValueError, 'no IDs shared'):
            md1.merge(md2)

    def test_duplicate_columns(self):
        md1 = qiime2.Metadata(pd.DataFrame(
            {'a': [1, 2], 'b': [3, 4]}, index=['id1', 'id2']))
        md2 = qiime2.Metadata(pd.DataFrame(
            {'c': [5, 6], 'b': [7, 8]}, index=['id1', 'id2']))

        with self.assertRaisesRegex(ValueError, "categories overlap: 'b'"):
            md1.merge(md2)

    def test_duplicate_columns_self_merge(self):
        md = qiime2.Metadata(pd.DataFrame(
            {'a': [1, 2], 'b': [3, 4]}, index=['id1', 'id2']))

        with self.assertRaisesRegex(ValueError,
                                    "categories overlap: 'a', 'b'"):
            md.merge(md)


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

    def test_query_by_index(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)

        actual = metadata.ids(where="id='S2' OR id='S1'")
        expected = {'S1', 'S2'}
        self.assertEqual(actual, expected)

    def test_no_columns(self):
        fp = pkg_resources.resource_filename(
            'qiime2.tests', 'data/metadata-no-columns.tsv')
        metadata = qiime2.Metadata.load(fp)

        obs = metadata.ids()

        exp = {'a', 'b', 'id'}
        self.assertEqual(obs, exp)


class TestEqualityOperators(unittest.TestCase, ReallyEqualMixin):
    def setUp(self):
        get_dummy_plugin()

    def test_type_mismatch(self):
        fp = pkg_resources.resource_filename(
            'qiime2.tests', 'data/metadata.tsv')
        md = qiime2.Metadata.load(fp)
        mdc = qiime2.MetadataCategory.load(fp, 'col1')

        self.assertIsInstance(md, qiime2.Metadata)
        self.assertIsInstance(mdc, qiime2.MetadataCategory)
        self.assertReallyNotEqual(md, mdc)

    def test_source_mismatch(self):
        # Metadata created from an artifact vs not shouldn't compare equal,
        # even if the data is the same.
        artifact = qiime2.Artifact.import_data('Mapping', {'a': '1', 'b': '2'})
        md_from_artifact = qiime2.Metadata.from_artifact(artifact)
        md_no_artifact = qiime2.Metadata(pd.DataFrame(
            {'a': '1', 'b': '2'}, index=['0']))

        pdt.assert_frame_equal(md_from_artifact.to_dataframe(),
                               md_no_artifact.to_dataframe())
        self.assertReallyNotEqual(md_from_artifact, md_no_artifact)

    def test_artifact_mismatch(self):
        # Metadata created from different artifacts shouldn't compare equal,
        # even if the data is the same.
        artifact1 = qiime2.Artifact.import_data('Mapping',
                                                {'a': '1', 'b': '2'})
        artifact2 = qiime2.Artifact.import_data('Mapping',
                                                {'a': '1', 'b': '2'})

        md1 = qiime2.Metadata.from_artifact(artifact1)
        md2 = qiime2.Metadata.from_artifact(artifact2)

        pdt.assert_frame_equal(md1.to_dataframe(), md2.to_dataframe())
        self.assertReallyNotEqual(md1, md2)

    def test_index_mismatch(self):
        md1 = qiime2.Metadata(pd.DataFrame({'a': '1', 'b': '2'}, index=['0']))
        md2 = qiime2.Metadata(pd.DataFrame({'a': '1', 'b': '2'}, index=['1']))

        self.assertReallyNotEqual(md1, md2)

    def test_column_mismatch(self):
        md1 = qiime2.Metadata(pd.DataFrame({'a': '1', 'b': '2'}, index=['0']))
        md2 = qiime2.Metadata(pd.DataFrame({'a': '1', 'c': '2'}, index=['0']))

        self.assertReallyNotEqual(md1, md2)

    def test_data_mismatch(self):
        md1 = qiime2.Metadata(pd.DataFrame({'a': '1', 'b': '3'}, index=['0']))
        md2 = qiime2.Metadata(pd.DataFrame({'a': '1', 'b': '2'}, index=['0']))

        self.assertReallyNotEqual(md1, md2)

    def test_equality_without_artifact(self):
        md1 = qiime2.Metadata(pd.DataFrame({'a': '1', 'b': '3'}, index=['0']))
        md2 = qiime2.Metadata(pd.DataFrame({'a': '1', 'b': '3'}, index=['0']))

        self.assertReallyEqual(md1, md2)

    def test_equality_with_artifact(self):
        artifact = qiime2.Artifact.import_data('Mapping', {'a': '1', 'b': '2'})
        md1 = qiime2.Metadata.from_artifact(artifact)
        md2 = qiime2.Metadata.from_artifact(artifact)

        self.assertReallyEqual(md1, md2)


if __name__ == '__main__':
    unittest.main()
