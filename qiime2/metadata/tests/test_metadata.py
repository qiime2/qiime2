# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources
import unittest

import pandas as pd
import pandas.util.testing as pdt
import numpy as np

from qiime2 import Artifact
from qiime2.metadata import (Metadata, CategoricalMetadataColumn,
                             NumericMetadataColumn, MetadataFileError)
from qiime2.core.testing.util import get_dummy_plugin, ReallyEqualMixin


class TestMetadata(unittest.TestCase):
    def setUp(self):
        self.illegal_chars = ['/', '\0', '\\', '*', '<', '>', '?', '|', '$']

    def test_valid_metadata(self):
        index = pd.Index(['a', 'b', 'c'], name='feature ID', dtype=object)
        df = pd.DataFrame({'col1': ['2', '1', '3']}, index=index, dtype=object)
        metadata = Metadata(df)

        obs_type = metadata.columns['col1'].type

        self.assertEqual(obs_type, 'categorical')

    def test_valid_metadata_str(self):
        index = pd.Index(['a', 'b', 'c'], name='sample id', dtype=str)
        df = pd.DataFrame({'col1': ['2', '1', '3']}, index=index, dtype=str)
        metadata = Metadata(df)

        obs_type = metadata.columns['col1'].type

        self.assertEqual(obs_type, 'categorical')

    def test_valid_metadata_id_column_only(self):
        index = pd.Index(['a', 'b', 'c'], name='ID', dtype=object)
        df = pd.DataFrame({}, index=index, dtype=object)
        metadata = Metadata(df)

        self.assertEqual(metadata.id_count, 3)
        self.assertEqual(metadata.column_count, 0)

    def test_case_insensitive_duplicate_ids(self):
        index = pd.Index(['a', 'b', 'A'], name='id')
        df = pd.DataFrame({'column': ['1', '2', '3']}, index=index)
        metadata = Metadata(df)

        self.assertEqual(metadata.ids, ('a', 'b', 'A'))

    def test_case_insensitive_duplicate_column_names(self):
        index = pd.Index(['a', 'b', 'c'], name='id')
        df = pd.DataFrame({'column': ['1', '2', '3'],
                           'Column': ['4', '5', '6']}, index=index)
        metadata = Metadata(df)

        self.assertEqual(set(metadata.columns), {'column', 'Column'})

    def test_artifacts(self):
        index = pd.Index(['a', 'b', 'c'], name='id', dtype=object)
        df = pd.DataFrame({'col1': ['2', '1', '3']}, index=index, dtype=object)

        metadata = Metadata(df)

        self.assertEqual(metadata.artifacts, ())

    def test_empty_metadata(self):
        # No index, no columns.
        df = pd.DataFrame([], index=pd.Index([], name='id'))

        with self.assertRaisesRegex(ValueError, 'Metadata.*empty'):
            Metadata(df)

        # No index, has columns.
        df = pd.DataFrame([], index=pd.Index([], name='id'),
                          columns=['a', 'b'])

        with self.assertRaisesRegex(ValueError, 'Metadata.*empty'):
            Metadata(df)

    def test_wrong_obj(self):
        with self.assertRaisesRegex(TypeError,
                                    'Metadata constructor.*pandas.DataFrame'):
            Metadata(pd.Series([1, 2, 3]))

        with self.assertRaisesRegex(TypeError,
                                    'Metadata constructor.*pandas.DataFrame'):
            Metadata({})

    def test_empty_column_name(self):
        with self.assertRaisesRegex(ValueError, 'empty metadata column name'):
            Metadata(pd.DataFrame({'': ['1']},
                     index=pd.Index(['a'], name='id')))

    def test_empty_id(self):
        with self.assertRaisesRegex(ValueError, 'empty metadata ID'):
            Metadata(pd.DataFrame({'a': ['1', '2']},
                     index=pd.Index(['id1', ''], name='id')))

    def test_invalid_column_dtype_w_null(self):
        columns = pd.Index(['a', float('nan')], dtype=object)
        with self.assertRaisesRegex(TypeError, 'non-string.*column name.*nan'):
            Metadata(pd.DataFrame([['val1', 'val2']],
                                  index=pd.Index(['x'], name='id'),
                                  columns=columns))

        columns = pd.Index(['a', None], dtype=object)
        with self.assertRaisesRegex(TypeError,
                                    'non-string.*column name.*None'):
            Metadata(pd.DataFrame([['val1', 'val2']],
                                  index=pd.Index(['x'], name='id'),
                                  columns=columns))

    def test_invalid_columns_dtype(self):
        with self.assertRaisesRegex(TypeError, 'non-string.*column name.*42'):
            Metadata(pd.DataFrame({'foo': ['a', 'b'], 42: ['c', 'd']},
                                  index=pd.Index(['0', '1'], name='id')))

    def test_invalid_index_dtype_w_null(self):
        index = pd.Index(['a', float('nan'), 'b'], name='id', dtype=object)
        with self.assertRaisesRegex(TypeError, 'non-string.*ID.*nan'):
            Metadata(pd.DataFrame({'x': [1, 2, 3], 'y': [4, 5, 6]},
                                  index=index))

        index = pd.Index(['a', None, 'c'], name='id', dtype=object)
        with self.assertRaisesRegex(TypeError, 'non-string.*ID.*None'):
            Metadata(pd.DataFrame({'x': [1, 2, 3], 'y': [4, 5, 6]},
                                  index=index))

    def test_invalid_index_dtype(self):
        with self.assertRaisesRegex(TypeError, 'non-string.*ID.*0'):
            Metadata(pd.DataFrame({'foo': ['a', 'b', 'c']},
                                  index=pd.Index([0, 1, 2], name='id')))

    def test_duplicate_columns(self):
        index = pd.Index(['a', 'b'], name='id', dtype=object)
        df = pd.DataFrame({'foo': [1, 2], 'bar': [3, 4]}, index=index)
        df.columns = ['foo', 'foo']

        with self.assertRaisesRegex(ValueError,
                                    "column names must be unique.*'foo'"):
            Metadata(df)

    def test_duplicate_indices(self):
        index = pd.Index(['a', 'b', 'b'], name='id', dtype=object)
        df = pd.DataFrame({'foo': [1, 2, 3]}, index=index)

        with self.assertRaisesRegex(ValueError,
                                    "IDs must be unique.*'b'"):
            Metadata(df)

    def test_to_dataframe(self):
        index = pd.Index(['a', 'b', 'c'], name='id', dtype=object)
        df = pd.DataFrame({'col1': [2, 1, 3],
                           'col2': ['2', '1', '3'],
                           'col3': ['4.0', '5.2', '6.9'],
                           'col4': [4.0, 5.2, 6.9],
                           'col5': ['a', 'b', 'c']},
                          index=index)
        metadata = Metadata(df)

        obs_df = metadata.to_dataframe()

        self.assertEqual(dict(obs_df.dtypes),
                         {'col1': np.float, 'col2': object, 'col3': object,
                          'col4': np.float, 'col5': object})


class TestMetadataLoad(unittest.TestCase):
    def test_comments_and_blank_lines(self):
        fp = pkg_resources.resource_filename(
            'qiime2.metadata.tests', 'data/comments-n-blanks.tsv')

        obs_md = Metadata.load(fp)

        exp_index = pd.Index(['id1', 'id2', 'id3'], name='ID',
                             dtype=object)
        exp_df = pd.DataFrame({'col1': [1.0, 2.0, 3.0],
                               'col2': ['a', 'b', 'c'],
                               'col3': ['foo', 'bar', '42']},
                              index=exp_index)
        exp_md = Metadata(exp_df)

        self.assertEqual(obs_md, exp_md)

    def test_empty_rows(self):
        fp = pkg_resources.resource_filename(
            'qiime2.metadata.tests', 'data/empty-rows.tsv')

        obs_md = Metadata.load(fp)

        exp_index = pd.Index(['id1', 'id2', 'id3'], name='id', dtype=object)
        exp_df = pd.DataFrame({'col1': [1.0, 2.0, 3.0],
                               'col2': ['a', 'b', 'c'],
                               'col3': ['foo', 'bar', '42']},
                              index=exp_index)
        exp_md = Metadata(exp_df)

        self.assertEqual(obs_md, exp_md)

    def test_qiime1_mapping_file(self):
        fp = pkg_resources.resource_filename(
            'qiime2.metadata.tests', 'data/qiime1.tsv')

        obs_md = Metadata.load(fp)

        exp_index = pd.Index(['id1', 'id2', 'id3'], name='#SampleID',
                             dtype=object)
        exp_df = pd.DataFrame({'col1': [1.0, 2.0, 3.0],
                               'col2': ['a', 'b', 'c'],
                               'col3': ['foo', 'bar', '42']},
                              index=exp_index)
        exp_md = Metadata(exp_df)

        self.assertEqual(obs_md, exp_md)

    def test_qiime1_empty_mapping_file(self):
        fp = pkg_resources.resource_filename(
            'qiime2.metadata.tests', 'data/qiime1-empty.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    'at least one ID.*empty'):
            Metadata.load(fp)

    def test_no_columns(self):
        fp = pkg_resources.resource_filename(
            'qiime2.metadata.tests', 'data/no-columns.tsv')

        obs_md = Metadata.load(fp)

        exp_index = pd.Index(['a', 'b', 'my-id'], name='id', dtype=object)
        exp_df = pd.DataFrame({}, index=exp_index, dtype=object)
        exp_md = Metadata(exp_df)

        self.assertEqual(obs_md, exp_md)

    def test_does_not_cast_ids(self):
        fp = pkg_resources.resource_filename(
            'qiime2.metadata.tests', 'data/no-type-cast.tsv')

        obs_md = Metadata.load(fp)

        exp_index = pd.Index(['0.000001', '0.004000', '0.000000'],
                             dtype=object, name='id')
        exp_df = pd.DataFrame({'col1': [2.0, 1.0, 3.0],
                               'col2': ['b', 'b', 'c'],
                               'col3': [2.5, 4.2, -9.999]},
                              index=exp_index)
        exp_md = Metadata(exp_df)

        self.assertEqual(obs_md, exp_md)

    def test_artifacts(self):
        fp = pkg_resources.resource_filename(
            'qiime2.metadata.tests', 'data/simple.tsv')

        metadata = Metadata.load(fp)

        self.assertEqual(metadata.artifacts, ())

    def test_empty_id(self):
        fp = pkg_resources.resource_filename(
            'qiime2.metadata.tests', 'data/comments-n-blanks-n-empty-id.tsv')

        with self.assertRaisesRegex(MetadataFileError, 'empty metadata ID'):
            Metadata.load(fp)

    def test_empty_file(self):
        fp = pkg_resources.resource_filename(
            'qiime2.metadata.tests', 'data/empty')

        with self.assertRaisesRegex(MetadataFileError,
                                    'locate header.*file may be empty'):
            Metadata.load(fp)

    def test_jagged_trailing_columns(self):
        # Test case based on https://github.com/qiime2/qiime2/issues/335
        fp = pkg_resources.resource_filename(
            'qiime2.metadata.tests', 'data/jagged-trailing-columns.tsv')

        obs_md = Metadata.load(fp)

        exp_index = pd.Index(['id1', 'id2', 'id3'], name='id', dtype=object)
        exp_df = pd.DataFrame({'col1': [1.0, 2.0, 3.0],
                               'col2': ['a', 'b', 'c'],
                               'col3': ['foo', 'bar', '42']},
                              index=exp_index)
        exp_md = Metadata(exp_df)

        self.assertEqual(obs_md, exp_md)


class TestGetColumn(unittest.TestCase):
    def setUp(self):
        get_dummy_plugin()

    def test_artifacts_are_propagated(self):
        A = Artifact.import_data('Mapping', {'a': '1', 'b': '3'})
        md = A.view(Metadata)

        obs = md.get_column('b')

        # TODO update to use MetadataColumn.__eq__
        self.assertEqual(obs.artifacts, (A,))
        pdt.assert_series_equal(
            obs.to_series(),
            pd.Series(['3'], index=pd.Index(['0'], name='id'), name='b'))


class TestMerge(unittest.TestCase):
    def setUp(self):
        get_dummy_plugin()

    def test_merging_one(self):
        md = Metadata(pd.DataFrame(
            {'a': [1, 2, 3], 'b': [4, 5, 6]},
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))

        obs = md.merge()

        self.assertIsNot(obs, md)
        self.assertEqual(obs, md)

    def test_merging_two(self):
        md1 = Metadata(pd.DataFrame(
            {'a': [1, 2, 3], 'b': [4, 5, 6]},
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))
        md2 = Metadata(pd.DataFrame(
            {'c': [7, 8, 9], 'd': [10, 11, 12]},
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))

        obs = md1.merge(md2)

        exp = Metadata(pd.DataFrame(
            {'a': [1, 2, 3], 'b': [4, 5, 6],
             'c': [7, 8, 9], 'd': [10, 11, 12]},
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))
        self.assertEqual(obs, exp)

    def test_merging_three(self):
        md1 = Metadata(pd.DataFrame(
            {'a': [1, 2, 3], 'b': [4, 5, 6]},
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))
        md2 = Metadata(pd.DataFrame(
            {'c': [7, 8, 9], 'd': [10, 11, 12]},
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))
        md3 = Metadata(pd.DataFrame(
            {'e': [13, 14, 15], 'f': [16, 17, 18]},
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))

        obs = md1.merge(md2, md3)

        exp = Metadata(pd.DataFrame(
            {'a': [1, 2, 3], 'b': [4, 5, 6],
             'c': [7, 8, 9], 'd': [10, 11, 12],
             'e': [13, 14, 15], 'f': [16, 17, 18]},
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))
        self.assertEqual(obs, exp)

    def test_merging_unaligned_indices(self):
        md1 = Metadata(pd.DataFrame(
            {'a': [1, 2, 3], 'b': [4, 5, 6]},
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))
        md2 = Metadata(pd.DataFrame(
            {'c': [9, 8, 7], 'd': [12, 11, 10]},
            index=pd.Index(['id3', 'id2', 'id1'], name='id')))
        md3 = Metadata(pd.DataFrame(
            {'e': [13, 15, 14], 'f': [16, 18, 17]},
            index=pd.Index(['id1', 'id3', 'id2'], name='id')))

        obs = md1.merge(md2, md3)

        exp = Metadata(pd.DataFrame(
            {'a': [1, 2, 3], 'b': [4, 5, 6],
             'c': [7, 8, 9], 'd': [10, 11, 12],
             'e': [13, 14, 15], 'f': [16, 17, 18]},
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))
        self.assertEqual(obs, exp)

    def test_inner_join(self):
        md1 = Metadata(pd.DataFrame(
            {'a': [1, 2, 3], 'b': [4, 5, 6]},
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))
        md2 = Metadata(pd.DataFrame(
            {'c': [7, 8, 9], 'd': [10, 11, 12]},
            index=pd.Index(['id2', 'X', 'Y'], name='id')))
        md3 = Metadata(pd.DataFrame(
            {'e': [13, 14, 15], 'f': [16, 17, 18]},
            index=pd.Index(['X', 'id3', 'id2'], name='id')))

        # Single shared ID.
        obs = md1.merge(md2, md3)

        exp = Metadata(pd.DataFrame(
            {'a': [2], 'b': [5], 'c': [7], 'd': [10], 'e': [15], 'f': [18]},
            index=pd.Index(['id2'], name='id')))
        self.assertEqual(obs, exp)

        # Multiple shared IDs.
        obs = md1.merge(md3)

        exp = Metadata(pd.DataFrame(
            {'a': [2, 3], 'b': [5, 6], 'e': [15, 14], 'f': [18, 17]},
            index=pd.Index(['id2', 'id3'], name='id')))
        self.assertEqual(obs, exp)

    def test_index_and_column_merge_order(self):
        md1 = Metadata(pd.DataFrame(
            [[1], [2], [3], [4]],
            index=pd.Index(['id1', 'id2', 'id3', 'id4'], name='id'),
            columns=['a']))
        md2 = Metadata(pd.DataFrame(
            [[5], [6], [7]], index=pd.Index(['id4', 'id3', 'id1'], name='id'),
            columns=['b']))
        md3 = Metadata(pd.DataFrame(
            [[8], [9], [10]], index=pd.Index(['id1', 'id4', 'id3'], name='id'),
            columns=['c']))

        obs = md1.merge(md2, md3)

        exp = Metadata(pd.DataFrame(
            [[1, 7, 8], [3, 6, 10], [4, 5, 9]],
            index=pd.Index(['id1', 'id3', 'id4'], name='id'),
            columns=['a', 'b', 'c']))
        self.assertEqual(obs, exp)

        # Merging in different order produces different ID/column order.
        obs = md2.merge(md1, md3)

        exp = Metadata(pd.DataFrame(
            [[5, 4, 9], [6, 3, 10], [7, 1, 8]],
            index=pd.Index(['id4', 'id3', 'id1'], name='id'),
            columns=['b', 'a', 'c']))
        self.assertEqual(obs, exp)

    def test_id_column_only(self):
        md1 = Metadata(pd.DataFrame({},
                       index=pd.Index(['id1', 'id2', 'id3'], name='id')))
        md2 = Metadata(pd.DataFrame({},
                       index=pd.Index(['id2', 'X', 'id1'], name='id')))
        md3 = Metadata(pd.DataFrame({},
                       index=pd.Index(['id1', 'id3', 'id2'], name='id')))

        obs = md1.merge(md2, md3)

        exp = Metadata(
            pd.DataFrame({}, index=pd.Index(['id1', 'id2'], name='id')))
        self.assertEqual(obs, exp)

    def test_merged_id_column_name(self):
        md1 = Metadata(pd.DataFrame(
            {'a': [1, 2]},
            index=pd.Index(['id1', 'id2'], name='sample ID')))
        md2 = Metadata(pd.DataFrame(
            {'b': [3, 4]},
            index=pd.Index(['id1', 'id2'], name='feature ID')))

        obs = md1.merge(md2)

        exp = Metadata(pd.DataFrame(
            {'a': [1, 2], 'b': [3, 4]},
            index=pd.Index(['id1', 'id2'], name='id')))
        self.assertEqual(obs, exp)

    def test_no_artifacts(self):
        md1 = Metadata(pd.DataFrame(
            {'a': [1, 2]}, index=pd.Index(['id1', 'id2'], name='id')))
        md2 = Metadata(pd.DataFrame(
            {'b': [3, 4]}, index=pd.Index(['id1', 'id2'], name='id')))

        metadata = md1.merge(md2)

        self.assertEqual(metadata.artifacts, ())

    def test_with_artifacts(self):
        artifact1 = Artifact.import_data('Mapping', {'a': '1', 'b': '2'})
        artifact2 = Artifact.import_data('Mapping', {'d': '4'})

        md_from_artifact1 = artifact1.view(Metadata)
        md_from_artifact2 = artifact2.view(Metadata)
        md_no_artifact = Metadata(pd.DataFrame(
            {'c': ['3', '42']}, index=pd.Index(['0', '1'], name='id')))

        # Merge three metadata objects -- the first has an artifact, the second
        # does not, and the third has an artifact.
        obs_md = md_from_artifact1.merge(md_no_artifact, md_from_artifact2)

        exp_df = pd.DataFrame(
            {'a': '1', 'b': '2', 'c': '3', 'd': '4'},
            index=pd.Index(['0'], name='id'))
        exp_md = Metadata(exp_df)
        exp_md._add_artifacts((artifact1, artifact2))

        self.assertEqual(obs_md, exp_md)
        self.assertEqual(obs_md.artifacts, (artifact1, artifact2))

    def test_disjoint_indices(self):
        md1 = Metadata(pd.DataFrame(
            {'a': [1, 2, 3], 'b': [4, 5, 6]},
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))
        md2 = Metadata(pd.DataFrame(
            {'c': [7, 8, 9], 'd': [10, 11, 12]},
            index=pd.Index(['X', 'Y', 'Z'], name='id')))

        with self.assertRaisesRegex(ValueError, 'no IDs shared'):
            md1.merge(md2)

    def test_duplicate_columns(self):
        md1 = Metadata(pd.DataFrame(
            {'a': [1, 2], 'b': [3, 4]},
            index=pd.Index(['id1', 'id2'], name='id')))
        md2 = Metadata(pd.DataFrame(
            {'c': [5, 6], 'b': [7, 8]},
            index=pd.Index(['id1', 'id2'], name='id')))

        with self.assertRaisesRegex(ValueError, "columns overlap: 'b'"):
            md1.merge(md2)

    def test_duplicate_columns_self_merge(self):
        md = Metadata(pd.DataFrame(
            {'a': [1, 2], 'b': [3, 4]},
            index=pd.Index(['id1', 'id2'], name='id')))

        with self.assertRaisesRegex(ValueError, "columns overlap: 'a', 'b'"):
            md.merge(md)


class TestGetIDs(unittest.TestCase):
    def test_default(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = Metadata(df)

        actual = metadata.get_ids()
        expected = {'S1', 'S2', 'S3'}
        self.assertEqual(actual, expected)

    def test_incomplete_where(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='sampleid'))
        metadata = Metadata(df)

        where = "Subject='subject-1' AND SampleType="
        with self.assertRaises(ValueError):
            metadata.get_ids(where)

        where = "Subject="
        with self.assertRaises(ValueError):
            metadata.get_ids(where)

    def test_invalid_where(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='sampleid'))
        metadata = Metadata(df)

        where = "not-a-column-name='subject-1'"
        with self.assertRaises(ValueError):
            metadata.get_ids(where)

    def test_empty_result(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = Metadata(df)

        where = "Subject='subject-3'"
        actual = metadata.get_ids(where)
        expected = set()
        self.assertEqual(actual, expected)

    def test_simple_expression(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = Metadata(df)

        where = "Subject='subject-1'"
        actual = metadata.get_ids(where)
        expected = {'S1', 'S2'}
        self.assertEqual(actual, expected)

        where = "Subject='subject-2'"
        actual = metadata.get_ids(where)
        expected = {'S3'}
        self.assertEqual(actual, expected)

        where = "Subject='subject-3'"
        actual = metadata.get_ids(where)
        expected = set()
        self.assertEqual(actual, expected)

        where = "SampleType='gut'"
        actual = metadata.get_ids(where)
        expected = {'S1', 'S3'}
        self.assertEqual(actual, expected)

        where = "SampleType='tongue'"
        actual = metadata.get_ids(where)
        expected = {'S2'}
        self.assertEqual(actual, expected)

    def test_more_complex_expressions(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = Metadata(df)

        where = "Subject='subject-1' OR Subject='subject-2'"
        actual = metadata.get_ids(where)
        expected = {'S1', 'S2', 'S3'}
        self.assertEqual(actual, expected)

        where = "Subject='subject-1' AND Subject='subject-2'"
        actual = metadata.get_ids(where)
        expected = set()
        self.assertEqual(actual, expected)

        where = "Subject='subject-1' AND SampleType='gut'"
        actual = metadata.get_ids(where)
        expected = {'S1'}
        self.assertEqual(actual, expected)

    def test_query_by_id(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = Metadata(df)

        actual = metadata.get_ids(where="id='S2' OR id='S1'")
        expected = {'S1', 'S2'}
        self.assertEqual(actual, expected)

    def test_no_columns(self):
        fp = pkg_resources.resource_filename(
            'qiime2.metadata.tests', 'data/no-columns.tsv')
        metadata = Metadata.load(fp)

        obs = metadata.get_ids()

        exp = {'a', 'b', 'my-id'}
        self.assertEqual(obs, exp)


# TODO test with different ID column names (shouldn't compare equal)
class TestEqualityOperators(unittest.TestCase, ReallyEqualMixin):
    def setUp(self):
        get_dummy_plugin()

    def test_type_mismatch(self):
        fp = pkg_resources.resource_filename(
            'qiime2.metadata.tests', 'data/simple.tsv')
        md = Metadata.load(fp)
        mdc = md.get_column('col1')

        self.assertIsInstance(md, Metadata)
        self.assertIsInstance(mdc, NumericMetadataColumn)
        self.assertReallyNotEqual(md, mdc)

    def test_source_mismatch(self):
        # Metadata created from an artifact vs not shouldn't compare equal,
        # even if the data is the same.
        artifact = Artifact.import_data('Mapping', {'a': '1', 'b': '2'})
        md_from_artifact = artifact.view(Metadata)

        md_no_artifact = Metadata(md_from_artifact.to_dataframe())

        pdt.assert_frame_equal(md_from_artifact.to_dataframe(),
                               md_no_artifact.to_dataframe())
        self.assertReallyNotEqual(md_from_artifact, md_no_artifact)

    def test_artifact_mismatch(self):
        # Metadata created from different artifacts shouldn't compare equal,
        # even if the data is the same.
        artifact1 = Artifact.import_data('Mapping', {'a': '1', 'b': '2'})
        artifact2 = Artifact.import_data('Mapping', {'a': '1', 'b': '2'})

        md1 = artifact1.view(Metadata)
        md2 = artifact2.view(Metadata)

        pdt.assert_frame_equal(md1.to_dataframe(), md2.to_dataframe())
        self.assertReallyNotEqual(md1, md2)

    def test_id_mismatch(self):
        md1 = Metadata(pd.DataFrame({'a': '1', 'b': '2'},
                                    index=pd.Index(['0'], name='id')))
        md2 = Metadata(pd.DataFrame({'a': '1', 'b': '2'},
                                    index=pd.Index(['1'], name='id')))

        self.assertReallyNotEqual(md1, md2)

    def test_column_mismatch(self):
        md1 = Metadata(pd.DataFrame({'a': '1', 'b': '2'},
                                    index=pd.Index(['0'], name='id')))
        md2 = Metadata(pd.DataFrame({'a': '1', 'c': '2'},
                                    index=pd.Index(['0'], name='id')))

        self.assertReallyNotEqual(md1, md2)

    def test_data_mismatch(self):
        md1 = Metadata(pd.DataFrame({'a': '1', 'b': '3'},
                       index=pd.Index(['0'], name='id')))
        md2 = Metadata(pd.DataFrame({'a': '1', 'b': '2'},
                       index=pd.Index(['0'], name='id')))

        self.assertReallyNotEqual(md1, md2)

    def test_equality_without_artifact(self):
        md1 = Metadata(pd.DataFrame({'a': '1', 'b': '3'},
                                    index=pd.Index(['0'], name='id')))
        md2 = Metadata(pd.DataFrame({'a': '1', 'b': '3'},
                                    index=pd.Index(['0'], name='id')))

        self.assertReallyEqual(md1, md2)

    def test_equality_with_artifact(self):
        artifact = Artifact.import_data('Mapping', {'a': '1', 'b': '2'})
        md1 = artifact.view(Metadata)
        md2 = artifact.view(Metadata)

        self.assertReallyEqual(md1, md2)


class TestMetadataColumn(unittest.TestCase):
    def test_wrong_obj(self):
        with self.assertRaisesRegex(
                TypeError,
                'NumericMetadataColumn constructor.*pandas.Series'):
            NumericMetadataColumn(pd.DataFrame([[1, 2, 3]]))

        with self.assertRaisesRegex(
                TypeError,
                'CategoricalMetadataColumn constructor.*pandas.Series'):
            CategoricalMetadataColumn({})


if __name__ == '__main__':
    unittest.main()
