# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import unittest
import warnings

import pandas as pd
import pandas.util.testing as pdt
import numpy as np

from qiime2 import Artifact
from qiime2.metadata import (Metadata, CategoricalMetadataColumn,
                             NumericMetadataColumn)
from qiime2.core.testing.util import get_dummy_plugin, ReallyEqualMixin


class TestInvalidMetadataConstruction(unittest.TestCase):
    def test_non_dataframe(self):
        with self.assertRaisesRegex(
                TypeError, 'Metadata constructor.*DataFrame.*not.*Series'):
            Metadata(pd.Series([1, 2, 3], name='col',
                               index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_no_ids(self):
        with self.assertRaisesRegex(ValueError, 'Metadata.*at least one ID'):
            Metadata(pd.DataFrame({}, index=pd.Index([], name='id')))

        with self.assertRaisesRegex(ValueError, 'Metadata.*at least one ID'):
            Metadata(pd.DataFrame({'column': []},
                                  index=pd.Index([], name='id')))

    def test_invalid_id_header(self):
        # default index name
        with self.assertRaisesRegex(ValueError, r'Index\.name.*None'):
            Metadata(pd.DataFrame(
                {'col': [1, 2, 3]}, index=pd.Index(['a', 'b', 'c'])))

        with self.assertRaisesRegex(ValueError, r'Index\.name.*my-id-header'):
            Metadata(pd.DataFrame(
                {'col': [1, 2, 3]},
                index=pd.Index(['a', 'b', 'c'], name='my-id-header')))

    def test_non_str_id(self):
        with self.assertRaisesRegex(
                TypeError, 'non-string metadata ID.*type.*float.*nan'):
            Metadata(pd.DataFrame(
                {'col': [1, 2, 3]},
                index=pd.Index(['a', np.nan, 'c'], name='id')))

    def test_non_str_column_name(self):
        with self.assertRaisesRegex(
                TypeError, 'non-string metadata column name.*type.*'
                           'float.*nan'):
            Metadata(pd.DataFrame(
                {'col': [1, 2, 3],
                 np.nan: [4, 5, 6]},
                index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_empty_id(self):
        with self.assertRaisesRegex(
                ValueError, 'empty metadata ID.*at least one character'):
            Metadata(pd.DataFrame(
                {'col': [1, 2, 3]}, index=pd.Index(['a', '', 'c'], name='id')))

    def test_empty_column_name(self):
        with self.assertRaisesRegex(
                ValueError, 'empty metadata column name.*'
                            'at least one character'):
            Metadata(pd.DataFrame(
                {'col': [1, 2, 3],
                 '': [4, 5, 6]}, index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_leading_trailing_whitespace_id(self):
        with self.assertRaisesRegex(ValueError, "metadata ID.*leading or "
                                                "trailing whitespace.*' b '"):
            Metadata(pd.DataFrame(
                {'col': [1, 2, 3]},
                index=pd.Index(['a', ' b ', 'c'], name='id')))

    def test_leading_trailing_whitespace_column_name(self):
        with self.assertRaisesRegex(
                ValueError, "metadata column name.*leading or trailing "
                            "whitespace.*' col2 '"):
            Metadata(pd.DataFrame(
                {'col1': [1, 2, 3],
                 ' col2 ': [4, 5, 6]},
                index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_pound_sign_id(self):
        with self.assertRaisesRegex(
                ValueError, "metadata ID.*begins with a pound sign.*'#b'"):
            Metadata(pd.DataFrame(
                {'col': [1, 2, 3]},
                index=pd.Index(['a', '#b', 'c'], name='id')))

    def test_id_conflicts_with_id_header(self):
        with self.assertRaisesRegex(
                ValueError, "metadata ID 'sample-id'.*conflicts.*reserved.*"
                            "ID header"):
            Metadata(pd.DataFrame(
                {'col': [1, 2, 3]},
                index=pd.Index(['a', 'sample-id', 'c'], name='id')))

    def test_column_name_conflicts_with_id_header(self):
        with self.assertRaisesRegex(
                ValueError, "metadata column name 'featureid'.*conflicts.*"
                            "reserved.*ID header"):
            Metadata(pd.DataFrame(
                {'col': [1, 2, 3],
                 'featureid': [4, 5, 6]},
                index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_duplicate_ids(self):
        with self.assertRaisesRegex(ValueError, "Metadata IDs.*unique.*'a'"):
            Metadata(pd.DataFrame(
                {'col': [1, 2, 3]},
                index=pd.Index(['a', 'b', 'a'], name='id')))

    def test_duplicate_column_names(self):
        data = [[1, 2, 3],
                [4, 5, 6],
                [7, 8, 9]]
        with self.assertRaisesRegex(ValueError,
                                    "Metadata column names.*unique.*'col1'"):
            Metadata(pd.DataFrame(data, columns=['col1', 'col2', 'col1'],
                                  index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_unsupported_column_dtype(self):
        with self.assertRaisesRegex(
                TypeError, "Metadata column 'col2'.*unsupported.*dtype.*bool"):
            Metadata(pd.DataFrame(
                {'col1': [1, 2, 3],
                 'col2': [True, False, True]},
                index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_categorical_column_unsupported_type(self):
        with self.assertRaisesRegex(
                TypeError, "CategoricalMetadataColumn.*strings or missing "
                           r"values.*42\.5.*float.*'col2'"):
            Metadata(pd.DataFrame(
                {'col1': [1, 2, 3],
                 'col2': ['foo', 'bar', 42.5]},
                index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_categorical_column_empty_str(self):
        with self.assertRaisesRegex(
                ValueError, "CategoricalMetadataColumn.*empty strings.*"
                            "column 'col2'"):
            Metadata(pd.DataFrame(
                {'col1': [1, 2, 3],
                 'col2': ['foo', '', 'bar']},
                index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_categorical_column_leading_trailing_whitespace_value(self):
        with self.assertRaisesRegex(
                ValueError, "CategoricalMetadataColumn.*leading or trailing "
                            "whitespace characters.*Column 'col2'.*' bar '"):
            Metadata(pd.DataFrame(
                {'col1': [1, 2, 3],
                 'col2': ['foo', ' bar ', 'baz']},
                index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_numeric_column_infinity(self):
        with self.assertRaisesRegex(
                ValueError, "NumericMetadataColumn.*positive or negative "
                            "infinity.*column 'col2'"):
            Metadata(pd.DataFrame(
                {'col1': ['foo', 'bar', 'baz'],
                 'col2': [42, float('+inf'), 4.3]},
                index=pd.Index(['a', 'b', 'c'], name='id')))


class TestMetadataConstructionAndProperties(unittest.TestCase):
    def assertEqualColumns(self, obs_columns, exp):
        obs = [(name, props.type) for name, props in obs_columns.items()]
        self.assertEqual(obs, exp)

    def test_minimal(self):
        md = Metadata(pd.DataFrame({}, index=pd.Index(['a'], name='id')))

        self.assertEqual(md.id_count, 1)
        self.assertEqual(md.column_count, 0)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(md.ids, ('a',))
        self.assertEqualColumns(md.columns, [])

    def test_single_id(self):
        index = pd.Index(['id1'], name='id')
        df = pd.DataFrame({'col1': [1.0], 'col2': ['a'], 'col3': ['foo']},
                          index=index)
        md = Metadata(df)

        self.assertEqual(md.id_count, 1)
        self.assertEqual(md.column_count, 3)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(md.ids, ('id1',))
        self.assertEqualColumns(md.columns,
                                [('col1', 'numeric'), ('col2', 'categorical'),
                                 ('col3', 'categorical')])

    def test_no_columns(self):
        index = pd.Index(['id1', 'id2', 'foo'], name='id')
        df = pd.DataFrame({}, index=index)
        md = Metadata(df)

        self.assertEqual(md.id_count, 3)
        self.assertEqual(md.column_count, 0)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(md.ids, ('id1', 'id2', 'foo'))
        self.assertEqualColumns(md.columns, [])

    def test_single_column(self):
        index = pd.Index(['id1', 'a', 'my-id'], name='id')
        df = pd.DataFrame({'column': ['foo', 'bar', 'baz']}, index=index)
        md = Metadata(df)

        self.assertEqual(md.id_count, 3)
        self.assertEqual(md.column_count, 1)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(md.ids, ('id1', 'a', 'my-id'))
        self.assertEqualColumns(md.columns, [('column', 'categorical')])

    def test_retains_column_order(self):
        # Supply DataFrame constructor with explicit column ordering instead of
        # a dict.
        index = pd.Index(['id1', 'id2', 'id3'], name='id')
        columns = ['z', 'a', 'ch']
        data = [
            [1.0, 'a', 'foo'],
            [2.0, 'b', 'bar'],
            [3.0, 'c', '42']
        ]
        df = pd.DataFrame(data, index=index, columns=columns)
        md = Metadata(df)

        self.assertEqual(md.id_count, 3)
        self.assertEqual(md.column_count, 3)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(md.ids, ('id1', 'id2', 'id3'))
        self.assertEqualColumns(md.columns,
                                [('z', 'numeric'), ('a', 'categorical'),
                                 ('ch', 'categorical')])

    def test_supported_id_headers(self):
        case_insensitive = {
            'id', 'sampleid', 'sample id', 'sample-id', 'featureid',
            'feature id', 'feature-id'
        }

        exact_match = {
            '#SampleID', '#Sample ID', '#OTUID', '#OTU ID', 'sample_name'
        }

        # Build a set of supported headers, including exact matches and headers
        # with different casing.
        headers = set()
        for header in case_insensitive:
            headers.add(header)
            headers.add(header.upper())
            headers.add(header.title())
        for header in exact_match:
            headers.add(header)

        count = 0
        for header in headers:
            index = pd.Index(['id1', 'id2'], name=header)
            df = pd.DataFrame({'column': ['foo', 'bar']}, index=index)
            md = Metadata(df)

            self.assertEqual(md.id_header, header)
            count += 1

        # Since this test case is a little complicated, make sure that the
        # expected number of comparisons are happening.
        self.assertEqual(count, 26)

    def test_recommended_ids(self):
        index = pd.Index(['c6ca034a-223f-40b4-a0e0-45942912a5ea', 'My.ID'],
                         name='id')
        df = pd.DataFrame({'col1': ['foo', 'bar']}, index=index)
        md = Metadata(df)

        self.assertEqual(md.id_count, 2)
        self.assertEqual(md.column_count, 1)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(md.ids,
                         ('c6ca034a-223f-40b4-a0e0-45942912a5ea', 'My.ID'))
        self.assertEqualColumns(md.columns, [('col1', 'categorical')])

    def test_non_standard_characters(self):
        index = pd.Index(['©id##1', '((id))2', "'id_3<>'", '"id#4"',
                          'i d\r\t\n5'], name='id')
        columns = ['↩c@l1™', 'col(#2)', "#col'3", '"<col_4>"', 'col\t  \r\n5']
        data = [
            ['ƒoo', '(foo)', '#f o #o', 'fo\ro', np.nan],
            ["''2''", 'b#r', 'ba\nr', np.nan, np.nan],
            ['b"ar', 'c\td', '4\r\n2', np.nan, np.nan],
            ['b__a_z', '<42>', '>42', np.nan, np.nan],
            ['baz', np.nan, '42']
        ]
        df = pd.DataFrame(data, index=index, columns=columns)
        md = Metadata(df)

        self.assertEqual(md.id_count, 5)
        self.assertEqual(md.column_count, 5)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(
            md.ids, ('©id##1', '((id))2', "'id_3<>'", '"id#4"', 'i d\r\t\n5'))
        self.assertEqualColumns(md.columns, [('↩c@l1™', 'categorical'),
                                             ('col(#2)', 'categorical'),
                                             ("#col'3", 'categorical'),
                                             ('"<col_4>"', 'categorical'),
                                             ('col\t  \r\n5', 'numeric')])

    def test_missing_data(self):
        index = pd.Index(['None', 'nan', 'NA', 'foo'], name='id')
        df = pd.DataFrame(collections.OrderedDict([
            ('col1', [1.0, np.nan, np.nan, np.nan]),
            ('NA', [np.nan, np.nan, np.nan, np.nan]),
            ('col3', ['null', 'N/A', np.nan, 'NA']),
            ('col4', np.array([np.nan, np.nan, np.nan, np.nan],
                              dtype=object))]),
            index=index)
        md = Metadata(df)

        self.assertEqual(md.id_count, 4)
        self.assertEqual(md.column_count, 4)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(md.ids, ('None', 'nan', 'NA', 'foo'))
        self.assertEqualColumns(md.columns, [('col1', 'numeric'),
                                             ('NA', 'numeric'),
                                             ('col3', 'categorical'),
                                             ('col4', 'categorical')])

    def test_does_not_cast_ids_or_column_names(self):
        index = pd.Index(['0.000001', '0.004000', '0.000000'], dtype=object,
                         name='id')
        columns = ['42.0', '1000', '-4.2']
        data = [
            [2.0, 'b', 2.5],
            [1.0, 'b', 4.2],
            [3.0, 'c', -9.999]
        ]
        df = pd.DataFrame(data, index=index, columns=columns)
        md = Metadata(df)

        self.assertEqual(md.id_count, 3)
        self.assertEqual(md.column_count, 3)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(md.ids, ('0.000001', '0.004000', '0.000000'))
        self.assertEqualColumns(md.columns, [('42.0', 'numeric'),
                                             ('1000', 'categorical'),
                                             ('-4.2', 'numeric')])

    def test_mixed_column_types(self):
        md = Metadata(
            pd.DataFrame({'col0': [1.0, 2.0, 3.0],
                          'col1': ['a', 'b', 'c'],
                          'col2': ['foo', 'bar', '42'],
                          'col3': ['1.0', '2.5', '-4.002'],
                          'col4': [1, 2, 3],
                          'col5': [1, 2, 3.5],
                          'col6': [1e-4, -0.0002, np.nan],
                          'col7': ['cat', np.nan, 'dog'],
                          'col8': ['a', 'a', 'a'],
                          'col9': [0, 0, 0]},
                         index=pd.Index(['id1', 'id2', 'id3'], name='id')))

        self.assertEqual(md.id_count, 3)
        self.assertEqual(md.column_count, 10)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(md.ids, ('id1', 'id2', 'id3'))
        self.assertEqualColumns(md.columns, [('col0', 'numeric'),
                                             ('col1', 'categorical'),
                                             ('col2', 'categorical'),
                                             ('col3', 'categorical'),
                                             ('col4', 'numeric'),
                                             ('col5', 'numeric'),
                                             ('col6', 'numeric'),
                                             ('col7', 'categorical'),
                                             ('col8', 'categorical'),
                                             ('col9', 'numeric')])

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


class TestSourceArtifacts(unittest.TestCase):
    def setUp(self):
        self.md = Metadata(pd.DataFrame(
            {'col': [1, 2, 3]}, index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_no_source_artifacts(self):
        self.assertEqual(self.md.artifacts, ())

    def test_add_zero_artifacts(self):
        self.md._source_artifact = None

        self.assertEqual(self.md.artifacts, ())

    def test_source_artifacts(self):
        artifact1 = Artifact.import_data('Mapping', {'a': '1', 'b': '2'})
        md1 = artifact1.view(Metadata)

        artifact2 = Artifact.import_data('Mapping', {'a': '3', 'b': '4'})
        md2 = artifact2.view(Metadata)

        artifact3 = Artifact.import_data('Mapping', {'c': '5', 'd': '6'})
        md3 = artifact3.view(Metadata)
        self.assertEqual(artifact3, md3._source_artifact)

        merged_md = md1.merge(md2, md3)
        # The artifacts do not come out in a completely deterministic order in
        # the tuple, it's alphabetical by column name which varies with the
        # ids
        self.assertEqual(len(merged_md.artifacts), 3)
        self.assertIn(artifact1, merged_md.artifacts)
        self.assertIn(artifact2, merged_md.artifacts)
        self.assertIn(artifact3, merged_md.artifacts)

    def test_add_non_artifact(self):
        with self.assertRaisesRegex(TypeError, "Artifact object.*42"):
            self.md._add_source_artifact(42)

        # Test that the object hasn't been mutated.
        self.assertEqual(self.md.artifacts, ())


class TestRepr(unittest.TestCase):
    def test_singular(self):
        md = Metadata(pd.DataFrame({'col1': [42]},
                                   index=pd.Index(['a'], name='id')))

        obs = repr(md)

        self.assertIn('Metadata', obs)
        self.assertIn('1 ID x 1 column', obs)
        self.assertIn("col1: ColumnProperties(type='numeric')", obs)

    def test_plural(self):
        md = Metadata(pd.DataFrame({'col1': [42, 42], 'col2': ['foo', 'bar']},
                                   index=pd.Index(['a', 'b'], name='id')))

        obs = repr(md)

        self.assertIn('Metadata', obs)
        self.assertIn('2 IDs x 2 columns', obs)
        self.assertIn("col1: ColumnProperties(type='numeric')", obs)
        self.assertIn("col2: ColumnProperties(type='categorical')", obs)

    def test_column_name_padding(self):
        data = [[0, 42, 'foo']]
        index = pd.Index(['my-id'], name='id')
        columns = ['col1', 'longer-column-name', 'c']
        md = Metadata(pd.DataFrame(data, index=index, columns=columns))

        obs = repr(md)

        self.assertIn('Metadata', obs)
        self.assertIn('1 ID x 3 columns', obs)
        self.assertIn(
            "col1:               ColumnProperties(type='numeric')", obs)
        self.assertIn(
            "longer-column-name: ColumnProperties(type='numeric')", obs)
        self.assertIn(
            "c:                  ColumnProperties(type='categorical')", obs)


class TestEqualityOperators(unittest.TestCase, ReallyEqualMixin):
    def setUp(self):
        get_dummy_plugin()

    def test_type_mismatch(self):
        md = Metadata(
            pd.DataFrame({'col1': [1.0, 2.0, 3.0],
                          'col2': ['a', 'b', 'c'],
                          'col3': ['foo', 'bar', '42']},
                         index=pd.Index(['id1', 'id2', 'id3'], name='id')))
        mdc = md.get_column('col1')

        self.assertIsInstance(md, Metadata)
        self.assertIsInstance(mdc, NumericMetadataColumn)
        self.assertReallyNotEqual(md, mdc)

    def test_id_header_mismatch(self):
        data = {'col1': ['foo', 'bar'], 'col2': [42, 43]}
        md1 = Metadata(pd.DataFrame(
            data, index=pd.Index(['id1', 'id2'], name='id')))
        md2 = Metadata(pd.DataFrame(
            data, index=pd.Index(['id1', 'id2'], name='ID')))

        self.assertReallyNotEqual(md1, md2)

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

    def test_column_name_mismatch(self):
        md1 = Metadata(pd.DataFrame({'a': '1', 'b': '2'},
                                    index=pd.Index(['0'], name='id')))
        md2 = Metadata(pd.DataFrame({'a': '1', 'c': '2'},
                                    index=pd.Index(['0'], name='id')))

        self.assertReallyNotEqual(md1, md2)

    def test_column_type_mismatch(self):
        md1 = Metadata(pd.DataFrame({'col1': ['42', '43']},
                                    index=pd.Index(['id1', 'id2'], name='id')))
        md2 = Metadata(pd.DataFrame({'col1': [42, 43]},
                                    index=pd.Index(['id1', 'id2'], name='id')))

        self.assertReallyNotEqual(md1, md2)

    def test_column_order_mismatch(self):
        index = pd.Index(['id1', 'id2'], name='id')
        md1 = Metadata(pd.DataFrame([[42, 'foo'], [43, 'bar']], index=index,
                                    columns=['z', 'a']))
        md2 = Metadata(pd.DataFrame([['foo', 42], ['bar', 43]], index=index,
                                    columns=['a', 'z']))

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

    def test_equality_with_missing_data(self):
        md1 = Metadata(pd.DataFrame(
            {'col1': [1, np.nan, 4.2],
             'col2': [np.nan, 'foo', np.nan]},
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))
        md2 = Metadata(pd.DataFrame(
            {'col1': [1, np.nan, 4.2],
             'col2': [np.nan, 'foo', np.nan]},
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))

        self.assertReallyEqual(md1, md2)


class TestToDataframe(unittest.TestCase):
    def test_minimal(self):
        df = pd.DataFrame({}, index=pd.Index(['id1'], name='id'))
        md = Metadata(df)

        obs = md.to_dataframe()

        pdt.assert_frame_equal(obs, df)

    def test_id_header_preserved(self):
        df = pd.DataFrame({'col1': [42, 2.5], 'col2': ['foo', 'bar']},
                          index=pd.Index(['id1', 'id2'], name='#SampleID'))
        md = Metadata(df)

        obs = md.to_dataframe()

        pdt.assert_frame_equal(obs, df)
        self.assertEqual(obs.index.name, '#SampleID')

    def test_dataframe_copy(self):
        df = pd.DataFrame({'col1': [42, 2.5], 'col2': ['foo', 'bar']},
                          index=pd.Index(['id1', 'id2'], name='id'))
        md = Metadata(df)

        obs = md.to_dataframe()

        pdt.assert_frame_equal(obs, df)
        self.assertIsNot(obs, df)

    def test_retains_column_order(self):
        index = pd.Index(['id1', 'id2'], name='id')
        columns = ['z', 'a', 'ch']
        data = [
            [1.0, 'a', 'foo'],
            [2.0, 'b', 'bar']
        ]
        df = pd.DataFrame(data, index=index, columns=columns)
        md = Metadata(df)

        obs = md.to_dataframe()

        pdt.assert_frame_equal(obs, df)
        self.assertEqual(obs.columns.tolist(), ['z', 'a', 'ch'])

    def test_missing_data(self):
        # Different missing data representations should be normalized to np.nan
        index = pd.Index(['None', 'nan', 'NA', 'id1'], name='id')
        df = pd.DataFrame(collections.OrderedDict([
            ('col1', [42.5, np.nan, float('nan'), 3]),
            ('NA', [np.nan, 'foo', float('nan'), None]),
            ('col3', ['null', 'N/A', np.nan, 'NA']),
            ('col4', np.array([np.nan, np.nan, np.nan, np.nan],
                              dtype=object))]),
            index=index)
        md = Metadata(df)

        obs = md.to_dataframe()

        exp = pd.DataFrame(collections.OrderedDict([
            ('col1', [42.5, np.nan, np.nan, 3.0]),
            ('NA', [np.nan, 'foo', np.nan, np.nan]),
            ('col3', ['null', 'N/A', np.nan, 'NA']),
            ('col4', np.array([np.nan, np.nan, np.nan, np.nan],
                              dtype=object))]),
            index=index)

        pdt.assert_frame_equal(obs, exp)
        self.assertEqual(obs.dtypes.to_dict(),
                         {'col1': np.float64, 'NA': object, 'col3': object,
                          'col4': object})
        self.assertTrue(np.isnan(obs['col1']['NA']))
        self.assertTrue(np.isnan(obs['NA']['NA']))
        self.assertTrue(np.isnan(obs['NA']['id1']))

    def test_dtype_int_normalized_to_dtype_float(self):
        index = pd.Index(['id1', 'id2', 'id3'], name='id')
        df = pd.DataFrame({'col1': [42, -43, 0],
                           'col2': [42.0, -43.0, 0.0],
                           'col3': [42, np.nan, 0]},
                          index=index)

        self.assertEqual(df.dtypes.to_dict(),
                         {'col1': np.int64, 'col2': np.float64,
                          'col3': np.float64})

        md = Metadata(df)
        obs = md.to_dataframe()

        exp = pd.DataFrame({'col1': [42.0, -43.0, 0.0],
                            'col2': [42.0, -43.0, 0.0],
                            'col3': [42.0, np.nan, 0.0]},
                           index=index)

        pdt.assert_frame_equal(obs, exp)
        self.assertEqual(obs.dtypes.to_dict(),
                         {'col1': np.float64, 'col2': np.float64,
                          'col3': np.float64})


class TestGetColumn(unittest.TestCase):
    def setUp(self):
        get_dummy_plugin()

    def test_column_name_not_found(self):
        df = pd.DataFrame({'col1': [42, 2.5], 'col2': ['foo', 'bar']},
                          index=pd.Index(['id1', 'id2'], name='id'))
        md = Metadata(df)

        with self.assertRaisesRegex(ValueError,
                                    "'col3'.*not a column.*'col1', 'col2'"):
            md.get_column('col3')

    def test_artifacts_are_propagated(self):
        A = Artifact.import_data('Mapping', {'a': '1', 'b': '3'})
        md = A.view(Metadata)

        obs = md.get_column('b')

        exp = CategoricalMetadataColumn(
            pd.Series(['3'], name='b', index=pd.Index(['0'], name='id')))
        exp._source_metadata = md

        self.assertEqual(obs, exp)
        self.assertEqual(obs.artifacts, (A,))

    def test_categorical_column(self):
        df = pd.DataFrame({'col1': [42, 2.5], 'col2': ['foo', 'bar']},
                          index=pd.Index(['id1', 'id2'], name='id'))
        md = Metadata(df)

        obs = md.get_column('col2')

        exp = CategoricalMetadataColumn(
            pd.Series(['foo', 'bar'], name='col2',
                      index=pd.Index(['id1', 'id2'], name='id')))

        self.assertEqual(obs, exp)

    def test_numeric_column(self):
        df = pd.DataFrame({'col1': [42, 2.5], 'col2': ['foo', 'bar']},
                          index=pd.Index(['id1', 'id2'], name='id'))
        md = Metadata(df)

        obs = md.get_column('col1')

        exp = NumericMetadataColumn(
            pd.Series([42, 2.5], name='col1',
                      index=pd.Index(['id1', 'id2'], name='id')))

        self.assertEqual(obs, exp)

    def test_id_header_preserved(self):
        df = pd.DataFrame({'col1': [42, 2.5], 'col2': ['foo', 'bar']},
                          index=pd.Index(['a', 'b'], name='#OTU ID'))
        md = Metadata(df)

        obs = md.get_column('col1')

        exp = NumericMetadataColumn(
            pd.Series([42, 2.5], name='col1',
                      index=pd.Index(['a', 'b'], name='#OTU ID')))

        self.assertEqual(obs, exp)
        self.assertEqual(obs.id_header, '#OTU ID')


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

    def test_query_by_alternate_id_header(self):
        metadata = Metadata(pd.DataFrame(
            {}, index=pd.Index(['id1', 'id2', 'id3'], name='#OTU ID')))

        obs = metadata.get_ids(where="\"#OTU ID\" IN ('id2', 'id3')")

        exp = {'id2', 'id3'}
        self.assertEqual(obs, exp)

    def test_no_columns(self):
        metadata = Metadata(
            pd.DataFrame({}, index=pd.Index(['a', 'b', 'my-id'], name='id')))

        obs = metadata.get_ids()

        exp = {'a', 'b', 'my-id'}
        self.assertEqual(obs, exp)

    def test_query_mixed_column_types(self):
        df = pd.DataFrame({'Name': ['Foo', 'Bar', 'Baz', 'Baaz'],
                           # numbers that would sort incorrectly as strings
                           'Age': [9, 10, 11, 101],
                           'Age_Str': ['9', '10', '11', '101'],
                           'Weight': [80.5, 85.3, np.nan, 120.0]},
                          index=pd.Index(['S1', 'S2', 'S3', 'S4'], name='id'))
        metadata = Metadata(df)

        # string pattern matching
        obs = metadata.get_ids(where="Name LIKE 'Ba_'")
        exp = {'S2', 'S3'}
        self.assertEqual(obs, exp)

        # string comparison
        obs = metadata.get_ids(where="Age_Str >= 11")
        exp = {'S1', 'S3'}
        self.assertEqual(obs, exp)

        # numeric comparison
        obs = metadata.get_ids(where="Age >= 11")
        exp = {'S3', 'S4'}
        self.assertEqual(obs, exp)

        # numeric comparison with missing data
        obs = metadata.get_ids(where="Weight < 100")
        exp = {'S1', 'S2'}
        self.assertEqual(obs, exp)

    def test_column_with_space_in_name(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'Sample Type': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = Metadata(df)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            metadata.get_ids()
            # The list of captured warnings should be empty
            self.assertFalse(w)


class TestMerge(unittest.TestCase):
    def setUp(self):
        get_dummy_plugin()

    def test_merging_nothing(self):
        md = Metadata(pd.DataFrame(
            {'a': [1, 2, 3], 'b': [4, 5, 6]},
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))

        with self.assertRaisesRegex(ValueError,
                                    'At least one Metadata.*nothing to merge'):
            md.merge()

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
                       index=pd.Index(['id1', 'A', 'id2'], name='id')))

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

    def test_merging_preserves_column_types(self):
        # Test that column types remain the same even if a categorical column
        # *could* be reinterpreted as numeric after the merge.
        md1 = Metadata(pd.DataFrame(
            {'a': [1, 2, 3],
             'b': [np.nan, np.nan, np.nan]},
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))
        md2 = Metadata(pd.DataFrame(
            {'c': ['1', 'foo', '3'],
             'd': np.array([np.nan, np.nan, np.nan], dtype=object)},
            index=pd.Index(['id1', 'id4', 'id3'], name='id')))

        obs = md1.merge(md2)

        exp = Metadata(pd.DataFrame(
            {'a': [1, 3], 'b': [np.nan, np.nan], 'c': ['1', '3'],
             'd': np.array([np.nan, np.nan], dtype=object)},
            index=pd.Index(['id1', 'id3'], name='id')))

        self.assertEqual(obs, exp)
        self.assertEqual(obs.columns['a'].type, 'numeric')
        self.assertEqual(obs.columns['b'].type, 'numeric')
        self.assertEqual(obs.columns['c'].type, 'categorical')
        self.assertEqual(obs.columns['d'].type, 'categorical')

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
        exp_md._column_sources = {'a': md_from_artifact1,
                                  'b': md_from_artifact1,
                                  'c': md_no_artifact,
                                  'd': md_from_artifact2}

        self.assertEqual(exp_md, obs_md)
        self.assertEqual(2, len(obs_md.artifacts))
        self.assertIn(artifact1, obs_md.artifacts)
        self.assertIn(artifact2, obs_md.artifacts)

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

        exp = Metadata(pd.DataFrame(
            {'a': [1, 2], 'b [31b345c7b509cf1a0bb0d8f56c8e08b2]': [3, 4],
             'c': [5, 6], 'b [af19e323d6b6f4feb8baaf5c5c6f1221]': [7, 8]},
            index=pd.Index(['id1', 'id2'], name='id')))

        obs = md1.merge(md2)
        self.assertEqual(obs, exp)


class TestFilterIDs(unittest.TestCase):
    def setUp(self):
        get_dummy_plugin()

    def test_supports_iterable(self):
        md = Metadata(pd.DataFrame(
            {'col1': [1, 2, 3], 'col2': ['foo', 'bar', 'baz']},
            index=pd.Index(['a', 'b', 'c'], name='id')))

        obs = md.filter_ids(iter({'a', 'c'}))

        exp = Metadata(pd.DataFrame(
            {'col1': [1, 3], 'col2': ['foo', 'baz']},
            index=pd.Index(['a', 'c'], name='id')))

        self.assertEqual(obs, exp)

    def test_keep_all(self):
        md = Metadata(pd.DataFrame(
            {'col1': [1, 2, 3], 'col2': ['foo', 'bar', 'baz']},
            index=pd.Index(['a', 'b', 'c'], name='id')))

        obs = md.filter_ids({'a', 'b', 'c'})

        self.assertEqual(obs, md)
        self.assertIsNot(obs, md)

    def test_keep_multiple(self):
        md = Metadata(pd.DataFrame(
            {'col1': [1, 2, 3], 'col2': ['foo', 'bar', 'baz']},
            index=pd.Index(['a', 'b', 'c'], name='id')))

        obs = md.filter_ids({'a', 'c'})

        exp = Metadata(pd.DataFrame(
            {'col1': [1, 3], 'col2': ['foo', 'baz']},
            index=pd.Index(['a', 'c'], name='id')))

        self.assertEqual(obs, exp)

    def test_keep_one(self):
        md = Metadata(pd.DataFrame(
            {'col1': [1, 2, 3], 'col2': ['foo', 'bar', 'baz']},
            index=pd.Index(['a', 'b', 'c'], name='id')))

        obs = md.filter_ids({'b'})

        exp = Metadata(pd.DataFrame(
            {'col1': [2], 'col2': ['bar']}, index=pd.Index(['b'], name='id')))

        self.assertEqual(obs, exp)

    def test_filtering_preserves_column_types(self):
        # Test that column types remain the same even if a categorical column
        # *could* be reinterpreted as numeric after the filter.
        md = Metadata(pd.DataFrame(
            {'a': [1, 2, 3],
             'b': [np.nan, np.nan, np.nan],
             'c': ['1', 'foo', '3'],
             'd': np.array([np.nan, np.nan, np.nan], dtype=object)},
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))

        obs = md.filter_ids({'id1', 'id3'})

        exp = Metadata(pd.DataFrame(
            {'a': [1, 3], 'b': [np.nan, np.nan], 'c': ['1', '3'],
             'd': np.array([np.nan, np.nan], dtype=object)},
            index=pd.Index(['id1', 'id3'], name='id')))

        self.assertEqual(obs, exp)
        self.assertEqual(obs.columns['a'].type, 'numeric')
        self.assertEqual(obs.columns['b'].type, 'numeric')
        self.assertEqual(obs.columns['c'].type, 'categorical')
        self.assertEqual(obs.columns['d'].type, 'categorical')

    def test_alternate_id_header(self):
        md = Metadata(pd.DataFrame(
            {'col1': [1, 2, 3, 4], 'col2': ['foo', 'bar', 'baz', 'bazz']},
            index=pd.Index(['a', 'b', 'c', 'd'], name='#Sample ID')))

        obs = md.filter_ids({'b', 'd'})

        exp = Metadata(pd.DataFrame(
            {'col1': [2, 4], 'col2': ['bar', 'bazz']},
            index=pd.Index(['b', 'd'], name='#Sample ID')))

        self.assertEqual(obs, exp)

    def test_retains_column_order(self):
        data = [[1, 'foo', 'cat'], [2, 'bar', 'dog'], [3, 'baz', 'bat']]
        md = Metadata(pd.DataFrame(
            data, columns=['z', 'a', 'ch'],
            index=pd.Index(['a', 'b', 'c'], name='id')))

        obs = md.filter_ids({'b', 'c'})

        exp_data = [[2, 'bar', 'dog'], [3, 'baz', 'bat']]
        exp = Metadata(pd.DataFrame(
            exp_data, columns=['z', 'a', 'ch'],
            index=pd.Index(['b', 'c'], name='id')))

        self.assertEqual(obs, exp)

    def test_no_artifacts(self):
        md = Metadata(pd.DataFrame(
            {'col1': [1, 2, 3], 'col2': ['foo', 'bar', 'baz']},
            index=pd.Index(['a', 'b', 'c'], name='id')))

        self.assertEqual(md.artifacts, ())

        filtered = md.filter_ids({'b'})

        self.assertEqual(filtered.artifacts, ())

    def test_empty_ids_to_keep(self):
        md = Metadata(pd.DataFrame(
            {'col1': [1, 2, 3], 'col2': ['foo', 'bar', 'baz']},
            index=pd.Index(['a', 'b', 'c'], name='id')))

        with self.assertRaisesRegex(ValueError,
                                    'ids_to_keep.*at least one ID'):
            md.filter_ids({})

    def test_duplicate_ids_to_keep(self):
        md = Metadata(pd.DataFrame(
            {'col1': [1, 2, 3], 'col2': ['foo', 'bar', 'baz']},
            index=pd.Index(['a', 'b', 'c'], name='id')))

        with self.assertRaisesRegex(ValueError,
                                    "ids_to_keep.*unique IDs.*'b'"):
            md.filter_ids(['b', 'c', 'b'])

    def test_missing_ids_to_keep(self):
        md = Metadata(pd.DataFrame(
            {'col1': [1, 2, 3], 'col2': ['foo', 'bar', 'baz']},
            index=pd.Index(['a', 'b', 'c'], name='id')))

        with self.assertRaisesRegex(ValueError,
                                    "IDs.*not present.*'d', 'id1'"):
            md.filter_ids({'b', 'id1', 'c', 'd'})


class TestFilterColumns(unittest.TestCase):
    def setUp(self):
        get_dummy_plugin()

        # This object can be reused in many of the tests because its columns
        # match various filtering criteria, allowing test cases to test
        # individual parameters or combinations of parameters.
        self.metadata = Metadata(pd.DataFrame(
            {'cat': ['foo', 'bar', np.nan, 'foo'],
             'num': [42, np.nan, -5.5, 42],
             'uniq-cat': ['foo', np.nan, 'bar', np.nan],
             'uniq-num': [np.nan, 9.9, np.nan, 42],
             'zvar-cat': ['foo', np.nan, 'foo', 'foo'],
             'zvar-num': [9.9, 9.9, np.nan, 9.9],
             'empty-cat': np.array([np.nan, np.nan, np.nan, np.nan],
                                   dtype=object),
             'empty-num': [np.nan, np.nan, np.nan, np.nan]},
            index=pd.Index(['a', 'b', 'c', 'd'], name='id')))

        # Basic sanity check to ensure column types are what we expect them to
        # be.
        obs = {n: p.type for n, p in self.metadata.columns.items()}

        exp = {'cat': 'categorical',
               'num': 'numeric',
               'uniq-cat': 'categorical',
               'uniq-num': 'numeric',
               'zvar-cat': 'categorical',
               'zvar-num': 'numeric',
               'empty-cat': 'categorical',
               'empty-num': 'numeric'}

        self.assertEqual(obs, exp)

    def test_unknown_column_type(self):
        with self.assertRaisesRegex(
                ValueError, "Unknown column type 'foo'.*categorical, numeric"):
            self.metadata.filter_columns(column_type='foo')

    def test_no_filters(self):
        obs = self.metadata.filter_columns()

        self.assertEqual(obs, self.metadata)
        self.assertIsNot(obs, self.metadata)

    def test_all_filters_no_columns(self):
        md = Metadata(pd.DataFrame(
            {}, index=pd.Index(['a', 'b', 'c'], name='id')))

        obs = md.filter_columns(
            column_type='categorical', drop_all_unique=True,
            drop_zero_variance=True, drop_all_missing=True)

        self.assertEqual(obs, md)
        self.assertIsNot(obs, md)

        obs = md.filter_columns(
            column_type='numeric', drop_all_unique=True,
            drop_zero_variance=True, drop_all_missing=True)

        self.assertEqual(obs, md)
        self.assertIsNot(obs, md)

    def test_all_filters(self):
        obs = self.metadata.filter_columns(
            column_type='categorical', drop_all_unique=True,
            drop_zero_variance=True, drop_all_missing=True)

        self.assertEqual(set(obs.columns), {'cat'})

        obs = self.metadata.filter_columns(
            column_type='numeric', drop_all_unique=True,
            drop_zero_variance=True, drop_all_missing=True)

        self.assertEqual(set(obs.columns), {'num'})

    def test_all_columns_filtered(self):
        categorical = self.metadata.filter_columns(column_type='categorical')

        obs = categorical.filter_columns(column_type='numeric')

        exp = Metadata(pd.DataFrame(
            {}, index=pd.Index(['a', 'b', 'c', 'd'], name='id')))

        self.assertEqual(obs, exp)

    def test_filter_to_categorical(self):
        obs = self.metadata.filter_columns(column_type='categorical')

        self.assertEqual(set(obs.columns),
                         {'cat', 'uniq-cat', 'zvar-cat', 'empty-cat'})

    def test_filter_to_numeric(self):
        obs = self.metadata.filter_columns(column_type='numeric')

        self.assertEqual(set(obs.columns),
                         {'num', 'uniq-num', 'zvar-num', 'empty-num'})

    def test_drop_all_unique(self):
        obs = self.metadata.filter_columns(drop_all_unique=True)

        self.assertEqual(set(obs.columns),
                         {'cat', 'num', 'zvar-cat', 'zvar-num'})

    def test_drop_zero_variance(self):
        obs = self.metadata.filter_columns(drop_zero_variance=True)

        self.assertEqual(set(obs.columns),
                         {'cat', 'num', 'uniq-cat', 'uniq-num'})

    def test_drop_all_missing(self):
        obs = self.metadata.filter_columns(drop_all_missing=True)

        self.assertEqual(
            set(obs.columns),
            {'cat', 'num', 'uniq-cat', 'uniq-num', 'zvar-cat', 'zvar-num'})

    def test_drop_all_unique_with_single_id(self):
        md = Metadata(pd.DataFrame(
            {'cat': ['foo'],
             'num': [-4.2],
             'empty-cat': np.array([np.nan], dtype=object),
             'empty-num': [np.nan]},
            index=pd.Index(['id1'], name='id')))

        obs = md.filter_columns(drop_all_unique=True)

        exp = Metadata(pd.DataFrame({}, index=pd.Index(['id1'], name='id')))

        self.assertEqual(obs, exp)

    def test_drop_zero_variance_with_single_id(self):
        md = Metadata(pd.DataFrame(
            {'cat': ['foo'],
             'num': [-4.2],
             'empty-cat': np.array([np.nan], dtype=object),
             'empty-num': [np.nan]},
            index=pd.Index(['id1'], name='id')))

        obs = md.filter_columns(drop_zero_variance=True)

        exp = Metadata(pd.DataFrame({}, index=pd.Index(['id1'], name='id')))

        self.assertEqual(obs, exp)

    def test_retains_column_order(self):
        data = [[42, 'foo', 2.5],
                [42, 'bar', 0.5],
                [11, 'foo', 0.0]]
        md = Metadata(pd.DataFrame(
            data, columns=['z', 'a', 'ch'],
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))

        obs = md.filter_columns(column_type='numeric')

        exp_data = [[42, 2.5],
                    [42, 0.5],
                    [11, 0.0]]
        exp = Metadata(pd.DataFrame(
            exp_data, columns=['z', 'ch'],
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))

        self.assertEqual(obs, exp)

    def test_alternate_id_header(self):
        md = Metadata(pd.DataFrame(
            {'col1': ['foo', 'bar'],
             'col2': [-4.2, -4.2],
             'col3': ['bar', 'baz']},
            index=pd.Index(['id1', 'id2'], name='feature-id')))

        obs = md.filter_columns(drop_zero_variance=True)

        exp = Metadata(pd.DataFrame(
            {'col1': ['foo', 'bar'],
             'col3': ['bar', 'baz']},
            index=pd.Index(['id1', 'id2'], name='feature-id')))

        self.assertEqual(obs, exp)

    def test_no_artifacts(self):
        md = Metadata(pd.DataFrame(
            {'col1': [1, 2, 3], 'col2': ['foo', 'bar', 'baz']},
            index=pd.Index(['a', 'b', 'c'], name='id')))

        self.assertEqual(md.artifacts, ())

        filtered = md.filter_columns(column_type='categorical')

        self.assertEqual(filtered.artifacts, ())


if __name__ == '__main__':
    unittest.main()
