# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os.path
import tempfile
import unittest

import pandas as pd
import numpy as np

from qiime2 import Artifact
from qiime2.metadata import (MetadataColumn, CategoricalMetadataColumn,
                             NumericMetadataColumn)
from qiime2.core.testing.util import get_dummy_plugin, ReallyEqualMixin


# Dummy class for testing MetadataColumn ABC
class DummyMetadataColumn(MetadataColumn):
    type = 'dummy'

    @classmethod
    def _is_supported_dtype(cls, dtype):
        return dtype == 'float' or dtype == 'int'

    @classmethod
    def _normalize_(cls, series):
        return series.astype(float, copy=True, errors='raise')


class TestInvalidMetadataColumnConstruction(unittest.TestCase):
    def test_non_series(self):
        with self.assertRaisesRegex(
                TypeError, 'DummyMetadataColumn constructor.*Series.*not.*'
                           'DataFrame'):
            DummyMetadataColumn(pd.DataFrame(
                {'col1': [1, 2, 3]},
                index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_no_ids(self):
        with self.assertRaisesRegex(ValueError,
                                    'DummyMetadataColumn.*at least one ID'):
            DummyMetadataColumn(pd.Series([], name='col',
                                          index=pd.Index([], name='id'),
                                          dtype=object))

    def test_invalid_id_header(self):
        # default index name
        with self.assertRaisesRegex(ValueError, r'Index\.name.*None'):
            DummyMetadataColumn(pd.Series([1, 2, 3], name='col',
                                          index=pd.Index(['a', 'b', 'c'],
                                          dtype=object)))

        with self.assertRaisesRegex(ValueError, r'Index\.name.*my-id-header'):
            DummyMetadataColumn(pd.Series(
                [1, 2, 3], name='col',
                index=pd.Index(['a', 'b', 'c'], name='my-id-header')))

    def test_non_str_id(self):
        with self.assertRaisesRegex(
                TypeError, 'non-string metadata ID.*type.*float.*nan'):
            DummyMetadataColumn(pd.Series(
                [1, 2, 3], name='col',
                index=pd.Index(['a', np.nan, 'c'], name='id')))

    def test_non_str_column_name(self):
        # default series name
        with self.assertRaisesRegex(
                TypeError, 'non-string metadata column name.*type.*'
                           'NoneType.*None'):
            DummyMetadataColumn(pd.Series(
                [1, 2, 3], index=pd.Index(['a', 'b', 'c'], name='id')))

        with self.assertRaisesRegex(
                TypeError, 'non-string metadata column name.*type.*'
                           'float.*nan'):
            DummyMetadataColumn(pd.Series(
                [1, 2, 3], name=np.nan,
                index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_empty_id(self):
        with self.assertRaisesRegex(
                ValueError, 'empty metadata ID.*at least one character'):
            DummyMetadataColumn(pd.Series(
                [1, 2, 3], name='col',
                index=pd.Index(['a', '', 'c'], name='id')))

    def test_empty_column_name(self):
        with self.assertRaisesRegex(
                ValueError, 'empty metadata column name.*'
                            'at least one character'):
            DummyMetadataColumn(pd.Series(
                [1, 2, 3], name='',
                index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_pound_sign_id(self):
        with self.assertRaisesRegex(
                ValueError, "metadata ID.*begins with a pound sign.*'#b'"):
            DummyMetadataColumn(pd.Series(
                [1, 2, 3], name='col',
                index=pd.Index(['a', '#b', 'c'], name='id')))

    def test_id_conflicts_with_id_header(self):
        with self.assertRaisesRegex(
                ValueError, "metadata ID 'sample-id'.*conflicts.*reserved.*"
                            "ID header"):
            DummyMetadataColumn(pd.Series(
                [1, 2, 3], name='col',
                index=pd.Index(['a', 'sample-id', 'c'], name='id')))

    def test_column_name_conflicts_with_id_header(self):
        with self.assertRaisesRegex(
                ValueError, "metadata column name 'featureid'.*conflicts.*"
                            "reserved.*ID header"):
            DummyMetadataColumn(pd.Series(
                [1, 2, 3], name='featureid',
                index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_duplicate_ids(self):
        with self.assertRaisesRegex(ValueError, "Metadata IDs.*unique.*'a'"):
            DummyMetadataColumn(pd.Series(
                [1, 2, 3], name='col',
                index=pd.Index(['a', 'b', 'a'], name='id')))

    def test_unsupported_column_dtype(self):
        with self.assertRaisesRegex(
                TypeError, "DummyMetadataColumn 'col1' does not support.*"
                           "Series.*dtype.*bool"):
            DummyMetadataColumn(pd.Series(
                [True, False, True], name='col1',
                index=pd.Index(['a', 'b', 'c'], name='id')))


class TestMetadataColumnConstructionAndProperties(unittest.TestCase):
    def test_single_id(self):
        index = pd.Index(['id1'], name='id')
        series = pd.Series([42], name='col1', index=index)
        mdc = DummyMetadataColumn(series)

        self.assertEqual(mdc.id_count, 1)
        self.assertEqual(mdc.id_header, 'id')
        self.assertEqual(mdc.ids, ('id1',))
        self.assertEqual(mdc.name, 'col1')

    def test_multiple_ids(self):
        index = pd.Index(['id1', 'a', 'my-id'], name='id')
        series = pd.Series([42, 4.2, -4.2], name='column', index=index)
        mdc = DummyMetadataColumn(series)

        self.assertEqual(mdc.id_count, 3)
        self.assertEqual(mdc.id_header, 'id')
        self.assertEqual(mdc.ids, ('id1', 'a', 'my-id'))
        self.assertEqual(mdc.name, 'column')

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
            series = pd.Series([0, 123], name='column', index=index)
            mdc = DummyMetadataColumn(series)

            self.assertEqual(mdc.id_header, header)
            count += 1

        # Since this test case is a little complicated, make sure that the
        # expected number of comparisons are happening.
        self.assertEqual(count, 26)

    def test_recommended_ids(self):
        index = pd.Index(['c6ca034a-223f-40b4-a0e0-45942912a5ea', 'My.ID'],
                         name='id')
        series = pd.Series([-1, -2], name='col1', index=index)
        mdc = DummyMetadataColumn(series)

        self.assertEqual(mdc.id_count, 2)
        self.assertEqual(mdc.id_header, 'id')
        self.assertEqual(mdc.ids,
                         ('c6ca034a-223f-40b4-a0e0-45942912a5ea', 'My.ID'))
        self.assertEqual(mdc.name, 'col1')

    def test_non_standard_characters(self):
        index = pd.Index(['©id##1', '((id))2', "'id_3<>'", '"id#4"',
                          'i d\r\t\n5'], name='id')
        series = pd.Series([0, 1, 2, 3, 4], name='↩c@l1™', index=index)
        mdc = DummyMetadataColumn(series)

        self.assertEqual(mdc.id_count, 5)
        self.assertEqual(mdc.id_header, 'id')
        self.assertEqual(
            mdc.ids, ('©id##1', '((id))2', "'id_3<>'", '"id#4"', 'i d\r\t\n5'))
        self.assertEqual(mdc.name, '↩c@l1™')

    def test_missing_data(self):
        index = pd.Index(['None', 'nan', 'NA'], name='id')
        series = pd.Series([np.nan, np.nan, np.nan], name='NA', index=index)
        mdc = DummyMetadataColumn(series)

        self.assertEqual(mdc.id_count, 3)
        self.assertEqual(mdc.id_header, 'id')
        self.assertEqual(mdc.ids, ('None', 'nan', 'NA'))
        self.assertEqual(mdc.name, 'NA')

    def test_does_not_cast_ids_or_column_name(self):
        index = pd.Index(['0.000001', '0.004000', '0.000000'], dtype=object,
                         name='id')
        series = pd.Series([2.0, 1.0, 3.0], name='42.0', index=index)
        mdc = DummyMetadataColumn(series)

        self.assertEqual(mdc.id_count, 3)
        self.assertEqual(mdc.id_header, 'id')
        self.assertEqual(mdc.ids, ('0.000001', '0.004000', '0.000000'))
        self.assertEqual(mdc.name, '42.0')

    def test_case_insensitive_duplicate_ids(self):
        index = pd.Index(['a', 'b', 'A'], name='id')
        series = pd.Series([1, 2, 3], name='column', index=index)
        mdc = DummyMetadataColumn(series)

        self.assertEqual(mdc.ids, ('a', 'b', 'A'))


class TestSourceArtifacts(unittest.TestCase):
    def setUp(self):
        self.mdc = DummyMetadataColumn(pd.Series(
            [1, 2, 3], name='col', index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_no_source_artifacts(self):
        self.assertEqual(self.mdc.artifacts, ())

    def test_add_zero_artifacts(self):
        self.mdc._add_artifacts([])

        self.assertEqual(self.mdc.artifacts, ())

    def test_add_artifacts(self):
        # First two artifacts have the same data but different UUIDs.
        artifact1 = Artifact.import_data('Mapping', {'a': '1', 'b': '3'})
        self.mdc._add_artifacts([artifact1])

        artifact2 = Artifact.import_data('Mapping', {'a': '1', 'b': '3'})
        artifact3 = Artifact.import_data('IntSequence1', [1, 2, 3, 4])
        self.mdc._add_artifacts([artifact2, artifact3])

        self.assertEqual(self.mdc.artifacts, (artifact1, artifact2, artifact3))

    def test_add_non_artifact(self):
        artifact = Artifact.import_data('Mapping', {'a': '1', 'b': '3'})

        with self.assertRaisesRegex(TypeError, "Artifact object.*42"):
            self.mdc._add_artifacts([artifact, 42])

        # Test that the object hasn't been mutated.
        self.assertEqual(self.mdc.artifacts, ())

    def test_add_duplicate_artifact(self):
        artifact1 = Artifact.import_data('Mapping', {'a': '1', 'b': '3'})
        artifact2 = Artifact.import_data('IntSequence1', [1, 2, 3, 4])
        self.mdc._add_artifacts([artifact1, artifact2])

        with self.assertRaisesRegex(
                ValueError, "Duplicate source artifacts.*DummyMetadataColumn.*"
                            "artifact: Mapping"):
            self.mdc._add_artifacts([artifact1])

        # Test that the object hasn't been mutated.
        self.assertEqual(self.mdc.artifacts, (artifact1, artifact2))


class TestRepr(unittest.TestCase):
    def test_single_id(self):
        mdc = DummyMetadataColumn(pd.Series(
            [42], name='foo', index=pd.Index(['id1'], name='id')))

        obs = repr(mdc)

        self.assertEqual(obs, "<DummyMetadataColumn name='foo' id_count=1>")

    def test_multiple_ids(self):
        mdc = DummyMetadataColumn(pd.Series(
            [42, 43, 44], name='my column',
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))

        obs = repr(mdc)

        self.assertEqual(
            obs, "<DummyMetadataColumn name='my column' id_count=3>")


class TestEqualityOperators(unittest.TestCase, ReallyEqualMixin):
    def setUp(self):
        get_dummy_plugin()

    def test_type_mismatch(self):
        dummy = DummyMetadataColumn(pd.Series(
            [1.0, 2.0, 3.0], name='col1',
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))
        numeric = NumericMetadataColumn(pd.Series(
            [1.0, 2.0, 3.0], name='col1',
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))
        categorical = CategoricalMetadataColumn(pd.Series(
            ['a', 'b', 'c'], name='col1',
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))

        self.assertReallyNotEqual(dummy, numeric)
        self.assertReallyNotEqual(dummy, categorical)

    def test_id_header_mismatch(self):
        mdc1 = DummyMetadataColumn(pd.Series(
            [42, 43], name='col1', index=pd.Index(['id1', 'id2'], name='id')))
        mdc2 = DummyMetadataColumn(pd.Series(
            [42, 43], name='col1', index=pd.Index(['id1', 'id2'], name='ID')))

        self.assertReallyNotEqual(mdc1, mdc2)

    def test_artifacts_mismatch(self):
        artifact1 = Artifact.import_data('Mapping', {'a': '1', 'b': '2'})
        artifact2 = Artifact.import_data('Mapping', {'a': '1', 'b': '2'})
        series = pd.Series([42, 43], name='col1',
                           index=pd.Index(['id1', 'id2'], name='id'))

        # No artifacts
        mdc1 = DummyMetadataColumn(series)

        # Has an artifact
        mdc2 = DummyMetadataColumn(series)
        mdc2._add_artifacts([artifact1])

        # Has a different artifact
        mdc3 = DummyMetadataColumn(series)
        mdc3._add_artifacts([artifact2])

        self.assertReallyNotEqual(mdc1, mdc2)
        self.assertReallyNotEqual(mdc2, mdc3)

    def test_id_mismatch(self):
        mdc1 = DummyMetadataColumn(pd.Series(
            [42, 43], name='col1', index=pd.Index(['id1', 'id2'], name='id')))
        mdc2 = DummyMetadataColumn(pd.Series(
            [42, 43], name='col1', index=pd.Index(['id1', 'id3'], name='id')))

        self.assertReallyNotEqual(mdc1, mdc2)

    def test_column_name_mismatch(self):
        mdc1 = DummyMetadataColumn(pd.Series(
            [42, 43], name='col1', index=pd.Index(['id1', 'id2'], name='id')))
        mdc2 = DummyMetadataColumn(pd.Series(
            [42, 43], name='col2', index=pd.Index(['id1', 'id2'], name='id')))

        self.assertReallyNotEqual(mdc1, mdc2)

    def test_data_mismatch(self):
        mdc1 = DummyMetadataColumn(pd.Series(
            [42, 43], name='col1', index=pd.Index(['id1', 'id2'], name='id')))
        mdc2 = DummyMetadataColumn(pd.Series(
            [42, 42], name='col1', index=pd.Index(['id1', 'id2'], name='id')))

        self.assertReallyNotEqual(mdc1, mdc2)

    def test_equality_without_artifact(self):
        mdc1 = DummyMetadataColumn(pd.Series(
            [42, 43], name='col1', index=pd.Index(['id1', 'id2'], name='id')))
        mdc2 = DummyMetadataColumn(pd.Series(
            [42, 43], name='col1', index=pd.Index(['id1', 'id2'], name='id')))

        self.assertReallyEqual(mdc1, mdc2)

    def test_equality_with_artifact(self):
        artifact = Artifact.import_data('Mapping', {'a': '1', 'b': '2'})

        mdc1 = DummyMetadataColumn(pd.Series(
            [42, 43], name='col1', index=pd.Index(['id1', 'id2'], name='id')))
        mdc1._add_artifacts([artifact])

        mdc2 = DummyMetadataColumn(pd.Series(
            [42, 43], name='col1', index=pd.Index(['id1', 'id2'], name='id')))
        mdc2._add_artifacts([artifact])

        self.assertReallyEqual(mdc1, mdc2)

    def test_equality_with_missing_data(self):
        mdc1 = DummyMetadataColumn(pd.Series(
            [42, np.nan, 43, np.nan], name='col1',
            index=pd.Index(['id1', 'id2', 'id3', 'id4'], name='id')))
        mdc2 = DummyMetadataColumn(pd.Series(
            [42, np.nan, 43, np.nan], name='col1',
            index=pd.Index(['id1', 'id2', 'id3', 'id4'], name='id')))

        self.assertReallyEqual(mdc1, mdc2)


# Extensive tests of the MetadataWriter are performed in test_io.py. This test
# is a sanity check that a new MetadataColumn subclass (DummyMetadataColumn)
# can be written to disk with its column type preserved. This test would have
# caught a bug in the original implementation of MetadataColumn.save(), which
# converted itself into a Metadata object, losing the "dummy" column type and
# replacing it with "numeric". In order for a MetadataColumn to turn itself
# into a Metadata object in a lossless/safe way, the Metadata constructor needs
# a `column_types` parameter to preserve column types.
class TestSave(unittest.TestCase):
    def setUp(self):
        self.temp_dir_obj = tempfile.TemporaryDirectory(
            prefix='qiime2-metadata-tests-temp-')
        self.temp_dir = self.temp_dir_obj.name

        self.filepath = os.path.join(self.temp_dir, 'metadata.tsv')

    def tearDown(self):
        self.temp_dir_obj.cleanup()

    def test_basic(self):
        mdc = DummyMetadataColumn(pd.Series(
            [42, 42.5, -999.123], name='dummy-column',
            index=pd.Index(['id1', 'id2', 'id3'], name='id')))

        mdc.save(self.filepath)

        with open(self.filepath, 'r') as fh:
            obs = fh.read()

        exp = (
            "id\tdummy-column\n"
            "#q2:types\tdummy\n"
            "id1\t42\n"
            "id2\t42.5\n"
            "id3\t-999.123\n"
        )

        self.assertEqual(obs, exp)


class TestToSeries(unittest.TestCase):
    def test_single_id(self):
        series = pd.Series([0.0], name='col',
                           index=pd.Index(['id1'], name='id'))
        mdc = DummyMetadataColumn(series)

        obs = mdc.to_series()

        pd.testing.assert_series_equal(obs, series)

    def test_multiple_ids(self):
        series = pd.Series([-1.5, np.nan, 42], name='col',
                           index=pd.Index(['id1', 'id2', 'id3'], name='id'))
        mdc = DummyMetadataColumn(series)

        obs = mdc.to_series()

        pd.testing.assert_series_equal(obs, series)

    def test_id_header_preserved(self):
        series = pd.Series(
            [-1.5, 0.0, 42], name='col',
            index=pd.Index(['id1', 'id2', 'id3'], name='#OTU ID'))
        mdc = DummyMetadataColumn(series)

        obs = mdc.to_series()

        pd.testing.assert_series_equal(obs, series)
        self.assertEqual(obs.index.name, '#OTU ID')

    def test_series_copy(self):
        series = pd.Series([1, 2.5, 3], name='col',
                           index=pd.Index(['id1', 'id2', 'id3'], name='id'))
        mdc = DummyMetadataColumn(series)

        obs = mdc.to_series()

        pd.testing.assert_series_equal(obs, series)
        self.assertIsNot(obs, series)


class TestToDataframe(unittest.TestCase):
    def test_single_id(self):
        series = pd.Series([0.0], name='col',
                           index=pd.Index(['id1'], name='id'))
        mdc = DummyMetadataColumn(series)

        obs = mdc.to_dataframe()

        exp = pd.DataFrame({'col': [0.0]}, index=pd.Index(['id1'], name='id'))

        pd.testing.assert_frame_equal(obs, exp)

    def test_multiple_ids(self):
        series = pd.Series([0.0, 4.2, np.nan], name='my column',
                           index=pd.Index(['a', 'b', 'c'], name='id'))
        mdc = DummyMetadataColumn(series)

        obs = mdc.to_dataframe()

        exp = pd.DataFrame({'my column': [0.0, 4.2, np.nan]},
                           index=pd.Index(['a', 'b', 'c'], name='id'))

        pd.testing.assert_frame_equal(obs, exp)

    def test_id_header_preserved(self):
        series = pd.Series([0.0, 4.2, 123], name='my column',
                           index=pd.Index(['a', 'b', 'c'], name='#Sample ID'))
        mdc = DummyMetadataColumn(series)

        obs = mdc.to_dataframe()

        exp = pd.DataFrame({'my column': [0.0, 4.2, 123]},
                           index=pd.Index(['a', 'b', 'c'], name='#Sample ID'))

        pd.testing.assert_frame_equal(obs, exp)
        self.assertEqual(obs.index.name, '#Sample ID')


class TestGetValue(unittest.TestCase):
    def test_id_not_found(self):
        series = pd.Series([1, 2, 3], name='col1',
                           index=pd.Index(['a', 'b', 'c'], name='id'))
        mdc = DummyMetadataColumn(series)

        with self.assertRaisesRegex(
                ValueError, "'d' is not present.*DummyMetadataColumn.*'col1'"):
            mdc.get_value('d')

    def test_get_value(self):
        series = pd.Series([1, 2, np.nan], name='col1',
                           index=pd.Index(['a', 'b', 'c'], name='id'))
        mdc = DummyMetadataColumn(series)

        obs = mdc.get_value('a')

        self.assertEqual(obs, 1.0)

        obs = mdc.get_value('b')

        self.assertEqual(obs, 2.0)

        obs = mdc.get_value('c')

        self.assertTrue(np.isnan(obs))


class TestHasMissingValues(unittest.TestCase):
    def test_no_missing_values(self):
        series = pd.Series([0.0, 2.2, 3.3], name='col1',
                           index=pd.Index(['a', 'b', 'c'], name='id'))
        mdc = DummyMetadataColumn(series)

        obs = mdc.has_missing_values()

        self.assertEqual(obs, False)

    def test_with_missing_values(self):
        series = pd.Series([0.0, np.nan, 3.3], name='col1',
                           index=pd.Index(['a', 'b', 'c'], name='id'))
        mdc = DummyMetadataColumn(series)

        obs = mdc.has_missing_values()

        self.assertEqual(obs, True)


class TestDropMissingValues(unittest.TestCase):
    def test_no_missing_values(self):
        series = pd.Series([0.0, 2.2, 3.3], name='col1',
                           index=pd.Index(['a', 'b', 'c'], name='id'))
        mdc = DummyMetadataColumn(series)

        obs = mdc.drop_missing_values()

        self.assertEqual(obs, mdc)
        self.assertIsNot(obs, mdc)

    def test_with_missing_values(self):
        series = pd.Series(
            [0.0, np.nan, 3.3, np.nan, np.nan, 4.4], name='col1',
            index=pd.Index(['a', 'b', 'c', 'd', 'e', 'f'], name='sampleid'))
        mdc = DummyMetadataColumn(series)

        obs = mdc.drop_missing_values()

        exp = DummyMetadataColumn(pd.Series(
            [0.0, 3.3, 4.4], name='col1',
            index=pd.Index(['a', 'c', 'f'], name='sampleid')))

        self.assertEqual(obs, exp)

    def test_artifacts_are_propagated(self):
        artifact = Artifact.import_data('Mapping', {'a': '1', 'b': '2'})

        series = pd.Series(
            [0.0, np.nan, 3.3, np.nan, np.nan, 4.4], name='col1',
            index=pd.Index(['a', 'b', 'c', 'd', 'e', 'f'], name='sampleid'))
        mdc = DummyMetadataColumn(series)
        mdc._add_artifacts([artifact])

        obs = mdc.drop_missing_values()

        exp = DummyMetadataColumn(pd.Series(
            [0.0, 3.3, 4.4], name='col1',
            index=pd.Index(['a', 'c', 'f'], name='sampleid')))
        exp._add_artifacts([artifact])

        self.assertEqual(obs, exp)
        self.assertEqual(obs.artifacts, (artifact,))


class TestGetIDs(unittest.TestCase):
    def test_single_id(self):
        series = pd.Series([1.234], name='col1',
                           index=pd.Index(['my id'], name='id'))
        mdc = DummyMetadataColumn(series)

        obs = mdc.get_ids()

        self.assertEqual(obs, {'my id'})

    def test_multiple_ids(self):
        series = pd.Series(
            [1.234, np.nan, 5.67, np.nan, 8.9], name='col1',
            index=pd.Index(['id1', 'id2', 'id3', 'id4', 'id5'], name='id'))
        mdc = DummyMetadataColumn(series)

        obs = mdc.get_ids()

        self.assertEqual(obs, {'id1', 'id2', 'id3', 'id4', 'id5'})

    def test_where_values_missing(self):
        series = pd.Series(
            [1.234, np.nan, 5.67, np.nan, 8.9], name='col1',
            index=pd.Index(['id1', 'id2', 'id3', 'id4', 'id5'], name='id'))
        mdc = DummyMetadataColumn(series)

        obs = mdc.get_ids(where_values_missing=True)

        self.assertEqual(obs, {'id2', 'id4'})

    def test_where_values_missing_all_missing(self):
        series = pd.Series(
            [np.nan, np.nan, np.nan], name='col1',
            index=pd.Index(['id1', 'id2', 'id3'], name='id'))
        mdc = DummyMetadataColumn(series)

        obs = mdc.get_ids(where_values_missing=True)

        self.assertEqual(obs, {'id1', 'id2', 'id3'})


class TestFilterIDs(unittest.TestCase):
    def setUp(self):
        get_dummy_plugin()

    def test_supports_iterable(self):
        mdc = DummyMetadataColumn(pd.Series(
            [1, 2, 3], name='col1',
            index=pd.Index(['a', 'b', 'c'], name='id')))

        obs = mdc.filter_ids(iter({'a', 'c'}))

        exp = DummyMetadataColumn(pd.Series(
            [1, 3], name='col1', index=pd.Index(['a', 'c'], name='id')))

        self.assertEqual(obs, exp)

    def test_keep_all(self):
        mdc = DummyMetadataColumn(pd.Series(
            [1, 2, 3], name='col1',
            index=pd.Index(['a', 'b', 'c'], name='id')))

        obs = mdc.filter_ids({'a', 'b', 'c'})

        self.assertEqual(obs, mdc)
        self.assertIsNot(obs, mdc)

    def test_keep_multiple(self):
        mdc = DummyMetadataColumn(pd.Series(
            [1, 2, 3], name='col1',
            index=pd.Index(['a', 'b', 'c'], name='id')))

        obs = mdc.filter_ids({'a', 'c'})

        exp = DummyMetadataColumn(pd.Series(
            [1, 3], name='col1',
            index=pd.Index(['a', 'c'], name='id')))

        self.assertEqual(obs, exp)

    def test_keep_one(self):
        mdc = DummyMetadataColumn(pd.Series(
            [1, 2, 3], name='col1',
            index=pd.Index(['a', 'b', 'c'], name='id')))

        obs = mdc.filter_ids({'b'})

        exp = DummyMetadataColumn(pd.Series(
            [2], name='col1', index=pd.Index(['b'], name='id')))

        self.assertEqual(obs, exp)

    def test_alternate_id_header(self):
        mdc = DummyMetadataColumn(pd.Series(
            [1, 2, 3], name='col1',
            index=pd.Index(['a', 'b', 'c'], name='#OTU ID')))

        obs = mdc.filter_ids({'b', 'c'})

        exp = DummyMetadataColumn(pd.Series(
            [2, 3], name='col1', index=pd.Index(['b', 'c'], name='#OTU ID')))

        self.assertEqual(obs, exp)

    def test_no_artifacts(self):
        mdc = DummyMetadataColumn(pd.Series(
            [1, 2, 3], name='col1',
            index=pd.Index(['a', 'b', 'c'], name='id')))

        self.assertEqual(mdc.artifacts, ())

        filtered = mdc.filter_ids({'b'})

        self.assertEqual(filtered.artifacts, ())

    def test_with_artifacts(self):
        artifact1 = Artifact.import_data('Mapping', {'a': '1', 'b': '2'})
        artifact2 = Artifact.import_data('Mapping', {'d': '4'})

        mdc = DummyMetadataColumn(pd.Series(
            [1, 2, 3], name='col1',
            index=pd.Index(['a', 'b', 'c'], name='id')))
        mdc._add_artifacts([artifact1, artifact2])

        obs = mdc.filter_ids({'a', 'c'})

        exp = DummyMetadataColumn(pd.Series(
            [1, 3], name='col1', index=pd.Index(['a', 'c'], name='id')))
        exp._add_artifacts([artifact1, artifact2])

        self.assertEqual(obs, exp)
        self.assertEqual(obs.artifacts, (artifact1, artifact2))

    def test_empty_ids_to_keep(self):
        mdc = DummyMetadataColumn(pd.Series(
            [1, 2, 3], name='col1',
            index=pd.Index(['a', 'b', 'c'], name='id')))

        with self.assertRaisesRegex(ValueError,
                                    'ids_to_keep.*at least one ID'):
            mdc.filter_ids({})

    def test_duplicate_ids_to_keep(self):
        mdc = DummyMetadataColumn(pd.Series(
            [1, 2, 3], name='col1',
            index=pd.Index(['a', 'b', 'c'], name='id')))

        with self.assertRaisesRegex(ValueError,
                                    "ids_to_keep.*unique IDs.*'b'"):
            mdc.filter_ids(['b', 'c', 'b'])

    def test_missing_ids_to_keep(self):
        mdc = DummyMetadataColumn(pd.Series(
            [1, 2, 3], name='col1',
            index=pd.Index(['a', 'b', 'c'], name='id')))

        with self.assertRaisesRegex(ValueError,
                                    "IDs.*not present.*'d', 'id1'"):
            mdc.filter_ids({'b', 'id1', 'c', 'd'})


# The tests for CategoricalMetadataColumn and NumericMetadataColumn only test
# behavior specific to these subclasses. More extensive tests of these objects
# are performed above by testing the MetadataColumn ABC in a generic way.
class TestCategoricalMetadataColumn(unittest.TestCase):
    def test_unsupported_dtype(self):
        with self.assertRaisesRegex(
                TypeError, "CategoricalMetadataColumn 'col1' does not support"
                           ".*Series.*dtype.*float64"):
            CategoricalMetadataColumn(pd.Series(
                [42.5, 42.6, 42.7], name='col1',
                index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_unsupported_type_value(self):
        with self.assertRaisesRegex(
                TypeError, "CategoricalMetadataColumn.*strings or missing "
                           r"values.*42\.5.*float.*'col1'"):
            CategoricalMetadataColumn(pd.Series(
                ['foo', 'bar', 42.5], name='col1',
                index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_empty_str_value(self):
        with self.assertRaisesRegex(
                ValueError, "CategoricalMetadataColumn.*empty strings.*"
                            "column 'col1'"):
            CategoricalMetadataColumn(pd.Series(
                ['foo', '', 'bar'], name='col1',
                index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_type_property(self):
        self.assertEqual(CategoricalMetadataColumn.type, 'categorical')

    def test_supported_dtype(self):
        series = pd.Series(
            ['foo', np.nan, 'bar', 'foo'], name='my column',
            index=pd.Index(['a', 'b', 'c', 'd'], name='id'))
        mdc = CategoricalMetadataColumn(series)

        self.assertEqual(mdc.id_count, 4)
        self.assertEqual(mdc.id_header, 'id')
        self.assertEqual(mdc.ids, ('a', 'b', 'c', 'd'))
        self.assertEqual(mdc.name, 'my column')

        obs_series = mdc.to_series()
        pd.testing.assert_series_equal(obs_series, series)
        self.assertEqual(obs_series.dtype, object)

    def test_numeric_strings_preserved_as_strings(self):
        series = pd.Series(
            ['1', np.nan, '2.5', '3.0'], name='my column',
            index=pd.Index(['a', 'b', 'c', 'd'], name='id'))
        mdc = CategoricalMetadataColumn(series)

        self.assertEqual(mdc.id_count, 4)
        self.assertEqual(mdc.id_header, 'id')
        self.assertEqual(mdc.ids, ('a', 'b', 'c', 'd'))
        self.assertEqual(mdc.name, 'my column')

        obs_series = mdc.to_series()
        pd.testing.assert_series_equal(obs_series, series)
        self.assertEqual(obs_series.dtype, object)

    def test_missing_data_normalized(self):
        # Different missing data representations should be normalized to np.nan
        mdc = CategoricalMetadataColumn(pd.Series(
            [np.nan, 'foo', float('nan'), None], name='col1',
            index=pd.Index(['a', 'b', 'c', 'd'], name='id')))

        obs = mdc.to_series()

        exp = pd.Series(
            [np.nan, 'foo', np.nan, np.nan], name='col1',
            index=pd.Index(['a', 'b', 'c', 'd'], name='id'))

        pd.testing.assert_series_equal(obs, exp)
        self.assertEqual(obs.dtype, object)
        self.assertTrue(np.isnan(obs['a']))
        self.assertTrue(np.isnan(obs['c']))
        self.assertTrue(np.isnan(obs['d']))

    def test_all_missing_data(self):
        mdc = CategoricalMetadataColumn(pd.Series(
            np.array([np.nan, np.nan, np.nan], dtype=object), name='col1',
            index=pd.Index(['a', 'b', 'c'], name='id')))

        obs = mdc.to_series()

        exp = pd.Series(
            np.array([np.nan, np.nan, np.nan], dtype=object), name='col1',
            index=pd.Index(['a', 'b', 'c'], name='id'))

        pd.testing.assert_series_equal(obs, exp)
        self.assertEqual(obs.dtype, object)

    def test_leading_trailing_whitespace_value(self):
        col1 = CategoricalMetadataColumn(pd.Series(
            ['foo', ' bar ', 'baz'], name='col1',
            index=pd.Index(['a', 'b', 'c'], name='id')))
        col2 = CategoricalMetadataColumn(pd.Series(
            ['foo', 'bar', 'baz'], name='col1',
            index=pd.Index(['a', 'b', 'c'], name='id')))

        self.assertEqual(col1, col2)

    def test_leading_trailing_whitespace_id(self):
        col1 = CategoricalMetadataColumn(pd.Series(
                ['foo', ' bar ', 'baz'], name='col',
                index=pd.Index(['a', ' b ', 'c'], name='id')))
        col2 = CategoricalMetadataColumn(pd.Series(
                ['foo', ' bar ', 'baz'], name='col',
                index=pd.Index(['a', 'b', 'c'], name='id')))

        self.assertEqual(col1, col2)

    def test_leading_trailing_whitespace_column_name(self):
        col1 = CategoricalMetadataColumn(pd.Series(
                ['foo', ' bar ', 'baz'], name=' col ',
                index=pd.Index(['a', 'b', 'c'], name='id')))
        col2 = CategoricalMetadataColumn(pd.Series(
                ['foo', ' bar ', 'baz'], name='col',
                index=pd.Index(['a', 'b', 'c'], name='id')))

        self.assertEqual(col1, col2)


class TestNumericMetadataColumn(unittest.TestCase):
    def test_unsupported_dtype(self):
        with self.assertRaisesRegex(
                TypeError, "NumericMetadataColumn 'col1' does not support"
                           ".*Series.*dtype.*bool"):
            NumericMetadataColumn(pd.Series(
                [True, False, True], name='col1',
                index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_infinity_value(self):
        with self.assertRaisesRegex(
                ValueError, "NumericMetadataColumn.*positive or negative "
                            "infinity.*column 'col1'"):
            NumericMetadataColumn(pd.Series(
                [42, float('+inf'), 4.3], name='col1',
                index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_type_property(self):
        self.assertEqual(NumericMetadataColumn.type, 'numeric')

    def test_supported_dtype_float(self):
        series = pd.Series(
            [1.23, np.nan, 4.56, -7.891], name='my column',
            index=pd.Index(['a', 'b', 'c', 'd'], name='id'))
        mdc = NumericMetadataColumn(series)

        self.assertEqual(mdc.id_count, 4)
        self.assertEqual(mdc.id_header, 'id')
        self.assertEqual(mdc.ids, ('a', 'b', 'c', 'd'))
        self.assertEqual(mdc.name, 'my column')

        obs_series = mdc.to_series()
        pd.testing.assert_series_equal(obs_series, series)
        self.assertEqual(obs_series.dtype, np.float64)

    def test_supported_dtype_int(self):
        series = pd.Series(
            [0, 1, 42, -2], name='my column',
            index=pd.Index(['a', 'b', 'c', 'd'], name='id'))
        mdc = NumericMetadataColumn(series)

        self.assertEqual(mdc.id_count, 4)
        self.assertEqual(mdc.id_header, 'id')
        self.assertEqual(mdc.ids, ('a', 'b', 'c', 'd'))
        self.assertEqual(mdc.name, 'my column')

        obs_series = mdc.to_series()

        exp_series = pd.Series(
            [0.0, 1.0, 42.0, -2.0], name='my column',
            index=pd.Index(['a', 'b', 'c', 'd'], name='id'))

        pd.testing.assert_series_equal(obs_series, exp_series)
        self.assertEqual(obs_series.dtype, np.float64)

    def test_missing_data_normalized(self):
        # Different missing data representations should be normalized to np.nan
        mdc = NumericMetadataColumn(pd.Series(
            [np.nan, 4.2, float('nan'), -5.678], name='col1',
            index=pd.Index(['a', 'b', 'c', 'd'], name='id')))

        obs = mdc.to_series()

        exp = pd.Series(
            [np.nan, 4.2, np.nan, -5.678], name='col1',
            index=pd.Index(['a', 'b', 'c', 'd'], name='id'))

        pd.testing.assert_series_equal(obs, exp)
        self.assertEqual(obs.dtype, np.float64)
        self.assertTrue(np.isnan(obs['a']))
        self.assertTrue(np.isnan(obs['c']))

    def test_all_missing_data(self):
        mdc = NumericMetadataColumn(pd.Series(
            [np.nan, np.nan, np.nan], name='col1',
            index=pd.Index(['a', 'b', 'c'], name='id')))

        obs = mdc.to_series()

        exp = pd.Series(
            [np.nan, np.nan, np.nan], name='col1',
            index=pd.Index(['a', 'b', 'c'], name='id'))

        pd.testing.assert_series_equal(obs, exp)
        self.assertEqual(obs.dtype, np.float64)


if __name__ == '__main__':
    unittest.main()
