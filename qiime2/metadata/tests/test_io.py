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

from qiime2.metadata import Metadata, MetadataFileError


class TestLoad(unittest.TestCase):
    def get_data_path(self, filename):
        return pkg_resources.resource_filename('qiime2.metadata.tests',
                                               'data/%s' % filename)

    def test_path_does_not_exist(self):
        with self.assertRaisesRegex(MetadataFileError,
                                    "Metadata file path doesn't exist"):
            Metadata.load(
                '/qiime2/unit/tests/hopefully/this/path/does/not/exist')

    def test_path_is_directory(self):
        fp = self.get_data_path('valid')

        with self.assertRaisesRegex(MetadataFileError,
                                    "path points to something other than a "
                                    "file"):
            Metadata.load(fp)

    def test_non_utf_8_file(self):
        fp = self.get_data_path('invalid/non-utf-8.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    'encoded as UTF-8 or ASCII'):
            Metadata.load(fp)

    def test_empty_file(self):
        fp = self.get_data_path('invalid/empty-file')

        with self.assertRaisesRegex(MetadataFileError,
                                    'locate header.*file may be empty'):
            Metadata.load(fp)

    def test_comments_and_empty_rows_only(self):
        fp = self.get_data_path('invalid/comments-and-empty-rows-only.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    'locate header.*only of comments or empty '
                                    'rows'):
            Metadata.load(fp)

    def test_header_only(self):
        fp = self.get_data_path('invalid/header-only.tsv')

        with self.assertRaisesRegex(MetadataFileError, 'at least one ID'):
            Metadata.load(fp)

    def test_header_only_with_comments_and_empty_rows(self):
        fp = self.get_data_path(
            'invalid/header-only-with-comments-and-empty-rows.tsv')

        with self.assertRaisesRegex(MetadataFileError, 'at least one ID'):
            Metadata.load(fp)

    def test_qiime1_empty_mapping_file(self):
        fp = self.get_data_path('invalid/qiime1-empty.tsv')

        with self.assertRaisesRegex(MetadataFileError, 'at least one ID'):
            Metadata.load(fp)

    def test_invalid_header(self):
        fp = self.get_data_path('invalid/invalid-header.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    'unrecognized ID column name.*'
                                    'invalid_id_header'):
            Metadata.load(fp)

    def test_empty_id(self):
        fp = self.get_data_path('invalid/empty-id.tsv')

        with self.assertRaisesRegex(MetadataFileError, 'empty metadata ID'):
            Metadata.load(fp)

    def test_whitespace_only_id(self):
        fp = self.get_data_path('invalid/whitespace-only-id.tsv')

        with self.assertRaisesRegex(MetadataFileError, 'empty metadata ID'):
            Metadata.load(fp)

    def test_empty_column_name(self):
        fp = self.get_data_path('invalid/empty-column-name.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    'column without a name'):
            Metadata.load(fp)

    def test_whitespace_only_column_name(self):
        fp = self.get_data_path('invalid/whitespace-only-column-name.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    'column without a name'):
            Metadata.load(fp)

    def test_duplicate_ids(self):
        fp = self.get_data_path('invalid/duplicate-ids.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    'IDs must be unique.*id1'):
            Metadata.load(fp)

    def test_duplicate_ids_with_whitespace(self):
        fp = self.get_data_path('invalid/duplicate-ids-with-whitespace.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    'IDs must be unique.*id1'):
            Metadata.load(fp)

    def test_duplicate_column_names(self):
        fp = self.get_data_path('invalid/duplicate-column-names.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    'Column names must be unique.*col1'):
            Metadata.load(fp)

    def test_duplicate_column_names_with_whitespace(self):
        fp = self.get_data_path(
            'invalid/duplicate-column-names-with-whitespace.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    'Column names must be unique.*col1'):
            Metadata.load(fp)

    def test_id_conflicts_with_id_header(self):
        fp = self.get_data_path('invalid/id-conflicts-with-id-header.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    "ID 'id' conflicts.*ID column header"):
            Metadata.load(fp)

    def test_column_name_conflicts_with_id_header(self):
        fp = self.get_data_path(
            'invalid/column-name-conflicts-with-id-header.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    "column name 'featureid' conflicts.*ID "
                                    "column header"):
            Metadata.load(fp)

    def test_column_types_unrecognized_column_name(self):
        fp = self.get_data_path('valid/simple.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    'not_a_column.*column_types.*not a column '
                                    'in the metadata file'):
            Metadata.load(fp, column_types={'not_a_column': 'numeric'})

    def test_column_types_unrecognized_column_type(self):
        fp = self.get_data_path('valid/simple.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    'col2.*column_types.*unrecognized column '
                                    'type.*CATEGORICAL'):
            Metadata.load(fp, column_types={'col1': 'numeric',
                                            'col2': 'CATEGORICAL'})

    def test_column_types_not_convertible_to_numeric(self):
        fp = self.get_data_path('valid/simple.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    "column 'col3' to numeric.*could not be "
                                    "interpreted as numeric: 'bar', 'foo'"):
            Metadata.load(fp, column_types={'col1': 'numeric',
                                            'col2': 'categorical',
                                            'col3': 'numeric'})

    def test_column_types_override_directive_not_convertible_to_numeric(self):
        fp = self.get_data_path('valid/simple-with-directive.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    "column 'col3' to numeric.*could not be "
                                    "interpreted as numeric: 'bar', 'foo'"):
            Metadata.load(fp, column_types={'col3': 'numeric'})

    def test_directive_before_header(self):
        fp = self.get_data_path('invalid/directive-before-header.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    'directive.*#q2:types.*searching for '
                                    'header'):
            Metadata.load(fp)

    def test_unrecognized_directive(self):
        fp = self.get_data_path('invalid/unrecognized-directive.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    'Unrecognized directive.*#q2:foo.*'
                                    '#q2:types directive is supported'):
            Metadata.load(fp)

    def test_duplicate_directives(self):
        fp = self.get_data_path('invalid/duplicate-directives.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    'duplicate directive.*#q2:types'):
            Metadata.load(fp)

    def test_unrecognized_column_type_in_directive(self):
        fp = self.get_data_path('invalid/unrecognized-column-type.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    'col2.*unrecognized column type.*foo.*'
                                    '#q2:types directive'):
            Metadata.load(fp)

    def test_column_types_directive_not_convertible_to_numeric(self):
        fp = self.get_data_path('invalid/types-directive-non-numeric.tsv')

        # This error message regex is intentionally verbose because we want to
        # assert that many different types of non-numeric strings aren't
        # interpreted as numbers. The error message displays a sorted list of
        # all values that couldn't be converted to numbers, making it possible
        # to test a variety of non-numeric strings in a single test case.
        msg = (r"column 'col2' to numeric.*could not be interpreted as "
               "numeric: '\$42', '\+inf', '-inf', '0xAF', '1,000', "
               "'1\.000\.0', '1_000_000', '1e3e4', 'Infinity', 'NA', 'NaN', "
               "'a', 'e3', 'foo', 'inf', 'nan', 'sample-1'")
        with self.assertRaisesRegex(MetadataFileError, msg):
            Metadata.load(fp)

    def test_directive_after_directives_section(self):
        fp = self.get_data_path(
            'invalid/directive-after-directives-section.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    '#q2:types.*outside of the directives '
                                    'section'):
            Metadata.load(fp)

    def test_directive_longer_than_header(self):
        fp = self.get_data_path('invalid/directive-longer-than-header.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    'row has 5 cells.*header declares 4 '
                                    'cells'):
            Metadata.load(fp)

    def test_data_longer_than_header(self):
        fp = self.get_data_path('invalid/data-longer-than-header.tsv')

        with self.assertRaisesRegex(MetadataFileError,
                                    'row has 5 cells.*header declares 4 '
                                    'cells'):
            Metadata.load(fp)

    def test_comments_and_blank_lines(self):
        fp = self.get_data_path('valid/comments-n-blanks.tsv')

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
        fp = self.get_data_path('valid/empty-rows.tsv')

        obs_md = Metadata.load(fp)

        exp_index = pd.Index(['id1', 'id2', 'id3'], name='id', dtype=object)
        exp_df = pd.DataFrame({'col1': [1.0, 2.0, 3.0],
                               'col2': ['a', 'b', 'c'],
                               'col3': ['foo', 'bar', '42']},
                              index=exp_index)
        exp_md = Metadata(exp_df)

        self.assertEqual(obs_md, exp_md)

    def test_qiime1_mapping_file(self):
        fp = self.get_data_path('valid/qiime1.tsv')

        obs_md = Metadata.load(fp)

        exp_index = pd.Index(['id1', 'id2', 'id3'], name='#SampleID',
                             dtype=object)
        exp_df = pd.DataFrame({'col1': [1.0, 2.0, 3.0],
                               'col2': ['a', 'b', 'c'],
                               'col3': ['foo', 'bar', '42']},
                              index=exp_index)
        exp_md = Metadata(exp_df)

        self.assertEqual(obs_md, exp_md)

    def test_no_columns(self):
        fp = self.get_data_path('valid/no-columns.tsv')

        obs_md = Metadata.load(fp)

        exp_index = pd.Index(['a', 'b', 'my-id'], name='id', dtype=object)
        exp_df = pd.DataFrame({}, index=exp_index, dtype=object)
        exp_md = Metadata(exp_df)

        self.assertEqual(obs_md, exp_md)

    def test_does_not_cast_ids(self):
        fp = self.get_data_path('valid/no-type-cast.tsv')

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
        fp = self.get_data_path('valid/simple.tsv')

        metadata = Metadata.load(fp)

        self.assertEqual(metadata.artifacts, ())

    def test_jagged_trailing_columns(self):
        # Test case based on https://github.com/qiime2/qiime2/issues/335
        fp = self.get_data_path('valid/jagged-trailing-columns.tsv')

        obs_md = Metadata.load(fp)

        exp_index = pd.Index(['id1', 'id2', 'id3'], name='id', dtype=object)
        exp_df = pd.DataFrame({'col1': [1.0, 2.0, 3.0],
                               'col2': ['a', 'b', 'c'],
                               'col3': ['foo', 'bar', '42']},
                              index=exp_index)
        exp_md = Metadata(exp_df)

        self.assertEqual(obs_md, exp_md)


if __name__ == '__main__':
    unittest.main()
