# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import collections
import tempfile
import unittest
import warnings

import qiime2.core.archive as archive

from qiime2.core.testing.util import get_dummy_plugin
from qiime2.plugin.testing import TestPluginBase

from qiime2.sdk import Artifact, Visualization
from qiime2.core.testing.type import IntSequence1, IntSequence2, SingleInt
from qiime2.core.testing.visualizer import most_common_viz
from qiime2 import Metadata
from qiime2.metadata.tests.test_io import get_data_path


# NOTE: This test suite exists for tests not easily split into
# test_method, test_visualizer, test_pipeline
# TestBadInputs tests type mismatches between Action signatures and passed args


class TestBadInputs(TestPluginBase):

    def make_provenance_capture(self):
        # importing visualizations is not supported, but we do that here to
        # simplify testing machinery
        return archive.ImportProvenanceCapture()

    def setUp(self):
        self.plugin = get_dummy_plugin()

        # TODO standardize temporary directories created by QIIME 2
        # create a temporary data_dir for sample Visualizations
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        self.data_dir = os.path.join(self.test_dir.name, 'viz-output')
        os.mkdir(self.data_dir)
        most_common_viz(self.data_dir, collections.Counter(range(42)))

    def tearDown(self):
        self.test_dir.cleanup()

    def test_viz_passed_as_input(self):
        saved_viz = Visualization._from_data_dir(
            self.data_dir, self.make_provenance_capture())
        method = self.plugin.methods['optional_artifacts_method']
        ints1 = Artifact.import_data(IntSequence1, [0, 42, 43])

        # tests Viz passed as primitive parameter
        with self.assertRaisesRegex(
                TypeError, 'Visualizations may not be used as inputs.'):
            method(saved_viz, 42)

        # tests Viz passed as Artifact input
        with self.assertRaisesRegex(
                TypeError, 'Visualizations may not be used as inputs.'):
            method(ints1, 42, optional1=saved_viz)

        # tests Viz passed as metadata
        method = self.plugin.methods['identity_with_optional_metadata']
        with self.assertRaisesRegex(
                TypeError, 'Visualizations may not be used as inputs.'):
            method(ints1, metadata=saved_viz)

    def test_artifact_passed_incorrectly(self):
        concatenate_ints = self.plugin.methods['concatenate_ints']
        identity_with_metadata = self.plugin.methods['identity_with_metadata']
        ints1 = Artifact.import_data(IntSequence1, [0, 42, 43])
        ints2 = Artifact.import_data(IntSequence1, [99, -22])
        ints3 = Artifact.import_data(IntSequence2, [12, 111])
        inappropriate_Artifact = Artifact.import_data(IntSequence1, [-9999999])
        int1 = 4
        int2 = 5

        # tests Artifact passed as integer
        with self.assertRaisesRegex(
                TypeError, 'int1.*type Int.*IntSequence1'):
            concatenate_ints(ints1, ints2, ints3, inappropriate_Artifact, int2)

        # tests Artifact passed as metadata
        with self.assertRaisesRegex(
                TypeError, '\'metadata\'.*type Metadata.*IntSequence1'):
            identity_with_metadata(ints1, inappropriate_Artifact)

        # tests wrong type of Artifact passed
        with self.assertRaisesRegex(
                TypeError, 'ints3.*IntSequence2.*IntSequence1'):
            concatenate_ints(ints1, ints2, inappropriate_Artifact, int1, int2)

    def test_primitive_passed_incorrectly(self):
        concatenate_ints = self.plugin.methods['concatenate_ints']
        identity_with_metadata = self.plugin.methods['identity_with_metadata']
        params_only_method = self.plugin.methods['params_only_method']

        md_fp = get_data_path('valid/simple.tsv')
        inappropriate_metadata = Metadata.load(md_fp)

        ints1 = Artifact.import_data(IntSequence1, [0, 42, 43])
        ints3 = Artifact.import_data(IntSequence1, [12, 111])
        int1 = 4
        int2 = 5
        arbitrary_int = 43

        # tests primitive int passed as IntSequence artifact
        with self.assertRaisesRegex(TypeError,
                                    'ints2.*43.*incompatible.*IntSequence1'):
            concatenate_ints(ints1, arbitrary_int, ints3, int1, int2)

        # tests primitive passed as metadata
        with self.assertRaisesRegex(TypeError,
                                    'metadata.*43.*incompatible.*Metadata'):
            identity_with_metadata(ints1, arbitrary_int)

        # tests wrong type of primitive passed
        with self.assertRaisesRegex(TypeError,
                                    'age.*arbitraryString.*incompatible.*Int'):
            params_only_method('key string', 'arbitraryString')

        # tests metadata passed as artifact
        with self.assertRaisesRegex(TypeError,
                                    '\'ints2\'.*Metadata.*IntSequence1'):
            concatenate_ints(ints1, inappropriate_metadata, ints3, int1, int2)

    def test_primitive_param_out_of_range(self):
        range_nested_in_list = self.plugin.methods['variadic_input_method']
        range_not_nested_in_list = self.plugin.visualizers['params_only_viz']
        ints_list = [Artifact.import_data(IntSequence1, [0, 42, 43]),
                     Artifact.import_data(IntSequence2, [4, 5, 6])]
        int_set = {Artifact.import_data(SingleInt, 7),
                   Artifact.import_data(SingleInt, 8)}
        nums = {9, 10}
        bad_range_val = [11, 12, -9999]
        invalid_age = -99999

        # Tests primitives of correct type but outside of Range...
        # ... in a list
        with self.assertRaisesRegex(
                TypeError, 'opt_nums.*-9999.*incompatible.*List'):
            range_nested_in_list(ints_list, int_set, nums, bad_range_val)

        # ... not in a list
        with self.assertRaisesRegex(
                TypeError,
                r'\'age\'.*-99999.*incompatible.*Int % Range\(0, None\)'):
            range_not_nested_in_list('John Doe', invalid_age)

    def test_primitive_param_not_valid_choice(self):
        pipeline = self.plugin.pipelines['failing_pipeline']
        int_sequence = Artifact.import_data(IntSequence1, [0, 42, 43])
        break_from = "invalid choice"

        # test String not a valid choice
        with self.assertRaisesRegex(
                TypeError, 'break_from.*\'invalid choice\''):
            pipeline(int_sequence, break_from)


class TestDeprecation(unittest.TestCase):
    def setUp(self):
        self.plugin = get_dummy_plugin()
        self.method = self.plugin.methods['deprecated_method']

    def test_successful_registration(self):
        self.assertTrue(self.method.deprecated)

    def test_deprecation_warning(self):
        with warnings.catch_warnings(record=True) as w:
            self.method()
            self.assertEqual(1, len(w))
            warning = w[0]
            self.assertEqual(warning.category, FutureWarning)
            self.assertTrue('Method is deprecated' in str(warning.message))

    def test_docstring(self):
        self.assertIn('Method is deprecated', self.method.__call__.__doc__)
