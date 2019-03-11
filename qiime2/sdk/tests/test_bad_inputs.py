# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import collections
import tempfile
import qiime2.core.archive as archive
# import qiime2.plugin
# from qiime2.core.type import MethodSignature, Int

from qiime2.core.testing.util import get_dummy_plugin
from qiime2.plugin.testing import TestPluginBase

from qiime2.sdk import Artifact, Visualization
from qiime2.core.testing.type import IntSequence1
from qiime2.core.testing.visualizer import most_common_viz

# from qiime2.sdk import Method, Results
# from qiime2.core.testing.method import (concatenate_ints, merge_mappings,
#                                         split_ints, params_only_method,
#                                         no_input_method)
# from qiime2.core.testing.type import (IntSequence2, SingleInt, Mapping)


class TestBadInputs(TestPluginBase):
    def make_provenance_capture(self):
        # You can't actually import a visualization, but I won't tell
        # visualization if you don't...
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
        # generate a sample viz and other params
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

    def test_artifact_passed_as_param(self):
        # TODO: implement
        # passed as params
        # passed as metadata
        # etc?
        pass

    def test_incorrect_artifact_type(self):
        # method = self.plugin.methods['optional_artifacts_method']
        # ints1 = Artifact.import_data(IntSequence1, [0, 42, 43])
        # testDoc = Artifact.import_data
        #
        #
        # with self.assertRaisesRegex(
        #         TypeError, 'Visualizations may not be used as inputs.'):
        #     method(saved_viz, 42)
        pass

    def test_incorrect_artifact_subtype(self):
        # TODO: implement
        pass

    def incorrect_primitive_type(self):
        # TODO: implement
        pass

    def test_metadata_passed_as_artifact(self):
        # TODO: implement
        pass

    def test_primitive_passed_as_input(self):
        # generate params
        concatenate_ints = self.plugin.methods['concatenate_ints']
        identity_with_metadata = self.plugin.methods['identity_with_metadata']
        params_only_method = self.plugin.methods['params_only_method']
        ints1 = Artifact.import_data(IntSequence1, [0, 42, 43])
        ints3 = Artifact.import_data(IntSequence1, [99, -22])
        int1 = 4
        int2 = 5
        arbitrary_int = 43

        # tests primitive passed as IntSequence artifact
        with self.assertRaisesRegex(TypeError,
                                    'type IntSequence1.*type \'int\''):
            concatenate_ints(ints1, arbitrary_int, ints3, int1, int2)

        # tests primitive passed as metadata
        with self.assertRaisesRegex(TypeError,
                                    'type Metadata.*type \'int\''):
            identity_with_metadata(ints1, arbitrary_int)

        # tests wrong type of primitive passed
        with self.assertRaisesRegex(TypeError,
                                    'type Int.*type \'str\''):
            params_only_method('key string', 'arbitrary string shuld be int')

    def test_primitive_param_out_of_range(self):
        # TODO: implement
        pass

    def test_primitive_param_not_valid_choice(self):
        # TODO: implement
        pass
