# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import unittest
import tempfile

from qiime2.core.testing.util import get_dummy_plugin
from qiime2.plugin.testing import TestPluginBase


class TestTesting(TestPluginBase):
    def setUp(self):
        self.plugin = get_dummy_plugin()

        # TODO standardize temporary directories created by QIIME 2
        # create a temporary data_dir for sample Visualizations
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        self.data_dir = os.path.join(self.test_dir.name, 'viz-output')
        os.mkdir(self.data_dir)

    def tearDown(self):
        self.test_dir.cleanup()

    def test_examples(self):
        def concatenate_ints_simple(use, scope):
            self.assertEqual(4, len(scope))

        def concatenate_ints_complex(use, scope):
            with self.subTest(test='scope length'):
                self.assertEqual(5, len(scope))

            with self.subTest(test='output type'):
                final_record = scope[list(scope.keys())[-1]]
                artifact = final_record.factory()
                self.assertEqual('IntSequence1', str(artifact.type))

        def identity_with_metadata_case_a(use, scope):
            self.assertEqual(3, len(scope))

        def most_common_viz_typical(use, scope):
            self.assertEqual(2, len(scope))

            final_record = scope[list(scope.keys())[-1]]
            self.assertEqual('visualization', final_record.type)

        callbacks = {
            'concatenate_ints_simple': concatenate_ints_simple,
            'concatenate_ints_complex': concatenate_ints_complex,
            'identity_with_metadata_case_a': identity_with_metadata_case_a,
            'most_common_viz_typical': most_common_viz_typical,
        }

        self.execute_examples(callbacks=callbacks)


if __name__ == '__main__':
    unittest.main()
