# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
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
        self.execute_examples()


if __name__ == '__main__':
    unittest.main()
