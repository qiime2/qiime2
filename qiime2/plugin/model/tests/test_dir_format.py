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
from qiime2.core.testing.plugin import IntSequenceDirectoryFormat


class TestBoundFile(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp')
        self.dummy_plugin = get_dummy_plugin()

    def test_filepath_expected(self):
        path = os.path.join(self.test_dir.name, 'file')
        with open(path, 'w') as fh:
            fh.write('1\n10')
        format = IntSequenceDirectoryFormat(path, mode='r')
        self.assertEqual(path, format.file.path)
