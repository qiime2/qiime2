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
from pathlib import Path

from qiime2.core.testing.plugin import IntSequenceDirectoryFormat


class TestBoundFile(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp')
        self.path = Path(os.path.join(self.test_dir.name, 'file'))
        with open(self.path, 'w') as fh:
            fh.write('1\n10')
        self.format = IntSequenceDirectoryFormat(self.path, mode='r')

    def tearDown(self):
        self.test_dir.cleanup()

    def test_filepath_expected(self):
        self.assertEqual(self.path, self.format.file.path)
