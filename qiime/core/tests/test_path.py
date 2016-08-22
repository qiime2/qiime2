# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import unittest

from qiime.core.path import TempPath


class TestTempPath(unittest.TestCase):
    def test_new_temppath(self):
        f = TempPath()
        self.assertIsInstance(f, TempPath)
        self.assertTrue(f.is_file())

        g = TempPath(dir=True)
        self.assertIsInstance(g, TempPath)
        self.assertTrue(g.is_dir())

    def test_new_temppath_context_mgr(self):
        with TempPath() as f:
            path = str(f)
            self.assertIsInstance(f, TempPath)
            self.assertTrue(os.path.isfile(path))
        self.assertFalse(os.path.isfile(path))

    def test_finalize(self):
        f = TempPath()
        path = str(f)

        self.assertTrue(os.path.isfile(path))
        f._finalize()
        self.assertFalse(os.path.isfile(path))
