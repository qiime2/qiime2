# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import unittest

from qiime2.core.path import OutPath


class TestOutPath(unittest.TestCase):
    def test_new_outpath(self):
        f = OutPath()
        self.assertIsInstance(f, OutPath)
        self.assertTrue(f.is_file())

        g = OutPath(dir=True)
        self.assertIsInstance(g, OutPath)
        self.assertTrue(g.is_dir())

    def test_new_outpath_context_mgr(self):
        with OutPath() as f:
            path = str(f)
            self.assertIsInstance(f, OutPath)
            self.assertTrue(os.path.isfile(path))
        self.assertFalse(os.path.isfile(path))

    def test_destructor(self):
        f = OutPath()
        path = str(f)

        self.assertTrue(os.path.isfile(path))
        f._destructor()
        self.assertFalse(os.path.isfile(path))


if __name__ == '__main__':
    unittest.main()
