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

import qiime2.plugin.model as model


class TestTextFileFormat(unittest.TestCase):
    PAYLOAD = "Somewhere over the rainbow."

    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')

    def tearDown(self):
        self.test_dir.cleanup()

    def test_open_read_good(self):
        path = os.path.join(self.test_dir.name, 'file')
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write(self.PAYLOAD)

        ff = model.TextFileFormat(path, mode='r')
        with ff.open() as fh:
            self.assertEqual(self.PAYLOAD, fh.read())

    def test_open_read_ignore_bom(self):
        path = os.path.join(self.test_dir.name, 'file')
        with open(path, 'w', encoding='utf-8-sig') as fh:
            fh.write(self.PAYLOAD)

        ff = model.TextFileFormat(path, mode='r')
        with ff.open() as fh:
            self.assertEqual(self.PAYLOAD, fh.read())

    def test_open_write_good(self):
        ff = model.TextFileFormat()

        with ff.open() as fh:
            fh.write(self.PAYLOAD)

        with open(str(ff), mode='r', encoding='utf-8') as fh:
            self.assertEqual(self.PAYLOAD, fh.read())

    def test_open_write_no_bom(self):
        ff = model.TextFileFormat()

        with ff.open() as fh:
            fh.write(self.PAYLOAD)

        with open(str(ff), mode='rb') as fh:
            self.assertEqual(b'S', fh.read(1))


if __name__ == '__main__':
    unittest.main()
