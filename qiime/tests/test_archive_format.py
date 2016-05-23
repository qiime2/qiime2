# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import tempfile
import os

import qiime.plugin


class DummyDataReader:
    # a temporary class for testing purposes

    def __init__(self, data_dir):
        self._data_dir = data_dir

    def get_paths(self):
        return [path for path in os.listdir(self._data_dir)
                if os.path.isfile(path)]

    def get_file(self, path):
        return open(os.path.join(self._data_dir, path))


class DummyDataWriter:
    # a temporary class for testing purposes

    def __init__(self, data_dir):
        self._data_dir = data_dir

    def create_file(self, path):
        return open(os.path.join(self._data_dir, path), 'w')


class TestFormat(unittest.TestCase):

    def test_init(self):
        def is_valid(data_reader):
            return True
        f = qiime.plugin.ArchiveFormat('my-format', '1.0.0', is_valid)
        self.assertEqual(f.name, 'my-format')
        self.assertEqual(f.version, '1.0.0')
        with tempfile.TemporaryDirectory() as data_dir:
            ddr = DummyDataReader(data_dir)
            # TODO: add function directly testing validate
            self.assertTrue(f.validate(ddr))
        self.assertEqual(f.get_reader_views(), set())
        self.assertEqual(f.get_writer_views(), set())

    def test_equality_operators(self):
        def is_valid(data_reader):
            return True
        f1 = qiime.plugin.ArchiveFormat('my-format', '1.0.0', is_valid)
        f2 = qiime.plugin.ArchiveFormat('my-format', '1.0.0', is_valid)
        f3 = qiime.plugin.ArchiveFormat('your-format', '1.0.0', is_valid)
        f4 = qiime.plugin.ArchiveFormat('my-format', '1.0.1', is_valid)
        self.assertEqual(f1, f1)
        self.assertEqual(f1, f2)
        self.assertNotEqual(f1, f3)
        self.assertNotEqual(f1, f4)
        self.assertNotEqual(f3, f4)

    def test_reader(self):
        def is_valid(data_reader):
            return True
        f = qiime.plugin.ArchiveFormat('my-format', '1.0.0', is_valid)

        @f.reader(float)
        def my_reader(data_reader):
            return 42.0
        # TODO: add functions specifically for testing get_reader_views
        self.assertEqual(f.get_reader_views(), set([float]))
        with tempfile.TemporaryDirectory() as data_dir:
            ddr = DummyDataReader(data_dir)
            self.assertEqual(f.readers[float](ddr), 42.0)

    def test_writer(self):
        def is_valid(data_reader):
            return True
        f = qiime.plugin.ArchiveFormat('my-format', '1.0.0', is_valid)

        @f.writer(float)
        def my_writer(data_writer):
            return None
        # TODO: add functions specifically for testing get_writer_views
        self.assertEqual(f.get_writer_views(), set([float]))
        with tempfile.TemporaryDirectory() as data_dir:
            ddw = DummyDataWriter(data_dir)
            self.assertEqual(f.writers[float](ddw), None)


if __name__ == "__main__":
    unittest.main()
