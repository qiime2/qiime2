# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import tempfile

import qiime.plugin


class TestFormat(unittest.TestCase):
    def test_init(self):
        def is_valid(data_dir):
            return True

        f = qiime.plugin.DataLayout('my-format', '1.0.0', is_valid)
        self.assertEqual(f.name, 'my-format')
        self.assertEqual(f.version, '1.0.0')
        with tempfile.TemporaryDirectory() as data_dir:
            # TODO: add function directly testing validate
            self.assertTrue(f.validate(data_dir))
        self.assertEqual(f.get_reader_views(), set())
        self.assertEqual(f.get_writer_views(), set())

    def test_equality_operators(self):
        def is_valid(data_dir):
            return True

        f1 = qiime.plugin.DataLayout('my-format', '1.0.0', is_valid)
        f2 = qiime.plugin.DataLayout('my-format', '1.0.0', is_valid)
        f3 = qiime.plugin.DataLayout('your-format', '1.0.0', is_valid)
        f4 = qiime.plugin.DataLayout('my-format', '1.0.1', is_valid)
        self.assertEqual(f1, f1)
        self.assertEqual(f1, f2)
        self.assertNotEqual(f1, f3)
        self.assertNotEqual(f1, f4)
        self.assertNotEqual(f3, f4)

    def test_reader(self):
        def is_valid(data_dir):
            return True

        f = qiime.plugin.DataLayout('my-format', '1.0.0', is_valid)

        @f.reader(float)
        def my_reader(data_dir):
            return 42.0

        # TODO: add functions specifically for testing get_reader_views
        self.assertEqual(f.get_reader_views(), set([float]))
        with tempfile.TemporaryDirectory() as data_dir:
            self.assertEqual(f.readers[float](data_dir), 42.0)

    def test_writer(self):
        def is_valid(data_dir):
            return True

        f = qiime.plugin.DataLayout('my-format', '1.0.0', is_valid)

        @f.writer(float)
        def my_writer(data_dir):
            return None

        # TODO: add functions specifically for testing get_writer_views
        self.assertEqual(f.get_writer_views(), set([float]))
        with tempfile.TemporaryDirectory() as data_dir:
            self.assertEqual(f.writers[float](data_dir), None)


if __name__ == "__main__":
    unittest.main()
