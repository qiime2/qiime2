# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import errno
import tempfile
import unittest
import unittest.mock as mock

import qiime2.util as util

EXDEV = OSError(errno.EXDEV, "Invalid cross-device link")
EPERM = PermissionError(errno.EPERM, "unsafe operation e.g. using link wrong")
EACCES = PermissionError(errno.EACCES, "insufficient r/w permissions")
SECRET = "this is a secret for testing, don't tell anyone!"


class TestDuplicate(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        self.dst1 = os.path.join(self.test_dir.name, 'dst1')
        self.dst2 = os.path.join(self.test_dir.name, 'dst2')
        self.dir = os.path.join(self.test_dir.name, 'dir')
        with open(self.dst2, 'w') as fh:
            fh.write("This is not the secret")
        os.mkdir(self.dir)

        self.src = os.path.join(self.test_dir.name, 'src')
        self.missing = os.path.join(self.test_dir.name, 'missing')
        with open(self.src, 'w') as fh:
            fh.write(SECRET)

    def tearDown(self):
        self.test_dir.cleanup()

    def test_src_not_exists(self):
        with self.assertRaisesRegex(FileNotFoundError, self.missing):
            util.duplicate(self.missing, self.dst1)

    def test_src_dir(self):
        with self.assertRaisesRegex(IsADirectoryError, self.dir):
            util.duplicate(self.dir, self.dst1)

    def test_dst_not_exists(self):
        util.duplicate(self.src, self.dst1)

        assert os.path.exists(self.dst1)
        with open(self.dst1) as fh:
            self.assertEqual(fh.read(), SECRET)

    def test_dst_exists(self):
        with self.assertRaisesRegex(FileExistsError, self.dst2):
            util.duplicate(self.src, self.dst2)

    def test_dst_dir(self):
        with self.assertRaisesRegex(IsADirectoryError, self.dir):
            util.duplicate(self.src, self.dir)

    @mock.patch('qiime2.util.os.link', side_effect=EACCES)
    def test_perm_error_EACCES(self, mocked_link):
        with self.assertRaisesRegex(
                PermissionError, "insufficient r/w permissions"):
            util.duplicate(self.src, self.dst1)

        assert mocked_link.called

    @mock.patch('qiime2.util.os.link', side_effect=EPERM)
    def test_perm_error_EPERM(self, mocked_link):
        util.duplicate(self.src, self.dst1)

        assert mocked_link.called
        assert os.path.exists(self.dst1)
        with open(self.dst1) as fh:
            self.assertEqual(fh.read(), SECRET)

    @mock.patch('qiime2.util.os.link', side_effect=EXDEV)
    def test_cross_device_src_not_exists(self, mocked_link):
        with self.assertRaisesRegex(FileNotFoundError, self.missing):
            util.duplicate(self.missing, self.dst1)

    @mock.patch('qiime2.util.os.link', side_effect=EXDEV)
    def test_cross_device_src_dir(self, mocked_link):
        with self.assertRaisesRegex(IsADirectoryError, self.dir):
            util.duplicate(self.dir, self.dst1)

    @mock.patch('qiime2.util.os.link', side_effect=EXDEV)
    def test_cross_device_dst_not_exists(self, mocked_link):
        util.duplicate(self.src, self.dst1)

        assert mocked_link.called
        assert os.path.exists(self.dst1)
        with open(self.dst1) as fh:
            self.assertEqual(fh.read(), SECRET)

    @mock.patch('qiime2.util.os.link', side_effect=EXDEV)
    def test_cross_device_dst_exists(self, mocked_link):
        with self.assertRaisesRegex(FileExistsError, self.dst2):
            util.duplicate(self.src, self.dst2)

    @mock.patch('qiime2.util.os.link', side_effect=EXDEV)
    def test_cross_device_dst_dir(self, mocked_link):
        with self.assertRaisesRegex(IsADirectoryError, self.dir):
            util.duplicate(self.src, self.dir)

    @mock.patch('qiime2.util.os.link', side_effect=EXDEV)
    @mock.patch('qiime2.util.shutil.copyfile', side_effect=EACCES)
    def test_cross_device_perm_error(self, mocked_link, mocked_copyfile):
        with self.assertRaisesRegex(
                PermissionError, "insufficient r/w permissions"):
            util.duplicate(self.src, self.dst1)

        assert mocked_link.called
        assert mocked_copyfile.called


if __name__ == '__main__':
    unittest.main()
