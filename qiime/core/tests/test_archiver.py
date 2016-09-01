# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
import os
import tempfile
import unittest
import zipfile

from qiime.core.archiver import Archiver


class TestArchiver(unittest.TestCase):
    def setUp(self):
        prefix = "qiime2-test-temp-"
        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)
        self.root_dir = self.temp_dir.name

        # FakeArchiver is for testing some of the private methods on Archiver
        # that don't require the full class.
        class FakeArchiver(Archiver):
            # Explicitly not calling super()
            def __init__(self, *args, **kwargs):
                self._temp_dir = tempfile.mkdtemp(prefix)
                self._data_dir = os.path.join(self._temp_dir,
                                              self.DATA_DIRNAME)
                self._pid = os.getpid()
                self._uuid = 'foo'
                self._type = 'bar'
                self._provenance = 'baz'
        self.FakeArchiver = FakeArchiver

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_load_version_file_decoding(self):
        fp = os.path.join(self.root_dir, "bad_data.qza")

        # Bypass Archiver.save to build a faux-qza file with a
        # non-UTF-8 encoded version file.
        with zipfile.ZipFile(fp, mode="w",
                             compression=zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zf:
            zf.writestr(os.path.join(self.root_dir,
                                     Archiver._VERSION_FILENAME),
                        "version info".encode("utf-16"))

            with self.assertRaises(UnicodeDecodeError):
                Archiver._load_version(zf, self.root_dir)

    def test_load_metadata_file_decoding(self):
        fp = os.path.join(self.root_dir, "bad_data.qza")

        # Bypass Archiver.save to build a faux-qza file with a
        # non-UTF-8 encoded metadata file.
        with zipfile.ZipFile(fp, mode="w",
                             compression=zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zf:
            zf.writestr(os.path.join(self.root_dir,
                                     Archiver._METADATA_FILENAME),
                        "metadata info".encode("utf-16"))

            with self.assertRaises(UnicodeDecodeError):
                Archiver._load_metadata(zf, self.root_dir)

    def test_save_version_file_encoding(self):
        archive = self.FakeArchiver()
        fp = os.path.join(self.root_dir, "archive.qza")

        with zipfile.ZipFile(fp, mode="w",
                             compression=zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zf:
            archive._save_version(zf, self.root_dir)

            # Manually load the saved file with strict errors on (this will
            # raise a UnicodeDecodeError if there is a mismatch
            # during decoding).
            version_path = os.path.join(self.root_dir,
                                        Archiver._VERSION_FILENAME)
            with zf.open(version_path) as bytes_fh:
                with io.TextIOWrapper(bytes_fh, newline=None,
                                      encoding='utf-8', errors='strict') as fh:
                    self.assertEqual(fh.read().rstrip('\n'), Archiver._VERSION)

    def test_save_readme_file_encoding(self):
        archive = self.FakeArchiver()
        fp = os.path.join(self.root_dir, "archive.qza")

        with zipfile.ZipFile(fp, mode="w",
                             compression=zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zf:
            archive._save_readme(zf, self.root_dir)

            # Manually load the saved file with strict errors on (this will
            # raise a UnicodeDecodeError if there is a mismatch during
            # decoding).
            readme_path = os.path.join(self.root_dir,
                                       Archiver._README_FILENAME)
            with zf.open(readme_path) as bytes_fh:
                with io.TextIOWrapper(bytes_fh, newline=None,
                                      encoding='utf-8', errors='strict') as fh:
                    self.assertTrue(fh.read().rstrip('\n'))

    def test_save_metadata_file_encoding(self):
        archive = self.FakeArchiver()
        fp = os.path.join(self.root_dir, "archive.qza")

        with zipfile.ZipFile(fp, mode="w",
                             compression=zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zf:
            archive._save_metadata(zf, self.root_dir)

            # Manually load the saved file with strict errors on (this will
            # raise a UnicodeDecodeError if there is a mismatch during
            # decoding).
            metadata_path = os.path.join(self.root_dir,
                                         Archiver._METADATA_FILENAME)
            with zf.open(metadata_path) as bytes_fh:
                with io.TextIOWrapper(bytes_fh, newline=None,
                                      encoding='utf-8', errors='strict') as fh:
                    self.assertTrue(fh.read().rstrip('\n'))

    def test_metadata_pprint_yaml(self):
        archive = self.FakeArchiver()
        fp = os.path.join(self.root_dir, "archive.qza")

        with zipfile.ZipFile(fp, mode="w",
                             compression=zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zf:
            archive._save_metadata(zf, self.root_dir)

            with zf.open(os.path.join(self.root_dir,
                                      Archiver._METADATA_FILENAME)) as zip_fh:
                with io.TextIOWrapper(zip_fh, newline=None,
                                      encoding='utf-8', errors='strict') as fh:
                    self.assertEqual(fh.read(),
                                     "uuid: foo\n"
                                     "type: '''bar'''\n"
                                     "provenance: baz\n")


if __name__ == '__main__':
    unittest.main()
