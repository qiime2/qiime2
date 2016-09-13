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
import uuid
import zipfile

from qiime.core.archiver import Archiver
from qiime.core.testing.type import IntSequence1


class TestArchiver(unittest.TestCase):
    def setUp(self):
        prefix = "qiime2-test-temp-"
        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)

        # Initialize an Archiver. The values passed to the constructor mostly
        # don't matter to the Archiver, but we'll pass valid Artifact test data
        # anyways in case Archiver's behavior changes in the future.
        def data_initializer(data_dir):
            fp = os.path.join(data_dir, 'ints.txt')
            with open(fp, 'w') as fh:
                fh.write('1\n')
                fh.write('2\n')
                fh.write('3\n')

        self.archiver = Archiver(uuid.uuid4(), IntSequence1,
                                 'IntSequenceDirectoryFormat', None,
                                 data_initializer=data_initializer)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_load_version_file_decoding(self):
        fp = os.path.join(self.temp_dir.name, "bad_data.qza")

        # Bypass Archiver.save to build a faux-qza file with a
        # non-UTF-8 encoded version file.
        with zipfile.ZipFile(fp, mode="w",
                             compression=zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zf:
            zf.writestr(os.path.join('bad_data', Archiver._VERSION_FILENAME),
                        "version info".encode("utf-16"))

            with self.assertRaises(UnicodeDecodeError):
                Archiver._load_version(zf, 'bad_data')

    def test_load_metadata_file_decoding(self):
        fp = os.path.join(self.temp_dir.name, "bad_data.qza")

        # Bypass Archiver.save to build a faux-qza file with a
        # non-UTF-8 encoded metadata file.
        with zipfile.ZipFile(fp, mode="w",
                             compression=zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zf:
            zf.writestr(os.path.join('bad_data', Archiver._METADATA_FILENAME),
                        "metadata info".encode("utf-16"))

            with self.assertRaises(UnicodeDecodeError):
                Archiver._load_metadata(zf, 'bad_data')

    def test_save_version_file_encoding(self):
        fp = os.path.join(self.temp_dir.name, "archive.qza")

        with zipfile.ZipFile(fp, mode="w",
                             compression=zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zf:
            self.archiver._save_version(zf, 'archive')

            # Manually load the saved file with strict errors on (this will
            # raise a UnicodeDecodeError if there is a mismatch
            # during decoding).
            version_path = os.path.join('archive', Archiver._VERSION_FILENAME)
            with zf.open(version_path) as bytes_fh:
                with io.TextIOWrapper(bytes_fh, newline=None,
                                      encoding='utf-8', errors='strict') as fh:
                    self.assertEqual(fh.read().rstrip('\n'), Archiver._VERSION)

    def test_save_readme_file_encoding(self):
        fp = os.path.join(self.temp_dir.name, "archive.qza")

        with zipfile.ZipFile(fp, mode="w",
                             compression=zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zf:
            self.archiver._save_readme(zf, 'archive')

            # Manually load the saved file with strict errors on (this will
            # raise a UnicodeDecodeError if there is a mismatch during
            # decoding).
            readme_path = os.path.join('archive', Archiver._README_FILENAME)
            with zf.open(readme_path) as bytes_fh:
                with io.TextIOWrapper(bytes_fh, newline=None,
                                      encoding='utf-8', errors='strict') as fh:
                    self.assertTrue(fh.read().rstrip('\n'))

    def test_save_metadata_file_encoding(self):
        fp = os.path.join(self.temp_dir.name, "archive.qza")

        with zipfile.ZipFile(fp, mode="w",
                             compression=zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zf:
            self.archiver._save_metadata(zf, 'archive')

            # Manually load the saved file with strict errors on (this will
            # raise a UnicodeDecodeError if there is a mismatch during
            # decoding).
            metadata_path = os.path.join('archive',
                                         Archiver._METADATA_FILENAME)
            with zf.open(metadata_path) as bytes_fh:
                with io.TextIOWrapper(bytes_fh, newline=None,
                                      encoding='utf-8', errors='strict') as fh:
                    self.assertTrue(fh.read().rstrip('\n'))

    def test_metadata_pprint_yaml(self):
        fp = os.path.join(self.temp_dir.name, "archive.qza")

        with zipfile.ZipFile(fp, mode="w",
                             compression=zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zf:
            self.archiver._save_metadata(zf, 'archive')

            with zf.open(os.path.join('archive',
                                      Archiver._METADATA_FILENAME)) as zip_fh:
                with io.TextIOWrapper(zip_fh, newline=None,
                                      encoding='utf-8', errors='strict') as fh:
                    self.assertEqual(fh.read(),
                                     "uuid: %s\n" % self.archiver.uuid +
                                     "type: IntSequence1\n"
                                     "format: IntSequenceDirectoryFormat\n"
                                     "provenance: None\n")

    def test_save_invalid_filepath(self):
        # Empty filepath.
        with self.assertRaisesRegex(ValueError, 'empty'):
            self.archiver.save('')

        # Directory.
        with self.assertRaisesRegex(ValueError, 'directory'):
            self.archiver.save(self.temp_dir.name)

        # Ends with path separator (no basename, e.g. /tmp/foo/).
        with self.assertRaisesRegex(ValueError, 'path separator'):
            self.archiver.save(os.path.join(self.temp_dir.name, 'foo', ''))

    def test_save_ignores_dotfiles(self):
        def data_initializer(data_dir):
            fp = os.path.join(data_dir, 'ints.txt')
            with open(fp, 'w') as fh:
                fh.write('1\n')
                fh.write('2\n')
                fh.write('3\n')

            hidden_fp = os.path.join(data_dir, '.hidden-file')
            with open(hidden_fp, 'w') as fh:
                fh.write("You can't see me if I can't see you\n")

            hidden_dir = os.path.join(data_dir, '.hidden-dir')
            os.mkdir(hidden_dir)
            with open(os.path.join(hidden_dir, 'ignored-file'), 'w') as fh:
                fh.write("I'm ignored because I live in a hidden dir :(\n")

        archiver = Archiver(uuid.uuid4(), IntSequence1,
                            'IntSequenceDirectoryFormat', None,
                            data_initializer=data_initializer)

        fp = os.path.join(self.temp_dir.name, 'archive.qza')
        archiver.save(fp)

        with zipfile.ZipFile(fp, mode='r') as zf:
            fps = set(zf.namelist())
            expected = {
                'archive/VERSION',
                'archive/metadata.yaml',
                'archive/README.md',
                'archive/data/ints.txt'
            }
            self.assertEqual(fps, expected)


if __name__ == '__main__':
    unittest.main()
