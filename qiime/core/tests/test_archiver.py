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
from qiime.core.testing.util import ArchiveTestingMixin


class TestArchiver(unittest.TestCase, ArchiveTestingMixin):
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

        self.archiver = Archiver(IntSequence1, 'IntSequenceDirectoryFormat',
                                 None, data_initializer=data_initializer)

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

    def test_save_excludes_dotfiles_in_data_dir(self):
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

        archiver = Archiver(IntSequence1, 'IntSequenceDirectoryFormat', None,
                            data_initializer=data_initializer)

        fp = os.path.join(self.temp_dir.name, 'archive.qza')
        archiver.save(fp)

        root_dir = str(archiver.uuid)
        expected = {
            'VERSION',
            'metadata.yaml',
            'README.md',
            'data/ints.txt'
        }

        self.assertArchiveMembers(fp, root_dir, expected)

    def test_save_archive_members(self):
        fp = os.path.join(self.temp_dir.name, 'archive.qza')

        self.archiver.save(fp)

        root_dir = str(self.archiver.uuid)
        expected = {
            'VERSION',
            'metadata.yaml',
            'README.md',
            'data/ints.txt'
        }

        self.assertArchiveMembers(fp, root_dir, expected)

    def test_load_archive(self):
        fp = os.path.join(self.temp_dir.name, 'archive.qza')
        self.archiver.save(fp)

        archiver = Archiver.load(fp)

        self.assertEqual(archiver.uuid, self.archiver.uuid)
        self.assertEqual(archiver.type, IntSequence1)
        self.assertEqual(archiver.format, 'IntSequenceDirectoryFormat')
        self.assertIsNone(archiver.provenance)
        self.assertEqual({e for (e, _) in archiver.get_data_paths()},
                         {'ints.txt'})

    def test_load_ignores_root_dotfiles(self):
        fp = os.path.join(self.temp_dir.name, 'archive.qza')
        self.archiver.save(fp)

        # Add some dotfiles to the archive.
        with zipfile.ZipFile(fp, mode='a') as zf:
            zf.writestr('.DS_Store', "The world's most beloved file\n")
            zf.writestr('.hidden-file',
                        "You can't see me if I can't see you\n")
            zf.writestr('.hidden-dir/ignored-file',
                        "I'm ignored because I live in a hidden dir :(\n")

        # Assert the expected files exist in the archive to verify this test
        # case is testing what we want it to.
        with zipfile.ZipFile(fp, mode='r') as zf:
            root_dir = str(self.archiver.uuid)
            expected = {
                '.DS_Store',
                '.hidden-file',
                '.hidden-dir/ignored-file',
                '%s/VERSION' % root_dir,
                '%s/metadata.yaml' % root_dir,
                '%s/README.md' % root_dir,
                '%s/data/ints.txt' % root_dir
            }

            observed = set(zf.namelist())

            # Not using self.assertArchiveMembers() because it accepts paths
            # relative to root_dir, and we have extra paths at the same level
            # as root_dir.
            self.assertEqual(observed, expected)

        archiver = Archiver.load(fp)

        self.assertEqual(archiver.uuid, self.archiver.uuid)
        self.assertEqual(archiver.type, IntSequence1)
        self.assertEqual(archiver.format, 'IntSequenceDirectoryFormat')
        self.assertIsNone(archiver.provenance)
        self.assertEqual({e for (e, _) in archiver.get_data_paths()},
                         {'ints.txt'})

    def test_load_ignores_directory_members(self):
        # Directory members aren't created by Python's zipfile module but can
        # be present if the archive is unzipped and then rezipped, for example,
        # using a command-line zip program.
        fp = os.path.join(self.temp_dir.name, 'archive.qza')
        self.archiver.save(fp)

        # Add directory entries to the archive.
        root_dir = str(self.archiver.uuid)
        with zipfile.ZipFile(fp, mode='a') as zf:
            zf.writestr('%s/' % root_dir, "")
            zf.writestr('%s/data/' % root_dir, "")

        # Assert the expected files exist in the archive to verify this test
        # case is testing what we want it to.
        expected = {
            '',  # Expected path: `root_dir`/
            'data/',
            'VERSION',
            'metadata.yaml',
            'README.md',
            'data/ints.txt'
        }

        self.assertArchiveMembers(fp, root_dir, expected)

        archiver = Archiver.load(fp)

        self.assertEqual(archiver.uuid, self.archiver.uuid)
        self.assertEqual(archiver.type, IntSequence1)
        self.assertEqual(archiver.format, 'IntSequenceDirectoryFormat')
        self.assertIsNone(archiver.provenance)
        self.assertEqual({e for (e, _) in archiver.get_data_paths()},
                         {'ints.txt'})

    def test_load_empty_archive(self):
        fp = os.path.join(self.temp_dir.name, 'empty.zip')

        with zipfile.ZipFile(fp, mode='w') as zf:
            pass

        with zipfile.ZipFile(fp, mode='r') as zf:
            expected = set()
            observed = set(zf.namelist())

            self.assertEqual(observed, expected)

        with self.assertRaisesRegex(ValueError, 'visible root directory'):
            Archiver.load(fp)

    def test_load_dotfile_only_archive(self):
        fp = os.path.join(self.temp_dir.name, 'dotfiles-only.zip')

        with zipfile.ZipFile(fp, mode='w') as zf:
            zf.writestr('.DS_Store', "The world's most beloved file\n")
            zf.writestr('.hidden-file',
                        "You can't see me if I can't see you\n")
            zf.writestr('.hidden-dir/ignored-file',
                        "I'm ignored because I live in a hidden dir :(\n")

        with zipfile.ZipFile(fp, mode='r') as zf:
            expected = {
                '.DS_Store',
                '.hidden-file',
                '.hidden-dir/ignored-file'
            }

            observed = set(zf.namelist())

            self.assertEqual(observed, expected)

        with self.assertRaisesRegex(ValueError, 'visible root directory'):
            Archiver.load(fp)

    def test_load_multiple_root_dirs(self):
        fp = os.path.join(self.temp_dir.name, 'multiple-root-dirs.zip')
        self.archiver.save(fp)

        # Add another semi-valid root dir.
        second_root_dir = str(uuid.uuid4())
        with zipfile.ZipFile(fp, mode='a') as zf:
            zf.writestr('%s/VERSION' % second_root_dir, "foo")

        with zipfile.ZipFile(fp, mode='r') as zf:
            root_dir = str(self.archiver.uuid)
            expected = {
                '%s/VERSION' % root_dir,
                '%s/metadata.yaml' % root_dir,
                '%s/README.md' % root_dir,
                '%s/data/ints.txt' % root_dir,
                '%s/VERSION' % second_root_dir
            }

            observed = set(zf.namelist())

            self.assertEqual(observed, expected)

        with self.assertRaisesRegex(ValueError, 'multiple root directories'):
            Archiver.load(fp)

    def test_load_invalid_uuid4_root_dir(self):
        fp = os.path.join(self.temp_dir.name, 'invalid-uuid4.zip')
        self.archiver.save(fp)

        # Invalid uuid4 taken from https://gist.github.com/ShawnMilo/7777304
        root_dir = '89eb3586-8a82-47a4-c911-758a62601cf7'
        with zipfile.ZipFile(fp, mode='w') as zf:
            self.archiver._save_version(zf, root_dir)
            self.archiver._save_readme(zf, root_dir)
            self.archiver._save_metadata(zf, root_dir)

            zf.writestr('%s/data/foo.txt' % root_dir, 'Hello, World!\n')

        with self.assertRaisesRegex(ValueError,
                                    'root directory.*valid version 4 UUID'):
            Archiver.load(fp)

    def test_load_root_dir_metadata_uuid_mismatch(self):
        fp = os.path.join(self.temp_dir.name, 'root-dir-metadata-mismatch.zip')
        self.archiver.save(fp)

        root_dir = str(uuid.uuid4())
        with zipfile.ZipFile(fp, mode='w') as zf:
            self.archiver._save_version(zf, root_dir)
            self.archiver._save_readme(zf, root_dir)
            self.archiver._save_metadata(zf, root_dir)

            zf.writestr('%s/data/foo.txt' % root_dir, 'Hello, World!\n')

        with self.assertRaisesRegex(
                ValueError, 'root directory must match UUID.*metadata'):
            Archiver.load(fp)

    def test_parse_uuid_valid(self):
        uuid_str = str(uuid.uuid4())

        obs = Archiver._parse_uuid(uuid_str)
        exp = uuid.UUID(hex=uuid_str, version=4)

        self.assertEqual(obs, exp)
        self.assertEqual(str(obs), uuid_str)

    def test_parse_uuid_invalid(self):
        # Invalid uuid4 taken from https://gist.github.com/ShawnMilo/7777304
        uuid_str = '89eb3586-8a82-47a4-c911-758a62601cf7'

        with self.assertRaisesRegex(ValueError, 'not a valid version 4 UUID'):
            Archiver._parse_uuid(uuid_str)

        # Not a UUID.
        uuid_str = 'abc123'
        with self.assertRaisesRegex(ValueError, 'not a valid version 4 UUID'):
            Archiver._parse_uuid(uuid_str)

        # Other UUID versions.
        for uuid_ in (uuid.uuid1(), uuid.uuid3(uuid.NAMESPACE_DNS, 'foo'),
                      uuid.uuid5(uuid.NAMESPACE_DNS, 'bar')):
            uuid_str = str(uuid_)
            with self.assertRaisesRegex(ValueError,
                                        'not a valid version 4 UUID'):
                Archiver._parse_uuid(uuid_str)

if __name__ == '__main__':
    unittest.main()
