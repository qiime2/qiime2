# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import unittest
import uuid
import zipfile
import pathlib

from qiime2.core.archive import Archiver
from qiime2.core.archive import ImportProvenanceCapture
from qiime2.core.archive.archiver import _ZipArchive
from qiime2.core.archive.format.util import artifact_version
from qiime2.core.testing.format import IntSequenceDirectoryFormat
from qiime2.core.testing.type import IntSequence1
from qiime2.core.testing.util import ArchiveTestingMixin


class TestArchiver(unittest.TestCase, ArchiveTestingMixin):
    def setUp(self):
        prefix = "qiime2-test-temp-"
        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)

        # Initialize an Archiver. The values passed to the constructor mostly
        # don't matter to the Archiver, but we'll pass valid Artifact test data
        # anyways in case Archiver's behavior changes in the future.
        def data_initializer(data_dir):
            fp = os.path.join(str(data_dir), 'ints.txt')
            with open(fp, 'w') as fh:
                fh.write('1\n')
                fh.write('2\n')
                fh.write('3\n')

        self.archiver = Archiver.from_data(
            IntSequence1, IntSequenceDirectoryFormat,
            data_initializer=data_initializer,
            provenance_capture=ImportProvenanceCapture())

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_save_invalid_filepath(self):
        # Empty filepath.
        with self.assertRaisesRegex(FileNotFoundError, 'No such file'):
            self.archiver.save('')

        # Directory.
        with self.assertRaisesRegex(IsADirectoryError, 'directory'):
            self.archiver.save(self.temp_dir.name)

        # Ends with path separator (no basename, e.g. /tmp/foo/).
        with self.assertRaises((IsADirectoryError, FileNotFoundError)):
            self.archiver.save(os.path.join(self.temp_dir.name, 'foo', ''))

    def test_save_excludes_dotfiles_in_data_dir(self):
        def data_initializer(data_dir):
            data_dir = str(data_dir)
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

        archiver = Archiver.from_data(
            IntSequence1, IntSequenceDirectoryFormat,
            data_initializer=data_initializer,
            provenance_capture=ImportProvenanceCapture())

        fp = os.path.join(self.temp_dir.name, 'archive.zip')
        archiver.save(fp)

        root_dir = str(archiver.uuid)
        expected = {
            'VERSION',
            'checksums.md5',
            'metadata.yaml',
            'data/ints.txt',
            'provenance/metadata.yaml',
            'provenance/VERSION',
            'provenance/citations.bib',
            'provenance/action/action.yaml'
        }

        self.assertArchiveMembers(fp, root_dir, expected)

    def test_save_archive_members(self):
        fp = os.path.join(self.temp_dir.name, 'archive.zip')

        self.archiver.save(fp)

        root_dir = str(self.archiver.uuid)
        expected = {
            'VERSION',
            'checksums.md5',
            'metadata.yaml',
            'data/ints.txt',
            'provenance/metadata.yaml',
            'provenance/VERSION',
            'provenance/citations.bib',
            'provenance/action/action.yaml'
        }

        self.assertArchiveMembers(fp, root_dir, expected)

    def test_load_archive(self):
        fp = os.path.join(self.temp_dir.name, 'archive.zip')
        self.archiver.save(fp)

        archiver = Archiver.load(fp)

        self.assertEqual(archiver.uuid, self.archiver.uuid)
        self.assertEqual(archiver.type, IntSequence1)
        self.assertEqual(archiver.format, IntSequenceDirectoryFormat)
        self.assertEqual({str(p.relative_to(archiver.data_dir))
                          for p in archiver.data_dir.iterdir()},
                         {'ints.txt'})

    def test_load_ignores_root_dotfiles(self):
        fp = os.path.join(self.temp_dir.name, 'archive.zip')
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
                '%s/checksums.md5' % root_dir,
                '%s/metadata.yaml' % root_dir,
                '%s/data/ints.txt' % root_dir,
                '%s/provenance/metadata.yaml' % root_dir,
                '%s/provenance/VERSION' % root_dir,
                '%s/provenance/citations.bib' % root_dir,
                '%s/provenance/action/action.yaml' % root_dir
            }

            observed = set(zf.namelist())

            # Not using self.assertArchiveMembers() because it accepts paths
            # relative to root_dir, and we have extra paths at the same level
            # as root_dir.
            self.assertEqual(observed, expected)

        archiver = Archiver.load(fp)

        self.assertEqual(archiver.uuid, self.archiver.uuid)
        self.assertEqual(archiver.type, IntSequence1)
        self.assertEqual(archiver.format, IntSequenceDirectoryFormat)
        self.assertEqual({str(p.relative_to(archiver.data_dir))
                          for p in archiver.data_dir.iterdir()},
                         {'ints.txt'})

    def test_load_ignores_directory_members(self):
        # Directory members aren't created by Python's zipfile module but can
        # be present if the archive is unzipped and then rezipped, for example,
        # using a command-line zip program.
        fp = os.path.join(self.temp_dir.name, 'archive.zip')
        self.archiver.save(fp)

        # Add directory entries to the archive.
        root_dir = str(self.archiver.uuid)
        with zipfile.ZipFile(fp, mode='a') as zf:
            zf.writestr('%s/' % root_dir, "")
            zf.writestr('%s/data/' % root_dir, "")
            zf.writestr('%s/data/nested/' % root_dir, "")
            zf.writestr('%s/data/nested/foo.txt' % root_dir, "bar")

        # Assert the expected files exist in the archive to verify this test
        # case is testing what we want it to.
        expected = {
            '',  # Expected path: `root_dir`/
            'data/',
            'data/nested/',
            'VERSION',
            'checksums.md5',
            'metadata.yaml',
            'data/ints.txt',
            'data/nested/foo.txt',
            'provenance/metadata.yaml',
            'provenance/VERSION',
            'provenance/citations.bib',
            'provenance/action/action.yaml'
        }

        self.assertArchiveMembers(fp, root_dir, expected)

        archiver = Archiver.load(fp)

        self.assertEqual(archiver.uuid, self.archiver.uuid)
        self.assertEqual(archiver.type, IntSequence1)
        self.assertEqual(archiver.format, IntSequenceDirectoryFormat)

        archiver.save(fp)

        root_dir = str(archiver.uuid)
        expected = {
            # Directory entries should not be present.
            'VERSION',
            'checksums.md5',
            'metadata.yaml',
            'data/ints.txt',
            'data/nested/foo.txt',
            'provenance/metadata.yaml',
            'provenance/VERSION',
            'provenance/citations.bib',
            'provenance/action/action.yaml'
        }

        self.assertArchiveMembers(fp, root_dir, expected)

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
                '%s/checksums.md5' % root_dir,
                '%s/metadata.yaml' % root_dir,
                '%s/data/ints.txt' % root_dir,
                '%s/provenance/metadata.yaml' % root_dir,
                '%s/provenance/VERSION' % root_dir,
                '%s/provenance/citations.bib' % root_dir,
                '%s/provenance/action/action.yaml' % root_dir,
                '%s/VERSION' % second_root_dir
            }

            observed = set(zf.namelist())

            self.assertEqual(observed, expected)

        with self.assertRaisesRegex(ValueError, 'multiple root directories'):
            Archiver.load(fp)

    def test_load_invalid_uuid4_root_dir(self):
        fp = pathlib.Path(self.temp_dir.name) / 'invalid-uuid4'
        zp = pathlib.Path(self.temp_dir.name) / 'bad.zip'
        fp.mkdir()
        # Invalid uuid4 taken from https://gist.github.com/ShawnMilo/7777304
        root_dir = '89eb3586-8a82-47a4-c911-758a62601cf7'

        record = _ZipArchive.setup(fp, 'foo', 'bar')
        (fp / str(record.uuid)).rename(fp / root_dir)
        _ZipArchive.save(fp, zp)

        with self.assertRaisesRegex(ValueError,
                                    'root directory.*valid version 4 UUID'):
            _ZipArchive(zp)

    def test_is_uuid4_valid(self):
        uuid_str = str(uuid.uuid4())

        self.assertTrue(_ZipArchive._is_uuid4(uuid_str))

    def test_parse_uuid_invalid(self):
        # Invalid uuid4 taken from https://gist.github.com/ShawnMilo/7777304
        uuid_str = '89eb3586-8a82-47a4-c911-758a62601cf7'
        self.assertFalse(_ZipArchive._is_uuid4(uuid_str))

        # Not a UUID.
        uuid_str = 'abc123'
        self.assertFalse(_ZipArchive._is_uuid4(uuid_str))

        # Other UUID versions.
        for uuid_ in (uuid.uuid1(), uuid.uuid3(uuid.NAMESPACE_DNS, 'foo'),
                      uuid.uuid5(uuid.NAMESPACE_DNS, 'bar')):
            uuid_str = str(uuid_)
            self.assertFalse(_ZipArchive._is_uuid4(uuid_str))

    def test_checksums_match(self):
        diff = self.archiver.validate_checksums()

        self.assertEqual(diff.added, {})
        self.assertEqual(diff.removed, {})
        self.assertEqual(diff.changed, {})

    def test_checksums_mismatch(self):
        with (self.archiver.root_dir / 'data' / 'ints.txt').open('w') as fh:
            fh.write('999\n')
        with (self.archiver.root_dir / 'tamper.txt').open('w') as fh:
            fh.write('extra file')

        (self.archiver.root_dir / 'VERSION').unlink()

        diff = self.archiver.validate_checksums()

        self.assertEqual(diff.added,
                         {'tamper.txt': '296583001b00d2b811b5871b19e0ad28'})
        # The contents of most files is either stochastic, or has the current
        # version (which is an unknown commit sha1), so just check name
        self.assertEqual(list(diff.removed.keys()), ['VERSION'])
        self.assertEqual(diff.changed,
                         {'data/ints.txt': ('c0710d6b4f15dfa88f600b0e6b624077',
                                            'f47bc36040d5c7db08e4b3a457dcfbb2')
                          })

    def test_checksum_backwards_compat(self):
        self.tearDown()
        with artifact_version(4):
            self.setUp()

        diff = self.archiver.validate_checksums()

        self.assertEqual(diff.added, {})
        self.assertEqual(diff.removed, {})
        self.assertEqual(diff.changed, {})


if __name__ == '__main__':
    unittest.main()
