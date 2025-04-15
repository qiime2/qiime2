# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest
import tempfile
import os
import zipfile

from qiime2.core.testing.type import FourInts
from qiime2.core.testing.util import ArchiveTestingMixin
import qiime2.core.archive as archive
from qiime2.core.archive.format.v7_0 import ArchiveFormat
from qiime2.core.archive.format.util import artifact_version
from qiime2.core.annotate import Note
from qiime2.sdk import Artifact


class TestArtifactVersion(unittest.TestCase, ArchiveTestingMixin):
    def setUp(self):
        prefix = "qiime2-test-temp-"
        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)
        self.provenance_capture = archive.ImportProvenanceCapture()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_nonexistent_archive_format(self):
        with self.assertRaisesRegex(ValueError, 'Version foo not supported'):
            with artifact_version('foo'):
                pass

    # ARCHIVE V0
    def test_write_v0_archive(self):
        fp = os.path.join(self.temp_dir.name, 'artifact_v0.qza')

        with artifact_version(0):
            artifact = Artifact._from_view(FourInts, [-1, 42, 0, 43], list,
                                           self.provenance_capture)
            artifact.save(fp)

        root_dir = str(artifact.uuid)
        # There should be no provenance
        expected = {
            'VERSION',
            'metadata.yaml',
            'data/file1.txt',
            'data/file2.txt',
            'data/nested/file3.txt',
            'data/nested/file4.txt',
        }
        self.assertArchiveMembers(fp, root_dir, expected)

        with zipfile.ZipFile(fp, mode='r') as zf:
            version = zf.read(os.path.join(root_dir, 'VERSION'))
        self.assertRegex(str(version), '^.*archive: 0.*$')

    # ARCHIVE V4
    def test_write_v4_archive(self):
        fp = os.path.join(self.temp_dir.name, 'artifact_v4.qza')

        with artifact_version(4):
            artifact = Artifact._from_view(FourInts, [-1, 42, 0, 43], list,
                                           self.provenance_capture)
            artifact.save(fp)

        root_dir = str(artifact.uuid)
        expected = {
            'VERSION',
            'metadata.yaml',
            'data/file1.txt',
            'data/file2.txt',
            'data/nested/file3.txt',
            'data/nested/file4.txt',
            'provenance/metadata.yaml',
            'provenance/VERSION',
            'provenance/citations.bib',
            'provenance/action/action.yaml',
        }
        self.assertArchiveMembers(fp, root_dir, expected)

        with zipfile.ZipFile(fp, mode='r') as zf:
            version = zf.read(os.path.join(root_dir, 'VERSION'))
        self.assertRegex(str(version), '^.*archive: 4.*$')

    # ARCHIVE V7.0
    def test_write_v7_0_archive_with_note_annotation(self):
        fp = os.path.join(self.temp_dir.name, 'artifact_v7_0.qza')

        with artifact_version(7.0):
            artifact = Artifact._from_view(FourInts, [-1, 42, 0, 43], list,
                                           self.provenance_capture)
            note = Note(name='mynote', text='my special text')
            artifact.add_annotation(note)
            artifact.save(fp)

        root_dir = str(artifact.uuid)
        expected = {
            'VERSION',
            'metadata.yaml',
            'checksums.md5',
            'data/file1.txt',
            'data/file2.txt',
            'data/nested/file3.txt',
            'data/nested/file4.txt',
            'provenance/metadata.yaml',
            'provenance/VERSION',
            'provenance/citations.bib',
            'provenance/conda-env.yaml',
            'provenance/action/action.yaml',
            f'provenance/annotations/{note.uuid}/metadata.yaml',
            f'provenance/annotations/{note.uuid}/note.txt'
        }
        self.assertArchiveMembers(fp, root_dir, expected)

        with zipfile.ZipFile(fp, mode='r') as zf:
            version = zf.read(os.path.join(root_dir, 'VERSION'))
            metadata = zf.read(os.path.join(root_dir, 'metadata.yaml'))
            conda_env = \
                zf.read(os.path.join(root_dir, 'provenance', 'conda-env.yaml'))
            annotation_metadata = \
                zf.read(os.path.join(root_dir, 'provenance', 'annotations',
                                     f'{note.uuid}', 'metadata.yaml'))
            note_contents = \
                zf.read(os.path.join(root_dir, 'provenance', 'annotations',
                                     f'{note.uuid}', 'note.txt'))
        self.assertRegex(str(version), '^.*archive: 7.0.*$')
        self.assertRegex(str(metadata), '^.*data-size: .*B.*$')
        self.assertRegex(str(conda_env), '^.*dependencies:.*- .*$')
        self.assertRegex(str(annotation_metadata),
                         f'^.*id: {note.uuid}.*name: mynote.*type: Note.*$')
        self.assertRegex(str(note_contents), 'my special text')

    # testing file size conversion helper
    def test_human_readable_size_util(self):
        cases = [
            (0, "0.0 B"),
            (500, "500.0 B"),
            (1023, "1023.0 B"),
            (1024, "1.0 KiB"),
            (1536, "1.5 KiB"),
            (1024**2, "1.0 MiB"),
        ]

        for num, expected in cases:
            result = ArchiveFormat._human_readable_size(num)
            self.assertEqual(result, expected)
