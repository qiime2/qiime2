# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
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
from qiime2.core.archive.format.util import artifact_version
from qiime2.sdk import Artifact


class TestArtifactVersion(unittest.TestCase, ArchiveTestingMixin):
    def setUp(self):
        prefix = "qiime2-test-temp-"
        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)
        self.provenance_capture = archive.ImportProvenanceCapture()

    def test_nonexistent_archive_format(self):
        with self.assertRaisesRegex(ValueError, 'Version foo not supported'):
            with artifact_version('foo'):
                pass

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

    def test_write_v4_archive(self):
        fp = os.path.join(self.temp_dir.name, 'artifact_v1.qza')

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
