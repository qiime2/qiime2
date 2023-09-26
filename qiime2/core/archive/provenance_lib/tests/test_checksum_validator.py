import os
import shutil
import tempfile
import unittest
import zipfile

import pytest

from qiime2 import Artifact
from qiime2.core.archive.archiver import ChecksumDiff
from qiime2.sdk.plugin_manager import PluginManager

from .._checksum_validator import validate_checksums, ValidationCode
from .testing_utilities import write_zip_archive


class ValidateChecksumTests(unittest.TestCase):
    def setUp(self):
        self.pm = PluginManager()
        self.dp = self.pm.plugins['dummy-plugin']
        self.tempdir = tempfile.mkdtemp(
            prefix='qiime2-test-checksum-validator-temp-'
        )

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_validate_checksums(self):
        int_seq = Artifact.import_data('IntSequence1', [1, 2, 3])
        fp = os.path.join(self.tempdir, 'int-seq.qza')
        int_seq.save(fp)

        with zipfile.ZipFile(fp) as zf:
            is_valid, diff = validate_checksums(zf)
        self.assertEqual(is_valid, ValidationCode.VALID)
        self.assertEqual(diff, ChecksumDiff({}, {}, {}))

    @pytest.mark.filterwarnings('ignore::UserWarning')
    def test_validate_checksums_invalid(self):
        '''
        Mangle an intact v5 Archive so that its checksums.md5 is invalid,
        and then confirm that we're catching all the changes we've made
        Specifically:
        - remove the root `<uuid>/metadata.yaml`
        - add a new file called '<uuid>/tamper.txt`
        - overwrite `<uuid>/provenance/citations.bib`
        '''
        int_seq = Artifact.import_data('IntSequence1', [1, 2, 3])
        fp = os.path.join(self.tempdir, 'int-seq-altered.qza')
        int_seq.save(fp)

        with tempfile.TemporaryDirectory() as tempdir:
            with zipfile.ZipFile(fp) as zf:
                zf.extractall(tempdir)

            uuid = os.listdir(tempdir)[0]
            root_dir = os.path.join(tempdir, uuid)
            print(os.listdir(root_dir))
            os.remove(os.path.join(root_dir, 'metadata.yaml'))
            print(os.listdir(root_dir))
            with open(os.path.join(root_dir, 'tamper.txt'), 'w') as fh:
                pass
            citations_path = \
                os.path.join(root_dir, 'provenance', 'citations.bib')
            with open(citations_path, 'w') as fh:
                fh.write('file overwritten\n')

            write_zip_archive(fp, tempdir)

        with zipfile.ZipFile(fp) as zf:
            is_valid, diff = validate_checksums(zf)

        self.assertEqual(is_valid, ValidationCode.INVALID)
        self.assertEqual(list(diff.added.keys()), ['tamper.txt'])
        self.assertEqual(list(diff.removed.keys()), ['metadata.yaml'])
        self.assertEqual(list(diff.changed.keys()),
                         ['provenance/citations.bib'])

    @pytest.mark.filterwarnings('ignore::UserWarning')
    def test_validate_checksums_checksums_missing(self):
        int_seq = Artifact.import_data('IntSequence1', [1, 2, 3])
        fp = os.path.join(self.tempdir, 'int-seq-missing-version.qza')
        int_seq.save(fp)

        with tempfile.TemporaryDirectory() as tempdir:
            with zipfile.ZipFile(fp) as zf:
                zf.extractall(tempdir)

            uuid = os.listdir(tempdir)[0]
            os.remove(os.path.join(tempdir, uuid, 'checksums.md5'))

            write_zip_archive(fp, tempdir)

        with zipfile.ZipFile(fp) as zf:
            is_valid, diff = validate_checksums(zf)

        self.assertEqual(is_valid, ValidationCode.INVALID)
        self.assertEqual(diff, None)
