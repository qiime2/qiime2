import codecs
import os
import shutil
import tempfile
import unittest
import zipfile

import qiime2
from qiime2 import Artifact
from qiime2.sdk.plugin_manager import PluginManager
from qiime2.core.archive.archiver import Archiver

from ..util import (
    _VERSION_MATCHER, parse_version, monkeypatch_archive_version,
    monkeypatch_framework_version, write_zip_archive
)


class TestVersionParser(unittest.TestCase):
    def setUp(self):
        self.pm = PluginManager()
        self.dp = self.pm.plugins['dummy-plugin']
        self.tempdir = tempfile.mkdtemp(
            prefix='qiime2-test-version-parser-temp-'
        )
        self.framework_version_exp = qiime2.__version__
        self.archive_version_exp = Archiver.CURRENT_FORMAT_VERSION

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_parse_version(self):
        int_seq = Artifact.import_data('IntSequence1', [1, 2, 3])
        int_seq.save(os.path.join(self.tempdir, 'int-seq.qza'))
        fp = os.path.join(self.tempdir, 'int-seq.qza')
        with zipfile.ZipFile(fp) as zf:
            actual = parse_version(zf)
            self.assertEqual(
                actual, (self.archive_version_exp, self.framework_version_exp)
            )

    def test_parse_version_old_archive_format(self):
        archive_version_exp = '2'
        fp = os.path.join(self.tempdir, 'int-seq-av2.qza')

        with monkeypatch_archive_version(archive_version_exp):
            int_seq = Artifact.import_data('IntSequence1', [1, 2, 3])
            int_seq.save(fp)

        with zipfile.ZipFile(fp) as zf:
            actual = parse_version(zf)
            self.assertEqual(
                actual, (archive_version_exp, self.framework_version_exp)
            )

    def test_artifact_with_commit_version(self):
        framework_version_exp = '2022.8.0+29.gb053440'
        fp = os.path.join(self.tempdir, 'int-seq-custom-fv.qza')

        with monkeypatch_framework_version(framework_version_exp):
            int_seq = Artifact.import_data('IntSequence1', [1, 2, 3])
            int_seq.save(fp)

        with zipfile.ZipFile(fp) as zf:
            actual = parse_version(zf)
            self.assertEqual(
                actual, (self.archive_version_exp, framework_version_exp)
            )

    def test_parse_version_no_VERSION_file(self):
        int_seq = Artifact.import_data('IntSequence1', [1, 2, 3])
        fp = os.path.join(self.tempdir, 'int-seq-no-v.qza')
        int_seq.save(fp)

        with tempfile.TemporaryDirectory() as tempdir:
            with zipfile.ZipFile(fp) as zf:
                zf.extractall(tempdir)

            uuid = os.listdir(tempdir)[0]
            os.remove(os.path.join(tempdir, uuid, 'VERSION'))
            write_zip_archive(fp, tempdir)

        with zipfile.ZipFile(fp) as zf:
            with self.assertRaisesRegex(ValueError,
                                        '(?s)VERSION.*nonexistent.*'):
                parse_version(zf)

    def test_parse_version_VERSION_file_missing_archive_field(self):
        int_seq = Artifact.import_data('IntSequence1', [1, 2, 3])
        fp = os.path.join(self.tempdir, 'int-seq-no-af.qza')
        int_seq.save(fp)

        with tempfile.TemporaryDirectory() as tempdir:
            with zipfile.ZipFile(fp) as zf:
                zf.extractall(tempdir)

            uuid = os.listdir(tempdir)[0]
            with open(os.path.join(tempdir, uuid, 'VERSION')) as fh:
                lines = fh.readlines()
                missing_archive_lines = [lines[0], lines[2]]

            with open(os.path.join(tempdir, uuid, 'VERSION'), 'w') as fh:
                for line in missing_archive_lines:
                    fh.write(line)

            write_zip_archive(fp, tempdir)

        with zipfile.ZipFile(fp) as zf:
            with self.assertRaisesRegex(ValueError, 'VERSION.*out of spec.*'):
                parse_version(zf)

    def test_parse_version_VERSION_file_extra_field(self):
        int_seq = Artifact.import_data('IntSequence1', [1, 2, 3])
        fp = os.path.join(self.tempdir, 'int-seq-extra-f.qza')
        int_seq.save(fp)

        with tempfile.TemporaryDirectory() as tempdir:
            with zipfile.ZipFile(fp) as zf:
                zf.extractall(tempdir)

            uuid = os.listdir(tempdir)[0]
            with open(os.path.join(tempdir, uuid, 'VERSION'), 'a') as fh:
                fh.write('fourth line\n')

            write_zip_archive(fp, tempdir)

        with zipfile.ZipFile(fp) as zf:
            with self.assertRaisesRegex(ValueError, 'VERSION.*out of spec.*'):
                parse_version(zf)

    '''
    Tests of the regex match itself below
    '''
    def test_version_too_short(self):
        short = (
            r'QIIME 2\n'
            r'archive: 4'
        )
        self.assertNotRegex(short, _VERSION_MATCHER)

    def test_version_too_long(self):
        long = (
            r'QIIME 2\n'
            r'archive: 4\n'
            r'framework: 2019.8.1.dev0\n'
            r'This line should not be here'
        )
        self.assertNotRegex(long, _VERSION_MATCHER)

    splitvm = codecs.decode(_VERSION_MATCHER.encode('utf-8'),
                            'unicode-escape').split(sep='\n')
    re_l1, re_l2, re_l3 = splitvm

    def test_line1_good(self):
        self.assertRegex('QIIME 2\n', self.re_l1)

    def test_line1_bad(self):
        self.assertNotRegex('SHIMMY 2\n', self.re_l1)

    def test_archive_version_1digit_numeric(self):
        self.assertRegex('archive: 1\n', self.re_l2)

    def test_archive_version_2digit_numeric(self):
        self.assertRegex('archive: 12\n', self.re_l2)

    def test_archive_version_bad(self):
        self.assertNotRegex('agama agama\n', self.re_l2)

    def test_archive_version_3digit_numeric(self):
        self.assertNotRegex('archive: 123\n', self.re_l2)

    def test_archive_version_nonnumeric(self):
        self.assertNotRegex('archive: 1a\n', self.re_l2)

    def test_fmwk_version_good_semver(self):
        self.assertRegex('framework: 2.0.6', self.re_l3)

    def test_fmwk_version_good_semver_dev(self):
        self.assertRegex('framework: 2.0.6.dev0', self.re_l3)

    def test_fmwk_version_good_year_month_patch(self):
        self.assertRegex('framework: 2020.2.0', self.re_l3)

    def test_fmwk_version_good_year_month_patch_2digit_month(self):
        self.assertRegex('framework: 2018.11.0', self.re_l3)

    def test_fmwk_version_good_year_month_patch_dev(self):
        self.assertRegex('framework: 2020.2.0.dev1', self.re_l3)

    def test_fmwk_version_good_ymp_2digit_month_dev(self):
        self.assertRegex('framework: 2020.11.0.dev0', self.re_l3)

    def test_fmwk_version_invalid_month(self):
        self.assertNotRegex('framework: 2020.13.0', self.re_l3)

    def test_fmwk_version_invalid_month_leading_zero(self):
        self.assertNotRegex('framework: 2020.03.0', self.re_l3)

    def test_fmwk_version_invalid_year(self):
        self.assertNotRegex('framework: 1953.3.0', self.re_l3)
