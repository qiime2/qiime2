import codecs
import os
import unittest
import warnings
import zipfile

from .test_parse import DATA_DIR, TEST_DATA
from ..version_parser import (
    _VERSION_MATCHER, parse_version, parse_version_from_fp
)


class ParseVersionFromFPTests(unittest.TestCase):
    """
    Simple tests for from_fp convenience function, which passes an fp through
    to `parse_version`. Comprehensive tests for that function below.
    """

    def test_parse_version_from_fp(self):
        fp = os.path.join(DATA_DIR, 'v5_uu_emperor.qzv')
        actual = parse_version_from_fp(fp)
        self.assertEqual(actual,
                         (TEST_DATA['5']['av'],
                          TEST_DATA['5']['fwv']))

    def test_qza_written_in_dev_state(self):
        """
        The test archive was produced by running feature-table merge in a dev
        version of QIIME 2.
        """
        fp = os.path.join(DATA_DIR, 'table_written_in_dev_version.qza')
        actual = parse_version_from_fp(fp)
        self.assertEqual(actual,
                         ('5', '2021.10.0.dev0'))

    def test_qza_written_in_dirty_dev_state(self):
        """
        The test archive was produced in a "dirty" dev version of QIIME 2.
        The regex handles these generously, matching a plus followed by
        any collection of . and word characters
        """
        fp = os.path.join(DATA_DIR, 'ns_collisions.qza')
        actual = parse_version_from_fp(fp)
        self.assertEqual(actual,
                         ('5', '2019.7.0.dev0+111.g40f412f'))


class GetVersionTests(unittest.TestCase):
    v5_no_version = os.path.join(DATA_DIR, 'VERSION_missing.qzv')
    v5_qzv_version_bad = os.path.join(DATA_DIR, 'VERSION_bad.qzv')
    v5_qzv_version_short = os.path.join(DATA_DIR, 'VERSION_short.qzv')
    v5_qzv_version_long = os.path.join(DATA_DIR, 'VERSION_long.qzv')
    uuid = '8854f06a-872f-4762-87b7-4541d0f283d4'

    # High-level checks only. Detailed tests of the VERSION_MATCHER regex are
    # in test_archive_formats.VersionMatcherTests to reduce overhead

    def test_parse_version_no_VERSION_file(self):
        with zipfile.ZipFile(self.v5_no_version) as zf:
            fn = 'VERSION_missing.qzv'
            with self.assertRaisesRegex(
                    ValueError, f'(?s)VERSION.*nonexistent.*{fn}'):
                parse_version(zf)

    def test_parse_version_VERSION_bad(self):
        with zipfile.ZipFile(self.v5_qzv_version_bad) as zf:
            fn = 'VERSION_bad.qzv'
            with self.assertRaisesRegex(
                    ValueError, f'VERSION.*out of spec.*{fn}'):
                parse_version(zf)

    def test_short_VERSION(self):
        with zipfile.ZipFile(self.v5_qzv_version_short) as zf:
            fn = 'VERSION_short.qzv'
            with self.assertRaisesRegex(
                    ValueError, f'VERSION.*out of spec.*{fn}'):
                parse_version(zf)

    def test_long_VERSION(self):
        fn = 'VERSION_long.qzv'
        with zipfile.ZipFile(self.v5_qzv_version_long) as zf:
            with self.assertRaisesRegex(
                    ValueError, f'VERSION.*out of spec.*{fn}'):
                parse_version(zf)

    def test_version_nums(self):
        for arch_ver in TEST_DATA:
            qzv = os.path.join(DATA_DIR, 'v' + arch_ver + '_uu_emperor.qzv')
            with zipfile.ZipFile(qzv) as zf:
                act_arch, act_frmwk = parse_version(zf)
                self.assertEqual(act_arch, TEST_DATA[arch_ver]['av'])
                self.assertEqual(act_frmwk, TEST_DATA[arch_ver]['fwv'])


class ArchiveVersionMatcherTests(unittest.TestCase):
    """Testing for the _VERSION_MATCHER regex itself"""

    def test_version_too_short(self):
        shorty = (
            r'QIIME 2\n'
            r'archive: 4'
        )
        self.assertNotRegex(shorty, _VERSION_MATCHER)

    def test_version_too_long(self):
        longy = (
            r'QIIME 2\n'
            r'archive: 4\n'
            r'framework: 2019.8.1.dev0\n'
            r'This line should not be here'
        )
        self.assertNotRegex(longy, _VERSION_MATCHER)

    warnings.filterwarnings('ignore', 'invalid escape sequence',
                            DeprecationWarning)
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
