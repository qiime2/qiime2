# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pathlib
import unittest
import zipfile

import pytest

from .testing_utilities import CustomAssertions, DummyArtifacts
from ..util import get_root_uuid, get_nonroot_uuid


class GetRootUUIDTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.das = DummyArtifacts()
        cls.tempdir = cls.das.tempdir

    @classmethod
    def tearDownClass(cls):
        cls.das.free()

    @pytest.mark.filterwarnings('ignore::UserWarning')
    def test_get_root_uuid(self):
        exp_root_uuids = {
            '0': '89af91c0-033d-4e30-8ac4-f29a3b407dc1',
            '1': '5b929500-e4d6-4d3f-8f5f-93fd95d1117d',
            '2': 'e01f0484-40d4-420e-adcf-ca9be58ed1ee',
            '3': 'aa960110-4069-4b7c-97a3-8a768875e515',
            '4': '856502cb-66f2-45aa-a86c-e484cc9bfd57',
            '5': '48af8384-2b0a-4b26-b85c-11b79c0d6ea6',
            '6': '6facaf61-1676-45eb-ada0-d530be678b27',
        }
        for artifact, exp_uuid in zip(
            self.das.all_artifact_versions, exp_root_uuids.values()
        ):
            with zipfile.ZipFile(artifact.filepath) as zfh:
                self.assertEqual(exp_uuid, get_root_uuid(zfh))


class GetNonRootUUIDTests(unittest.TestCase):
    def test_get_nonroot_uuid(self):
        md_example = pathlib.Path(
            'arch_root/provenance/artifacts/uuid123/metadata.yaml')
        action_example = pathlib.Path(
            'arch_root/provenance/artifacts/uuid123/action/action.yaml')
        exp = 'uuid123'

        self.assertEqual(get_nonroot_uuid(md_example), exp)
        self.assertEqual(get_nonroot_uuid(action_example), exp)


class CustomAssertionsTests(CustomAssertions):
    def test_assert_re_appears_only_once(self):
        t = ("Lick an orange. It tastes like an orange.\n"
             "The strawberries taste like strawberries!\n"
             "The snozzberries taste like snozzberries!")
        self.assertREAppearsOnlyOnce(t, 'Lick an orange')
        self.assertREAppearsOnlyOnce(t, 'tastes like')
        with self.assertRaisesRegex(AssertionError, 'Regex.*match.*orange'):
            self.assertREAppearsOnlyOnce(t, 'orange')
        with self.assertRaisesRegex(AssertionError, 'Regex.*taste like'):
            self.assertREAppearsOnlyOnce(t, 'taste like')
        with self.assertRaisesRegex(AssertionError, 'Regex.*snozzberries'):
            self.assertREAppearsOnlyOnce(t, 'snozzberries')
        with self.assertRaisesRegex(AssertionError, 'Regex.*!'):
            self.assertREAppearsOnlyOnce(t, '!')
