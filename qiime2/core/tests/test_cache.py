# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import unittest

import qiime2
from qiime2.core.cache import Cache
from qiime2.core.testing.type import IntSequence1
from qiime2.sdk.result import Artifact


class TestCache(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')

    def tearDown(self):
        self.test_dir.cleanup()

    def test_existing_cache(self):
        # Assert path exists
        # self.assertTrue()
        # Assert path is cache
        pass

    # This test was written for version 1 of the cache
    def test_create_cache(self):
        cache_contents = set(('data', 'keys', 'pools', 'VERSION'))
        path = os.path.join(self.test_dir.name, 'new_cache')

        # Assert path does not exist
        self.assertFalse(os.path.exists(path))
        # Create cache at path
        cache = Cache(path)

        # Assert cache looks how we expect
        self.assertTrue(os.path.exists(path))
        contents = os.listdir(path)
        dedup = set(contents)

        self.assertEqual(len(contents), len(dedup))
        self.assertEqual(cache_contents, dedup)

        # Assert version file looks how we want
        with open(cache.version) as fh:
            lines = fh.readlines()

            self.assertEqual(lines[0], 'QIIME 2\n')
            self.assertEqual(lines[1],
                             f'cache: {cache.CURRENT_FORMAT_VERSION}\n')
            self.assertEqual(lines[2], f'framework: {qiime2.__version__}\n')

    def test_not_a_cache(self):
        # Assert path exists
        # Assert path is not cache
        pass

    def test_roundtrip(self):
        # Create artifact and cache
        art = Artifact.import_data(IntSequence1, [0, 1, 2])
        cache = Cache(os.path.join(self.test_dir.name, 'new_cache'))

        # Save artifact to cache
        cache.save(art, 'foo')

        # Load artifact from cache
        art2 = cache.load('foo')

        # Ensure our data is correct
        self.assertEqual(art.view(list), art2.view(list))
