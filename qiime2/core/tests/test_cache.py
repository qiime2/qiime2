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
    base_cache_contents = set(('data', 'keys', 'pools', 'VERSION'))

    def setUp(self):
        # Create temp test dir
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')

        # Create artifact and cache
        self.art1 = Artifact.import_data(IntSequence1, [0, 1, 2])
        self.art2 = Artifact.import_data(IntSequence1, [3, 4, 5])
        self.art3 = Artifact.import_data(IntSequence1, [6, 7, 8])
        self.art4 = Artifact.import_data(IntSequence1, [9, 10, 11])
        self.cache = Cache(os.path.join(self.test_dir.name, 'new_cache'))

        self.not_cache_path = os.path.join(self.test_dir.name, 'not_cache')
        os.mkdir(self.not_cache_path)

    def tearDown(self):
        # Remove our cache and all that from last test
        self.test_dir.cleanup()

    # Verifies that is_cache is identifying a cache
    def test_is_cache(self):
        # Assert path is cache
        self.assertTrue(Cache.is_cache(self.cache.path))

    # Verifies that is_cache is identifying things aren't caches
    def test_is_not_cache(self):
        self.assertFalse(Cache.is_cache(self.not_cache_path))

    # This test manually asserts the cache created by the constructor looks
    # exactly as expected.
    def test_cache_manually_V1(self):
        self.assertTrue(os.path.exists(self.cache.path))
        contents = set(os.listdir(self.cache.path))

        self.assertEqual(self.base_cache_contents, contents)

        # Assert version file looks how we want
        with open(self.cache.version) as fh:
            lines = fh.readlines()

            self.assertEqual(lines[0], 'QIIME 2\n')
            self.assertEqual(lines[1],
                             f'cache: {self.cache.CURRENT_FORMAT_VERSION}\n')
            self.assertEqual(lines[2], f'framework: {qiime2.__version__}\n')

    def test_roundtrip(self):
        # Save artifact to cache
        self.cache.save(self.art1, 'foo')

        # Load artifact from cache
        art2 = self.cache.load('foo')

        # Ensure our data is correct
        self.assertEqual(self.art1.view(list), art2.view(list))

    def test_delete(self):
        # Save our artifact
        self.cache.save(self.art1, 'foo')

        # Show that we can load our artifact
        self.cache.load('foo')

        # delete our artifact
        self.cache.delete('foo')

        # Show that we can no longer load our artifact
        with self.assertRaisesRegex(FileNotFoundError,
                                    'No such file or directory'):
            self.cache.load('foo')

    def test_garbage_collection(self):
        # Data referenced directly by key
        self.cache.save(self.art1, 'foo')
        # Data referenced by pool that is referenced by key
        self.cache.save(self.art2, 'bar', complete=False)
        # We will be manually deleting the keys that back these two
        self.cache.save(self.art3, 'baz')
        self.cache.save(self.art4, 'qux', complete=False)

        # What we expect to see before and after gc
        expected_pre_gc_contents = \
            set(('./VERSION', 'keys/foo', 'keys/bar',
                 'keys/baz', 'keys/qux',
                 f'pools/{self.art2.uuid}/{self.art2.uuid}',
                 f'pools/{self.art4.uuid}/{self.art4.uuid}',
                 f'data/{self.art1.uuid}.qza', f'data/{self.art2.uuid}.qza',
                 f'data/{self.art3.uuid}.qza', f'data/{self.art4.uuid}.qza'))

        expected_post_gc_contents = \
            set(('./VERSION', 'keys/foo', 'keys/bar',
                 f'pools/{self.art2.uuid}/{self.art2.uuid}',
                 f'data/{self.art1.uuid}.qza', f'data/{self.art2.uuid}.qza'))

        # Assert cache looks how we want pre gc
        pre_gc_contents = self.get_cache_contents()
        self.assertEqual(expected_pre_gc_contents, pre_gc_contents)

        # Delete keys
        os.remove(self.cache.keys / 'baz')
        os.remove(self.cache.keys / 'qux')

        # Run gc
        self.cache.garbage_collection()

        # Assert cache looks how we want post gc
        post_gc_contents = self.get_cache_contents()
        self.assertEqual(expected_post_gc_contents, post_gc_contents)

    def get_cache_contents(self):
        cache_contents = set()

        for dir, _, files in os.walk(self.cache.path):
            for file in files:
                rel_dir = os.path.relpath(dir, self.cache.path)
                rel_file = os.path.join(rel_dir, file)
                cache_contents.add(rel_file)

        return cache_contents
