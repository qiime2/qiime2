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

from flufl.lock import LockState

import qiime2
from qiime2.core.cache import Cache
from qiime2.core.testing.type import IntSequence1
from qiime2.sdk.result import Artifact


class TestCache(unittest.TestCase):
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

        self.assertEqual(Cache.base_cache_contents, contents)

        # Assert version file looks how we want
        with open(self.cache.version) as fh:
            lines = fh.readlines()

            self.assertEqual(lines[0], 'QIIME 2\n')
            self.assertEqual(lines[1],
                             f'cache: {self.cache.CURRENT_FORMAT_VERSION}\n')
            self.assertEqual(lines[2], f'framework: {qiime2.__version__}\n')

    def test_roundtrip(self):
        # Save artifact to cache
        art1 = Artifact.import_data(IntSequence1, [0, 1, 2])
        expected = art1.view(list)
        self.cache.save(art1, 'foo')

        # Delete artifact
        art1._archiver.path._destructor()
        del art1

        # Load artifact from cache
        art2 = self.cache.load('foo')

        # Ensure our data is correct
        self.assertEqual(expected, art2.view(list))

    def test_remove(self):
        # Save our artifact
        self.cache.save(self.art1, 'foo')

        # Show that we can load our artifact
        self.cache.load('foo')

        # remove our artifact
        self.cache.remove('foo')

        # Show that we can no longer load our artifact
        with self.assertRaisesRegex(FileNotFoundError,
                                    'No such file or directory'):
            self.cache.load('foo')

    def test_invalid_keys(self):
        # Invalid data key
        with self.assertRaisesRegex(ValueError, 'valid Python identifier'):
            self.cache.save(self.art1, '1')

        # Invalid pool key
        with self.assertRaisesRegex(ValueError, 'valid Python identifier'):
            self.cache.create_pool('1')

    def test_remove_locks(self):
        """Create some locks then see if we can remove them
        """
        test_pool = self.cache.create_pool('test')

        self.cache.lock.lock()
        test_pool.lock.lock()

        self.assertEqual(self.cache.lock.state, LockState.ours)
        self.assertEqual(test_pool.lock.state, LockState.ours)

        self.cache.clear_locks()

        self.assertEqual(self.cache.lock.state, LockState.unlocked)
        self.assertEqual(test_pool.lock.state, LockState.unlocked)

    # Might create another class for garbage collection tests to test more
    # cases with shared boilerplate
    def test_garbage_collection(self):
        # Data referenced directly by key
        self.cache.save(self.art1, 'foo')
        # Data referenced by pool that is referenced by key
        pool = self.cache.create_pool(['bar'])
        pool.save(self.art2)
        # We will be manually deleting the keys that back these two
        self.cache.save(self.art3, 'baz')
        pool = self.cache.create_pool(['qux'])
        pool.save(self.art4)

        # What we expect to see before and after gc
        expected_pre_gc_contents = \
            set(('./VERSION', 'keys/foo', 'keys/bar',
                 'keys/baz', 'keys/qux',
                 f'pools/bar/{self.art2.uuid}',
                 f'pools/qux/{self.art4.uuid}',
                 f'data/{self.art1.uuid}', f'data/{self.art2.uuid}',
                 f'data/{self.art3.uuid}', f'data/{self.art4.uuid}'))

        expected_post_gc_contents = \
            set(('./VERSION', 'keys/foo', 'keys/bar',
                 f'pools/bar/{self.art2.uuid}',
                 f'data/{self.art1.uuid}', f'data/{self.art2.uuid}'))

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
        """ Gets contents of cache not including contents of the artifacts
            themselves relative to the root of the cache
        """
        cache_contents = set()

        rel_keys = os.path.relpath(self.cache.keys, self.cache.path)
        rel_data = os.path.relpath(self.cache.data, self.cache.path)
        rel_pools = os.path.relpath(self.cache.pools, self.cache.path)
        rel_cache = os.path.relpath(self.cache.path, self.cache.path)

        for key in os.listdir(self.cache.keys):
            cache_contents.add(os.path.join(rel_keys, key))

        for art in os.listdir(self.cache.data):
            cache_contents.add(os.path.join(rel_data, art))

        for pool in os.listdir(self.cache.pools):
            for link in os.listdir(os.path.join(self.cache.pools, pool)):
                cache_contents.add(os.path.join(rel_pools, pool, link))

        for elem in os.listdir(self.cache.path):
            if os.path.isfile(os.path.join(self.cache.path, elem)):
                cache_contents.add(os.path.join(rel_cache, elem))

        return cache_contents
