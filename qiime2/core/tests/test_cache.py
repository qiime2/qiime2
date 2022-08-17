# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import atexit
import tempfile
import unittest

from flufl.lock import LockState

import qiime2
from qiime2.core.cache import Cache, _exit_cleanup, get_cache
from qiime2.core.testing.type import IntSequence1, IntSequence2
from qiime2.core.testing.util import get_dummy_plugin
from qiime2.sdk.result import Artifact

TEST_POOL = '__TEST_POOL__'


# TODO: Check process contents too
def _get_cache_contents(cache):
    """Gets contents of cache not including contents of the artifacts
    themselves relative to the root of the cache
    """
    cache_contents = set()

    rel_keys = os.path.relpath(cache.keys, cache.path)
    rel_data = os.path.relpath(cache.data, cache.path)
    rel_pools = os.path.relpath(cache.pools, cache.path)
    rel_cache = os.path.relpath(cache.path, cache.path)

    for key in os.listdir(cache.keys):
        cache_contents.add(os.path.join(rel_keys, key))

    for art in os.listdir(cache.data):
        cache_contents.add(os.path.join(rel_data, art))

    for pool in os.listdir(cache.pools):
        for link in os.listdir(os.path.join(cache.pools, pool)):
            cache_contents.add(os.path.join(rel_pools, pool, link))

    for elem in os.listdir(cache.path):
        if os.path.isfile(os.path.join(cache.path, elem)):
            cache_contents.add(os.path.join(rel_cache, elem))

    return cache_contents


def _on_exit_validate(cache, expected):
    observed = _get_cache_contents(cache)
    cache.remove(TEST_POOL)
    assert expected.issubset(observed)


class TestCache(unittest.TestCase):
    def setUp(self):
        # Create temp test dir
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')

        # Create artifact and cache
        self.art1 = Artifact.import_data(IntSequence1, [0, 1, 2])
        self.art2 = Artifact.import_data(IntSequence1, [3, 4, 5])
        self.art3 = Artifact.import_data(IntSequence1, [6, 7, 8])
        self.art4 = Artifact.import_data(IntSequence2, [9, 10, 11])
        self.cache = Cache(os.path.join(self.test_dir.name, 'new_cache'))

        self.not_cache_path = os.path.join(self.test_dir.name, 'not_cache')
        os.mkdir(self.not_cache_path)

    def tearDown(self):
        """Remove our cache and all that from last test
        """
        self.test_dir.cleanup()

    def test_is_cache(self):
        """Verifies that is_cache is identifying a cache
        """
        self.assertTrue(Cache.is_cache(self.cache.path))

    def test_is_not_cache(self):
        """Verifies that is_cache is identifying when things aren't caches
        """
        self.assertFalse(Cache.is_cache(self.not_cache_path))

    def test_cache_manually_V1(self):
        """This test manually asserts the cache created by the constructor
        looks exactly as expected.
        """
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
        self.cache.lock.lock()

        self.assertEqual(self.cache.lock.state, LockState.ours)

        self.cache.clear_lock()

        self.assertEqual(self.cache.lock.state, LockState.unlocked)

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
                 f'data/{self.art1.uuid}', f'data/{self.art2.uuid}',
                 f'data/{self.art3.uuid}', f'data/{self.art4.uuid}'))

        # Assert cache looks how we want pre gc
        pre_gc_contents = _get_cache_contents(self.cache)
        self.assertEqual(expected_pre_gc_contents, pre_gc_contents)

        # Delete keys
        self.cache.remove(self.cache.keys / 'baz')
        self.cache.remove(self.cache.keys / 'qux')

        # Assert cache looks how we want post gc
        post_gc_contents = _get_cache_contents(self.cache)
        self.assertEqual(expected_post_gc_contents, post_gc_contents)

    def test_asynchronous(self):
        plugin = get_dummy_plugin()
        concatenate_ints = plugin.methods['concatenate_ints']

        with self.cache:
            future = concatenate_ints.asynchronous(self.art1, self.art2,
                                                   self.art4, 4, 5)
            result = future.result()

        result = result[0]

        expected = set(('./VERSION', f'data/{result._archiver.uuid}'))

        observed = _get_cache_contents(self.cache)
        self.assertEqual(expected, observed)

    def test_asynchronous_pool(self):
        plugin = get_dummy_plugin()
        concatenate_ints = plugin.methods['concatenate_ints']
        test_pool = self.cache.create_pool(keys=[TEST_POOL])

        with self.cache:
            with test_pool:
                future = concatenate_ints.asynchronous(self.art1, self.art2,
                                                       self.art4, 4, 5)
                result = future.result()

        result = result[0]

        expected = set((
            './VERSION', f'data/{result._archiver.uuid}', 'keys/__TEST_POOL__',
            f'pools/__TEST_POOL__/{result._archiver.uuid}'
        ))

        observed = _get_cache_contents(self.cache)
        self.assertEqual(expected, observed)

    def test_no_dangling_ref(self):
        ref = self.cache.save(self.art1, 'foo')
        ref.validate()

        # This would create a dangling ref if we were not properly saving
        # things to the process pool when we save them with a cache key
        self.cache.remove('foo')
        ref.validate()

    def test_no_dangling_ref_pool(self):
        pool = self.cache.create_pool(keys=['pool'])
        ref = pool.save(self.art1)
        ref.validate()

        # This would create a dangling ref if we were not properly saving
        # things to the process pool when we save them to a named pool
        self.cache.remove('pool')
        ref.validate()

    def test_pool(self):
        pool = self.cache.create_pool(keys=['pool'])

        # Create an artifact in the cache and the pool
        with self.cache:
            with pool:
                ref = Artifact.import_data(IntSequence1, [0, 1, 2])

        uuid = str(ref.uuid)
        self.assertIn(uuid, os.listdir(self.cache.data))
        self.assertIn(uuid, os.listdir(self.cache.pools / 'pool'))

    def test_pool_no_cache_set(self):
        pool = self.cache.create_pool(keys=['pool'])

        with pool:
            ref = Artifact.import_data(IntSequence1, [0, 1, 2])

        uuid = str(ref.uuid)
        self.assertIn(uuid, os.listdir(self.cache.data))
        self.assertIn(uuid, os.listdir(self.cache.pools / 'pool'))

    def test_pool_wrong_cache_set(self):
        cache = Cache(os.path.join(self.test_dir.name, 'cache'))
        pool = self.cache.create_pool(keys=['pool'])

        with cache:
            with self.assertRaisesRegex(ValueError,
                                        'pool that is not on cache'):
                with pool:
                    Artifact.import_data(IntSequence1, [0, 1, 2])

    def test_enter_multiple_caches(self):
        cache = Cache(os.path.join(self.test_dir.name, 'cache'))

        with self.cache:
            with self.assertRaisesRegex(ValueError,
                                        'cannot enter multiple caches.*'
                                        f'{self.cache.path}'):
                with cache:
                    pass

    def test_enter_multiple_pools(self):
        pool1 = self.cache.create_pool(keys=['pool1'])
        pool2 = self.cache.create_pool(keys=['pool2'])

        with pool1:
            with self.assertRaisesRegex(ValueError,
                                        'cannot enter multiple pools.*'
                                        f'{pool1.path}'):
                with pool2:
                    pass

    # TODO: This test may need to go at some point
    def test_asynchronous_pool_post_exit(self):
        """This test determines if all of the data is still in the cache when
        we are getting ready to exit. This was put here when ensuring we do not
        destroy our data when running asynchronous actions, and it can probably
        be removed once Archiver is reworked
        """
        plugin = get_dummy_plugin()
        concatenate_ints = plugin.methods['concatenate_ints']

        # This test needs to use a cache that exists past the lifespan of the
        # function
        cache = get_cache()
        test_pool = cache.create_pool(keys=[TEST_POOL], reuse=True)

        with test_pool:
            future = concatenate_ints.asynchronous(self.art1, self.art2,
                                                   self.art4, 4, 5)
            result = future.result()

        result = result[0]

        expected = set((
            './VERSION', f'data/{result._archiver.uuid}', 'keys/__TEST_POOL__',
            f'pools/__TEST_POOL__/{result._archiver.uuid}'
        ))

        atexit.unregister(_exit_cleanup)
        atexit.register(_on_exit_validate, cache, expected)
        atexit.register(_exit_cleanup)
