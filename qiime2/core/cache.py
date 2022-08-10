# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
import os
import stat
import yaml
import time
import atexit
import psutil
import shutil
import getpass
import pathlib
import tempfile
import threading
from datetime import timedelta

from flufl.lock import Lock

import qiime2
from qiime2.sdk.result import Result
from qiime2.core.archive.archiver import Archiver
from qiime2.core.util import is_uuid4

_VERSION_TEMPLATE = """\
QIIME 2
cache: %s
framework: %s
"""

_KEY_TEMPLATE = """\
origin:
 %s
data:
 %s
pool:
 %s
"""

# Thread local indicating the cache to use
__CACHE__ = threading.local()
__CACHE__.cache = None

# TODO: Do we want this on the threadlocal? I feel like maybe we do
# Keep track of every cache used by this process for cleanup later
USED_CACHES = set()


def get_cache():
    """ Gets our cache if we have one and creates one in temp if we don't
    """
    # If we are on a new thread we may in fact not have a cache attribute here
    # at all
    if not hasattr(__CACHE__, 'cache') or __CACHE__.cache is None:
        __CACHE__.cache = Cache(None)

    return __CACHE__.cache


def _get_process_pool_name():
    """Gets the necessary info to create/identify a process pool for this
    process
    """
    pid = os.getpid()
    user = getpass.getuser()

    process = psutil.Process(pid)
    time = process.create_time()

    return f'{pid}-{time}@{user}'


@atexit.register
def _exit_cleanup():
    """For each cache used by this process we remove the process pool created
    by this process then run garbage collection
    """
    for cache in USED_CACHES:
        process_pool = _get_process_pool_name()
        target = cache.process / process_pool

        # There are several legitimate reasons the path could not exist. It
        # happens during our cache tests when the entire cache is nuked in the
        # end. It also happens in asycnhronous runs where the worker process
        # does not create a process pool (on mac this atexit is invoked on
        # workers). It could also happen if someone deleted the process pool
        # but... They probably shouldn't do that
        if os.path.exists(target):
            shutil.rmtree(target)
            cache.garbage_collection()


class Cache:
    """General structure of the cache (tmp optional)
    artifact_cache/
    ├── data
    │   ├── uuid1
    │   ├── uuid2
    │   ├── uuid3
    │   └── uuid4
    ├── keys
    │   ├── bar.yaml
    │   ├── baz.yaml
    │   └── foo.yaml
    ├── pools
    │   └── puuid1
    │       ├── uuid2 -> ../../data/uuid2/
    │       └── uuid3 -> ../../data/uuid3/
    ├── process
    ├── tmp
    └── VERSION
    Process folder contains pid-created_at@host some kinda reference to
    anonymous pools
    Create anonymous pools backing all final outputs. We have a named pool that
    tracks all intermediate and final results. We have an anonymous pool that
    tracks only final results (pid pool). Ensures that we don't gc the final
    results if we gc the named pool because they didn't pre register all keys
    """
    CURRENT_FORMAT_VERSION = '1'

    # The files and folders you expect to see at the top level of a cache
    base_cache_contents = set(('data', 'keys', 'pools', 'process', 'VERSION'))

    def __init__(self, path, process_pool_lifespan=1):
        """Creates a cache object backed by the directory specified by path. If
        no path is provided it gets a path to a temp cache.
        """
        if path is not None:
            self.path = pathlib.Path(path)
        else:
            self.path = pathlib.Path(self._get_temp_path())

        # Do we want a more rigorous check for whether or not we've been
        # pointed at an existing cache?
        if not os.path.exists(self.path):
            self._create_cache()
        elif not self.is_cache(self.path):
            raise ValueError(f"Path: \'{path}\' already exists and is not a"
                             " cache")

        self.lock = Lock(str(self.lockfile), lifetime=timedelta(minutes=10))
        # Make our process pool.
        # TODO: We currently will only end up with a process pool for the
        # process that originally launched QIIME 2 (not seperate ones for
        # parsl workers). We might want to change this in the future
        self.process_pool = self._create_process_pool()
        # Lifespan is supplied in days and converted to seconds for internal
        # use
        self.process_pool_lifespan = process_pool_lifespan * 3600 * 24
        # This is set if a named pool is created on this cache and withed in
        self.named_pool = None

        # We were used by this process
        USED_CACHES.add(self)

    def __enter__(self):
        """Set this cache on the thread local
        """
        self.backup = __CACHE__.cache
        __CACHE__.cache = self

    def __exit__(self, *args):
        """Set the thread local back to whatever cache it was using before this
        one
        """
        __CACHE__.cache = self.backup

    # TODO: This is only going to need to be a thing if we want worker
    # processes to have their own process pools. I'm leaving it here as a
    # skeleton for now, but if we decide we're cool with just having a process
    # pool for the process that started QIIME 2 then we probably don't need
    # anything like this
    def get_process_pool(self):
        """Before we save things into a process pool, we want to make sure we
        actually have a pool for the process and use it
        """
        pass

    def _create_cache(self):
        """Create the cache directory, all sub directories, and the version
        file
        """
        # Construct the cache root recursively
        os.makedirs(self.path)
        os.mkdir(self.path / 'data')
        os.mkdir(self.path / 'keys')
        os.mkdir(self.path / 'pools')
        os.mkdir(self.path / 'process')
        # Do we want this right off the bat? How exactly is setting tmp in the
        # cache going to work? tmp is never going to be managed by the cache,
        # it's just so they're both on the same disk, so they'll probably just
        # set the tmp location in the config or something. I feel like if we're
        # going to manage the cache, we should manage the cache which means if
        # they're going to create_poolput tmp in the cache it should have to be
        # in a set directory within the cache like tmp not just whatever they
        # want it to be in the cache. Not sure how we would really enforce
        # that, but we can just... Heavily encourage it I guess
        # os.mkdir('tmp')

        self.version.write_text(
            _VERSION_TEMPLATE % (self.CURRENT_FORMAT_VERSION,
                                 qiime2.__version__))

    def _get_temp_path(self):
        """ Get path to temp cache if the user did not specify a named cache.
        """
        TMPDIR = tempfile.gettempdir()
        USER = getpass.getuser()

        cache_dir = os.path.join(TMPDIR, 'qiime2')
        if not os.path.exists(cache_dir):
            os.mkdir(cache_dir)

        # Make sure the sticky bit is set on the cache directory. Documentation
        # on what a sitcky bit is found here
        # https://docs.python.org/3/library/stat.html#stat.S_ISVTX
        permissions = os.stat(cache_dir).st_mode
        sticky_permissions = permissions | stat.S_ISVTX
        os.chmod(cache_dir, sticky_permissions)

        user_path = os.path.join(TMPDIR, 'qiime2', USER)
        return user_path

    def _create_process_pool(self):
        """Creates a process pool which is identical in function to a named
        pool, but it lives in the process subdirectory not the pools
        subdirectory, and is handled differently by garbage collection due to
        being unkeyed
        """
        return Pool(self, reuse=True)

    def create_pool(self, keys=[], reuse=False):
        """ Used internally to create the process pool and externally to create
        named pools.
        keys: A list of keys to point to a named pool. The pool name will be a
        concatenation of all the keys
        reuse: If False, we will error if the pool we are trying to create
        already exists (if the path to it already exists)
        """
        pool_name = '_'.join(keys)
        pool = Pool(self, name=pool_name, reuse=reuse)

        self._create_pool_keys(pool_name, keys)

        return pool

    def _create_pool_keys(self, pool_name, keys):
        """A pool can have many keys refering to it.
        """
        for key in keys:
            self._register_key(key, pool_name, pool=True)

    @classmethod
    def is_cache(cls, path):
        """Tells us if the path we were given is a cache
        """
        contents = set(os.listdir(path))
        if not contents.issuperset(cls.base_cache_contents):
            return False

        regex = \
            re.compile(
                r"QIIME 2\ncache: \d+\nframework: 20\d\d\.")
        with open(path / 'VERSION') as fh:
            version_file = fh.read()
            return regex.match(version_file) is not None

    def garbage_collection(self):
        """Runs garbage collection on the cache. We log all data and pools
        pointed to by keys. Then we go through all pools and delete any that
        were not referred to by a key while logging all data in pools that are
        refered to by keys. Then we go through all process pools and log all
        data they point to. Then we go through the data and remove any that was
        not logged. This only destroys data and named pools.
        """
        referenced_pools = set()
        referenced_data = set()

        # Walk over keys and track all pools and data referenced
        # This needs to be locked so we ensure that we don't have other threads
        # or processes writing refs that we don't see leading to us deleting
        # their data
        with self.lock:
            for key in os.listdir(self.keys):
                loaded_key = yaml.safe_load(open(self.keys / key))
                referenced_pools.add(loaded_key['pool'])
                referenced_data.add(loaded_key['data'])

            # Since each key has at most a pool or data, we will end up with a
            # None in at least one of these sets. We don't want it
            referenced_pools.discard(None)
            referenced_data.discard(None)

            # Walk over pools and remove any that were not refered to by keys
            # while tracking all data within those that were referenced
            for pool in os.listdir(self.pools):
                if pool not in referenced_pools:
                    shutil.rmtree(self.pools / pool)
                else:
                    for data in os.listdir(self.pools / pool):
                        referenced_data.add(data)

            # Add references to data in process pools
            for process_pool in os.listdir(self.process):
                # Pick the creation time out of the pool name of format
                # {pid}-time@user
                create_time = float(process_pool.split('-')[1].split('@')[0])

                if time.time() - create_time >= self.process_pool_lifespan:
                    shutil.rmtree(self.process / process_pool)
                else:
                    for data in os.listdir(self.process / process_pool):
                        referenced_data.add(data)

            # Walk over all data and remove any that was not referenced
            for data in os.listdir(self.data):
                # If this assert is ever tripped something real bad happened
                assert is_uuid4(data)

                if data not in referenced_data:
                    shutil.rmtree(self.data / data, ignore_errors=True)

    def save(self, ref, key):
        """Create our key then create our data. Returns a version of the data
        backed by the key in the cache
        """
        # Create the key before the data, this is so that if another thread or
        # process is running garbage collection it doesn't see our unkeyed data
        # and remove it leaving us with a dangling reference and no data
        with self.lock:
            self._register_key(key, str(ref.uuid))

        if not os.path.exists(self.data / str(ref.uuid)):
            shutil.copytree(ref._archiver.path, self.data, dirs_exist_ok=True)

        # Give back an instance of the Artifact they can use if they want
        return self.load(key)

    def _register_key(self, key, value, pool=False):
        """Create a new key pointing at data or a named pool
        """
        if not key.isidentifier():
            raise ValueError('Key must be a valid Python identifier. Python '
                             'identifier rules may be found here '
                             'https://www.askpython.com/python/'
                             'python-identifiers-rules-best-practices')

        key_fp = self.keys / key

        if pool:
            key_fp.write_text(_KEY_TEMPLATE % (key, '', value))
        else:
            key_fp.write_text(_KEY_TEMPLATE % (key, value, ''))

    def load(self, key):
        """Load the data pointed to by a key. Only works on a key that refers
        to a data item will error on a key that points to a pool
        """
        archiver = \
            Archiver.load(
                self.data / yaml.safe_load(open(self.keys / key))['data'],
                allow_no_op=True)
        return Result._from_archiver(archiver)

    def remove(self, key):
        """Remove a key from the cache then run garbage collection to remove
        anything it was referencing and any other loose data
        """
        os.remove(self.keys / key)
        self.garbage_collection()

    def clear_locks(self):
        """Clears all locks in this cache and in pools under this cache
        NOTE: Forcibly removes the locks outside of the locking library's API
        """
        if os.path.exists(self.lockfile):
            os.remove(self.lockfile)

    @property
    def data(self):
        return self.path / 'data'

    @property
    def keys(self):
        return self.path / 'keys'

    @property
    def lockfile(self):
        return self.path / 'LOCK'

    @property
    def pools(self):
        return self.path / 'pools'

    @property
    def process(self):
        return self.path / 'process'

    @property
    def version(self):
        return self.path / 'VERSION'


# Assume we will make this its own class for now
class Pool:
    def __init__(self, cache, name=None, reuse=False):
        """ Used with name=None and reuse=True to create a process pool. Used
        with a name to create named pools
        """
        # The pool keeps track of the cache it belongs to
        self.cache = cache

        # If they are creating a named pool, we already have this info
        if name:
            self.name = name
            self.path = cache.pools / name
        # The alternative is that we have a process pool. We want this pool to
        # exist in the process directory under the cache not the pools
        # directory. The name follows the scheme
        # pid-process_start_time@user
        else:
            self.name = _get_process_pool_name()
            self.path = cache.process / self.name

        # Raise a value error if we thought we were making a new pool but
        # actually are not
        if not reuse and os.path.exists(self.path):
            raise ValueError("Pool already exists, please use reuse=True to "
                             "reuse existing pool, or remove all keys "
                             "indicating this pool to remove the pool")

        if not os.path.exists(self.path):
            os.mkdir(self.path)

    def __enter__(self):
        """Set this pool to be our named pool on the current cache
        """
        self.old_pool = self.cache.named_pool
        self.cache.named_pool = self

    def __exit__(self, *args):
        """Set the named pool on the current cache back to whatever it was
        before this one
        """
        self.cache.named_pool = self.old_pool

    def save(self, ref):
        """Save the data into the pool then load a new ref backed by the data
        in the pool
        """
        # Create the key before the data, this is so that if another thread or
        # process is running garbage collection it doesn't see our unkeyed data
        # and remove it leaving us with a dangling reference and no data
        with self.cache.lock:
            os.symlink(self.cache.data / str(ref.uuid),
                       self.path / str(ref.uuid))

        if not (self.cache.data / str(ref.uuid)).exists():
            shutil.copytree(ref._archiver.path, self.cache.data,
                            dirs_exist_ok=True)

        return self.load(ref)

    def load(self, ref):
        """Load a reference to an element in the pool
        """
        archiver = Archiver.load(self.cache.data / str(ref.uuid),
                                 allow_no_op=True)
        return Result._from_archiver(archiver)

    def remove(self, ref):
        """Remove an element from the pool
        """
        # TODO: This guard should be removed when we rework the logic
        if (self.path / str(ref.uuid)).exists():
            os.remove(self.path / str(ref.uuid))
            self.cache.garbage_collection()
