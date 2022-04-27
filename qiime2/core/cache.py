# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
import os
import yaml
import psutil
import shutil
import pathlib

import qiime2
from qiime2.sdk.result import Result
from qiime2.sdk.cache_config import CACHE_CONFIG
from qiime2.core.archive.archiver import Archiver

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

    def __init__(self, path, named=True):
        if path is not None:
            self.path = pathlib.Path(path)
        else:
            self.path = pathlib.Path(self.get_temp_path())

        self.named = named

        # Do we want a more rigorous check for whether or not we've been
        # pointed at an existing cache?
        if not os.path.exists(self.path):
            self.create_cache()
        elif not self.is_cache(self.path):
            raise ValueError(f"Path: \'{path}\' already exists and is not a"
                             " cache")

        # Make our process pool, we probably want this to happen somewhere else
        # CACHE_CONFIG.process_pool = self.create_pool(process_pool=True,
        #                                              reuse=True)

    def __enter__(self):
        """If you with a cache you are using it as a named cache
        """
        self.backup = CACHE_CONFIG.cache
        CACHE_CONFIG.cache = self

    def __exit__(self, *args):
        CACHE_CONFIG.cache = self.backup

    def get_process_pool(self):
        """ Before we save things into a process pool, we want to make sure
        we actually have a pool for the process and use it
        """
        pass

    # Surely this needs to be a thing? I suppose if they hand us a path that
    # doesn't exist we just create a cache there? do we want to create it at
    # that exact path do we slap a 'cache' sub-directory at that location and
    # use that?
    def create_cache(self):
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

    def get_temp_path(self):
        """ Get path to temp cache if the user did not specify a named cache.
        """
        # Get location of tmp
        TMPDIR = os.getenv("TMPDIR")
        # Get current username
        USER = os.getenv("USER")

        # If not set default to root tmp
        if TMPDIR is None:
            TMPDIR = '/tmp'

        # TODO: Gonna need to figure out setting our sticky bit, doing so seems
        # to make things explode
        user_path = os.path.join(TMPDIR, 'qiime2', USER)

        # return our path
        return user_path

    # Maybe this is create named pool specifically and we just auto create a
    # process pool
    def create_pool(self, keys=[], reuse=False, process_pool=False):
        # if reuse, look for an existing pool that matches keys
        # otherwise create a new pool matching keys and overwrite if one exists
        # Always create an anonymous pool keyed on pid-created_at@host in the
        # process folder
        # Need some kinda default name

        # If we are making a process pool just do it
        if process_pool:
            pool = Pool(self.process, self, reuse=reuse)
        else:
            pool_name = '_'.join(keys)
            pool_fp = self.pools / pool_name
            pool = Pool(pool_fp, self, name=pool_name, reuse=reuse)

            self.create_pool_keys(pool_name, keys)

        return pool

    def create_pool_keys(self, pool_name, keys):
        for key in keys:
            self._register_key(key, pool_name, pool=True)

    # Tell us if the path is a cache or not
    # NOTE: maybe we want this to be raising errors and whatnot instead of just
    # returning false?
    @classmethod
    def is_cache(cls, path):
        base_cache_contents = set(('data', 'keys', 'pools', 'process',
                                   'VERSION'))

        contents = set(os.listdir(path))
        if not contents.issuperset(base_cache_contents):
            return False

        regex = \
            re.compile(
                r"QIIME 2\ncache: \d+\nframework: 20\d\d\.")
        with open(path / 'VERSION') as fh:
            version_file = fh.read()
            return regex.match(version_file) is not None

    # Run the garbage collection algorithm
    # TODO: Needs to account for process pools
    def garbage_collection(self):
        referenced_pools = set()
        referenced_data = set()

        # Walk over keys and track all pools and data referenced
        for key in os.listdir(self.keys):
            loaded_key = yaml.safe_load(open(self.keys / key))
            referenced_pools.add(loaded_key['pool'])
            referenced_data.add(loaded_key['data'])

        # Since each key has at most a pool or data, we will end up with a None
        # in at least one of these sets. We don't want it
        referenced_pools.discard(None)
        referenced_data.discard(None)

        # Walk over pools and remove any that were not refered to by keys while
        # tracking all data within those that were referenced
        for pool in os.listdir(self.pools):
            if pool not in referenced_pools:
                shutil.rmtree(self.pools / pool)
            else:
                for data in os.listdir(self.pools / pool):
                    referenced_data.add(data)

        # Add references to data in process pools
        # TODO: Remove all process pools that are older than some configured
        # amount of time
        for process_pool in os.listdir(self.process):
            for data in os.listdir(self.process / process_pool):
                referenced_data.add(data)

        # Walk over all data and remove any that was not referenced
        for data in os.listdir(self.data):
            if data not in referenced_data:
                shutil.rmtree(self.data / data)

    # Save data and create key or pool entry
    def save(self, ref, key):
        # Move the data into cache under key
        shutil.copytree(ref._archiver.path, self.data, dirs_exist_ok=True)
        self._register_key(key, str(ref.uuid))

        # Collect garbage after a save
        self.garbage_collection()
        # Give back an instance of the Artifact they can use if they want
        return self.load(key)

    # Load the data pointed to by the key. Does not work on pools. Only works
    # if you have data
    def load(self, key):
        archiver = \
            Archiver.load(
                self.data / yaml.safe_load(open(self.keys / key))['data'],
                allow_no_op=True)
        return Result._from_archiver(archiver)

    # Remove key from cache
    def delete(self, key):
        os.remove(self.keys / key)
        self.garbage_collection()

    # Create a new key pointing at data or a pool
    def _register_key(self, key, value, pool=False):
        key_fp = self.keys / key

        if pool:
            key_fp.write_text(_KEY_TEMPLATE % (key, '', value))
        else:
            key_fp.write_text(_KEY_TEMPLATE % (key, value, ''))

    # Not entirely clear how this will work yet. We are assuming multiple
    # processes from multiple systems will be interacting with the cache. This
    # means we can't even safely assume unique PIDs. We will probably create
    # some kind of lock file to lock the entire cache or to list locked
    # elements of the cache or something
    def lock(self):
        pass

    def unlock(self):
        pass

    # Let's think about this a bit more than not at all. What do we do to check
    # for a lock, and if there is one, then what do we do? We need to in some
    # way wait for it to unlock, which suggests we need some way of managing
    # locks across multiple processes. We're working on a locking transactional
    # server in 565, and Dr. Otte told us to write a lock manager that exists
    # entirely to manage locks and nothing else. Surely we have to wait in a
    # queue until unlock or something right? We can't just drop what we're
    # doing because we have a lock. Could we conceivably use built in Python
    # locking?
    def check_lock(self):
        pass

    @property
    def data(self):
        return self.path / 'data'

    @property
    def keys(self):
        return self.path / 'keys'

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
    def __init__(self, path, cache, name=None, reuse=False):
        self.cache = cache

        if name:
            self.path = path
            self.name = name
        else:
            pid = os.getpid()
            user = os.getlogin()

            process = psutil.Process(pid)
            time = process.create_time()
            self.name = f'{pid}-{time}@{user}'
            self.path = path / self.name

        if not reuse and os.path.exists(self.path):
            raise ValueError("Pool already exists, please use reuse=True to "
                             "reuse existing pool, or remove all keys "
                             "indicating this pool to remove the pool")

        if not os.path.exists(self.path):
            os.mkdir(self.path)

    def __enter__(self):
        """If you with a pool you are using it as your named pool
        """
        self.old_pool = CACHE_CONFIG.named_pool
        CACHE_CONFIG.named_pool = self

    def __exit__(self, type, value, tb):
        CACHE_CONFIG.named_pool = self.old_pool

    def save(self, ref):
        shutil.copytree(ref._archiver.path, self.cache.data,
                        dirs_exist_ok=True)
        os.symlink(self.cache.data / str(ref.uuid), self.path / str(ref.uuid))
        self.cache.garbage_collection()

        return self.load(ref)

    # Load a reference to an element in the pool
    def load(self, ref):
        archiver = Archiver.load(self.cache.data / str(ref.uuid),
                                 allow_no_op=True)
        return Result._from_archiver(archiver)

    # Remove an element from the pool
    def remove(self, ref):
        os.remove(self.path / str(ref.uuid))
