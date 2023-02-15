# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
"""
The cache is used to store unzipped data on disk in a predictable and user
controlled location. This allows us to skip constantly zipping and unzipping
large amounts of data and taking up CPU time when storage space is not an
issue. It also allows us to to know exactly what data has been created and
where.

By default, a cache will be created under $TMPDIR/qiime2/$USER and all
intermediate data created by QIIME 2 as it executes will be written into that
directory. This means QIIME 2 reserves usage of the $TMPDIR/qiime2 directory.
The user may also specify a new location to be used in place of this default
directory. This location must meet a few criteria.

**1.** It must be writable from any and all locations the QIIME 2 command
intending to use it will be running. This means that in an HPC context, the
location specified for the cache must be writable from the node QIIME 2 will be
executing on.

**2.** It must either not exist or already be a cache. The first time a
directory is specified to be used as a cache, it should not exist. QIIME 2 will
create a cache structure on disk at that location. Any existing directory you
attempt to use as a cache should have been created as a cache by QIIME 2.

"""
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
import weakref
import tempfile
import warnings
import threading
from sys import maxsize
from random import randint
from datetime import timedelta

import flufl.lock

import qiime2
from .path import ArchivePath
from qiime2.sdk.result import Result
from qiime2.core.util import (is_uuid4, set_permissions, touch_under_path,
                              READ_ONLY_FILE, READ_ONLY_DIR, USER_GROUP_RWX)
from qiime2.core.archive.archiver import Archiver

_VERSION_TEMPLATE = """\
QIIME 2
cache: %s
framework: %s
"""

# Thread local indicating the cache to use
_CACHE = threading.local()
_CACHE.cache = None
_CACHE.temp_cache = None

# TODO: Do we want this on the thread local? I feel like maybe we do
# Keep track of every cache used by this process for cleanup later
USED_CACHES = set()

# These permissions are directory with sticky bit and rwx for all set
EXPECTED_PERMISSIONS = 0o41777


def get_cache():
    """Gets the cache we have instructed QIIME 2 to use in this invocation.
    By default this is a cache located at $TMPDIR/qiime2/$USER, but if the
    user has set a cache it is the cache they set. This is used by various
    parts of the framework to determine what cache they should be saving
    to/loading from.

    Returns
    -------
    Cache
        The cache QIIME 2 is using for the current invocation.

    Examples
    --------
    >>> test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
    >>> cache_path = os.path.join(test_dir.name, 'cache')
    >>> cache = Cache(cache_path)
    >>> # get_cache() will return the temp cache, not the one we just made.
    >>> get_cache() == cache
    False
    >>> # After withing in the cache we just made, get_cache() will return it.
    >>> with cache:
    ...     get_cache() == cache
    True
    >>> # Now that we have exited our cache, we will get the temp cache again.
    >>> get_cache() == cache
    False
    >>> test_dir.cleanup()

    """
    # If we are on a new thread we may in fact not have a cache attribute here
    # at all
    if not hasattr(_CACHE, 'cache') or _CACHE.cache is None:
        if not hasattr(_CACHE, 'temp_cache') or _CACHE.temp_cache is None:
            _CACHE.temp_cache = Cache()
        return _CACHE.temp_cache

    return _CACHE.cache


def _get_temp_path():
    """Get path to temp cache if the user did not specify a named cache.
    This function will create the path if it does not exist and ensure it
    is suitable for use as a cache if it does.

    Returns
    -------
    str
        The path created for the temp cache.
    """
    tmpdir = tempfile.gettempdir()

    cache_dir = os.path.join(tmpdir, 'qiime2')

    # Make sure the sticky bit is set on the cache directory. Documentation on
    # what a sticky bit is can be found here
    # https://docs.python.org/3/library/stat.html#stat.S_ISVTX We also set
    # read/write/execute permissions for everyone on this directory. We only do
    # this if we are the owner of the /tmp/qiime2  directory or in other words
    # the first person to run QIIME 2 with this /tmp since the /tmp was wiped
    if not os.path.exists(cache_dir):
        try:
            os.mkdir(cache_dir)
        except FileExistsError:
            # we know that it didn't exist a moment ago, so we're probably
            # about to set it up in a different process
            time.sleep(0.5)
            # this sleep is to give the first process enough time to create
            # a cache object which we will then re-use. Ideally this would
            # be handled with a lock, but we don't have anywhere to put it
            # yet. Since this is the kind of thing that can only happen when
            # QIIME 2 has to create a new temp cache and there's a race for it
            # this small hack seems not too bad.
        else:
            # skip this if there was an error we ignored
            sticky_permissions = stat.S_ISVTX | stat.S_IRWXU | stat.S_IRWXG \
                | stat.S_IRWXO
            os.chmod(cache_dir, sticky_permissions)
    elif os.stat(cache_dir).st_mode != EXPECTED_PERMISSIONS:
        raise ValueError(f"Directory '{cache_dir}' already exists without "
                         f"proper permissions '{oct(EXPECTED_PERMISSIONS)}' "
                         "set. Current permissions are "
                         f"'{oct(os.stat(cache_dir).st_mode)}.' This most "
                         "likely means something other than QIIME 2 created "
                         f"the directory '{cache_dir}' or QIIME 2 failed "
                         f"between creating '{cache_dir}' and setting "
                         "permissions on it.")

    user = _get_user()
    user_dir = os.path.join(cache_dir, user)

    # It is conceivable that we already have a path matching this username that
    # belongs to another uid, if we do then we want to create a garbage name
    # for the temp cache that will be used by this user
    if os.path.exists(user_dir) and os.stat(user_dir).st_uid != os.getuid():
        uid_name = _get_uid_cache_name()
        # This really shouldn't happen
        if user == uid_name:
            raise ValueError(f'Temp cache for uid path {user} already exists '
                             'but does not belong to us.')

        user_dir = os.path.join(cache_dir, uid_name)

    return user_dir


def _get_user():
    """Get the uname for our default cache. Internally getpass.getuser is
    getting the uid then looking up the username associated with it. This could
    fail it we are running inside a container because the container is looking
    for its parent's uid in its own /etc/passwd which is unlikely to contain a
    user associated with that uid. If that failure does occur, we create an
    alternate default username.

    Returns
    -------
    str
        The value we will be using as uname for our default cache.
    """
    try:
        return getpass.getuser()
    except KeyError:
        return _get_uid_cache_name()


def _get_uid_cache_name():
    """Create an esoteric name that is unlikely to be the name of a real user
    in cases were getpass.getuser fails. This name is of the form
    uid=#<uid> which should be consistent across invocations of this
    function by the same user.

    Returns
    -------
    str
        The aforementioned stand in name.
    """
    return f'uid=#{os.getuid()}'


@atexit.register
def _exit_cleanup():
    """When the process ends, for each cache used by this process we remove the
    process pool created by this process then run garbage collection.
    """
    for cache in USED_CACHES:
        target = cache.processes / os.path.basename(cache.process_pool.path)

        # There are several legitimate reasons the path could not exist. It
        # happens during our cache tests when the entire cache is nuked in the
        # end. It also happens in asynchronous runs where the worker process
        # does not create a process pool (on Mac this atexit is invoked on
        # workers). It could also happen if someone deleted the process pool
        # but... They probably shouldn't do that
        if os.path.exists(target):
            shutil.rmtree(target)
            cache.garbage_collection()


def monitor_thread(cache_dir, is_done):
    """MacOS reaps temp files that are three days old or older. This function
    will be running in a separate daemon and making sure MacOS doesn't cull
    anything still needed by a long running process by periodically updating
    the last accessed times on all files in the cache by touching them every
    six hours. The daemon running this function will be terminated when the
    process that invoked it ends.

    Parameters
    ----------
    cache_dir : str or PathLike object
        The path to the cache that invoked this daemon.
    is_done : threading.Event
        The process that invoked this daemon sets this flag on exit to notify
        this daemon to terminate.
    """
    while not is_done.is_set():
        touch_under_path(cache_dir)
        time.sleep(60 * 60 * 6)


# This is very important to our trademark
tm = object


class MEGALock(tm):
    """ We need to lock out other processes with flufl, but we also need to
    lock out other threads with a Python thread lock (because parsl
    threadpools), so we put them together in one MEGALock(tm)
    """

    def __init__(self, flufl_fp, lifetime):
        self.flufl_fp = flufl_fp
        self.lifetime = lifetime
        self.re_entries = 0

        self.thread_lock = threading.Lock()
        self.flufl_lock = flufl.lock.Lock(flufl_fp, lifetime=lifetime)

    def __enter__(self):
        """ We acquire the thread lock first because the flufl lock isn't
        thread-safe which is why we need both locks in the first place
        """
        if self.re_entries == 0:
            self.thread_lock.acquire()
            self.flufl_lock.lock()

        self.re_entries += 1

    def __exit__(self, *args):
        if self.re_entries > 0:
            self.re_entries -= 1

        if self.re_entries == 0:
            self.flufl_lock.unlock()
            self.thread_lock.release()

    def __getstate__(self):
        lockless_dict = self.__dict__.copy()

        del lockless_dict['thread_lock']
        del lockless_dict['flufl_lock']

        return lockless_dict

    def __setstate__(self, state):
        self.__dict__.update(state)

        self.thread_lock = threading.Lock()
        self.flufl_lock = \
            flufl.lock.Lock(self.flufl_fp, lifetime=self.lifetime)


class Cache:
    """General structure of the cache:

    ::

        artifact_cache/
        ├── data/
        │   ├── uuid1/
        │   ├── uuid2/
        │   ├── uuid3/
        │   └── uuid4/
        ├── keys/
        │   ├── bar.yaml
        │   ├── baz.yaml
        │   └── foo.yaml
        ├── pools/
        │   └── puuid1/
        │       ├── uuid1 -> ../../data/uuid1/
        │       └── uuid2 -> ../../data/uuid2/
        ├── processes/
        │   └── <process-id>-<process-create-time>@<user>/
        │       ├── uuid3 -> ../../data/uuid3/
        │       └── uuid4 -> ../../data/uuid4/
        └── VERSION

    **Data:** The data directory contains all of the artifacts in the cache in
    unzipped form.

    **Keys:** The keys directory contains yaml files that refer to either a
    piece of data or a pool. The data/pool referenced by the key will be kept
    as long as the key exists.

    **Pools:** The pools directory contains all named (keyed) pools in the
    cache. Each pool contains symlinks to all of the data it contains.

    **Processes:** The processes directory contains process pools of the format
    <process-id>-<process-create-time>@<user> for each process that has used
    this cache. Each pool contains symlinks to each element in the data
    directory the process that created the pool has used in some way (created,
    loaded, etc.). These symlinks are ephemeral and have lifetimes <= the
    lifetime of the process that created them. More permanent storage is done
    using keys.

    **VERSION:** This file contains some information QIIME 2 uses to determine
    what version of QIIME 2 was used to create the cache and what version of
    cache it is (if we make breaking changes in the future this version number
    will allow for backwards compatibility).
    """
    CURRENT_FORMAT_VERSION = '1'

    # The files and folders you expect to see at the top level of a cache
    base_cache_contents = \
        set(('data', 'keys', 'pools', 'processes', 'VERSION'))

    def __new__(cls, path=None):
        if path is None:
            path = _get_temp_path()

        # We have to ensure we really have the same path here because otherwise
        # something as simple as path='/tmp/qiime2/x' and path='/tmp/qiime2/x/'
        # would create two different Cache objects
        for cache in USED_CACHES:
            if os.path.exists(path) and os.path.exists(cache.path) and \
                    os.path.samefile(path, cache.path):
                return cache

        return super(Cache, cls).__new__(cls)

    def __init__(self, path=None, process_pool_lifespan=45):
        """Creates a Cache object backed by the directory specified by path. If
        no path is provided, it gets a path to a temp cache.

        Warning
        -------
        If no path is provided and the path $TMPDIR/qiime2/$USER exists but is
        not a valid cache, we remove the directory and create a cache there.

        Parameters
        ----------
        path : str or PathLike object
            Should point either to a non-existent writable directory to be
            created as a cache or to an existing writable cache. Defaults to
            None which creates the cache at $TMPDIR/qiime2/$USER.
        process_pool_lifespan : int
            The number of days we should allow process pools to exist for
            before culling them.
        """
        # If this is a new cache or if the cache somehow got invalidated
        # (MacOS culling) then we need to re-init the cache. This could
        # theoretically cause us to end up with two Cache instances pointing at
        # the same path again should a cache be in some way invalidated during
        # the lifetime of a process with an existing Cache instance pointing to
        # it, but if that happens you're probably in trouble anyway.
        if self not in USED_CACHES or not self.is_cache(self.path):
            self.__init(path=path, process_pool_lifespan=process_pool_lifespan)

    def __init(self, path=None, process_pool_lifespan=45):
        if path is not None:
            self.path = pathlib.Path(path)
        else:
            self.path = pathlib.Path(_get_temp_path())

        if not os.path.exists(self.path):
            os.makedirs(self.path)

        self.lock = \
            MEGALock(str(self.lockfile), lifetime=timedelta(minutes=10))

        # We need to lock here to ensure that if we have multiple processes
        # trying to create the same cache one of them can actually succeed at
        # creating the cache without interference from the other processes.
        with self.lock:
            if not Cache.is_cache(self.path):
                try:
                    self._create_cache_contents()
                except FileExistsError as e:
                    if path is None:
                        warnings.warn(
                            "Your temporary cache was found to be in an "
                            "inconsistent state. It has been recreated.")
                        set_permissions(self.path, USER_GROUP_RWX,
                                        USER_GROUP_RWX, skip_root=True)
                        self._remove_cache_contents()
                        self._create_cache_contents()
                    else:
                        raise ValueError(
                            f"Path: \'{self.path}\' already exists and is not "
                            "a cache.") from e

        # Make our process pool.
        self.process_pool = self._create_process_pool()
        # Lifespan is supplied in days and converted to seconds for internal
        # use
        self.process_pool_lifespan = process_pool_lifespan * 3600 * 24
        # This is set if a named pool is created on this cache and withed in
        self.named_pool = None

        # We were used by this process
        USED_CACHES.add(self)

        # Start thread that pokes things in the cache to ensure they aren't
        # culled for being too old (only if we are in a temp cache)
        if path is None:
            self._thread_is_done = threading.Event()
            self._thread_destructor = \
                weakref.finalize(self, self._thread_is_done.set)

            self._thread = threading.Thread(
                target=monitor_thread, args=(self.path, self._thread_is_done),
                daemon=True)

            self._thread.start()

    def __enter__(self):
        """Tell QIIME 2 to use this cache in its current invocation (see
        get_cache).
        """
        if _CACHE.cache is not None and _CACHE.cache.path != self.path:
            raise ValueError("You cannot enter multiple caches at once, "
                             "currently entered cache is located at: "
                             f"'{_CACHE.cache.path}'")

        _CACHE.cache = self

    def __exit__(self, *args):
        """Tell QIIME 2 to go back to using the default cache.
        """
        _CACHE.cache = None

    def __getstate__(self):
        """Tell the cache not to pickle anything related to the daemon that
        keeps files around on MacOS because it can't pickle, and we don't need
        it after pickling and rehydrating. It will already be managed by the
        original process.
        """
        threadless_dict = self.__dict__.copy()

        # This will only even exist if we are a temp cache not a named cache.
        # If _thread exists the others should as well
        if '_thread' in threadless_dict:
            del threadless_dict['_thread_is_done']
            del threadless_dict['_thread_destructor']
            del threadless_dict['_thread']

        return threadless_dict

    @classmethod
    def is_cache(cls, path):
        """Tells us if the path we were given is a cache.

        Parameters
        ----------
        path : str or PathLike object
            The path to the cache we are checking.

        Returns
        -------
        bool
            Whether the path we were given is a cache or not.

        Examples
        --------
        >>> test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        >>> cache_path = os.path.join(test_dir.name, 'cache')
        >>> cache = Cache(cache_path)
        >>> Cache.is_cache(cache_path)
        True
        >>> test_dir.cleanup()
        """
        path = pathlib.Path(path)
        contents = set(os.listdir(path))
        if not contents.issuperset(cls.base_cache_contents):
            return False

        regex = \
            re.compile(
                r"QIIME 2\ncache: \d+\nframework: 20\d\d\.")
        with open(path / 'VERSION') as fh:
            version_file = fh.read()
            return regex.match(version_file) is not None

    def _create_cache_contents(self):
        """Create the cache directory, all sub directories, and the version
        file.
        """
        os.mkdir(self.data)
        os.mkdir(self.keys)
        os.mkdir(self.pools)
        os.mkdir(self.processes)

        self.version.write_text(
            _VERSION_TEMPLATE % (self.CURRENT_FORMAT_VERSION,
                                 qiime2.__version__))

    def _remove_cache_contents(self):
        """Removes everything in a cache that isn't a lock file. If you want to
        completely remove a cache, just use shutil.rmtree (make sure you have
        permissions).

        Note
        ----
        We ignore lock files because we want the process that is running this
        method to maintain its lock on the cache.
        """
        for elem in os.listdir(self.path):
            if 'LOCK' not in elem:
                fp = os.path.join(self.path, elem)
                if os.path.isdir(fp):
                    shutil.rmtree(os.path.join(self.path, fp))
                else:
                    os.unlink(fp)

    def _create_process_pool(self):
        """Creates a process pool which is identical in function to a named
        pool, but it lives in the processes subdirectory not the pools
        subdirectory, and is handled differently by garbage collection due to
        being un-keyed. Process pools are used to keep track of results for
        currently running processes and are removed when the process that
        created them ends.

        Returns
        -------
        Pool
            The pool we created.
        """
        return Pool(self, reuse=True)

    def _create_collection_pool(self, ref_collection, key):
        pool_name = f'{key}_collection'
        pool = Pool(self, name=pool_name, reuse=False)
        self._register_key(
            key, pool_name, pool=True, collection=ref_collection)

        return pool

    def create_pool(self, keys=[], reuse=False):
        """Used to create named pools. A named pool's name is all of the keys
        given for it separated by underscores. All of the given keys are
        created individually and refer to the named pool as opposed to saving a
        single piece of data where a single key is created referring to that
        data.

        Named pools can be used by pipelines to store all intermediate results
        created by the pipeline and prevent it from being reaped. This allows
        us to resume failed pipelines by collecting all of the data the
        pipeline saved to the named pool before it crashed and reusing it so we
        don't need to run the steps that created it again and can instead rerun
        the pipeline from where it failed.

        Once the pipeline completes, all of its final results will be saved to
        the pool as well with the idea being that the user can then reuse the
        pool keys to refer to the final data and get rid of the pool now that
        the pipeline that created it has completed.

        Parameters
        ----------
        keys : List[str]
            A list of keys to use to reference the pool.
        reuse : bool
            Whether to reuse a pool if a pool with the given keys already
            exists.

        Returns
        -------
        Pool
            The pool we created.

        Examples
        --------
        >>> test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        >>> cache_path = os.path.join(test_dir.name, 'cache')
        >>> cache = Cache(cache_path)
        >>> pool = cache.create_pool(keys=['some', 'kinda', 'keys'])
        >>> cache.get_keys() == set(['some', 'kinda', 'keys'])
        True
        >>> cache.get_pools() == set(['some_kinda_keys'])
        True
        >>> test_dir.cleanup()
        """
        pool_name = '_'.join(keys)
        pool = Pool(self, name=pool_name, reuse=reuse)

        self._create_pool_keys(pool_name, keys)

        return pool

    def _create_pool_keys(self, pool_name, keys):
        """A pool can have many keys referring to it. This function creates all
        of the keys referring to the pool.

        Parameters
        ----------
        pool_name : str
            The name of the pool we are keying.
        keys : List[str]
            A list of all the keys to create referring to the pool.
        """
        for key in keys:
            self._register_key(key, pool_name, pool=True)

    def garbage_collection(self):
        """Runs garbage collection on the cache in the following steps:

        **1.** Iterate over all keys and log all data and pools referenced by
        the keys.

        **2.** Iterate over all named pools and delete any that were not
        referred to by a key while logging all data in pools that were referred
        to by keys.

        **3.** Iterate over all process pools and log all data they refer to.

        **4.** Iterate over all data and remove any that was not referenced.

        This process destroys data and named pools that do not have keys along
        with process pools older than the process_pool_lifespan on the cache
        which defaults to 45 days. It never removes keys.

        We lock out other processes and threads from accessing the cache while
        garbage collecting to ensure the cache remains in a consistent state.
        """
        referenced_pools = set()
        referenced_data = set()

        # Walk over keys and track all pools and data referenced
        # This needs to be locked so we ensure that we don't have other threads
        # or processes writing refs that we don't see leading to us deleting
        # their data
        with self.lock:
            for key in self.get_keys():
                with open(self.keys / key) as fh:
                    loaded_key = yaml.safe_load(fh)

                if 'data' in loaded_key:
                    referenced_data.add(loaded_key['data'])
                elif 'pool' in loaded_key:
                    referenced_pools.add(loaded_key['pool'])
                # This really should never be happening unless someone messes
                # with things manually
                else:
                    raise ValueError(f"The key '{key}' in the cache"
                                     f"'{self.path}' does not point to"
                                     " anything")

            # Walk over pools and remove any that were not referred to by keys
            # while tracking all data within those that were referenced
            for pool in self.get_pools():
                if pool not in referenced_pools:
                    shutil.rmtree(self.pools / pool)
                else:
                    for data in os.listdir(self.pools / pool):
                        referenced_data.add(data)

            # Add references to data in process pools
            for process_pool in self.get_processes():
                # Pick the creation time out of the pool name of format
                # <process-id>-<process-create-time>@<user>
                create_time = float(process_pool.split('-')[1].split('@')[0])

                if time.time() - create_time >= self.process_pool_lifespan:
                    shutil.rmtree(self.processes / process_pool)
                else:
                    for data in os.listdir(self.processes / process_pool):
                        referenced_data.add(data.split('.')[0])

            # Walk over all data and remove any that was not referenced
            for data in self.get_data():
                # If this assert is ever tripped something real bad happened
                assert is_uuid4(data)

                if data not in referenced_data:
                    target = self.data / data

                    set_permissions(target, None, USER_GROUP_RWX)
                    shutil.rmtree(target)

    def save(self, ref, key):
        """Saves data into the cache by creating a key referring to the data
        then copying the data if it is not already in the cache.

        Parameters
        ----------
        ref : Result
            The QIIME 2 result we are saving into the cache.
        key : str
            The key we are saving the result under.

        Returns
        -------
        Result
            A Result backed by the data in the cache.

        Examples
        --------
        >>> from qiime2.sdk.result import Artifact
        >>> from qiime2.core.testing.type import IntSequence1
        >>> test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        >>> cache_path = os.path.join(test_dir.name, 'cache')
        >>> cache = Cache(cache_path)
        >>> artifact = Artifact.import_data(IntSequence1, [0, 1, 2])
        >>> saved_artifact = cache.save(artifact, 'key')
        >>> # save returned an artifact that is backed by the data in the cache
        >>> str(saved_artifact._archiver.path) == \
                str(cache.data / str(artifact.uuid))
        True
        >>> cache.get_keys() == set(['key'])
        True
        >>> test_dir.cleanup()
        """
        # Create the key before the data, this is so that if another thread or
        # process is running garbage collection it doesn't see our un-keyed
        # data and remove it leaving us with a dangling reference and no data
        with self.lock:
            self._register_key(key, str(ref.uuid))
            self._copy_to_data(ref)

        return self.load(key)

    def save_collection(self, ref_collection, key):
        """Saves a Collection to a pool in the cache with the given key. This
        pool's key file will keep track of the order of the Collection.
        """
        with self.lock:
            pool = self._create_collection_pool(ref_collection, key)
            # self._register_key(key, ref_collection, pool=True, ordered=True)
            for ref in list(ref_collection.values()):
                pool.save(ref)

        return self.load_collection(key)

    def _register_key(self, key, value, pool=False, collection=None):
        """Creates a key file pointing at the specified data or pool.

        Parameters
        ----------
        key : str
            The name of the key to create.
        value : str
            The path to the data or pool we are keying.
        pool : bool
            Whether we are keying a pool or not.

        Raises
        ------
        ValueError
            If the key passed in is not a valid Python identifier. We enforce
            this to ensure no one creates keys that cause issues when we try to
            load them.
        """
        if not key.isidentifier():
            raise ValueError('Key must be a valid Python identifier. Python '
                             'identifier rules may be found here '
                             'https://www.askpython.com/python/'
                             'python-identifiers-rules-best-practices')

        key_fp = self.keys / key

        key = {}
        key['origin'] = key

        if pool:
            key['pool'] = value

            if collection is not None:
                key['order'] = \
                    [{k: str(v.uuid)} for k, v in collection.items()]
        else:
            key['data'] = value

            if collection is not None:
                raise ValueError('An ordered Collection key can only be made'
                                 ' for a pool.')

        with open(key_fp, 'w') as fh:
            yaml.dump(key, fh)

    def read_key(self, key):
        """Reads the contents of a given key.

        Parameters
        ----------
        key : str
            The name of the key to read

        Returns
        -------
        dict
            Maps 'data' -> the data referenced or 'pool' -> the pool
            referenced. Only 'data' or 'pool' will have a value the other will
            be none.

        Raises
        ------
        KeyError
            If the key does not exists in the cache.
        """
        with self.lock:
            try:
                with open(self.keys / key) as fh:
                    return yaml.safe_load(fh)
            except FileNotFoundError as e:
                raise KeyError(f"The cache '{self.path}' does not contain the "
                               f"key '{key}'") from e

    def load(self, key):
        """Loads the data pointed to by a key. Only works on keys that refer to
        data items and will error on keys that refer to pools.

        Parameters
        ----------
        key : str
            The key to the data we are loading.

        Returns
        -------
        Result
            The loaded data pointed to by the key.

        Raises
        ------
        ValueError
            If the key does not reference any data meaning you probably tried
            to load a pool.

        Examples
        --------
        >>> from qiime2.sdk.result import Artifact
        >>> from qiime2.core.testing.type import IntSequence1
        >>> test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        >>> cache_path = os.path.join(test_dir.name, 'cache')
        >>> cache = Cache(cache_path)
        >>> artifact = Artifact.import_data(IntSequence1, [0, 1, 2])
        >>> saved_artifact = cache.save(artifact, 'key')
        >>> loaded_artifact = cache.load('key')
        >>> loaded_artifact == saved_artifact == artifact
        True
        >>> str(loaded_artifact._archiver.path) == \
                str(cache.data / str(artifact.uuid))
        True
        >>> test_dir.cleanup()
        """
        with self.lock:
            key_values = self.read_key(key)

            if 'data' not in key_values:
                raise ValueError(f"The key file '{key}' does not point to any "
                                 "data. This most likely occurred because you "
                                 "tried to load a pool which is not "
                                 "supported.")

            path = self.data / key_values['data']
            archiver = Archiver.load_raw(path, self)

        return Result._from_archiver(archiver)

    def load_collection(self, key):
        """Loads a pool referenced by a given key as a Collection. The pool
        loaded must have been created using Cache.save_collection.
        """
        pass

    def _load_uuid(self, uuid):
        """Load raw from the cache if the uuid is in the cache. Return None if
        it isn't. This is done so if someone already has an artifact in the
        cache then tries to use their qza for the artifact we can use the
        already cached version instead.
        """
        path = self.data / str(uuid)

        with self.lock:
            if os.path.exists(path):
                return Archiver.load_raw(path, self)
            else:
                return None

    def remove(self, key):
        """Removes a key from the cache then runs garbage collection to remove
        anything the removed key was referencing and any other loose data.

        Parameters
        ----------
        key : str
            The key we are removing.

        Raises
        ------
        KeyError
            If the key does not exist in the cache.

        Examples
        --------
        >>> from qiime2.sdk.result import Artifact
        >>> from qiime2.core.testing.type import IntSequence1
        >>> test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        >>> cache_path = os.path.join(test_dir.name, 'cache')
        >>> cache = Cache(cache_path)
        >>> artifact = Artifact.import_data(IntSequence1, [0, 1, 2])
        >>> saved_artifact = cache.save(artifact, 'key')
        >>> cache.get_keys() == set(['key'])
        True
        >>> cache.remove('key')
        >>> cache.get_keys() == set()
        True
        >>> # Note that the data is still in the cache due to our
        >>> # saved_artifact causing the process pool to keep a reference to it
        >>> cache.get_data() == set([str(saved_artifact.uuid)])
        True
        >>> del saved_artifact
        >>> # The data is still there even though the reference is gone because
        >>> # the cache has not run its own garbage collection yet. For various
        >>> # reasons, it is not feasible for us to safely garbage collect the
        >>> # cache when a reference in memory is deleted. Note also that
        >>> # "artifact" is not backed by the data in the cache, it only lives
        >>> # in memory, but it does have the same uuid as "saved_artifact."
        >>> cache.get_data() == set([str(artifact.uuid)])
        True
        >>> cache.garbage_collection()
        >>> # Now it is gone
        >>> cache.get_data() == set()
        True
        >>> test_dir.cleanup()
        """
        with self.lock:
            try:
                os.remove(self.keys / key)
            except FileNotFoundError as e:
                raise KeyError(f"The cache '{self.path}' does not contain the"
                               f" key '{key}'") from e

            self.garbage_collection()

    def clear_lock(self):
        """Clears the flufl lock on the cache. This exists in case something
        goes horribly wrong and we end up in an unrecoverable state. It's
        easy to tell the user "Recreate the failed cache (use the same path)
        and run this method on it."

        Note
        ----
        Forcibly removes the lock outside of the locking library's API.
        """
        if os.path.exists(self.lockfile):
            os.remove(self.lockfile)

    def _copy_to_data(self, ref):
        """If the data does not already exist in the cache, it will copy the
        data into the cache's data directory and set the appropriate
        permissions on the data. If the data does already exist in the cache,
        it will do nothing. This is generally used to copy data from outside
        the cache into the cache.

        Parameters
        ----------
        ref : Result
            The data we are copying into the cache's data directory.
        """
        destination = self.data / str(ref.uuid)

        with self.lock:
            if not os.path.exists(destination):
                # We need to actually create the cache/data/uuid directory
                # manually because the uuid isn't a part of the ArchivePath
                if not isinstance(ref._archiver.path, ArchivePath):
                    os.mkdir(destination)
                    shutil.copytree(
                        ref._archiver.path, destination, dirs_exist_ok=True)
                # Otherwise, the path we are copying should already contain the
                # uuid, so we don't need to manually create the uuid directory
                else:
                    shutil.copytree(
                        ref._archiver.path, self.data, dirs_exist_ok=True)

                set_permissions(destination, READ_ONLY_FILE, READ_ONLY_DIR)

    def _rename_to_data(self, uuid, src):
        """Takes some data in src and renames it into the cache's data dir. It
        then ensures there are symlinks for this data in the process pool and
        the named pool if one exists. This is generally used to move data from
        temporary per thread mount points in the process pool into the cache's
        data directory in one atomic action.

        Parameters
        ----------
        uuid : str or uuid4
            The uuid of the artifact whose data we are renaming into self.data
        src : str or Pathlike object
            The location of the data we are renaming into self.data.

        Returns
        -------
        str
            The alias we created for the artifact in the cache's process pool.
        pathlib.Path
            The location we renamed the data into.
        """
        uuid = str(uuid)

        dest = self.data / uuid
        alias = os.path.split(src)[0]
        with self.lock:
            # Rename errors if the destination already exists
            if not os.path.exists(dest):
                os.rename(src, dest)
                set_permissions(dest, READ_ONLY_FILE, READ_ONLY_DIR)

            # Create a new alias whether we renamed or not because this is
            # still loading a new reference to the data even if the data is
            # already there
            process_alias = self._alias(uuid)

        # Remove the aliased directory above the one we renamed. We need to do
        # this whether we renamed or not because we aren't renaming this
        # directory but the one beneath it
        shutil.rmtree(alias)
        return process_alias, dest

    def _alias(self, uuid):
        """Creates an alias and a symlink for the artifact with the given uuid
        in both the cache's process pool and its named pool if it has one.

        Parameters
        ----------
        uuid : str or uuid4
            The uuid of the artifact we are aliasing.

        Returns
        -------
        str
            The alias we created for the artifact.
        """
        with self.lock:
            process_alias = self.process_pool._alias(uuid)
            self.process_pool._make_symlink(uuid, process_alias)

            # Named pool links are not aliased
            if self.named_pool is not None:
                self.named_pool._make_symlink(uuid, uuid)

        return process_alias

    def _deallocate(self, symlink):
        """Removes a specific symlink from the process pool. This happens when
        an archiver goes out of scope. We remove that archiver's reference to
        the data from the process pool. We do this to prevent the cache from
        growing wildly during long running processes.

        Parameters
        ----------
        symlink : str
            The basename of the symlink we are going to be removing from the
            process pool.
        """
        # NOTE: Beware locking inside of this method. This method is called by
        # Python's garbage collector and that seems to cause deadlocks when
        # acquiring the thread lock
        target = self.process_pool.path / symlink

        if target.exists():
            os.remove(target)

    @property
    def data(self):
        """The directory in the cache that stores the data.
        """
        return self.path / 'data'

    def get_data(self):
        """Returns a set of all data in the cache.

        Returns
        -------
        set[str]
            All of the data in the cache in the form of the top level
            directories which will be the uuids of the artifacts.
        """
        with self.lock:
            return set(os.listdir(self.data))

    @property
    def keys(self):
        """The directory in the cache that stores the keys.
        """
        return self.path / 'keys'

    def get_keys(self):
        """Returns a set of all keys in the cache.

        Returns
        -------
        set[str]
            All of the keys in the cache. Just the names now what they refer
            to.
        """
        with self.lock:
            return set(os.listdir(self.keys))

    @property
    def lockfile(self):
        """The path to the flufl lock file.
        """
        return self.path / 'LOCK'

    @property
    def pools(self):
        """The directory in the cache that stores the named pools.
        """
        return self.path / 'pools'

    def get_pools(self):
        """Returns a set of all pools in the cache.

        Returns
        -------
        set[str]
            The names of all of the named pools in the cache.
        """
        with self.lock:
            return set(os.listdir(self.pools))

    @property
    def processes(self):
        """The directory in the cache that stores the process pools.
        """
        return self.path / 'processes'

    def get_processes(self):
        """Returns a set of all process pools in the cache.

        Returns
        -------
        set[str]
            The names of all of the process pools in the cache.
        """
        with self.lock:
            return set(os.listdir(self.processes))

    @property
    def version(self):
        """The path to the version file.
        """
        return self.path / 'VERSION'


class Pool:
    """Pools are folders in the cache that contain many symlinks to many
    different piece of data. There are two types of pool:

    **Process Pools:** These pools have names of the form
    <process-id>-<process-create-time>@<user> based on the process that created
    them. They only exist for the length of the process that created them and
    ensure data that process is using stays in the cache.

    **Named Pools:** Named pools are keyed just like individual pieces of data.
    They exist for as long as they have a key, and all of the data they symlink
    to is retained in the cache for as long as the pool exists.
    """

    def __init__(self, cache, name=None, reuse=False):
        """Used with name=None and reuse=True to create a process pool. Used
        with a name to create named pools.

        Note
        ----
        In general, you should not invoke this constructor directly and should
        instead use qiime2.core.cache.Cache.create_pool to create a pool
        properly on a given cache.

        Parameters
        ----------
        cache : Cache
            The cache this pool will be created under.
        named : str
            The name of the pool we are creating if it is a named pool.
        reuse : bool
            Whether we will be reusing this pool if it already exists.

        Raises
        ------
        ValueError
            If the pool already exists and reuse is False.
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
        # <process-id>-<process-start-time>@<user>
        else:
            self.name = self._get_process_pool_name()
            self.path = cache.processes / self.name

        # Raise a value error if we thought we were making a new pool but
        # actually are not
        if not reuse and os.path.exists(self.path):
            raise ValueError("Pool already exists, please use reuse=True to "
                             "reuse existing pool, or remove all keys "
                             "indicating this pool to remove the pool")

        if not os.path.exists(self.path):
            os.mkdir(self.path)

    def __enter__(self):
        """Tells the currently set cache to use this named pool. If there is no
        cache set then set the cache this named pool is on as well.

        Note
        ----
        If you have already set a cache then you cannot set a named pool that
        belongs to a different cache.

        Raises
        ------
        ValueError
            If you try to set a pool that is not on the currently set cache.
        ValueError
            If you have already set a pool and try to set another.

        Examples
        --------
        >>> test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        >>> cache_path = os.path.join(test_dir.name, 'cache')
        >>> cache = Cache(cache_path)
        >>> pool = cache.create_pool(keys=['pool'])
        >>> # When we with in the pool the set cache will be the cache the pool
        >>> # belongs to, and the named pool on that cache will be the pool
        >>> # we withed in
        >>> with pool:
        ...     current_cache = get_cache()
        ...     cache.named_pool == pool
        True
        >>> current_cache == cache
        True
        >>> # Now that we have exited the with, both cache and pool are unset
        >>> get_cache() == cache
        False
        >>> cache.named_pool == pool
        False
        >>> test_dir.cleanup()
        """
        if _CACHE.cache is not None and _CACHE.cache.path != self.cache.path:
            raise ValueError('Cannot enter a pool that is not on the '
                             'currently set cache. The current cache is '
                             f'located at: {_CACHE.cache.path}')
        else:
            self.previously_entered_cache = _CACHE.cache
            _CACHE.cache = self.cache

        if self.cache.named_pool is not None:
            raise ValueError("You cannot enter multiple pools at once, "
                             "currently entered pool is located at: "
                             f"'{self.cache.named_pool.path}'")

        self.cache.named_pool = self

    def __exit__(self, *args):
        """Unsets the named pool on the currently set cache. If there was no
        cache set before setting this named pool then unset the cache as well.

        Note
        ----
        self.previously_entered_cache will either be None or the cache this
        named pool belongs to. It will be None if there was no cache set when
        we set this named pool. It will be this named pool's cache if that
        cache was already set when we set this named pool. If there was a
        different cache set when we set this named pool, we would have errored
        in __enter__.
        """
        _CACHE.cache = self.previously_entered_cache
        self.cache.named_pool = None

    def _get_process_pool_name(self):
        """Creates a process pool name of the format
        <process-id>-<process-create-time>@<user> for the process that invoked
        this function.

        Returns
        -------
        str
            The name of this process pool.
        """
        pid = os.getpid()
        user = _get_user()

        process = psutil.Process(pid)
        time = process.create_time()

        return f'{pid}-{time}@{user}'

    def save(self, ref):
        """Saves the data into the pool then loads a new ref backed by the data
        in the pool.

        Parameters
        ----------
        ref : Result
            The QIIME 2 result we are saving into this pool.

        Returns
        -------
        Result
            A QIIME 2 result backed by the data in the cache the pool belongs
            to.

        Examples
        --------
        >>> from qiime2.sdk.result import Artifact
        >>> from qiime2.core.testing.type import IntSequence1
        >>> test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        >>> cache_path = os.path.join(test_dir.name, 'cache')
        >>> cache = Cache(cache_path)
        >>> pool = cache.create_pool(keys=['pool'])
        >>> artifact = Artifact.import_data(IntSequence1, [0, 1, 2])
        >>> pool_artifact = pool.save(artifact)
        >>> # The data itself resides in the cache this pool belongs to
        >>> str(pool_artifact._archiver.path) == \
                str(cache.data / str(artifact.uuid))
        True
        >>> # The pool now contains a symlink to the data. The symlink is named
        >>> # after the uuid of the data.
        >>> pool.get_data() == set([str(artifact.uuid)])
        True
        >>> test_dir.cleanup()
        """
        uuid = str(ref.uuid)
        if self.path == self.cache.process_pool.path:
            alias = self._alias(uuid)
        else:
            alias = uuid

        self._make_symlink(uuid, alias)

        self.cache._copy_to_data(ref)
        return self.load(ref)

    def _alias(self, uuid):
        """We may want to create multiple references to a single artifact in a
        process pool, but we cannot create multiple symlinks with the same
        name, so we take the uuid and add a random number to the end of it and
        use uuid.random_number as the name of the symlink. This means when you
        look in a process pool you may see the same uuid multiple times with
        different random numbers appended. This means there are multiple
        references to the artifact with that uuid in the process poole because
        it was loaded multiple times in the process.

        Parameters
        ----------
        uuid : str or uuid4
            The uuid we are creating an alias for.

        Returns
        -------
        str
            The aliased uuid.

        """
        MAX_RETRIES = 5

        uuid = str(uuid)
        with self.cache.lock:
            for _ in range(MAX_RETRIES):
                alias = uuid + '.' + str(randint(0, maxsize))
                path = self.path / alias

                # os.path.exists returns false on broken symlinks
                if not os.path.exists(path) and not os.path.islink(path):
                    break
            else:
                raise ValueError(f'Too many collisions ({MAX_RETRIES}) '
                                 'occurred while trying to save artifact '
                                 f'<{uuid}> to process pool {self.path}. It '
                                 'is likely you have attempted to load the '
                                 'same artifact a very large number of times.')
        return alias

    def _allocate(self, uuid):
        """Allocate an empty directory under the process pool to extract to.
        This directory is of the form alias / uuid and provides a per thread
        mount location for artifacts.

        Parameters
        ----------
        uuid : str or uuid4
            The uuid of the artifact we are creating an extract location for.

        Returns
        -------
        pathlib.Path
            The path we allocated to extract the artifact into.
        """
        uuid = str(uuid)

        # We want to extract artifacts to this thread unique location in the
        # process pool before using Cache.rename to put them into Cache.data.
        # We need to do this in order to ensure that if a uuid exists in
        # Cache.data, it is actually populated with data and is usable as an
        # artifact. Otherwise it could just be an empty directory (or only
        # contain part of the artifact) when another thread/process tries to
        # access it.
        with self.cache.lock:
            alias = self._alias(uuid)
            allocated_path = self.path / alias / uuid
            os.makedirs(allocated_path)

        return allocated_path

    def _make_symlink(self, uuid, alias):
        """Symlinks self.path / alias to self.cache.data / uuid. This creates a
        reference to the artifact with the given uuid in the cache.

        Parameters
        ----------
        uuid : str or uuid4
            The uuid of the artifact we are creating a symlink reference for.
        alias : str
            The alias we are using as the actual name of the symlink.
        """
        uuid = str(uuid)
        src = self.cache.data / uuid
        dest = self.path / alias

        # Symlink will error if the location we are creating the link at
        # already exists. This could happen legitimately from trying to save
        # the same thing to a named pool several times.
        with self.cache.lock:
            if not os.path.exists(dest):
                os.symlink(src, dest)

    def load(self, ref):
        """Loads a reference to an element in the pool.

        Parameters
        ----------
        ref : str or Result
            The result we are loading out of this pool, or just its uuid as a
            string.

        Returns
        -------
        Result
            A result backed by the data in the cache that this pool belongs to.

        Examples
        --------
        >>> from qiime2.sdk.result import Artifact
        >>> from qiime2.core.testing.type import IntSequence1
        >>> test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        >>> cache_path = os.path.join(test_dir.name, 'cache')
        >>> cache = Cache(cache_path)
        >>> pool = cache.create_pool(keys=['pool'])
        >>> artifact = Artifact.import_data(IntSequence1, [0, 1, 2])
        >>> pool_artifact = pool.save(artifact)
        >>> loaded_artifact = pool.load(str(artifact.uuid))
        >>> artifact == pool_artifact == loaded_artifact
        True
        >>> str(loaded_artifact._archiver.path) == \
                str(cache.data / str(artifact.uuid))
        True
        >>> test_dir.cleanup()
        """
        # Could receive an artifact or just a string uuid
        if isinstance(ref, str):
            uuid = ref
        else:
            uuid = str(ref.uuid)

        path = self.cache.data / uuid

        archiver = Archiver.load_raw(path, self.cache)
        return Result._from_archiver(archiver)

    def remove(self, ref):
        """Removes an element from the pool. The element can be just the uuid
        of the data as a string, or it can be a Result object referencing the
        data we are trying to remove.

        Parameters
        ----------
        ref : str or Result
            The result we are removing from this pool, or just its uuid as a
            string.

        Examples
        --------
        >>> from qiime2.sdk.result import Artifact
        >>> from qiime2.core.testing.type import IntSequence1
        >>> test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        >>> cache_path = os.path.join(test_dir.name, 'cache')
        >>> cache = Cache(cache_path)
        >>> pool = cache.create_pool(keys=['pool'])
        >>> artifact = Artifact.import_data(IntSequence1, [0, 1, 2])
        >>> pool_artifact = pool.save(artifact)
        >>> pool.get_data() == set([str(artifact.uuid)])
        True
        >>> pool.remove(str(artifact.uuid))
        >>> pool.get_data() == set()
        True
        >>> # Note that the data is still in the cache due to our
        >>> # pool_artifact causing the process pool to keep a reference to it
        >>> cache.get_data() == set([str(pool_artifact.uuid)])
        True
        >>> del pool_artifact
        >>> # The data is still there even though the reference is gone because
        >>> # the cache has not run its own garbage collection yet. For various
        >>> # reasons, it is not feasible for us to safely garbage collect the
        >>> # cache when a reference in memory is deleted. Note also that
        >>> # "artifact" is not backed by the data in the cache, it only lives
        >>> # in memory, but it does have the same uuid as "pool_artifact."
        >>> cache.get_data() == set([str(artifact.uuid)])
        True
        >>> cache.garbage_collection()
        >>> # Now it is gone
        >>> cache.get_data() == set()
        True
        >>> test_dir.cleanup()
        """
        # Could receive an artifact or just a string uuid
        if isinstance(ref, str):
            uuid = ref
        else:
            uuid = str(ref.uuid)

        target = self.path / uuid
        with self.cache.lock:
            if target.exists():
                if os.path.islink(target):
                    os.remove(target)
                else:
                    shutil.rmtree(target)
                self.cache.garbage_collection()

    def get_data(self):
        """Returns a set of all data in the pool.

        Returns
        -------
        set[str]
            The uuids of all of the data in the pool.
        """
        return set(os.listdir(self.path))
