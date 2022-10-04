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

By default, a cache will be created under tmpdir/qiime2/<uname> and all
intermediate data that was previously being written all over tmpdir will be
written into that specific directory. That means QIIME 2 reserves usage of the
tmpdir/qiime2 directory. The user may also specify a new location to be used
in place of this default directory. This location must meet a few criteria

1. It must be writable from any and all locations the QIIME 2 command intending
to use it will be running. This means that in an HPC context, the location
specified for the cache must be writable from the node QIIME 2 will be
executing on

2. It must either not exist or already be a cache. The first time a directory
is specified to be used as a cache, it should not exist. QIIME 2 will create a
cache structure on disk at that location. Any existing directory you attempt to
use as a cache should have been created as a cache by QIIME 2.

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

from flufl.lock import Lock

import qiime2
from .path import ArchivePath
from qiime2.sdk.result import Result
from qiime2.core.util import (is_uuid4, set_permissions, touch_under_path,
                              READ_ONLY_FILE, READ_ONLY_DIR, ALL_PERMISSIONS)
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

# Thread local indicating the cache to use
_CACHE = threading.local()
_CACHE.cache = None
_CACHE.temp_cache = None

# TODO: Do we want this on the threadlocal? I feel like maybe we do
# Keep track of every cache used by this process for cleanup later
USED_CACHES = set()

# These permissions are directory with sticky bit and rwx for all set
EXPECTED_PERMISSIONS = 0o41777


def get_cache():
    """Gets the cache we have instructed QIIME 2 to use in this invocation.
    By default this is a cache located at tmpdir/qiime2/<uname>, but if the
    user has set a cache it is the cache they set.

    Returns
    -------
    Cache
        The cache QIIME 2 is using for the current invocation.
    """
    # If we are on a new thread we may in fact not have a cache attribute here
    # at all
    if not hasattr(_CACHE, 'cache') or _CACHE.cache is None:
        if not hasattr(_CACHE, 'temp_cache') or _CACHE.temp_cache is None:
            _CACHE.temp_cache = Cache()
        return _CACHE.temp_cache

    return _CACHE.cache


# TODO: maybe hand shutil.copytree qiime2.util.duplicate
def _copy_to_data(cache, ref):
    """If the data does not already exist in the cache, it will copy the data
    into the cache's data directory and set the appropriate permissions on the
    data. If the data does already exist in the cache, it will do nothing.

    Parameters
    ----------
    cache : Cache
        The cache whose data directory we are moving data into.
    ref : Result
        The data we are copying into the cache's data directory.
    """
    destination = cache.data / str(ref.uuid)

    if not os.path.exists(destination):
        if not isinstance(ref._archiver.path, ArchivePath):
            os.mkdir(destination)
            shutil.copytree(ref._archiver.path, destination,
                            dirs_exist_ok=True)
        else:
            shutil.copytree(ref._archiver.path, cache.data, dirs_exist_ok=True)

        set_permissions(destination, READ_ONLY_FILE, READ_ONLY_DIR)


def _get_user():
    """Get the <uname> for our default cache. Internally getpass.getuser is
    getting the uid then looking up the username associated with it. This could
    fail it we are running inside a container because the container is looking
    for its parent's uid in its own /etc/passwd which is unlikely to contain a
    user associated with that uid. If that failure does occur, we create an
    alternate default username.

    Returns
    -------
    str
        The value we will be using as <uname> for our default cache.
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
    ├── processes
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
    base_cache_contents = \
        set(('data', 'keys', 'pools', 'processes', 'VERSION'))

    def __init__(self, path=None, process_pool_lifespan=45):
        """Creates a cache object backed by the directory specified by path. If
        no path is provided it gets a path to a temp cache.

        Parameters
        ----------
        path : str or PathLike object
            Should point either to a non-existent writable directory to be
            created as a cache or to an existing writable cache. Defaults to
            None which creates the cache at tmpdir/qiime2/<uname>
        process_pool_lifespan : int
            The number of days we should allow process_pools to exist for
            before culling them.

        Examples
        --------
        >>> test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        >>> cache_path = os.path.join(test_dir.name, 'cache')
        >>> cache = Cache(cache_path)
        >>> Cache.is_cache(cache_path)
        True
        >>> test_dir.cleanup()
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
            # MacOS culls files in the temp dir that haven't been used for a
            # few days. This can lead to the VERSION file being deleted while
            # we still have a cache dir, so we see the directory but don't
            # think it's a cache. Our solution is to just kill this directory
            # and recreate it. We only do this on the temp cache which is not
            # storing anything long term anyway
            if path is None:
                warnings.warn("Your temporary cache was found to be in an "
                              "inconsistent state. It has been recreated.")
                set_permissions(self.path, ALL_PERMISSIONS, ALL_PERMISSIONS)
                shutil.rmtree(self.path)
                self._create_cache()
            else:
                raise ValueError(
                    f"Path: \'{self.path}\' already exists and is not a "
                    "cache.")

        self.lock = Lock(str(self.lockfile), lifetime=timedelta(minutes=10))
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
        """Tell QIIME 2 to use this cache on its current invocation.
        """
        if _CACHE.cache is not None and _CACHE.cache.path != self.path:
            raise ValueError("You cannot enter multiple caches at once, "
                             "currently entered cache is located at: "
                             f"'{_CACHE.cache.path}'")

        _CACHE.cache = self

    def __exit__(self, *args):
        """Tell QIIME 2 to stop using this cache and go back to using whatever
        cache it was using before.
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

    def _create_cache(self):
        """Create the cache directory, all sub directories, and the version
        file.
        """
        # Construct the cache root recursively
        os.makedirs(self.path)
        os.mkdir(self.data)
        os.mkdir(self.keys)
        os.mkdir(self.pools)
        os.mkdir(self.processes)
        # Do we want this right off the bat? How exactly is setting tmp in the
        # cache going to work? tmp is never going to be managed by the cache,
        # it's just so they're both on the same disk, so they'll probably just
        # set the tmp location in the config or something. I feel like if we're
        # going to manage the cache, we should manage the cache which means if
        # they're going to create_pool put tmp in the cache it should have to
        # be in a set directory within the cache like tmp not just whatever
        # they want it to be in the cache. Not sure how we would really enforce
        # that, but we can just... Heavily encourage it I guess
        # os.mkdir('tmp')

        self.version.write_text(
            _VERSION_TEMPLATE % (self.CURRENT_FORMAT_VERSION,
                                 qiime2.__version__))

    def _get_temp_path(self):
        """Get path to temp cache if the user did not specify a named cache.
        This function will create the path if it does not exist and ensure it
        is suitable for use as a cache if it does.

        Warning
        -------
        If the path tmpdir/qiime2/<uname> exists but is not a valid cache, we
        remove the directory and create a cache there.

        Returns
        -------
        str
            The path created for the temp cache.
        """
        tmpdir = tempfile.gettempdir()

        cache_dir = os.path.join(tmpdir, 'qiime2')

        # Make sure the sticky bit is set on the cache directory. Documentation
        # on what a sitcky bit is found here
        # https://docs.python.org/3/library/stat.html#stat.S_ISVTX
        # We also set read/write/execute permissions for everyone on this
        # directory. We only do this if we are the owner of the /tmp/qiime2
        # directory or in other words the first person to run QIIME 2 with this
        # /tmp since the /tmp was wiped
        if not os.path.exists(cache_dir):
            os.mkdir(cache_dir)
            sticky_permissions = stat.S_ISVTX | stat.S_IRWXU | stat.S_IRWXG \
                | stat.S_IRWXO
            os.chmod(cache_dir, sticky_permissions)
        elif os.stat(cache_dir).st_mode != EXPECTED_PERMISSIONS:
            raise ValueError(f"Directory '{cache_dir}' already exists without "
                             f"proper permissions "
                             f"'{oct(EXPECTED_PERMISSIONS)}' set. Current "
                             "permissions are "
                             f"'{oct(os.stat(cache_dir).st_mode)}.' This most "
                             "likely means something other than QIIME 2 "
                             f"created the directory '{cache_dir}' or QIIME 2 "
                             f"failed between creating '{cache_dir}' and "
                             "setting permissions on it.")

        user = _get_user()
        user_dir = os.path.join(cache_dir, user)

        # It is conceivable that we already have a path matching this username
        # that belongs to another uid, if we do then we want to create a
        # garbage name for the temp cache that will be used by this user
        if os.path.exists(user_dir) and \
                os.stat(user_dir).st_uid != os.getuid():
            uid_name = _get_uid_cache_name()
            # This really shouldn't happen
            if user == uid_name:
                raise ValueError(f'Temp cache for uid path {user} already '
                                 'exists but does not belong to us.')

            user_dir = os.path.join(cache_dir, uid_name)

        return user_dir

    def _create_process_pool(self):
        """Creates a process pool which is identical in function to a named
        pool, but it lives in the processes subdirectory not the pools
        subdirectory, and is handled differently by garbage collection due to
        being unkeyed. Process pools are used to keep track of results for
        currently running processes and are removed when the process that
        created them ends.

        Returns
        -------
        Pool
            The pool we created.
        """
        return Pool(self, reuse=True)

    def create_pool(self, keys=[], reuse=False):
        """Used to create named pools. A named pool's name is all of the keys
        given for it separated by underscores.

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
            The pool we created
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
            The name of the pool we are keying
        keys : List[str]
            A list of all the keys to create referring to the pool
        """
        for key in keys:
            self._register_key(key, pool_name, pool=True)

    def garbage_collection(self):
        """Runs garbage collection on the cache. We log all data and pools
        pointed to by keys. Then we go through all pools and delete any that
        were not referred to by a key while logging all data in pools that are
        referred to by keys. Then we go through all process pools and log all
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
                with open(self.keys / key) as fh:
                    loaded_key = yaml.safe_load(fh)
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
            for process_pool in os.listdir(self.processes):
                # Pick the creation time out of the pool name of format
                # {pid}-time@user
                create_time = float(process_pool.split('-')[1].split('@')[0])

                if time.time() - create_time >= self.process_pool_lifespan:
                    shutil.rmtree(self.processes / process_pool)
                else:
                    for data in os.listdir(self.processes / process_pool):
                        referenced_data.add(data.split('.')[0])

            # Walk over all data and remove any that was not referenced
            for data in os.listdir(self.data):
                # If this assert is ever tripped something real bad happened
                assert is_uuid4(data)

                if data not in referenced_data:
                    target = self.data / data

                    set_permissions(target, None, ALL_PERMISSIONS)
                    shutil.rmtree(target)

    def save(self, ref, key):
        """Save data into the cache by creating a key referring to the data
        then copying the data if it is not already in the cache. We create the
        key first because if we created the data first it would be unkeyed for
        a brief period of time and if someone else were garbage collecting the
        cache they would remove our unkeyed data between its creation and the
        creation of its key.

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
        """
        # Create the key before the data, this is so that if another thread or
        # process is running garbage collection it doesn't see our unkeyed data
        # and remove it leaving us with a dangling reference and no data
        with self.lock:
            self._register_key(key, str(ref.uuid))

        _copy_to_data(self, ref)
        return self.load(key)

    def _register_key(self, key, value, pool=False):
        """Create a key file pointing at the specified data or pool.

        Parameters
        ----------
        key : str
            The name of the key to create
        value : str
            The path to the data or pool we are keying
        pool : bool
            Whether we are keying a pool or not

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

        if pool:
            key_fp.write_text(_KEY_TEMPLATE % (key, '', value))
        else:
            key_fp.write_text(_KEY_TEMPLATE % (key, value, ''))

    def load(self, key):
        """Load the data pointed to by a key. Only works on a key that refers
        to a data item will error on a key that points to a pool

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
        ValueError
            If the cache does not contain the specified key.
        """
        try:
            with open(self.keys / key) as fh:
                path = self.data / yaml.safe_load(fh)['data']
        except TypeError as e:
            raise ValueError(f"The key file '{key}' does not point to any "
                             "data. This most likely occurred because you "
                             "tried to load a pool which is not supported.") \
                from e
        except FileNotFoundError as e:
            raise ValueError(f"The cache '{self.path}' does not contain the "
                             f"key '{key}'") from e

        archiver = Archiver.load_raw(path, self)
        return Result._from_archiver(archiver)

    def remove(self, key):
        """Remove a key from the cache then run garbage collection to remove
        anything it was referencing and any other loose data.

        Parameters
        ----------
        key : str
            The key we are removing.
        """
        os.remove(self.keys / key)
        self.garbage_collection()

    def clear_lock(self):
        """Clears the flufl lock on the cache.

        Note
        ----
        Forcibly removes the lock outside of the locking library's API.
        """
        if os.path.exists(self.lockfile):
            os.remove(self.lockfile)

    def _allocate(self, uuid):
        """Creates a space in the cache for the archiver to put its data. Also
        creates symlinks in the process pool and named pool (if there is one)
        to keep track of the data. These symlinks will contain the uuid of the
        data with some random bits appended if the symlinks were already
        present to create a new reference without a duplicate name.

        Parameters
        ----------
        uuid : str
            The uuid of the data we are going to be saving.

        Returns
        -------
        pathlib.Path
            The path to the allocated space in the data directory.
        pathlib.Path
            The name of the symlink to the data in the process pool.
        """
        process_alias = self._alias(uuid)
        path = self.data / uuid

        if not os.path.exists(path):
            os.mkdir(path)

        return path, process_alias

    def _alias(self, uuid):
        """Aliases a uuid with some extra random characters at the end to allow
        us to have multiple references to the same data in one process pool.

        Parameters
        ----------
        uuid : str
            The uuid of the data we are going to be saving.

        Returns
        -------
        str
            The aliased version of the uuid
        """
        process_alias = self.process_pool._make_symlink(uuid)

        if self.named_pool is not None:
            self.named_pool._make_symlink(uuid)

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
            process pool
        """
        target = self.process_pool.path / symlink

        if target.exists():
            os.remove(target)

    @property
    def data(self):
        """The directory in the cache that stores the data.
        """
        return self.path / 'data'

    @property
    def keys(self):
        """The directory in the cache that stores the keys.
        """
        return self.path / 'keys'

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

    @property
    def processes(self):
        """The directory in the cache that stores the process pools.
        """
        return self.path / 'processes'

    @property
    def version(self):
        """The path to the version file.
        """
        return self.path / 'VERSION'


# Assume we will make this its own class for now
class Pool:
    def __init__(self, cache, name=None, reuse=False):
        """Used with name=None and reuse=True to create a process pool. Used
        with a name to create named pools.

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
        # pid-process_start_time@user
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
        """Tell the currently set cache to use this named pool. If there is no
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
        """Unset the named pool on the currently set cache. If there was no
        cache set before setting this named pool then unset the cache as well.

        Note
        ----
        self.previously_entered_cache will either be None or the cache this
        named pool belongs to. It will be None if there was no cache set when
        we set this named pool. It will be this named pool's cache if that
        cache was already set when we set this named pool. If there was a
        different cache set when we set this named pool, we would have errored
        in __enter__
        """
        _CACHE.cache = self.previously_entered_cache
        self.cache.named_pool = None

    def _get_process_pool_name(self):
        """Creates a process pool name of the format <pid>-<start_time>@<user>
        for the process that invoked this function.

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
        """Save the data into the pool then load a new ref backed by the data
        in the pool.

        Parameters
        ----------
        ref : Result
            The QIIME 2 result we are saving into this pool.

        Returns
        -------
        Result
            A QIIME 2 result backed by the data in the pool.
        """
        self._make_symlink(str(ref.uuid))

        _copy_to_data(self.cache, ref)
        return self.load(ref)

    def _make_symlink(self, uuid):
        """Creates a symlink in this pool to the data referred to by the given
        uuid in the data path of the cache this pool belongs to. It is possible
        (especially if this is a process pool) that we already contain a
        symlink to the given data and therefore already contain a symlink named
        after the given uuid. If this happens, we will create a new symlink
        in the form of <uuid>.<random_characters>.

        Parameters
        ----------
        uuid : str
            The uuid of the result we are saving into the pool.

        Returns
        -------
        pathlib.Path
            The path to the symlink within this pool.

        Raises
        ------
        ValueError
            It took too many attempts to find a unique name for the symlink.
            This should really never happen the odds are so small.
        """
        MAX_RETRIES = 5

        src = self.cache.data / uuid
        dest = self.path / uuid

        with self.cache.lock:
            if self._guarded_symlink(src, dest):
                pass
            elif self.path == self.cache.process_pool.path:
                for _ in range(MAX_RETRIES):
                    new_uuid = uuid + '.' + str(randint(0, maxsize))
                    dest = self.path / new_uuid

                    if self._guarded_symlink(src, dest):
                        break
                else:
                    raise ValueError(f'Too many collisions ({MAX_RETRIES}) '
                                     'occurred while trying to save artifact '
                                     f'<{uuid}> to process pool {self.path}.'
                                     'It is likely you have attempted to load '
                                     'the same artifact a very large number '
                                     'of times.')
        return dest

    def _guarded_symlink(self, src, dest):
        """Creates a symlink at dest pointing to src only if the symlink does
        not already exist.

        Parameters
        ----------
        src : pathlib.Path
            The location of the data we are symlinking too.
        dest : pathlib.Path
            The location of the symlink we are creating.

        Returns
        -------
        bool
            Whether we were able to create the symlink or not.
        """
        if not os.path.exists(dest):
            os.symlink(src, dest)
            return True

        return False

    def load(self, ref):
        """Load a reference to an element in the pool.

        Parameters
        ----------
        ref : Result
            The result we are loading out of this pool.

        Returns
        -------
        Result
            A result backed by the data in the cache that this pool belongs to.
        """
        path = self.cache.data / str(ref.uuid)

        archiver = Archiver.load_raw(path, self.cache)
        return Result._from_archiver(archiver)

    def remove(self, ref):
        """Remove an element from the pool. The element can be just the uuid of
        the data as a string, or it can be a Result object referencing the data
        we are trying to remove.

        Parameters
        ----------
        ref : str or Result
            The result we are removing from this pool.
        """
        # Could receive an artifact or just a uuid
        if isinstance(ref, str):
            uuid = ref
        else:
            uuid = str(ref.uuid)

        # TODO: This guard should be removed when we rework the logic
        target = self.path / uuid
        if target.exists():
            os.remove(target)
            self.cache.garbage_collection()
