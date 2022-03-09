# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import yaml
import pathlib

import qiime2
from qiime2.sdk.result import Artifact

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
    ├── tmp
    └── VERSION
    """
    VERSION_FILE = 'VERSION'
    CURRENT_FORMAT_VERSION = '1'

    def __init__(self, path):
        self.path = pathlib.Path(path)

        # Do we want a more rigorous check for whether or not we've been
        # pointed at an existing cache?
        if not os.path.exists(self.path):
            self.create_cache()
        elif not self.is_cache():
            raise ValueError(f"Path: \'{path}\' already exists and is not a"
                             " cache")

    # Surely this needs to be a thing? I suppose if they hand us a path that
    # doesn't exist we just create a cache there? do we want to create it at
    # that exact path do we slap a 'cache' sub-directory at that location and
    # use that?
    def create_cache(self):
        os.mkdir(self.path)
        os.mkdir(self.path / 'data')
        os.mkdir(self.path / 'keys')
        os.mkdir(self.path / 'pools')
        # Do we want this right off the bat? How exactly is setting tmp in the
        # cache going to work? tmp is never going to be managed by the cache,
        # it's just so they're both on the same disk, so they'll probably just
        # set the tmp location in the config or something. I feel like if we're
        # going to manage the cache, we should manage the cache which means if
        # they're going to put tmp in the cache it should have to be in a set
        # directory within the cache like tmp not just whatever they want it to
        # be in the cache. Not sure how we would really enforce that, but we
        # can just... Heavily encourage it I guess
        # os.mkdir('tmp')

        self.version.write_text(
            _VERSION_TEMPLATE % (self.CURRENT_FORMAT_VERSION,
                                 qiime2.__version__))

    # Tell us if the path is a cache or not
    def is_cache(self):
        pass

    # Run the garbage collection algorithm
    def garbage_collection(self):
        # Walk over keys and track all pools and data referenced
        #
        # Walk over pools and remove any that were not refered to by keys while
        # tracking all data within those that were referenced
        #
        # Walk over all data and remove any that was not referenced
        pass

    # Export artifact to zip
    def export(self, key):
        pass

    # Save artifact to key in cache
    def save(self, artifact, key):
        data_fp = str(self.data / str(artifact.uuid))
        artifact.save(data_fp)

        # Right now we aren't worrying about pools at all
        key_fp = self.keys / key
        key_fp.write_text(
            _KEY_TEMPLATE % (key, data_fp + artifact.extension, ''))

    # Artifact the load the data pointed to by the key. Does not work on pools.
    # Only works if you have data
    def load(self, key):
        return Artifact.load(yaml.safe_load(open(self.keys / key))['data'])

    # Remove key from cache
    def delete(self, key):
        pass

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
    # doing because we have a lock. Could we conceivably used built in Python
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
    def version(self):
        return self.path / self.VERSION_FILE
