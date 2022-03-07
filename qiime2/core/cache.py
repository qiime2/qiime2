# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

class Cache:
    """General structure of the cache
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
    └── version.txt
    """
    VERSION_FILE = 'VERSION'

    def __init__(self, path):
        self.path = path

    # Surely this needs to be a thing? I suppose if they hand us a path that
    # doesn't exist we just create a cache there? do we want to create it at
    # that exact path do we slap a 'cache' sub-directory at that location and
    # use that?
    def create_cache(self):
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

    # Remove key and backing from cache
    def remove(self, key):
        pass

    # Not entirely clear how this will work yet. We are assuming multi
    # processes from multiple systems will be interacting with the cache. This
    # means we can't even safely assume unique PIDs. We will probably create
    # some kind of lock file to lock the entire cache or to list locked
    # elements of the cache or something
    def lock(self):
        pass
