# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import stat
import threading

from qiime2.core.cache import Cache


CACHE_CONFIG = threading.local()

# Cache might be named or default
CACHE_CONFIG.cache = None
# The named pool that may or may not be set
CACHE_CONFIG.named_pool = None


def get_cache(path=None, named=True):
    """ TODO: Create some system for loading a cache config similar to a parsl
    config?
    """
    if not hasattr(CACHE_CONFIG.cache):
        CACHE_CONFIG.cache = Cache(path, named)

    return CACHE_CONFIG.cache


class CacheConfig():
    def __init__(self, path, named):
        """Do I even want this?
        """
        pass
