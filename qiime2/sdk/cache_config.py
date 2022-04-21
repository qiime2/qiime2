# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import threading

CACHE_CONFIG = threading.local()

# Cache might be named or default
CACHE_CONFIG.cache = None
# The named pool that may or may not be set
CACHE_CONFIG.named_pool = None
