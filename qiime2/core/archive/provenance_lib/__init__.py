# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
"""
Software to support scientific reproducibility, attribution, and
collaboration on the QIIME 2 platform.

Core objects:
- ProvDAG: A directed, acyclic graph (DAG) describing QIIME 2 provenance
- ProvNode: Parsed data about a single QIIME 2 Result, generated internally by
  a provenance parser and available through the ProvDAG
"""

from .parse import ProvDAG, archive_not_parsed
from .replay import (
    replay_provenance, replay_citations, replay_supplement,
)
from .util import get_root_uuid, get_nonroot_uuid
from .usage_drivers import ReplayPythonUsage
from .tests.testing_utilities import TestArtifacts

__all__ = [
    'ProvDAG', 'archive_not_parsed', 'get_root_uuid', 'get_nonroot_uuid',
    'replay_provenance', 'replay_citations', 'replay_supplement',
    'ReplayPythonUsage', 'TestArtifacts'
]
