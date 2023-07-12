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

Core commands:
- replay_supplement
- replay_provenance
- replay_citations

These take a payload including QIIME 2 provenance information,
and write replay documentation to disk.

Core objects:
- ProvDAG: A directed, acyclic graph (DAG) describing QIIME 2 provenance
- ProvNode: Parsed data about a single QIIME 2 Result, generated internally by
  a provenance parser and available through the ProvDAG
"""

from .parse import ProvDAG, archive_not_parsed, UnparseableDataError
from .replay import (
    replay_provenance, replay_citations, replay_supplement,
)
from .util import get_root_uuid, get_nonroot_uuid, camel_to_snake
from .version_parser import parse_version

__all__ = [
    'ProvDAG', 'archive_not_parsed', 'UnparseableDataError',
    'get_root_uuid', 'get_nonroot_uuid', 'camel_to_snake',
    'replay_provenance', 'replay_citations', 'replay_supplement',
    'parse_version',
]
