# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.sdk import Artifact, Visualization, ResultCollection
from qiime2.metadata import (Metadata, MetadataColumn,
                             CategoricalMetadataColumn, NumericMetadataColumn)
from qiime2.plugin import Citations
from qiime2.core.cache import Cache, Pool
from ._version import get_versions

__version__ = get_versions()['version']
del get_versions

# "Train release" version includes <year>.<month> and excludes patch numbers
# and pre/post-release tags. All versions within a train release are expected
# to be compatible.
__release__ = '.'.join(__version__.split('.')[:2])
__citations__ = tuple(Citations.load('citations.bib', package='qiime2'))
__website__ = 'https://qiime2.org'

__all__ = ['Artifact', 'Visualization', 'ResultCollection', 'Metadata',
           'MetadataColumn', 'CategoricalMetadataColumn',
           'NumericMetadataColumn', 'Cache', 'Pool']


# Used by `jupyter serverextension enable`
def _jupyter_server_extension_paths():
    return [{"module": "qiime2.jupyter"}]
