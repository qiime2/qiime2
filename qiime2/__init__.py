# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.sdk import Artifact, Visualization
from qiime2.metadata import (Metadata, MetadataColumn,
                             CategoricalMetadataColumn, NumericMetadataColumn)
from ._version import get_versions

__version__ = get_versions()['version']
del get_versions

# "Train release" version includes <year>.<month> and excludes patch numbers
# and pre/post-release tags. All versions within a train release are expected
# to be compatible.
__release__ = '.'.join(__version__.split('.')[:2])

__all__ = ['Artifact', 'Visualization', 'Metadata', 'MetadataColumn',
           'CategoricalMetadataColumn', 'NumericMetadataColumn']
