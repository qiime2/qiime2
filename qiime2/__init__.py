# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources

import qiime2.core  # noqa
from qiime2.metadata import Metadata, MetadataCategory
import qiime2.sdk as _sdk

__all__ = ['Metadata', 'MetadataCategory', 'Visualization', 'Artifact']

__version__ = pkg_resources.get_distribution('qiime2').version

# "Train release" version includes <year>.<month> and excludes patch numbers
# and pre/post-release tags. All versions within a train release are expected
# to be compatible.
__release__ = '.'.join(__version__.split('.')[:2])

# `from qiime2 import Artifact` fails if `from qiime2.sdk` is used above so
# import and alias instead
Visualization = _sdk.Visualization
Artifact = _sdk.Artifact
