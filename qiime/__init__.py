# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import qiime.core  # noqa
from qiime.citation import Citation
from qiime.metadata import Metadata, MetadataCategory
import qiime.sdk as _sdk

__all__ = ['Citation', 'Metadata', 'MetadataCategory', 'Visualization',
           'Artifact']

__version__ = '2.0.3.dev'

# `from qiime import Artifact` fails if `from qiime.sdk` is used above so
# import and alias instead
Visualization = _sdk.Visualization
Artifact = _sdk.Artifact
