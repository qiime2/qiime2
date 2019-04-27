# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from .collection import List, Set
from .semantic import SemanticType, Properties
from .primitive import (Str, Int, Float, Metadata, Bool, MetadataColumn,
                        Categorical, Numeric, Range, Start, End, Choices)
from .visualization import Visualization
from .signature import PipelineSignature, MethodSignature, VisualizerSignature
from .meta import TypeMap, TypeMatch
from .util import (is_primitive_type, is_semantic_type, is_metadata_type,
                   is_collection_type, is_visualization_type,
                   interrogate_collection_type, parse_primitive)

__all__ = [
    # Type Helpers
    'is_semantic_type', 'is_visualization_type', 'is_primitive_type',
    'is_metadata_type', 'is_collection_type', 'interrogate_collection_type',
    'parse_primitive',
    # Collection Types
    'Set', 'List',
    # Semantic Types
    'SemanticType',
    'Properties',
    # Primitive Types
    'Str', 'Int', 'Float', 'Bool', 'Metadata', 'MetadataColumn',
    'Categorical', 'Numeric', 'Range', 'Start', 'End', 'Choices',
    # Visualization Type
    'Visualization',
    # Signatures
    'PipelineSignature', 'MethodSignature', 'VisualizerSignature',
    # Variables
    'TypeMap', 'TypeMatch'
]
