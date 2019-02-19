# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from .collection import List, Set, is_collection_type
from .semantic import SemanticType, is_semantic_type, Properties
from .primitive import (Str, Int, Float, Metadata, Bool, MetadataColumn,
                        Categorical, Numeric, Range, Start, End, Choices,
                        is_primitive_type)
from .visualization import Visualization
from .signature import PipelineSignature, MethodSignature, VisualizerSignature

__all__ = [
    # Type Helpers
    'is_semantic_type', 'is_primitive_type', 'is_collection_type',
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
    'PipelineSignature', 'MethodSignature', 'VisualizerSignature'
]
