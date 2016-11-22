# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from importlib import import_module

import qiime
import qiime.plugin

from .format import (
    IntSequenceFormat,
    MappingFormat,
    SingleIntFormat,
    IntSequenceDirectoryFormat,
    MappingDirectoryFormat,
    FourIntsDirectoryFormat
)

from .type import IntSequence1, IntSequence2, Mapping, FourInts
from .method import (concatenate_ints, split_ints, merge_mappings,
                     identity_with_metadata, identity_with_metadata_category)
from .visualizer import most_common_viz, mapping_viz

dummy_plugin = qiime.plugin.Plugin(
    name='dummy-plugin',
    version='0.0.0-dev',
    website='https://github.com/qiime2/qiime2',
    package='qiime.core.testing',
    citation_text='No relevant citation.',
    user_support_text='For help, see https://qiime2.org'
)

import_module('qiime.core.testing.transformer')

# Register semantic types
dummy_plugin.register_semantic_types(IntSequence1, IntSequence2, Mapping,
                                     FourInts)

# Register formats
dummy_plugin.register_formats(
    IntSequenceFormat, MappingFormat, SingleIntFormat,
    IntSequenceDirectoryFormat, MappingDirectoryFormat,
    FourIntsDirectoryFormat
)


dummy_plugin.register_semantic_type_to_format(
    IntSequence1,
    artifact_format=IntSequenceDirectoryFormat
)
dummy_plugin.register_semantic_type_to_format(
    IntSequence2,
    artifact_format=IntSequenceDirectoryFormat
)
dummy_plugin.register_semantic_type_to_format(
    Mapping,
    artifact_format=MappingDirectoryFormat
)
dummy_plugin.register_semantic_type_to_format(
    FourInts,
    artifact_format=FourIntsDirectoryFormat
)

# TODO add an optional parameter to this method when they are supported
dummy_plugin.methods.register_function(
    function=concatenate_ints,
    inputs={
        'ints1': IntSequence1 | IntSequence2,
        'ints2': IntSequence1,
        'ints3': IntSequence2
    },
    parameters={
        'int1': qiime.plugin.Int,
        'int2': qiime.plugin.Int
    },
    outputs=[
        ('concatenated_ints', IntSequence1)
    ],
    name='Concatenate integers',
    description='This method concatenates integers into a single sequence in '
                'the order they are provided.'
)

# TODO update to use TypeMap so IntSequence1 | IntSequence2 are accepted, and
# the return type is IntSequence1 or IntSequence2.
dummy_plugin.methods.register_function(
    function=split_ints,
    inputs={
        'ints': IntSequence1
    },
    parameters={},
    outputs=[
        ('left', IntSequence1),
        ('right', IntSequence1)
    ],
    name='Split sequence of integers in half',
    description='This method splits a sequence of integers in half, returning '
                'the two halves (left and right). If the input sequence\'s '
                'length is not evenly divisible by 2, the right half will '
                'have one more element than the left.'
)

dummy_plugin.methods.register_function(
    function=merge_mappings,
    inputs={
        'mapping1': Mapping,
        'mapping2': Mapping
    },
    parameters={},
    outputs=[
        ('merged_mapping', Mapping)
    ],
    name='Merge mappings',
    description='This method merges two mappings into a single new mapping. '
                'If a key is shared between mappings and the values differ, '
                'an error will be raised.'
)

dummy_plugin.methods.register_function(
    function=identity_with_metadata,
    inputs={
        'ints': IntSequence1 | IntSequence2
    },
    parameters={
        'metadata': qiime.plugin.Metadata
    },
    outputs=[
        ('out', IntSequence1)
    ],
    name='Identity',
    description='This method does nothing, but takes metadata'
)

dummy_plugin.methods.register_function(
    function=identity_with_metadata_category,
    inputs={
        'ints': IntSequence1 | IntSequence2
    },
    parameters={
        'metadata': qiime.plugin.MetadataCategory
    },
    outputs=[
        ('out', IntSequence1)
    ],
    name='Identity',
    description='This method does nothing, but takes a metadata category'
)

dummy_plugin.visualizers.register_function(
    function=most_common_viz,
    inputs={
        'ints': IntSequence1 | IntSequence2
    },
    parameters={},
    name='Visualize most common integers',
    description='This visualizer produces HTML and TSV outputs containing the '
                'input sequence of integers ordered from most- to '
                'least-frequently occurring, along with their respective '
                'frequencies.'
)

# TODO add optional parameters to this method when they are supported
dummy_plugin.visualizers.register_function(
    function=mapping_viz,
    inputs={
        'mapping1': Mapping,
        'mapping2': Mapping
    },
    parameters={
        'key_label': qiime.plugin.Str,
        'value_label': qiime.plugin.Str
    },
    name='Visualize two mappings',
    description='This visualizer produces an HTML visualization of two '
                'key-value mappings, each sorted in alphabetical order by key.'
)
