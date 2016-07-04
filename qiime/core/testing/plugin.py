# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections

import qiime
import qiime.plugin

from .data_layout import (
    int_sequence_validator,
    int_sequence_to_list,
    int_sequence_to_counter,
    list_to_int_sequence,

    mapping_validator,
    mapping_to_dict,
    dict_to_mapping,

    four_ints_validator,
    four_ints_to_list,
    list_to_four_ints
)
from .type import IntSequence1, IntSequence2, Mapping, FourInts
from .method import concatenate_ints, split_ints, merge_mappings
from .visualizer import most_common_viz, mapping_viz


dummy_plugin = qiime.plugin.Plugin(
    name='dummy-plugin',
    version='0.0.0-dev',
    website='https://github.com/qiime2/qiime2',
    package='qiime.core.testing'
)

dummy_plugin.register_data_layout('int-sequence', 1, int_sequence_validator)

dummy_plugin.register_data_layout_reader('int-sequence', 1, list,
                                         int_sequence_to_list)

dummy_plugin.register_data_layout_writer('int-sequence', 1, list,
                                         list_to_int_sequence)

dummy_plugin.register_data_layout_reader(
    'int-sequence', 1, collections.Counter, int_sequence_to_counter)

dummy_plugin.register_data_layout('mapping', 1, mapping_validator)

dummy_plugin.register_data_layout_reader('mapping', 1, dict, mapping_to_dict)

dummy_plugin.register_data_layout_writer('mapping', 1, dict, dict_to_mapping)

dummy_plugin.register_data_layout('four-ints', 1, four_ints_validator)

dummy_plugin.register_data_layout_reader(
    'four-ints', 1, list, four_ints_to_list)

dummy_plugin.register_data_layout_writer(
    'four-ints', 1, list, list_to_four_ints)

dummy_plugin.register_semantic_type(IntSequence1)
dummy_plugin.register_semantic_type(IntSequence2)
dummy_plugin.register_semantic_type(Mapping)
dummy_plugin.register_semantic_type(FourInts)

dummy_plugin.register_type_to_data_layout(IntSequence1, 'int-sequence', 1)
dummy_plugin.register_type_to_data_layout(IntSequence2, 'int-sequence', 1)
dummy_plugin.register_type_to_data_layout(Mapping, 'mapping', 1)
dummy_plugin.register_type_to_data_layout(FourInts, 'four-ints', 1)

# TODO add an optional parameter to this method when they are supported
dummy_plugin.methods.register_markdown('markdown/concatenate_ints_markdown.md')
dummy_plugin.methods.register_markdown('markdown/split_ints_markdown.md')

# This method is equivalent to its markdown version above, they just have
# different IDs.
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
