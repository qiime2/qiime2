# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from .type import SingleInt, Mapping


def parameter_only_pipeline(ctx, int1, int2=2, metadata=None):
    identity_with_optional_metadata = ctx.get_action(
        'dummy_plugin', 'identity_with_optional_metadata')
    concatenate_ints = ctx.get_action('dummy_plugin', 'concatenate_ints')

    ints1 = ctx.make_artifact('IntSequence2', [int1, int2, 3])
    ints2, = identity_with_optional_metadata(ints1, metadata)
    ints3, = identity_with_optional_metadata(ints1, metadata)
    more_ints, = concatenate_ints(ints3, ints2, ints1, int1=int1, int2=int2)

    return ints1, more_ints


def typical_pipeline(ctx, int_sequence, mapping, do_extra_thing, add=1):
    split_ints = ctx.get_action('dummy_plugin', 'split_ints')
    most_common_viz = ctx.get_action('dummy_plugin', 'most_common_viz')

    left, right = split_ints(int_sequence)
    if do_extra_thing:
        left = ctx.make_artifact(
            'IntSequence1', [i + add for i in left.view(list)])

    val, = mapping.view(dict).values()
    # Some kind of runtime failure
    if val != '42':
        raise ValueError("Bad mapping")

    left_viz, = most_common_viz(left)
    right_viz, = most_common_viz(right)

    return mapping, left, right, left_viz, right_viz


def optional_artifact_pipeline(ctx, int_sequence, single_int=None):
    optional_artifact_method = ctx.get_action(
        'dummy_plugin', 'optional_artifacts_method')

    if single_int is None:
        # not a nested pipeline, just sharing the ctx object
        single_int = pointless_pipeline(ctx)

    num1 = single_int.view(int)
    ints, = optional_artifact_method(int_sequence, num1)
    return ints


def visualizer_only_pipeline(ctx, mapping):
    no_input_viz = ctx.get_action('dummy_plugin', 'no_input_viz')
    mapping_viz = ctx.get_action('dummy_plugin', 'mapping_viz')

    viz1, = no_input_viz()
    viz2, = mapping_viz(mapping, mapping, 'foo', 'bar')

    return viz1, viz2


def pipelines_in_pipeline(ctx, int_sequence, mapping):
    pointless_pipeline = ctx.get_action('dummy_plugin', 'pointless_pipeline')
    typical_pipeline = ctx.get_action('dummy_plugin', 'typical_pipeline')
    visualizer_only_pipeline = ctx.get_action(
        'dummy_plugin', 'visualizer_only_pipeline')

    results = []
    results += pointless_pipeline()
    typical_results = typical_pipeline(int_sequence, mapping, True)
    results += typical_results
    results += visualizer_only_pipeline(typical_results[0])

    return tuple(results)


def pointless_pipeline(ctx):
    # Use a real type expression instead of a string.
    return ctx.make_artifact(SingleInt, 4)


def failing_pipeline(ctx, int_sequence, break_from='arity'):
    merge_mappings = ctx.get_action('dummy_plugin', 'merge_mappings')

    list_ = int_sequence.view(list)
    if list_:
        integer = list_[0]
    else:
        integer = 0

    # Made here so that we can make sure it gets cleaned up
    wrong_output = ctx.make_artifact(SingleInt, integer)

    if break_from == 'arity':
        return int_sequence, int_sequence, int_sequence
    elif break_from == 'return-view':
        return None
    elif break_from == 'type':
        return wrong_output
    elif break_from == 'method':
        a = ctx.make_artifact(Mapping, {'foo': 'a'})
        b = ctx.make_artifact(Mapping, {'foo': 'b'})
        # has the same key
        merge_mappings(a, b)
    elif break_from == 'no-plugin':
        ctx.get_action('not%a$plugin', 'foo')
    elif break_from == 'no-action':
        ctx.get_action('dummy_plugin', 'not%a$method')
    else:
        raise ValueError('this never works')
