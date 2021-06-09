# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from qiime2 import Artifact, Metadata

from .type import IntSequence1, IntSequence2, Mapping, SingleInt


def ints1_factory():
    return Artifact.import_data(IntSequence1, [0, 1, 2])


def ints2_factory():
    return Artifact.import_data(IntSequence1, [3, 4, 5])


def ints3_factory():
    return Artifact.import_data(IntSequence2, [6, 7, 8])


def mapping1_factory():
    return Artifact.import_data(Mapping, {'a': 42})


def md1_factory():
    return Metadata(pd.DataFrame({'a': ['1', '2', '3']},
                                 index=pd.Index(['0', '1', '2'],
                                                name='id')))


def md2_factory():
    return Metadata(pd.DataFrame({'b': ['4', '5', '6']},
                                 index=pd.Index(['0', '1', '2'],
                                                name='id')))


def single_int1_factory():
    return Artifact.import_data(SingleInt, 10)


def single_int2_factory():
    return Artifact.import_data(SingleInt, 11)


def concatenate_ints_simple(use):
    ints_a = use.init_data('ints_a', ints1_factory)
    ints_b = use.init_data('ints_b', ints2_factory)
    ints_c = use.init_data('ints_c', ints3_factory)

    use.comment('This example demonstrates basic usage.')
    use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='concatenate_ints'),
        use.UsageInputs(ints1=ints_a, ints2=ints_b, ints3=ints_c, int1=4,
                        int2=2),
        use.UsageOutputNames(concatenated_ints='ints_d'),
    )


def concatenate_ints_complex(use):
    ints_a = use.init_data('ints_a', ints1_factory)
    ints_b = use.init_data('ints_b', ints2_factory)
    ints_c = use.init_data('ints_c', ints3_factory)

    use.comment('This example demonstrates chained usage (pt 1).')
    use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='concatenate_ints'),
        use.UsageInputs(ints1=ints_a, ints2=ints_b, ints3=ints_c, int1=4,
                        int2=2),
        use.UsageOutputNames(concatenated_ints='ints_d'),
    )

    ints_d = use.get_result('ints_d')
    use.comment('This example demonstrates chained usage (pt 2).')
    use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='concatenate_ints'),
        use.UsageInputs(ints1=ints_d, ints2=ints_b, ints3=ints_c, int1=41,
                        int2=0),
        use.UsageOutputNames(concatenated_ints='concatenated_ints'),
    )


def typical_pipeline_simple(use):
    ints = use.init_data('ints', ints1_factory)
    mapper = use.init_data('mapper', mapping1_factory)

    use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='typical_pipeline'),
        use.UsageInputs(int_sequence=ints, mapping=mapper,
                        do_extra_thing=True),
        use.UsageOutputNames(out_map='out_map', left='left', right='right',
                             left_viz='left_viz', right_viz='right_viz')
    )


def typical_pipeline_complex(use):
    ints1 = use.init_data('ints1', ints1_factory)
    mapper1 = use.init_data('mapper1', mapping1_factory)

    use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='typical_pipeline'),
        use.UsageInputs(int_sequence=ints1, mapping=mapper1,
                        do_extra_thing=True),
        use.UsageOutputNames(out_map='out_map1', left='left1', right='right1',
                             left_viz='left_viz1', right_viz='right_viz1')
    )

    ints2 = use.get_result('left1')
    mapper2 = use.get_result('out_map1')

    use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='typical_pipeline'),
        use.UsageInputs(int_sequence=ints2, mapping=mapper2,
                        do_extra_thing=False),
        use.UsageOutputNames(out_map='out_map2', left='left2', right='right2',
                             left_viz='left_viz2', right_viz='right_viz2')
    )

    right2 = use.get_result('right2')
    right2.assert_has_line_matching(
        label='a nice label about this assertion',
        path='ints.txt',
        expression='1',
    )


def comments_only(use):
    use.comment('comment 1')
    use.comment('comment 2')


def comments_only_factory():
    def comments_only_closure(use):
        use.comment('comment 1')
        use.comment('comment 2')

    return comments_only_closure


def identity_with_metadata_simple(use):
    ints = use.init_data('ints', ints1_factory)
    md = use.init_metadata('md', md1_factory)

    use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='identity_with_metadata'),
        use.UsageInputs(ints=ints, metadata=md),
        use.UsageOutputNames(out='out'),
    )


def identity_with_metadata_merging(use):
    ints = use.init_data('ints', ints1_factory)
    md1 = use.init_metadata('md1', md1_factory)
    md2 = use.init_metadata('md2', md2_factory)

    md3 = use.merge_metadata('md3', md1, md2)

    use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='identity_with_metadata'),
        use.UsageInputs(ints=ints, metadata=md3),
        use.UsageOutputNames(out='out'),
    )


def identity_with_metadata_column_get_mdc(use):
    ints = use.init_data('ints', ints1_factory)
    md = use.init_metadata('md', md1_factory)
    mdc = use.get_metadata_column('a', md)

    use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='identity_with_metadata_column'),
        use.UsageInputs(ints=ints, metadata=mdc),
        use.UsageOutputNames(out='out'),
    )


def variadic_input_simple(use):
    ints_a = use.init_data('ints_a', ints1_factory)
    ints_b = use.init_data('ints_b', ints2_factory)
    ints = use.init_data_collection('ints', list, ints_a, ints_b)

    single_int1 = use.init_data('single_int1', single_int1_factory)
    single_int2 = use.init_data('single_int2', single_int2_factory)
    int_set = use.init_data_collection('int_set', set, single_int1,
                                       single_int2)

    use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='variadic_input_method'),
        use.UsageInputs(ints=ints, int_set=int_set, nums={7, 8, 9}),
        use.UsageOutputNames(output='out'),
    )


def optional_inputs(use):
    ints_a = use.init_data('ints', ints1_factory)

    use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='optional_artifacts_method'),
        use.UsageInputs(ints=ints_a, num1=1),
        use.UsageOutputNames(output='output'),
    )

    use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='optional_artifacts_method'),
        use.UsageInputs(ints=ints_a, num1=1, num2=2),
        use.UsageOutputNames(output='output'),
    )

    use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='optional_artifacts_method'),
        use.UsageInputs(ints=ints_a, num1=1, num2=None),
        use.UsageOutputNames(output='ints_b'),
    )

    ints_b = use.get_result('ints_b')

    use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='optional_artifacts_method'),
        use.UsageInputs(ints=ints_a, optional1=ints_b, num1=3, num2=4),
        use.UsageOutputNames(output='output'),
    )
