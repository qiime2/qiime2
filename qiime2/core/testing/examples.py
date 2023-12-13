# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from qiime2 import Artifact, Metadata, ResultCollection

from .type import IntSequence1, IntSequence2, Mapping, SingleInt


def ints1_factory():
    return Artifact.import_data(IntSequence1, [0, 1, 2])


def ints2_factory():
    return Artifact.import_data(IntSequence1, [3, 4, 5])


def ints3_factory():
    return Artifact.import_data(IntSequence2, [6, 7, 8])


def artifact_collection_factory():
    return ResultCollection({'Foo': Artifact.import_data(SingleInt, 1),
                             'Bar': Artifact.import_data(SingleInt, 2)})


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
    ints_a = use.init_artifact('ints_a', ints1_factory)
    ints_b = use.init_artifact('ints_b', ints2_factory)
    ints_c = use.init_artifact('ints_c', ints3_factory)

    use.comment('This example demonstrates basic usage.')
    ints_d, = use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='concatenate_ints'),
        use.UsageInputs(ints1=ints_a, ints2=ints_b, ints3=ints_c, int1=4,
                        int2=2),
        use.UsageOutputNames(concatenated_ints='ints_d'),
    )


def concatenate_ints_complex(use):
    ints_a = use.init_artifact('ints_a', ints1_factory)
    ints_b = use.init_artifact('ints_b', ints2_factory)
    ints_c = use.init_artifact('ints_c', ints3_factory)

    use.comment('This example demonstrates chained usage (pt 1).')
    ints_d, = use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='concatenate_ints'),
        use.UsageInputs(ints1=ints_a, ints2=ints_b, ints3=ints_c, int1=4,
                        int2=2),
        use.UsageOutputNames(concatenated_ints='ints_d'),
    )

    use.comment('This example demonstrates chained usage (pt 2).')
    concatenated_ints, = use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='concatenate_ints'),
        use.UsageInputs(ints1=ints_d, ints2=ints_b, ints3=ints_c, int1=41,
                        int2=0),
        use.UsageOutputNames(concatenated_ints='concatenated_ints'),
    )


def typical_pipeline_simple(use):
    ints = use.init_artifact('ints', ints1_factory)
    mapper = use.init_artifact('mapper', mapping1_factory)

    use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='typical_pipeline'),
        use.UsageInputs(int_sequence=ints, mapping=mapper,
                        do_extra_thing=True),
        use.UsageOutputNames(out_map='out_map', left='left', right='right',
                             left_viz='left_viz', right_viz='right_viz')
    )


def typical_pipeline_complex(use):
    ints1 = use.init_artifact('ints1', ints1_factory)
    mapper1 = use.init_artifact('mapper1', mapping1_factory)

    mapper2, ints2, *_ = use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='typical_pipeline'),
        use.UsageInputs(int_sequence=ints1, mapping=mapper1,
                        do_extra_thing=True),
        use.UsageOutputNames(out_map='out_map1', left='left1', right='right1',
                             left_viz='left_viz1', right_viz='right_viz1')
    )

    _, _, right2, *_ = use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='typical_pipeline'),
        use.UsageInputs(int_sequence=ints2, mapping=mapper2,
                        do_extra_thing=False),
        use.UsageOutputNames(out_map='out_map2', left='left2', right='right2',
                             left_viz='left_viz2', right_viz='right_viz2')
    )

    right2.assert_has_line_matching(
        path='ints.txt',
        expression='1',
    )

    # test that the non-string type works
    right2.assert_output_type(semantic_type=IntSequence1)
    # test that the string type works
    mapper2.assert_output_type(semantic_type='Mapping')


def comments_only(use):
    use.comment('comment 1')
    use.comment('comment 2')


def comments_only_factory():
    def comments_only_closure(use):
        use.comment('comment 1')
        use.comment('comment 2')

    return comments_only_closure


def identity_with_metadata_simple(use):
    ints = use.init_artifact('ints', ints1_factory)
    md = use.init_metadata('md', md1_factory)

    use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='identity_with_metadata'),
        use.UsageInputs(ints=ints, metadata=md),
        use.UsageOutputNames(out='out'),
    )


def identity_with_metadata_merging(use):
    ints = use.init_artifact('ints', ints1_factory)
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
    ints = use.init_artifact('ints', ints1_factory)
    md = use.init_metadata('md', md1_factory)
    mdc = use.get_metadata_column('mdc', 'a', md)

    out, = use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='identity_with_metadata_column'),
        use.UsageInputs(ints=ints, metadata=mdc),
        use.UsageOutputNames(out='out'),
    )


def variadic_input_simple(use):
    ints_a = use.init_artifact('ints_a', ints1_factory)
    ints_b = use.init_artifact('ints_b', ints2_factory)

    single_int1 = use.init_artifact('single_int1', single_int1_factory)
    single_int2 = use.init_artifact('single_int2', single_int2_factory)

    use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='variadic_input_method'),
        use.UsageInputs(ints=[ints_a, ints_b],
                        int_set={single_int1, single_int2},
                        nums={7, 8, 9}),
        use.UsageOutputNames(output='out'),
    )


def optional_inputs(use):
    ints = use.init_artifact('ints', ints1_factory)

    output1, = use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='optional_artifacts_method'),
        use.UsageInputs(ints=ints, num1=1),
        use.UsageOutputNames(output='output1'),
    )

    output2, = use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='optional_artifacts_method'),
        use.UsageInputs(ints=ints, num1=1, num2=2),
        use.UsageOutputNames(output='output2'),
    )

    output3, = use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='optional_artifacts_method'),
        use.UsageInputs(ints=ints, num1=1, num2=None),
        use.UsageOutputNames(output='output3'),
    )

    output4, = use.action(
        use.UsageAction(plugin_id='dummy_plugin',
                        action_id='optional_artifacts_method'),
        use.UsageInputs(ints=ints, optional1=output3, num1=3, num2=4),
        use.UsageOutputNames(output='output4'),
    )


def collection_list_of_ints(use):
    ints = use.init_artifact_collection('ints', artifact_collection_factory)

    out, = use.action(
                use.UsageAction(plugin_id='dummy_plugin',
                                action_id='list_of_ints'),
                use.UsageInputs(ints=ints),
                use.UsageOutputNames(output='out'),
    )


def collection_dict_of_ints(use):
    ints = use.init_artifact_collection('ints', artifact_collection_factory)

    out, = use.action(
                use.UsageAction(plugin_id='dummy_plugin',
                                action_id='dict_of_ints'),
                use.UsageInputs(ints=ints),
                use.UsageOutputNames(output='out'),
    )

    out.assert_output_type(semantic_type='SingleInt', key='Foo')


def construct_and_access_collection(use):
    ints_a = use.init_artifact('ints_a', single_int1_factory)
    ints_b = use.init_artifact('ints_b', single_int2_factory)

    rc_in = use.construct_artifact_collection(
        'rc_in', {'a': ints_a, 'b': ints_b}
    )

    rc_out, = use.action(
        use.UsageAction(plugin_id='dummy_plugin', action_id='dict_of_ints'),
        use.UsageInputs(ints=rc_in),
        use.UsageOutputNames(output='rc_out')
    )

    ints_b_from_collection = use.get_artifact_collection_member(  # noqa: F841
        'ints_b_from_collection', rc_out, 'b'
    )
