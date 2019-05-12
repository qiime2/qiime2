# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2 import Artifact

from .type import IntSequence1, IntSequence2


def ints1_factory():
    return Artifact.import_data(IntSequence1, [0, 1, 2])


def ints2_factory():
    return Artifact.import_data(IntSequence1, [3, 4, 5])


def ints3_factory():
    return Artifact.import_data(IntSequence2, [6, 7, 8])


def concatenate_ints_simple(use):
    concatenate_ints = use.get_action('dummy-plugin', 'concatenate_ints')
    use.scope.add_artifact('byod', ints1_factory)
    use.scope.add_artifact('ints2', ints2_factory)
    use.scope.add_artifact('this_one_is_important', ints3_factory)

    use.comment('big data == concatenating ints')
    use.action(
        concatenate_ints,
        {
            'ints1': 'byod',
            'ints2': 'ints2',
            'ints3': 'this_one_is_important',
            'int1': 4,
            'int2': 2,
        },
        {
            'concatenated_ints': 'youre_just_a_copy_of_an_imitation',
        },
    )
    use.comment('as you can clearly see, p == np')


def concatenate_ints_complex(use):
    concatenate_ints = use.get_action('dummy-plugin', 'concatenate_ints')
    use.scope.add_artifact('byod', ints1_factory)
    use.scope.add_artifact('ints2', ints2_factory)
    use.scope.add_artifact('this_one_is_important', ints3_factory)

    use.comment('big data == concatenating ints')
    use.action(
        concatenate_ints,
        {
            'ints1': 'byod',
            'ints2': 'ints2',
            'ints3': 'this_one_is_important',
            'int1': 4,
            'int2': 2,
        },
        {
            'concatenated_ints': 'youre_just_a_copy_of_an_imitation',
        },
    )
    use.comment('as you can clearly see, p == np')
    use.action(
        concatenate_ints,
        {
            'ints1': 'youre_just_a_copy_of_an_imitation',
            'ints2': 'youre_just_a_copy_of_an_imitation',
            'ints3': 'this_one_is_important',
            'int1': 6,
            'int2': 7,
        },
        {
            'concatenated_ints': 'well_well_well_what_do_we_have_here',
        },
    )
    use.comment('fin')


def most_common_viz_typical(use):
    most_common_viz = use.get_action('dummy-plugin', 'most_common_viz')
    use.scope.add_artifact('int', ints1_factory)

    use.comment('doing things')
    use.action(
        most_common_viz,
        dict(),
        dict(),
    )
