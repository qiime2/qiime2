# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os.path

import qiime
import qiime.plugin

# TODO centralize "dummy" functions and workflows that are currently duplicated
# across unit tests.
plugin = qiime.plugin.Plugin(
    name='test-plugin',
    version=qiime.__version__,
    website='https://github.com/qiime2/qiime2',
    package='qiime'
)


def validator(data_dir):
    raise NotImplementedError()


plugin.register_archive_format('example-archive-format', 1, validator)


def example_archive_format_to_list(data_dir):
    with open(os.path.join(data_dir, 'data.txt'), 'r') as fh:
        model = []
        for line in fh:
            model.append(int(line.rstrip()))
        return model


plugin.register_archive_format_reader('example-archive-format', 1, list,
                                      example_archive_format_to_list)


def list_to_example_archive_format(view, data_dir):
    with open(os.path.join(data_dir, 'data.txt'), 'w') as fh:
        for num in view:
            fh.write('%d\n' % num)


plugin.register_archive_format_writer('example-archive-format', 1, list,
                                      list_to_example_archive_format)

TestType = qiime.plugin.SemanticType('TestType')

plugin.register_semantic_type(TestType)

plugin.register_type_to_archive_format(TestType, 'example-archive-format', 1)
