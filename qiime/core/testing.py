# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os.path

import nose

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

plugin.register_archive_format('test-archive-format', 1, validator)

# Not using `qiime.plugin.get_archive_format` to avoid circular ImportError.
test_archive_format = plugin.archive_formats['test-archive-format', 1]


@test_archive_format.reader(list)
@nose.tools.nottest
def test_archive_format_to_list(data_dir):
    with open(os.path.join(data_dir, 'data.txt'), 'r') as fh:
        model = []
        for line in fh:
            model.append(int(line.rstrip()))
        return model


@test_archive_format.writer(list)
@nose.tools.nottest
def list_to_test_archive_format(view, data_dir):
    with open(os.path.join(data_dir, 'data.txt'), 'w') as fh:
        for num in view:
            fh.write('%d\n' % num)


class TestType(qiime.plugin.Type, variant_of=qiime.plugin.Type.Artifact):
    pass
