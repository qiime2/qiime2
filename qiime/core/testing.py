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
import qiime.core.type

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


plugin.register_data_layout('example-data-layout', 1, validator)


def example_data_layout_to_list(data_dir):
    with open(os.path.join(data_dir, 'data.txt'), 'r') as fh:
        model = []
        for line in fh:
            model.append(int(line.rstrip()))
        return model


plugin.register_data_layout_reader('example-data-layout', 1, list,
                                   example_data_layout_to_list)


def list_to_example_data_layout(view, data_dir):
    with open(os.path.join(data_dir, 'data.txt'), 'w') as fh:
        for num in view:
            fh.write('%d\n' % num)


plugin.register_data_layout_writer('example-data-layout', 1, list,
                                   list_to_example_data_layout)

TestType = qiime.plugin.SemanticType('TestType')

plugin.register_semantic_type(TestType)

plugin.register_type_to_data_layout(TestType, 'example-data-layout', 1)


def visualizer1(output_dir: str, ints: list) -> None:
    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write('<html><body>\n')
        for i in ints:
            fh.write('%d<br>\n' % i)
        fh.write('</body></html>')

plugin.register_visualization(
    function=visualizer1,
    inputs={'ints': TestType},
    parameters={},
    name='Visualizer #1',
    doc="Let's write some integers to an html file!"
)
