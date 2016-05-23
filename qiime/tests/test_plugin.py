# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import os
import pkg_resources
import tempfile
import unittest
import unittest.mock

import frontmatter

from qiime.plugin import Plugin, Int, ArchiveFormat
from qiime.sdk import Workflow, Signature


def dummy_function() -> int:
    return 42


def other_dummy_function() -> int:
    return 42


class TestPlugin(unittest.TestCase):
    def setUp(self):
        # TODO standardize temporary directories created by QIIME
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-temp-')

        self.markdown_fp = os.path.join(self.test_dir.name,
                                        'dummy_markdown_workflow.md')
        with open(self.markdown_fp, 'w') as markdown_fh:
            markdown_fh.write(markdown_template)

        self.plugin = Plugin(
            name='dummy-plugin',
            version='0.4.2',
            website='www.dummy-plugin-hub.com',
            package='dummy_plugin'
        )

    def tearDown(self):
        self.test_dir.cleanup()

    def test_constructor(self):
        self.assertEqual(self.plugin.name, 'dummy-plugin')
        self.assertEqual(self.plugin.version, '0.4.2')
        self.assertEqual(self.plugin.website, 'www.dummy-plugin-hub.com')
        self.assertEqual(self.plugin.package, 'dummy_plugin')
        self.assertEqual(self.plugin.workflows, {})
        self.assertEqual(self.plugin.archive_formats, {})

    def test_register_archive_format(self):
        self.assertEqual(self.plugin.archive_formats, {})

        def is_valid(data_reader):
            return True
        self.plugin.register_archive_format('my-format', '1.0.0', is_valid)
        exp_format = ArchiveFormat('my-format', '1.0.0', is_valid)
        self.assertEqual(self.plugin.archive_formats,
                         {('my-format', '1.0.0'): exp_format})
        # TODO: test execution of archive format reader and writer

    def test_register_function(self):
        self.assertEqual(self.plugin.workflows, {})

        self.plugin.register_function(
            name='Dummy function',
            function=dummy_function,
            inputs={},
            outputs=[('answer', Int)],
            doc='Computes the answer to life, the universe, and everything'
        )

        self.plugin.register_function(
            name='Dummy function',
            function=other_dummy_function,
            inputs={},
            outputs=[('answer', Int)],
            doc='Computes the answer to life, the universe, and everything'
        )

        workflows = {
            'dummy_function':
                Workflow(
                    signature=Signature(
                        name='Dummy function',
                        inputs={},
                        outputs=collections.OrderedDict([('answer', Int)])),
                    template=expected_dummy_function_template,
                    id_='dummy_function'
                ),
            'other_dummy_function':
                Workflow(
                    signature=Signature(
                        name='Dummy function',
                        inputs={},
                        outputs=collections.OrderedDict([('answer', Int)])),
                    template=expected_other_dummy_function_template,
                    id_='other_dummy_function'
                )
        }

        self.assertEqual(self.plugin.workflows, workflows)

    def test_register_workflow(self):
        self.assertEqual(self.plugin.workflows, {})

        with unittest.mock.patch.object(pkg_resources, 'resource_filename',
                                        return_value=self.markdown_fp):
            self.plugin.register_workflow(self.markdown_fp)

        workflows = {
            'dummy_markdown_workflow':
                Workflow(
                    signature=Signature(
                        name='Dummy markdown workflow',
                        inputs={'param1': Int, 'param2': Int},
                        outputs=collections.OrderedDict([('the_sum', Int)])),
                    template=frontmatter.parse(markdown_template)[1],
                    id_='dummy_markdown_workflow'
                )
        }

        self.assertEqual(self.plugin.workflows, workflows)

    def test_register_function_and_workflow(self):
        self.assertEqual(self.plugin.workflows, {})

        self.plugin.register_function(
            name='Dummy function',
            function=dummy_function,
            inputs={},
            outputs=[('answer', Int)],
            doc='Computes the answer to life, the universe, and everything'
        )

        with unittest.mock.patch.object(pkg_resources, 'resource_filename',
                                        return_value=self.markdown_fp):
            self.plugin.register_workflow(self.markdown_fp)

        workflows = {
            'dummy_function':
                Workflow(
                    signature=Signature(
                        name='Dummy function',
                        inputs={},
                        outputs=collections.OrderedDict([('answer', Int)])),
                    template=expected_dummy_function_template,
                    id_='dummy_function'
                ),
            'dummy_markdown_workflow':
                Workflow(
                    signature=Signature(
                        name='Dummy markdown workflow',
                        inputs={'param1': Int, 'param2': Int},
                        outputs=collections.OrderedDict([('the_sum', Int)])),
                    template=frontmatter.parse(markdown_template)[1],
                    id_='dummy_markdown_workflow'
                )
        }

        self.assertEqual(self.plugin.workflows, workflows)


markdown_template = """---
name: Dummy markdown workflow
type-imports:
    - qiime.plugin:Int
inputs:
    param1:
        - Int
        - int
    param2:
        - Int
        - int
outputs:
    - the_sum:
        - Int
        - int
---
## Sum some integers and return their summation

This workflow sums integers in the following way:

```python
>>> the_sum = param1 + param2
```
"""

expected_dummy_function_template = """Computes the answer to life, the universe, and everything

```python
>>> from qiime.tests.test_plugin import dummy_function
>>> answer = dummy_function()
```
"""

expected_other_dummy_function_template = """Computes the answer to life, the universe, and everything

```python
>>> from qiime.tests.test_plugin import other_dummy_function
>>> answer = other_dummy_function()
```
"""


if __name__ == '__main__':
    unittest.main()
