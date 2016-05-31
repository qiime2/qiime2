# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import os
import tempfile
import unittest
import uuid

import frontmatter

import qiime.plugin
import qiime.core.testing
from qiime.sdk import Artifact, Signature, Workflow


class TestWorkflow(unittest.TestCase):
    def setUp(self):
        # TODO standardize temporary directories created by QIIME
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-temp-')

        self.markdown_fp = os.path.join(self.test_dir.name,
                                        'dummy_markdown_workflow.md')
        with open(self.markdown_fp, 'w') as markdown_fh:
            markdown_fh.write(markdown_template)

    def tearDown(self):
        self.test_dir.cleanup()

    def test_signature(self):
        workflow = Workflow.from_markdown(self.markdown_fp)

        signature = Signature(
            name='Dummy markdown workflow',
            inputs={
                'input1': (qiime.core.testing.TestType, list),
                'input2': (qiime.core.testing.TestType, list),
                'param1': (qiime.plugin.Int, int),
                'param2': (qiime.plugin.Int, int),
            },
            outputs=collections.OrderedDict([
                ('concatenated_inputs', (qiime.core.testing.TestType, list))
            ])
        )

        self.assertEqual(workflow.signature, signature)

    def test_template(self):
        workflow = Workflow.from_markdown(self.markdown_fp)

        self.assertEqual(workflow.template,
                         frontmatter.parse(markdown_template)[1])

    def test_id(self):
        workflow = Workflow.from_markdown(self.markdown_fp)

        self.assertEqual(workflow.id, 'dummy_markdown_workflow')

    def test_name(self):
        workflow = Workflow.from_markdown(self.markdown_fp)

        self.assertEqual(workflow.name, 'Dummy markdown workflow')

    def test_reference(self):
        workflow = Workflow.from_markdown(self.markdown_fp)

        self.assertTrue(workflow.reference.startswith(workflow.id))

    def test_from_markdown(self):
        workflow = Workflow.from_markdown(self.markdown_fp)

        expected = Workflow(
            signature=Signature(
                name='Dummy markdown workflow',
                inputs={
                    'input1': (qiime.core.testing.TestType, list),
                    'input2': (qiime.core.testing.TestType, list),
                    'param1': (qiime.plugin.Int, int),
                    'param2': (qiime.plugin.Int, int),
                },
                outputs=collections.OrderedDict([
                    ('concatenated_inputs',
                     (qiime.core.testing.TestType, list))
                ])
            ),
            template=frontmatter.parse(markdown_template)[1],
            id_='dummy_markdown_workflow'
        )

        self.assertEqual(workflow, expected)

    def test_from_function(self):
        def dummy_function(input1: list, input2: list,
                           param1: int, param2: int) -> list:
            return input1 + input2 + [param1] + [param2]

        workflow = Workflow.from_function(
            dummy_function,
            inputs={
                'input1': qiime.core.testing.TestType,
                'input2': qiime.core.testing.TestType,
                'param1': qiime.plugin.Int,
                'param2': qiime.plugin.Int
            },
            outputs=collections.OrderedDict([
                ('concatenated_inputs', qiime.core.testing.TestType)
            ]),
            name='Concatenate things',
            doc="Let's concatenate some things!"
        )

        expected = Workflow(
            signature=Signature(
                name='Concatenate things',
                inputs={
                    'input1': (qiime.core.testing.TestType, list),
                    'input2': (qiime.core.testing.TestType, list),
                    'param1': (qiime.plugin.Int, int),
                    'param2': (qiime.plugin.Int, int),
                },
                outputs=collections.OrderedDict([
                    ('concatenated_inputs',
                     (qiime.core.testing.TestType, list))
                ])
            ),
            template=expected_template,
            id_='dummy_function'
        )

        self.assertEqual(workflow, expected)

    def test_to_script(self):
        # Python script should have non-code lines commented out.
        template_lines = [
            "# ### Join lists together",
            "concatenated_inputs = input1 + input2",
            "concatenated_inputs.append(param2)"
        ]

        self._test_to_script_or_to_markdown(Workflow.to_script, template_lines)

    def test_to_markdown(self):
        # Should not comment markdown text.
        template_lines = [
            "### Join lists together",
            "concatenated_inputs = input1 + input2",
            "concatenated_inputs.append(param2)"
        ]

        self._test_to_script_or_to_markdown(Workflow.to_markdown,
                                            template_lines)

    def _test_to_script_or_to_markdown(self, to_method, template_lines):
        # These methods are so similar that it makes sense to have a helper
        # that can test either one instead of duplicating a bunch of code.
        workflow = Workflow.from_markdown(self.markdown_fp)

        artifact_fp1 = os.path.join(self.test_dir.name, 'artifact1.qzf')
        artifact = Artifact._from_view(
            [-1, 42, 0, 43, 43], qiime.core.testing.TestType, None)
        artifact.save(artifact_fp1)

        artifact_fp2 = os.path.join(self.test_dir.name, 'artifact2.qzf')
        artifact = Artifact._from_view(
            [1, 2, 100], qiime.core.testing.TestType, None)
        artifact.save(artifact_fp2)

        artifact_fp3 = os.path.join(self.test_dir.name, 'artifact3.qzf')

        job = to_method(
            workflow,
            input_artifact_filepaths={
                'input1': artifact_fp1,
                'input2': artifact_fp2
            },
            parameter_references={
                'param1': 99,
                'param2': -999,
            },
            output_artifact_filepaths={
                'concatenated_inputs': artifact_fp3
            }
        )

        provenance_lines = [
            "provenance = Provenance(",
            "parameters={",
            "'param2': -999"
        ]

        setup_lines = [
            "input1 = Artifact.load(%r).view(list)" % artifact_fp1,
            "input2 = Artifact.load(%r).view(list)" % artifact_fp2,
            "param1 = 99",
            "param2 = -999"
        ]

        teardown_lines = [
            "artifact = Artifact._from_view(concatenated_inputs, TestType, "
            "provenance)",
            "artifact.save(%r)" % artifact_fp3
        ]

        for expected_lines in (provenance_lines, setup_lines, template_lines,
                               teardown_lines):
            for expected_line in expected_lines:
                self.assertIn(expected_line, job.code)

        self.assertIsInstance(job.uuid, uuid.UUID)
        self.assertEqual(job.uuid.version, 4)

        self.assertEqual(
            job.input_artifact_filepaths,
            {'input1': artifact_fp1, 'input2': artifact_fp2})

        self.assertEqual(
            job.parameter_references,
            {'param1': 99, 'param2': -999})

        self.assertEqual(
            job.output_artifact_filepaths,
            {'concatenated_inputs': artifact_fp3})


markdown_template = """---
name: Dummy markdown workflow
type-imports:
    - qiime.core.testing:TestType
    - qiime.plugin:Int
inputs:
    input1:
        - TestType
        - list
    input2:
        - TestType
        - list
    param1:
        - Int
        - int
    param2:
        - Int
        - int
outputs:
    - concatenated_inputs:
        - TestType
        - list
---
## Concatenate some things together

This workflow concatenates inputs together in the following steps.

### Join lists together

```python
>>> concatenated_inputs = input1 + input2
```

### Append integers

```python
>>> concatenated_inputs.append(param1)
>>> concatenated_inputs.append(param2)
```
"""

expected_template = """Let's concatenate some things!

```python
>>> from qiime.sdk.tests.test_workflow import dummy_function
>>> concatenated_inputs = dummy_function(input1=input1, input2=input2, \
param1=param1, param2=param2)
```
"""


if __name__ == '__main__':
    unittest.main()
