# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import inspect
import frontmatter

from qiime.sdk.type import Type


class Signature:
    def __init__(self, name, inputs, outputs):
        """
           Parameters
           ----------
           name : str, name for this workflow
           inputs : dict, names mapped to Types
           outputs : list of tuples, names mapped to artifact Types
        """
        self.name = name
        self.input_artifacts = {}
        self.input_parameters = {}
        for input_name, input_type in inputs.items():
            if Type.Primitive.is_member(input_type):
                self.input_parameters[input_name] = input_type
            elif Type.Artifact.is_member(input_type):
                self.input_artifacts[input_name] = input_type
            else:
                raise TypeError("Unrecognized input type: %r" % input_type)
        self.output_artifacts = outputs

    def __call__(self, artifact, parameters):
        # TODO implement me!
        return self.output_artifacts


class WorkflowTemplate:
    def __init__(self, signature, template):
        """
           Parameters
           ----------
           signature : Signature
           template : str, markdown text in ipymd style
        """
        self.signature = signature
        self.template = template

    def create_job(self, setup_lines, teardown_lines):
        setup_str = '\n'.join(['>>> %s' % line for line in setup_lines])
        teardown_str = '\n'.join(['>>> %s' % line for line in teardown_lines])
        setup_str = _markdown_code_cell_template.format(content=setup_str)
        teardown_str = _markdown_code_cell_template.format(
                content=teardown_str)
        return "\n\n".join([setup_str, self.template, teardown_str])

    @classmethod
    def from_markdown(cls, markdown):
        """
           Parameters
           ----------
           markdown : filepath
        """
        with open(markdown) as fh:
            metadata, template = frontmatter.parse(fh.read())
        inputs = metadata['inputs']
        outputs = metadata['outputs'].items()
        name = metadata['name']
        signature = Signature(name, inputs, outputs)
        return cls(signature, template)

    @classmethod
    def from_function(cls, function, inputs, outputs, name, doc):
        # TODO: the import paths will not necessarily be a public
        # import, but might be a full (private) import. fix that. 
        import_path = inspect.getmodule(function).__name__
        function_name = function.__name__
        parameters = ', '.join(['%s=%s' % (k, k) for k in inputs])
        results = ', '.join(outputs.keys())
        template = _markdown_template.format(
            doc=doc, import_path=import_path, function_name=function_name,
            parameters=parameters, results=results)
        signature = Signature(name, inputs, outputs)
        return cls(signature, template)

_markdown_template = """{doc}

```python
>>> from {import_path} import {function_name}
>>> {results} = {function_name}({parameters})
```"""

_markdown_code_cell_template = """```python
{content}
```"""
