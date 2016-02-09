# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import inspect


class Workflow(object):

    def __init__(self, inputs, outputs, workflow_template, doc=""):
        """
           Parameters
           ----------
           inputs : dict, names mapped to ArtifactTypes
           outputs : list of tuples, names mapped to ArtifactTypes
           workflow_template : str, markdown text in ipymd style
           doc : str
        """
        self.inputs = inputs
        self.outputs = outputs
        self.workflow_template = workflow_template
        self.doc = doc

    @classmethod
    def from_markdown(cls, markdown):
        """
           Parameters
           ----------
           markdown : file handle
        """
        inputs, outputs, workflow_template, doc = [], [], "", ""
        return cls(inputs, outputs, workflow_template, doc)

    @classmethod
    def from_function(cls, function, inputs, outputs, doc):
        import_path = inspect.getmodule(function).__name__
        function_name = function.__name__
        parameters = ', '.join(['%s=%s' % (k, k) for k in inputs])
        results = ', '.join([output[0] for output in outputs])
        workflow_template = _markdown_template.format(
            doc=doc, import_path=import_path, function_name=function_name,
            parameters=parameters, results=results)
        return cls(inputs, outputs, workflow_template, doc)

_markdown_template = """
{doc}

```python
>>> from {import_path} import {function_name}
>>> {results} = {function_name}({parameters})
```
"""
