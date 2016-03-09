# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import importlib
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
                self.input_parameters[input_name] = input_type()
            elif Type.Artifact.is_member(input_type):
                self.input_artifacts[input_name] = input_type()
            else:
                raise TypeError("Unrecognized input type: %r" % input_type)
        self.output_artifacts = outputs

    def __call__(self, artifacts, parameters):
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

    @property
    def name(self):
        return self.signature.name

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

        type_imports = metadata['type-imports']
        input_types = {}
        for name, type_expr in metadata['inputs'].items():
            input_types[name] = cls._parse_type(type_imports, type_expr)

        output_types = collections.OrderedDict()
        for output in metadata['outputs']:
            name, type_expr = next(output.items())
            output_types[name] = cls._parse_type(type_imports, type_expr)

        name = metadata['name']
        signature = Signature(name, input_types, output_types)
        return cls(signature, template)

    # TODO this is duplicated from Artifact._parse_type. Refactor!
    @classmethod
    def _parse_type(cls, imports, type_exp):
        type_exp = type_exp.split('\n')
        if len(type_exp) != 1:
            raise TypeError("Multiple lines in type expression of"
                            " artifact. Will not load to avoid arbitrary"
                            " code execution.")
        type_exp, = type_exp

        if ';' in type_exp:
            raise TypeError("Invalid type expression in artifact. Will not"
                            " load to avoid arbitrary code execution.")

        locals_ = {}
        for import_ in imports:
            path, class_ = import_.split(":")
            try:
                module = importlib.import_module(path)
            except ImportError:
                raise ImportError("The plugin which defines: %r is not"
                                  " installed." % path)
            class_ = getattr(module, class_)
            if not issubclass(class_, Type):
                raise TypeError("Non-Type artifact. Will not load to avoid"
                                " arbitrary code execution.")
            if class_.__name__ in locals_:
                raise TypeError("Duplicate type name (%r) in expression."
                                % class_.__name__)
            locals_[class_.__name__] = class_
        type_ = eval(type_exp, {'__builtins__': {}}, locals_)
        return type_

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
