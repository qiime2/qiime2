# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import uuid
import collections
import importlib
import inspect
import os.path
import textwrap

import frontmatter

import qiime.sdk
import qiime.sdk.job
import qiime.core.type
import qiime.plugin


# TODO `name` is being used as a human-friendly name. I don't think this should
# be stored in Signature, it makes more sense in Workflow.
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
            if qiime.core.type.BaseType.Primitive.is_member(input_type):
                self.input_parameters[input_name] = input_type()
            elif qiime.plugin.Type.Artifact.is_member(input_type):
                self.input_artifacts[input_name] = input_type()
            else:
                raise TypeError("Unrecognized input type: %r" % input_type)
        self.output_artifacts = outputs

    def __call__(self, artifacts, parameters):
        # TODO implement me!
        return self.output_artifacts


class Workflow:
    @classmethod
    def from_markdown(cls, markdown):
        """
           Parameters
           ----------
           markdown : filepath
        """
        id_ = os.path.splitext(os.path.split(markdown)[1])[0]
        # TODO: verify that `id_` is machine-friendly
        with open(markdown) as fh:
            metadata, template = frontmatter.parse(fh.read())

        type_imports = metadata['type-imports']
        input_types = {}
        for name, type_expr in metadata['inputs'].items():
            input_types[name] = cls._parse_type(type_imports, type_expr)

        output_types = collections.OrderedDict()
        for output in metadata['outputs']:
            # TODO validate each nested dict has exactly one item
            name, type_expr = list(output.items())[0]
            output_types[name] = cls._parse_type(type_imports, type_expr)

        name = metadata['name']
        signature = Signature(name, input_types, output_types)
        return cls(signature, template, id_)

    # TODO this is mostly duplicated from Artifact._parse_type. Refactor!
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
            if not issubclass(class_, qiime.core.type.BaseType):
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
        return cls(signature, template, function.__name__)

    def __init__(self, signature, template, id_):
        """
           Parameters
           ----------
           signature : Signature
           template : str, markdown text in ipymd style
        """
        self.signature = signature
        self.template = template
        self.id = id_

    @property
    def name(self):
        return self.signature.name

    @property
    def reference(self):
        # TODO: fill this in with more information for tracking the source of
        # the template. This will come from the plugin system, through a
        # reverse lookup on the workflow name (or better identifier)
        #
        # Should we also/instead store the workflow template markdown, maybe
        # under a `provenance` directory?
        return ("%s. Details on plugin, version, and website will also be "
                "included, see https://github.com/biocore/qiime2/issues/26 "
                % self.id)

    def to_markdown(self, input_artifact_filepaths, parameter_references,
                    output_artifact_filepaths):
        artifacts = self._load_artifacts(input_artifact_filepaths)
        output_artifact_types = self.signature(artifacts, parameter_references)

        parameters = {}
        for name, ref in parameter_references.items():
            parameters[name] = \
                    self.signature.input_parameters[name].from_string(ref)

        job_uuid = uuid.uuid4()
        provenance_lines = self._provenance_lines(job_uuid, artifacts,
                                                  parameters)

        setup_lines = []
        setup_lines.append('from qiime.sdk.artifact import Artifact')

        for name, filepath in input_artifact_filepaths.items():
            setup_lines.append('%s = Artifact(%r).data' % (name, filepath))

        for name, value in parameters.items():
            setup_lines.append('%s = %r' % (name, value))

        teardown_lines = []
        teardown_lines.extend(
            self._teardown_import_lines(output_artifact_types))

        # TODO make sure output order is respected
        for name, output_filepath in output_artifact_filepaths.items():
            output_type = output_artifact_types[name]
            teardown_lines.append(
                'Artifact.save(%s, %r, provenance, %r)' % (name, output_type,
                                                           output_filepath))

        markdown = self._format_markdown(
            provenance_lines, setup_lines, teardown_lines)
        return qiime.sdk.job.Job(markdown, job_uuid, input_artifact_filepaths,
                                 parameter_references,
                                 output_artifact_filepaths)

    def _load_artifacts(self, artifact_filepaths):
        artifacts = {}
        for name, filepath in artifact_filepaths.items():
            artifacts[name] = qiime.sdk.Artifact(filepath)
        return artifacts

    def _provenance_lines(self, job_uuid, artifacts, parameters):
        # TODO find a better way to format provenance, this is painful...
        provenance_lines = [
            "from qiime.sdk.provenance import Provenance",
            "provenance = Provenance(",
            "    job_uuid=%r," % str(job_uuid),
            "    artifact_uuids={"
        ]

        for name, artifact in artifacts.items():
            provenance_lines.append(
                "        %r: %r," % (name, str(artifact.uuid)))
        provenance_lines.append("    },")

        provenance_lines.append("    parameters={")
        for name, value in parameters.items():
            provenance_lines.append("        %r: %r," % (name, value))
        provenance_lines.append("    },")

        provenance_lines.append("    workflow_reference=(")
        indent = ' ' * 8
        workflow_reference_lines = textwrap.wrap(
            self.reference, width=79 - len(indent), drop_whitespace=False)
        for line in workflow_reference_lines:
            provenance_lines.append('%s%r' % (indent, line))
        provenance_lines.append("    )")

        provenance_lines.append(")")
        return provenance_lines

    def _teardown_import_lines(self, output_artifact_types):
        # TODO collapse imports with common prefix
        imports = set()
        for output_type in output_artifact_types.values():
            imports = imports.union(output_type().get_imports())
        return ['from %s import %s' % (path, name) for name, path in imports]

    def _format_markdown(self, provenance_lines, setup_lines, teardown_lines):
        provenance_str = self._format_markdown_code_cell(provenance_lines)
        setup_str = self._format_markdown_code_cell(setup_lines)
        teardown_str = self._format_markdown_code_cell(teardown_lines)
        return "\n\n".join(
            [provenance_str, setup_str, self.template, teardown_str])

    def _format_markdown_code_cell(self, code_lines):
        code_str = '\n'.join(['>>> %s' % line for line in code_lines])
        return _markdown_code_cell_template.format(content=code_str)

    # TODO implement to return Python script that can be executed
    def to_script(self):
        pass

    # TODO implement to return Python function that can be executed
    def to_function(self):
        pass


_markdown_template = """{doc}

```python
>>> from {import_path} import {function_name}
>>> {results} = {function_name}({parameters})
```"""

_markdown_code_cell_template = """```python
{content}
```"""
