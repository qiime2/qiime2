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
import ipymd

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
        name : str
            Human-readable name for this workflow.
        inputs : dict
            Parameter name to semantic type.
        outputs : collections.OrderedDict
            Named output to semantic type.


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

    def __eq__(self, other):
        return (
            self.name == other.name and
            self.input_artifacts == other.input_artifacts and
            self.input_parameters == other.input_parameters and
            self.output_artifacts == other.output_artifacts
        )


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
        input_views = {}
        for name, type_expr in metadata['inputs'].items():
            semantic_type, view_type = cls._parse_type(type_imports, type_expr)
            input_types[name] = semantic_type
            input_views[name] = view_type

        output_types = collections.OrderedDict()
        output_views = collections.OrderedDict()
        for output in metadata['outputs']:
            # TODO validate each nested dict has exactly one item
            name, type_expr = list(output.items())[0]
            semantic_type, view_type = cls._parse_type(type_imports, type_expr)
            output_types[name] = semantic_type
            output_views[name] = view_type

        name = metadata['name']
        signature = Signature(name, input_types, output_types)
        return cls(signature, template, id_)

    # TODO this is mostly duplicated from Artifact._parse_type. Refactor!
    @classmethod
    def _parse_type(cls, imports, type_exp):
        # Split the type expression into its components: the semantic_type_exp
        # and the view_type. Note that this differs from the type definitions
        # in Artifact._parse_type, as those won't have view types.
        semantic_type_exp, view_type = type_exp

        semantic_type_exp = semantic_type_exp.split('\n')
        if len(semantic_type_exp) != 1:
            raise TypeError("Multiple lines in type expression of"
                            " artifact. Will not load to avoid arbitrary"
                            " code execution.")
        semantic_type_exp, = semantic_type_exp

        if ';' in semantic_type_exp:
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
                raise TypeError("Non-Type artifact (%r). Will not load to"
                                " avoid arbitrary code execution."
                                % class_.__name__)
            if class_.__name__ in locals_:
                raise TypeError("Duplicate type name (%r) in expression."
                                % class_.__name__)
            locals_[class_.__name__] = class_
        type_ = eval(semantic_type_exp, {'__builtins__': {}}, locals_)
        return type_, view_type

    # TODO can we drop the names from `outputs`?
    @classmethod
    def from_function(cls, function, inputs, outputs, name, doc):
        """

        Parameters
        ----------
        function : Python function
            Function to wrap as a workflow.
        inputs : dict
            Parameter name to semantic type.
        outputs : collections.OrderedDict
            Named output to semantic type.
        name : str
            Human-readable name for this workflow.
        doc : str
            Description of the function.

        """
        # TODO: the import paths will not necessarily be a public
        # import, but might be a full (private) import. fix that.
        import_path = inspect.getmodule(function).__name__
        function_name = function.__name__
        # TODO sorting by parameter name for reproducible output necessary for
        # testing. Should devs be able to define an ordering to parameters?
        parameters = ', '.join(['%s=%s' % (k, k) for k in sorted(inputs)])
        results = ', '.join(outputs.keys())
        template = _markdown_template.format(
            doc=doc, import_path=import_path, function_name=function_name,
            parameters=parameters, results=results)
        signature = Signature(name, inputs, outputs)

        # just to show that we've got the views
        input_views = {input: function.__annotations__[input]
                       for input in inputs}
        output_view_types = qiime.core.type.util.tuplize(
                           function.__annotations__['return'])
        output_views = dict(zip(outputs, output_view_types))
        # and to keep flake8 happy
        input_views, output_views

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

    def __eq__(self, other):
        return (
            self.signature == other.signature and
            self.template == other.template and
            self.id == other.id
        )

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

    def _context_lines(self, input_artifact_filepaths,
                       parameter_references, output_artifact_filepaths):
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

        return provenance_lines, setup_lines, teardown_lines, job_uuid

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

    def _format_script(self, provenance_lines, setup_lines, teardown_lines):
        provenance_str = "\n".join(provenance_lines)
        setup_str = "\n".join(setup_lines)
        teardown_str = "\n".join(teardown_lines)
        py_template = ipymd.convert(self.template, from_="markdown",
                                    to='python')
        return "\n\n".join(
            [provenance_str, setup_str, py_template, teardown_str])

    def _format_markdown_code_cell(self, code_lines):
        code_str = '\n'.join(['>>> %s' % line for line in code_lines])
        return _markdown_code_cell_template.format(content=code_str)

    def to_markdown(self, input_artifact_filepaths, parameter_references,
                    output_artifact_filepaths):
        """

        Parameters
        ----------
        input_artifact_filepaths : dict
            Input artifact names to artifact filepaths.
        parameter_references : dict
            Parameter names to values.
        output_artifact_filepaths : dict
            Output artifact names to artifact filepaths.

        Returns
        -------
        qiime.sdk.job.Job
            Job with executable ipymd markdown as `code`.

        """
        provenance_lines, setup_lines, teardown_lines, job_uuid = \
            self._context_lines(input_artifact_filepaths,
                                parameter_references,
                                output_artifact_filepaths)
        markdown = self._format_markdown(provenance_lines, setup_lines,
                                         teardown_lines)
        return qiime.sdk.job.Job(markdown, job_uuid, input_artifact_filepaths,
                                 parameter_references,
                                 output_artifact_filepaths)

    def to_script(self, input_artifact_filepaths, parameter_references,
                  output_artifact_filepaths):
        """

        Parameters
        ----------
        input_artifact_filepaths : dict
            Input artifact names to artifact filepaths.
        parameter_references : dict
            Parameter names to values.
        output_artifact_filepaths : dict
            Output artifact names to artifact filepaths.

        Returns
        -------
        qiime.sdk.job.Job
            Job with executable Python script as `code`.

        """
        provenance_lines, setup_lines, teardown_lines, job_uuid = \
            self._context_lines(input_artifact_filepaths,
                                parameter_references,
                                output_artifact_filepaths)
        script = self._format_script(provenance_lines, setup_lines,
                                     teardown_lines)
        return qiime.sdk.job.Job(script, job_uuid, input_artifact_filepaths,
                                 parameter_references,
                                 output_artifact_filepaths)


_markdown_template = """{doc}

```python
>>> from {import_path} import {function_name}
>>> {results} = {function_name}({parameters})
```
"""

_markdown_code_cell_template = """```python
{content}
```"""
