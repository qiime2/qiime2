# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import uuid
import collections
import inspect
import os.path
import pydoc
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
    def __init__(self, name, inputs, parameters, outputs):
        """

        Parameters
        ----------
        name : str
            Human-readable name for this workflow.
        inputs : dict
            Parameter name to tuple of semantic type and view type.
        parameters : dict
            Parameter name to tuple of primitive type and view type.
        outputs : collections.OrderedDict
            Named output to tuple of semantic type and view type.

        """
        for input_name, (semantic_type, input_view_type) in \
                inputs.items():
            if not qiime.core.type.is_semantic_type(semantic_type):
                raise TypeError("%r for %r is not a semantic qiime type." %
                                (semantic_type, input_name))

        for param_name, (primitive_type, param_view_type) in \
                parameters.items():
            if not qiime.core.type.is_primitive_type(primitive_type):
                raise TypeError("%r for %r is not a primitive qiime type."
                                % (primitive_type, param_name))

        for output_name, (output_semantic_type, output_view_type) in \
                outputs.items():
            if not qiime.core.type.is_semantic_type(output_semantic_type) and \
               output_semantic_type is not qiime.core.type.Visualization:
                raise TypeError("%r for %r is not a semantic qiime type or "
                                "visualization." %
                                (output_semantic_type, output_name))

        self.name = name
        self.inputs = inputs
        self.parameters = parameters
        self.outputs = outputs

    def __call__(self, artifacts, arguments):
        # TODO implement me!
        return self.outputs

    def __eq__(self, other):
        return (
            type(self) is type(other) and
            self.name == other.name and
            self.inputs == other.inputs and
            self.parameters == other.parameters and
            self.outputs == other.outputs
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

        input_types = {}
        for input_ in metadata['inputs']:
            # TODO validate each nested dict has exactly two items
            name, type_expr = list(input_.items())[0]
            input_types[name] = cls._parse_semantic_type(type_expr)

        parameter_types = {}
        for parameter in metadata['parameters']:
            # TODO validate each nested dict has exactly two items
            name, type_expr = list(parameter.items())[0]
            parameter_types[name] = cls._parse_primitive_type(type_expr)

        output_types = collections.OrderedDict()
        for output in metadata['outputs']:
            # TODO validate each nested dict has exactly two items
            name, type_expr = list(output.items())[0]
            output_types[name] = cls._parse_semantic_type(type_expr)

        name = metadata['name']
        signature = Signature(name, input_types, parameter_types, output_types)
        return cls(signature, template, id_)

    @classmethod
    def _parse_primitive_type(cls, type_expr):
        # TODO: what is factoring?
        primitive_type_expr, view_type = type_expr
        view_type = cls._parse_view_type(view_type)

        # TODO: Get these primitives from somewhere else.
        locals_ = {
            t.name: t for t in
            {qiime.plugin.Int, qiime.plugin.Str, qiime.plugin.Float}}

        type_ = eval(primitive_type_expr, {'__builtins__': {}}, locals_)
        return type_, view_type

    # TODO this is mostly duplicated from Artifact._parse_type. Refactor!
    @classmethod
    def _parse_semantic_type(cls, type_exp):
        # Split the type expression into its components: the semantic_type_exp
        # and the view_type. Note that this differs from the type definitions
        # in Artifact._parse_type, as those won't have view types.
        semantic_type_exp, view_type = type_exp
        view_type = cls._parse_view_type(view_type)

        semantic_type_exp = semantic_type_exp.split('\n')
        if len(semantic_type_exp) != 1:
            raise TypeError("Multiple lines in type expression of"
                            " artifact. Will not load to avoid arbitrary"
                            " code execution.")
        semantic_type_exp, = semantic_type_exp

        if ';' in semantic_type_exp:
            raise TypeError("Invalid type expression in artifact. Will not"
                            " load to avoid arbitrary code execution.")

        pm = qiime.sdk.PluginManager()
        locals_ = {k: v[1] for k, v in pm.semantic_types.items()}
        # Set up all of the types we know about in local scope of the eval
        # so that complicated type expressions are evaluated.
        type_ = eval(semantic_type_exp, {'__builtins__': {}}, locals_)
        return type_, view_type

    # TODO this isn't safe nor robust, and not how we want to ultimately handle
    # importing a view type (`_setup_import_lines` has a similar problem).
    # Refactor to use a similar import mechanism as semantic types when that
    # part is finalized. Hack taken from:
    #     http://stackoverflow.com/a/29831586/3776794
    @classmethod
    def _parse_view_type(cls, view_type_str):
        view_type = pydoc.locate(view_type_str)
        if view_type is None:
            raise ImportError("Could not import view type %r" % view_type_str)
        return view_type

    # TODO can we drop the names from `outputs`?
    @classmethod
    def from_function(cls, function, inputs, parameters, outputs, name, doc):
        """

        Parameters
        ----------
        function : Python function
            Function to wrap as a workflow.
        inputs : dict
            Parameter name to semantic type.
        parameters : dict
            Parameter name to primitive type.
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
        # testing. This should match `function` signature
        md_params = ', '.join([
            '%s=%s' % (k, k) for k in sorted(inputs) + sorted(parameters)])

        input_types = {}
        for param_name, semantic_type in inputs.items():
            view_type = function.__annotations__[param_name]
            input_types[param_name] = (semantic_type, view_type)

        param_types = {}
        for param_name, primitive_type in parameters.items():
            view_type = function.__annotations__[param_name]
            param_types[param_name] = (primitive_type, view_type)

        output_view_types = qiime.core.type.util.tuplize(
                           function.__annotations__['return'])
        output_view_types = dict(zip(outputs, output_view_types))

        output_types = collections.OrderedDict()
        for output_name, semantic_type in outputs.items():
            view_type = output_view_types[output_name]
            output_types[output_name] = (semantic_type, view_type)

        if function.__annotations__['return'] is None:
            template = _markdown_template_without_outputs.format(
                doc=doc, import_path=import_path, function_name=function_name,
                parameters=md_params)
        else:
            results = ', '.join(outputs.keys())
            template = _markdown_template_with_outputs.format(
                doc=doc, import_path=import_path, function_name=function_name,
                parameters=md_params, results=results)

        signature = Signature(name, input_types, param_types, output_types)

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
                    self.signature.parameters[name][0].decode(ref)

        job_uuid = uuid.uuid4()
        provenance_lines = self._provenance_lines(job_uuid, artifacts,
                                                  parameters)

        setup_lines = self._setup_import_lines()

        for name, filepath in input_artifact_filepaths.items():
            view_type = self.signature.inputs[name][1]
            # TODO explicitly clean up Artifact object instead of relying on
            # gc?
            full_view_type_name = '.'.join([_get_import_path(view_type),
                                            view_type.__name__])
            setup_lines.append(
                '%s = Artifact.load(%r).view(%s)' % (name, filepath,
                                                     full_view_type_name))

        for name, value in parameters.items():
            setup_lines.append('%s = %r' % (name, value))

        teardown_lines = []

        # TODO make sure output order is respected
        for name, output_filepath in output_artifact_filepaths.items():
            output_semantic_type, _ = output_artifact_types[name]
            if output_semantic_type == qiime.core.type.Visualization:
                teardown_lines.append(
                    'from qiime.sdk.visualization import Visualization')
                # visualizations will by definition have an output_dir
                # parameter, so that is hard-coded here.
                teardown_lines.append(
                    'visualization = Visualization._from_data_dir('
                    'output_dir, provenance)'
                )
                teardown_lines.append(
                    'visualization.save(%r)' % output_filepath
                )
            else:
                teardown_lines.append(
                    # TODO make sure `artifact` is a reserved local variable
                    # name so we don't shadow a plugin's local vars. Are there
                    # other reserved names the framework uses?
                    'artifact = Artifact._from_view(%s, %r, provenance)' %
                    (name, str(output_semantic_type)))
                # TODO explicitly clean up Artifact object instead of relying
                # on gc?
                teardown_lines.append(
                    'artifact.save(%r)' % output_filepath)

        return provenance_lines, setup_lines, teardown_lines, job_uuid

    # TODO this isn't safe nor robust, and not how we want to ultimately handle
    # importing a view type (`_parse_view_type` has a similar problem).
    # Refactor to use a similar import mechanism as semantic types when that
    # part is finalized.
    def _setup_import_lines(self):
        lines = set()

        for _, input_view_type in self.signature.inputs.values():
            lines.add('import %s' % _get_import_path(input_view_type))

        lines.add('from qiime.sdk.artifact import Artifact')
        lines = list(lines)
        lines.sort(key=lambda x: len(x))

        return lines

    def _load_artifacts(self, artifact_filepaths):
        artifacts = {}
        for name, filepath in artifact_filepaths.items():
            artifacts[name] = qiime.sdk.Artifact.load(filepath)
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


def _get_import_path(cls):
    # wrapping this in its own function now so we can do more complex things
    # later, like find the shortest import path
    return inspect.getmodule(cls).__name__

_markdown_template_with_outputs = """{doc}

```python
>>> from {import_path} import {function_name}
>>> {results} = {function_name}({parameters})
```
"""

_markdown_template_without_outputs = """{doc}

```python
>>> from {import_path} import {function_name}
>>> {function_name}({parameters})
```
"""

_markdown_code_cell_template = """```python
{content}
```"""
