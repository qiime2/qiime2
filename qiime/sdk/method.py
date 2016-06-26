# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import concurrent.futures
import inspect
import os.path
import pydoc
import textwrap
import uuid

import decorator
import frontmatter
import ipymd

import qiime.sdk


# Descriptor protocol for methods with dynamic signatures built from the
# callables they wrap.
class DispatchableSignature:
    def __init__(self, method):
        self._method = method

    def __get__(self, obj, cls=None):
        return staticmethod(getattr(obj, self._method)).__get__(obj, cls)


class Method:
    __call__ = DispatchableSignature('_dynamic_call')
    async = DispatchableSignature('_dynamic_async')

    def __new__(cls):
        self = super().__new__(cls)
        self._pid = os.getpid()
        return self

    # TODO can we drop the names from `outputs`?
    @classmethod
    def from_function(cls, function, inputs, parameters, outputs, name,
                      description):
        """

        Parameters
        ----------
        function : Python function
            Function defining method computation.
        inputs : dict
            Parameter name to semantic type.
        parameters : dict
            Parameter name to primitive type.
        outputs : collections.OrderedDict consumable
            Named output to semantic type.
        name : str
            Human-readable name for this method.
        description : str
            Human-readable description for this method.

        """
        input_types = {}
        for param_name, semantic_type in inputs.items():
            view_type = function.__annotations__[param_name]
            input_types[param_name] = (semantic_type, view_type)

        param_types = {}
        for param_name, primitive_type in parameters.items():
            view_type = function.__annotations__[param_name]
            param_types[param_name] = (primitive_type, view_type)

        outputs = collections.OrderedDict(outputs)
        output_view_types = qiime.core.util.tuplize(
                           function.__annotations__['return'])
        output_view_types = dict(zip(outputs, output_view_types))

        output_types = collections.OrderedDict()
        for output_name, semantic_type in outputs.items():
            view_type = output_view_types[output_name]
            output_types[output_name] = (semantic_type, view_type)

        signature = qiime.sdk.Signature(input_types, param_types, output_types)

        id_ = function.__name__
        try:
            source = inspect.getsource(function)
        except OSError:
            raise TypeError(
                "Cannot retrieve source code for function %r" %
                function.__name__)
        markdown_source = markdown_source_template % {'source': source}

        self = cls.__new__(cls)
        self._init(id_, signature, function, ('function', function), name,
                   description, markdown_source)
        return self

    @classmethod
    def from_markdown(cls, markdown_filepath):
        self = cls.__new__(cls)
        self._from_markdown(markdown_filepath)
        return self

    def _from_markdown(self, markdown_filepath):
        with open(markdown_filepath) as fh:
            metadata, template = frontmatter.parse(fh.read())

        input_types = collections.OrderedDict()
        for input_ in metadata['inputs']:
            # TODO validate each nested dict has exactly two items
            name, type_expr = list(input_.items())[0]
            input_types[name] = self._parse_semantic_type(type_expr)

        param_types = collections.OrderedDict()
        for parameter in metadata['parameters']:
            # TODO validate each nested dict has exactly two items
            name, type_expr = list(parameter.items())[0]
            param_types[name] = self._parse_primitive_type(type_expr)

        output_types = collections.OrderedDict()
        for output in metadata['outputs']:
            # TODO validate each nested dict has exactly two items
            name, type_expr = list(output.items())[0]
            output_types[name] = self._parse_semantic_type(type_expr)

        signature = qiime.sdk.Signature(input_types, param_types, output_types)

        # TODO: verify that `id_` is a valid Python identifier
        id_ = os.path.splitext(os.path.basename(markdown_filepath))[0]

        # TODO handle default values for optional parameters when that's
        # supported
        function_def_line = 'def %s(%s, %s):' % (id_,
                                                 ', '.join(input_types),
                                                 ', '.join(param_types))
        indent = ' ' * 4
        function_body = ipymd.convert(template, from_='markdown', to='python')
        function_body = textwrap.indent(function_body, indent)
        function_return_line = '%sreturn %s' % (indent,
                                                ', '.join(output_types))

        function_str = '\n'.join([function_def_line,
                                  function_body,
                                  function_return_line])

        scope = {}
        exec(function_str, scope)
        function = scope[id_]

        name = metadata['name']
        description = metadata['description']

        self._init(id_, signature, function, ('markdown', markdown_filepath),
                   name, description, template)

    # TODO this is mostly duplicated from
    # qiime.core.archiver.Archiver._parse_type. Refactor!
    @classmethod
    def _parse_semantic_type(cls, type_exp):
        # Split the type expression into its components: the semantic_type_exp
        # and the view_type. Note that this differs from the type definitions
        # in Archiver._parse_type, as those won't have view types.
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

    # TODO this isn't safe nor robust, and not how we want to ultimately handle
    # importing a view type. Refactor to use a similar import mechanism as
    # semantic types when that part is finalized. Hack taken from:
    #     http://stackoverflow.com/a/29831586/3776794
    @classmethod
    def _parse_view_type(cls, view_type_str):
        view_type = pydoc.locate(view_type_str)
        if view_type is None:
            raise ImportError("Could not import view type %r" % view_type_str)
        return view_type

    def __init__(self):
        raise NotImplementedError(
            "%(clsname)s constructor is private. Use "
            "%(clsname)s.from_function or %(clsname)s.from_markdown." %
            {'clsname': self.__class__.__name__})

    def _init(self, id, signature, callable, callable_ref, name,
              description, source):
        """

        Parameters
        ----------
        id_ : str
        signature : Signature
        callable_ : callable
        callable_ref: (str, object)
        name : str
            Human-readable name for this method.
        description : str
            Human-readable description for this method.
        source : str
            Markdown text defining/describing this method's computation.

        """
        self.id = id
        self.signature = signature
        self._callable = callable
        self._callable_ref = callable_ref
        self.name = name
        self.description = description
        self.source = source

        self._bind_executors()

    def _bind_executors(self):
        callable_wrapper = self._get_callable_wrapper()

        # TODO drop function annotations as they only make sense for the
        # "view API". Necessary for `async` too. Simply setting __annotations__
        # to {} or None doesn't work for some reason.
        __call__ = decorator.decorator(callable_wrapper, self._callable)
        __call__.__name__ = '__call__'
        self._dynamic_call = __call__

        # TODO `async` execution has some problems with garbage-collection in
        # the subprocess given certain inputs. Needs further investigation.
        def async_wrapper(*args, **kwargs):
            args = args[1:]
            pool = concurrent.futures.ProcessPoolExecutor(max_workers=1)
            future = pool.submit(self, *args, **kwargs)
            pool.shutdown(wait=False)
            return future

        async = decorator.decorator(async_wrapper, self._callable)
        async.__name__ = 'async'
        self._dynamic_async = async

    def _get_callable_wrapper(self):
        def callable_wrapper(*args, **kwargs):
            args = args[1:]
            callable_sig = inspect.Signature.from_callable(self._callable)
            callable_params = callable_sig.parameters

            args = {name: value for value, name in zip(args, callable_params)}
            args.update(kwargs)

            artifacts = {}
            for name in self.signature.inputs:
                artifacts[name] = args[name]

            parameters = {}
            for name in self.signature.parameters:
                parameters[name] = args[name]

            output_types = self.signature(artifacts, parameters)

            view_args = parameters.copy()
            for name, (_, view_type) in self.signature.inputs.items():
                view_args[name] = artifacts[name].view(view_type)

            output_views = self._callable(**view_args)

            if type(output_views) is not tuple:
                output_views = (output_views,)

            if len(output_views) != len(output_types):
                raise TypeError(
                    "Number of output views must match number of output "
                    "semantic types: %d != %d"
                    % (len(output_views), len(output_types)))

            execution_uuid = uuid.uuid4()
            method_reference = (
                "%s. Details on plugin, version, website, etc. will also be "
                "included, see https://github.com/biocore/qiime2/issues/26 "
                % self.id)
            artifact_uuids = {name: str(artifact.uuid) for name, artifact in
                              artifacts.items()}
            provenance = qiime.sdk.Provenance(execution_uuid, method_reference,
                                              artifact_uuids, parameters)

            output_artifacts = []
            for output_view, (semantic_type, view_type) in \
                    zip(output_views, output_types.values()):
                if not type(output_view) is view_type:
                    raise TypeError(
                        "Expected output view type %r, received %r" %
                        (view_type.__name__, type(output_view).__name__))
                artifact = qiime.sdk.Artifact._from_view(
                    output_view, semantic_type, provenance)
                artifact._orphan(self._pid)
                output_artifacts.append(artifact)

            if len(output_artifacts) == 0:
                return None
            elif len(output_artifacts) == 1:
                return output_artifacts[0]
            else:
                return tuple(output_artifacts)

        return callable_wrapper

    def __getstate__(self):
        return {
            'type': self._callable_ref[0],
            'src': self._callable_ref[1],
            'attrs': {
                'id': self.id,
                'signature': self.signature,
                'name': self.name,
                'description': self.description,
                'source': self.source
            },
            'pid': self._pid,
        }

    def __setstate__(self, state):
        if state['type'] == 'function':
            attrs = state['attrs']
            self._init(callable=state['src'],
                       callable_ref=(state['type'], state['src']),
                       **attrs)
        elif state['type'] == 'markdown':
            markdown_filepath = state['src']
            self._from_markdown(markdown_filepath)
        else:
            raise NotImplementedError
        self._pid = state['pid']

    def __eq__(self, other):
        return (
            isinstance(other, Method) and
            self.id == other.id and
            self.signature == other.signature and
            self._callable_ref == other._callable_ref and
            self.name == other.name and
            self.description == other.description and
            self.source == other.source
        )


markdown_source_template = """
```python
%(source)s
```
"""
