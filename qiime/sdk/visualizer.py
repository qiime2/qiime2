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
import tempfile
import uuid

import decorator

import qiime.sdk
import qiime.core.type as qtype


# Descriptor protocol for methods with dynamic signatures built from the
# callables they wrap.
class DispatchableSignature:
    def __init__(self, method):
        self._method = method

    def __get__(self, obj, cls=None):
        return staticmethod(getattr(obj, self._method)).__get__(obj, cls)


class DropFirstParameter(decorator.FunctionMaker):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.signature = self._remove_first_arg(self.signature)
        self.shortsignature = self._remove_first_arg(self.shortsignature)

    def _remove_first_arg(self, string):
        return ",".join(string.split(',')[1:])[1:]

    @classmethod
    def from_function(cls, function):
        return cls.create(function, "return None", {})


# TODO a ton of code is shared between Method and Visualizer, refactor!
class Visualizer:
    __call__ = DispatchableSignature('_dynamic_call')
    async = DispatchableSignature('_dynamic_async')

    def __new__(cls):
        self = super().__new__(cls)
        self._pid = os.getpid()
        return self

    @classmethod
    def from_function(cls, function, inputs, parameters, name, description,
                      plugin=None):
        if 'output_dir' in inputs or 'output_dir' in parameters:
            raise TypeError(
                "`output_dir` is a reserved parameter name and cannot be used "
                "in `inputs` or `parameters`")

        function_parameters = \
            list(inspect.signature(function).parameters.keys())
        if len(function_parameters) > 0:
            first_parameter = function_parameters[0]
            if first_parameter != 'output_dir':
                raise TypeError(
                    "Visualizer function must have `output_dir` as its first "
                    "argument, not %r" % first_parameter)
        else:
            raise TypeError(
                "Visualizer function must have at least one argument")

        if ('return' in function.__annotations__ and
                function.__annotations__['return'] is not None):
            raise TypeError(
                "Visualizer function %r cannot return anything. Write output "
                "to `output_dir`." % function.__name__)

        input_types = {}
        for param_name, semantic_type in inputs.items():
            view_type = function.__annotations__[param_name]
            input_types[param_name] = (semantic_type, view_type)

        param_types = {}
        for param_name, primitive_type in parameters.items():
            view_type = function.__annotations__[param_name]
            param_types[param_name] = (primitive_type, view_type)

        output_types = collections.OrderedDict(
            [('visualization', (qtype.Visualization, None))])
        signature = qtype.Signature(input_types, param_types, output_types)

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
                   description, markdown_source, plugin)
        return self

    def __init__(self):
        raise NotImplementedError(
            "%(clsname)s constructor is private. Use "
            "%(clsname)s.from_function." %
            {'clsname': self.__class__.__name__})

    def _init(self, id, signature, callable, callable_ref, name,
              description, source, plugin):
        """

        Parameters
        ----------
        id_ : str
        signature : Signature
        callable_ : callable
        callable_ref: (str, object)
        name : str
            Human-readable name for this visualizer.
        description : str
            Human-readable description for this visualizer.
        source : str
            Markdown text defining/describing this visualizer's computation.
        plugin : str
            Name of the plugin this visualizer is registered to.
        """
        self.id = id
        self.signature = signature
        self._callable = callable
        self._callable_ref = callable_ref
        self.name = name
        self.description = description
        self.source = source
        self._plugin = plugin

        self._bind_executors()

    def _bind_executors(self):
        callable_wrapper = self._get_callable_wrapper()

        __call__ = decorator.decorator(
            callable_wrapper,
            DropFirstParameter.from_function(self._callable))
        __call__.__name__ = __call__.__qualname__ = '__call__'
        __call__.__module__ = self._get_import_path()
        del __call__.__annotations__
        del __call__.__wrapped__
        self._dynamic_call = __call__

        # TODO `async` execution has some problems with garbage-collection in
        # the subprocess given certain inputs. Needs further investigation.
        def async_wrapper(*args, **kwargs):
            args = args[1:]
            pool = concurrent.futures.ProcessPoolExecutor(max_workers=1)
            future = pool.submit(self, *args, **kwargs)
            pool.shutdown(wait=False)
            return future

        async = decorator.decorator(
            async_wrapper,
            DropFirstParameter.from_function(self._callable))
        async.__name__ = async.__qualname__ = 'async'
        async.__module__ = self._get_import_path()
        del async.__annotations__
        del async.__wrapped__
        self._dynamic_async = async

    def _get_callable_wrapper(self):
        def callable_wrapper(*args, **kwargs):
            args = args[1:]
            callable_sig = inspect.Signature.from_callable(self._callable)
            callable_params = collections.OrderedDict(
                callable_sig.parameters.items())
            del callable_params['output_dir']

            args = {name: value for value, name in zip(args, callable_params)}
            args.update(kwargs)

            self.signature.check_types(**args)

            artifacts = {}
            for name in self.signature.inputs:
                artifacts[name] = args[name]

            parameters = {}
            for name in self.signature.parameters:
                parameters[name] = args[name]

            view_args = parameters.copy()
            for name, (_, view_type) in self.signature.inputs.items():
                view_args[name] = artifacts[name].view(view_type)

            execution_uuid = uuid.uuid4()
            visualizer_reference = (
                "%s. Details on plugin, version, website, etc. will also be "
                "included, see https://github.com/biocore/qiime2/issues/26"
                % self.id)
            artifact_uuids = {name: artifact.uuid for name, artifact in
                              artifacts.items()}
            parameter_references = {name: str(param) for name, param in
                                    parameters.items()}
            provenance = qiime.sdk.Provenance(
                execution_uuid, visualizer_reference, artifact_uuids,
                parameter_references)

            # TODO use user-configured temp dir
            with tempfile.TemporaryDirectory('qiime2-temp-') as temp_dir:
                # TODO make sure _callable doesn't return anything
                self._callable(output_dir=temp_dir, **view_args)
                visualization = qiime.sdk.Visualization._from_data_dir(
                    temp_dir, provenance)
                visualization._orphan(self._pid)
                return visualization

        return callable_wrapper

    def __repr__(self):
        return "<visualizer %s>" % self._get_import_path()

    def _get_import_path(self):
        plugin_path = ''
        if self._plugin is not None:
            plugin_path = ('qiime.plugin.%s.visualizers.'
                           % self._plugin.replace('-', '_'))
        return plugin_path + self.id

    def __getstate__(self):
        return {
            'type': self._callable_ref[0],
            'src': self._callable_ref[1],
            'attrs': {
                'id': self.id,
                'signature': self.signature,
                'name': self.name,
                'description': self.description,
                'source': self.source,
                'plugin': self._plugin
            },
            'pid': self._pid,
        }

    def __setstate__(self, state):
        if state['type'] == 'function':
            attrs = state['attrs']
            self._init(callable=state['src'],
                       callable_ref=(state['type'], state['src']),
                       **attrs)
        else:
            raise NotImplementedError
        self._pid = state['pid']


markdown_source_template = """
```python
%(source)s
```
"""
