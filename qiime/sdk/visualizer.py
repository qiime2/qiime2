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
import qiime.core.type


# Descriptor protocol for methods with dynamic signatures built from the
# callables they wrap.
class DispatchableSignature:
    def __init__(self, method):
        self._method = method

    def __get__(self, obj, cls=None):
        return staticmethod(getattr(obj, self._method)).__get__(obj, cls)


class FunctionMakerDropFirstArg(decorator.FunctionMaker):
    def __init__(self, *args, **kwargs):
        super(FunctionMakerDropFirstArg, self).__init__(*args, **kwargs)
        self.signature = self._remove_first_arg(self.signature)
        self.shortsignature = self._remove_first_arg(self.shortsignature)

    def _remove_first_arg(self, string):
        return ",".join(string.split(',')[1:])[1:]

    @classmethod
    def clone_sig_without_first_arg(cls, function):
        evaldict = {}
        fun = cls.create(function, "return None", evaldict, 
                         __wrapped__=function)
        fun.__func__ = function  # Doctests need the orginal function
        fun.__annotations__ = {}
        return fun


# TODO a ton of code is shared between Method and Visualizer, refactor!
class Visualizer:
    __call__ = DispatchableSignature('_dynamic_call')
    async = DispatchableSignature('_dynamic_async')

    def __new__(cls):
        self = super().__new__(cls)
        self._pid = os.getpid()
        return self

    @classmethod
    def from_function(cls, function, inputs, parameters, name, description):
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
                    "Visualization function must have `output_dir` as its "
                    "first argument, not %r" % first_parameter)
        else:
            raise TypeError(
                "Visualization function must have at least one argument")

        if function.__annotations__['return'] is not None:
            raise TypeError(
                "Visualization function cannot return anything. Write "
                "visualization output to `output_dir`.")

        input_types = {}
        for param_name, semantic_type in inputs.items():
            view_type = function.__annotations__[param_name]
            input_types[param_name] = (semantic_type, view_type)

        param_types = {}
        for param_name, primitive_type in parameters.items():
            view_type = function.__annotations__[param_name]
            param_types[param_name] = (primitive_type, view_type)

        output_types = collections.OrderedDict(
            [('visualization', (qiime.core.type.Visualization, None))])
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

    def __init__(self):
        raise NotImplementedError(
            "%(clsname)s constructor is private. Use "
            "%(clsname)s.from_function." %
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
            Human-readable name for this visualizer.
        description : str
            Human-readable description for this visualizer.
        source : str
            Markdown text defining/describing this visualizer's computation.

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

        # TODO the signature (as `inspect` and IPython sees it) isn't correct,
        # it displays `output_dir` when there shouldn't be one
        __call__ = decorator.decorator(
            callable_wrapper,
            FunctionMakerDropFirstArg.clone_sig_without_first_arg(self._callable))
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

        async = decorator.decorator(async_wrapper,
            FunctionMakerDropFirstArg.clone_sig_without_first_arg(self._callable))
        async.__name__ = 'async'
        self._dynamic_async = async

    def _get_callable_wrapper(self):
        def callable_wrapper(*args, **kwargs):
            args = args[1:]
            callable_sig = inspect.Signature.from_callable(self._callable)
            callable_params = dict(callable_sig.parameters)
            del callable_params['output_dir']

            args = {name: value for value, name in zip(args, callable_params)}
            args.update(kwargs)

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
            method_reference = (
                "%s. Details on plugin, version, website, etc. will also be "
                "included, see https://github.com/biocore/qiime2/issues/26 "
                % self.id)
            artifact_uuids = {name: str(artifact.uuid) for name, artifact in
                              artifacts.items()}
            provenance = qiime.sdk.Provenance(execution_uuid, method_reference,
                                              artifact_uuids, parameters)

            # TODO use sane prefix
            with tempfile.TemporaryDirectory() as temp_dir:
                # TODO make sure _callable doesn't return anything
                self._callable(output_dir=temp_dir, **view_args)
                visualization = qiime.sdk.Visualization._from_data_dir(temp_dir, provenance)
                visualization._orphan(self._pid)
                return visualization

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
            isinstance(other, Visualizer) and
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

