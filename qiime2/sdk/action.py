# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import abc
import concurrent.futures
import inspect
import os.path
import tempfile
import textwrap

import decorator

import qiime2.sdk
import qiime2.core.type as qtype
import qiime2.core.archive as archive
from qiime2.core.util import LateBindingAttribute, DropFirstParameter, tuplize


# These aren't globals as much as a process locals. This is necessary because
# atexit handlers aren't invoked during a subprocess exit (`_exit` is called)
_FAILURE_PROCESS_CLEANUP = []
_ALWAYS_PROCESS_CLEANUP = []


def _async_action(action, args, kwargs):
    """Helper to cleanup because atexit destructors are not called"""
    try:
        return action(*args, **kwargs)
    except:  # This is cleanup, even KeyboardInterrupt should be caught
        for garbage in _FAILURE_PROCESS_CLEANUP:
            garbage._destructor()
        raise
    finally:
        for garbage in _ALWAYS_PROCESS_CLEANUP:
            garbage._destructor()


class Action(metaclass=abc.ABCMeta):
    """QIIME 2 Action"""

    type = 'action'

    __call__ = LateBindingAttribute('_dynamic_call')
    async = LateBindingAttribute('_dynamic_async')

    # Converts a callable's signature into its wrapper's signature (i.e.
    # converts the "view API" signature into the "artifact API" signature).
    # Accepts a callable as input and returns a callable as output with
    # converted signature.
    @abc.abstractmethod
    def _callable_sig_converter_(self, callable):
        raise NotImplementedError

    # Executes a callable on the provided `view_args`, wrapping and returning
    # the callable's outputs. In other words, executes the "view API", wrapping
    # and returning the outputs as the "artifact API". `view_args` is a dict
    # mapping parameter name to unwrapped value (i.e. view). `view_args`
    # contains an entry for each parameter accepted by the wrapper. It is the
    # executor's responsibility to perform any additional transformations on
    # these parameters, or provide extra parameters, in order to execute the
    # callable. `output_types` is an OrderedDict mapping output name to QIIME
    # type (e.g. semantic type).
    @abc.abstractmethod
    def _callable_executor_(self, callable, view_args, output_types):
        raise NotImplementedError

    # Private constructor
    @classmethod
    def _init(cls, callable, signature, package, name, description):
        """

        Parameters
        ----------
        callable : callable
        signature : qiime2.core.type.Signature
        package : str
        name : str
            Human-readable name for this action.
        description : str
            Human-readable description for this action.

        """
        self = cls.__new__(cls)
        self.__init(callable, signature, package, name, description)
        return self

    # This "extra private" constructor is necessary because `Action` objects
    # can be initialized from a static (classmethod) context or on an
    # existing instance (see `_init` and `__setstate__`, respectively).
    def __init(self, callable, signature, package, name, description,
               pid=None):
        self._callable = callable
        self.signature = signature
        self.package = package
        self.name = name
        self.description = description

        self.id = callable.__name__
        if pid is None:
            pid = os.getpid()
        self._pid = pid

        self._dynamic_call = self._get_callable_wrapper()
        self._dynamic_async = self._get_async_wrapper()

    def __init__(self):
        raise NotImplementedError(
            "%s constructor is private." % self.__class__.__name__)

    @property
    def source(self):
        """

        Returns
        -------
        str
            The source code of this action's callable formatted as Markdown
            text.

        """
        try:
            source = inspect.getsource(self._callable)
        except OSError:
            raise TypeError(
                "Cannot retrieve source code for callable %r" %
                self._callable.__name__)
        return markdown_source_template % {'source': source}

    def get_import_path(self):
        return self.package + '.' + self.id

    def __repr__(self):
        return "<%s %s>" % (self.type, self.get_import_path())

    def __getstate__(self):
        return {
            'callable': self._callable,
            'signature': self.signature,
            'package': self.package,
            'name': self.name,
            'description': self.description,
            'pid': self._pid
        }

    def __setstate__(self, state):
        self.__init(**state)

    def _get_callable_wrapper(self):
        def callable_wrapper(*args, **kwargs):
            provenance = archive.ActionProvenanceCapture(
                self.type, self.package, self.id)
            if self._is_subprocess():
                _ALWAYS_PROCESS_CLEANUP.append(provenance)

            # This function's signature is rewritten below using
            # `decorator.decorator`. When the signature is rewritten, args[0]
            # is the function whose signature was used to rewrite this
            # function's signature.
            args = args[1:]

            user_input = {name: value for value, name in
                          zip(args, self.signature.ordered_parameters)}
            user_input.update(kwargs)

            self.signature.check_types(**user_input)
            output_types = self.signature.solve_output(**user_input)

            artifacts = {}
            for name in self.signature.inputs:
                artifact = artifacts[name] = user_input[name]
                provenance.add_input(name, artifact)
                if self._is_subprocess() and artifact is not None:
                    # Cleanup shouldn't be handled in the subprocess, it
                    # doesn't own any of these inputs, they were just provided.
                    # We also can't rely on the subprocess preventing atexit
                    # hooks as the destructor is also called when the artifact
                    # goes out of scope (which happens).
                    artifact._orphan()

            parameters = {}
            for name, spec in self.signature.parameters.items():
                parameter = parameters[name] = user_input[name]
                provenance.add_parameter(name, spec.qiime_type, parameter)

            view_args = parameters.copy()
            for name, spec in self.signature.inputs.items():
                recorder = provenance.transformation_recorder(name)
                artifact = artifacts[name]
                if artifact is None:
                    view_args[name] = artifact
                else:
                    view_args[name] = artifact._view(spec.view_type, recorder)

            outputs = self._callable_executor_(self._callable, view_args,
                                               output_types, provenance)
            # The outputs don't need to be orphaned, because their destructors
            # aren't invoked in atexit for a subprocess, instead the
            # `_async_action` helper will detect failure and cleanup if needed.
            # These are meant to be owned by the parent process.

            if len(outputs) != len(self.signature.outputs):
                raise ValueError(
                    "Number of callable outputs must match number of outputs "
                    "defined in signature: %d != %d" %
                    (len(outputs), len(self.signature.outputs)))

            # Wrap in a Results object mapping output name to value so users
            # have access to outputs by name or position.
            return qiime2.sdk.Results(self.signature.outputs.keys(),
                                      outputs)

        callable_wrapper = self._rewrite_wrapper_signature(callable_wrapper)
        self._set_wrapper_properties(callable_wrapper, '__call__')
        return callable_wrapper

    def _get_async_wrapper(self):
        def async_wrapper(*args, **kwargs):
            # TODO handle this better in the future, but stop the massive error
            # caused by MacOSX async runs for now.
            try:
                import matplotlib as plt
                if plt.rcParams['backend'].lower() == 'macosx':
                    raise EnvironmentError(backend_error_template %
                                           plt.matplotlib_fname())
            except ImportError:
                pass

            # This function's signature is rewritten below using
            # `decorator.decorator`. When the signature is rewritten, args[0]
            # is the function whose signature was used to rewrite this
            # function's signature.
            args = args[1:]

            pool = concurrent.futures.ProcessPoolExecutor(max_workers=1)
            future = pool.submit(_async_action, self, args, kwargs)
            pool.shutdown(wait=False)
            return future

        async_wrapper = self._rewrite_wrapper_signature(async_wrapper)
        self._set_wrapper_properties(async_wrapper, 'async')
        return async_wrapper

    def _rewrite_wrapper_signature(self, wrapper):
        # Convert the callable's signature into the wrapper's signature and set
        # it on the wrapper.
        return decorator.decorator(
            wrapper, self._callable_sig_converter_(self._callable))

    def _set_wrapper_properties(self, wrapper, name):
        wrapper.__name__ = wrapper.__qualname__ = name
        wrapper.__module__ = self.package
        wrapper.__doc__ = "{}\n\n{}".format(
            self.name, textwrap.fill(self.description, width=79))
        wrapper.__annotations__ = self._build_annotations()
        # This is necessary so that `inspect` doesn't display the wrapped
        # function's annotations (the annotations apply to the "view API" and
        # not the "artifact API").
        del wrapper.__wrapped__

    def _build_annotations(self):
        annotations = {}
        for name, spec in self.signature.ordered_parameters.items():
            annotations[name] = spec.qiime_type

        output = []
        for spec in self.signature.outputs.values():
            output.append(spec.qiime_type)
        output = tuple(output)

        annotations["return"] = output

        return annotations

    def _is_subprocess(self):
        return self._pid != os.getpid()


class Method(Action):
    """QIIME 2 Method"""

    type = 'method'

    # Abstract method implementations:

    def _callable_sig_converter_(self, callable):
        # No conversion necessary.
        return callable

    def _callable_executor_(self, callable, view_args, output_types,
                            provenance):
        is_subprocess = self._is_subprocess()
        output_views = callable(**view_args)
        output_views = tuplize(output_views)

        # TODO this won't work if the user has annotated their "view API" to
        # return a `typing.Tuple` with some number of components. Python will
        # return a tuple when there are multiple return values, and this length
        # check will fail because the tuple as a whole should be matched up to
        # a single output type instead of its components. This is an edgecase
        # due to how Python handles multiple returns, and can be worked around
        # by using something like `typing.List` instead.
        if len(output_views) != len(output_types):
            raise TypeError(
                "Number of output views must match number of output "
                "semantic types: %d != %d"
                % (len(output_views), len(output_types)))

        output_artifacts = []
        for output_view, spec in zip(output_views, output_types.values()):
            if type(output_view) is not spec.view_type:
                raise TypeError(
                    "Expected output view type %r, received %r" %
                    (spec.view_type.__name__, type(output_view).__name__))
            artifact = qiime2.sdk.Artifact._from_view(
                spec.qiime_type, output_view, spec.view_type,
                provenance.fork())
            output_artifacts.append(artifact)
            if is_subprocess:
                _FAILURE_PROCESS_CLEANUP.append(artifact._archiver)

        return tuple(output_artifacts)

    @classmethod
    def _init(cls, callable, inputs, parameters, outputs, package, name,
              description, input_descriptions, parameter_descriptions,
              output_descriptions):
        signature = qtype.MethodSignature(callable, inputs, parameters,
                                          outputs, input_descriptions,
                                          parameter_descriptions,
                                          output_descriptions)
        return super()._init(callable, signature, package, name, description)


class Visualizer(Action):
    """QIIME 2 Visualizer"""

    type = 'visualizer'

    # Abstract method implementations:

    def _callable_sig_converter_(self, callable):
        return DropFirstParameter.from_function(callable)

    def _callable_executor_(self, callable, view_args, output_types,
                            provenance):
        # TODO use qiime2.plugin.OutPath when it exists, and update visualizers
        # to work with OutPath instead of str. Visualization._from_data_dir
        # will also need to be updated to support OutPath instead of str.
        with tempfile.TemporaryDirectory(prefix='qiime2-temp-') as temp_dir:
            ret_val = callable(output_dir=temp_dir, **view_args)
            if ret_val is not None:
                raise TypeError(
                    "Visualizer %r should not return anything. "
                    "Received %r as a return value." % (self, ret_val))
            viz = qiime2.sdk.Visualization._from_data_dir(temp_dir,
                                                          provenance)
            if self._is_subprocess():
                _FAILURE_PROCESS_CLEANUP.append(viz._archiver)

            return (viz,)

    @classmethod
    def _init(cls, callable, inputs, parameters, package, name, description,
              input_descriptions, parameter_descriptions):
        signature = qtype.VisualizerSignature(callable, inputs, parameters,
                                              input_descriptions,
                                              parameter_descriptions)
        return super()._init(callable, signature, package, name, description)


markdown_source_template = """
```python
%(source)s
```
"""

# TODO add unit test for callables raising this
backend_error_template = """
Your current matplotlib backend (MacOSX) does not work with async calls.
A recommended backend is Agg, and can be changed by modifying your
matplotlibrc "backend" parameter, which can be found at: \n\n %s
"""
