# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import abc
import concurrent.futures
import inspect
import tempfile
import textwrap
import itertools

import decorator
import dill
import parsl
from parsl.app.app import python_app, join_app

import qiime2.sdk
import qiime2.core.type as qtype
import qiime2.core.archive as archive
from qiime2.core.util import LateBindingAttribute, DropFirstParameter, tuplize
from qiime2.sdk.parsl_config import PARSL_CONFIG, get_config
from qiime2.sdk.context import Context
from qiime2.sdk.results import Results
from qiime2.sdk.util import ProxyArtifact
from qiime2.sdk.cache_config import CACHE_CONFIG


def _subprocess_apply(action, args, kwargs):
    # Preprocess input artifacts as we've got pickled clones which shouldn't
    # self-destruct.
    for arg in itertools.chain(args, kwargs.values()):
        if isinstance(arg, qiime2.sdk.Artifact):
            # We can't rely on the subprocess preventing atexit hooks as the
            # destructor is also called when the artifact goes out of scope
            # (which happens).
            arg._destructor.detach()

    results = action(*args, **kwargs)
    for r in results:
        # The destructor doesn't keep its detatched state when sent back to the
        # main process. Something about the context-manager from ctx seems to
        # cause a GC of the artifacts before the process actually ends, so we
        # do need to detach these. The specifics are not understood.
        r._destructor.detach()
    return results


def run_parsl_action(action, ctx, args, kwargs, inputs=[]):
    import qiime2.sdk.context

    remapped_kwargs = {}
    for key, value in kwargs.items():
        if isinstance(value, qiime2.sdk.util.ProxyArtifact):
            remapped_kwargs[key] = value.get_element(inputs[value.future])
        else:
            remapped_kwargs[key] = value

    remapped_args = []
    for arg in args:
        if isinstance(arg, qiime2.sdk.util.ProxyArtifact):
            remapped_args.append(arg.get_element(inputs[arg.future]))
        else:
            remapped_args.append(arg)

    exe = action._bind(lambda: Context(parent=ctx))
    return exe(*remapped_args, **remapped_kwargs)


class Action(metaclass=abc.ABCMeta):
    """QIIME 2 Action"""
    type = 'action'
    _ProvCaptureCls = archive.ActionProvenanceCapture

    __call__ = LateBindingAttribute('_dynamic_call')
    asynchronous = LateBindingAttribute('_dynamic_async')
    parsl = LateBindingAttribute('_dynamic_parsl')

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
    def _callable_executor_(self, scope, view_args, output_types):
        raise NotImplementedError

    # Private constructor
    @classmethod
    def _init(cls, callable, signature, plugin_id, name, description,
              citations, deprecated, examples):
        """

        Parameters
        ----------
        callable : callable
        signature : qiime2.core.type.Signature
        plugin_id : str
        name : str
            Human-readable name for this action.
        description : str
            Human-readable description for this action.

        """
        self = cls.__new__(cls)
        self.__init(callable, signature, plugin_id, name, description,
                    citations, deprecated, examples)
        return self

    # This "extra private" constructor is necessary because `Action` objects
    # can be initialized from a static (classmethod) context or on an
    # existing instance (see `_init` and `__setstate__`, respectively).
    def __init(self, callable, signature, plugin_id, name, description,
               citations, deprecated, examples):
        self._callable = callable
        self.signature = signature
        self.plugin_id = plugin_id
        self.name = name
        self.description = description
        self.citations = citations
        self.deprecated = deprecated
        self.examples = examples

        self.id = callable.__name__
        self._dynamic_call = self._get_callable_wrapper()
        self._dynamic_async = self._get_async_wrapper()
        # This a temp thing to play with parsl before integrating more deeply
        self._dynamic_parsl = self._get_parsl_wrapper()

    def __init__(self):
        raise NotImplementedError(
            "%s constructor is private." % self.__class__.__name__)

    @property
    def source(self):
        """
        The source code for the action's callable.

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

    def get_import_path(self, include_self=True):
        path = f'qiime2.plugins.{self.plugin_id}.{self.type}s'
        if include_self:
            path += f'.{self.id}'
        return path

    def __repr__(self):
        return "<%s %s>" % (self.type, self.get_import_path())

    def __getstate__(self):
        return dill.dumps({
            'callable': self._callable,
            'signature': self.signature,
            'plugin_id': self.plugin_id,
            'name': self.name,
            'description': self.description,
            'citations': self.citations,
            'deprecated': self.deprecated,
            'examples': self.examples,
        })

    def __setstate__(self, state):
        self.__init(**dill.loads(state))

    # Use bind to bind parsl config to action when action is being sent to
    # executor
    def _bind(self, context_factory):
        """Bind an action to a Context factory, returning a decorated function.

        This is a very primitive API and should be used primarily by the
        framework and very advanced interfaces which need deep control over
        the calling semantics of pipelines and garbage collection.

        The basic idea behind this is outlined as follows:

        Every action is defined as an *instance* that a plugin constructs.
        This means that `self` represents the internal details as to what
        the action is. If you need to associate additional state with the
        *application* of an action, you cannot mutate `self` without
        changing all future applications. So there needs to be an
        additional instance variable that can serve as the state of a given
        application. We call this a Context object. It is also important
        that each application of an action has *independent* state, so
        providing an instance of Context won't work. We need a factory.

        Parameterizing the context is necessary because it is possible for
        an action to call other actions. The details need to be coordinated
        behind the scenes to the user, so we can parameterize the behavior
        by providing different context factories to `bind` at different
        points in the "call stack".

        """
        def bound_callable(*args, **kwargs):
            # This function's signature is rewritten below using
            # `decorator.decorator`. When the signature is rewritten,
            # args[0] is the function whose signature was used to rewrite
            # this function's signature.
            args = args[1:]
            ctx = context_factory()
            # Set up a scope under which we can track destructable references
            # if something goes wrong, the __exit__ handler of this context
            # manager will clean up. (It also cleans up when things go right)
            with ctx as scope:
                provenance = self._ProvCaptureCls(
                    self.type, self.plugin_id, self.id)
                scope.add_reference(provenance)

                # Collate user arguments
                user_input = {name: value for value, name in
                              zip(args, self.signature.signature_order)}
                user_input.update(kwargs)

                # Type management
                self.signature.check_types(**user_input)
                output_types = self.signature.solve_output(**user_input)
                callable_args = {}

                # Record parameters
                for name, spec in self.signature.parameters.items():
                    parameter = callable_args[name] = user_input[name]
                    provenance.add_parameter(name, spec.qiime_type, parameter)

                # Record and transform inputs
                for name, spec in self.signature.inputs.items():
                    artifact = user_input[name]
                    provenance.add_input(name, artifact)
                    if artifact is None:
                        callable_args[name] = None
                    elif spec.has_view_type():
                        recorder = provenance.transformation_recorder(name)
                        if qtype.is_collection_type(spec.qiime_type):
                            # Always put in a list. Sometimes the view isn't
                            # hashable, which isn't relevant, but would break
                            # a Set[SomeType].
                            callable_args[name] = [
                                a._view(spec.view_type, recorder)
                                for a in user_input[name]]
                        else:
                            callable_args[name] = artifact._view(
                                spec.view_type, recorder)
                    else:
                        callable_args[name] = artifact

                if self.deprecated:
                    with qiime2.core.util.warning() as warn:
                        warn(self._build_deprecation_message(),
                             FutureWarning)

                # Wrap in a Results object mapping output name to value so
                # users have access to outputs by name or position.
                if isinstance(self, Pipeline) and ctx.parsl:
                    return self._parsl_callable_executor_(
                        scope, callable_args, output_types, provenance)

                return self._callable_executor_(
                    scope, callable_args, output_types, provenance)

        bound_callable = self._rewrite_wrapper_signature(bound_callable)
        self._set_wrapper_properties(bound_callable)
        self._set_wrapper_name(bound_callable, self.id)
        return bound_callable

    def _get_callable_wrapper(self):
        # This is a "root" level invocation (not a nested call within a
        # pipeline), so no special factory is needed.
        callable_wrapper = self._bind(qiime2.sdk.Context)
        self._set_wrapper_name(callable_wrapper, '__call__')
        return callable_wrapper

    def _get_async_wrapper(self):
        def async_wrapper(*args, **kwargs):
            # TODO handle this better in the future, but stop the massive error
            # caused by MacOSX asynchronous runs for now.
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
            future = pool.submit(_subprocess_apply, self, args, kwargs)
            # TODO: pool.shutdown(wait=False) caused the child process to
            # hang unrecoverably. This seems to be a bug in Python 3.7
            # It's probably best to gut concurrent.futures entirely, so we're
            # ignoring the resource leakage for the moment.
            return future

        async_wrapper = self._rewrite_wrapper_signature(async_wrapper)
        self._set_wrapper_properties(async_wrapper)
        self._set_wrapper_name(async_wrapper, 'asynchronous')
        return async_wrapper

    def _bind_parsl(self, ctx, *args, **kwargs):
        # If you find a good way to determine if a parsl config is loaded.
        # Use it here
        try:
            parsl.load(get_parsl_config())
        except RuntimeError:
            pass

        futures = []
        remapped_args = []
        for arg in args[1:]:
            if isinstance(arg, ProxyArtifact):
                futures.append(arg.future)
                remapped_args.append(ProxyArtifact(len(futures) - 1,
                                     arg.selector))
            else:
                remapped_args.append(arg)

        remapped_kwargs = {}
        for key, value in kwargs.items():
            if isinstance(value, ProxyArtifact):
                futures.append(value.future)
                remapped_kwargs[key] = ProxyArtifact(len(futures) - 1,
                                                     value.selector)
            else:
                remapped_kwargs[key] = value

        if self.id in ctx.action_executor_mapping:
            executor = ctx.action_executor_mapping[self.id]
        else:
            executor = 'default'

        # TODO: CREATE ISSUE IN PARSL ABOUT PYTHON_APP(JOIN=TRUE) SELECTING
        # EXECUTOR THAT IS NOT LOCAL TO THE MAIN PROCESS (EX: HTEX) BLOWING UP
        if isinstance(self, qiime2.sdk.action.Pipeline):
            # NOTE: Do not make this a python_app(join=True). We need it to run
            # in the parsl main thread
            future = join_app()(
                    run_parsl_action)(self, ctx, remapped_args,
                                      remapped_kwargs, inputs=futures)
        else:
            future = python_app(
                executors=[executor])(
                    run_parsl_action)(self, ctx, remapped_args,
                                      remapped_kwargs, inputs=futures)

        return qiime2.sdk.util.ProxyResults(future, self.signature.outputs)

    def _get_parsl_wrapper(self):
        def parsl_wrapper(*args, **kwargs):
            return self._bind_parsl(qiime2.sdk.Context(parsl=True), *args,
                                    **kwargs)

        parsl_wrapper = self._rewrite_wrapper_signature(parsl_wrapper)
        self._set_wrapper_properties(parsl_wrapper)
        self._set_wrapper_name(parsl_wrapper, 'parsl')
        return parsl_wrapper

    def _rewrite_wrapper_signature(self, wrapper):
        # Convert the callable's signature into the wrapper's signature and set
        # it on the wrapper.
        return decorator.decorator(
            wrapper, self._callable_sig_converter_(self._callable))

    def _set_wrapper_name(self, wrapper, name):
        wrapper.__name__ = wrapper.__qualname__ = name

    def _set_wrapper_properties(self, wrapper):
        wrapper.__module__ = self.get_import_path(include_self=False)
        wrapper.__doc__ = self._build_numpydoc()
        wrapper.__annotations__ = self._build_annotations()
        # This is necessary so that `inspect` doesn't display the wrapped
        # function's annotations (the annotations apply to the "view API" and
        # not the "artifact API").
        del wrapper.__wrapped__

    def _build_annotations(self):
        annotations = {}
        for name, spec in self.signature.signature_order.items():
            annotations[name] = spec.qiime_type

        output = []
        for spec in self.signature.outputs.values():
            output.append(spec.qiime_type)
        output = tuple(output)

        annotations["return"] = output

        return annotations

    def _build_numpydoc(self):
        numpydoc = []
        numpydoc.append(textwrap.fill(self.name, width=75))
        if self.deprecated:
            base_msg = textwrap.indent(
                textwrap.fill(self._build_deprecation_message(), width=72),
                '   ')
            numpydoc.append('.. deprecated::\n' + base_msg)
        numpydoc.append(textwrap.fill(self.description, width=75))

        sig = self.signature
        parameters = self._build_section("Parameters", sig.signature_order)
        returns = self._build_section("Returns", sig.outputs)

        # TODO: include Usage-rendered examples here

        for section in (parameters, returns):
            if section:
                numpydoc.append(section)

        return '\n\n'.join(numpydoc) + '\n'

    def _build_section(self, header, iterable):
        section = []

        if iterable:
            section.append(header)
            section.append('-'*len(header))
            for key, value in iterable.items():
                variable_line = (
                    "{item} : {type}".format(item=key, type=value.qiime_type))
                if value.has_default():
                    variable_line += ", optional"
                section.append(variable_line)
                if value.has_description():
                    section.append(textwrap.indent(textwrap.fill(
                        str(value.description), width=71), '    '))

        return '\n'.join(section).strip()

    def _build_deprecation_message(self):
        return (f'This {self.type.title()} is deprecated and will be removed '
                'in a future version of this plugin.')


class Method(Action):
    """QIIME 2 Method"""

    type = 'method'

    # Abstract method implementations:

    def _callable_sig_converter_(self, callable):
        # No conversion necessary.
        return callable

    def _callable_executor_(self, scope, view_args, output_types, provenance):
        output_views = self._callable(**view_args)
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
        for output_view, (name, spec) in zip(output_views,
                                             output_types.items()):
            if type(output_view) is not spec.view_type:
                raise TypeError(
                    "Expected output view type %r, received %r" %
                    (spec.view_type.__name__, type(output_view).__name__))

            prov = provenance.fork(name)
            scope.add_reference(prov)

            artifact = qiime2.sdk.Artifact._from_view(
                spec.qiime_type, output_view, spec.view_type, prov)
            artifact = scope.add_parent_reference(artifact)

            output_artifacts.append(artifact)

        results = Results(self.signature.outputs.keys(), output_artifacts)

        if len(results) != len(self.signature.outputs):
            raise ValueError(
                "Number of callable outputs must match number of "
                "outputs defined in signature: %d != %d" %
                (len(results), len(self.signature.outputs)))

        return results

    @classmethod
    def _init(cls, callable, inputs, parameters, outputs, plugin_id, name,
              description, input_descriptions, parameter_descriptions,
              output_descriptions, citations, deprecated, examples):
        signature = qtype.MethodSignature(callable, inputs, parameters,
                                          outputs, input_descriptions,
                                          parameter_descriptions,
                                          output_descriptions)
        return super()._init(callable, signature, plugin_id, name, description,
                             citations, deprecated, examples)


class Visualizer(Action):
    """QIIME 2 Visualizer"""

    type = 'visualizer'

    # Abstract method implementations:

    def _callable_sig_converter_(self, callable):
        return DropFirstParameter.from_function(callable)

    def _callable_executor_(self, scope, view_args, output_types, provenance):
        # TODO use qiime2.plugin.OutPath when it exists, and update visualizers
        # to work with OutPath instead of str. Visualization._from_data_dir
        # will also need to be updated to support OutPath instead of str.
        with tempfile.TemporaryDirectory(prefix='qiime2-temp-') as temp_dir:
            ret_val = self._callable(output_dir=temp_dir, **view_args)
            if ret_val is not None:
                raise TypeError(
                    "Visualizer %r should not return anything. "
                    "Received %r as a return value." % (self, ret_val))
            provenance.output_name = 'visualization'
            viz = qiime2.sdk.Visualization._from_data_dir(temp_dir,
                                                          provenance)
            viz = scope.add_parent_reference(viz)

            results = Results(self.signature.outputs.keys(), (viz,))

            if len(results) != len(self.signature.outputs):
                raise ValueError(
                    "Number of callable outputs must match number of "
                    "outputs defined in signature: %d != %d" %
                    (len(results), len(self.signature.outputs)))

            return results

    @classmethod
    def _init(cls, callable, inputs, parameters, plugin_id, name, description,
              input_descriptions, parameter_descriptions, citations,
              deprecated, examples):
        signature = qtype.VisualizerSignature(callable, inputs, parameters,
                                              input_descriptions,
                                              parameter_descriptions)
        return super()._init(callable, signature, plugin_id, name, description,
                             citations, deprecated, examples)


class Pipeline(Action):
    """QIIME 2 Pipeline"""
    type = 'pipeline'
    _ProvCaptureCls = archive.PipelineProvenanceCapture

    def _callable_sig_converter_(self, callable):
        return DropFirstParameter.from_function(callable)

    def _parsl_callable_executor_(self, scope, view_args, output_types,
                                  provenance):
        outputs = self._callable(scope.ctx, **view_args)
        outputs = tuple(output.get_element(output.future.result())
                        for output in tuplize(outputs))

        for output in outputs:
            if not isinstance(output, qiime2.sdk.Result):
                raise TypeError("Pipelines must return `Result` objects, "
                                "not %s" % (type(output), ))

        # This condition *is* tested by the caller of _callable_executor_, but
        # the kinds of errors a plugin developer see will make more sense if
        # this check happens before the subtype check. Otherwise forgetting an
        # output would more likely error as a wrong type, which while correct,
        # isn't root of the problem.
        if len(outputs) != len(output_types):
            raise TypeError(
                "Number of outputs must match number of output "
                "semantic types: %d != %d"
                % (len(outputs), len(output_types)))

        results = []
        for output, (name, spec) in zip(outputs, output_types.items()):
            if not (output.type <= spec.qiime_type):
                raise TypeError(
                    "Expected output type %r, received %r" %
                    (spec.qiime_type, output.type))
            prov = provenance.fork(name, output)
            scope.add_reference(prov)

            aliased_result = output._alias(prov)
            scope.add_parent_reference(aliased_result)

            results.append(aliased_result)

        if len(results) != len(self.signature.outputs):
            raise ValueError(
                "Number of callable outputs must match number of "
                "outputs defined in signature: %d != %d" %
                (len(results), len(self.signature.outputs)))

        results = Results(self.signature.outputs.keys(), tuple(results))
        return _create_future(results)

    def _callable_executor_(self, scope, view_args, output_types, provenance):
        outputs = self._callable(scope.ctx, **view_args)
        outputs = tuplize(outputs)

        for output in outputs:
            if not isinstance(output, qiime2.sdk.Result):
                raise TypeError("Pipelines must return `Result` objects, "
                                "not %s" % (type(output), ))

        # This condition *is* tested by the caller of _callable_executor_, but
        # the kinds of errors a plugin developer see will make more sense if
        # this check happens before the subtype check. Otherwise forgetting an
        # output would more likely error as a wrong type, which while correct,
        # isn't root of the problem.
        if len(outputs) != len(output_types):
            raise TypeError(
                "Number of outputs must match number of output "
                "semantic types: %d != %d"
                % (len(outputs), len(output_types)))

        results = []
        for output, (name, spec) in zip(outputs, output_types.items()):
            if not (output.type <= spec.qiime_type):
                raise TypeError(
                    "Expected output type %r, received %r" %
                    (spec.qiime_type, output.type))
            prov = provenance.fork(name, output)
            scope.add_reference(prov)

            aliased_result = output._alias(prov)
            aliased_result = scope.add_parent_reference(aliased_result)

            results.append(aliased_result)

        if len(results) != len(self.signature.outputs):
            raise ValueError(
                "Number of callable outputs must match number of "
                "outputs defined in signature: %d != %d" %
                (len(results), len(self.signature.outputs)))

        results = Results(self.signature.outputs.keys(), tuple(results))
        return results

    @classmethod
    def _init(cls, callable, inputs, parameters, outputs, plugin_id, name,
              description, input_descriptions, parameter_descriptions,
              output_descriptions, citations, deprecated, examples):
        signature = qtype.PipelineSignature(callable, inputs, parameters,
                                            outputs, input_descriptions,
                                            parameter_descriptions,
                                            output_descriptions)
        return super()._init(callable, signature, plugin_id, name, description,
                             citations, deprecated, examples)


@python_app
def _create_future(results):
    return results


markdown_source_template = """
```python
%(source)s
```
"""

# TODO add unit test for callables raising this
backend_error_template = """
Your current matplotlib backend (MacOSX) does not work with asynchronous calls.
A recommended backend is Agg, and can be changed by modifying your
matplotlibrc "backend" parameter, which can be found at: \n\n %s
"""
