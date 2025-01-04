# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import abc
import inspect
import tempfile
import textwrap

import decorator
import dill

import qiime2.sdk
import qiime2.core.type as qtype
import qiime2.core.archive as archive
from qiime2.core.util import (LateBindingAttribute, DropFirstParameter,
                              tuplize, create_collection_name)


def _coerce_pipeline_outputs(ctx, outputs):
    """Ensure all futures are resolved and all collections are of type
       ResultCollection
    """
    coerced_outputs = []

    for output in outputs:
        # Ensure collection outputs are ResultCollection
        if isinstance(output, dict) or isinstance(output, list):
            output = qiime2.sdk.ResultCollection(output)

        # Handle proxy outputs if root
        if ctx._parent is None and output is not None:
            output = output.result()

            # Handle proxies as elements of collections
            if isinstance(output, qiime2.sdk.ResultCollection):
                for key, value in output.items():
                    output[key] = value.result()

        coerced_outputs.append(output)

    return tuple(coerced_outputs)


class Action(metaclass=abc.ABCMeta):
    """QIIME 2 Action"""
    type = 'action'
    _ProvCaptureCls = archive.ActionProvenanceCapture

    __call__ = LateBindingAttribute('_dynamic_call')
    asynchronous = LateBindingAttribute('_dynamic_async')
    parallel = LateBindingAttribute('_dynamic_parsl')

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
    def _callable_executor_(self, ctx, view_args, output_types):
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
        self._dynamic_call = self._callable_action_wrapper()
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

    def _bind(self, context_factory, execution_ctx={'type': 'synchronous'}):
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
            provenance = self._ProvCaptureCls(
                self.type, self.plugin_id, self.id, execution_ctx)

            if self.deprecated:
                with qiime2.core.util.warning() as warn:
                    warn(self._build_deprecation_message(), FutureWarning)

            # Type management
            collated_inputs = self.signature.collate_inputs(
                *args, **kwargs)
            self.signature.check_types(**collated_inputs)
            output_types = self.signature.solve_output(**collated_inputs)
            callable_args = self.signature.coerce_user_input(
                **collated_inputs)

            callable_args = \
                self.signature.transform_and_add_callable_args_to_prov(
                    provenance, **callable_args)

            outputs = self._callable_executor_(
                ctx, callable_args, output_types, provenance)

            if len(outputs) != len(self.signature.outputs):
                raise ValueError(
                    "Number of callable outputs must match number of "
                    "outputs defined in signature: %d != %d" %
                    (len(outputs), len(self.signature.outputs)))

            # Wrap in a Results object mapping output name to value so
            # users have access to outputs by name or position.
            results = qiime2.sdk.Results(
                self.signature.outputs.keys(), outputs)

            return results

        bound_callable = self._rewrite_wrapper_signature(bound_callable)
        self._set_wrapper_properties(bound_callable)
        self._set_wrapper_name(bound_callable, self.id)
        return bound_callable

    def _callable_action_wrapper(self):
        # This is a "root" level invocation (not a nested call within a
        # pipeline), so no special factory is needed.
        callable_wrapper = self._bind(lambda: qiime2.sdk.Context(self))
        self._set_wrapper_name(callable_wrapper, '__call__')
        return callable_wrapper

    def _get_async_wrapper(self):
        async_wrapper = self._rewrite_wrapper_signature(
            qiime2.sdk.AsynchronousContext(self)._dispatch_)
        self._set_wrapper_properties(async_wrapper)
        self._set_wrapper_name(async_wrapper, 'asynchronous')
        return async_wrapper

    def _get_parsl_wrapper(self):
        def parsl_wrapper(*args, **kwargs):
            # TODO: Maybe make this a warning instead?
            if not isinstance(self, Pipeline):
                raise ValueError('Only pipelines may be run in parallel')

            # TODO: could call callable_action here instead and do index check
            return qiime2.sdk.ParallelContext(self)._dispatch_(*args, **kwargs)

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

    def _callable_executor_(self, ctx, view_args, output_types, provenance):
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

        output_artifacts = \
            self.signature.coerce_given_outputs(output_views, output_types,
                                                ctx, provenance)

        return tuple(output_artifacts)

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

    def _callable_executor_(self, ctx, view_args, output_types, provenance):
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
            viz = ctx.add_reference(viz)

            return (viz, )

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

    def _callable_executor_(self, ctx, view_args, output_types, provenance):
        outputs = self._callable(ctx, **view_args)
        # Just make sure we have an iterable even if there was only one output
        outputs = tuplize(outputs)

        # Make sure any collections returned are in the form of
        # ResultCollections and resolve proxies if this is the root pipeline
        outputs = _coerce_pipeline_outputs(ctx, outputs)

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

        message = "Pipelines must return `Result` objects, not %s"
        for output in outputs:
            if isinstance(output, qiime2.sdk.ResultCollection):
                for elem in output.values():
                    if not isinstance(elem, qiime2.sdk.IResult):
                        raise TypeError(message % type(elem))
            elif not isinstance(output, qiime2.sdk.IResult):
                raise TypeError(message % type(output))

        results = []

        # If we don't have a Collection, we should have a Result, if we
        # have neither, or our types just don't match up, something bad
        # happened
        for output, (name, spec) in zip(outputs, output_types.items()):
            if spec.qiime_type.name == 'Collection' and \
                    output.collection in spec.qiime_type:
                size = len(output)
                aliased_output = qiime2.sdk.ResultCollection()

                for idx, (key, value) in enumerate(output.items()):
                    collection_name = create_collection_name(
                        name=name, key=key, idx=idx, size=size)
                    aliased_result = \
                        value._alias(collection_name, provenance, ctx)

                    aliased_output[str(key)] = aliased_result
                results.append(aliased_output)
            elif output.type <= spec.qiime_type:
                aliased_result = output._alias(name, provenance, ctx)

                results.append(aliased_result)
            else:
                _type = output.type if hasattr(output, 'type') \
                    else type(output)
                raise TypeError(
                    "Expected output type %r, received %r" %
                    (spec.qiime_type, _type))

        if len(results) != len(self.signature.outputs):
            raise ValueError(
                "Number of callable outputs must match number of "
                "outputs defined in signature: %d != %d" %
                (len(results), len(self.signature.outputs)))

        return tuple(results)

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


markdown_source_template = """
```python
%(source)s
```
"""
