# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import abc
import concurrent.futures
import inspect
import os.path
import tempfile
import textwrap

import decorator
import frontmatter
import ipymd

import qiime.sdk
import qiime.core.type as qtype
from qiime.core.util import LateBindingAttribute, DropFirstParameter, tuplize


# This class provides an interface similar to a `namedtuple` type. We can't use
# `namedtuple` directly because the `Results` type must be defined dynamically
# on a per-instance basis of `Callable` (`Callable.signature` determines the
# structure of the `namedtuple`). This prevents the dynamic `namedtuple` types
# from being pickled (necessary for `async`) because `Results` must be
# accessible as a module global, but this global type will redefined each time
# a `Callable` is instantiated.
class Results(tuple):
    """Tuple class representing the named results of a `Callable`.

    Instances of this class may be returned by callables that have multiple
    return values, and are compatible with how Python handles multiple return
    values from callables. Provides an interface similar to a `namedtuple`
    type (e.g. fields are accessible as attributes).

    Users should not instantiate this class directly.

    """

    # Subclassing `tuple` requires `__new__` override.
    def __new__(cls, fields, values):
        fields = tuple(fields)
        values = tuple(values)

        if len(fields) != len(values):
            raise ValueError(
                "`fields` and `values` must have matching length: %d != %d" %
                (len(fields), len(values)))

        # Create tuple instance, store fields, and create read-only attributes
        # for each field name. Fields must be stored for pickling/copying (see
        # `__getnewargs__`).
        #
        # Note: setting field names as attributes allows for tab-completion in
        # interactive contexts! Using `__getattr__` does not support this.
        self = super().__new__(cls, values)

        # Must set attributes this way because `__setattr__` prevents
        # setting directly (necessary for immutability).
        object.__setattr__(self, '_fields', fields)

        # Attach field names as instance attributes.
        for field, value in zip(fields, values):
            object.__setattr__(self, field, value)

        return self

    def __getnewargs__(self):
        """Arguments to pass to `__new__`. Used by copy and pickle."""
        # `tuple(self)` returns `values`.
        return self._fields, tuple(self)

    # `__setattr__` and `__delattr__` must be defined to prevent users from
    # creating or deleting attributes after this class has been instantiated.
    # `tuple` and `namedtuple` do not have this problem because they are
    # immutable (`__slots__ = ()`). We cannot make this class immutable because
    # we cannot define nonempty `__slots__` when subclassing `tuple`, and we
    # need the `_fields` attribute. We work around this issue by disallowing
    # setting and deleting attributes. The error messages here match those
    # raised by `namedtuple` in Python 3.5.1.
    def __setattr__(self, name, value):
        raise AttributeError("can't set attribute")

    def __delattr__(self, name):
        raise AttributeError("can't delete attribute")

    def __eq__(self, other):
        # Results with different field names should not compare equal, even if
        # their values are equal.
        return (
            isinstance(other, Results) and
            self._fields == other._fields and
            tuple(self) == tuple(other)
        )

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):
        # It is possible to provide an evalable repr but this type of repr does
        # not make the field/value pairs apparent. If the constructor accepted
        # **kwargs, the order of field/value pairs would be lost.
        lines = []
        lines.append('%s (name = value)' % self.__class__.__name__)
        lines.append('')

        max_len = -1
        for field in self._fields:
            if len(field) > max_len:
                max_len = len(field)

        for field, value in zip(self._fields, self):
            field_padding = ' ' * (max_len - len(field))
            lines.append('%s%s = %r' % (field, field_padding, value))

        max_len = -1
        for line in lines:
            if len(line) > max_len:
                max_len = len(line)
        lines[1] = '-' * max_len

        return '\n'.join(lines)


class Callable(metaclass=abc.ABCMeta):
    action_type = 'callable'
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

    @classmethod
    def from_function(cls, function, signature, package):
        id = function.__name__
        return cls(package, id, signature, function, ('function', function))

    def __init__(self, package, id, signature, callable, callable_source,
                 pid=None):
        """

        Parameters
        ----------
        package : str
        id : str
        signature : qiime.core.type.Signature
        callable : callable
        callable_source : (str, object)
        pid : int, optional

        """
        self.package = package
        self.id = id
        self.signature = signature
        self._callable = callable
        self._callable_source = callable_source

        if pid is None:
            pid = os.getpid()
        self._pid = pid

        self._dynamic_call = self._get_callable_wrapper()
        self._dynamic_async = self._get_async_wrapper()

    def _get_callable_wrapper(self):
        def callable_wrapper(*args, **kwargs):
            # This function's signature is rewritten below using
            # `decorator.decorator`. When the signature is rewritten, args[0]
            # is the function whose signature was used to rewrite this
            # function's signature.
            args = args[1:]

            # TODO this may be able to be simplified once Signature respects
            # order.
            wrapper_sig = self._callable_sig_converter_(self._callable)
            wrapper_sig = inspect.Signature.from_callable(wrapper_sig)
            wrapper_params = wrapper_sig.parameters

            user_input = {name: value for value, name in
                          zip(args, wrapper_params)}
            user_input.update(kwargs)

            self.signature.check_types(**user_input)
            output_types = self.signature.solve_output(**user_input)

            artifacts = {}
            for name in self.signature.inputs:
                artifacts[name] = user_input[name]

            parameters = {}
            for name in self.signature.parameters:
                parameters[name] = user_input[name]

            view_args = parameters.copy()
            for name, (_, view_type) in self.signature.inputs.items():
                view_args[name] = artifacts[name].view(view_type)

            outputs = self._callable_executor_(self._callable, view_args,
                                               output_types)
            # `outputs` matches a Python function's return: either a single
            # value is returned, or it is a tuple of return values. Treat both
            # cases uniformly.
            outputs_tuple = tuplize(outputs)
            for output in outputs_tuple:
                output._orphan(self._pid)

            if len(outputs_tuple) != len(self.signature.outputs):
                raise ValueError(
                    "Number of callable outputs must match number of outputs "
                    "defined in signature: %d != %d" %
                    (len(outputs_tuple), len(self.signature.outputs)))

            # Wrap in a Results object mapping output name to value so users
            # have access to outputs by name or position.
            return Results(self.signature.outputs.keys(), outputs_tuple)

        callable_wrapper = self._rewrite_wrapper_signature(callable_wrapper)
        self._set_wrapper_properties(callable_wrapper, '__call__')
        return callable_wrapper

    def _get_async_wrapper(self):
        def async_wrapper(*args, **kwargs):
            # TODO handle this better in the future, but stop the massive error
            # caused by MacOSX async runs for now.
            try:
                import matplotlib as plt
                if plt.rcParams['backend'] == 'MacOSX':
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
            future = pool.submit(self, *args, **kwargs)
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
        del wrapper.__annotations__
        # This is necessary so that `inspect` doesn't display the wrapped
        # function's annotations (the annotations apply to the "view API" and
        # not the "artifact API").
        del wrapper.__wrapped__

    def get_import_path(self):
        return self.package + '.' + self.id

    def __repr__(self):
        return "<%s %s>" % (self.action_type, self.get_import_path())

    def __getstate__(self):
        return {
            'type': self._callable_source[0],
            'src': self._callable_source[1],
            'attrs': {
                'id': self.id,
                'signature': self.signature,
                'package': self.package,
                'pid': self._pid
            }
        }

    def __setstate__(self, state):
        if state['type'] == 'function':
            self.__init__(callable=state['src'],
                          callable_source=(state['type'], state['src']),
                          **state['attrs'])
        else:
            raise NotImplementedError


class MethodCallable(Callable):
    action_type = 'method'
    # Abstract method implementations:

    def _callable_sig_converter_(self, callable):
        # No conversion necessary.
        return callable

    def _callable_executor_(self, callable, view_args, output_types):
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
        for output_view, (semantic_type, view_type) in \
                zip(output_views, output_types.values()):
            if type(output_view) is not view_type:
                raise TypeError(
                    "Expected output view type %r, received %r" %
                    (view_type.__name__, type(output_view).__name__))
            artifact = qiime.sdk.Artifact._from_view(
                semantic_type, output_view, view_type)
            output_artifacts.append(artifact)

        if len(output_artifacts) == 1:
            return output_artifacts[0]
        else:
            return tuple(output_artifacts)

    @classmethod
    def from_function(cls, function, inputs, parameters, outputs, package):
        signature = qtype.MethodSignature.from_function(
            function, inputs, parameters, outputs)
        return super().from_function(function, signature, package)

    @classmethod
    def from_markdown(cls, markdown_filepath, package):
        self = cls.__new__(cls)
        self._from_markdown(markdown_filepath, package)
        return self

    def _from_markdown(self, markdown_filepath, package, pid=None):
        signature = qtype.MethodSignature.from_markdown(markdown_filepath)

        with open(markdown_filepath) as fh:
            _, template = frontmatter.parse(fh.read())

        # TODO: verify that `id` is a valid Python identifier
        id = os.path.splitext(os.path.basename(markdown_filepath))[0]

        # TODO handle default values for optional parameters when that's
        # supported
        function_def_line = 'def %s(%s, %s):' % (
            id, ', '.join(signature.inputs), ', '.join(signature.parameters))
        indent = ' ' * 4
        function_body = ipymd.convert(template, from_='markdown', to='python')
        function_body = textwrap.indent(function_body, indent)
        function_return_line = '%sreturn %s' % (
            indent, ', '.join(signature.outputs))

        function_str = '\n'.join([function_def_line,
                                  function_body,
                                  function_return_line])

        scope = {}
        exec(function_str, scope)
        function = scope[id]

        self.__init__(package, id, signature, function,
                      ('markdown', markdown_filepath), pid=pid)

    def __setstate__(self, state):
        try:
            super().__setstate__(state)
        except NotImplementedError:
            if state['type'] == 'markdown':
                self._from_markdown(state['src'], state['attrs']['package'],
                                    state['attrs']['pid'])
            else:
                raise NotImplementedError


class VisualizerCallable(Callable):
    action_type = 'visualizer'
    # Abstract method implementations:

    def _callable_sig_converter_(self, callable):
        return DropFirstParameter.from_function(callable)

    def _callable_executor_(self, callable, view_args, output_types):
        # TODO use qiime.plugin.OutPath when it exists, and update visualizers
        # to work with OutPath instead of str. Visualization._from_data_dir
        # will also need to be updated to support OutPath instead of str.
        with tempfile.TemporaryDirectory(prefix='qiime2-temp-') as temp_dir:
            ret_val = callable(output_dir=temp_dir, **view_args)
            if ret_val is not None:
                raise TypeError(
                    "Visualizer %r should not return anything. "
                    "Received %r as a return value." % (self, ret_val))
            return qiime.sdk.Visualization._from_data_dir(temp_dir)

    @classmethod
    def from_function(cls, function, inputs, parameters, package):
        signature = qtype.VisualizerSignature.from_function(
            function, inputs, parameters)
        return super().from_function(function, signature, package)


# TODO add unit test for callables raising this
backend_error_template = """
Your current matplotlib backend (MacOSX) does not work with async calls.
A recommended backend is Agg, and can be changed by modifying your
matplotlibrc "backend" parameter, which can be found at: \n\n %s
"""
