# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import dataclasses
import functools
import re

from qiime2 import sdk
from qiime2.core.type import is_semantic_type, is_visualization_type


def assert_usage_var_type(usage_variable, *valid_var_types):
    if usage_variable.var_type not in valid_var_types:
        tmpl = (
            usage_variable.name, valid_var_types, usage_variable.var_type,
        )
        raise ValueError('Incorrect var_type for %s, need %s got %s' % tmpl)


class UsageAction:
    def __init__(self, *, plugin_id: str, action_id: str):
        if plugin_id == '':
            raise ValueError('Must specify a value for plugin_id.')

        if action_id == '':
            raise ValueError('Must specify a value for action_id.')

        self.plugin_id = plugin_id
        self.action_id = action_id
        try:
            self._plugin_manager = sdk.PluginManager.reuse_existing()
        except sdk.UninitializedPluginManagerError:
            raise sdk.UninitializedPluginManagerError(
                'Please create an instance of sdk.PluginManager'
            )

    def __repr__(self):
        return 'UsageAction(plugin_id=%r, action_id=%r)' %\
            (self.plugin_id, self.action_id)

    def get_action(self):
        plugin = self._plugin_manager.get_plugin(id=self.plugin_id)
        try:
            action_f = plugin.actions[self.action_id]
        except KeyError:
            raise KeyError('No action currently registered with '
                           'id: "%s".' % (self.action_id,))
        return action_f


class UsageInputs:
    def __init__(self, **kwargs):
        self.values = kwargs

    def __repr__(self):
        return 'UsageInputs(**%r)' % (self.values,)

    def __getitem__(self, key):
        return self.values[key]

    def __contains__(self, key):
        return key in self.values

    def items(self):
        return self.values.items()

    def keys(self):
        return self.values.keys()

    def values(self):
        return self.values.values()

    def map_variables(self, function):
        result = {}

        def mapped(v):
            if isinstance(v, UsageVariable):
                assert_usage_var_type(v, 'artifact', 'metadata', 'column')
                v = function(v)
            return v

        for name, value in self.items():
            if isinstance(value, (list, set)):
                collection_type = type(value)
                value = [mapped(v) for v in value]
                value = collection_type(value)
            else:
                value = mapped(value)

            result[name] = value

        return result


class UsageOutputNames:
    def __init__(self, **kwargs):
        for key, val in kwargs.items():
            if not isinstance(val, str):
                raise TypeError(
                    'Name provided for key %r must be a string, not a %r.' %
                    (key, type(val)))

        self.values = kwargs

    def __repr__(self):
        return 'UsageOutputNames(**%r)' % (self.values, )

    def __getitem__(self, key):
        return self.values[key]

    def __contains__(self, key):
        return key in self.values

    def items(self):
        return self.values.items()

    def keys(self):
        return self.values.keys()

    def values(self):
        return self.values.values()


class UsageOutputs(sdk.Results):
    pass


class UsageVariable:
    DEFERRED = object()
    VAR_TYPES = (
        'artifact',
        'visualization',
        'metadata',
        'column',
        'format',
    )

    def __init__(self, name: str, factory: callable, var_type: str, usage):
        if not callable(factory):
            raise TypeError('value for `factory` should be a `callable`, '
                            'recieved %s' % (type(factory),))

        if var_type not in self.VAR_TYPES:
            raise ValueError('value for `var_type` should be one of %r, '
                             'received %s' % (self.VAR_TYPES, var_type))

        self.name = name
        self.factory = factory
        self.var_type = var_type
        self.value = self.DEFERRED
        self.use = usage

    def __repr__(self):
        return '%s<name=%r, var_type=%r>' % (self.__class__.__name__,
                                             self.name, self.var_type)

    @property
    def is_deferred(self):
        return self.value is self.DEFERRED

    def execute(self):
        if self.is_deferred:
            self.value = self.factory()
        return self.value

    def save(self, filepath, ext=None):
        value = self.execute()
        return value.save(filepath, ext=ext)

    def to_interface_name(self):
        return self.name

    def assert_has_line_matching(self, path, expression):
        pass

    def assert_output_type(self, semantic_type):
        pass


class Usage:
    # these are here for namespace/import convenience
    UsageAction = UsageAction
    UsageInputs = UsageInputs
    UsageOutputNames = UsageOutputNames
    # NOTE: not exporting UsageOutputs here because example writers shouldn't
    # be instantiating those on their own, anyway.

    def __init__(self, asynchronous=False):
        self.asynchronous = asynchronous
        self.namespace = set()

    def _usage_variable(self, name, factory, var_type):
        variable = self.usage_variable(name, factory, var_type)
        var_name = variable.to_interface_name()
        if var_name in self.namespace:
            raise ValueError(
                '%r namespace collision (%r)' % (variable, var_name))
        self.namespace.add(var_name)
        return variable

    def usage_variable(self, name, factory, var_type):
        return UsageVariable(name, factory, var_type, self)

    def render(self, flush=False):
        raise NotImplementedError

    def init_artifact(self, name, factory):
        return self._usage_variable(name, factory, 'artifact')

    def init_metadata(self, name, factory):
        return self._usage_variable(name, factory, 'metadata')

    def init_format(self, name, factory, ext=None):
        return self._usage_variable(name, factory, 'format')

    def import_from_format(self, name, semantic_type,
                           variable, view_type=None):
        assert_usage_var_type(variable, 'format')

        def factory():
            from qiime2 import Artifact

            fmt = variable.execute()
            artifact = Artifact.import_data(
                semantic_type, str(fmt), view_type=view_type)

            return artifact
        return self._usage_variable(name, factory, 'artifact')

    def merge_metadata(self, name, *variables):
        if len(variables) < 2:
            raise ValueError('Must provide two or more Metadata inputs.')

        for variable in variables:
            assert_usage_var_type(variable, 'metadata')

        def factory():
            mds = [v.execute() for v in variables]
            return mds[0].merge(*mds[1:])
        return self._usage_variable(name, factory, 'metadata')

    def get_metadata_column(self, name, column_name, variable):
        assert_usage_var_type(variable, 'metadata')

        def factory():
            return variable.execute().get_column(column_name)
        return self._usage_variable(name, factory, 'column')

    def view_as_metadata(self, name, variable):
        assert_usage_var_type(variable, 'artifact')

        def factory():
            from qiime2 import Metadata
            return variable.execute().view(Metadata)
        return self._usage_variable(name, factory, 'metadata')

    def peek(self, variable):
        assert_usage_var_type(variable, 'artifact', 'visualization')

    def comment(self, text):
        pass

    def help(self, action: UsageAction):
        pass

    def action(self, action: UsageAction, inputs: UsageInputs,
               outputs: UsageOutputNames):
        if not isinstance(action, UsageAction):
            raise ValueError('Invalid value for `action`: expected %r, '
                             'received %r.' % (UsageAction, type(action)))

        if not isinstance(inputs, UsageInputs):
            raise ValueError('Invalid value for `inputs`: expected %r, '
                             'received %r.' % (UsageInputs, type(inputs)))

        if not isinstance(outputs, UsageOutputNames):
            raise ValueError('Invalid value for `outputs`: expected %r, '
                             'received %r.' % (UsageOutputNames,
                                               type(outputs)))

        action_f = action.get_action()

        @functools.lru_cache(maxsize=None)
        def memoized_action():
            execed_inputs = inputs.map_variables(lambda v: v.execute())
            if self.asynchronous:
                return action_f.asynchronous(**execed_inputs).result()
            return action_f(**execed_inputs)

        usage_results = []
        # outputs will be ordered by the `UsageOutputNames` order, not the
        # signature order - this makes it so that the example writer doesn't
        # need to be explicitly aware of the signature order
        for param_name, var_name in outputs.items():
            qiime_type = action_f.signature.outputs[param_name].qiime_type
            if is_visualization_type(qiime_type):
                var_type = 'visualization'
            elif is_semantic_type(qiime_type):
                var_type = 'artifact'
            else:
                raise ValueError('unknown output type: %r' % (qiime_type,))

            def factory(name=param_name):
                results = memoized_action()
                result = getattr(results, name)
                return result

            variable = self._usage_variable(var_name, factory, var_type)
            usage_results.append(variable)

        results = UsageOutputs(outputs.keys(), usage_results)
        cache_info = memoized_action.cache_info
        cache_clear = memoized_action.cache_clear
        # manually graft on cache operations
        object.__setattr__(results, '_cache_info', cache_info)
        object.__setattr__(results, '_cache_reset', cache_clear)
        return results


class DiagnosticUsage(Usage):
    @dataclasses.dataclass(frozen=True)
    class DiagnosticUsageRecord:
        source: str
        variable: UsageVariable

    def __init__(self):
        super().__init__()
        self.recorder = []

    def _append_record(self, source, variable):
        self.recorder.append(self.DiagnosticUsageRecord(source, variable))

    def init_artifact(self, name, factory):
        variable = super().init_artifact(name, factory)
        self._append_record('init_artifact', variable)
        return variable

    def init_metadata(self, name, factory):
        variable = super().init_metadata(name, factory)
        self._append_record('init_metadata', variable)
        return variable

    def init_format(self, name, factory, ext=None):
        variable = super().init_format(name, factory, ext=ext)
        self._append_record('init_format', variable)
        return variable

    def import_from_format(self, name, semantic_type,
                           fmt_variable, view_type=None):
        variable = super().import_from_format(
            name, semantic_type, fmt_variable, view_type=view_type)
        self._append_record('import_from_format', variable)
        return variable

    def merge_metadata(self, name, *variables):
        variable = super().merge_metadata(name, *variables)
        self._append_record('merge_metadata', variable)
        return variable

    def get_metadata_column(self, name, column_name, variable):
        variable = super().get_metadata_column(name, column_name, variable)
        self._append_record('get_metadata_column', variable)
        return variable

    def view_as_metadata(self, name, artifact_variable):
        variable = super().view_as_metadata(name, artifact_variable)
        self._append_record('view_as_metadata', variable)
        return variable

    def peek(self, variable):
        self._append_record('peek', variable)

    def comment(self, text):
        self._append_record('comment', text)

    def help(self, action):
        self._append_record('help', action)

    def action(self, action, input_opts, output_opts):
        variables = super().action(action, input_opts, output_opts)
        self._append_record('action', variables)
        return variables


class ExecutionUsageVariable(UsageVariable):
    def assert_has_line_matching(self, path, expression):
        assert_usage_var_type(self, 'artifact', 'visualization')

        data = self.value

        hits = sorted(data._archiver.data_dir.glob(path))
        if len(hits) != 1:
            raise ValueError('Value provided for path (%s) did not produce '
                             'exactly one hit: %s' % (path, hits))

        target = hits[0].read_text()
        match = re.search(expression, target, flags=re.MULTILINE)
        if match is None:
            raise AssertionError('Expression %r not found in %s.' %
                                 (expression, path))

    def assert_output_type(self, semantic_type):
        data = self.value

        if str(data.type) != str(semantic_type):
            raise AssertionError("Output %r has type %s, which does not match"
                                 " expected output type of %s"
                                 % (self.name, data.type, semantic_type))


class ExecutionUsage(Usage):
    def __init__(self, asynchronous=False):
        super().__init__(asynchronous)
        # This is here for testing-purposes
        self.recorder = dict()

    def usage_variable(self, name, factory, var_type):
        return ExecutionUsageVariable(name, factory, var_type, self)

    def init_artifact(self, name, factory):
        variable = super().init_artifact(name, factory)

        variable.execute()
        self.recorder[variable.name] = variable

        return variable

    def init_metadata(self, name, factory):
        variable = super().init_metadata(name, factory)

        variable.execute()
        self.recorder[variable.name] = variable

        return variable

    def init_format(self, name, factory, ext=None):
        variable = super().init_format(name, factory, ext=ext)

        variable.execute()
        self.recorder[variable.name] = variable

        return variable

    def import_from_format(self, name, semantic_type,
                           fmt_variable, view_type=None):
        variable = super().import_from_format(
            name, semantic_type, fmt_variable, view_type=view_type)

        variable.execute()
        self.recorder[variable.name] = variable

        return variable

    def merge_metadata(self, name, *variables):
        variable = super().merge_metadata(name, *variables)

        variable.execute()
        self.recorder[variable.name] = variable

        return variable

    def get_metadata_column(self, name, column_name, variable):
        variable = super().get_metadata_column(name, column_name, variable)

        variable.execute()
        self.recorder[variable.name] = variable

        return variable

    def view_as_metadata(self, name, artifact_variable):
        variable = super().view_as_metadata(name, artifact_variable)

        variable.execute()
        self.recorder[variable.name] = variable

        return variable

    def action(self, action, input_opts, output_opts):
        variables = super().action(action, input_opts, output_opts)

        for variable in variables:
            variable.execute()
            self.recorder[variable.name] = variable

        return variables
