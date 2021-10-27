# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from dataclasses import dataclass
import re

from qiime2 import sdk
from qiime2.core.type import is_semantic_type, is_visualization_type


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

    def validate(self, inputs: 'UsageInputs', outputs: 'UsageOutputNames'):
        if not isinstance(inputs, UsageInputs):
            raise TypeError('Must provide an instance of UsageInputs.')
        if not isinstance(outputs, UsageOutputNames):
            raise TypeError('Must provide an instance of UsageOutputNames.')

        action_f = self.get_action()

        inputs.validate(action_f.signature)
        outputs.validate(action_f.signature)


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

    def validate(self, signature):
        provided = set(self.keys())
        inputs, params = signature.inputs, signature.parameters

        exp_inputs, optional_inputs = set(), set()
        for name, sig in inputs.items():
            if sig.has_default():
                optional_inputs.add(name)
            else:
                exp_inputs.add(name)

        exp_params, optional_params = set(), set()
        for name, sig in params.items():
            if sig.has_default():
                optional_params.add(name)
            else:
                exp_params.add(name)

        missing = exp_inputs - provided
        if len(missing) > 0:
            raise ValueError('Missing input(s): %r' % (missing, ))

        missing = exp_params - provided
        if len(missing) > 0:
            raise ValueError('Missing parameter(s): %r' % (missing, ))

        all_vals = exp_inputs | optional_inputs | exp_params | optional_params
        extra = provided - all_vals
        if len(extra) > 0:
            raise ValueError('Extra input(s) or parameter(s): %r' %
                             (extra, ))

    def map_variables(self, function):
        result = {}

        def mapped(v):
            if isinstance(v, UsageVariable):
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

    def validate(self, signature):
        provided = set(self.values.keys())
        exp_outputs = set(signature.outputs)

        missing = exp_outputs - provided
        if len(missing) > 0:
            raise ValueError('Missing output(s): %r' % (missing, ))

        extra = provided - exp_outputs
        if len(extra) > 0:
            raise ValueError('Extra output(s): %r' % (extra, ))


class UsageOutputs(sdk.Results):
    pass


class UsageVariable:
    DEFERRED = object()
    VAR_TYPES = (
        'artifact',
        'visualization',
        'metadata',
        'column',
        'file',
    )

    def __init__(self, name: str, factory: callable, var_type: str, usage):
        if not name.isidentifier():
            raise ValueError('invalid identifier: %r' % (name,))

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
        if self.value is self.DEFERRED:
            self.value = self.factory()
        return self.value

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

    def variable_factory(self, name, factory, var_type):
        return UsageVariable(
            name,
            factory,
            var_type,
            self,
        )

    def init_artifact(self, name, factory):
        return self.variable_factory(name, factory, 'artifact')

    def init_metadata(self, name, factory):
        return self.variable_factory(name, factory, 'metadata')

    def import_data(*args, **kwargs):
        # TODO
        ...

    def import_file(*args, **kwargs):
        # TODO
        ...

    def merge_metadata(self, name, *variables):
        if len(variables) < 2:
            raise ValueError('Must provide two or more Metadata inputs.')

        for variable in variables:
            if variable.var_type != 'metadata':
                raise ValueError('Incorrect input type for %s, need '
                                 '`metadata` got %s' %
                                 (variable.name, variable.var_type))

        def factory():
            mds = [v.execute() for v in variables]
            return mds[0].merge(*mds[1:])

        return self.variable_factory(name, factory, 'metadata')

    def get_metadata_column(self, name, column_name, variable):
        def factory():
            return variable.execute().get_column(column_name)

        return self.variable_factory(name, factory, 'column')

    def comment(self, text: str):
        pass

    def action(self, action: UsageAction, inputs: UsageInputs,
               outputs: UsageOutputNames):
        if not isinstance(action, UsageAction):
            raise TypeError('Must provide an instance of UsageAction.')
        action.validate(inputs, outputs)
        action_f = action.get_action()

        def memoized_action(_result=[]):
            if len(_result) == 1:
                return _result[0]

            execed_inputs = inputs.map_variables(lambda v: v.execute())

            if self.asynchronous:
                results = action_f.asynchronous(**execed_inputs).result()
            else:
                results = action_f(**execed_inputs)

            _result.append(results)
            return results

        usage_results = []
        for idx, (param_name, var_name) in enumerate(outputs.items()):
            qiime_type = action_f.signature.outputs[param_name].qiime_type
            if is_visualization_type(qiime_type):
                var_type = 'visualization'
            elif is_semantic_type(qiime_type):
                var_type = 'artifact'
            else:
                raise

            usage_results.append(
                self.variable_factory(
                    var_name,
                    lambda i=idx: memoized_action()[i],
                    var_type,
                )
            )

        return UsageOutputs(outputs.keys(), usage_results)


class DiagnosticUsage(Usage):
    def __init__(self):
        super().__init__()
        self.records = []

        @dataclass(frozen=True)
        class DiagnosticUsageRecord:
            source: str
            variable: UsageVariable

        self.record_class = DiagnosticUsageRecord

    def append_record(self, source, variable):
        self.records.append(
            self.record_class(source, variable)
        )

    def init_artifact(self, name, factory):
        variable = super().init_artifact(name, factory)
        self.append_record('init_artifact', variable)
        return variable

    def init_metadata(self, name, factory):
        variable = super().init_metadata(name, factory)
        self.append_record('init_metadata', variable)
        return variable

    def merge_metadata(self, name, *variables):
        variable = super().merge_metadata(name, *variables)
        self.append_record('merge_metadata', variable)
        return variable

    def get_metadata_column(self, name, column_name, variable):
        variable = super().get_metadata_column(name, column_name, variable)
        self.append_record('get_metadata_column', variable)
        return variable

    def comment(self, text):
        self.append_record('comment', text)

    def action(self, action, input_opts, output_opts):
        variables = super().action(action, input_opts, output_opts)
        self.append_record('action', variables)
        return variables


class ExecutionUsageVariable(UsageVariable):
    def __repr__(self):
        return 'ExecutionUsageVariable<name=%r, var_type=%r>' % \
                (self.name, self.var_type)

    def assert_has_line_matching(self, path, expression):
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
        self.records = dict()

    def variable_factory(self, name, factory, var_type):
        return ExecutionUsageVariable(
            name,
            factory,
            var_type,
            self,
        )

    def init_artifact(self, name, factory):
        variable = super().init_artifact(name, factory)

        variable.execute()
        self.records[variable.name] = variable

        return variable

    def init_metadata(self, name, factory):
        variable = super().init_metadata(name, factory)

        variable.execute()
        self.records[variable.name] = variable

        return variable

    def merge_metadata(self, name, *variables):
        variable = super().merge_metadata(name, *variables)

        variable.execute()
        self.records[variable.name] = variable

        return variable

    def get_metadata_column(self, name, column_name, variable):
        variable = super().get_metadata_column(name, column_name, variable)

        variable.execute()
        self.records[variable.name] = variable

        return variable

    def action(self, action, input_opts, output_opts):
        variables = super().action(action, input_opts, output_opts)

        for variable in variables:
            variable.execute()
            self.records[variable.name] = variable

        return variables
