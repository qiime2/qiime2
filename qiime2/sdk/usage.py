# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import abc
import os
import re
import types
import tempfile
import zipfile

from qiime2 import sdk

# TODO: docstrings


class UsageAction:
    def __init__(self, *, plugin_id: str, action_id: str):
        if plugin_id == '':
            raise ValueError('Must specify a value for plugin_id.')

        if action_id == '':
            raise ValueError('Must specify a value for action_id.')

        self.plugin_id = plugin_id
        self.action_id = action_id
        self._plugin_manager = sdk.PluginManager()

    def __repr__(self):
        return 'UsageAction(plugin_id=%r, action_id=%r)' %\
            (self.plugin_id, self.action_id)

    def get_action(self):
        plugin = self._plugin_manager.get_plugin(id=self.plugin_id)
        # TODO: should this validation be pushed up into
        # plugin.py or action.py?
        try:
            action_f = plugin.actions[self.action_id]
        except KeyError:
            raise KeyError('No action currently registered with '
                           'id: "%s".' % (self.action_id,))
        return (action_f, action_f.signature)

    def validate(self, inputs, outputs):
        if not isinstance(inputs, UsageInputs):
            raise TypeError('Must provide an instance of UsageInputs.')
        if not isinstance(outputs, UsageOutputNames):
            raise TypeError('Must provide an instance of UsageOutputNames.')

        _, sig = self.get_action()

        inputs.validate(sig)
        outputs.validate(sig)


class UsageInputs:
    def __init__(self, **kwargs):
        self.values = kwargs

    def __repr__(self):
        return 'UsageInputs(**%r)' % (self.values,)

    def validate(self, signature):
        provided = set(self.values.keys())
        inputs, params = signature.inputs, signature.parameters
        exp_inputs = {k for k, v in inputs.items() if not v.has_default()}
        exp_params = {k for k, v in params.items() if not v.has_default()}

        missing = exp_inputs - provided
        if len(missing) > 0:
            raise ValueError('Missing input(s): %r' % (missing, ))

        missing = exp_params - provided
        if len(missing) > 0:
            raise ValueError('Missing parameter(s): %r' % (missing, ))

        extra = provided - exp_inputs - exp_params
        if len(extra) > 0:
            raise ValueError('Extra input(s) or parameter(s): %r' %
                             (extra, ))

    def build_opts(self, signature, scope):
        opts = {}

        for name, signature in signature.signature_order.items():
            if name in self.values:
                if isinstance(self.values[name], ScopeRecord) \
                        and self.values[name].id in scope.records:
                    value = self.values[name].result
                else:
                    value = self.values[name]
                opts[name] = value

        return opts


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

    def get(self, key):
        return self.values[key]

    def validate(self, signature):
        provided = set(self.values.keys())
        exp_outputs = set(signature.outputs)

        missing = exp_outputs - provided
        if len(missing) > 0:
            raise ValueError('Missing output(s): %r' % (missing, ))

        extra = provided - exp_outputs
        if len(extra) > 0:
            raise ValueError('Extra output(s): %r' % (extra, ))

    def validate_derived(self, derived_outputs):
        provided = set(derived_outputs.keys())
        exp_outputs = set(self.values.values())

        missing = exp_outputs - provided
        if len(missing) > 0:
            raise ValueError('SDK implementation is missing output(s): %r' %
                             (missing, ))

        extra = provided - exp_outputs
        if len(extra) > 0:
            raise ValueError('SDK implementation has specified extra '
                             'output(s): %r' % (extra, ))

    def build_opts(self, action_signature, scope):
        opts = {}

        for output in action_signature.outputs.keys():
            opts[self.get(output)] = output

        return opts


class ScopeRecord:
    STATE_PENDING = 'pending'
    STATE_FINISHED = 'finished'

    def __init__(self, id: str, factory: callable = None,
                 value: object = None,
                 assert_has_line_matching: callable = None):
        if factory is not None and not callable(factory):
            # TODO
            raise Exception

        if assert_has_line_matching is not None and \
                not callable(assert_has_line_matching):
            # TODO
            raise Exception

        self.id = id
        if value is None:
            self.state = self.STATE_PENDING
        else:
            self.state = self.STATE_FINISHED
        self._factory = factory
        self._result = value
        self._assert_has_line_matching_ = assert_has_line_matching

    def __repr__(self):
        if self.state == self.STATE_PENDING:
            return 'ScopeRecord<id=%s, state=%s>' % (self.id, self.state)
        else:
            return 'ScopeRecord<id=%s, result=%r>' % (self.id, self.result)

    @property
    def result(self):
        if self.is_pending:
            # TODO
            raise Exception
        return self._result

    @property
    def is_pending(self):
        return self.state == self.STATE_PENDING

    def resolve(self):
        if self.is_pending:
            self._result = self._factory()
            self.state = self.STATE_FINISHED

    def assert_has_line_matching(self, label, path, expression):
        return self._assert_has_line_matching_(self.id, label, path,
                                               expression)


class Scope:
    def __init__(self):
        self._records = dict()

    def __repr__(self):
        return '%r' % (self._records, )

    @property
    def records(self):
        return types.MappingProxyType(self._records)

    def push_record(self, ref, factory=None, value=None,
                    assert_has_line_matching=None):
        record = ScopeRecord(id=ref, factory=factory, value=value,
                             assert_has_line_matching=assert_has_line_matching)
        self._records[ref] = record
        return record

    def get_record(self, id):
        try:
            return self.records[id]
        except KeyError:
            raise KeyError('No record with id: "%s" in scope.' % (id,))


class Usage(metaclass=abc.ABCMeta):
    def __init__(self):
        self._scope = Scope()

    def init_data(self, ref, factory):
        value = self._init_data_(ref, factory)
        return self._push_record(ref=ref, value=value)

    def merge_metadata(self, ref, *records):
        if len(records) < 2:
            raise ValueError('Must provide two or more Metadata inputs.')

        value = self._merge_metadata_(ref, records)
        return self._push_record(ref=ref, value=value)

    def get_result(self, ref):
        return self._get_record(ref)

    def comment(self, text: str):
        return self._comment_(text)

    def action(self, action: UsageAction, inputs: UsageInputs,
               outputs: UsageOutputNames):

        if not isinstance(action, UsageAction):
            raise TypeError('Must provide an instance of UsageAction.')
        action.validate(inputs, outputs)

        _, action_signature = action.get_action()

        input_opts = inputs.build_opts(action_signature, self._scope)
        output_opts = outputs.build_opts(action_signature, self._scope)

        derived_outputs = self._action_(action, input_opts, output_opts)
        self._add_outputs_to_scope(outputs, derived_outputs)

    def _action_(self, action: UsageAction,
                 input_opts: dict, output_opts: list):
        raise NotImplementedError

    def _comment_(self, text: str):
        raise NotImplementedError

    def _assert_has_line_matching_(self, scope_id, label, path, expression):
        raise NotImplementedError

    def _merge_metadata_(self, ref, records):
        raise NotImplementedError

    def _init_data_(self, ref, factory):
        raise NotImplementedError

    def _add_outputs_to_scope(self, outputs, derived_outputs):
        outputs.validate_derived(derived_outputs)
        for ref, result in derived_outputs.items():
            self._push_record(ref, value=result)

    def _push_record(self, ref, factory=None, value=None):
        return self._scope.push_record(
            ref=ref, factory=factory, value=value,
            assert_has_line_matching=self._assert_has_line_matching_)

    def _get_record(self, id):
        return self._scope.get_record(id)

    def _get_records(self):
        return self._scope.records


class DiagnosticUsage(Usage):
    def __init__(self):
        super().__init__()
        self._recorder = []

    def _action_(self, action, input_opts, output_opts):
        self._recorder.append({
            'type': 'action',
            'action': action,
            'input_opts': input_opts,
            'output_opts': output_opts,
        })
        return output_opts

    def _comment_(self, text):
        self._recorder.append({
            'type': 'comment',
            'text': text,
        })

    def _assert_has_line_matching_(self, scope_id, label, path, expression):
        self._recorder.append({
            'type': 'assert_has_line_matching',
            'scope_id': scope_id,
            'label': label,
            'path': path,
            'expression': expression,
        })

    def _init_data_(self, ref, factory):
        return ref

    def _merge_metadata_(self, ref, records):
        return ref


class ExecutionUsage(Usage):
    def _init_data_(self, ref, factory):
        return factory()

    def _comment_(self, text):
        pass

    def _action_(self, action: UsageAction,
                 input_opts: dict, output_opts: dict):
        action_f, _ = action.get_action()
        results = action_f(**input_opts)
        return {k: getattr(results, v) for k, v in output_opts.items()}

    def _assert_has_line_matching_(self, scope_id, label, path, expression):
        data = self._get_record(scope_id).result
        with tempfile.TemporaryDirectory(prefix='q2-exc-usage-') as temp_dir:
            fp = data.save(os.path.join(temp_dir, str(scope_id)))
            with zipfile.ZipFile(fp, 'r') as zip_temp:
                for fn in zip_temp.namelist():
                    if re.match(path, fn):
                        path = fn
                        break  # Will only match on first hit

                with zip_temp.open(path) as file_temp:
                    target = file_temp.read().decode('utf-8')
                    match = re.search(expression, target,
                                      flags=re.MULTILINE)
                    if match is None:
                        raise AssertionError(
                            'Expression %r not found in %s.' %
                            (expression, path))

    def _merge_metadata_(self, ref, records):
        mds = [r.result for r in records]
        return mds[0].merge(*mds[1:])
