# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import abc
import re
import types

from qiime2 import sdk, metadata

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

    def build_opts(self, signature, scope):
        opts = {}

        for name, signature in signature.signature_order.items():
            if name in self.values:
                v = self.values[name]
                if isinstance(v, ScopeRecord) and v.ref in scope.records:
                    value = self.values[name].result
                elif (isinstance(v, list)
                      and all(isinstance(n, ScopeRecord) for n in v)):
                    value = [item.result for item in v]
                elif (isinstance(v, set)
                      and all(isinstance(n, ScopeRecord) for n in v)):
                    value = {item.result for item in v}
                else:
                    value = v
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

    def validate_computed(self, computed_outputs):
        provided = set(computed_outputs.keys())
        exp_outputs = set(self.values.keys())

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
            opts[output] = self.get(output)

        return opts


class ScopeRecord:
    def __init__(self, ref: str, value: object, source: str,
                 assert_has_line_matching: callable = None):
        if assert_has_line_matching is not None and \
                not callable(assert_has_line_matching):
            raise TypeError('Value for `assert_has_line_matching` should be a '
                            '`callable`.')

        self.ref = ref
        self._result = value
        self._source = source
        self._assert_has_line_matching_ = assert_has_line_matching

    def __repr__(self):
        return 'ScopeRecord<ref=%s, result=%r, source=%s>' % (self.ref,
                                                              self.result,
                                                              self.source)

    @property
    def result(self):
        return self._result

    @property
    def source(self):
        return self._source

    def assert_has_line_matching(self, label, path, expression):
        return self._assert_has_line_matching_(self.ref, label, path,
                                               expression)


class Scope:
    def __init__(self):
        self._records = dict()

    def __repr__(self):
        return '%r' % (self._records, )

    @property
    def records(self):
        return types.MappingProxyType(self._records)

    def push_record(self, ref, value, source, assert_has_line_matching=None):
        record = ScopeRecord(ref=ref, value=value, source=source,
                             assert_has_line_matching=assert_has_line_matching)
        self._records[ref] = record
        return record

    def get_record(self, ref):
        try:
            return self.records[ref]
        except KeyError:
            raise KeyError('No record with ref id: "%s" in scope.' % (ref, ))


class Usage(metaclass=abc.ABCMeta):
    def __init__(self):
        self._scope = Scope()

    def init_data(self, ref, factory):
        value = self._init_data_(ref, factory)
        return self._push_record(ref, value, 'init_data')

    def _init_data_(self, ref, factory):
        raise NotImplementedError

    def init_metadata(self, ref, factory):
        value = self._init_metadata_(ref, factory)
        return self._push_record(ref, value, 'init_metadata')

    def _init_metadata_(self, ref, factory):
        raise NotImplementedError

    def init_data_collection(self, ref, records):
        if len(records) < 2:
            raise ValueError('Must provide two or more Artifact inputs.')

        value = self._init_data_collection_(ref, records)
        return self._push_record(ref, value, 'init_data_collection')

    def _init_data_collection_(self, ref, records):
        raise NotImplementedError

    def merge_metadata(self, ref, *records):
        if len(records) < 2:
            raise ValueError('Must provide two or more Metadata inputs.')

        value = self._merge_metadata_(ref, records)
        return self._push_record(ref, value, 'merge_metadata')

    def _merge_metadata_(self, ref, records):
        raise NotImplementedError

    def get_metadata_column(self, column_name, record):
        value = self._get_metadata_column_(column_name, record)
        return self._push_record(column_name, value, 'get_metadata_column')

    def _get_metadata_column_(self, column_name, record):
        raise NotImplementedError

    def comment(self, text: str):
        return self._comment_(text)

    def _comment_(self, text: str):
        raise NotImplementedError

    def action(self, action: UsageAction, inputs: UsageInputs,
               outputs: UsageOutputNames):

        if not isinstance(action, UsageAction):
            raise TypeError('Must provide an instance of UsageAction.')
        action.validate(inputs, outputs)

        _, action_signature = action.get_action()

        input_opts = inputs.build_opts(action_signature, self._scope)
        output_opts = outputs.build_opts(action_signature, self._scope)

        computed_outputs = self._action_(action, input_opts, output_opts)
        self._add_outputs_to_scope(outputs, computed_outputs)

    def _action_(self, action: UsageAction,
                 input_opts: dict, output_opts: dict):
        raise NotImplementedError

    def _assert_has_line_matching_(self, ref, label, path, expression):
        raise NotImplementedError

    def get_result(self, ref):
        record = self._get_record(ref)
        source = record.source
        if source != 'action':
            raise TypeError('source == %s but must be "action"' % source)
        return record

    def _add_outputs_to_scope(self, outputs, computed_outputs):
        outputs.validate_computed(computed_outputs)
        for output, result in computed_outputs.items():
            ref = outputs.get(output)
            self._push_record(ref, result, 'action')

    def _push_record(self, ref, value, source):
        return self._scope.push_record(
            ref=ref, value=value, source=source,
            assert_has_line_matching=self._assert_has_line_matching_)

    def _get_record(self, ref):
        return self._scope.get_record(ref)

    def _get_records(self):
        return self._scope.records


class DiagnosticUsage(Usage):
    def __init__(self):
        super().__init__()
        self.recorder = []

    def _init_data_(self, ref, factory):
        self.recorder.append({
            'source': 'init_data',
            'ref': ref,
        })
        return ref

    def _init_metadata_(self, ref, factory):
        self.recorder.append({
            'source': 'init_data',
            'ref': ref,
        })
        return ref

    def _init_data_collection_(self, ref, records):
        self.recorder.append({
            'source': 'init_data_collection',
            'ref': ref,
        })
        return ref

    def _merge_metadata_(self, ref, records):
        self.recorder.append({
            'source': 'merge_metadata',
            'ref': ref,
            'records_refs': [r.ref for r in records],
        })
        return ref

    def _get_metadata_column_(self, column_name, record):
        self.recorder.append({
            'source': 'get_metadata_column',
            'ref': column_name,
            'record_ref': record.ref,
            'column_name': column_name,
        })
        return column_name

    def _comment_(self, text):
        self.recorder.append({
            'source': 'comment',
            'text': text,
        })

    def _action_(self, action, input_opts, output_opts):
        self.recorder.append({
            'source': 'action',
            'action': action,
            'input_opts': input_opts,
            'output_opts': output_opts,
        })
        return output_opts

    def _assert_has_line_matching_(self, ref, label, path, expression):
        self.recorder.append({
            'source': 'assert_has_line_matching',
            'ref': ref,
            'label': label,
            'path': path,
            'expression': expression,
        })


class ExecutionUsage(Usage):
    def _init_data_(self, ref, factory):
        result = factory()
        result_type = type(result)

        if result_type not in (list, set, sdk.Artifact):
            raise ValueError('Factory (%r) returned a %s, expected an '
                             'Artifact.' % (factory, result_type))

        if result_type in (list, set):
            if not all([isinstance(i, sdk.Artifact) for i in result]):
                raise ValueError('Factory (%r) returned a %s where not all '
                                 'elements were Artifacts.' %
                                 (factory, result_type))

        return result

    def _init_metadata_(self, ref, factory):
        result = factory()
        result_type = type(result)

        if not isinstance(result, metadata.Metadata):
            raise TypeError('Factory (%r) returned a %s, but expected '
                            'Metadata.' % (factory, result_type))

        return result

    def _init_data_collection_(self, ref, records):
        for factory in records:
            result = factory()
            result_type = type(result)

            if not isinstance(result_type, sdk.Artifact):
                raise ValueError('Factory (%r) returned a %s, expected an '
                                 'Artifact.' % (factory, result_type))

    def _merge_metadata_(self, ref, records):
        mds = [r.result for r in records]
        return mds[0].merge(*mds[1:])

    def _get_metadata_column_(self, column_name, record):
        return record.result.get_column(column_name)

    def _comment_(self, text):
        pass

    def _action_(self, action: UsageAction,
                 input_opts: dict, output_opts: dict):
        action_f, _ = action.get_action()
        results = action_f(**input_opts)
        return {k: getattr(results, k) for k in output_opts.keys()}

    def _assert_has_line_matching_(self, ref, label, path, expression):
        data = self._get_record(ref).result

        hits = sorted(data._archiver.data_dir.glob(path))
        if len(hits) != 1:
            raise ValueError('Value provided for path (%s) did not produce '
                             'exactly one hit: %s' % (path, hits))

        target = hits[0].read_text()
        match = re.search(expression, target, flags=re.MULTILINE)
        if match is None:
            raise AssertionError('Expression %r not found in %s.' %
                                 (expression, path))
