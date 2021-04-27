# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""
The Usage module enables generating examples for multiple QIIME 2 interfaces.
"""

import abc
import re
import types
import typing

from qiime2 import core, sdk, metadata


class UsageAction:
    """
    Parameters
    ----------
    plugin_id : str
        Plugin ID.
    action_id : str
        Action ID.
    """

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

    def get_action(self) -> typing.Tuple[
            typing.Union['sdk.Method', 'sdk.Pipeline'],
            typing.Type['core.PipelineSignature']]:
        """
        Get this example's action and signature.

        Returns
        -------
        action_f : QIIME 2 Method, Visualizer, or Pipeline
            The plugin action.
        action_f.signature : QIIME 2 Method, Visualizer, or Pipeline Signature
            The method signature for the plugin action.
        """

        plugin = self._plugin_manager.get_plugin(id=self.plugin_id)
        # TODO: should this validation be pushed up into
        # plugin.py or action.py?
        try:
            action_f = plugin.actions[self.action_id]
        except KeyError:
            raise KeyError('No action currently registered with '
                           'id: "%s".' % (self.action_id,))
        return (action_f, action_f.signature)

    def validate(self, inputs: 'UsageInputs',
                 outputs: 'UsageOutputNames') -> None:
        """
        Ensure that all required inputs and outputs are declared and reference
        appropriate scope records, if necessary.

        Parameters
        ----------
        inputs : UsageInputs
        outputs : UsageOutputNames
        """
        if not isinstance(inputs, UsageInputs):
            raise TypeError('Must provide an instance of UsageInputs.')
        if not isinstance(outputs, UsageOutputNames):
            raise TypeError('Must provide an instance of UsageOutputNames.')

        _, sig = self.get_action()

        inputs.validate(sig)
        outputs.validate(sig)


class UsageInputs:
    """
    A convenience object for assisting with construction of input options
    and performing signature validation.

    Parameters
    ----------
    kwargs : Example inputs
        Inputs to be passed in to ``Usage.action``.
    """

    def __init__(self, **kwargs: typing.Any):
        self.values = kwargs

    def __repr__(self):
        return 'UsageInputs(**%r)' % (self.values,)

    def validate(self,
                 signature: typing.Type['core.PipelineSignature']) -> None:
        """
        Ensure that all required inputs are accounted for and that there are no
        unnecessary or extra parameters provided.

        Parameters
        ----------
        signature : QIIME 2 Method, Visualizer, or Pipeline Signature
            The plugin action's signature.

        Raises
        ------
        ValueError
            If there are missing or extra inputs or parameters.
        """
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

    def build_opts(self, signature: typing.Type['core.PipelineSignature'],
                   scope: 'Scope') -> dict:
        """
        Build a dictionary mapping example ``ScopeRecord``s and primitives to
        action inputs and parameters, respectively. Values are derived from
        either an input's ``ScopeRecord`` (``ScopeRecord.value``), or the value
        a keyword argument passed into the ``UsageInputs`` constructor.

        Parameters
        ----------
        signature : QIIME 2 Method, Visualizer, or Pipeline Signature
            The plugin action's signature.
        scope : Scope
            A Usage example's current scope.

        Returns
        -------
        dict
            Input identifiers and their example values.
        """
        opts = {}

        for name, signature in signature.signature_order.items():
            if name in self.values:
                v = self.values[name]
                if isinstance(v, ScopeRecord) and v.ref in scope.records:
                    value = self.values[name].result
                else:
                    value = v
                opts[name] = value

        return opts


class UsageOutputNames:
    """
    Parameters
    ----------
    kwargs : str
        A mapping between the Action's output name and the desired scope record
        identifier in which to "store" the results.
    """

    def __init__(self, **kwargs: str):
        for key, val in kwargs.items():
            if not isinstance(val, str):
                raise TypeError(
                    'Name provided for key %r must be a string, not a %r.' %
                    (key, type(val)))

        self.values = kwargs

    def __repr__(self):
        return 'UsageOutputNames(**%r)' % (self.values, )

    def get(self, key) -> str:
        """
        Get an example output's identifier.

        Returns
        -------
        str
            The identifier for an example output.
        """
        return self.values[key]

    def validate(self,
                 signature: typing.Type['core.PipelineSignature']) -> None:
        """
        Ensure that all required outputs are accounted for and that there are
        no unnecessary or extra outputs provided.

        Parameters
        ----------
        signature
            Action signature.

        Raises
        ------
        ValueError
            If the example has missing or extra outputs as per the action
            signature.
        """
        provided = set(self.values.keys())
        exp_outputs = set(signature.outputs)

        missing = exp_outputs - provided
        if len(missing) > 0:
            raise ValueError('Missing output(s): %r' % (missing, ))

        extra = provided - exp_outputs
        if len(extra) > 0:
            raise ValueError('Extra output(s): %r' % (extra, ))

    def validate_computed(self,
                          computed_outputs: typing.Dict[
                              str,
                              typing.Union['sdk.Artifact', 'sdk.Visualization',
                                           'metadata.Metadata']]
                          ) -> None:
        """
        Check that outputs are still valid after being processed by a Usage
        driver's ``_action_`` method.

        Parameters
        ----------
        computed_outputs : dict of outputs
            Outputs returned by the Usage driver's ``_action_`` method.

        Raises
        ------
        ValueError
            If there are missing or extra outputs as per the action signature.
        """
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

    def build_opts(self,
                   action_signature: typing.Type['core.PipelineSignature'],
                   scope: 'Scope') -> dict:
        """
        Build a dictionary mapping action output identifiers to example output
        value.

        Parameters
        ----------
        action_signature : QIIME 2 Method, Visualizer, or Pipeline Signature
            The plugin action's signature.
        scope : Scope
            A Usage example's current scope.

        Returns
        -------
        dict
            Output identifiers and their example values.
        """
        opts = {}

        for output in action_signature.outputs.keys():
            opts[output] = self.get(output)

        return opts


class ScopeRecord:
    """
    Represents a discrete step in a Usage example.

    Provides information needed by Usage drivers to render Usage examples.

    Parameters
    ----------
    ref : str
        A unique identifier for referring to the record.
    value : Artifact, Visualization, or Metadata
        The value referred to by ``ref``.
    source : str
        The Usage method called to initialize example data.
    assert_has_line_matching : callable
        A function for asserting something about rendered example data.

    Notes
    -----
    ``ScopeRecord`` is an internal implementation and need not be
    instantiated manually.
    """

    def __init__(self,
                 ref: str,
                 value: typing.Union['sdk.Artifact', 'sdk.Visualization',
                                     'metadata.Metadata'],
                 source: str,
                 assert_has_line_matching: typing.Optional[
                     typing.Callable] = None):

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
    def result(self) -> typing.Union[
            'sdk.Artifact', 'sdk.Visualization', 'metadata.Metadata']:
        """
        Artifact, Visualization, or Metadata value referred to by ``self.ref``.
        """
        return self._result

    @property
    def source(self) -> str:
        """
        The origin of a ``ScopeRecord``. Possible values for this property
        are ``init_data``, ``init_metadata``, ``init_data_collection``,
        ``merge_metadata``, ``get_metadata_column``, or ``action``.  This
        information is required by some Usage drivers in order to render
        examples properly.

        See Also
        --------
        Scope.push_record

        q2cli.q2cli.core.usage.CLIRenderer
        """
        return self._source

    def assert_has_line_matching(self, label: str, path: str,
                                 expression: str) -> None:
        """
        A proxy for a Usage driver's reference implementation for verifying
        that the file at ``path`` contains a line matching ``expression``.

        Parameters
        ----------
        label : str
            A label for describing this assertion. Interface drivers may
            choose to omit this information in the rendered output.
        path : str
            Path to example data file.
        expression : str
            A regex pattern to be passed as the first argument to
            ``re.search``.

        Raises
        ______
        AssertionError
            If ``expression`` is not found in ``path``.

        See Also
        --------
        See ``ExecutionUsage`` for an example implementation.
        """
        return self._assert_has_line_matching_(self.ref, label, path,
                                               expression)


class Scope:
    """
    Sequentially track all ``ScopeRecord``s for a Usage example. The Scope
    provides a detailed, start-to-finish accounting of all of the steps
    necessary for translating the Usage example into an interface-specific
    representation.

    Notes
    -----
    ``Scope`` is an internal implementation and need not be instantiated
    manually.
    """

    def __init__(self):
        self._records: typing.Dict[str, 'ScopeRecord'] = dict()

    def __repr__(self):
        return '%r' % (self._records, )

    @property
    def records(self) -> types.MappingProxyType:
        """
        A read-only view of ``ScopeRecords`` in the current scope.
        """
        return types.MappingProxyType(self._records)

    def push_record(self,
                    ref: str,
                    value: typing.Union['sdk.Artifact', 'sdk.Visualization',
                                        'metadata.Metadata'],
                    source: str,
                    assert_has_line_matching: typing.Callable = None,
                    ) -> 'ScopeRecord':
        """
        Appends a new ``ScopeRecord`` to the Usage Example's scope.

        Parameters
        ----------
        ref : str
        value : Artifact, Visualization, or Metadata
            Data from a Usage data initialization method.
        source : str
            The origin of Usage example data.
        assert_has_line_matching : callable
            Verify that the file at ``path`` contains a line matching
            ``expression`` within an Artifact. See
            ``ScopeRecord.assert_has_line_matching``. This is a proxy for a
            Usage driver's reference implementation.

        Returns
        -------
        record : ScopeRecord
        """
        record = ScopeRecord(ref=ref, value=value, source=source,
                             assert_has_line_matching=assert_has_line_matching)
        self._records[ref] = record
        return record

    def get_record(self, ref: str) -> ScopeRecord:
        """
        Look up a ``ScopeRecord`` from the current scope by identifier.

        Parameters
        ----------
        ref : str
            The identifier for a ``ScopeRecord``.

        Raises
        ------
        KeyError
            If the record identifier isn't in the scope.

        Returns
        -------
        record : ScopeRecord
        """
        try:
            return self.records[ref]
        except KeyError:
            raise KeyError('No record with ref id: "%s" in scope.' % (ref, ))


class Usage(metaclass=abc.ABCMeta):
    """
    ``Usage`` is the base class for Usage driver implementations and should be
    customized for the interface or implementation.
    """

    UsageAction = UsageAction
    UsageInputs = UsageInputs
    UsageOutputNames = UsageOutputNames

    def __init__(self):
        self._scope = Scope()

    def init_data(self, ref: str,
                  factory: typing.Callable[[], 'sdk.Artifact']) \
            -> 'ScopeRecord':
        """
        Initialize example Artifact data from a factory. Whether or not the
        example data is actually created is dependent on the driver executing
        the example.

        Parameters
        ----------
        ref : str
            Unique identifier for example data.
        factory : callable
            A factory that returns an example Artifact.

        Returns
        -------
        record : ScopeRecord
            A record with information about example data.
        """
        value = self._init_data_(ref, factory)
        return self._push_record(ref, value, 'init_data')

    def _init_data_(self, ref, factory):
        raise NotImplementedError

    def init_metadata(self, ref: str,
                      factory: typing.Callable[[], 'metadata.Metadata']) \
            -> 'ScopeRecord':
        """
        Initialize Metadata for a Usage example. Whether or not the example
        metadata is actually created is dependent on the driver executing the
        example.

        Parameters
        ----------
        ref : str
            Unique identifier for example metadata.
        factory : callable
            A factory that returns example Metadata.

        Returns
        -------
        record : ScopeRecord
            A record with information about example Metadata.
        """
        value = self._init_metadata_(ref, factory)
        return self._push_record(ref, value, 'init_metadata')

    def _init_metadata_(self, ref, factory):
        raise NotImplementedError

    def init_data_collection(self,
                             ref: str,
                             collection_type: typing.Union[list, set],
                             *records: ScopeRecord) -> 'ScopeRecord':
        """
        Initialize a collection of data for a Usage example.

        Parameters
        ----------
        ref : str
            Unique identifier for example data collection.
        collection_type : list or set
            The type of collection required by an action.
        records : ScopeRecords belonging to the collection
            The records associated with data to be initialized in the
            collection.

        Returns
        -------
        record : ScopeRecord
            A record with information about an example of collection data.
        """
        if len(records) < 1:
            raise ValueError('Must provide at least one ScopeRecord input.')
        for record in records:
            if not isinstance(record, ScopeRecord):
                raise ValueError('Record (%r) returned a %s, expected a '
                                 'ScopeRecord.' % (record, type(record)))

        value = self._init_data_collection_(ref, collection_type, records)
        return self._push_record(ref, value, 'init_data_collection')

    def _init_data_collection_(self, ref, collection_type, records):
        raise NotImplementedError

    def merge_metadata(self, ref: str,
                       *records: typing.List['ScopeRecord']) -> 'ScopeRecord':
        """
        Merge previously initialized example Metadata.

        Parameters
        ----------
        ref : str
            Unique identifier for merged Metadata.
        records : ScopeRecords
            Records for the example Metadata to be merged.

        Returns
        -------
        record : ScopeRecord
            A new record with information about the merged example Metadata.
        """

        if len(records) < 2:
            raise ValueError('Must provide two or more Metadata inputs.')

        value = self._merge_metadata_(ref, records)
        return self._push_record(ref, value, 'merge_metadata')

    def _merge_metadata_(self, ref, records):
        raise NotImplementedError

    def get_metadata_column(self, column_name: str,
                            record: ScopeRecord) -> 'ScopeRecord':
        """
        Get a Metadata column from previously initialized example Metadata.

        Parameters
        ----------
        column_name : str
            The name of a column in example Metadata.
        record : ScopeRecord
            The record associated with example Metadata.

        Returns
        -------
        record : ScopeRecord
            A new scope record for example Metadata column ``column_name``.
        """
        value = self._get_metadata_column_(column_name, record)
        return self._push_record(column_name, value, 'get_metadata_column')

    def _get_metadata_column_(self, column_name, record):
        raise NotImplementedError

    def comment(self, text: str):
        comment = self._comment_(text)
        return self._push_record(str(comment), comment, 'comment')

    def _comment_(self, text: str):
        raise NotImplementedError

    def action(self, action: UsageAction, inputs: UsageInputs,
               outputs: UsageOutputNames) -> None:
        """
        This method is a proxy for invoking a QIIME 2 Action. Whether or not
        the Action is actually invoked is dependent on the driver executing the
        example.

        Parameters
        ----------
        action : UsageAction
            Specifies the Plugin and Action for a Usage example.
        inputs : UsageInputs
            Specifies the inputs to an Action for a Usage example.
        outputs : UsageOutputNames
            Species the outputs of an Action for a Usage example.

        Examples
        --------
        qiime2.core.testing.examples : Usage examples
        """

        if not isinstance(action, UsageAction):
            raise TypeError('Must provide an instance of UsageAction.')
        action.validate(inputs, outputs)

        _, action_signature = action.get_action()

        input_opts = inputs.build_opts(action_signature, self._scope)
        output_opts = outputs.build_opts(action_signature, self._scope)

        computed_outputs = self._action_(action, input_opts, output_opts)
        self._add_outputs_to_scope(outputs, computed_outputs)

    def _action_(self, action: UsageAction, input_opts: dict,
                 output_opts: dict) -> dict:
        raise NotImplementedError

    def _assert_has_line_matching_(self, ref, label, path, expression):
        raise NotImplementedError

    def get_result(self, ref: str) -> 'ScopeRecord':
        """
        Get the record for a Usage example output. This is a convenience
        method used to access records generated after running ``Usage.action``.

        Parameters
        ----------
        ref : str
            Output identifier.

        Raises
        ------
        KeyError
            If ``ref`` is not associated with a record generated by
            ``Usage.action``.

        TypeError
            If the source type is not an Action.
        """
        record = self._get_record(ref)
        source = record.source
        if source != 'action':
            raise TypeError('source == %s but must be "action"' % source)
        return record

    def _add_outputs_to_scope(self, outputs: UsageOutputNames,
                              computed_outputs):
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

    def _destructure_signature(self, action_sig):
        # In the future this could return a more robust spec subset,
        # if necessary.
        def distill_spec(spec):
            return str(spec.qiime_type)

        inputs = {k: distill_spec(v) for k, v in action_sig.inputs.items()}
        outputs = {k: distill_spec(v) for k, v in action_sig.outputs.items()}
        params, mds = {}, {}

        for param_name, spec in action_sig.parameters.items():
            if sdk.util.is_metadata_type(spec.qiime_type):
                mds[param_name] = distill_spec(spec)
            else:
                params[param_name] = distill_spec(spec)

        return {'inputs': inputs, 'params': params,
                'mds': mds, 'outputs': outputs}

    def _destructure_opts(self, signature, input_opts, output_opts):
        inputs, params, mds, outputs = {}, {}, {}, {}

        for opt_name, val in input_opts.items():
            if opt_name in signature['inputs'].keys():
                inputs[opt_name] = (val, signature['inputs'][opt_name])
            elif opt_name in signature['params'].keys():
                params[opt_name] = (val, signature['params'][opt_name])
            elif opt_name in signature['mds'].keys():
                mds[opt_name] = (val, signature['mds'][opt_name])

        for opt_name, val in output_opts.items():
            outputs[opt_name] = (val, signature['outputs'][opt_name])

        return inputs, params, mds, outputs


class DiagnosticUsage(Usage):
    """
    Reference implementation of a Usage driver.

    Used to generate information for testing the Usage API.

    See Also
    --------
    qiime2.sdk.tests.test_usage.TestUsage : Unit tests using this driver.
    """
    def __init__(self):
        super().__init__()

    def _init_data_(self, ref, factory):
        return {
            'source': 'init_data',
            'ref': ref,
        }

    def _init_metadata_(self, ref, factory):
        return {
            'source': 'init_metadata',
            'ref': ref,
        }

    def _init_data_collection_(self, ref, collection_type, records):
        return {
            'source': 'init_data_collection',
            'ref': ref,
        },  collection_type([i.ref for i in records])

    def _merge_metadata_(self, ref, records):
        return {
            'source': 'merge_metadata',
            'ref': ref,
            'records_refs': [r.ref for r in records],
        }

    def _get_metadata_column_(self, column_name, record):
        return {
            'source': 'get_metadata_column',
            'ref': column_name,
            'record_ref': record.ref,
            'column_name': column_name,
        }

    def _comment_(self, text):
        return {
            'source': 'comment',
            'text': text,
        }

    def _action_(self, action, input_opts, output_opts):
        results = dict()
        for output_opt in output_opts.keys():
            results[output_opt] = {
                'source': 'action',
                'plugin_id': action.plugin_id,
                'action_id': action.action_id,
                'input_opts': input_opts,
                'output_opts': output_opts,
                'output_opt': output_opt,
            }
        return results

    def _assert_has_line_matching_(self, ref, label, path, expression):
        return {
            'source': 'assert_has_line_matching',
            'ref': ref,
            'label': label,
            'path': path,
            'expression': expression,
        }


class ExecutionUsage(Usage):
    """
    Execute and test rendered examples.

    See Also
    --------
    qiime2.sdk.tests.test_usage.TestExecutionUsage : Tests using this driver.
    qiime2.plugin.testing.TestPluginBase.execute_examples : Executes examples.
    """
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

    def _init_data_collection_(self, ref, collection_type, records):
        collection = []
        for record in records:
            collection.append(record.result)

        return collection_type(collection)

    def _merge_metadata_(self, ref, records):
        mds = [r.result for r in records]
        return mds[0].merge(*mds[1:])

    def _get_metadata_column_(self, column_name, record):
        return record.result.get_column(column_name)

    def _comment_(self, text):
        pass

    def _action_(self,
                 action: UsageAction,
                 input_opts: dict,
                 output_opts: dict):
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
