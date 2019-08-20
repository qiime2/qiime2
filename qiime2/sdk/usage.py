# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import abc
import collections
import contextlib
import os
import re
import tempfile
import zipfile

import qiime2
from qiime2.core.type.util import is_metadata_column_type


RefAction = collections.namedtuple('RefAction', 'plugin_id action_name')
class RefInputs(dict): pass
class RefOutputs(dict): pass


class ScopeRecord:
    def __init__(self, name, factory=None, is_bound=False,
            assert_has_line_matching=None):
        self.name = name
        self.factory = factory
        self.is_bound = is_bound
        self._assert_has_line_matching_ = assert_has_line_matching

    def __repr__(self):
        return 'ScopeRecord(%s, %s, %s)' % (self.name, self.factory,
                                            self.is_bound)

    def assert_has_line_matching(self ,label, path, expression):
        # TODO: did somebody say "curry"?
        return self._assert_has_line_matching_(label, self.name, path,
                                              expression)


class Scope:
    """
    TODO: docstring
    """

    # TODO: when would i ever seed scope with records?
    def __init__(self, records=None):
        if records is None:
            self.records = dict()
        else:
            self.records = dict(records)

    def keys(self):
        return self.records.keys()

    def values(self):
        return self.records.values()

    def __len__(self):
        return len(self.records)

    def __str__(self):
        return str(self.records)

    def __iter__(self):
        yield from self.records.values()

    # TODO: this might not work quite like we would hope, since staring at it
    # makes something pop into existence
    def __getitem__(self, key):
        if key not in self.records:
            self.records[key] = ScopeRecord(
                name=key,
                assert_has_line_matching=self._assert_has_line_matching_)
        return self.records[key]

    def __setitem__(self, key, value):
        self.records[key] = ScopeRecord(
            name=key, factory=value, is_bound=True,
            assert_has_line_matching=self._assert_has_line_matching_)

    def __contains__(self, key):
        return key in self.records

    def _bind(self, name, obj):
        setattr(self, name, obj)


class Usage(metaclass=abc.ABCMeta):
    """
    TODO: docstring
    """

    RefAction = RefAction
    RefInputs = RefInputs
    RefOutputs = RefOutputs

    def __init__(self, scope=None):
        self.plugin_manager = None

        if scope is None:
            self.scope = Scope()

        self.scope._bind('_assert_has_line_matching_',
                         self._assert_has_line_matching_)

    def get_action(self, action: RefAction):
        if self.plugin_manager is None:
            self.plugin_manager = qiime2.sdk.PluginManager()

        try:
            plugin = self.plugin_manager.get_plugin(action.plugin_id)
        except KeyError:
            raise KeyError('No plugin currently registered with id: "%s".' %
                           (action.plugin_id, ))
        return plugin.actions[action.action_name]

    @abc.abstractmethod
    def comment(self, text):
        raise NotImplementedError

    @abc.abstractmethod
    def action(self, action: RefAction, inputs: RefInputs,
               outputs: RefOutputs=None):
        raise NotImplementedError

    @abc.abstractmethod
    def import_file(self, type, input_name, output_name=None, format=None):
        raise NotImplementedError

    @abc.abstractmethod
    def export_file(self, input_name, output_name, format=None):
        raise NotImplementedError

    @abc.abstractmethod
    def _assert_has_line_matching_(self, label, result, path, expression):
        raise NotImplementedError


class NoOpUsage(Usage):
    """
    TODO: docstring
    """

    def get_action(self, action: RefAction):
        return None

    def comment(self, text):
        return None

    def action(self, action: RefAction, inputs: RefInputs,
               outputs: RefOutputs=None):
        return None

    def import_file(self, type, input_name, output_name=None, format=None):
        return None

    def export_file(self, input_name, output_name, format=None):
        return None

    def _assert_has_line_matching_(self, label, result, path, expression):
        return None


class ExecutionUsage(Usage):
    """
    TODO: docstring
    """

    def comment(self, text):
        return None

    def import_file(self, type, input_name, output_name=None, format=None):
        return None

    def export_file(self, input_name, output_name, format=None):
        return None

    def _assert_has_line_matching_(self, label, result, path, expression):
        data = self.scope[result].factory()
        with tempfile.TemporaryDirectory(prefix='q2-exc-usage-') as temp_dir:
            fp = data.save(os.path.join(temp_dir, result))
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
                        raise ValueError('omz')

        return None

    def action(self, action: RefAction, inputs: RefInputs,
               outputs: RefOutputs=None):
        action_f = self.get_action(action)
        sig = action_f.signature
        opts = {}

        for name, signature in sig.signature_order.items():
            if name in inputs:
                if is_metadata_column_type(signature.qiime_type):
                    ref, col = inputs[name]
                    value = self.scope[ref].factory()
                    value = value.get_column(col)
                elif inputs[name] in self.scope.values() \
                        and inputs[name].is_bound:
                    value = inputs[name].factory()
                else:
                    value = inputs[name]
                opts[name] = value

        results = action_f(**opts)

        for output in sig.outputs.keys():
            scope_record = outputs[output] if output in outputs else output
            result = getattr(results, output)
            self.scope[scope_record.name] = lambda: result
