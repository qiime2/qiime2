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


ScopeRecord = collections.namedtuple('ScopeRecord',
                                     ['name', 'type', 'factory'])


class Scope:
    """
    TODO: docstring
    """

    def __init__(self, records=None):
        if records is None:
            self.records = dict()
        else:
            self.records = dict(records)

    def _add_record(self, name, factory, type):
        if name in self.records:
            record = self.records[name]
            if type != record.type or factory != record.factory:
                raise KeyError('A(n) %s record with the name "%s" is already '
                               'registered.' % (type, name))

        self.records[name] = ScopeRecord(name=name, type=type, factory=factory)

    def add_artifact(self, name, factory):
        self._add_record(name, factory, 'artifact')

    def add_visualization(self, name, factory):
        self._add_record(name, factory, 'visualization')

    def add_file(self, name, factory):
        self._add_record(name, factory, 'file')

    def add_metadata(self, name, factory):
        self._add_record(name, factory, 'metadata')

    def keys(self):
        return self.records.keys()

    def __len__(self):
        return len(self.records)

    def __str__(self):
        return str(self.records)

    def __iter__(self):
        yield from self.records.values()

    def __getitem__(self, key):
        # TODO: Better error here?
        return self.records[key]

    def __contains__(self, key):
        return key in self.records

    def assert_has_line_matching(self, label, result, path, expression):
        pass


class Usage(metaclass=abc.ABCMeta):
    """
    TODO: docstring
    """

    def __init__(self, scope=None):
        self.plugin_manager = None

        if scope is None:
            self.scope = Scope()

        self.scope.assert_has_line_matching = self._assert_has_line_matching_

    def get_action(self, plugin_id, action):
        if self.plugin_manager is None:
            self.plugin_manager = qiime2.sdk.PluginManager()

        try:
            plugin = self.plugin_manager.get_plugin(plugin_id)
        except KeyError:
            raise KeyError('No plugin currently registered with id: "%s".' %
                           (plugin_id, ))
        return plugin.actions[action]

    @abc.abstractmethod
    def comment(self, text):
        raise NotImplementedError

    @abc.abstractmethod
    def action(self, action, inputs, outputs=None):
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

    def get_action(self, plugin_id, action):
        return None

    def comment(self, text):
        return None

    def action(self, action, inputs, outputs=None):
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

    def action(self, action, inputs, outputs=None):
        sig = action.signature
        opts = {}

        for name, signature in sig.signature_order.items():
            if name in inputs:
                if is_metadata_column_type(signature.qiime_type):
                    ref, col = inputs[name]
                    value = self.scope[ref].factory()
                    value = value.get_column(col)
                elif inputs[name] in self.scope:
                    value = self.scope[inputs[name]].factory()
                else:
                    value = inputs[name]
                opts[name] = value

        results = action(**opts)

        for output in sig.outputs.keys():
            scope_name = outputs[output] if output in outputs else output
            result = getattr(results, output)
            if isinstance(result, qiime2.Artifact):
                self.scope.add_artifact(scope_name, lambda: result)
            else:
                self.scope.add_visualization(scope_name, lambda: result)
