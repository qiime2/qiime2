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

import qiime2


ScopeRecord = collections.namedtuple('ScopeRecord',
                                     ['name', 'type', 'factory'])


class Scope:
    def __init__(self, records=None):
        if records is None:
            self.records = dict()
        else:
            self.records = dict(records)

    def _add_record(self, name, factory, type):
        # TODO: do i need this?
        # if name in self.records:
        #     raise ValueError

        self.records[name] = ScopeRecord(name=name, type=type, factory=factory)

    def add_artifact(self, name, factory):
        self._add_record(name, factory, 'artifact')

    def add_visualization(self, name, factory):
        self._add_record(name, factory, 'visualization')

    def add_file(self, name, factory):
        self._add_record(name, factory, 'file')

    def add_metadata(self, name, factory):
        self._add_record(name, factory, 'metadata')

    def __iter__(self):
        yield from self.records.values()

    def __getitem__(self, key):
        return self.records[key]


class Usage(metaclass=abc.ABCMeta):
    def __init__(self):
        self.scope = None
        self.plugin_manager = None

    @contextlib.contextmanager
    def bind(self, scope):
        self.scope = scope
        yield self
        self.scope = None

    def get_action(self, plugin, action):
        if self.plugin_manager is None:
            self.plugin_manager = qiime2.sdk.PluginManager()

        plugin = self.plugin_manager.plugins[plugin]
        return plugin.actions[action]

    @abc.abstractmethod
    def comment(self, comment):
        raise NotImplementedError

    @abc.abstractmethod
    def action(self, action, inputs, outputs):
        raise NotImplementedError


class NoOpUsage(Usage):
    def get_action(self, plugin, action):
        return None

    def comment(self, comment):
        return None

    def action(self, action, inputs, outputs):
        return None


# TODO: think about this one
# class TypeCheckUsage(Usage):
#     pass


# TODO: this doesn't really work yet, needs to backtrace input refs, and
# also update the scope as it goes with the stuff that it creates
class ExecutionUsage(Usage):
    def __init__(self):
        super().__init__()
        # TODO: not sure about this
        self.inputs = dict()
        self.outputs = dict()

    def comment(self, comment):
        return None

    def action(self, action, inputs, outputs):
        sig_inputs, sig_params = dict(), dict()

        # TODO: deal with MD
        for param, _ in action.signature.parameters.items():
            p = inputs.pop(param, None)
            if p is not None:
                sig_params[param] = p

        for input, _ in action.signature.inputs.items():
            try:
                record = self.scope[input]
            except KeyError:
                raise ValueError('Missing data registration for input %s' %
                                 (input, ))

            if record.type == 'artifact':
                data = record.factory()
            else:
                # TODO: add in the other branches here
                data = None

            sig_inputs[input] = data
            self.inputs[data] = input

        outputs = action(**sig_inputs, **sig_params)
        for output, _ in action.signature.outputs.items():
            self.outputs[getattr(outputs, output)] = output
