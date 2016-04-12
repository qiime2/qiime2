# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import pkg_resources

import qiime.sdk
from qiime.core.type.type import BaseType
from qiime.core.type.primitive import Int, Str, Float

__all__ = ['Plugin', 'Type', 'Int', 'Str', 'Float']


class Plugin:
    def __init__(self, name, version, website, package):
        self.package = package
        self.name = name
        self.version = version
        self.website = website
        self.workflows = {}

    def register_workflow(self, workflow):
        fn = pkg_resources.resource_filename(self.package, workflow)
        w = qiime.sdk.Workflow.from_markdown(fn)
        self.workflows[w.id] = w

    def register_function(self, name, function, inputs, outputs, doc=""):
        # TODO where is the best place to convert outputs as a list of tuples
        # into an OrderedDict?
        outputs = collections.OrderedDict(outputs)
        w = qiime.sdk.Workflow.from_function(function, inputs, outputs, name,
                                             doc)
        self.workflows[w.id] = w

    def __eq__(self, other):
        return (
            self.package == other.package and
            self.name == other.name and
            self.version == other.version and
            self.website == other.website and
            self.workflows == other.workflows
        )


class Type(BaseType, fields=('Artifact', 'Metadata')):
    class Artifact:
        def save(self, data, data_writer):
            pass

        def load(self, data_reader):
            pass

    class Metadata:
        def get_columns(self, data):
            pass

        def get_series(self, data, column):
            pass
