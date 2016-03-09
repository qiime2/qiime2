# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import pkg_resources

from qiime.sdk.workflow import WorkflowTemplate


class Plugin(object):

    def __init__(self, name, version, website, package):
        """

        """
        self.package = package
        self.name = name
        self.version = version
        self.website = website
        self.workflows = {}

    def register_workflow(self, workflow):
        fn = pkg_resources.resource_filename(self.package, workflow)
        w = WorkflowTemplate.from_markdown(fn)
        self.workflows[w.name] = w

    def register_function(self, name, function, inputs, outputs, doc=""):
        # TODO where is the best place to convert outputs as a list of tuples
        # into an OrderedDict?
        outputs = collections.OrderedDict(outputs)
        w = WorkflowTemplate.from_function(function, inputs, outputs, name, doc)
        self.workflows[w.name] = w
