# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources

from qiime.sdk.workflow import Workflow


class Plugin(object):

    def __init__(self, name, version, website, package):
        """

        """
        self.package = package
        self.name = name
        self.version = version
        self.website = website
        self.workflows = {}
        self.artifact_types = {}
        self.traits = {}

    def register_workflow(self, workflow):
        fn = pkg_resources.resource_filename(self.package, workflow)
        w = Workflow.from_markdown(fn)
        self.workflows[w.name] = w

    def register_function(self, name, function, inputs, outputs, doc=""):
        w = Workflow.from_function(function, inputs, outputs, name, doc)
        self.workflows[w.name] = w

    def register_artifact_type(self, artifact_type):
        self.artifact_types[artifact_type.__name__] = artifact_type

    def register_trait(self, trait):
        self.traits[trait.__name__] = trait
