# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
import collections
import os

from qiime.sdk import PluginManager
from qiime.sdk import Artifact
from qiime.interface.q2d3.context import Q2D3Context

pm = PluginManager()


def get_import_types():
    return [('FeatureTable[Frequency]', 'Feature Table of Frequencies'),
            ('Phylogeny', 'Phylogenetic Tree')]


def get_artifacts(named=False):
    context_manager = Q2D3Context(os.getcwd())
    return sorted([(context_manager.names[Artifact(fp).uuid], Artifact(fp))
                   if named else Artifact(fp)
                   for fp in context_manager.data.values()],
                  key=lambda art: repr(art[1].type) if named else art.type)


def get_workflows():
    artifacts = get_artifacts()
    available = collections.defaultdict(dict)
    unavailable = collections.defaultdict(dict)
    global pm
    for pname, plugin in pm.plugins.items():
        for wname, workflow in plugin.workflows.items():
            completed = True
            for input_name in workflow.signature.input_artifacts:
                if not get_input_artifacts(workflow.signature, input_name,
                                           artifacts):
                    completed = False
            if completed:
                available[pname][wname] = workflow
            else:
                unavailable[pname][wname] = workflow
    return dict(available), dict(unavailable)


def get_input_artifacts(signature, name, artifacts):
    inputs = []
    target_type = signature.input_artifacts[name]
    for artifact in artifacts:
        if artifact.type < target_type:
            inputs.append(artifact)
    return inputs


def get_plugin_workflow(plugin, workflow):
    plugin = pm.plugins[plugin]
    workflow = plugin.workflows[workflow]
    return plugin, workflow


def get_artifact_name(uuid):
    context_manager = Q2D3Context(os.getcwd())
    return context_manager.names[uuid]
