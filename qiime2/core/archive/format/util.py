# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import uuid
from contextlib import contextmanager
from collections import OrderedDict

import yaml
import networkx as nx

from qiime2.core.archive import Archiver
from qiime2.core.archive.provenance import NoProvenance


@contextmanager
def artifact_version(version):
    version = str(version)
    if version not in Archiver._FORMAT_REGISTRY:
        raise ValueError("Version %s not supported" % version)
    original_version = Archiver.CURRENT_FORMAT_VERSION
    try:
        Archiver.CURRENT_FORMAT_VERSION = version
        yield
    finally:
        Archiver.CURRENT_FORMAT_VERSION = original_version


# parse_actions is a private utility helper that accepts v0- and v1- format
# Artifacts, and returns a NetworkX DiGraph where the nodes represent artifacts
# and the edges represent actions. The nodes and edges each carry the necessary
# attributes (as extracted from the provenance) to recreate the execution
# pipeline necessary to derive the  artifact being processed.
# This utility is intended to be a stop-gap measure for now, and we anticipate
# replacing it with a more full-featured provenance object.
def parse_actions(artifact):
    prov_dir = artifact._archiver.provenance_dir
    prov = nx.DiGraph()

    # v0-format artifact, nothing we can do here.
    if prov_dir is None or not prov_dir.exists():
        prov.add_node(artifact.uuid)
        return prov

    def dict_list_to_ordereddict(l):
        return OrderedDict([tuple(x.items())[0] for x in l])

    def load_yaml(prov_dir):
        action_yaml = prov_dir / 'action' / 'action.yaml'
        metadata_yaml = prov_dir / 'metadata.yaml'
        with action_yaml.open() as fh_a, metadata_yaml.open() as fh_b:
            return yaml.load(fh_a), yaml.load(fh_b)

    def add_edges(action, metadata):
        id = uuid.UUID(metadata['uuid'])
        action = action['action']  # pardon the stutter
        if action.get('inputs') is not None:
            # Non-import action
            for i in action['inputs']:
                for k, v in i.items():
                    if type(v) is not NoProvenance:
                        # We want to keep the `NoProvenance` type, it seems
                        # useful for downstream processing. Otherwise,
                        # nodes are a `uuid.UUID`.
                        v = uuid.UUID(v)
                    params = dict_list_to_ordereddict(action['parameters'])
                    prov.add_edge(v, id,
                                  type=action['type'],
                                  plugin=action['plugin'],
                                  action=action['action'],
                                  parameters=params)
            for attr in ['type', 'format']:
                nx.set_node_attributes(prov, attr, {id: metadata[attr]})
        else:
            # Import
            prov.add_node(id, type=metadata['type'],
                          format=metadata['format'],
                          manifest=action.get('manifest'))  # not always there
    base_action, base_metadata = load_yaml(prov_dir)
    add_edges(base_action, base_metadata)
    artifacts_dir = prov_dir / 'artifacts'
    # It doesn't matter what order we walk over things, NetworkX will handle
    # wiring everything up.
    for artifact_dir in [a for a in artifacts_dir.iterdir() if a.is_dir()]:
        action, metadata = load_yaml(artifact_dir)
        add_edges(action, metadata)
    return prov
