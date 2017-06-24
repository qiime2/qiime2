# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import time
import collections
import pkg_resources
import uuid
import copy
import importlib
import shutil
import sys
from datetime import datetime

import distutils
import yaml
import tzlocal
import dateutil.relativedelta as relativedelta

import qiime2
import qiime2.core.util as util


# Used to give PyYAML something to recognize for custom tags
ForwardRef = collections.namedtuple('ForwardRef', ['reference'])
NoProvenance = collections.namedtuple('NoProvenance', ['uuid'])
MetadataPath = collections.namedtuple('MetadataPath', ['path'])
ColorPrimitive = collections.namedtuple('ColorPrimitive', ['hex'])
LiteralString = collections.namedtuple('LiteralString', ['string'])


class OrderedKeyValue(collections.OrderedDict):
    pass


# Used for yaml that looks like:
#   - key1: value1
#   - key2: value2
yaml.add_representer(OrderedKeyValue, lambda dumper, data:
                     dumper.represent_list([
                        {k: v} for k, v in data.items()]))


# Controlling the order of dictionaries (even if semantically irrelevant) is
# important to making it look nice.
yaml.add_representer(collections.OrderedDict, lambda dumper, data:
                     dumper.represent_dict(data.items()))


# LiteralString uses the | character and has literal newlines
yaml.add_representer(LiteralString, lambda dumper, data:
                     dumper.represent_scalar('tag:yaml.org,2002:str',
                                             data.string, style='|'))


# Make our timestamps pretty (unquoted).
yaml.add_representer(datetime, lambda dumper, data:
                     dumper.represent_scalar('tag:yaml.org,2002:timestamp',
                                             data.isoformat()))


# Forward reference to something else in the document, namespaces are
# delimited by colons (:).
yaml.add_representer(ForwardRef, lambda dumper, data:
                     dumper.represent_scalar('!ref', data.reference))


# This tag represents an artifact without provenance, this is to support
# archive format v0. Ideally this won't be seen in the wild in practice.
yaml.add_representer(NoProvenance, lambda dumper, data:
                     dumper.represent_scalar('!no-provenance', str(data.uuid)))


# A reference to Metadata and MetadataCategory who's data can be found at the
# relative path indicated as its value
yaml.add_representer(MetadataPath, lambda dumper, data:
                     dumper.represent_scalar('!metadata', data.path))

# A color primitive.
yaml.add_representer(ColorPrimitive, lambda dumper, data:
                     dumper.represent_scalar('!color', data.hex))


class ProvenanceCapture:
    ANCESTOR_DIR = 'artifacts'
    ACTION_DIR = 'action'
    ACTION_FILE = 'action.yaml'

    def __init__(self):
        self.start = time.time()
        self.uuid = uuid.uuid4()
        self.end = None
        self.plugins = collections.OrderedDict()

        # For the purposes of this dict, `return` is a special case for output
        # we expect to transform this later when serializing, but this lets
        # us treat all transformations uniformly.
        self.transformers = collections.OrderedDict()

        self._build_paths()

    def _destructor(self):
        self.path._destructor()

    def _build_paths(self):
        self.path = qiime2.core.path.ProvenancePath()

        self.ancestor_dir = self.path / self.ANCESTOR_DIR
        self.ancestor_dir.mkdir()

        self.action_dir = self.path / self.ACTION_DIR
        self.action_dir.mkdir()

    def add_ancestor(self, artifact):
        other_path = artifact._archiver.provenance_dir
        if other_path is None:
            # The artifact doesn't have provenance (e.g. version 0)
            # it would be possible to invent a metadata.yaml, but we won't know
            # the framework version for the VERSION file. Even if we did
            # it won't accomplish a lot and there shouldn't be enough
            # version 0 artifacts in the wild to be important in practice.
            # NOTE: this implies that it is possible for an action.yaml file to
            # contain an artifact UUID that is not in the artifacts/ directory.
            return NotImplemented

        destination = self.ancestor_dir / str(artifact.uuid)
        if destination.exists():
            # This artifact is already in the provenance (and so are its
            # ancestors)
            return

        # Handle root node of ancestor
        shutil.copytree(
            str(other_path), str(destination),
            ignore=shutil.ignore_patterns(self.ANCESTOR_DIR + '*'))

        # Handle ancestral nodes of ancestor
        grandcestor_path = other_path / self.ANCESTOR_DIR
        if grandcestor_path.exists():
            for grandcestor in grandcestor_path.iterdir():
                destination = self.ancestor_dir / grandcestor.name
                if not destination.exists():
                    shutil.copytree(str(grandcestor), str(destination))

    def reference_plugin(self, plugin):
        self.plugins[plugin.name] = plugin
        return ForwardRef('environment:plugins:' + plugin.name)

    def capture_env(self):
        return collections.OrderedDict(
            (d.project_name, d.version) for d in pkg_resources.working_set)

    def transformation_recorder(self, name):
        # TODO: this is currently stubbed, but not used.
        record = self.transformers[name] = []
        return record.append

    def _ts_to_date(self, ts):
        return datetime.fromtimestamp(ts, tzlocal.get_localzone())

    def make_execution_section(self):
        execution = collections.OrderedDict()
        execution['uuid'] = str(self.uuid)
        execution['runtime'] = runtime = collections.OrderedDict()
        runtime['start'] = start = self._ts_to_date(self.start)
        runtime['end'] = end = self._ts_to_date(self.end)
        runtime['duration'] = \
            util.duration_time(relativedelta.relativedelta(end, start))

        return execution

    def make_env_section(self):
        env = collections.OrderedDict()
        env['platform'] = pkg_resources.get_build_platform()
        # There is a trailing whitespace in sys.version, strip so that YAML can
        # use literal formatting.
        env['python'] = LiteralString('\n'.join(line.strip() for line in
                                      sys.version.split('\n')))
        env['framework'] = qiime2.__version__
        env['plugins'] = self.plugins
        env['python-packages'] = self.capture_env()

        return env

    def write_action_yaml(self):
        settings = dict(default_flow_style=False, indent=4)
        with (self.action_dir / self.ACTION_FILE).open(mode='w') as fh:
            fh.write(yaml.dump({'execution': self.make_execution_section()},
                               **settings))
            fh.write('\n')
            fh.write(yaml.dump({'action': self.make_action_section()},
                               **settings))
            fh.write('\n')
            fh.write(yaml.dump({'environment': self.make_env_section()},
                               **settings))

    def finalize(self, final_path, node_members):
        self.end = time.time()

        for member in node_members:
            shutil.copy(str(member), str(self.path))

        self.write_action_yaml()

        self.path.rename(final_path)

    def fork(self):
        forked = copy.copy(self)
        # Unique `result` key for each output of an action
        forked.transformers = forked.transformers.copy()
        # create a copy of the backing dir so factory (the hard stuff is
        # mostly done by this point)
        forked._build_paths()
        distutils.dir_util.copy_tree(str(self.path), str(forked.path))

        return forked


class ImportProvenanceCapture(ProvenanceCapture):
    def __init__(self, format=None, checksums=None):
        super().__init__()
        self.format_name = format.__name__ if format is not None else None
        self.checksums = checksums

    def make_action_section(self):
        action = collections.OrderedDict()
        action['type'] = 'import'
        if self.format_name is not None:
            action['format'] = self.format_name
        if self.checksums is not None:
            action['manifest'] = [
                collections.OrderedDict([('name', name), ('md5sum', md5sum)])
                for name, md5sum in self.checksums.items()]

        return action


class ActionProvenanceCapture(ProvenanceCapture):
    def __init__(self, action_type, import_path, action_id):
        super().__init__()
        self._plugin = importlib.import_module(import_path).__plugin__
        self.action = self._plugin.actions[action_id]
        self.action_type = action_type
        self.inputs = OrderedKeyValue()
        self.parameters = OrderedKeyValue()

    def handle_metadata(self, name, value):
        if value is None:
            return None

        if isinstance(value, qiime2.MetadataCategory):
            pandas_obj = value.to_series()
        elif isinstance(value, qiime2.Metadata):
            pandas_obj = value.to_dataframe()
        else:
            raise NotImplementedError

        uuid_ref = ""
        if value.artifacts:
            uuids = []
            for artifact in value.artifacts:
                uuids.append(str(artifact.uuid))
                self.add_ancestor(artifact)
            uuid_ref = ",".join(uuids) + ":"

        relpath = name + '.tsv'
        pandas_obj.to_csv(str(self.action_dir / relpath), sep='\t')

        return MetadataPath(uuid_ref + relpath)

    def add_parameter(self, name, type_expr, parameter):
        type_map = {
            'Color': ColorPrimitive,
            'Metadata': lambda x: self.handle_metadata(name, x),
            'MetadataCategory': lambda x: self.handle_metadata(name, x)
            # TODO: handle collection primitives (not currently used)
        }

        handler = type_map.get(type_expr.to_ast()['name'], lambda x: x)
        self.parameters[name] = handler(parameter)

    def add_input(self, name, artifact):
        if artifact is None:
            self.inputs[name] = artifact
        else:
            ancestral_provenance = self.add_ancestor(artifact)
            if ancestral_provenance is NotImplemented:
                self.inputs[name] = NoProvenance(artifact.uuid)
            else:
                self.inputs[name] = str(artifact.uuid)

    def make_action_section(self):
        action = collections.OrderedDict()
        action['type'] = self.action_type
        action['plugin'] = self.reference_plugin(self._plugin)
        action['action'] = self.action.id
        action['inputs'] = self.inputs
        action['parameters'] = self.parameters

        return action
