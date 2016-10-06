import time
import collections
import pkg_resources
import uuid
import copy
import importlib
import pathlib
import shutil
import tempfile
import os
import sys
from datetime import datetime

import distutils
import yaml
import tzlocal
import dateutil.relativedelta as relativedelta

import qiime
import qiime.core.util as util


ForwardRef = collections.namedtuple('ForwardRef', ['reference'])
MetadataPath = collections.namedtuple('MetadataPath', ['path'])
ColorPrimitive = collections.namedtuple('ColorPrimitive', ['hex'])


yaml.add_representer(collections.OrderedDict, lambda dumper, data:
                     dumper.represent_dict(data.items()))

yaml.add_representer(datetime, lambda dumper, data:
                     dumper.represent_scalar('tag:yaml.org,2002:timestamp',
                                             data.isoformat()))

yaml.add_representer(ForwardRef, lambda dumper, data:
                     dumper.represent_scalar('!ref', data.reference))

yaml.add_representer(MetadataPath, lambda dumper, data:
                     dumper.represent_scalar('!metadata', data.path))

yaml.add_representer(ColorPrimitive, lambda dumper, data:
                     dumper.represent_scalar('!color', data.hex))


class ProvenanceCapture:
    ANCESTOR_DIR = 'artifacts'
    ACTION_DIR = 'action'
    ACTION_FILE = 'action.yml'

    def __init__(self):
        self.start = time.time()
        self.uuid = uuid.uuid4()
        self.end = None
        self.plugins = collections.OrderedDict()

        # For the purposes of this dict, `return` is a special case for output
        # we expect to transform this later when serializing, but this lets
        # us treat all transformations uniformly.
        self.transformers = collections.OrderedDict()

        # TODO: normalize `mkdtemp` when we have framework temp locations.
        self.path = pathlib.Path(tempfile.mkdtemp(prefix='qiime2-prov-'))
        self._build_paths()

    def _build_paths(self):
        self.ancestor_dir = self.path / self.ANCESTOR_DIR
        self.ancestor_dir.mkdir()

        self.action_dir = self.path / self.ACTION_DIR
        self.action_dir.mkdir()

    def add_ancestor(self, artifact):
        other_path = artifact._archiver.provenance_dir

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
        relative_delta = relativedelta.relativedelta(end, start)
        runtime['duration'] = util.human_readable(relative_delta)

        return execution

    def make_env_section(self):
        env = collections.OrderedDict()
        env['platform'] = pkg_resources.get_build_platform()
        env['python'] = ' '.join(sys.version.split('\n'))
        env['framework'] = qiime.__version__
        env['plugins'] = self.plugins
        env['pip'] = self.capture_env()

        return env

    def write_action_yaml(self):
        action_data = collections.OrderedDict()
        action_data['execution'] = self.make_execution_section()
        action_data['action'] = self.make_action_section()
        action_data['environment'] = self.make_env_section()

        with (self.action_dir / self.ACTION_FILE).open(mode='w') as fh:
            fh.write(yaml.dump(action_data, default_flow_style=False,
                               width=float("inf")))


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
        forked.path = pathlib.Path(tempfile.mkdtemp(prefix='qiime2-prov-'))
        forked._build_paths()

        distutils.dir_util.copy_tree(str(self.path), str(forked.path))

        return forked

    def __del__(self):
        # Used to delete the original source when "forking"
        if self.path.exists():
            shutil.rmtree(str(self.path))


class ImportProvenanceCapture(ProvenanceCapture):
    def __init__(self, format=None, checksums=None):
        super().__init__()
        self.format_name = format.__name__ if format is not None else None
        self.checksums = checksums

    def make_action_section(self):
        action = collections.OrderedDict()
        action['type'] = 'import'
        # if self.citations:
        #     action['citations'] = self.citations
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
        self.inputs = collections.OrderedDict()
        self.parameters = collections.OrderedDict()

    def add_parameter(self, name, type_expr, parameter):
        pass

    def make_action_section(self):
        action = collections.OrderedDict()
        action['type'] = self.action_type
        action['plugin'] = self.reference_plugin(self._plugin)
        action['action'] = self.action.id
        # TODO: citations
        action['inputs'] = [
            {key: str(uuid)} for key, uuid in self.inputs.items()]

        return action

    def add_input(self, name, artifact):
        self.add_ancestor(artifact)
        self.inputs[name] = artifact.uuid
