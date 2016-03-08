# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from collections import OrderedDict
import abc

from .artifact import Artifact

class SystemContext(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def _output_to_filepath(self, name, output_type, workflow_template,
                            artifacts, parameters):
        pass

    @abc.abstractmethod
    def _uuid_to_filepath(self, uuid):
        pass

    def __call__(self, workflow_template, artifact_uuids, parameters):
        artifacts = self._load_artifacts(artifact_uuids)
        outputs = worflow_template.signature(artifacts, parameters)
        # create setup_lines
        setup_lines = self._system_context_setup_lines()
        for name, artifact_uuid in artifact_uuids.items():
            setup_lines.append(
                self._create_artifact_setup_lines(name, artifact_uuid))
        for name, value in parameters.items():
            setup_lines.append(self._create_parameter_setup_lines(name, value))

        # create teardown lines
        teardown_lines = []
        output_filepaths = OrderedDict()
        for name, output_type in outputs:
            provenance = None
            output_filepath = self._output_to_filepath(
                name, output_type, workflow_template, artifacts, parameters)
            output_filepaths[name] = output_filepath
            teardown_lines.append(self._create_output_teardown_lines(
                name, output_type, provenance, output_filepath))
        teardown_lines.extend(self._system_context_teardown_lines())

        # create and return job
        return workflow_template.create_job(setup_lines, teardown_lines)

    def _load_artifacts(self, artifact_uuids):
        result = {}
        for name, uuid in artifact_uuids.items():
            result[name] = Artifact(self._uuid_to_filepath(uuid))
        return result

    def _system_context_setup_lines(self):
        return ['from qiime.sdk.artifact import Artifact']

    def _create_artifact_setup_lines(self, name, artifact_uuid):
        artifact_filepath = self._uuid_to_filepath(artifact_uuid)
        return ['%s = Artifact(%r)' % (name, artifact_filepath)]

    def _create_parameter_setup_lines(self, name, value):
        result = []
        for name, value:
            result.append('%s = %r' % (name, value))
        return result

    def _create_output_teardown_lines(self, name, output_type, provenance, output_reference):
        return ['Artifact.save(%r, %r, %r, %r)' %
                (name, output_type, provenance, output_reference)]

    def _system_context_teardown_lines(self):
        return []
