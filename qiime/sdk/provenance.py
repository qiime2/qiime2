# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


class Provenance:
    def __init__(self, execution_uuid, method_reference, artifact_uuids,
                 parameters):
        self._execution_uuid = execution_uuid
        self._method_reference = method_reference
        self._artifact_uuids = artifact_uuids
        self._parameters = parameters

    @property
    def execution_uuid(self):
        return self._execution_uuid

    @property
    def method_reference(self):
        return self._method_reference

    @property
    def artifact_uuids(self):
        return self._artifact_uuids

    @property
    def parameters(self):
        return self._parameters

    def __eq__(self, other):
        return (
            isinstance(other, Provenance) and
            (self.execution_uuid == other.execution_uuid) and
            (self.method_reference == other.method_reference) and
            (self.artifact_uuids == other.artifact_uuids) and
            (self.parameters == other.parameters)
        )
