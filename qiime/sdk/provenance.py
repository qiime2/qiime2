# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import uuid


# TODO Provenance is stubbed for now. Type-checking here is pretty restrictive
# to ensure these objects are easily serializable and can be round-tripped.
class Provenance:
    def __init__(self, execution_uuid, executor_reference, artifact_uuids,
                 parameter_references):
        if not isinstance(execution_uuid, uuid.UUID):
            raise TypeError(
                "`execution_uuid` must be type %r, not %r"
                % (uuid.UUID.__name__, type(execution_uuid).__name__))
        self._execution_uuid = execution_uuid

        if not isinstance(executor_reference, str):
            raise TypeError(
                "`executor_reference` must be type %r, not %r"
                % (str.__name__, type(executor_reference).__name__))
        self._executor_reference = executor_reference

        for uuid_ in artifact_uuids.values():
            if not isinstance(uuid_, uuid.UUID):
                raise TypeError(
                    "Values of `artifact_uuids` must be type %r, not %r"
                    % (uuid.UUID.__name__, type(uuid_).__name__))
        self._artifact_uuids = artifact_uuids

        for param_ref in parameter_references.values():
            if not isinstance(param_ref, str):
                raise TypeError(
                    "Values of `parameter_references` must be type %r, not %r"
                    % (str.__name__, type(param_ref).__name__))
        self._parameter_references = parameter_references

    def __repr__(self):
        return repr(dict(execution_uuid=self.execution_uuid,
                         executor_reference=self.executor_reference,
                         artifact_uuids=self.artifact_uuids,
                         parameter_references=self.parameter_references))

    @property
    def execution_uuid(self):
        return self._execution_uuid

    @property
    def executor_reference(self):
        return self._executor_reference

    @property
    def artifact_uuids(self):
        return self._artifact_uuids

    @property
    def parameter_references(self):
        return self._parameter_references

    def __eq__(self, other):
        return (
            isinstance(other, Provenance) and
            (self.execution_uuid == other.execution_uuid) and
            (self.executor_reference == other.executor_reference) and
            (self.artifact_uuids == other.artifact_uuids) and
            (self.parameter_references == other.parameter_references)
        )

    def __ne__(self, other):
        return not (self == other)
