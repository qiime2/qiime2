# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
import os
import uuid
import tarfile
import tempfile

from qiime.sdk.artifact_data import ArtifactDataReader, ArtifactDataWriter


# TODO use real type when it is ready
class DummyType:
    # model is simply some `bytes`

    @classmethod
    def load(cls, data_reader):
        with data_reader.get_file('model.dat') as fh:
            return fh.read()

    @classmethod
    def save(cls, model, data_writer):
        with data_writer.create_file('model.dat') as fh:
            fh.write(model)


class ArtifactError(Exception):
    pass


class ArtifactLoadError(ArtifactError):
    pass


# TODO make sure files are closed appropriately
# TODO write README
# TODO check file extension on load/save and warn if not .qtf
# TODO make sure path separator code is correct for Unix and Windows systems
class Artifact:
    _metadata_path = 'metadata.txt'

    @classmethod
    def load(cls, tarfilepath):
        if not tarfile.is_tarfile(tarfilepath):
            raise ArtifactLoadError(
                "%r is not a readable tar archive file" % tarfilepath)

        with tarfile.open(tarfilepath, mode='r') as tar:
            root_dir = cls._get_root_dir(tarfilepath)
            type_, uuid_, provenance = cls._load_metadata(tar, root_dir)

            data_reader = ArtifactDataReader(tar, root_dir)
            model = type_.load(data_reader)
            return cls(model, type_, provenance, uuid_=uuid_)

    @classmethod
    def _get_root_dir(cls, tarfilepath):
        return os.path.splitext(os.path.basename(tarfilepath))[0]

    @classmethod
    def _load_metadata(cls, tar, root_dir):
        metadata_path = os.path.join(root_dir, cls._metadata_path)

        try:
            filehandle = tar.extractfile(metadata_path)
        except KeyError:
            raise ArtifactLoadError(
                "Artifact metadata (%r) does not exist" % metadata_path)

        if filehandle is None:
            raise ArtifactLoadError(
                "Artifact metadata (%r) is not a file" % metadata_path)

        with io.TextIOWrapper(filehandle) as fh:
            for line in fh:
                if not line.strip().startswith('#'):
                    break
            type_ = cls._parse_type(line)
            uuid_ = cls._parse_uuid(next(fh))
            provenance = cls._parse_provenance(next(fh))

        filehandle.close()
        return type_, uuid_, provenance

    @classmethod
    def _parse_type(cls, line):
        type_name = line.strip()

        if type_name == 'DummyType':
            return DummyType
        else:
            raise ArtifactLoadError(
                "Unrecognized artifact type %r" % type_name)

    @classmethod
    def _parse_uuid(cls, line):
        return uuid.UUID(hex=line.strip())

    @classmethod
    def _parse_provenance(cls, line):
        return line.strip()

    def __init__(self, model, type_, provenance, uuid_=None):
        self._model = model
        self._type = type_
        self._provenance = provenance

        if uuid_ is None:
            uuid_ = uuid.uuid4()
        elif not isinstance(uuid_, uuid.UUID):
            raise TypeError("`uuid_` must be of type %r, not %r" %
                    (uuid.UUID.__name__, type(uuid_).__name__))
        self._uuid = uuid_

    def __eq__(self, other):
        # TODO this method is mainly for sanity checking right now; revisit
        # semantics if we decide to keep the method
        if not isinstance(other, Artifact):
            return False

        if self._model != other._model:
            return False

        if self._type is not other._type:
            return False

        if self._uuid != other._uuid:
            return False

        if self._provenance != other._provenance:
            return False

        return True

    def save(self, tarfilepath):
        # TODO support compression?
        with tarfile.open(tarfilepath, mode='w') as tar:
            root_dir = self._get_root_dir(tarfilepath)
            self._save_metadata(tar, root_dir)

            data_writer = ArtifactDataWriter(root_dir)
            self._type.save(self._model, data_writer)
            data_writer.save(tar)
            # TODO fix cleanup
            del data_writer

    def _save_metadata(self, tar, root_dir):
        # TODO pass temp dir specified by user in config
        with tempfile.TemporaryDirectory(
                prefix='qiime2-temp-artifact-metadata-') as tempdir:
            fp = os.path.join(tempdir, self._metadata_path)
            with open(fp, mode='w') as fh:
                fh.write('# QIIME 2 artifact metadata\n')
                fh.write(
                    '# WARNING: This is a temporary file format used for '
                    'prototyping QIIME 2. Do not rely on this format for any '
                    'reason as it will not be supported past the prototyping '
                    'stage.\n')

                fh.write(self._formatted_type())
                fh.write(self._formatted_uuid())
                fh.write(self._formatted_provenance())

            # TODO set file metadata appropriately (e.g., owner, permissions)
            tar.add(fp, arcname=os.path.join(root_dir, self._metadata_path), recursive=False)

    def _formatted_type(self):
        return '%s\n' % self._type.__name__

    def _formatted_uuid(self):
        return '%s\n' % self._uuid

    def _formatted_provenance(self):
        return '%s\n' % self._provenance
