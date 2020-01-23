# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import uuid as _uuid

import yaml

import qiime2.sdk as sdk

# Allow OrderedDict to be serialized for YAML representation
yaml.add_representer(collections.OrderedDict, lambda dumper, data:
                     dumper.represent_dict(data.items()))


class ArchiveFormat:
    DATA_DIR = 'data'
    METADATA_FILE = 'metadata.yaml'

    @classmethod
    def _parse_metadata(self, fh, expected_uuid):
        metadata = yaml.safe_load(fh)
        if metadata['uuid'] != str(expected_uuid):
            raise ValueError(
                "Archive root directory must match UUID present in archive's"
                " metadata: %s != %s" % (expected_uuid,  metadata['uuid']))

        return metadata['uuid'], metadata['type'], metadata['format']

    @classmethod
    def _format_metadata(self, fh, uuid, type, format):
        metadata = collections.OrderedDict()
        metadata['uuid'] = str(uuid)
        metadata['type'] = repr(type)
        metadata['format'] = None
        if format is not None:
            metadata['format'] = format.__name__

        fh.write(yaml.dump(metadata, default_flow_style=False))

    @classmethod
    def load_metadata(self, archive):
        with archive.open(self.METADATA_FILE) as fh:
            return self._parse_metadata(fh, expected_uuid=archive.uuid)

    @classmethod
    def write(cls, archive_record, type, format, data_initializer, _):
        root = archive_record.root
        metadata_fp = root / cls.METADATA_FILE

        with metadata_fp.open(mode='w') as fh:
            cls._format_metadata(fh, archive_record.uuid, type, format)

        data_dir = root / cls.DATA_DIR
        data_dir.mkdir()

        data_initializer(data_dir)

    def __init__(self, archive_record):
        path = archive_record.root

        with (path / self.METADATA_FILE).open() as fh:
            uuid, type, format = \
                self._parse_metadata(fh, expected_uuid=archive_record.uuid)

        self.uuid = _uuid.UUID(uuid)
        self.type = sdk.parse_type(type)
        self.format = sdk.parse_format(format)

        self.path = path
        self.data_dir = path / self.DATA_DIR
