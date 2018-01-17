# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import csv
import itertools
import os.path
import re

import numpy as np
import pandas as pd

from qiime2.core.util import find_duplicates


class MetadataFileError(Exception):
    _suffix = (
        "There may be more errors present in the metadata file. To get a full "
        "report, sample/feature metadata files can be validated with Keemei: "
        "http://keemei.qiime.org\n\nFind details on QIIME 2 metadata "
        "requirements here: https://docs.qiime2.org/%s/tutorials/metadata/")

    def __init__(self, message):
        # Lazy import because `qiime2.__release__` is available at runtime but
        # not at import time (otherwise the release value could be interpolated
        # into `_suffix` in the class definition above).
        import qiime2

        super().__init__(message + '\n\n' + self._suffix % qiime2.__release__)


class MetadataReader:
    _supported_column_types = {'categorical', 'numeric'}

    def __init__(self, filepath):
        if not os.path.isfile(filepath):
            raise MetadataFileError(
                "Metadata file path doesn't exist, or the path points to "
                "something other than a file. Please check that the path "
                "exists, has read permissions, and points to a regular file "
                "(not a directory): %s" % filepath)

        self._filepath = filepath

        # Used by `read()` to store an iterator yielding records with
        # leading/trailing whitespace stripped from their cells (this is a
        # preprocessing step that should happen with *every* record). The
        # iterator protocol is the only guaranteed API on this object.
        self._reader = None

    def read(self, into, column_types=None):
        if column_types is None:
            column_types = {}

        try:
            # Newline settings based on recommendation from csv docs:
            #     https://docs.python.org/3/library/csv.html#id3
            with open(self._filepath, 'r', newline='', encoding='utf-8') as fh:
                tsv_reader = csv.reader(fh, dialect='excel-tab', strict=True)
                self._reader = (self._strip_cell_whitespace(record)
                                for record in tsv_reader)
                header = self._read_header()
                directives = self._read_directives(header)
                ids, data = self._read_data(header)
        finally:
            self._reader = None

        index = pd.Index(ids, name=header[0], dtype=object)
        df = pd.DataFrame(data, columns=header[1:], index=index, dtype=object)

        df.replace('', np.nan, inplace=True)

        for name, type in column_types.items():
            if name not in df.columns:
                raise MetadataFileError(
                    "Unrecognized column name %r specified in `column_types`."
                    % name)
            if type not in self._supported_column_types:
                raise MetadataFileError(
                    "Unrecognized column type %r specified in `column_types`. "
                    "Supported column types: %s" %
                    (type, ', '.join(sorted(self._supported_column_types))))

        resolved_column_types = directives.get('types', {})
        resolved_column_types.update(column_types)

        # Cast each column to the appropriate dtype based on column type.
        df = df.apply(self._cast_column, axis='index',
                      column_types=resolved_column_types)

        try:
            return into(df)
        except Exception as e:
            raise MetadataFileError(
                "There was an issue with loading the metadata file:\n\n%s" % e)

    def _read_header(self):
        header = None
        for record in self._reader:
            if self._is_header(record):
                header = record
                break
            elif self._is_comment(record):
                continue
            elif self._is_empty(record):
                continue
            elif self._is_directive(record):
                raise MetadataFileError(
                    "Found directive %r when searching for header. Directives "
                    "may only appear immediately after the header."
                    % record[0])
            else:
                # TODO better error message to hint at what to do
                raise MetadataFileError("Invalid header: %r" % record)

        if header is None:
            raise MetadataFileError(
                "Failed to locate header. The metadata file may be empty, or "
                "consists only of comments or empty records.")

        # Trim trailing empty cells from header.
        data_extent = None
        for idx, cell in enumerate(header):
            if cell != '':
                data_extent = idx
        header = header[:data_extent+1]

        # Basic validation to 1) fail early before processing entire file; and
        # 2) make some basic guarantees about the header for things in this
        # class that use the header as part of reading the file.
        column_names = set(header)
        if '' in column_names:
            raise MetadataFileError(
                "Found at least one column without a name in the header. Each "
                "column must be named.")
        elif len(header) != len(column_names):
            duplicates = find_duplicates(header)
            raise MetadataFileError(
                "Column names must be unique. The following column name(s) "
                "are duplicated: %s" %
                (', '.join(repr(e) for e in sorted(duplicates))))

        return header

    def _read_directives(self, header):
        directives = {}
        for record in self._reader:
            if not self._is_directive(record):
                self._reader = itertools.chain([record], self._reader)
                break

            if not self._is_column_types_directive(record):
                raise MetadataFileError(
                        "Unrecognized directive %r. Only the #q2:types "
                        "directive is supported at this time." % record[0])
            if 'types' in directives:
                raise MetadataFileError(
                    "Found duplicate directive %r. Each directive may "
                    "only be specified a single time." % record[0])

            record = self._match_header_len(record, header)

            column_types = {}
            for column_name, column_type in zip(header[1:], record[1:]):
                if column_type:
                    type_nocase = column_type.lower()
                    if type_nocase in self._supported_column_types:
                        column_types[column_name] = type_nocase
                    else:
                        supported_column_types = ', '.join(
                                sorted(self._supported_column_types))
                        raise MetadataFileError(
                            "Column %r has unrecognized column type %r "
                            "specified in its #q2:types directive. "
                            "Supported column types (case-insensitive): %s"
                            % (column_name, column_type,
                               supported_column_types))
            directives['types'] = column_types
        return directives

    def _read_data(self, header):
        ids = []
        data = []
        for record in self._reader:
            if self._is_comment(record):
                continue
            elif self._is_empty(record):
                continue
            elif self._is_directive(record):
                raise MetadataFileError(
                    "Found directive %r outside of the directives section of "
                    "the file. Directives may only appear immediately after "
                    "the header." % record[0])
            elif self._is_header(record):
                raise MetadataFileError(
                    "Detected metadata ID that conflicts with a name reserved "
                    "for headers: %r" % record[0])

            record = self._match_header_len(record, header)
            ids.append(record[0])
            data.append(record[1:])
        return ids, data

    def _strip_cell_whitespace(self, record):
        return [cell.strip() for cell in record]

    def _match_header_len(self, record, header):
        record_len = len(record)
        header_len = len(header)

        if record_len < header_len:
            # Pad record with empty cells to match header length.
            record = record + [''] * (header_len - record_len)
        elif record_len > header_len:
            trailing_record = record[header_len:]
            if not self._is_empty_record(trailing_record):
                raise MetadataFileError(
                    "Metadata record contains more fields than are declared "
                    "by the header. The record has %d field(s), and the "
                    "header declares %d field(s)."
                    % (record_len, header_len))
            record = record[:header_len]
        return record

    def _is_empty(self, record):
        # `all` returns True for an empty iterable, so this check works for a
        # record of zero elements (corresponds to a blank line in the file).
        return all((cell == '' for cell in record))

    def _is_comment(self, record):
        return (
            len(record) > 0 and
            record[0].startswith('#') and
            not self._is_directive(record) and
            not self._is_header(record)
        )

    def _is_header(self, record):
        if len(record) == 0:
            return False

        # TODO call predicate when it exists
        try:
            self._assert_valid_id_header(record[0])
        except MetadataFileError:
            return False
        else:
            return True

    def _is_directive(self, record):
        return len(record) > 0 and record[0].startswith('#q2:')

    def _is_column_types_directive(self, record):
        return len(record) > 0 and record[0] == '#q2:types'

    # TODO add a predicate version of this, and then have this method call the
    # predicate
    # TODO this code is duplicated from metadata.py, refactor to avoid that
    def _assert_valid_id_header(self, name):
        case_insensitive = {'id', 'sampleid', 'sample id', 'sample-id',
                            'featureid', 'feature id', 'feature-id'}

        # For backwards-compatibility with existing formats.
        exact_match = {
            # QIIME 1 mapping files. "#Sample ID" was never supported, but
            # we're including it here for symmetry with the other supported
            # headers that allow a space between words.
            '#SampleID', '#Sample ID',

            # biom-format: observation metadata and "classic" (TSV) OTU tables.
            '#OTUID', '#OTU ID',

            # Qiita sample/prep information files.
            'sample_name'
        }

        if not (name and
                (name in exact_match or name.lower() in case_insensitive)):
            raise MetadataFileError(
                "Metadata ID header must be one of the following values, not "
                "%r:\n\ncase-insensitive: %s\n\nexact match: %s" %
                (name,
                 ', '.join(sorted(case_insensitive)),
                 ', '.join(sorted(exact_match))))

    def _cast_column(self, series, column_types):
        if series.name in column_types:
            if column_types[series.name] == 'numeric':
                return self._to_numeric(series)
            else:  # 'categorical'
                return series
        else:
            # Infer type
            try:
                return self._to_numeric(series)
            except MetadataFileError:
                return series

    def _to_numeric(self, series):
        is_numeric = series.apply(self._is_numeric)
        if is_numeric.all():
            return pd.to_numeric(series, errors='raise')
        else:
            non_numerics = series[~is_numeric].unique()
            raise MetadataFileError(
                "Cannot convert column %r to numeric. The following value(s) "
                "could not be interpreted as numeric: %s" %
                (series.name,
                 ', '.join(repr(e) for e in sorted(non_numerics))))

    def _is_numeric(self, value):
        return (isinstance(value, float) or
                len(_numeric_regex.findall(value)) == 1)


class MetadataWriter:
    def __init__(self, metadata):
        self._metadata = metadata

    def write(self, filepath):
        # Newline settings based on recommendation from csv docs:
        #     https://docs.python.org/3/library/csv.html#id3
        with open(filepath, 'w', newline='', encoding='utf-8') as fh:
            tsv_writer = csv.writer(fh, dialect='excel-tab', strict=True)

            md = self._metadata
            header = [md.id_header]
            types_directive = ['#q2:types']
            for name, props in md.columns.items():
                header.append(name)
                types_directive.append(props.type)
            tsv_writer.writerow(header)
            tsv_writer.writerow(types_directive)

            df = md.to_dataframe()
            df.fillna('', inplace=True)
            # TODO if writing floats, don't include trailing `.0` so that those
            # numbers look like ints
            df = df.applymap(lambda e: e if isinstance(e, str) else repr(e))
            tsv_writer.writerows(df.itertuples(index=True))


# Credit: https://stackoverflow.com/a/4703508/3776794
_numeric_pattern = r"""
    ^[-+]? # optional sign
    (?:
        (?: \d* \. \d+ ) # .1 .12 .123 etc 9.1 etc 98.1 etc
        |
        (?: \d+ \.? ) # 1. 12. 123. etc 1 12 123 etc
    )
    # followed by optional exponent part if desired
    (?: [Ee] [+-]? \d+ ) ?$
"""

_numeric_regex = re.compile(_numeric_pattern, re.VERBOSE)
