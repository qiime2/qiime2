# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import abc
import collections
import itertools
import sqlite3
import types

import pandas as pd
import numpy as np

from qiime2.core.util import find_duplicates
from .io import MetadataReader, MetadataWriter


class _MetadataBase:
    @property
    def id_header(self):
        return self._id_header

    @property
    def ids(self):
        return self._ids

    @property
    def id_count(self):
        return len(self._ids)

    @property
    def artifacts(self):
        return tuple(self._artifacts)

    def __init__(self, index):
        if index.empty:
            raise ValueError(
                "%s must contain at least one ID and cannot be empty." %
                self.__class__.__name__)

        id_header = index.name
        self._assert_valid_id_header(id_header)
        self._id_header = id_header

        self._validate_pandas_index(index, 'ID')
        self._ids = tuple(index)

        self._artifacts = []

    def _add_artifacts(self, artifacts):
        import qiime2

        deduped = set(self._artifacts)
        for artifact in artifacts:
            if not isinstance(artifact, qiime2.Artifact):
                raise TypeError(
                    "Expected Artifact object, received %r" % artifact)
            if artifact in deduped:
                raise ValueError(
                    "Duplicate source artifacts are not supported on %s "
                    "objects. The following artifact is a duplicate of "
                    "another source artifact: %r" %
                    (self.__class__.__name__, artifact))
            deduped.add(artifact)
        self._artifacts.extend(artifacts)

    @classmethod
    def _assert_valid_id_header(cls, name):
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
            # TODO reference `df.index.name` since this error will only happen
            # at the API level (loading the metadata file should provide its
            # own error message)
            # TODO this code overlaps with similar checks in MetadataReader --
            # refactor to avoid duplication
            raise ValueError(
                "Metadata ID header must be one of the following values, not "
                "%r:\n\ncase-insensitive: %s\n\nexact match: %s" %
                (name,
                 ', '.join(sorted(case_insensitive)),
                 ', '.join(sorted(exact_match))))

    # TODO this only requires an iterable, rename
    def _validate_pandas_index(self, index, label):
        for value in index:
            # TODO raise a better error message for "missing values"
            # (e.g. np.nan, None), right now users will get a "non-string
            # metadata ID" error message, which isn't the most intuitive.
            if not isinstance(value, str):
                raise TypeError(
                    "Detected non-string metadata %s: %r" % (label, value))

            if not value:
                raise ValueError(
                    "Detected empty metadata %s: %r" % (label, value))

            if value != value.strip():
                raise ValueError(
                    "Detected metadata %s with leading or trailing "
                    "whitespace characters: %r" % (label, value))

            # HACK: don't use label as a conditional here
            if label == 'ID' and value.startswith('#'):
                raise ValueError(
                    "Detected metadata %s that begins with the pound sign "
                    "(#): %r" % (label, value))

            try:
                self._assert_valid_id_header(value)
            except ValueError:
                pass
            else:
                raise ValueError(
                    "Detected metadata %s that conflicts with a name reserved "
                    "for ID headers: %r" % (label, value))

        if len(index) != len(set(index)):
            duplicates = find_duplicates(index)
            raise ValueError(
                "Metadata %ss must be unique. The following %ss are "
                "duplicated: %s" %
                (label, label, ', '.join(repr(e) for e in sorted(duplicates))))


# Other properties such as units can be included here in the future!
ColumnProperties = collections.namedtuple('ColumnProperties', ['type'])


class Metadata(_MetadataBase):
    _supported_column_types = {'categorical', 'numeric'}

    @classmethod
    def load(cls, filepath, column_types=None):
        """Load a TSV metadata file.

        Parameters
        ----------
        filepath : str
            Path to TSV metadata file to be loaded.
        column_types : dict, optional
            Override metadata column types specified or inferred in the file.
            This is a dict mapping column names (str) to column types (str).
            Valid column types are 'categorical' and 'numeric'. Column names
            may be omitted from this dict to use the column types read from the
            file.

        Returns
        -------
        Metadata
            Metadata object loaded from `filepath`.

        Raises
        ------
        MetadataFileError
            If the metadata file is invalid in any way (e.g. doesn't meet the
            file format's requirements).

        """
        return MetadataReader(filepath).read(into=cls,
                                             column_types=column_types)

    @classmethod
    def from_artifact(cls, artifact):
        """
        Parameters
        ----------
        artifact : qiime2.Artifact
           Loaded artifact object.

        Returns
        -------
        Metadata

        """
        if not artifact.has_metadata():
            raise ValueError('Artifact has no metadata.')

        return artifact.view(cls)

    @property
    def columns(self):
        """Ordered mapping of column names to ColumnProperties.

        Mapping is read-only.

        """
        # Read-only proxy to the OrderedDict mapping column names to
        # ColumnProperties.
        return types.MappingProxyType(self._columns)

    @property
    def column_count(self):
        """Number of metadata columns.

        Notes
        -----
        Zero metadata columns are allowed.

        """
        return len(self._columns)

    def __init__(self, dataframe):
        if not isinstance(dataframe, pd.DataFrame):
            raise TypeError(
                "%s constructor requires a pandas.DataFrame object, not "
                "%r" % (self.__class__.__name__, type(dataframe)))

        super().__init__(dataframe.index)

        self._dataframe, self._columns = self._normalize_dataframe(dataframe)

    def _normalize_dataframe(self, dataframe):
        self._validate_pandas_index(dataframe.columns, 'column name')

        norm_df = dataframe.copy()
        columns = collections.OrderedDict()
        for column_name, series in dataframe.items():
            metadata_column = self._metadata_column_factory(series)
            norm_df[column_name] = metadata_column.to_series()
            properties = ColumnProperties(type=metadata_column.type)
            columns[column_name] = properties

        return norm_df, columns

    def _metadata_column_factory(self, series):
        dtype = series.dtype
        if NumericMetadataColumn._is_supported_dtype(dtype):
            column = NumericMetadataColumn(series)
        elif CategoricalMetadataColumn._is_supported_dtype(dtype):
            column = CategoricalMetadataColumn(series)
        else:
            # TODO list supported dtypes
            raise TypeError(
                "Metadata column %r has an unsupported pandas dtype %r" %
                (series.name, dtype))

        column._add_artifacts(self.artifacts)
        return column

    def __repr__(self):
        lines = []

        # Header
        lines.append(self.__class__.__name__)
        lines.append('-' * len(self.__class__.__name__))

        # Dimensions
        lines.append('%d ID%s x %d column%s' % (
            self.id_count,
            '' if self.id_count == 1 else 's',
            self.column_count,
            '' if self.column_count == 1 else 's',
        ))

        # Column properties
        max_name_len = max((len(name) for name in self.columns))
        for name, props in self.columns.items():
            padding = ' ' * ((max_name_len - len(name)) + 1)
            lines.append('%s:%s%r' % (name, padding, props))

        # Epilogue
        lines.append('')
        lines.append('Call to_dataframe() for a tabular representation.')

        return '\n'.join(lines)

    # TODO defer some of this to base class __eq__?
    def __eq__(self, other):
        return (
            isinstance(other, self.__class__) and
            self._artifacts == other._artifacts and
            self._columns == other._columns and
            self._dataframe.equals(other._dataframe)
        )

    def __ne__(self, other):
        return not (self == other)

    def get_column(self, name):
        """Retrieve metadata column based on column name.

        Parameters
        ----------
        name : str
            Metadata column name to retrieve.

        Returns
        -------
        MetadataColumn
            Requested metadata column.

        """
        try:
            series = self._dataframe[name]
        except KeyError:
            raise ValueError(
                '%s is not a column in the metadata. Available columns: '
                '%s' % (name, ', '.join(self.columns)))

        return self._metadata_column_factory(series)

    def save(self, filepath):
        MetadataWriter(self).write(filepath)

    def to_dataframe(self):
        return self._dataframe.copy()

    def merge(self, *others):
        """Merge this ``Metadata`` object with other ``Metadata`` objects.

        Returns a new ``Metadata`` object containing the merged contents of
        this ``Metadata`` object and `others`. The merge is not in-place and
        will always return a **new** merged ``Metadata`` object.

        The merge will include only those IDs that are shared across **all**
        ``Metadata`` objects being merged (i.e. the merge is an *inner join*).

        Each metadata column being merged must be unique; merging metadata with
        overlapping columns will result in an error.

        Parameters
        ----------
        others : tuple
            Zero or more ``Metadata`` objects to merge with this ``Metadata``
            object.

        Returns
        -------
        Metadata
            New object containing merged metadata. The merged IDs will be in
            the same relative order as the IDs in this ``Metadata`` object
            after performing the inner join. The merged column order will match
            the column order of ``Metadata`` objects being merged from left to
            right.

        Notes
        -----
        The merged metadata object tracks all source artifacts that it was
        built from to preserve provenance (i.e. the ``.artifacts`` property
        on all ``Metadata`` objects is merged).

        """
        dfs = []
        columns = []
        artifacts = []
        for md in itertools.chain([self], others):
            df = md._dataframe
            dfs.append(df)
            columns.extend(df.columns.tolist())
            artifacts.extend(md.artifacts)

        columns = pd.Index(columns)
        if columns.has_duplicates:
            raise ValueError(
                "Cannot merge metadata with overlapping columns. The "
                "following columns overlap: %s" %
                ', '.join([repr(e) for e in columns.get_duplicates()]))

        merged_df = dfs[0].join(dfs[1:], how='inner')

        # Not using DataFrame.empty because empty columns are allowed in
        # Metadata.
        if merged_df.index.empty:
            raise ValueError(
                "Cannot merge because there are no IDs shared across metadata "
                "objects.")

        merged_df.index.name = 'id'
        merged_md = self.__class__(merged_df)
        merged_md._add_artifacts(artifacts)
        return merged_md

    def get_ids(self, where=None):
        """Retrieve IDs matching search criteria.

        Parameters
        ----------
        where : str, optional
            SQLite WHERE clause specifying criteria IDs must meet to be
            included in the results. All IDs are included by default.

        Returns
        -------
        set
            IDs matching search criteria specified in `where`.

        See Also
        --------
        ids

        """
        if where is None:
            return set(self._ids)

        conn = sqlite3.connect(':memory:')
        conn.row_factory = lambda cursor, row: row[0]

        self._dataframe.to_sql('metadata', conn, index=True,
                               index_label=self.id_header)
        c = conn.cursor()

        # In general we wouldn't want to format our query in this way because
        # it leaves us open to sql injection, but it seems acceptable here for
        # a few reasons:
        # 1) This is a throw-away database which we're just creating to have
        #    access to the query language, so any malicious behavior wouldn't
        #    impact any data that isn't temporary
        # 2) The substitution syntax recommended in the docs doesn't allow
        #    us to specify complex `where` statements, which is what we need to
        #    do here. For example, we need to specify things like:
        #        WHERE Subject='subject-1' AND SampleType='gut'
        #    but their qmark/named-style syntaxes only supports substition of
        #    variables, such as:
        #        WHERE Subject=?
        # 3) sqlite3.Cursor.execute will only execute a single statement so
        #    inserting multiple statements
        #    (e.g., "Subject='subject-1'; DROP...") will result in an
        #    OperationalError being raised.
        query = ('SELECT "{0}" FROM metadata WHERE {1} GROUP BY "{0}" '
                 'ORDER BY "{0}";'.format(self.id_header, where))

        try:
            c.execute(query)
        except sqlite3.OperationalError:
            conn.close()
            raise ValueError("Selection of IDs failed with query:\n %s"
                             % query)

        ids = set(c.fetchall())
        conn.close()
        return ids

    def filter_ids(self, ids_to_keep):
        """Filter metadata by IDs.

        Parameters
        ----------
        ids_to_keep : iterable of str
            IDs that should be retained in the filtered ``Metadata`` object. If
            any IDs in `ids_to_keep` are not contained in this ``Metadata``
            object, a ``ValueError`` will be raised. The filtered ``Metadata``
            object will retain the same relative ordering of IDs in this
            ``Metadata`` object. Thus, the ordering of IDs in `ids_to_keep`
            does not determine the ordering of IDs in the filtered ``Metadata``
            object.

        Returns
        -------
        Metadata
            The metadata filtered by IDs.

        """
        ids_to_keep_set = set(ids_to_keep)
        if len(ids_to_keep) != len(ids_to_keep_set):
            duplicates = find_duplicates(ids_to_keep)
            raise ValueError(
                "`ids_to_keep` must consist of unique IDs. The following "
                "ID(s) are duplicated: %s"
                % (', '.join(repr(e) for e in sorted(duplicates))))

        missing_ids = ids_to_keep_set - self.get_ids()
        if missing_ids:
            raise ValueError(
                "The following ID(s) are not present in the Metadata: %s"
                % (', '.join(repr(e) for e in sorted(missing_ids))))

        # While preserving order, get rid of any IDs not contained in
        # `ids_to_keep`.
        ids_to_discard = self.get_ids() - ids_to_keep_set
        filtered_df = self._dataframe.drop(labels=ids_to_discard, axis='index',
                                           inplace=False, errors='raise')

        # Not using DataFrame.empty because empty columns are allowed in
        # Metadata.
        # TODO instead of erroring here, just check that `ids_to_keep` isn't
        # empty at the start of this method.
        if filtered_df.index.empty:
            raise ValueError(
                "All IDs were filtered out of the Metadata, resulting in an "
                "empty Metadata object.")

        filtered_md = self.__class__(filtered_df)
        filtered_md._add_artifacts(self.artifacts)
        return filtered_md

    def filter_columns(self, *, column_type=None, drop_all_unique=False,
                       drop_zero_variance=False, drop_all_missing=False):
        """Filter metadata by columns.

        Parameters
        ----------
        column_type : str, optional
            If supplied, will retain only columns of this type. The currently
            supported column types are 'numeric' and 'categorical'.
        drop_all_unique : bool, optional
            If True, columns that contain a unique value for every ID will be
            dropped.
        drop_zero_variance : bool, optional
            If True, columns that contain the same value for every ID will be
            dropped.
        drop_all_missing : bool, optional
            If True, columns that have a missing value for every ID will be
            dropped.

        Returns
        -------
        Metadata
            The metadata filtered by columns.

        """
        df = self._dataframe

        if column_type is not None:
            df = self._filter_columns_by_type(df, column_type)

        if drop_all_unique or drop_zero_variance:
            df = self._filter_columns_by_variance(df, drop_all_unique,
                                                  drop_zero_variance)
        if drop_all_missing:
            df = self._filter_empty_columns(df)

        filtered_md = self.__class__(df)
        filtered_md._add_artifacts(self.artifacts)
        return filtered_md

    def _filter_columns_by_type(self, df, column_type):
        if column_type not in self._supported_column_types:
            raise ValueError(
                "Unknown column type %r. Supported column types: %s" %
                (column_type, ', '.join(sorted(self._supported_column_types))))

        columns_to_keep = [name for name in self.columns
                           if self.columns[name].type == column_type]
        # TODO use `drop` instead, or does this preserve relative sort order?
        # What about filtering to nothing, or filtering nothing?
        return df.filter(items=columns_to_keep, axis=1)

    def _filter_columns_by_variance(self, df, drop_all_unique=False,
                                    drop_zero_variance=False):
        num_ids = df.shape[0]
        all_unique = []
        zero_variance = []
        for column in df.columns:
            # TODO handle when all nans (nunique=0)
            num_unique = df[column].nunique()
            if num_unique == num_ids:
                all_unique.append(column)
            elif num_unique == 1:
                zero_variance.append(column)

        # TODO handle case where nunique==1 and nunique==num_ids
        if drop_all_unique:
            df = df.drop(all_unique, axis=1)

        if drop_zero_variance:
            df = df.drop(zero_variance, axis=1)

        return df

    def _filter_empty_columns(self, df):
        return df.dropna(axis='columns', how='all', inplace=False)


class MetadataColumn(_MetadataBase, metaclass=abc.ABCMeta):
    # Abstract, must be defined by subclasses.
    type = None

    @classmethod
    @abc.abstractmethod
    def _is_supported_dtype(cls, dtype):
        pass

    @classmethod
    @abc.abstractmethod
    def _normalize_(cls, series):
        """

        Contract: Return a copy of `series` that has been converted to the
        appropriate internal dtype and has any other necessary normalization or
        validation applied (e.g. missing value representations, disallowing
        certain values, etc). Raise an error with a detailed error message if
        the operation cannot be completed.

        """
        pass

    @property
    def name(self):
        return self._series.name

    def __init__(self, series):
        if not isinstance(series, pd.Series):
            raise TypeError(
                "%s constructor requires a pandas.Series object, not %r" %
                (self.__class__.__name__, type(series)))

        super().__init__(series.index)

        # TODO make this error msg clearer, i.e. similar to setting
        # df.index.name, the Series.name must be set
        self._validate_pandas_index([series.name], 'column name')

        if not self._is_supported_dtype(series.dtype):
            raise TypeError(
                "%s does not support pandas dtype %r" %
                (self.__class__.__name__, series.dtype))

        self._series = self._normalize_(series)

    def __repr__(self):
        return '<%s name=%r id_count=%d>' % (self.__class__.__name__,
                                             self.name, self.id_count)

    def save(self, filepath):
        # TODO consider adding a `to_metadata()` method, which this method
        # could call instead of `to_dataframe()`. Converting a MetadataColumn
        # into a Metadata may become more complicated in the future (e.g.
        # passing other column directives through, source artifacts, etc.).
        metadata = Metadata(self.to_dataframe())
        MetadataWriter(metadata).write(filepath)

    def to_series(self):
        return self._series.copy()

    def to_dataframe(self):
        return self._series.to_frame()

    def get_value(self, id):
        if id not in self._series.index:
            raise ValueError(
                "ID %r is not present in the MetadataColumn." % id)
        return self._series.loc[id]

    def has_missing_values(self):
        return len(self.get_ids(where_values_missing=True)) > 0

    def drop_missing_values(self):
        missing = self.get_ids(where_values_missing=True)
        present = self.get_ids() - missing
        return self.filter_ids(present)

    # TODO get_ids() and filter_ids() are mostly duplicated from Metadata,
    # refactor to avoid this.
    def get_ids(self, where_values_missing=False):
        """Retrieve IDs matching search criteria.

        Parameters
        ----------
        where_values_missing : bool, optional
            If ``True``, only return IDs that have missing values. If ``False``
            (the default), return all IDs.

        Returns
        -------
        set
            IDs matching search criteria.

        See Also
        --------
        ids

        """
        if where_values_missing:
            ids = self._series[self._series.isnull()].index
        else:
            ids = self._ids
        return set(ids)

    def filter_ids(self, ids_to_keep):
        ids_to_keep_set = set(ids_to_keep)
        if len(ids_to_keep) != len(ids_to_keep_set):
            duplicates = find_duplicates(ids_to_keep)
            raise ValueError(
                "`ids_to_keep` must consist of unique IDs. The following "
                "ID(s) are duplicated: %s"
                % (', '.join(repr(e) for e in sorted(duplicates))))

        missing_ids = ids_to_keep_set - self.get_ids()
        if missing_ids:
            raise ValueError(
                "The following ID(s) are not present in the MetadataColumn: %s"
                % (', '.join(repr(e) for e in sorted(missing_ids))))

        # While preserving order, get rid of any IDs not contained in
        # `ids_to_keep`.
        ids_to_discard = self.get_ids() - ids_to_keep_set
        filtered_series = self._series.drop(
            labels=ids_to_discard, axis='index', inplace=False, errors='raise')

        # Not using Series.empty because empty columns are allowed in
        # Metadata.
        # TODO instead of erroring here, just check that `ids_to_keep` isn't
        # empty at the start of this method.
        if filtered_series.index.empty:
            raise ValueError(
                "All IDs were filtered out of the MetadataColumn, resulting "
                "in an empty MetadataColumn object.")

        filtered_mdc = self.__class__(filtered_series)
        filtered_mdc._add_artifacts(self.artifacts)
        return filtered_mdc


class CategoricalMetadataColumn(MetadataColumn):
    type = 'categorical'

    @classmethod
    def _is_supported_dtype(cls, dtype):
        return dtype == 'object'

    @classmethod
    def _normalize_(cls, series):
        def normalize(value):
            if isinstance(value, str):
                if value == '':
                    raise ValueError(
                        "%s does not support empty strings. Use an "
                        "appropriate pandas missing value type "
                        "(e.g. `numpy.nan`) or supply a non-empty string as "
                        "the value." % cls.__name__)
                elif value != value.strip():
                    raise ValueError(
                        "%s does not support strings with leading or trailing "
                        "whitespace characters: %r" % (cls.__name__, value))
                else:
                    return value
            elif pd.isnull(value):  # permits np.nan, Python float nan, None
                return np.nan
            else:
                raise TypeError(
                    "%s only supports strings or missing values. Found value "
                    "%r with type %r" % (cls.__name__, value, type(value)))

        return series.apply(normalize, convert_dtype=False)


class NumericMetadataColumn(MetadataColumn):
    type = 'numeric'

    @classmethod
    def _is_supported_dtype(cls, dtype):
        return dtype == 'float' or dtype == 'int'

    @classmethod
    def _normalize_(cls, series):
        series = series.astype(float, copy=True, errors='raise')
        if np.isinf(series).any():
            raise ValueError(
                "%s does not support positive or negative infinity as a "
                "floating point value." % cls.__name__)
        return series
