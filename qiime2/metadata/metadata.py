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

import qiime2
from qiime2.core.util import find_duplicates
from .base import SUPPORTED_COLUMN_TYPES, FORMATTED_ID_HEADERS, is_id_header


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
                "%s must contain at least one ID." % self.__class__.__name__)

        id_header = index.name
        self._assert_valid_id_header(id_header)
        self._id_header = id_header

        self._validate_index(index, axis='id')
        self._ids = tuple(index)

        self._artifacts = []

    def __eq__(self, other):
        return (
            isinstance(other, self.__class__) and
            self._id_header == other._id_header and
            self._artifacts == other._artifacts
        )

    def __ne__(self, other):
        return not (self == other)

    def _add_artifacts(self, artifacts):
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

    # Static helpers below for code reuse in Metadata and MetadataColumn

    @classmethod
    def _assert_valid_id_header(cls, name):
        if not is_id_header(name):
            raise ValueError(
                "pandas index name (`Index.name`) must be one of the "
                "following values, not %r:\n\n%s" %
                (name, FORMATTED_ID_HEADERS))

    @classmethod
    def _validate_index(cls, index, *, axis):
        if axis == 'id':
            label = 'ID'
        elif axis == 'column':
            label = 'column name'
        else:
            raise NotImplementedError

        for value in index:
            if not isinstance(value, str):
                raise TypeError(
                    "Detected non-string metadata %s of type %r: %r" %
                    (label, type(value), value))

            if not value:
                raise ValueError(
                    "Detected empty metadata %s. %ss must consist of at least "
                    "one character." % (label, label))

            if value != value.strip():
                raise ValueError(
                    "Detected metadata %s with leading or trailing "
                    "whitespace characters: %r" % (label, value))

            if axis == 'id' and value.startswith('#'):
                raise ValueError(
                    "Detected metadata %s that begins with a pound sign "
                    "(#): %r" % (label, value))

            if is_id_header(value):
                raise ValueError(
                    "Detected metadata %s %r that conflicts with a name "
                    "reserved for the ID header. Reserved ID headers:\n\n%s" %
                    (label, value, FORMATTED_ID_HEADERS))

        if len(index) != len(set(index)):
            duplicates = find_duplicates(index)
            raise ValueError(
                "Metadata %ss must be unique. The following %ss are "
                "duplicated: %s" %
                (label, label, ', '.join(repr(e) for e in sorted(duplicates))))

    @classmethod
    def _filter_ids_helper(cls, df_or_series, ids, ids_to_keep):
        # `ids_to_keep` can be any iterable, so turn it into a list so that it
        # can be iterated over multiple times below (and length-checked).
        ids_to_keep = list(ids_to_keep)

        if len(ids_to_keep) == 0:
            raise ValueError("`ids_to_keep` must contain at least one ID.")

        duplicates = find_duplicates(ids_to_keep)
        if duplicates:
            raise ValueError(
                "`ids_to_keep` must contain unique IDs. The following IDs are "
                "duplicated: %s" %
                (', '.join(repr(e) for e in sorted(duplicates))))

        ids_to_keep = set(ids_to_keep)
        missing_ids = ids_to_keep - ids
        if missing_ids:
            raise ValueError(
                "The following IDs are not present in the metadata: %s"
                % (', '.join(repr(e) for e in sorted(missing_ids))))

        # While preserving order, get rid of any IDs not contained in
        # `ids_to_keep`.
        ids_to_discard = ids - ids_to_keep
        return df_or_series.drop(labels=ids_to_discard, axis='index',
                                 inplace=False, errors='raise')


# Other properties such as units can be included here in the future!
ColumnProperties = collections.namedtuple('ColumnProperties', ['type'])


class Metadata(_MetadataBase):
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
        from .io import MetadataReader
        return MetadataReader(filepath).read(into=cls,
                                             column_types=column_types)

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
        self._validate_index(dataframe.columns, axis='column')

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
            raise TypeError(
                "Metadata column %r has an unsupported pandas dtype of %s. "
                "Supported dtypes: float, int, object" %
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

    def __eq__(self, other):
        return (
            super().__eq__(other) and
            self._columns == other._columns and
            self._dataframe.equals(other._dataframe)
        )

    def __ne__(self, other):
        return not (self == other)

    def save(self, filepath):
        from .io import MetadataWriter
        MetadataWriter(self).write(filepath)

    def to_dataframe(self):
        return self._dataframe.copy()

    def get_column(self, name):
        """Retrieve metadata column based on column name.

        Parameters
        ----------
        name : str
            Name of the metadata column to retrieve.

        Returns
        -------
        MetadataColumn
            Requested metadata column.

        """
        try:
            series = self._dataframe[name]
        except KeyError:
            raise ValueError(
                '%r is not a column in the metadata. Available columns: '
                '%s' % (name, ', '.join(repr(c) for c in self.columns)))

        return self._metadata_column_factory(series)

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
        except sqlite3.OperationalError as e:
            conn.close()
            raise ValueError("Selection of IDs failed with query:\n %s\n\n"
                             "If one of the metadata column names specified "
                             "in the `where` statement is on this list "
                             "of reserved keywords "
                             "(http://www.sqlite.org/lang_keywords.html), "
                             "please ensure it is quoted appropriately in the "
                             "`where` statement." % query) from e

        ids = set(c.fetchall())
        conn.close()
        return ids

    def merge(self, *others):
        """Merge this ``Metadata`` object with other ``Metadata`` objects.

        Returns a new ``Metadata`` object containing the merged contents of
        this ``Metadata`` object and `others`. The merge is not in-place and
        will always return a **new** merged ``Metadata`` object.

        The merge will include only those IDs that are shared across **all**
        ``Metadata`` objects being merged (i.e. the merge is an *inner join*).

        Each metadata column being merged must have a unique name; merging
        metadata with overlapping column names will result in an error.

        Parameters
        ----------
        others : tuple
            One or more ``Metadata`` objects to merge with this ``Metadata``
            object.

        Returns
        -------
        Metadata
            New object containing merged metadata. The merged IDs will be in
            the same relative order as the IDs in this ``Metadata`` object
            after performing the inner join. The merged column order will match
            the column order of ``Metadata`` objects being merged from left to
            right.

        Raises
        ------
        ValueError
            If zero ``Metadata`` objects are provided in `others` (there is
            nothing to merge in this case).

        Notes
        -----
        The merged ``Metadata`` object will always have its ``id_header``
        property set to ``'id'``, regardless of the ``id_header`` values on the
        ``Metadata`` objects being merged.

        The merged ``Metadata`` object tracks all source artifacts that it was
        built from to preserve provenance (i.e. the ``.artifacts`` property
        on all ``Metadata`` objects is merged).

        """
        if len(others) < 1:
            raise ValueError(
                "At least one Metadata object must be provided to merge into "
                "this Metadata object (otherwise there is nothing to merge).")

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
        filtered_df = self._filter_ids_helper(self._dataframe, self.get_ids(),
                                              ids_to_keep)
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
        if (column_type is not None and
                column_type not in SUPPORTED_COLUMN_TYPES):
            raise ValueError(
                "Unknown column type %r. Supported column types: %s" %
                (column_type, ', '.join(sorted(SUPPORTED_COLUMN_TYPES))))

        # Build up a set of columns to drop. Short-circuit as soon as we know a
        # given column can be dropped (no need to apply further filters to it).
        columns_to_drop = set()
        for column, props in self.columns.items():
            if column_type is not None and props.type != column_type:
                columns_to_drop.add(column)
                continue

            series = self._dataframe[column]
            if drop_all_unique or drop_zero_variance:
                # Ignore nans in the unique count, and compare to the number of
                # non-nan values in the series.
                num_unique = series.nunique(dropna=True)
                if drop_all_unique and num_unique == series.count():
                    columns_to_drop.add(column)
                    continue

                # If num_unique == 0, the series was empty (all nans). If
                # num_unique == 1, the series contained only a single unique
                # value (ignoring nans).
                if drop_zero_variance and num_unique < 2:
                    columns_to_drop.add(column)
                    continue

            if drop_all_missing and series.isnull().all():
                columns_to_drop.add(column)
                continue

        filtered_df = self._dataframe.drop(columns_to_drop, axis=1,
                                           inplace=False)
        filtered_md = self.__class__(filtered_df)
        filtered_md._add_artifacts(self.artifacts)
        return filtered_md


class MetadataColumn(_MetadataBase, metaclass=abc.ABCMeta):
    # Abstract, must be defined by subclasses.
    type = None

    @classmethod
    @abc.abstractmethod
    def _is_supported_dtype(cls, dtype):
        raise NotImplementedError

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
        raise NotImplementedError

    @property
    def name(self):
        return self._series.name

    def __init__(self, series):
        if not isinstance(series, pd.Series):
            raise TypeError(
                "%s constructor requires a pandas.Series object, not %r" %
                (self.__class__.__name__, type(series)))

        super().__init__(series.index)

        self._validate_index([series.name], axis='column')

        if not self._is_supported_dtype(series.dtype):
            raise TypeError(
                "%s %r does not support a pandas.Series object with dtype %s" %
                (self.__class__.__name__, series.name, series.dtype))

        self._series = self._normalize_(series)

    def __repr__(self):
        return '<%s name=%r id_count=%d>' % (self.__class__.__name__,
                                             self.name, self.id_count)

    def __eq__(self, other):
        return (
            super().__eq__(other) and
            self.name == other.name and
            self._series.equals(other._series)
        )

    def __ne__(self, other):
        return not (self == other)

    def save(self, filepath):
        from .io import MetadataWriter
        MetadataWriter(self).write(filepath)

    def to_series(self):
        return self._series.copy()

    def to_dataframe(self):
        return self._series.to_frame()

    def get_value(self, id):
        if id not in self._series.index:
            raise ValueError("ID %r is not present in %r" % (id, self))
        return self._series.loc[id]

    def has_missing_values(self):
        return len(self.get_ids(where_values_missing=True)) > 0

    def drop_missing_values(self):
        missing = self.get_ids(where_values_missing=True)
        present = self.get_ids() - missing
        return self.filter_ids(present)

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
        filtered_series = self._filter_ids_helper(self._series, self.get_ids(),
                                                  ids_to_keep)
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
                        "%s does not support empty strings as values. Use an "
                        "appropriate pandas missing value type "
                        "(e.g. `numpy.nan`) or supply a non-empty string as "
                        "the value in column %r." %
                        (cls.__name__, series.name))
                elif value != value.strip():
                    raise ValueError(
                        "%s does not support values with leading or trailing "
                        "whitespace characters. Column %r has the following "
                        "value: %r" % (cls.__name__, series.name, value))
                else:
                    return value
            elif pd.isnull(value):  # permits np.nan, Python float nan, None
                return np.nan
            else:
                raise TypeError(
                    "%s only supports strings or missing values. Found value "
                    "%r of type %r in column %r." %
                    (cls.__name__, value, type(value), series.name))

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
                "floating point value in column %r." %
                (cls.__name__, series.name))
        return series
