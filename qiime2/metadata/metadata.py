# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import abc
import collections
import uuid
import tempfile
import copy

import sqlite3
import types
import warnings

import pandas as pd
import numpy as np

import qiime2
from qiime2.core.util import find_duplicates, md5sum
from .base import SUPPORTED_COLUMN_TYPES, FORMATTED_ID_HEADERS, is_id_header


class _MetadataBase:
    """Base class for functionality shared between Metadata and MetadataColumn.

    Parameters
    ----------
    index : pandas.Index
        IDs associated with the metadata.

    """

    @property
    def id_header(self):
        """Name identifying the IDs associated with the metadata.

        This property is read-only.

        Returns
        -------
        str
            Name of IDs associated with the metadata.

        """
        return self._id_header

    @property
    def ids(self):
        """IDs associated with the metadata.

        This property is read-only.

        Returns
        -------
        tuple of str
            Metadata IDs.

        """
        return self._ids

    @property
    def id_count(self):
        """Number of metadata IDs.

        This property is read-only.

        Returns
        -------
        int
            Number of metadata IDs.

        """
        return len(self._ids)

    def __init__(self, index):
        if index.empty:
            raise ValueError(
                "%s must contain at least one ID." % self.__class__.__name__)

        id_header = index.name
        self._assert_valid_id_header(id_header)
        self._id_header = id_header

        self._validate_index(index, axis='id')
        self._ids = tuple(index)

    def __eq__(self, other):
        return (
            isinstance(other, self.__class__) and
            self._id_header == other._id_header
        )

    def __ne__(self, other):
        return not (self == other)

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
    """Store metadata associated with identifiers in a study.

    Metadata is tabular in nature, mapping study identifiers (e.g. sample or
    feature IDs) to columns of metadata associated with each ID.

    For more details about metadata in QIIME 2, including the TSV metadata file
    format, see the Metadata Tutorial at https://docs.qiime2.org.

    The following text focuses on design and considerations when working with
    ``Metadata`` objects at the API level.

    A ``Metadata`` object is composed of zero or more ``MetadataColumn``
    objects. A ``Metadata`` object always contains at least one ID, regardless
    of the number of columns. Each column in the ``Metadata`` object has an
    associated column type representing either *categorical* or *numeric*
    data. Each metadata column is represented by an object corresponding to the
    column's type: ``CategoricalMetadataColumn`` or ``NumericMetadataColumn``,
    respectively.

    A ``Metadata`` object is closely linked to its corresponding TSV metadata
    file format described at https://docs.qiime2.org. Therefore, certain
    requirements present in the file format are also enforced on the in-memory
    object in order to make serialized ``Metadata`` objects roundtrippable when
    loaded from disk again. For example, IDs cannot begin with a pound
    character (``#``) because those IDs would be interpreted as comment rows
    when written to disk as TSV. See the metadata file format spec for more
    details about data formatting requirements.

    In addition to being loaded from or saved to disk, a ``Metadata`` object
    can be constructed from a ``pandas.DataFrame`` object. See the *Parameters*
    section below for details on how to construct ``Metadata`` objects from
    dataframes.

    ``Metadata`` objects have various methods to access, filter, and merge
    data. A dataframe can be retrieved from the ``Metadata`` object for further
    data manipulation using the pandas API. Individual ``MetadataColumn``
    objects can be retrieved to gain access to APIs applicable to a single
    metadata column.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        Dataframe containing metadata. The dataframe's index defines the IDs,
        and the index name (``Index.name``) must match one of the required ID
        headers described in the metadata file format spec. Each column in the
        dataframe defines a metadata column, and the metadata column's type
        (i.e. *categorical* or *numeric*) is determined based on the column's
        dtype. If a column has ``dtype=object``, it may contain strings or
        pandas missing values (e.g. ``np.nan``, ``None``). Columns matching
        this requirement are assumed to be *categorical*. If a column in the
        dataframe has ``dtype=float`` or ``dtype=int``, it may contain floating
        point numbers or integers, as well as pandas missing values
        (e.g. ``np.nan``). Columns matching this requirement are assumed to be
        *numeric*. Regardless of column type (categorical vs numeric), the
        dataframe stored within the ``Metadata`` object will have any missing
        values normalized to ``np.nan``. Columns with ``dtype=int`` will be
        cast to ``dtype=float``. To obtain a dataframe from the ``Metadata``
            object containing these normalized data types and values, use
        ``Metadata.to_dataframe()``.

    """

    @classmethod
    def load(cls, filepath, column_types=None):
        """Load a TSV metadata file.

        The TSV metadata file format is described at https://docs.qiime2.org in
        the Metadata Tutorial.

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

        See Also
        --------
        save

        """
        from .io import MetadataReader
        return MetadataReader(filepath).read(into=cls,
                                             column_types=column_types)

    @property
    def columns(self):
        """Ordered mapping of column names to ColumnProperties.

        The mapping that is returned is read-only. This property is also
        read-only.

        Returns
        -------
        types.MappingProxyType
            Ordered mapping of column names to ColumnProperties.

        """
        # Read-only proxy to the OrderedDict mapping column names to
        # ColumnProperties.
        return types.MappingProxyType(self._columns)

    @property
    def column_count(self):
        """Number of metadata columns.

        This property is read-only.

        Returns
        -------
        int
            Number of metadata columns.

        Notes
        -----
        Zero metadata columns are allowed.

        See Also
        --------
        id_count

        """
        return len(self._columns)

    @property
    def artifacts(self):
        """Artifacts that are the source of the metadata.

        This property is read-only.

        Returns
        -------
        tuple of qiime2.Artifact
            Source artifacts of the metadata.

        """
        if self._source_artifact is not None:
            return (self._source_artifact,)

        artifacts = set()
        for source in self._column_sources.values():
            if source is not None:
                artifacts.update(artifact for artifact in source.artifacts)
        return tuple(artifacts)

    @property
    def source_artifact(self):
        """Artifact that is the source of the metadata.

        This property is read-only.

        Returns
        -------
        qiime2.Artifact
            Source artifact of the metadata.

        """
        return self._source_artifact

    # Try to ensure the source artifact is actually an artifact
    def _add_source_artifact(self, artifact):
        if not isinstance(artifact, qiime2.Artifact):
            raise TypeError('Source Artifact must an Artifact object not '
                            f'{artifact}.')
        self._source_artifact = artifact

    def __init__(self, dataframe):
        if not isinstance(dataframe, pd.DataFrame):
            raise TypeError(
                "%s constructor requires a pandas.DataFrame object, not "
                "%r" % (self.__class__.__name__, type(dataframe)))

        super().__init__(dataframe.index)

        self._dataframe, self._columns = self._normalize_dataframe(dataframe)
        self.contains_renamed_columns = False

        self._column_names = {}
        self._column_sources = {}
        for column in self._columns:
            self._column_names[column] = column
            self._column_sources[column] = None

        # This is set post facto at the point of the creation of this metadata
        # from an artifact
        self._source_artifact = None

    def _init(self, column_names, column_sources):
        self._column_names = column_names
        self._column_sources = column_sources

    @property
    def id(self):
        if not isinstance(self._id, uuid.UUID):
            return f'{self._id:x}'
        else:
            return self._id

    @property
    def _id(self):
        if self._source_artifact is not None:
            return self._source_artifact.uuid
        with tempfile.NamedTemporaryFile(prefix='md5-') as fh:
            self.save(fh.name)
            return int(md5sum(fh.name), 16)

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
            column._init(self)
        elif CategoricalMetadataColumn._is_supported_dtype(dtype):
            column = CategoricalMetadataColumn(series)
            column._init(self)
        else:
            raise TypeError(
                "Metadata column %r has an unsupported pandas dtype of %s. "
                "Supported dtypes: float, int, object" %
                (series.name, dtype))

        return column

    def __repr__(self):
        """String summary of the metadata and its columns."""
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
        if self.column_count != 0:
            max_name_len = max((len(name) for name in self.columns))
            for name, props in self.columns.items():
                padding = ' ' * ((max_name_len - len(name)) + 1)
                lines.append('%s:%s%r' % (name, padding, props))

        # Epilogue
        lines.append('')
        lines.append('Call to_dataframe() for a tabular representation.')

        return '\n'.join(lines)

    def __eq__(self, other):
        """Determine if this metadata is equal to another.

        ``Metadata`` objects are equal if their IDs, columns (including column
        names, types, and ordering), ID headers, source artifacts, and metadata
        values are equal.

        Parameters
        ----------
        other : Metadata
            Metadata to test for equality.

        Returns
        -------
        bool
            Indicates whether this ``Metadata`` object is equal to `other`.

        See Also
        --------
        __ne__

        """
        return (
            super().__eq__(other) and
            self._columns == other._columns and
            self._dataframe.equals(other._dataframe) and
            self._id == other._id
        )

    def __ne__(self, other):
        """Determine if this metadata is not equal to another.

        ``Metadata`` objects are not equal if their IDs, columns (including
        column names, types, or ordering), ID headers, source artifacts, or
        metadata values are not equal.

        Parameters
        ----------
        other : Metadata
            Metadata to test for inequality.

        Returns
        -------
        bool
            Indicates whether this ``Metadata`` object is not equal to `other`.

        See Also
        --------
        __eq__

        """
        return not (self == other)

    def save(self, filepath):
        """Save a TSV metadata file.

        The TSV metadata file format is described at https://docs.qiime2.org in
        the Metadata Tutorial.

        The file will always include the ``#q2:types`` directive in order to
        make the file roundtrippable without relying on column type inference.

        Parameters
        ----------
        filepath : str
            Path to save TSV metadata file at.

        See Also
        --------
        load

        """
        from .io import MetadataWriter
        MetadataWriter(self).write(filepath)

    def to_dataframe(self):
        """Create a pandas dataframe from the metadata.

        The dataframe's index name (``Index.name``) will match this metadata
        object's ``id_header``, and the index will contain this metadata
        object's IDs. The dataframe's column names will match the column names
        in this metadata. Categorical columns will be stored as
        ``dtype=object`` (containing strings), and numeric columns will be
        stored as ``dtype=float``.

        Returns
        -------
        pandas.DataFrame
            Dataframe constructed from the metadata.

        """
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
            Requested metadata column (``CategoricalMetadataColumn`` or
            ``NumericMetadataColumn``).

        See Also
        --------
        get_ids

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
        filter_ids
        get_column

        Notes
        -----
        The ID header (``Metadata.id_header``) may be used in the `where`
        clause to query the table's ID column.

        """
        if where is None:
            return set(self._ids)

        conn = sqlite3.connect(':memory:')
        conn.row_factory = lambda cursor, row: row[0]

        # https://github.com/pandas-dev/pandas/blob/
        # 7c7bd569ce8e0f117c618d068e3d2798134dbc73/pandas/io/sql.py#L1306
        with warnings.catch_warnings():
            warnings.filterwarnings(
                'ignore', 'The spaces in these column names will not.*')
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
        column_names = []
        column_sources = []
        columns = {}
        df_changes = {}

        mds = [self]
        mds.extend(others)
        dupes = set()
        for md in mds:
            if md._id in dupes:
                raise ValueError(
                    'Your input contained duplicate metadata files. Merging '
                    'the same file multiple times is not allowed.')
            dupes.add(md._id)

        for i, md in enumerate(mds):
            df = md._dataframe
            dfs.append(df)
            names = copy.deepcopy(md._column_names)
            column_names.append(names)
            sources = copy.deepcopy(md._column_sources)
            column_sources.append(sources)
            items = names.items()

            for new_name, old_name in items:
                if old_name not in columns:
                    columns[old_name] = i

                else:
                    if old_name in column_names[columns[old_name]]:
                        old_md_index = columns[old_name]
                        modified_name = str(old_name) + \
                            f' [{mds[old_md_index].id}]'

                        if columns[old_name] not in df_changes:
                            df_changes[old_md_index] = {}

                        df_changes[old_md_index][old_name] = modified_name
                        column_names[old_md_index][modified_name] = \
                            column_names[old_md_index].pop(old_name)

                        column_sources[old_md_index].pop(old_name)
                        column_sources[old_md_index][modified_name] = \
                            mds[old_md_index]

                    if old_name == new_name:
                        if i not in df_changes:
                            df_changes[i] = {}
                        df_changes[i][old_name] = str(old_name) + \
                            f' [{mds[i].id}]'
                        new_name = str(old_name) + f' [{mds[i].id}]'
                        names[new_name] = names.pop(old_name)

                if old_name in sources:
                    sources.pop(old_name)
                sources[new_name] = md

        for change in df_changes:
            dfs[change] = dfs[change].rename(columns=df_changes[change])

        merged_df = dfs[0].join(dfs[1:], how='inner')
        # Not using DataFrame.empty because empty columns are allowed in
        # Metadata.
        if merged_df.index.empty:
            raise ValueError(
                "Cannot merge because there are no IDs shared across metadata "
                "objects.")

        column_names = {k: v for d in column_names for k, v in d.items()}
        column_sources = {k: v for d in column_sources for k, v in d.items()}
        merged_df.index.name = 'id'
        merged_md = self.__class__(merged_df)
        merged_md._init(column_names, column_sources)
        if df_changes:
            merged_md.contains_renamed_columns = True
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

        See Also
        --------
        get_ids
        filter_columns

        """
        filtered_df = self._filter_ids_helper(self._dataframe, self.get_ids(),
                                              ids_to_keep)
        filtered_md = self.__class__(filtered_df)
        filtered_md._init(self._column_names, self._column_sources)
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
            If ``True``, columns that contain a unique value for every ID will
            be dropped. Missing data (``np.nan``) are ignored when determining
            unique values. If a column consists solely of missing data, it will
            be dropped.
        drop_zero_variance : bool, optional
            If ``True``, columns that contain the same value for every ID will
            be dropped. Missing data (``np.nan``) are ignored when determining
            variance. If a column consists solely of missing data, it will be
            dropped.
        drop_all_missing : bool, optional
            If ``True``, columns that have a missing value (``np.nan``) for
            every ID will be dropped.

        Returns
        -------
        Metadata
            The metadata filtered by columns.

        See Also
        --------
        filter_ids

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
        filtered_md._init(self._column_names, self._column_sources)
        return filtered_md


class MetadataColumn(_MetadataBase, metaclass=abc.ABCMeta):
    """Abstract base class representing a single metadata column.

    Concrete subclasses represent specific metadata column types, e.g.
    ``CategoricalMetadataColumn`` and ``NumericMetadataColumn``.

    See the ``Metadata`` class docstring for details about ``Metadata`` and
    ``MetadataColumn`` objects, including a description of column types.

    The main difference in constructing ``MetadataColumn`` vs ``Metadata``
    objects is that ``MetadataColumn`` objects are constructed from a
    ``pandas.Series`` object instead of a ``pandas.DataFrame``. Otherwise, the
    same restrictions, considerations, and data normalization are applied as
    with ``Metadata`` objects.

    """
    # Abstract, must be defined by subclasses.
    type = None

    @classmethod
    @abc.abstractmethod
    def _is_supported_dtype(cls, dtype):
        """

        Contract: Return ``True`` if the series `dtype` is supported by this
        object and can be handled appropriately by ``_normalize_``. Return
        ``False`` otherwise.

        """
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
        """Metadata column name.

        This property is read-only.

        Returns
        -------
        str
            Metadata column name.

        """
        return self._series.name

    @property
    def artifacts(self):
        """Artifact that is the source of the column.

        This property is read-only.

        Returns
        -------
        tuple of qiime2.Artifact
            Source artifact of the column.

        """
        if self._source_metadata is None:
            return ()
        if self._source_metadata._source_artifact is not None:
            return (self._source_metadata._source_artifact,)

        name = self._series.name
        old_md = self._source_metadata
        new_md = self._source_metadata._column_sources[name]
        while True:
            if new_md is None:
                if old_md._source_artifact is not None:
                    return (old_md._source_artifact,)
                return ()
            if name in new_md._column_sources:
                old_md = new_md
                new_md = new_md._column_sources[name]
            else:
                name = old_md._column_names[name]
                old_md = new_md
                new_md = new_md._column_sources[name]

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
        self._source_metadata = None

    def _init(self, md):
        self._source_metadata = md

    def __repr__(self):
        """String summary of the metadata column."""
        return '<%s name=%r id_count=%d>' % (self.__class__.__name__,
                                             self.name, self.id_count)

    def __eq__(self, other):
        """Determine if this metadata column is equal to another.

        ``MetadataColumn`` objects are equal if their IDs, column names, column
        types, ID headers, source artifacts, and metadata values are equal.

        Parameters
        ----------
        other : MetadataColumn
            Metadata column to test for equality.

        Returns
        -------
        bool
            Indicates whether this ``MetadataColumn`` object is equal to
            `other`.

        See Also
        --------
        __ne__

        """
        return (
            super().__eq__(other) and
            self.name == other.name and
            self._series.equals(other._series) and
            # TODO: Consider how to handle this equality
            # My initial inclination was to make this test for equalty of
            # _source_metadata, but this seriously messes with some of the
            # rountrip tests in test_io, so it looks like this for now
            self.artifacts == other.artifacts
        )

    def __ne__(self, other):
        """Determine if this metadata column is not equal to another.

        ``MetadataColumn`` objects are not equal if their IDs, column names,
        column types, ID headers, source artifacts, or metadata values are not
        equal.

        Parameters
        ----------
        other : MetadataColumn
            Metadata column to test for inequality.

        Returns
        -------
        bool
            Indicates whether this ``MetadataColumn`` object is not equal to
            `other`.

        See Also
        --------
        __eq__

        """
        return not (self == other)

    def save(self, filepath):
        """Save a TSV metadata file containing this metadata column.

        The TSV metadata file format is described at https://docs.qiime2.org in
        the Metadata Tutorial.

        The file will always include the ``#q2:types`` directive in order to
        make the file roundtrippable without relying on column type inference.

        Parameters
        ----------
        filepath : str
            Path to save TSV metadata file at.

        """
        from .io import MetadataWriter
        MetadataWriter(self).write(filepath)

    def to_series(self):
        """Create a pandas series from the metadata column.

        The series index name (``Index.name``) will match this metadata
        column's ``id_header``, and the index will contain this metadata
        column's IDs. The series name will match this metadata column's name.

        Returns
        -------
        pandas.Series
            Series constructed from the metadata column.

        See Also
        --------
        to_dataframe

        """
        return self._series.copy()

    def to_dataframe(self):
        """Create a pandas dataframe from the metadata column.

        The dataframe will contain exactly one column. The dataframe's index
        name (``Index.name``) will match this metadata column's ``id_header``,
        and the index will contain this metadata column's IDs. The dataframe's
        column name will match this metadata column's name.

        Returns
        -------
        pandas.DataFrame
            Dataframe constructed from the metadata column.

        See Also
        --------
        to_series

        """
        return self._series.to_frame()

    def get_value(self, id):
        """Retrieve metadata column value associated with an ID.

        Parameters
        ----------
        id : str
            ID corresponding to the metadata column value to retrieve.

        Returns
        -------
        object
            Value associated with the provided `id`.

        """
        if id not in self._series.index:
            raise ValueError("ID %r is not present in %r" % (id, self))
        return self._series.loc[id]

    def has_missing_values(self):
        """Determine if the metadata column has one or more missing values.

        Returns
        -------
        bool
            ``True`` if the metadata column has one or more missing values
            (``np.nan``), ``False`` otherwise.

        See Also
        --------
        drop_missing_values
        get_ids

        """
        return len(self.get_ids(where_values_missing=True)) > 0

    def drop_missing_values(self):
        """Filter out missing values from the metadata column.

        Returns
        -------
        MetadataColumn
            Metadata column with missing values removed.

        See Also
        --------
        has_missing_values
        get_ids

        """
        missing = self.get_ids(where_values_missing=True)
        present = self.get_ids() - missing
        return self.filter_ids(present)

    def get_ids(self, where_values_missing=False):
        """Retrieve IDs matching search criteria.

        Parameters
        ----------
        where_values_missing : bool, optional
            If ``True``, only return IDs that are associated with missing
            values (``np.nan``). If ``False`` (the default), return all IDs in
            the metadata column.

        Returns
        -------
        set
            IDs matching search criteria.

        See Also
        --------
        ids
        filter_ids
        has_missing_values
        drop_missing_values

        """
        if where_values_missing:
            ids = self._series[self._series.isnull()].index
        else:
            ids = self._ids
        return set(ids)

    def filter_ids(self, ids_to_keep):
        """Filter metadata column by IDs.

        Parameters
        ----------
        ids_to_keep : iterable of str
            IDs that should be retained in the filtered ``MetadataColumn``
            object. If any IDs in `ids_to_keep` are not contained in this
            ``MetadataColumn`` object, a ``ValueError`` will be raised. The
            filtered ``MetadataColumn`` object will retain the same relative
            ordering of IDs in this ``MetadataColumn`` object. Thus, the
            ordering of IDs in `ids_to_keep` does not determine the ordering of
            IDs in the filtered ``MetadataColumn`` object.

        Returns
        -------
        MetadataColumn
            The metadata column filtered by IDs.

        See Also
        --------
        get_ids

        """
        filtered_series = self._filter_ids_helper(self._series, self.get_ids(),
                                                  ids_to_keep)
        filtered_mdc = self.__class__(filtered_series)
        filtered_mdc._init(self._source_metadata)
        return filtered_mdc


class CategoricalMetadataColumn(MetadataColumn):
    """A single metadata column containing categorical data.

    See the ``Metadata`` class docstring for details about ``Metadata`` and
    ``MetadataColumn`` objects, including a description of column types and
    supported data formats.

    """
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
    """A single metadata column containing numeric data.

    See the ``Metadata`` class docstring for details about ``Metadata`` and
    ``MetadataColumn`` objects, including a description of column types and
    supported data formats.

    """
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
