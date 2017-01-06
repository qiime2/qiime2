# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sqlite3

import pandas as pd


class Metadata:
    def __init__(self, dataframe):
        self._dataframe = dataframe

    def __repr__(self):
        return repr(self._dataframe)

    def _repr_html_(self):
        return self._dataframe._repr_html_()

    @classmethod
    def load(cls, path):
        try:
            df = pd.read_csv(path, sep='\t', dtype=object)
            df.set_index(df.columns[0], drop=True, append=False, inplace=True)
        except OSError:
            raise OSError(
                "Metadata file %s doesn't exist or isn't accessible (e.g., "
                "due to incompatible file permissions)." % path)
        except (pd.io.common.CParserError, KeyError):
            raise ValueError(
                'Metadata file format is invalid for file %s. Currently only '
                'QIIME 1 sample/feature metadata mapping files are officially '
                'supported. Sample metadata mapping files can be validated '
                'using Keemei: http://keemei.qiime.org' % path)
        return cls(df)

    def get_category(self, *names):
        if len(names) != 1:
            # TODO: Make this work with multiple columns as a single series
            raise NotImplementedError("Extracting multiple columns is not yet"
                                      " supported.")
        try:
            result = MetadataCategory(self._dataframe[names[0]])
        except KeyError:
            raise KeyError(
                '%s is not a category in metadata file. Available '
                'categories are %s.' %
                (names[0], ', '.join(self._dataframe.columns)))
        return result

    def to_dataframe(self):
        return self._dataframe.copy()

    def ids(self, where=None):
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

        """
        if where is None:
            return set(self._dataframe.index)

        conn = sqlite3.connect(':memory:')
        conn.row_factory = lambda cursor, row: row[0]

        self._dataframe.to_sql('metadata', conn)
        id_column = self._dataframe.index.name

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
                 'ORDER BY "{0}";'.format(id_column, where))

        try:
            c.execute(query)
        except sqlite3.OperationalError:
            conn.close()
            raise ValueError("Selection of IDs failed with query:\n %s"
                             % query)

        ids = set(c.fetchall())
        conn.close()
        return ids


class MetadataCategory:
    def __init__(self, series):
        self._series = series

    def __repr__(self):
        return repr(self._series)

    @classmethod
    def load(cls, path, category):
        return Metadata.load(path).get_category(category)

    def to_series(self):
        return self._series.copy()
