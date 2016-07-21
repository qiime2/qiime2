# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

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
            df = pd.read_csv(path, sep='\t', dtype=object, index_col=0)
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
