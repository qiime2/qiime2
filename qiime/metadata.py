# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


class Metadata:
    def __init__(self, dataframe):
        self._dataframe = dataframe

    def __repr__(self):
        return self._dataframe.to_json()

    def get_category(self, *names):
        if len(names) != 1:
            # TODO: Make this work with multiple columns as a single series
            raise NotImplementedError("Extracting multiple columns is not yet"
                                      " supported.")
        return MetadataCategory(self._dataframe[names[0]])

    def to_dataframe(self):
        return self._dataframe.copy()


class MetadataCategory:
    def __init__(self, series):
        self._series = series

    def __repr__(self):
        return self._series.to_json()

    def to_series(self):
        return self._series.copy()
