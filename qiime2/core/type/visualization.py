# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.core.type.template import TypeTemplate, instantiate

@instantiate
class Visualization(TypeTemplate):
    def get_kind(self):
        return "visualization"
    def __eq__(self, other):
        return type(self) is type(other)

    def get_field_names(self):
        return []

    def get_name(self):
        return self.__class__.__name__

    def is_element(self, value, expr):
        import qiime2.sdk
        return isinstance(value, qiime2.sdk.Visualization)

    def validate_fields(self, fields, expr):
        raise TypeError

    def validate_intersection(self, other):
        pass

    def collapse_intersection(self, other):
        return self

    def validate_union(self, other):
        raise TypeError

    def validate_predicate(self, predicate, expr):
        raise TypeError
