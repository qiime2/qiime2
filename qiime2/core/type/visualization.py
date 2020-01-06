# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.core.type.template import TypeTemplate


class _Visualization(TypeTemplate):
    def get_kind(self):
        return "visualization"

    def __eq__(self, other):
        return type(self) is type(other)

    def get_field_names(self):
        return []

    def get_name(self):
        return "Visualization"

    def is_element(self, value):
        import qiime2.sdk

        return isinstance(value, qiime2.sdk.Visualization)

    def validate_field(self, name, field):
        raise TypeError

    def get_union_membership_expr(self, self_expr):
        return None

    def validate_predicate(self, predicate, expr):
        raise TypeError


Visualization = _Visualization()
