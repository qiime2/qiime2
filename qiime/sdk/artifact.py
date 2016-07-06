# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import functools
import uuid

import qiime.core.archiver
import qiime.core.result_base
import qiime.core.type
import qiime.sdk


class Artifact(qiime.core.result_base.ResultBase):
    @classmethod
    def _assert_valid_type(cls, type_):
        if not (qiime.core.type.is_semantic_type(type_) and
                type_.is_concrete()):
            error_suffix = ""
            if type_ == qiime.core.type.Visualization:
                error_suffix = " Use qiime.sdk.Visualization with this type."
            raise TypeError(
                "An artifact requires a concrete semantic type, not type %r.%s"
                % (type_, error_suffix))

    @classmethod
    def _from_view(cls, view, type_, provenance):
        """

        Parameters
        ----------
        view : Python object
            View to serialize.
        type_ : qiime.plugin.Type
            Concrete semantic type of the artifact.
        provenance : qiime.sdk.Provenance
            Artifact provenance.

        """
        cls._assert_valid_type(type_)
        data_layout = cls._get_data_layout(type_)
        view_type = type(view)
        # TODO better error handling for when `view` cannot be written to
        # `type_` data layout.
        writer = data_layout.writers[view_type]
        writer = functools.partial(writer, view)

        artifact = cls.__new__(cls)
        artifact._archiver = qiime.core.archiver.Archiver(
            uuid.uuid4(), type_, provenance, data_initializer=writer)
        return artifact

    @classmethod
    def _get_data_layout(cls, type_):
        pm = qiime.sdk.PluginManager()

        data_layout = None
        for semantic_type, datalayout in \
                pm.semantic_type_to_data_layouts.items():
            if type_ <= semantic_type:
                data_layout = datalayout
                break

        if data_layout is None:
            raise TypeError(
                "Artifact semantic type %r does not have a compatible data "
                "layout." % type_)

        return data_layout

    def view(self, view_type):
        data_layout = self._get_data_layout(self.type)
        reader = data_layout.readers[view_type]
        return self._archiver.load_data(reader)
