# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

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
            if type_ is qiime.core.type.Visualization:
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
        type_ : qiime.plugin.Type or str
            Concrete semantic type of the artifact.
        provenance : qiime.sdk.Provenance
            Artifact provenance.

        """
        artifact = cls.__new__(cls)
        artifact._archiver = qiime.core.archiver.Archiver(uuid.uuid4(), type_,
                                                          provenance)
        artifact._assert_valid_type(artifact.type)
        # TODO do some validation, such as making sure `view` can be
        # written to `type_` data layout.
        artifact._save_view(view)
        return artifact

    def _save_view(self, view):
        data_layout = self._get_data_layout()
        view_type = type(view)
        writer = data_layout.writers[view_type]
        self._archiver.save_data(view, writer)

    def _get_data_layout(self):
        pm = qiime.sdk.PluginManager()

        data_layout = None
        for semantic_type, datalayout in \
                pm.semantic_type_to_data_layouts.items():
            if self.type <= semantic_type:
                data_layout = datalayout
                break

        if data_layout is None:
            raise TypeError(
                "Artifact semantic type %r does not have a compatible data "
                "layout." % self.type)

        return data_layout

    def _orphan(self, pid):
        self._archiver.orphan(pid)

    def view(self, view_type):
        data_layout = self._get_data_layout()
        reader = data_layout.readers[view_type]
        return self._archiver.load_data(reader)
