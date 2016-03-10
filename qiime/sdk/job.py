# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


# TODO figure out what the public API should be. For now it's a simple struct.
class Job:
    def __init__(self, markdown, uuid, input_artifact_filepaths,
                 parameter_references, output_artifact_filepaths):
        self.markdown = markdown
        self._uuid = uuid
        self.input_artifact_filepaths = input_artifact_filepaths
        self.parameter_references = parameter_references
        self.output_artifact_filepaths = output_artifact_filepaths

    @property
    def uuid(self):
        return self._uuid
