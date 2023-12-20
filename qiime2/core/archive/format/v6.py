# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.core.archive.format.v5 as v5


class ArchiveFormat(v5.ArchiveFormat):
    # - Adds execution_context to the execution section of action.yaml.
    #   This looks like:
    #
    #   execution_context:
    #       type: parsl/synchronous/asynchronous
    #       parsl_type (if type is parsl): Type of executor
    #
    #   NOTE: Import actions will not have an execution_context section
    #
    # - Adds support for output collections.
    #   This looks like:
    #
    #   output-name:
    #   - The name of the entire output collection (the qiime output name)
    #   - The key of this element in the collection
    #   - The index of this element in the collection ex. '5/10' for the 5th
    #     out of 10 elements in the collection
    #
    #   Input collections now look like:
    #
    #   input-name:
    #   - key: value
    #   - key: value
    #   etc. for n elements
    pass
