# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.core.archive.format.v5 as v5


class ArchiveFormat(v5.ArchiveFormat):
    # Exactly the same as v5, but we will have an execution_context describing
    # what executor was used and can have output collections
    pass
