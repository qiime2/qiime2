# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.core.archive.format.v1 as v1


class ArchiveFormat(v1.ArchiveFormat):
    # Exactly the same as v1, but in provenance, when the action type isn't
    # import, there is an `output-name` key in the action section with that
    # node's output name according to the action's signature object. Also has
    # pipeline action types.
    pass
