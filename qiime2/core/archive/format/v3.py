# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.core.archive.format.v2 as v2


class ArchiveFormat(v2.ArchiveFormat):
    # Exactly the same as v2, but inputs may be variadic where the UUIDs are in
    # a YAML sequence. Additionally `Set` is now represented as a sequence
    # with a custom !set tag.
    pass
