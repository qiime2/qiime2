# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import qiime


# Artifacts and parameters.
def concatenate_ints(ints1: list, ints2: list, ints3: list, int1: int,
                     int2: int) -> list:
    return ints1 + ints2 + ints3 + [int1] + [int2]


# Multiple output artifacts.
def split_ints(ints: list) -> (list, list):
    middle = int(len(ints) / 2)
    left = ints[:middle]
    right = ints[middle:]
    return left, right


# No parameters, only artifacts.
def merge_mappings(mapping1: dict, mapping2: dict) -> dict:
    merged = mapping1.copy()
    for key, value in mapping2.items():
        if key in merged and merged[key] != value:
            raise ValueError(
                "Key %r exists in `mapping1` and `mapping2` with conflicting "
                "values: %r != %r" % (key, merged[key], value))
        merged[key] = value
    return merged


def identity_with_metadata(ints: list, metadata: qiime.Metadata) -> list:
    return ints


def identity_with_metadata_category(ints: list,
                                    metadata: qiime.MetadataCategory) -> list:
    return ints
