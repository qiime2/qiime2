# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2


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


# No input artifacts, only parameters.
def params_only_method(name: str, age: int) -> dict:
    return {name: age}


# No input artifacts or parameters.
def no_input_method() -> dict:
    return {'foo': 42}


def long_description_method(mapping1: dict, name: str, age: int) -> dict:
    return {name: age}


def identity_with_metadata(ints: list, metadata: qiime2.Metadata) -> list:
    return ints


def identity_with_metadata_category(ints: list,
                                    metadata: qiime2.MetadataCategory) -> list:
    return ints


def identity_with_optional_metadata(ints: list,
                                    metadata: qiime2.Metadata=None) -> list:
    return ints


def identity_with_optional_metadata_category(
        ints: list, metadata: qiime2.MetadataCategory=None) -> list:
    return ints


def optional_artifacts_method(ints: list, num1: int, optional1: list=None,
                              optional2: list=None, num2: int=None) -> list:
    result = ints + [num1]
    if optional1 is not None:
        result += optional1
    if optional2 is not None:
        result += optional2
    if num2 is not None:
        result += [num2]
    return result


def variadic_input_method(ints: list, int_set: int, nums: int,
                          opt_nums: int=None) -> list:
    results = []

    for int_list in ints:
        results += int_list
    results += sorted(int_set)
    results += nums
    if opt_nums:
        results += opt_nums

    return results
