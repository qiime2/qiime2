# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
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


# Unioned primitives
def unioned_primitives(foo: int, bar: str = 'auto_bar') -> dict:
    return {'foo': foo, 'bar': bar}


# No input artifacts or parameters.
def no_input_method() -> dict:
    return {'foo': 42}


def deprecated_method() -> dict:
    return {'foo': 43}


def long_description_method(mapping1: dict, name: str, age: int) -> dict:
    return {name: age}


def docstring_order_method(req_input: dict, req_param: str,
                           opt_input: dict = None,
                           opt_param: int = None) -> dict:
    return {req_param: opt_param}


def identity_with_metadata(ints: list, metadata: qiime2.Metadata) -> list:
    assert isinstance(metadata, qiime2.Metadata)
    return ints


# TODO unit tests (test_method.py) for 3 variations of MetadataColumn methods
# below
def identity_with_metadata_column(ints: list,
                                  metadata: qiime2.MetadataColumn) -> list:
    assert isinstance(metadata, (qiime2.CategoricalMetadataColumn,
                                 qiime2.NumericMetadataColumn))
    return ints


def identity_with_categorical_metadata_column(
        ints: list, metadata: qiime2.CategoricalMetadataColumn) -> list:
    assert isinstance(metadata, qiime2.CategoricalMetadataColumn)
    return ints


def identity_with_numeric_metadata_column(
        ints: list, metadata: qiime2.NumericMetadataColumn) -> list:
    assert isinstance(metadata, qiime2.NumericMetadataColumn)
    return ints


def identity_with_optional_metadata(ints: list,
                                    metadata: qiime2.Metadata = None) -> list:
    assert isinstance(metadata, (qiime2.Metadata, type(None)))
    return ints


def identity_with_optional_metadata_column(
        ints: list, metadata: qiime2.MetadataColumn = None) -> list:
    assert isinstance(metadata, (qiime2.CategoricalMetadataColumn,
                                 qiime2.NumericMetadataColumn,
                                 type(None)))
    return ints


def optional_artifacts_method(ints: list, num1: int, optional1: list = None,
                              optional2: list = None,
                              num2: int = None) -> list:
    result = ints + [num1]
    if optional1 is not None:
        result += optional1
    if optional2 is not None:
        result += optional2
    if num2 is not None:
        result += [num2]
    return result


def variadic_input_method(ints: list, int_set: int, nums: int,
                          opt_nums: int = None) -> list:
    results = []

    for int_list in ints:
        results += int_list
    results += sorted(int_set)
    results += nums
    if opt_nums:
        results += opt_nums

    return results
