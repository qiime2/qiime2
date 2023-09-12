# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re

import qiime2.sdk
import qiime2.core.type as qtype
import qiime2.core.type.parse as _parse
from qiime2.core.type import (
    is_semantic_type, is_primitive_type, is_collection_type, is_metadata_type,
    is_visualization_type, interrogate_collection_type, parse_primitive,
    is_union, is_metadata_column_type)


__all__ = [
    'is_semantic_type', 'is_primitive_type', 'is_collection_type',
    'is_metadata_type', 'is_visualization_type', 'interrogate_collection_type',
    'type_from_ast', 'parse_primitive', 'parse_type', 'parse_format',
    'actions_by_input_type', 'is_union', 'is_metadata_column_type',
]


def camel_to_snake(name: str) -> str:
    """
    There are more comprehensive and faster ways of doing this (incl compiling)
    but it handles acronyms in semantic types nicely
    e.g. EMPSingleEndSequences -> emp_single_end_sequences
    c/o https://stackoverflow.com/a/1176023/9872253
    """
    # this will frequently be called on QIIME type expressions, so drop [ and ]
    name = re.sub(r'[\[\]]', '', name)
    # camel to snake
    name = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', name).lower()


def type_from_ast(ast, scope=None):
    """Convert a type ast (from `.to_ast()`) to a type expression.

    Parameters
    ----------
    ast : json compatible object
        The abstract syntax tree produced by `to_ast` on a type.
    scope : dict
        A dictionary to use between multiple calls to share scope between
        different types. This allows type variables from the same type map
        to be constructed from an equivalent type map. Scope should be shared
        within a given call signature, but not between call signatures.

    Returns
    -------
    type expression

    """
    return _parse.ast_to_type(ast, scope=scope)


def parse_type(string, expect=None):
    """Convert a string into a type expression

    Parameters
    ----------
    string : str
        The string type expression to convert into a TypeExpression
    expect : {'semantic', 'primitive', 'visualization'}, optional
        Will raise a TypeError if the resulting TypeExpression is not a member
        of `expect`.

    Returns
    -------
    type expression

    """
    if expect is not None and expect not in {'semantic', 'primitive',
                                             'visualization'}:
        raise ValueError("`expect` got %r, must be 'semantic', 'primitive',"
                         " 'visualization', or None." % (expect,))

    type_expr = _parse.ast_to_type(_parse.string_to_ast(string))

    if expect is None:
        pass
    elif expect == 'semantic' and qtype.is_semantic_type(type_expr):
        pass
    elif expect == 'primitive' and qtype.is_primitive_type(type_expr):
        pass
    elif expect == 'visualization' and type_expr == qtype.Visualization:
        pass
    else:
        raise TypeError("Type expression %r is not a %s type."
                        % (type_expr, expect))
    return type_expr


def parse_format(format_str):
    if format_str is None:
        return None

    pm = qiime2.sdk.PluginManager()
    try:
        format_record = pm.formats[format_str]
    except KeyError:
        raise TypeError("No format: %s" % format_str)
    return format_record.format


def actions_by_input_type(string):
    """Plugins and actions that have as input the artifact type (string)

    Parameters
    ----------
    string : str
        QIIME2 artifact type

    Returns
    -------
    list of tuples: [(q2.plugin, [q2.actions, ...]), ...]
    """
    commands = []
    if string is not None:
        query_type = qiime2.sdk.util.parse_type(string)

        pm = qiime2.sdk.PluginManager()
        for pgn, pg in pm.plugins.items():
            actions = list({a for an, a in pg.actions.items()
                            for iname, i in a.signature.inputs.items()
                            if i.qiime_type >= query_type})
            if actions:
                commands.append((pg, actions))

    return commands


def validate_result_collection_keys(*args):
    """Validate one or more strings intended for use as ResultCollection keys.

    This can be called on one or more keys provided as arguments:
    qiime2.sdk.util.validate_result_collection_keys('@', 'a1@')

    Or on a list, by unpacking it in the call:
    l = ['@', 'a1@']
    qiime2.sdk.util.validate_result_collection_keys(*l)
    """
    invalid_keys = []
    for key in args:
        if not isinstance(key, str) or bool(re.search(r'[^\w+-.]', key)):
            invalid_keys.append(key)

    if len(invalid_keys) > 0:
        raise KeyError('Invalid key(s) provided for ResultCollection. '
                       'ResultCollection keys must be strings and may only '
                       'contain the following characters: A-Z, a-z, 0-9, +, '
                       '-, ., and _. Offending keys include: '
                       f'{", ".join(map(str, invalid_keys))}')


def view_collection(collection, view_type):
    return {k: v.view(view_type) for k, v in collection.items()}
