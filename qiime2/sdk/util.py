# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

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


# NOTE: Can put all artifact API methods on this class because we know what
# type it will be (maybe not uuid)
class ProxyArtifact:
    """This represents a future artifact that is being returned by a Parsl app
    """
    def __init__(self, future, selector, signature=None):
        """We have a future that represents the results of some QIIME 2 action,
        and we have a selector indicating specifically which result we want
        """
        self.future = future
        self.selector = selector
        self.signature = signature

    def __repr__(self):
        if self.signature is None:
            return f'Unknown Type: {object.__repr__(self)}'
        else:
            return f'<artifact: {self.signature[self.selector].qiime_type}>'

    def get_element(self, results):
        """Get the result we want off of the future we have
        """
        return getattr(results, self.selector)

    def view(self, type):
        """If we want to view the result we need the future to be resolved
        """
        return self.get_element(self.future.result()).view(type)


class ProxyResults:
    """This represents future results that are being returned by a Parsl app
    """
    def __init__(self, future, signature):
        """We have the future results and the outputs portion of the signature
        of the action creating the results
        """
        self.future = future
        self.signature = signature

    def __iter__(self):
        """Give us a ProxyArtifact for each result in the future
        """
        for s in self.signature:
            yield ProxyArtifact(self.future, s, self.signature)

    def __getattr__(self, attr):
        """Get a particular ProxyArtifact out of the future
        """
        return ProxyArtifact(self.future, attr, self.signature)

    def __getitem__(self, index):
        return ProxyArtifact(
            self.future, list(self.signature.keys())[index], None)

    def result(self):
        return self.future.result()
