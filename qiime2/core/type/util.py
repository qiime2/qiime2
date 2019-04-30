# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import collections

from qiime2.core.util import tuplize
from qiime2.core.type.collection import List, Set
from qiime2.core.type.primitive import Int, Float, Bool, Str
from qiime2.core.type.grammar import UnionExp, _ExpBase
from qiime2.core.type.parse import ast_to_type


def _booler(v):
    '''
    This is a psuedo-type --- we don't want Python's usual str->bool
    coercion business, just simple string matching (I think)
    '''
    if isinstance(v, bool):
        return v
    # TODO: make case insensitive
    elif v == 'True':
        return True
    elif v == 'False':
        return False
    else:
        raise ValueError('nuh uh uh, you didnt say the magic word')


def _inter(v):
    if isinstance(v, bool):
        raise ValueError('uh no')
    elif isinstance(v, int):
        return v
    else:
        return int(v)


_VARIADIC = {'List': list, 'Set': set}
_VARIADIC2 = {list: List, set: Set}
# Order matters here:
_SEMANTIC_TO_PYTHON_TYPE_FUNC = {Int: _inter, Float: float, Bool: _booler, Str: str}
_SEMANTIC_TO_PYTHON_TYPE = {Int: int, Float: float, Bool: bool,    Str: str}
_COERCE_ERROR = ValueError(
    'Could not coerce value based on expression provided.')
CollectionStyle = collections.namedtuple(
    'CollectionStyle', ['style', 'members', 'view', 'expr'])


def _norm_input(t):
    if type(t) is dict:
        return ast_to_type(t)
    elif not isinstance(t, _ExpBase):
        raise TypeError("%r is not a QIIME 2 type" % (t,))
    return t


def is_qiime_type(t):
    try:
        _norm_input(t)
    except Exception:
        return False
    else:
        return True


def is_primitive_type(t):
    expr = _norm_input(t)
    return hasattr(expr, 'kind') and expr.kind == 'primitive'


def is_metadata_type(t):
    expr = _norm_input(t)
    return is_primitive_type(t) and expr.name.startswith('Metadata')


def is_semantic_type(t):
    expr = _norm_input(t)
    return hasattr(expr, 'kind') and expr.kind == 'semantic-type'


def is_visualization_type(t):
    expr = _norm_input(t)
    return hasattr(expr, 'kind') and expr.kind == 'visualization'


def is_collection_type(t):
    expr = _norm_input(t)

    if expr.name in _VARIADIC:
        return True

    if isinstance(expr, UnionExp):
        for m in expr.members:
            if m.name in _VARIADIC:
                return True

    return False


def interrogate_collection_type(t):
    expr = _norm_input(t)
    style = None    # simple, monomorphic, composite, complex
    members = None  # T     , [T1, T2]   , [T1, T2],  [[T1], [T2, T3]]
    view = None  # set, list

    if expr.name in _VARIADIC:
        view = _VARIADIC[expr.name]
        field, = expr.fields
        if isinstance(field, UnionExp):
            style = 'composite'
            members = list(field.members)
        else:
            style = 'simple'
            members = field
    elif isinstance(expr, UnionExp):
        if expr.members[0].name in _VARIADIC:
            members = []
            for member in expr.members:
                field, = member.fields
                if isinstance(field, UnionExp):
                    style = 'complex'
                    members.append(list(field.members))
                else:
                    members.append([field])
                    if style != 'complex':
                        style = 'monomorphic'

            # use last iteration
            view = _VARIADIC[member.name]
            if style == 'monomorphic':
                members = [m[0] for m in members]

    return CollectionStyle(style=style, members=members, view=view, expr=expr)


def _ordered_coercion(types):
    types = tuplize(types)
    return tuple(k for k in _SEMANTIC_TO_PYTHON_TYPE.keys() if k in types)


def _interrogate_types(allowed, value):
    for coerce_type in (_SEMANTIC_TO_PYTHON_TYPE_FUNC[x] for x in allowed):
        try:
            return coerce_type(value)
        except ValueError:
            pass
    raise _COERCE_ERROR


# Value might be a boolean
def parse_primitive(t, value):
    expr = _norm_input(t)
    collection_style = interrogate_collection_type(expr)
    result = []
    allowed = tuple()
    homogeneous = True

    if is_metadata_type(expr):
        raise ValueError('what on earth were you thinking??')

    if collection_style.style == 'simple':
        allowed = _ordered_coercion(collection_style.members)
    elif collection_style.style == 'monomorphic':
        allowed = _ordered_coercion(tuple(collection_style.members))
    elif collection_style.style == 'composite':
        allowed = _ordered_coercion(tuple(collection_style.members))
        homogeneous = False
    elif collection_style.style == 'complex':
        # Sort by smallest group of members first, go with the simplest
        # explanation of values first
        for subexpr in sorted(collection_style.members, key=len):
            expr = _VARIADIC2[collection_style.view][UnionExp(subexpr)]
            try:
                return parse_primitive(expr, value)
            except ValueError:
                # TODO: is it possible to hit this branch?
                pass
        # TODO: is it possible to hit this branch?
        raise _COERCE_ERROR
    elif collection_style.style is None:
        value = tuplize(value)
        if expr in (Int, Float, Bool, Str):
            # No sense in walking over all options when we know
            # what it should be
            allowed = tuplize(expr)
        else:
            allowed = _ordered_coercion(tuple(_SEMANTIC_TO_PYTHON_TYPE.keys()))

    assert allowed

    for v in value:
        result.append(_interrogate_types(allowed, v))

    if not result:
        raise ValueError('bar')

    # Some exprs require homogeneous values, make it so
    if homogeneous:
        all_matching = False
        for member in allowed:
            if all(type(x) == _SEMANTIC_TO_PYTHON_TYPE[member]
                   for x in result):
                all_matching = True
                break
        if not all_matching:
            # Should simple collections be included here, too?
            if collection_style and collection_style.style == 'monomorphic':
                for subexpr in allowed:
                    expr = _VARIADIC2[collection_style.view][subexpr]
                    try:
                        return parse_primitive(expr, value)
                    except ValueError:
                        # TODO: is it possible to hit this branch?
                        pass
                raise _COERCE_ERROR
            else:
                # TODO: is it possible to hit this branch?
                raise ValueError('not all matching')

    if collection_style.view is None:
        return result[0]
    else:
        return collection_style.view(result)
