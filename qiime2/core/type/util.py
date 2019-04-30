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


def val_to_bool(value):
    if type(value) == bool:
        return value
    elif str(value).lower() == 'true':
        return True
    elif str(value).lower() == 'false':
        return False
    else:
        raise ValueError('Could not cast to bool')


def val_to_int(v):
    if type(v) == bool:
        raise ValueError('Could not cast to int')
    elif type(v) == int:
        return v
    else:
        return int(v)


VariadicRecord = collections.namedtuple('VariadicRecord', ['pytype', 'q2type'])
_VARIADIC = {
    'List': VariadicRecord(pytype=list, q2type=List),
    'Set': VariadicRecord(pytype=set, q2type=Set),
}

CoercionRecord = collections.namedtuple('CoercionRecord', ['func', 'pytype'])
# Beware visitor, order matters in this here mapper
_COERCION_MAPPER = {
    Int: CoercionRecord(pytype=int, func=val_to_int),
    Float: CoercionRecord(pytype=float, func=float),
    Bool: CoercionRecord(pytype=bool, func=val_to_bool),
    Str: CoercionRecord(pytype=str, func=str),
}
_COERCE_ERROR = ValueError(
    'Could not coerce value based on expression provided.')

CollectionStyle = collections.namedtuple(
    'CollectionStyle', ['style', 'members', 'view', 'expr', 'base'])


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
    base = None

    if expr.name in _VARIADIC:
        view, base = _VARIADIC[expr.name]
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
            view, base = _VARIADIC[member.name]
            if style == 'monomorphic':
                members = [m[0] for m in members]

    return CollectionStyle(style=style, members=members, view=view,
                           expr=expr, base=base)


def _ordered_coercion(types):
    types = tuple(types)
    return tuple(k for k in _COERCION_MAPPER.keys() if k in types)


def _interrogate_types(allowed, value):
    ordered_allowed = _ordered_coercion(allowed)
    for coerce_type in (_COERCION_MAPPER[x].func for x in ordered_allowed):
        try:
            return coerce_type(value)
        except ValueError:
            pass
    raise _COERCE_ERROR


def parse_primitive(t, value):
    expr = _norm_input(t)
    collection_style = interrogate_collection_type(expr)
    result = []
    allowed = None
    homogeneous = True

    if is_metadata_type(expr):
        raise ValueError('%r may not be parsed with this util.' % (expr,))

    if collection_style.style in ('simple', 'monomorphic', 'composite'):
        allowed = collection_style.members

    if collection_style.style == 'composite':
        homogeneous = False
    elif collection_style.style == 'complex':
        # Sort here so that we can start with any simple lists in the memberset
        for subexpr in sorted(collection_style.members, key=len):
            expr = collection_style.base[UnionExp(subexpr)]
            try:
                return parse_primitive(expr, value)
            except ValueError:
                pass
        raise _COERCE_ERROR
    elif collection_style.style is None:
        value = tuplize(value)
        if expr in (Int, Float, Bool, Str):
            # No sense in walking over all options when we know
            # what it should be
            allowed = expr
        else:
            allowed = _COERCION_MAPPER.keys()
    else:
        pass

    assert allowed is not None

    for v in value:
        result.append(_interrogate_types(allowed, v))

    # Some exprs require homogeneous values, make it so
    if homogeneous:
        all_matching = False
        for member in allowed:
            if all(type(x) == _COERCION_MAPPER[member].pytype
                   for x in result):
                all_matching = True
                break
        if not all_matching and collection_style.style == 'monomorphic':
            for subexpr in allowed:
                expr = collection_style.base[subexpr]
                try:
                    return parse_primitive(expr, value)
                except ValueError:
                    pass
            raise _COERCE_ERROR

    if collection_style.view is None:
        return result[0]
    else:
        return collection_style.view(result)
