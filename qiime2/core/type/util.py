# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import collections

from qiime2.core.util import tuplize
from qiime2.core.type.primitive import Int, Float, Bool, Str
from qiime2.core.type.grammar import UnionExp, _ExpBase
from qiime2.core.type.parse import ast_to_type


def _booler(v):
    '''
    This is a psuedo-type --- we don't want Python's usual str->bool
    coercion business, just simple string matching (I think)
    '''
    # TODO: make case insensitive
    if v == 'True':
        return True
    elif v == 'False':
        return False
    else:
        raise ValueError('nuh uh uh, you didnt say the magic word')


_VARIADIC = {'List': list, 'Set': set}
# Order matters here:
_SEMANTIC_TO_PYTHON_TYPE_FUNC = {Int: int, Float: float, Bool: _booler, Str: str}
_SEMANTIC_TO_PYTHON_TYPE =      {Int: int, Float: float, Bool: bool,    Str: str}
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
    raise ValueError('Could not coerce value based on expression provided.')


# Value might be a boolean
def parse_primitive(t, value):
    expr = _norm_input(t)
    result = []
    collection_style = None
    allowed = tuple()
    homogeneous = True

    if is_metadata_type(expr):
        raise ValueError('what on earth were you thinking??')

    if expr in (Int, Float, Bool, Str):
        # No sense in walking over all options when we know what it should be
        allowed = tuplize(expr)
        value = tuplize(value)
    elif is_collection_type(expr):
        collection_style = interrogate_collection_type(expr)

        if collection_style.style == 'simple':
            allowed = _ordered_coercion(collection_style.members)
        elif collection_style.style == 'monomorphic':
            allowed = _ordered_coercion(tuple(collection_style.members))
        elif collection_style.style == 'composite':
            allowed = _ordered_coercion(tuple(collection_style.members))
            homogeneous = False
        elif collection_style.style == 'complex':
            pass
        else:
            raise ValueError('yikes, what are you doing here?')
    else:
        allowed = _ordered_coercion(tuple(_SEMANTIC_TO_PYTHON_TYPE.keys()))
        value = tuplize(value)

    for v in value:
        result.append(_interrogate_types(allowed, v))

    if not result:
        raise ValueError('bar')

    # Post-coerce validation --- some collection exprs require
    # homogeneous values
    if homogeneous:
        all_matching = False
        for member in allowed:
            if all(isinstance(x, _SEMANTIC_TO_PYTHON_TYPE[member])
                   for x in result):
                all_matching = True
        if not all_matching:
            allowed_ = tuple(map(lambda x: _SEMANTIC_TO_PYTHON_TYPE[x],
                                 allowed))
            # Should simple collections be included here, too?
            if collection_style \
                    and collection_style.style == 'monomorphic' \
                    and all(isinstance(x, allowed_) for x in result):
                # last allowed is the least common denominator for all vals
                # TODO: probably should _interrogate_types again...
                result = map(allowed_[-1], result)
            else:
                # TODO: is it even possible to hit this branch?
                raise ValueError('not all matching')

    if collection_style:
        return collection_style.view(result)
    else:
        return result[0]
    return result
