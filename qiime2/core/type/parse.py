# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import ast

from . import grammar, meta, collection, primitive, semantic, visualization


def string_to_ast(type_expr):
    try:
        parsed = ast.parse(type_expr)
    except SyntaxError:
        raise ValueError("%r could not be parsed, it may not be a QIIME 2 type"
                         " or it may not be an atomic type. Use"
                         " `ast_to_type` instead." % (type_expr,))

    if type(parsed) is not ast.Module:
        # I don't think this branch *can* be hit
        raise ValueError("%r is not a type expression." % (type_expr,))

    try:
        expr, = parsed.body
    except ValueError:
        raise ValueError("Only one type expression may be parse at a time, got"
                         ": %r" % (type_expr,))

    return _expr(expr.value)


def _expr(expr):
    node = type(expr)

    if node is ast.Name:
        return _build_atomic(expr.id)

    if node is ast.Call:
        args = _parse_args(expr.args)
        kwargs = _parse_kwargs(expr.keywords)
        return _build_predicate(expr.func.id, args, kwargs)

    if node is ast.Subscript:
        field_expr = expr.slice.value

        if type(field_expr) is ast.Tuple:
            field_expr = field_expr.elts
        else:
            field_expr = (field_expr,)

        base = _expr(expr.value)
        base['fields'] = [_expr(e) for e in field_expr]
        return base

    if node is ast.BinOp:
        op = type(expr.op)
        left = _expr(expr.left)
        right = _expr(expr.right)

        if op is ast.Mod:
            left['predicate'] = right
            return left
        if op is ast.BitOr:
            return _build_union(left, right)
        if op is ast.BitAnd:
            return _build_intersection(left, right)

    raise ValueError("Unknown expression: %r" % node)


def _convert_literals(expr):
    node = type(expr)

    if node is ast.List:
        return [_convert_literals(e) for e in expr.elts]

    if node is ast.Set:
        return {_convert_literals(e) for e in expr.elts}

    if node is ast.Tuple:
        return tuple(_convert_literals(e) for e in expr.elts)

    if node is ast.Dict:
        return {_convert_literals(k): _convert_literals(v)
                for k, v in zip(expr.keys, expr.values)}

    if node is ast.NameConstant:
        return expr.value

    if node is ast.Name and expr.id == 'inf':
        return float('inf')

    if node is ast.Num:
        return expr.n

    if node is ast.Str:
        return expr.s

    raise ValueError("Unknown literal: %r" % node)


def _parse_args(args):
    return tuple(_convert_literals(e) for e in args)


def _parse_kwargs(kwargs):
    return {e.arg: _convert_literals(e.value) for e in kwargs}


def _build_predicate(name, args, kwargs):
    base = {
        'type': 'predicate',
        'name': name
    }

    if name == 'Properties':
        return _build_properties(base, args, kwargs)
    if name == 'Range':
        return _build_range(base, args, kwargs)
    if name == 'Choices':
        return _build_choices(base, args, kwargs)


def _normalize_input_collection(args):
    if len(args) == 1 and isinstance(args[0], (list, set, tuple)):
        return tuple(args[0])
    return args


def _build_choices(base, args, kwargs):
    if 'choices' in kwargs:
        args = (kwargs['choices'],)
    args = _normalize_input_collection(args)
    base['choices'] = list(args)
    return base


def _build_range(base, args, kwargs):
    inclusive_start = kwargs.get('inclusive_start', True)
    inclusive_end = kwargs.get('inclusive_end', False)
    start = None
    end = None

    if len(args) == 1:
        end = args[0]
    elif len(args) != 0:
        start, end = args

    if start == float('-inf'):
        start = None
    if end == float('inf'):
        end = None

    base['range'] = [start, end]
    base['inclusive'] = [inclusive_start, inclusive_end]
    return base


def _build_properties(base, args, kwargs):
    exclude = kwargs.get('exclude', [])
    if 'include' in kwargs:
        args = (kwargs['include'],)

    args = _normalize_input_collection(args)
    base['include'] = list(args)
    base['exclude'] = list(exclude)
    return base


def _build_atomic(name):
    return {
        'type': 'expression',
        'builtin': name in {'Str', 'Int', 'Float', 'Bool',
                            'List', 'Set', 'Tuple', 'Visualization',
                            'Metadata', 'MetadataColumn', 'Numeric',
                            'Categorical'},
        'name': name,
        'predicate': None,
        'fields': []
    }


def _build_union(left, right):
    return _build_ident(left, right, 'union')


def _build_intersection(left, right):
    return _build_ident(left, right, 'intersection')


def _build_ident(left, right, type):
    members = []
    if left['type'] == type:
        members.extend(left['members'])
    else:
        members.append(left)

    if right['type'] == type:
        members.extend(right['members'])
    else:
        members.append(right)

    return {
        'type': type,
        'members': members
    }


def ast_to_type(json_ast, scope=None):
    if scope is None:
        scope = {}
    type_ = json_ast['type']

    if type_ == 'expression':
        predicate = json_ast['predicate']
        if predicate is not None:
            predicate = ast_to_type(predicate, scope=scope)

        fields = json_ast['fields']
        if len(fields) > 0:
            fields = [ast_to_type(f, scope=scope) for f in fields]

        name = json_ast['name']
        if not json_ast['builtin']:
            base_template = semantic.SemanticType(name).template
        elif name == 'Visualization':
            return visualization.Visualization
        elif name in {'List', 'Set', 'Tuple'}:
            base_template = getattr(collection, name).template
        else:
            base_template = getattr(primitive, name).template

        return grammar.TypeExp(base_template,
                               fields=fields, predicate=predicate)

    if type_ == 'predicate':
        name = json_ast['name']
        if name == 'Choices':
            return primitive.Choices(json_ast['choices'])
        if name == 'Range':
            return primitive.Range(*json_ast['range'],
                                   inclusive_start=json_ast['inclusive'][0],
                                   inclusive_end=json_ast['inclusive'][1])
        if name == 'Properties':
            return semantic.Properties(json_ast['include'],
                                       exclude=json_ast['exclude'])

    if type_ == 'union':
        members = [ast_to_type(m, scope=scope) for m in json_ast['members']]
        return grammar.UnionExp(members)

    if type_ == 'intersection':
        members = [ast_to_type(m, scope=scope) for m in json_ast['members']]
        return grammar.IntersectionExp(members)

    if type_ == 'variable':
        var_group = json_ast['group']
        if var_group not in scope:
            mapping = {}
            out_idx = json_ast['outputs']
            for entry in json_ast['mapping']:
                entry = [ast_to_type(e) for e in entry]
                mapping[tuple(entry[:out_idx])] = tuple(entry[out_idx:])
            scope[var_group] = list(meta.TypeMap(mapping))

        return scope[var_group][json_ast['index']]
