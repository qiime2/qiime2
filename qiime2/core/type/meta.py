# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import itertools
from types import MappingProxyType

from ..util import superscript, tuplize, ImmutableBase
from .grammar import UnionExp, TypeExp
from .collection import Tuple


class TypeVarExp(UnionExp):
    def __init__(self, members, tmap, input=False, output=False, index=None):
        self.mapping = tmap
        self.input = input
        self.output = output
        self.index = index

        super().__init__(members)

    def __repr__(self):
        numbers = {}
        for idx, m in enumerate(self.members, 1):
            if m in numbers:
                numbers[m] += superscript(',' + str(idx))
            else:
                numbers[m] = superscript(idx)
        return " | ".join([repr(k) + v for k, v in numbers.items()])

    def uniq_upto_sub(self, a_expr, b_expr):
        """
        Two elements are unique up to a subtype if they are indistinguishable
        with respect to that subtype. In the case of a type var, that means
        the same branches must be "available" in the type map.

        This means that A or B may have additional refinements (or may even be
        subtypes of each other), so long as that does not change the branch
        chosen by the type map.

        """
        a_branches = [m for m in self.members if a_expr <= m]
        b_branches = [m for m in self.members if b_expr <= m]
        return a_branches == b_branches

    def __eq__(self, other):
        return (type(self) is type(other)
                and self.index == other.index
                and self.mapping == other.mapping)

    def __hash__(self):
        return hash(self.index) ^ hash(self.mapping)

    def is_concrete(self):
        return False

    def can_intersect(self):
        return False

    def get_union_membership_expr(self, self_expr):
        return None

    def _is_subtype_(self, other):
        return all(m <= other for m in self.members)

    def _is_supertype_(self, other):
        return any(m >= other for m in self.members)

    def __iter__(self):
        yield from self.members

    def unpack_union(self):
        yield self

    def to_ast(self):
        return {
            "type": "variable",
            "index": self.index,
            "group": id(self.mapping),
            "outputs": self.mapping.input_width(),
            "mapping": [
                ([k.to_ast() for k in key.fields]
                 + [v.to_ast() for v in value.fields])
                for key, value in self.mapping.lifted.items()]
        }


class TypeMap(ImmutableBase):
    def __init__(self, mapping):
        mapping = {Tuple[tuplize(k)]: Tuple[tuplize(v)]
                   for k, v in mapping.items()}
        branches = list(mapping)
        for i, a in enumerate(branches):
            for j in range(i, len(branches)):
                b = branches[j]
                try:
                    intersection = a & b
                except TypeError:
                    raise ValueError("Cannot place %r and %r in the same "
                                     "type variable." % (a, b))
                if (intersection.is_bottom()
                        or intersection is a or intersection is b):
                    continue

                for k in range(i):
                    if intersection <= branches[k]:
                        break
                else:
                    raise ValueError(
                        "Ambiguous resolution for invocations with type %r."
                        " Could match %r or %r, add a new branch ABOVE these"
                        " two (or modify these branches) to correct this."
                        % (intersection.fields, a.fields, b.fields))
        self.__lifted = mapping
        super()._freeze_()

    @property
    def lifted(self):
        return MappingProxyType(self.__lifted)

    def __eq__(self, other):
        return self is other

    def __hash__(self):
        return hash(id(self))

    def __iter__(self):
        for idx, members in enumerate(
                zip(*(k.fields for k in self.lifted.keys()))):
            yield TypeVarExp(members, self, input=True, index=idx)

        yield from self.iter_outputs()

    def solve(self, *inputs):
        inputs = Tuple[inputs]
        for branch, outputs in self.lifted.items():
            if inputs <= branch:
                return outputs.fields

    def input_width(self):
        return len(next(iter(self.lifted.keys())).fields)

    def iter_outputs(self, *, _double_as_input=False):
        start = self.input_width()
        for idx, members in enumerate(
                zip(*(v.fields for v in self.lifted.values())), start):
            yield TypeVarExp(members, self, output=True, index=idx,
                             input=_double_as_input)


def _get_intersections(listing):
    intersections = []
    for a, b in itertools.combinations(listing, 2):
        i = a & b
        if i.is_bottom() or i is a or i is b:
            continue
        intersections.append(i)
    return intersections


def TypeMatch(listing):
    listing = list(listing)
    intersections = _get_intersections(listing)
    to_add = []
    while intersections:
        to_add.extend(intersections)
        intersections = _get_intersections(intersections)
    mapping = TypeMap({l: l for l in list(reversed(to_add)) + listing})
    # TypeMatch only produces a single variable
    # iter_outputs is used by match for solving, so the index must match
    return next(iter(mapping.iter_outputs(_double_as_input=True)))


def select_variables(expr):
    """When called on an expression, will yield selectors to the variable.

    A selector will either return the variable (or equivalent fragment) in
    an expression, or will return an entirely new expression with the
    fragment replaced with the value of `swap`.

    e.g.
    >>> select_u, select_t = select_variables(Example[T] % U)
    >>> t = select_t(Example[T] % U)
    >>> assert T is t
    >>> u = select_u(Example[T] % U)
    >>> assert U is u

    >>> frag = select_t(Example[Foo] % Bar)
    >>> assert frag is Foo
    >>> new_expr = select_t(Example[T] % U, swap=frag)
    >>> assert new_expr == Example[Foo] % U

    """
    if type(expr) is TypeVarExp:
        def select(x, swap=None):
            if swap is not None:
                return swap
            return x

        yield select
        return

    if type(expr) is not TypeExp:
        return

    if type(expr.full_predicate) is TypeVarExp:
        def select(x, swap=None):
            if swap is not None:
                return x.duplicate(predicate=swap)
            return x.full_predicate

        yield select

    for idx, field in enumerate(expr.fields):
        for sel in select_variables(field):
            # Without this closure, the idx in select will be the last
            # value of the enumerate, same for sel
            # (Same problem as JS with callbacks inside a loop)
            def closure(idx, sel):
                def select(x, swap=None):
                    if swap is not None:
                        new_fields = list(x.fields)
                        new_fields[idx] = sel(x.fields[idx], swap)
                        return x.duplicate(fields=tuple(new_fields))
                    return sel(x.fields[idx])
                return select
            yield closure(idx, sel)


def match(provided, inputs, outputs):
    provided_binding = {}
    error_map = {}
    for key, expr in inputs.items():
        for selector in select_variables(expr):
            var = selector(expr)
            provided_fragment = selector(provided[key])
            try:
                current_binding = provided_binding[var]
            except KeyError:
                provided_binding[var] = provided_fragment
                error_map[var] = provided[key]
            else:
                if not var.uniq_upto_sub(current_binding, provided_fragment):
                    raise ValueError("Received %r and %r, but expected %r"
                                     " and %r to match (or to select the same"
                                     " output)."
                                     % (error_map[var], provided[key],
                                        current_binding, provided_fragment))

    # provided_binding now maps TypeVarExp instances to a TypeExp instance
    # which is the relevent fragment from the provided input types

    grouped_maps = {}
    for item in provided_binding.items():
        var = item[0]
        if var.mapping not in grouped_maps:
            grouped_maps[var.mapping] = [item]
        else:
            grouped_maps[var.mapping].append(item)

    # grouped_maps now maps a TypeMap instance to tuples of
    # (TypeVarExp, TypeExp) which are the items of provided_binding
    # i.e. all of the bindings are now grouped under their shared type maps

    output_fragments = {}
    for mapping, group in grouped_maps.items():
        if len(group) != mapping.input_width():
            raise ValueError("Missing input variables")

        inputs = [x[1] for x in sorted(group, key=lambda x: x[0].index)]
        solved = mapping.solve(*inputs)
        if solved is None:
            provided = tuple(error_map[x[0]]
                             for x in sorted(group, key=lambda x: x[0].index))
            raise ValueError("No solution for inputs: %r, check the signature "
                             "to see valid combinations." % (provided,))

        # type vars share identity by instance of map and index, so we will
        # be able to see the "same" vars again when looking up the outputs
        for var, out in zip(mapping.iter_outputs(), solved):
            output_fragments[var] = out

    # output_fragments now maps a TypeVarExp to a TypeExp which is the solved
    # fragment for the given output type variable

    results = {}
    for key, expr in outputs.items():
        r = expr  # output may not have a typevar, so default is the expr
        for selector in select_variables(expr):
            var = selector(expr)
            r = selector(r, swap=output_fragments[var])
        results[key] = r

    # results now maps a key to a full TypeExp as solved by the inputs
    return results
