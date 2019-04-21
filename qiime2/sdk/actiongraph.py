# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


from itertools import combinations
import networkx as nx
import copy
import qiime2


def get_next_arguments(action, type=0):
    """
    Get a tuple of required/nonrequired inputs or outputs for each method

    Parameters
    ----------
    action : Qiime2.action

    type : {0,1,2}
        Delineates if getting the action input, param, or output types

    Returns
    -------
    Tuple - list of required and optional semantic types

    """

    req = []
    non_req = []

    if type == 0:
        for k, v in action.signature.inputs.items():
            if not v.has_default():
                req.append(v.qiime_type)
            else:
                non_req.append(v.qiime_type)
    elif type == 1:
        for k, v in action.signature.parameters.items():
            if not v.has_default():
                req.append(v.qiime_type)
            else:
                non_req.append(v.qiime_type)
    else:
        for k, v in action.signature.outputs.items():
            if not v.has_default():
                req.append(v.qiime_type)
            else:
                non_req.append(v.qiime_type)

    return req, non_req


def get_combinations(no_req):
    """
    Get all permutations of list of semantic types

    Parameters
    ----------
    non_req : list of semantic types

    Returns
    -------
    List of lists - each list contains an instance of a permutation

    """

    if not no_req:
        return no_req
    no_req_comb = []
    for i in range(len(no_req)+1):
        combs = combinations(no_req, i)
        for c in combs:
            no_req_comb.append(list(c))
    return no_req_comb


def generate_nodes_by_action(action):
    """
    Given a method, generates all combinations of inputs and
    outputs for that particular method and and stores the combinations
    as dictionaries in a resulting list.

    Parameters
    ----------
    method : Qiime2.action

    Returns
    -------
    list of dictionaries - each dictionary is a combination inputs and
        outputs for particular node

    """

    results = []

    input, input_nr = get_next_arguments(action, 0)
    param, param_nr = get_next_arguments(action, 1)
    output, _ = get_next_arguments(action, 2)

    non_req = []

    # unravel potential unions
    input = unravel_union(input)
    param = unravel_union(param)

    for x in input_nr + param_nr:
        non_req += unravel_union(x)

    for i in input:
        for p in param:
            if non_req:
                non_req = get_combinations(non_req)
                for nr in non_req:
                    d = {}
                    d['inputs'] = i
                    d['params'] = p
                    d['outputs'] = output
                    d['non_req'] = nr
                    results.append(d)
            else:
                d = {}
                d['inputs'] = i
                d['params'] = p
                d['outputs'] = output
                d['non_req'] = []
                results.append(d)

    return results


def unravel_union(l):
    """
    Unravel Union node to get all permutations of types for each action

    Parameters
    ----------
    list : list of Qiime2.types

    Returns
    -------
    list of lists - list of permuations of types for each action

    """

    result = [l]
    for i, x in enumerate(l):
        if 'members' in x.__dict__ and len(x.__dict__['members']) > 0:
            members = list(x.__dict__['members'])
            temp = copy.deepcopy(result)

            # update result with first element of types in member
            for each_list in result:
                each_list[i] = members[0]

            # add in other permutations of types in member
            for n in range(1, len(members)):
                copy_result = copy.deepcopy(temp)
                for each_list in copy_result:
                    each_list[i] = members[n]
                result += copy_result
    return result


def build_graph(sigs=[]):
    """
    Constructs a networkx graph with different semantic types
    and actions as nodes

    Parameters
    ----------
    sigs : list of strings

    Returns
    -------
    nx.DiGraph - networkx graph connected based on all or specified methods

    """

    G = nx.DiGraph()
    action_list = []

    # get all actions or specifc actions if specified in sigs
    pm = qiime2.sdk.PluginManager()
    if not sigs:
        for pgn, pg in pm.plugins.items():
            action_list += list(pg.actions.values())
    else:
        for pgn, pg in pm.plugins.items():
            for action in pg.actions.keys():
                if str(action) in sigs:
                    action_list.append(pg.actions[str(action)])

    for action in action_list:
        node_combs = generate_nodes_by_action(action)

        # each action node is a str(dict)
        # each Qiime2.type node is the type itself
        for dict in node_combs:
            if not G.has_node(str(dict)):
                G.add_node(str(dict), value=action)

            for k in dict.keys():
                for type in dict[k]:
                    if not G.has_node(type):
                        G.add_node(type, value=type)
                    if k == 'outputs':
                        G.add_edge(str(dict), type)
                    elif k == 'inputs':
                        G.add_edge(type, str(dict))
                    elif k == 'params':
                        G.add_edge(type, str(dict))
    return G
