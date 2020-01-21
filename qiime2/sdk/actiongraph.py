# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from itertools import product, chain
import networkx as nx
import copy
import qiime2


def get_next_arguments(action, type="input"):
    """
    Get a tuple of required/nonrequired inputs or outputs for each method

    Parameters
    ----------
    action : Qiime2.action

    type : {"input", "param", "output"}
        Delineates if getting the action input, param, or output types

    Returns
    -------
    List of tuples containing name and required semantic types
    List of tuples containing name and optional semantic types
    """

    req = []
    non_req = []

    if type == "input":
        for k, v in action.signature.inputs.items():
            if not v.has_default():
                req.append([k, v.qiime_type])
            else:
                non_req.append(["."+k, v.qiime_type])
    elif type == "param":
        for k, v in action.signature.parameters.items():
            if not v.has_default():
                req.append([k, v.qiime_type])
            else:
                non_req.append(["."+k, v.qiime_type])
    else:
        for k, v in action.signature.outputs.items():
            if not v.has_default():
                req.append([k, v.qiime_type])
            else:
                non_req.append(["."+k, v.qiime_type])

    return req, non_req


def unravel(l):
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
        if len(list(x[1])) > 1:
            members = list(x[1])
            temp = copy.deepcopy(result)

            # update result with first element of types in member
            for each_list in result:
                each_list[i][1] = members[0]

            # add in other permutations of types in member
            for n in range(1, len(members)):
                copy_result = copy.deepcopy(temp)
                for each_list in copy_result:
                    each_list[i][1] = members[n]
                result += copy_result
    return result


def generate_nodes_by_action(action, opt=False):
    """
    Given a method, generates all combinations of inputs and
    outputs for that particular method and and stores the combinations
    as dictionaries in a resulting list.

    Parameters
    ----------
    method : Qiime2.action
    opt : {True, False}
        Delineates if optional types should be included

    Returns
    -------
    list of dictionaries - each dictionary is a combination inputs and
        outputs for particular node

    """

    input, input_nr = get_next_arguments(action, "input")
    param, param_nr = get_next_arguments(action, "param")
    output, output_nr = get_next_arguments(action, "output")

    input = unravel(input)
    param = unravel(param)

    opt_in_list = []

    if opt:
        opt_in_list += input_nr
        opt_in_list += param_nr
        opt_in_list = unravel(opt_in_list)
        ins = [dict(x) for x in
               [list(chain.from_iterable(i)) for i in
                list(product(input, param, opt_in_list))]]
        outs = dict(output + output_nr)
        results = [{'inputs': i, 'outputs': outs} for i in ins]
        return results

    ins = [dict(x) for x in
           [list(chain.from_iterable(i)) for i in list(product(input, param))]]
    outs = dict(output)
    results = [{'inputs': i, 'outputs': outs} for i in ins]
    return results


def build_graph(action_list=[], opt=False):
    """
    Constructs a networkx graph with different semantic types
    and actions as nodes

    Parameters
    ----------
    action_list : list of Qiime2.action
        If list is empty, will pull from all methods in the Qiime2 plugin
    opt : {True, False}
        Delineates if optional types should be included in the graph

    Returns
    -------
    nx.DiGraph - networkx graph connected based on all or specified methods

    """

    G = nx.DiGraph()
    G.edges(data=True)

    # get all actions or specifc actions if specified in sigs
    pm = qiime2.sdk.PluginManager()
    if not action_list:
        for _, pg in pm.plugins.items():
            action_list += list(pg.actions.values())

    for action in action_list:
        results = generate_nodes_by_action(action, opt)
        for dict_ in results:
            for k, v in dict_.items():
                if not v:
                    continue

                # renaming dictionary to remove '.'
                action_node = {}
                for x, y in v.items():
                    if x[0] == '.':
                        action_node[x[1:]] = y
                    else:
                        action_node[x] = y
                dict_[k] = action_node

                if not G.has_node(str(dict_)):
                    G.add_node(str(dict_), value=action, node='action')

                if k == 'inputs':
                    for in_k, in_v in v.items():
                        if not in_v:
                            continue
                        if in_k[0] == '.':
                            name = "opt_"+str(in_v)
                            G.add_edge(name, str(dict_))
                            G[name][str(dict_)]['name'] = in_k[1:]
                            G.nodes[name]['type'] = in_v
                            G.nodes[name]['optional'] = True
                            G.nodes[name]['node'] = 'type'
                        else:
                            G.add_edge(in_v, str(dict_))
                            G[in_v][str(dict_)]['name'] = in_k
                            G.nodes[in_v]['type'] = in_v
                            G.nodes[in_v]['optional'] = False
                            G.nodes[in_v]['node'] = 'type'
                else:
                    for out_k, out_v in v.items():
                        if not out_v:
                            continue
                        if out_k[0] == '.':
                            name = "opt_"+str(out_v)
                            G.add_edge("opt_"+str(out_v), str(dict_))
                            G[str(dict_)][name]['name'] = out_k[1:]
                            G.nodes[name]['type'] = in_v
                            G.nodes[name]['optional'] = True
                            G.nodes[name]['node'] = 'type'
                        else:
                            G.add_edge(str(dict_), out_v)
                            G[str(dict_)][out_v]['name'] = out_k
                            G.nodes[out_v]['type'] = out_v
                            G.nodes[out_v]['optional'] = False
                            G.nodes[out_v]['node'] = 'type'
    return G
