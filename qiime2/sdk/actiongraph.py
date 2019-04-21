from itertools import combinations
import networkx as nx
import qiime2


def get_next_param(method, input=True):
    """
    Get a tuple of required/nonrequired inputs or outputs for each method

    Parameters
    ----------
    method : Qiime2.method (check the right signature for this)
        The method to get the semantic types from
    input : {True,False}
        Delineates if getting the method input or output types

    Returns
    -------
    Tuple - list of required and optional semantic types

    """

    req = []
    non_req = []

    if input:
        for k, v in method.signature.inputs.items():
            if not v.has_default():
                req.append(v.qiime_type)
            else:
                non_req.append(v.qiime_type)
        for k, v in method.signature.parameters.items():
            if not v.has_default():
                print("parameter has default!")
                req.append(v.qiime_type)
            else:
                print("parameter doesn't have default!")
                non_req.append(v.qiime_type)
    else:
        for k, v in method.signature.outputs.items():
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


def build_graph(sigs=[]):
    """
    Constructs a networkx graph with different semantic types
    and methods as nodes

    Parameters
    ----------
    sigs : list of strings

    Returns
    -------
    nx.DiGraph - networkx graph connected based on all or specified methods
    """
    G = nx.DiGraph()
    method_list = []

    pm = qiime2.sdk.PluginManager()
    if not sigs:
        for pgn, pg in pm.plugins.items():
            method_list += list(pg.actions.values())
    else:
        for pgn, pg in pm.plugins.items():
            for method in pg.actions.keys():
                if str(method) in sigs:
                    method_list.append(pg.actions[str(method)])

    for method in method_list:
        req_in, non_req_in = get_next_param(method, 1)
        req_out, non_req_out = get_next_param(method, 0)
        for key in req_in:
            if not G.has_node(key):
                G.add_node(key, value=key, color='blue')
        for key in non_req_in:
            if not G.has_node(key):
                G.add_node(key, value=key, color='blue')
        for key in req_out:
            if not G.has_node(key):
                G.add_node(key, value=key, color='blue')
        for key in non_req_out:
            if not G.has_node(key):
                G.add_node(key, value=key, color='blue')
        if non_req_in:
            combs = get_combinations(non_req_in)
            for non_req in combs:
                new_method_key = str(method)
                str_non_req = [str(i) for i in non_req]
                if str_non_req:
                    non_req_key = ";".join(str_non_req)
                    new_method_key += ";"+non_req_key
                if not G.has_node(new_method_key):
                    G.add_node(new_method_key, value=method, color='red')
                    for key in req_in:
                        G.add_edge(key, new_method_key)
                    #for key in non_req:
                    #    G.add_edge(key, new_method_key)
                    for key in req_out:
                        G.add_edge(new_method_key, key)
                    #for key in non_req_out:
                    #    G.add_edge(new_method_key, key)
        else:
            if not G.has_node(method):
                m = str(method)
                G.add_node(m, value=method, color='red')
                for key in req_in:
                    G.add_edge(key, m)
                for key in req_out:
                    G.add_edge(m, key)
                #for key in non_req_out:
                #    G.add_edge(m, key)
    return G
