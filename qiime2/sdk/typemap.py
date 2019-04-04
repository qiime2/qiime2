import random
from itertools import combinations
import networkx as nx
import matplotlib.pyplot as plt

import qiime2 

#q2methods = qiime2.sdk.util.actions_by_input_type('FeatureTable[Frequency]')

#############testing##################
q2p = qiime2.sdk.PluginManager()

#for plugin, value in q2p.plugins.items():
#    print(value.methods.keys())
##############testing#################

def getNextParam(method, input=True):
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
        for k,v in method.signature.inputs.items():
            if not v.has_default():
                req.append(v.qiime_type)
            else:
                non_req.append(v.qiime_type)
    else:
        for k,v in method.signature.outputs.items():
            if not v.has_default():
                req.append(v.qiime_type)
            else:
                non_req.append(v.qiime_type)
    return req, non_req

def getCombinations(no_req):
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



class TypeMap:
    """
    This class is is used to generate a graph of all methods in the Qiime2
    plugin and to create random testing paths from the methods.

    """
    def __init__(self, types):
        """
        Constructor method for TypeMap

        Parameters
        ----------
        types : list of strings of all semantic types to add to the graph

        """
        self.G = nx.DiGraph()
        self.methods = []
        for type in types:
            self.methods += qiime2.sdk.util.actions_by_input_type(type)
 

    def buildGraph(self):
        """Constructs a networkx graph with different semantic types 
        and methods as nodes
    
        """

        for plugin in self.methods:
            method_list = plugin[1]
            for method in method_list:
                req_in, non_req_in = getNextParam(method,1)
                req_out, non_req_out = getNextParam(method,0)
            
                for key in req_in:
                    if not self.G.has_node(key):
                        self.G.add_node(key, value=key, color='blue')
                for key in non_req_in:
                    if not self.G.has_node(key):
                        self.G.add_node(key, value=key, color='blue')
                for key in req_out:
                    if not self.G.has_node(key):
                        self.G.add_node(key, value=key, color='blue')
                for key in non_req_out:
                    if not self.G.has_node(key):
                        self.G.add_node(key, value=key, color='blue')

                #there are combinations
                if non_req_in:
                    #construct combinations
                    #print("COMBINATIONS: ")
                    #print(getCombinations(non_req_in))
                    combs = getCombinations(non_req_in)
                
                    for non_req in combs:
                        #print("DEBUG")
                        #print(non_req)
                        new_method_key= str(method)
                        str_non_req = [str(i) for i in non_req]
                        if str_non_req:
                            non_req_key = ";".join(str_non_req)
                            new_method_key += ";"+non_req_key
                    
                    #check if there exists a method node with same in-edges
                        if not self.G.has_node(new_method_key):
                            self.G.add_node(new_method_key,value=method,color='red')
                            for key in req_in:
                                self.G.add_edge(key, new_method_key)
                            for key in non_req:
                                self.G.add_edge(key, new_method_key)
                            for key in req_out:
                                self.G.add_edge(new_method_key, key)
                            for key in non_req_out:
                                self.G.add_edge(new_method_key, key)
                else:
                    if not self.G.has_node(method):
                        self.G.add_node(method,value=method,color='red')
                        for key in req_in:
                            self.G.add_edge(key, method)
                        for key in req_out:
                            self.G.add_edge(method, key)
                        for key in non_req_out:
                            self.G.add_edge(method, key)
                            








