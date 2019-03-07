import qiime2
import random
from itertools import combinations
import networkx as nx
import matplotlib.pyplot as plt


#q2methods = qiime2.sdk.util.actions_by_input_type('PCoAResults')
q2methods = qiime2.sdk.util.actions_by_input_type('FeatureTable[Frequency]')

#############testing##################
q2p = qiime2.sdk.PluginManager()

#for plugin, value in q2p.plugins.items():
#    print(value.methods.keys())
##############testing#################


def getNextParam(method, b):
    '''
        description: given an input method, return a tuple of parameter types
        input: method object
        output: tuple of parameter types
    '''
    req = []
    non_req = []
    if b:
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
    #print("req: "+str(req))
    #print("non_req: "+str(non_req))
    return req, non_req

def getCombinations(no_req):
    if not no_req:
        return no_req
    no_req_comb = []
    for i in range(len(no_req)+1):
        combs = combinations(no_req, i)
        for c in combs:
            no_req_comb.append(list(c))
    return no_req_comb


def buildGraph(G):
    for plugin in q2methods:
        method_list = plugin[1]
        for method in method_list:
            #check if there are any combinations in parameters
            req_in, non_req_in = getNextParam(method,1)
            req_out, non_req_out = getNextParam(method,0)

            #add keys to graph
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
                    #print("new method"+new_method_key)

                    #check if there exists a method node with same in-edges
                    if not G.has_node(new_method_key):
                        G.add_node(new_method_key,value=method,color='red')
                        for key in req_in:
                            G.add_edge(key, new_method_key)
                        for key in non_req:
                            G.add_edge(key, new_method_key)
                        for key in req_out:
                            G.add_edge(new_method_key, key)
                        for key in non_req_out:
                            G.add_edge(new_method_key, key)
            else:
                if not G.has_node(method):
                    G.add_node(method,value=method,color='red')
                    for key in req_in:
                        G.add_edge(key, method)
                    for key in req_out:
                        G.add_edge(method, key)
                    for key in non_req_out:
                        G.add_edge(method, key)


'''
if __name__ == '__main__':
    G=nx.DiGraph()
    c_map = [] #color map
    buildGraph(G)
    c_map = [n[1]['color'] for n in list(G.nodes(data=True))]
    print(c_map)
    print()
    nx.draw(G, node_color=c_map, with_labels=True)
    plt.draw()
    plt.show()
    plt.savefig('graph2.png')


    for i in G.nodes:
        print(i)
        print(G.in_edges(nbunch=i))
        print()
'''
