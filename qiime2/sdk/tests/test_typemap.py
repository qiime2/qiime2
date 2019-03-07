import qiime2
import qiime2.sdk
import networkx as nx

q2methods = qiime2.sdk.util.actions_by_input_type('FeatureTable[Frequency]')
G=nx.DiGraph()
buildGraph(G)
cmap=[n[1]['color'] for n in list(G.nodes(data=True))]

def test_get_input_keys():
    print("Testing Input Keys")
    for plugin in q2methods:
        method_list = plugin[1]
        for method in method_list:
            req_in, _ = getNextParam(method,1)
            for input in req_in:
                assert G.has_node(input), "Should have this semantic node"

def test_get_input_combination_keys():
    print("Testing Input Combination Keys")
    for plugin in q2methods:
        method_list = plugin[1]
        for method in method_list:
            _, non_req_in = getNextParam(method,1)
            combs= getCombinations(non_req_in)
            for non_req in combs:
                new_method_key = str(method)
                str_non_req = [str(i) for i in non_req]
                if str_non_req:
                    non_req_key = ";".join(str_non_req)
                    new_method_key += ";"+non_req_key
                assert G.has_node(new_method_key), "Should have this combination node"

def test_get_output_keys():
    print("Testing Input Keys")
    for plugin in q2methods:
        method_list = plugin[1]
        for method in method_list:
            req_out, _ = getNextParam(method,0)
            for output in req_out:
                assert G.has_node(output), "Should have this semantic node"

def test_color_keys():
    print("Testing Color Keys")
    all_keys = {}
    for plugin in q2methods:
        method_list = plugin[1]
        for method in method_list:
            req_in, non_req_in = getNextParam(method,1)
            req_out, _ = getNextParam(method,0)

            for input in req_in:
                if input not in all_keys:
                    all_keys[input] = 1

            for output in req_out:
                if output not in all_keys:
                    all_keys[output] = 1

    #count number of blues in cmap
    num_blues = len([x for x in cmap if x=='blue'])
    num_keys = len(all_keys)
    assert num_blues == num_keys, "Number of keys incorrect"



def test_color_methods():
    print("Testing Color Methods")
    all_methods = {}
    for plugin in q2methods:
        method_list = plugin[1]
        for method in method_list:

            _, non_req_in = getNextParam(method,1)

            combs=getCombinations(non_req_in)
            for non_req in combs:
                new_method_key = str(method)
                str_non_req = [str(i) for i in non_req]
                if str_non_req:
                    non_req_key = ";".join(str_non_req)
                    new_method_key += ";"+non_req_key
                if new_method_key not in all_methods:
                    all_methods[new_method_key] = 1

            if str(method) not in all_methods:
                all_methods[method] = 1


    #count number of red in cmap
    num_reds = len([x for x in cmap if x=='red'])
    num_methods = len(all_methods)
    m_num = []
    for plugin in q2methods:
        m_num.append(len(plugin[1]))
    assert num_reds == num_methods, "Number of keys incorrect"

if __name__ == "__main__":
    test_get_input_keys()
    test_get_input_combination_keys()
    test_get_output_keys()
    test_color_keys()
    test_color_methods()
    print("Everything passed")
