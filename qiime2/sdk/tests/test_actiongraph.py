import unittest
import qiime2.sdk
from qiime2.core.testing.method import (split_ints, params_only_method,
                                        no_input_method, long_description_method)

from qiime2.core.testing.type import (Mapping)
from qiime2.core.testing.util import get_dummy_plugin
from actiongraph import *

class TestUtil(unittest.TestCase):
    def setUp(self):
        self.plugin = get_dummy_plugin()
        self.g = None


    def test_build_graph(self):
        #set up graph to contain multiple methods
        methods = ['long_description_method', 'no_input_method', 
                   'params_only_method', 'docstring_order_method'] 
    
        all_methods = list(self.plugin.methods.keys())
        for method in all_methods:
            if method not in methods:
                del self.plugin.methods[method]

        #test 
        self.g = build_graph()
        obs = list(self.g.nodes)
        exp = [Mapping, str(self.plugin.methods['no_input_method']), 
               str(self.plugin.methods['long_description_method']),
               str(self.plugin.methods['params_only_method']),
               str(self.plugin.methods['docstring_order_method']), 
               str(self.plugin.methods['docstring_order_method'])+';'+str(Mapping)]

        for item in exp:
            assert item in obs
        
        docstring_node = str(self.plugin.methods['docstring_order_method'])
        docstring_node_2 = str(self.plugin.methods['docstring_order_method'])+';'+str(Mapping)
        long_node = str(self.plugin.methods['long_description_method'])
        params_node = str(self.plugin.methods['params_only_method'])
        no_input_node = str(self.plugin.methods['no_input_method'])
        type_node = Mapping

        assert self.g.has_edge(docstring_node, type_node)
        assert self.g.has_edge(type_node, docstring_node)
        assert self.g.has_edge(docstring_node_2, type_node)
        assert self.g.has_edge(type_node, docstring_node)
        assert self.g.has_edge(long_node, type_node)
        assert self.g.has_edge(type_node, long_node)
        assert self.g.has_edge(params_node, type_node)
        assert self.g.has_edge(no_input_node, type_node)
    

        #set up graph to contain multiple methods
        methods = ['long_description_method', 'no_input_method', 
                   'params_only_method'] 
    
        all_methods = list(self.plugin.methods.keys())
        for method in all_methods:
            if method not in methods:
                del self.plugin.methods[method]

        #test 1
        self.g = build_graph()
        obs = list(self.g.nodes)
        exp = [Mapping, str(self.plugin.methods['no_input_method']), 
               str(self.plugin.methods['long_description_method']),
               str(self.plugin.methods['params_only_method'])]

        for item in exp:
            assert item in obs
        
        long_node = str(self.plugin.methods['long_description_method'])
        params_node = str(self.plugin.methods['params_only_method'])
        no_input_node = str(self.plugin.methods['no_input_method'])
        type_node = Mapping

        assert self.g.has_edge(long_node, type_node)
        assert self.g.has_edge(type_node, long_node)
        assert self.g.has_edge(params_node, type_node)
        assert self.g.has_edge(no_input_node, type_node)
    
        #set up graph to only contain a base method
        methods = ['no_input_method']
        all_methods = list(self.plugin.methods.keys())
        for method in all_methods:
            if method not in methods:
                del self.plugin.methods[method]

        #test 2
        self.g = build_graph()
        obs = list(self.g.nodes)
        exp = [Mapping, str(self.plugin.methods['no_input_method'])]
        for item in exp:
            assert item in obs

        no_input_node = str(self.plugin.methods['no_input_method'])
        type_node = Mapping
        assert self.g.has_edge(no_input_node, type_node)

        #set up graph to contain no methods
        methods = []
        all_methods = list(self.plugin.methods.keys())
        for method in all_methods:
            if method not in methods:
                del self.plugin.methods[method]
        
        #test 3
        self.g = build_graph()
        obs = list(self.g.nodes)
        self.assertEqual(obs, [])
        

if __name__ == '__main__':
    unittest.main()
