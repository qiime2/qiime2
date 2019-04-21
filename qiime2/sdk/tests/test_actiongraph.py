import unittest
import qiime2.sdk
from qiime2.core.testing.method import (split_ints, params_only_method,
                                        no_input_method, long_description_method)

from qiime2.core.testing.type import (Mapping)
from qiime2.core.type.primitive import (Int, Str)
from qiime2.core.testing.util import get_dummy_plugin
from actiongraph import *

class TestActiongraph(unittest.TestCase):
    def setUp(self):
        self.plugin = get_dummy_plugin()
        self.g = None

    def test_no_graph(self):
        methods = ['not_a_valid_method']
        self.g = build_graph(methods)
        obs = list(self.g.nodes)
        self.assertEqual(obs, [])
      
    def test_simple_graph(self):
        methods = ['no_input_method']
        self.g = build_graph(methods)
        obs = list(self.g.nodes)
        exp = [Mapping, str(self.plugin.methods['no_input_method'])]
        for item in obs:
            assert item in exp

        no_input_node = str(self.plugin.methods['no_input_method'])
        type_node = Mapping
        assert self.g.has_edge(no_input_node, type_node)

    def test_cycle_in_graph(self):
        methods = ['long_description_method','docstring_order_method'] 
        self.g = build_graph(methods)
        obs = list(self.g.nodes)
        exp = [Mapping, Int, Str, 'oInt', 'oMapping', 
            str(self.plugin.methods['long_description_method']),
            str(self.plugin.methods['docstring_order_method']), 
            str(self.plugin.methods['docstring_order_method'])+';'+str(Mapping), 
            str(self.plugin.methods['docstring_order_method'])+';'+str(Mapping)+';'+'Int',
            str(self.plugin.methods['docstring_order_method'])+';Int']

        for item in obs:
            assert item in exp
        
        ds = str(self.plugin.methods['docstring_order_method'])
        ds_map = str(self.plugin.methods['docstring_order_method'])+';'+str(Mapping)
        ds_map_int = str(self.plugin.methods['docstring_order_method'])+';'+\
        str(Mapping)+';Int'
        ds_int = str(self.plugin.methods['docstring_order_method'])+';Int'
        long_node = str(self.plugin.methods['long_description_method'])
        map_node = Mapping
        int_node = Int
        str_node = Str

        assert self.g.has_edge(int_node, long_node)
        assert self.g.has_edge(str_node, long_node)
        assert self.g.has_edge(map_node, long_node)
        assert self.g.has_edge(long_node, map_node)
        
        assert self.g.in_degree(ds) == 2
        assert self.g.out_degree(ds) == 1
    
        assert self.g.in_degree(ds_map) == 2
        assert self.g.out_degree(ds_map) == 1
        
        assert self.g.in_degree(ds_int) == 2
        assert self.g.out_degree(ds_int) == 1
        
        assert self.g.in_degree(ds_map_int) == 2
        assert self.g.out_degree(ds_map_int) == 1

if __name__ == '__main__':
    unittest.main()
