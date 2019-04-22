# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
from qiime2.core.testing.type import (Mapping, IntSequence1, IntSequence2)
from qiime2.core.type.primitive import (Int, Str, Metadata)
from qiime2.core.type.visualization import (Visualization)
from qiime2.core.testing.util import get_dummy_plugin
from qiime2.sdk.actiongraph import build_graph


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

        exp_node = str({
            "inputs": [],
            "params": [],
            "outputs": [Mapping],
            "non_req": []
        })

        type_node = Mapping
        exp = [type_node, exp_node]

        for item in obs:
            assert item in exp

        assert self.g.has_edge(exp_node, type_node)

    def test_cycle_in_graph(self):
        methods = ['docstring_order_method']
        self.g = build_graph(methods)
        obs = list(self.g.nodes)
        exp = [Mapping, Int, Str]

        exp_node_1 = str({
            "inputs": [Mapping],
            "params": [Str],
            "outputs": [Mapping],
            "non_req": []
        })

        exp_node_2 = str({
            "inputs": [Mapping],
            "params": [Str],
            "outputs": [Mapping],
            "non_req": [Int]
        })

        exp_node_3 = str({
            "inputs": [Mapping],
            "params": [Str],
            "outputs": [Mapping],
            "non_req": [Mapping]
        })

        exp_node_4 = str({
            "inputs": [Mapping],
            "params": [Str],
            "outputs": [Mapping],
            "non_req": [Mapping, Int]
        })

        exp += [exp_node_1, exp_node_2, exp_node_3, exp_node_4]

        for item in obs:
            assert item in exp

        assert self.g.in_degree(exp_node_1) == 2
        assert self.g.out_degree(exp_node_1) == 1

        assert self.g.in_degree(exp_node_2) == 2
        assert self.g.out_degree(exp_node_2) == 1

        assert self.g.in_degree(exp_node_3) == 2
        assert self.g.out_degree(exp_node_3) == 1

        assert self.g.in_degree(exp_node_4) == 2
        assert self.g.out_degree(exp_node_4) == 1

    def test_union(self):
        vis = ['most_common_viz']
        self.g = build_graph(vis)
        obs = list(self.g.nodes)
        exp = [Visualization, IntSequence1, IntSequence2]

        exp_node_1 = str({
            "inputs": [IntSequence1],
            "params": [],
            "outputs": [Visualization],
            "non_req": []
        })

        exp_node_2 = str({
            "inputs": [IntSequence2],
            "params": [],
            "outputs": [Visualization],
            "non_req": []
        })

        exp += [exp_node_1, exp_node_2]

        for item in obs:
            assert item in exp

        assert self.g.in_degree(exp_node_1) == 1
        assert self.g.out_degree(exp_node_1) == 1

        assert self.g.in_degree(exp_node_2) == 1
        assert self.g.out_degree(exp_node_2) == 1

        assert self.g.in_degree(Visualization) == 2
        assert self.g.out_degree(Visualization) == 0

    def test_multiple_outputs(self):
        actions = ['visualizer_only_pipeline']
        self.g = build_graph(actions)
        obs = list(self.g.nodes)
        exp = [Visualization, Mapping]

        exp_node_1 = str({
            "inputs": [Mapping],
            "params": [],
            "outputs": [Visualization, Visualization],
            "non_req": []
        })

        exp += [exp_node_1]

        for item in obs:
            assert item in exp

        assert self.g.in_degree(exp_node_1) == 1
        assert self.g.out_degree(exp_node_1) == 1

    def test_metadata(self):
        actions = ['identity_with_metadata']
        self.g = build_graph(actions)
        obs = list(self.g.nodes)
        exp = [Metadata, IntSequence1, IntSequence2]

        exp_node_1 = str({
            "inputs": [IntSequence1],
            "params": [Metadata],
            "outputs": [IntSequence1],
            "non_req": []
        })

        exp_node_2 = str({
            "inputs": [IntSequence2],
            "params": [Metadata],
            "outputs": [IntSequence1],
            "non_req": []
        })

        exp += [exp_node_1, exp_node_2]

        for item in obs:
            assert item in exp

        assert self.g.in_degree(exp_node_1) == 2
        assert self.g.out_degree(exp_node_1) == 1

        assert self.g.in_degree(exp_node_1) == 2
        assert self.g.out_degree(exp_node_1) == 1

        assert self.g.in_degree(IntSequence1) == 2
        assert self.g.out_degree(IntSequence1) == 1


if __name__ == '__main__':
    unittest.main()
