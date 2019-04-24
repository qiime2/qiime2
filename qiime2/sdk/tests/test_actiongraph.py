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
from qiime2.core.testing.method import (no_input_method,
                                        docstring_order_method,
                                        identity_with_metadata)
from qiime2.core.testing.visualizer import (params_only_viz, most_common_viz)
from qiime2.core.testing.pipeline import (visualizer_only_pipeline)
from qiime2.sdk.actiongraph import build_graph


class TestActiongraph(unittest.TestCase):
    def setUp(self):
        self.plugin = get_dummy_plugin()
        self.g = None

    def test_simple_graph(self):
        methods = [no_input_method]
        self.g = build_graph(methods)
        obs = list(self.g.nodes)

        exp_node = str({
            'inputs': {},
            'outputs': {
                'out': Mapping
            },
        })

        type_node = Mapping
        exp = [type_node, exp_node]

        for item in obs:
            assert item in exp

        assert self.g.has_edge(str(exp_node), type_node)

    def test_cycle_in_graph(self):
        methods = [docstring_order_method]
        self.g = build_graph(methods)
        obs = list(self.g.nodes)
        exp = [Mapping, Int, Str]

        exp_node_1 = str({
            'inputs': {
                'req_input': Mapping,
                'opt_input': Mapping,
                'req_param': Str,
                'opt_param': Int
            },
            'outputs': {
                'out': Mapping
            },
        })

        exp_node_2 = str({
            'inputs': {
                'req_input': Mapping,
                'opt_input': Mapping,
                'req_param': Str,
                'opt_param': None
            },
            'outputs': {
                'out': Mapping
            },
        })

        exp_node_3 = str({
            'inputs': {
                'req_input': Mapping,
                'opt_input': None,
                'req_param': Str,
                'opt_param': Int
            },
            'outputs': {
                'out': Mapping
            },
        })

        exp_node_4 = str({
            'inputs': {
                'req_input': Mapping,
                'opt_input': None,
                'req_param': Str,
                'opt_param': None
            },
            'outputs': {
                'out': Mapping
            },
        })

        exp += [exp_node_1, exp_node_2, exp_node_3, exp_node_4]

        for item in obs:
            assert item in exp

        assert self.g.in_degree(exp_node_1) == 3
        assert self.g.out_degree(exp_node_1) == 1

        assert self.g.in_degree(exp_node_2) == 2
        assert self.g.out_degree(exp_node_2) == 1

        assert self.g.in_degree(exp_node_3) == 3
        assert self.g.out_degree(exp_node_3) == 1

        assert self.g.in_degree(exp_node_4) == 2
        assert self.g.out_degree(exp_node_4) == 1

    def test_visualizer(self):
        vis = [params_only_viz]

        self.g = build_graph(vis)
        obs = list(self.g.nodes)
        exp = [Visualization, Int, Str]

        exp_node_1 = str({
            'inputs': {
                'name': Str,
                'age': Int
            },
            'outputs': {
                'visualization': Visualization
            },
        })

        exp_node_2 = str({
            'inputs': {
                'name': None,
                'age': Int
            },
            'outputs': {
                'visualization': Visualization
            },
        })

        exp_node_3 = str({
            'inputs': {
                'name': Str,
                'age': None
            },
            'outputs': {
                'visualization': Visualization
            },
        })

        exp_node_4 = str({
            'inputs': {
                'name': None,
                'age': None
            },
            'outputs': {
                'visualization': Visualization
            },
        })

        exp += [exp_node_1, exp_node_2, exp_node_3, exp_node_4]
        for item in obs:
            assert item in exp

        assert self.g.in_degree(exp_node_1) == 2
        assert self.g.out_degree(exp_node_1) == 1

        assert self.g.in_degree(exp_node_2) == 1
        assert self.g.out_degree(exp_node_2) == 1

        assert self.g.in_degree(exp_node_3) == 1
        assert self.g.out_degree(exp_node_3) == 1

        assert self.g.in_degree(exp_node_4) == 0
        assert self.g.out_degree(exp_node_4) == 1

        assert self.g.in_degree(Visualization) == 4
        assert self.g.out_degree(Visualization) == 0

    def test_union(self):
        vis = [most_common_viz]
        self.g = build_graph(vis)
        obs = list(self.g.nodes)
        exp = [Visualization, IntSequence1, IntSequence2]

        exp_node_1 = str({
            'inputs': {
                'ints': IntSequence1,
            },
            'outputs': {
                'visualization': Visualization
            },
        })
        exp_node_2 = str({
            'inputs': {
                'ints': IntSequence2,
            },
            'outputs': {
                'visualization': Visualization
            },
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
        actions = [visualizer_only_pipeline]
        self.g = build_graph(actions)
        obs = list(self.g.nodes)
        exp = [Visualization, Mapping]

        exp_node = str({
            'inputs': {
                'mapping': Mapping
            },
            'outputs': {
                'viz1': Visualization,
                'viz2': Visualization
            },
        })

        exp += [exp_node]

        for item in obs:
            assert item in exp

        assert self.g.in_degree(exp_node) == 1
        assert self.g.out_degree(exp_node) == 1

    def test_metadata(self):
        actions = [identity_with_metadata]
        self.g = build_graph(actions)
        obs = list(self.g.nodes)
        exp = [Metadata, IntSequence1, IntSequence2]

        exp_node_1 = str({
            'inputs': {
                'ints': IntSequence1,
                'metadata': Metadata
            },
            'outputs': {
                'out': IntSequence1
            },
        })

        exp_node_2 = str({
            'inputs': {
                'ints': IntSequence2,
                'metadata': Metadata
            },
            'outputs': {
                'out': IntSequence1
            },
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
