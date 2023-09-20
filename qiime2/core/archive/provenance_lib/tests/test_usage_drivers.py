import os
import shutil
import tempfile
import unittest

from qiime2.sdk.plugin_manager import PluginManager
from qiime2.core.testing.type import IntSequence1

from ..replay import replay_provenance


class ReplayPythonUsageTests(unittest.TestCase):
    def setUp(self):
        self.pm = PluginManager()
        self.dp = self.pm.plugins['dummy-plugin']
        self.tempdir = tempfile.mkdtemp(
            prefix='qiime2-test-usage-drivers-temp-'
        )

        def return_many_ints() -> (list, list, list, list, list, list):
            return ([1, 2, 3], [4, 5, 6], [7], [4, 4], [0], [9, 8])

        self.dp.methods.register_function(
            function=return_many_ints,
            inputs={},
            parameters={},
            outputs=[
                ('ints1', IntSequence1),
                ('ints2', IntSequence1),
                ('ints3', IntSequence1),
                ('ints4', IntSequence1),
                ('ints5', IntSequence1),
                ('ints6', IntSequence1),
            ],
            output_descriptions={
                'ints1': 'ints',
                'ints2': 'ints',
                'ints3': 'ints',
                'ints4': 'ints',
                'ints5': 'ints',
                'ints6': 'ints',
            },
            name='return_many_ints',
            description=''
        )

        def return_four_ints() -> (list, list, list, list):
            return ([1, 2, 3], [4, 5, 6], [7, 8, 9], [4, 4])

        self.dp.methods.register_function(
            function=return_four_ints,
            inputs={},
            parameters={},
            outputs=[
                ('ints1', IntSequence1),
                ('ints2', IntSequence1),
                ('ints3', IntSequence1),
                ('ints4', IntSequence1),
            ],
            output_descriptions={
                'ints1': 'ints',
                'ints2': 'ints',
                'ints3': 'ints',
                'ints4': 'ints',
            },
            name='return_four_ints',
            description=''
        )

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_template_action_lumps_many_outputs(self):
        """
        ReplayPythonUsage._template_action should "lump" multiple outputs from
        one command into a single Results-like object when the total number of
        outputs from a single command > 5

        In these cases, our rendering should look like:
        `action_results = plugin_actions.action()...`
        instead of:
        `_, _, thing3, _, _, _ = plugin_actions.action()...`
        """
        ints = self.dp.actions['return_many_ints']()
        first_ints = ints[0]
        first_ints.save(os.path.join(self.tempdir, 'int-seq.qza'))
        fp = os.path.join(self.tempdir, 'int-seq.qza')
        out_fp = os.path.join(self.tempdir, 'action_collection.txt')
        replay_provenance(fp, out_fp, 'python3')

        exp = 'action_results = dummy_plugin_actions.return_many_ints'
        with open(out_fp) as fh:
            rendered = fh.read()
        self.assertRegex(rendered, exp)

    def test_template_action_does_not_lump_four_outputs(self):
        """
        ReplayPythonUsage._template_action should not "lump" multiple outputs
        one command into a single Results-like object when the total number of
        outputs from a single command <= 5, unless the total number of results
        is high (see above).

        In these cases, our rendering should look like:
        `_, _, thing3, _ = plugin_actions.action()...`
        instead of:
        `action_results = plugin_actions.action()...`
        """
        ints = self.dp.actions['return_four_ints']()
        first_ints = ints[0]
        first_ints.save(os.path.join(self.tempdir, 'int-seq.qza'))
        fp = os.path.join(self.tempdir, 'int-seq.qza')
        out_fp = os.path.join(self.tempdir, 'action_collection.txt')
        replay_provenance(fp, out_fp, 'python3')

        exp = 'ints1_0, _, _, _ = dummy_plugin_actions.return_four_ints'
        with open(out_fp) as fh:
            rendered = fh.read()
        self.assertRegex(rendered, exp)

    def test_template_action_lumps_three_variables(self):
        """
        ReplayPythonUsage._template_action should "lump" multiple outputs from
        one command into a single Results-like object when there are more than
        two usage variables (i.e. replay of 3+ results from a single command)

        In these cases, our rendering should look like:
        ```
        action_results = plugin_actions.action(...)
        thing1 = action_results.thinga
        etc.
        ```
        instead of:
        `thing1, _, thing3, _, thing5, _ = plugin_actions.action()...`
        """
        ints = self.dp.actions['return_four_ints']()
        os.mkdir(os.path.join(self.tempdir, 'three-ints-dir'))
        for i in range(3):
            out_path = os.path.join(self.tempdir, 'three-ints-dir',
                                    f'int-seq-{i}.qza')
            ints[i].save(out_path)

        fp = os.path.join(self.tempdir, 'three-ints-dir')
        out_fp = os.path.join(self.tempdir, 'action_collection.txt')
        replay_provenance(fp, out_fp, 'python3')

        exp = (
            'action_results = dummy_plugin_actions.return_four_ints',
            'ints1_0 = action_results.ints1',
            'ints2_0 = action_results.ints2',
            'ints3_0 = action_results.ints3'
        )
        with open(out_fp) as fh:
            rendered = fh.read()
        for pattern in exp:
            self.assertRegex(rendered, pattern)

    def test_template_action_does_not_lump_two_vars(self):
        """
        ReplayPythonUsage._template_action should not "lump" multiple outputs
        from one command into a single Results-like object when the total count
        of usage variables (i.e. replayed outputs) from a single command < 3,
        unless the total number of outputs is high (see above).

        In these cases, our rendering should look like:
        `thing1, _, thing3, _ = plugin_actions.action()...`
        instead of:
        `action_results = plugin_actions.action()...`
        """
        ints1, ints2, _, _ = self.dp.actions['return_four_ints']()
        os.mkdir(os.path.join(self.tempdir, 'two-ints-dir'))
        ints1.save(os.path.join(self.tempdir, 'two-ints-dir', 'int-seq-1.qza'))
        ints2.save(os.path.join(self.tempdir, 'two-ints-dir', 'int-seq-2.qza'))
        fp = os.path.join(self.tempdir, 'two-ints-dir')
        out_fp = os.path.join(self.tempdir, 'action_collection.txt')
        replay_provenance(fp, out_fp, 'python3')

        exp = 'ints1_0, ints2_0, _, _ = dummy_plugin_actions.return_four_ints'
        with open(out_fp) as fh:
            rendered = fh.read()
        self.assertRegex(rendered, exp)
