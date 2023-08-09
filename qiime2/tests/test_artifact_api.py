# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib
import sys
import tempfile
import types
import unittest

import qiime2.sdk
from qiime2.core.testing.util import get_dummy_plugin
from qiime2.plugins import ArtifactAPIUsage


class TestImports(unittest.TestCase):
    def setUp(self):
        self._sys_modules = sys.modules.copy()
        # Ignore the returned dummy plugin object, just run this to verify the
        # plugin exists as the tests rely on it being loaded.
        get_dummy_plugin()

    def tearDown(self):
        # This allows us to reset the state of our imports to test different
        # expressions to make sure they work independently.
        remove = []
        for key in sys.modules:
            if key not in self._sys_modules:
                remove.append(key)  # can't delete in place while we iterate

        for key in remove:
            del sys.modules[key]

    def _check_spec(self, module, name, is_package):
        self.assertIsInstance(module, types.ModuleType)
        self.assertEqual(module.__name__, name)
        self.assertEqual(module.__spec__.name, name)
        self.assertEqual(module.__spec__.submodule_search_locations,
                         [] if is_package else None)
        self.assertFalse(module.__spec__.has_location)
        self.assertIn('generated QIIME 2 API', repr(module))

    def _check_plugin(self, module):
        self._check_spec(module, 'qiime2.plugins.dummy_plugin', True)
        self._check_methods(module.methods)
        self._check_visualizers(module.visualizers)
        self.assertEqual(set(x for x in dir(module) if not x.startswith('_')),
                         {'visualizers', 'methods', 'actions', 'pipelines'})

    def _check_methods(self, module):
        self._check_spec(module, 'qiime2.plugins.dummy_plugin.methods', False)
        self.assertTrue(hasattr(module, 'concatenate_ints'))
        self.assertFalse(hasattr(module, 'most_common_viz'))

        self.assertIsInstance(module.concatenate_ints, qiime2.sdk.Action)

    def _check_visualizers(self, module):
        self._check_spec(
            module, 'qiime2.plugins.dummy_plugin.visualizers', False)
        self.assertTrue(hasattr(module, 'most_common_viz'))
        self.assertFalse(hasattr(module, 'concatenate_ints'))
        self.assertIsInstance(module.most_common_viz, qiime2.sdk.Action)

    def test_import_root(self):
        import qiime2.plugins.dummy_plugin
        self._check_plugin(qiime2.plugins.dummy_plugin)

    def test_import_root_from(self):
        from qiime2.plugins import dummy_plugin
        self._check_plugin(dummy_plugin)

    def test_import_methods(self):
        import qiime2.plugins.dummy_plugin.methods
        self._check_methods(qiime2.plugins.dummy_plugin.methods)

    def test_import_visualizers(self):
        import qiime2.plugins.dummy_plugin.visualizers
        self._check_visualizers(qiime2.plugins.dummy_plugin.visualizers)

    def test_import_methods_from(self):
        from qiime2.plugins.dummy_plugin import methods
        self._check_methods(methods)

    def test_import_visualizers_from(self):
        from qiime2.plugins.dummy_plugin import visualizers
        self._check_visualizers(visualizers)

    def test_import_non_plugin(self):
        with self.assertRaises(ImportError):
            import qiime2.plugins.dummy_not_plugin  # noqa

    def test_import_non_action(self):
        with self.assertRaises(ImportError):
            import qiime2.plugins.dummy_plugin.non_action  # noqa

    def test_import_side_module(self):
        # Certain implementations of __PATH__ can cause a module to load
        # siblings (__PATH__ = ['.'] for example)
        import qiime2.metadata
        self.assertIsInstance(qiime2.metadata, types.ModuleType)
        with self.assertRaises(ImportError):
            import qiime2.plugins.metadata  # noqa

    def test_import_too_deep(self):
        with self.assertRaises(ImportError):
            import qiime2.plugins.dummy_plugin.methods.too_deep  # noqa

    def test_import_non_module(self):
        with self.assertRaises(ImportError):
            import qiime2.plugins.dummy_plugin.methods.concatenate_ints  # noqa

    def test_reload_fails(self):
        import qiime2.plugins.dummy_plugin
        with self.assertRaises(ImportError):
            importlib.reload(qiime2.plugins.dummy_plugin)


class TestArtifactAPIUsage(unittest.TestCase):
    def setUp(self):
        # TODO standardize temporary directories created by QIIME 2
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        self.plugin = get_dummy_plugin()

    def tearDown(self):
        self.test_dir.cleanup()

    def test_basic(self):
        action = self.plugin.actions['concatenate_ints']
        use = ArtifactAPIUsage()
        action.examples['concatenate_ints_simple'](use)
        exp = """\
import qiime2.plugins.dummy_plugin.actions as dummy_plugin_actions

# This example demonstrates basic usage.
ints_d, = dummy_plugin_actions.concatenate_ints(
    ints1=ints_a,
    ints2=ints_b,
    ints3=ints_c,
    int1=4,
    int2=2,
)"""
        self.assertEqual(exp, use.render())

    def test_chained(self):
        action = self.plugin.actions['concatenate_ints']
        use = ArtifactAPIUsage()
        action.examples['concatenate_ints_complex'](use)
        exp = """\
import qiime2.plugins.dummy_plugin.actions as dummy_plugin_actions

# This example demonstrates chained usage (pt 1).
ints_d, = dummy_plugin_actions.concatenate_ints(
    ints1=ints_a,
    ints2=ints_b,
    ints3=ints_c,
    int1=4,
    int2=2,
)
# This example demonstrates chained usage (pt 2).
concatenated_ints, = dummy_plugin_actions.concatenate_ints(
    ints1=ints_d,
    ints2=ints_b,
    ints3=ints_c,
    int1=41,
    int2=0,
)"""
        self.assertEqual(exp, use.render())

    def test_dereferencing(self):
        action = self.plugin.actions['typical_pipeline']
        use = ArtifactAPIUsage()
        action.examples['typical_pipeline_simple'](use)
        exp = """\
import qiime2.plugins.dummy_plugin.actions as dummy_plugin_actions

action_results = dummy_plugin_actions.typical_pipeline(
    int_sequence=ints,
    mapping=mapper,
    do_extra_thing=True,
)
out_map = action_results.out_map
left = action_results.left
right = action_results.right
left_viz_viz = action_results.left_viz
right_viz_viz = action_results.right_viz"""
        self.assertEqual(exp, use.render())

    def test_chained_dereferencing(self):
        action = self.plugin.actions['typical_pipeline']
        use = ArtifactAPIUsage()
        action.examples['typical_pipeline_complex'](use)
        exp = """\
import qiime2.plugins.dummy_plugin.actions as dummy_plugin_actions

action_results = dummy_plugin_actions.typical_pipeline(
    int_sequence=ints1,
    mapping=mapper1,
    do_extra_thing=True,
)
out_map1 = action_results.out_map
left1 = action_results.left
right1 = action_results.right
left_viz1_viz = action_results.left_viz
right_viz1_viz = action_results.right_viz
action_results = dummy_plugin_actions.typical_pipeline(
    int_sequence=left1,
    mapping=out_map1,
    do_extra_thing=False,
)
out_map2 = action_results.out_map
left2 = action_results.left
right2 = action_results.right
left_viz2_viz = action_results.left_viz
right_viz2_viz = action_results.right_viz"""
        self.assertEqual(exp, use.render())

    def test_metadata_merging(self):
        action = self.plugin.actions['identity_with_metadata']
        use = ArtifactAPIUsage()
        action.examples['identity_with_metadata_merging'](use)
        exp = """\
import qiime2.plugins.dummy_plugin.actions as dummy_plugin_actions

md3_md = md1_md.merge(md2_md)
out, = dummy_plugin_actions.identity_with_metadata(
    ints=ints,
    metadata=md3_md,
)"""
        self.assertEqual(exp, use.render())

    def test_metadata_column_from_helper(self):
        action = self.plugin.actions['identity_with_metadata_column']
        use = ArtifactAPIUsage()
        action.examples['identity_with_metadata_column_get_mdc'](use)
        exp = """\
import qiime2.plugins.dummy_plugin.actions as dummy_plugin_actions

mdc_mdc = md_md.get_column('a')
out, = dummy_plugin_actions.identity_with_metadata_column(
    ints=ints,
    metadata=mdc_mdc,
)"""
        self.assertEqual(exp, use.render())

    def test_optional_inputs(self):
        self.maxDiff = None
        action = self.plugin.actions['optional_artifacts_method']
        use = ArtifactAPIUsage()
        action.examples['optional_inputs'](use)
        exp = """\
import qiime2.plugins.dummy_plugin.actions as dummy_plugin_actions

output1, = dummy_plugin_actions.optional_artifacts_method(
    ints=ints,
    num1=1,
)
output2, = dummy_plugin_actions.optional_artifacts_method(
    ints=ints,
    num1=1,
    num2=2,
)
output3, = dummy_plugin_actions.optional_artifacts_method(
    ints=ints,
    num1=1,
    num2=None,
)
output4, = dummy_plugin_actions.optional_artifacts_method(
    ints=ints,
    optional1=output3,
    num1=3,
    num2=4,
)"""

        self.assertEqual(exp, use.render())

    def test_artifact_collection_dict_of_ints(self):
        action = self.plugin.actions['dict_of_ints']
        use = ArtifactAPIUsage()
        action.examples['collection_dict_of_ints'](use)
        exp = """\
import qiime2.plugins.dummy_plugin.actions as dummy_plugin_actions

out_artifact_collection, = dummy_plugin_actions.dict_of_ints(
    ints=ints_artifact_collection,
)"""

        self.assertEqual(exp, use.render())


if __name__ == '__main__':
    unittest.main()
