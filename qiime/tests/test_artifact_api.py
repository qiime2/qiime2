# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import sys
import types
import importlib

import qiime.sdk
from qiime.core.testing.util import get_dummy_plugin


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
        self.assertIn('generated QIIME API', repr(module))

    def _check_plugin(self, module):
        self._check_spec(module, 'qiime.plugins.dummy_plugin', True)
        self._check_methods(module.methods)
        self._check_visualizers(module.visualizers)
        self.assertEqual(set(x for x in dir(module) if not x.startswith('_')),
                         {'visualizers', 'methods', 'actions'})

    def _check_methods(self, module):
        self._check_spec(module, 'qiime.plugins.dummy_plugin.methods', False)
        self.assertTrue(hasattr(module, 'concatenate_ints'))
        self.assertFalse(hasattr(module, 'most_common_viz'))

        self.assertIsInstance(module.concatenate_ints, qiime.sdk.Action)

    def _check_visualizers(self, module):
        self._check_spec(
            module, 'qiime.plugins.dummy_plugin.visualizers', False)
        self.assertTrue(hasattr(module, 'most_common_viz'))
        self.assertFalse(hasattr(module, 'concatenate_ints'))
        self.assertIsInstance(module.most_common_viz, qiime.sdk.Action)

    def test_import_root(self):
        import qiime.plugins.dummy_plugin
        self._check_plugin(qiime.plugins.dummy_plugin)

    def test_import_root_from(self):
        from qiime.plugins import dummy_plugin
        self._check_plugin(dummy_plugin)

    def test_import_methods(self):
        import qiime.plugins.dummy_plugin.methods
        self._check_methods(qiime.plugins.dummy_plugin.methods)

    def test_import_visualizers(self):
        import qiime.plugins.dummy_plugin.visualizers
        self._check_visualizers(qiime.plugins.dummy_plugin.visualizers)

    def test_import_methods_from(self):
        from qiime.plugins.dummy_plugin import methods
        self._check_methods(methods)

    def test_import_visualizers_from(self):
        from qiime.plugins.dummy_plugin import visualizers
        self._check_visualizers(visualizers)

    def test_import_non_plugin(self):
        with self.assertRaises(ImportError):
            import qiime.plugins.dummy_not_plugin  # noqa

    def test_import_non_action(self):
        with self.assertRaises(ImportError):
            import qiime.plugins.dummy_plugin.non_action  # noqa

    def test_import_side_module(self):
        # Certain implementations of __PATH__ can cause a module to load
        # siblings (__PATH__ = ['.'] for example)
        import qiime.metadata
        self.assertIsInstance(qiime.metadata, types.ModuleType)
        with self.assertRaises(ImportError):
            import qiime.plugins.metadata  # noqa

    def test_import_too_deep(self):
        with self.assertRaises(ImportError):
            import qiime.plugins.dummy_plugin.methods.too_deep  # noqa

    def test_import_non_module(self):
        with self.assertRaises(ImportError):
            import qiime.plugins.dummy_plugin.methods.concatenate_ints  # noqa

    def test_reload_fails(self):
        import qiime.plugins.dummy_plugin
        with self.assertRaises(ImportError):
            importlib.reload(qiime.plugins.dummy_plugin)


if __name__ == '__main__':
    unittest.main()
