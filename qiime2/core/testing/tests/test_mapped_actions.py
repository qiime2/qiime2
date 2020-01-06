# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from qiime2 import Artifact
from qiime2.core.testing.util import get_dummy_plugin


class ActionTester(unittest.TestCase):
    ACTION = 'N/A'

    def setUp(self):
        plugin = get_dummy_plugin()
        self.action = plugin.actions[self.ACTION]

    def run_action(self, **inputs):
        results = self.action(**inputs)

        future = self.action.asynchronous(**inputs)
        async_results = future.result()

        for a, b in zip(async_results, results):
            self.assertEqual(a.type, b.type)

        return results


class TestConstrainedInputVisualization(ActionTester):
    ACTION = 'constrained_input_visualization'

    def test_match_foo(self):
        a = Artifact.import_data('Foo', "element 1", view_type=str)
        b = Artifact.import_data('Foo', "element 2", view_type=str)

        viz, = self.run_action(a=a, b=b)

        contents = (viz._archiver.data_dir / 'index.html').read_text()
        self.assertIn('element 1', contents)
        self.assertIn('element 2', contents)

    def test_match_nested(self):
        a = Artifact.import_data('C1[Baz]', "element 1", view_type=str)
        b = Artifact.import_data('C1[Baz]', "element 2", view_type=str)

        viz, = self.run_action(a=a, b=b)

        contents = (viz._archiver.data_dir / 'index.html').read_text()
        self.assertIn('element 1', contents)
        self.assertIn('element 2', contents)

    def test_mismatch_foo_bar(self):
        a = Artifact.import_data('Foo', "element 1", view_type=str)
        b = Artifact.import_data('Bar', "element 2", view_type=str)

        with self.assertRaisesRegex(ValueError, 'No solution.*Foo'):
            viz, = self.run_action(a=a, b=b)

    def test_mismatch_nested(self):
        a = Artifact.import_data('C1[Foo]', "element 1", view_type=str)
        b = Artifact.import_data('Foo', "element 2", view_type=str)

        with self.assertRaisesRegex(ValueError, 'No solution.*C1'):
            viz, = self.run_action(a=a, b=b)


class TestCombinatoricallyMappedMethod(ActionTester):
    ACTION = 'combinatorically_mapped_method'

    def test_match_foo(self):
        a = Artifact.import_data('C1[Foo]', 'element 1', view_type=str)
        b = Artifact.import_data('C3[Foo, Foo, Foo]',
                                 'element 2', view_type=str)

        x, y = self.run_action(a=a, b=b)

        self.assertEqual(repr(x.type), 'C2[Bar, Bar]')
        self.assertEqual(repr(y.type), 'Foo')

    def test_match_bar_foo(self):
        a = Artifact.import_data('C1[Bar]', 'element 1', view_type=str)
        b = Artifact.import_data('C3[Foo, Foo, Foo]',
                                 'element 2', view_type=str)

        x, y = self.run_action(a=a, b=b)

        self.assertEqual(repr(x.type), 'C2[Baz, Baz]')
        self.assertEqual(repr(y.type), 'Foo')

    def test_match_baz_misc(self):
        a = Artifact.import_data('C1[Baz]', 'element 1', view_type=str)
        b = Artifact.import_data('C3[Foo, Bar, Baz]',
                                 'element 2', view_type=str)

        x, y = self.run_action(a=a, b=b)

        self.assertEqual(repr(x.type), 'C2[Foo, Foo]')
        self.assertEqual(repr(y.type), 'Baz')

    def test_mismatch(self):
        a = Artifact.import_data('Bar', 'element 1', view_type=str)
        b = Artifact.import_data('C3[Foo, Foo, Foo]',
                                 'element 2', view_type=str)

        with self.assertRaises(TypeError):
            self.run_action(a=a, b=b)


class TestDoubleBoundVariableMethod(ActionTester):
    ACTION = 'double_bound_variable_method'

    def test_predicate_on_second(self):
        a = Artifact.import_data('Bar', 'element 1', view_type=str)
        b = Artifact.import_data('Bar % Properties("A")',
                                 'element 2', view_type=str)
        extra = Artifact.import_data('Foo', 'always foo', view_type=str)

        x, = self.run_action(a=a, b=b, extra=extra)

        self.assertEqual(repr(x.type), 'Baz')

    def test_mismatch(self):
        a = Artifact.import_data('Foo', 'element 1', view_type=str)
        b = Artifact.import_data('Bar', 'element 2', view_type=str)
        extra = Artifact.import_data('Foo', 'always foo', view_type=str)

        with self.assertRaisesRegex(ValueError, 'match.*same output'):
            self.run_action(a=a, b=b, extra=extra)


class TestBoolFlagSwapsOutputMethod(ActionTester):
    ACTION = 'bool_flag_swaps_output_method'

    def test_true(self):
        a = Artifact.import_data('Bar', 'element', view_type=str)

        x, = self.run_action(a=a, b=True)

        self.assertEqual(repr(x.type), 'C1[Foo]')

    def test_false(self):
        a = Artifact.import_data('Bar', 'element', view_type=str)

        x, = self.run_action(a=a, b=False)

        self.assertEqual(repr(x.type), 'Foo')


class TestPredicatesPreservedMethod(ActionTester):
    ACTION = 'predicates_preserved_method'

    def test_simple(self):
        a = Artifact.import_data("Foo % Properties('A')",
                                 'element 1', view_type=str)

        x, = self.run_action(a=a)

        self.assertEqual(repr(x.type), "Foo % Properties('A')")

    def test_mismatch(self):
        a = Artifact.import_data("Foo % Properties('X')",
                                 'element 1', view_type=str)

        with self.assertRaises(TypeError):
            self.run_action(a=a)

    def test_combinations_preserved(self):
        a = Artifact.import_data("Foo % Properties('A', 'B')",
                                 'element 1', view_type=str)

        x, = self.run_action(a=a)

        self.assertEqual(repr(x.type), "Foo % Properties('A', 'B')")

    def test_extra_dropped(self):
        a = Artifact.import_data("Foo % Properties('Extra', 'A', 'B')",
                                 'element 1', view_type=str)

        x, = self.run_action(a=a)

        self.assertEqual(repr(x.type), "Foo % Properties('A', 'B')")


if __name__ == '__main__':
    unittest.main()
