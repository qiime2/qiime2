# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources
import tempfile
import unittest
import shutil
import pathlib

import qiime2

from qiime2.sdk import usage
from qiime2.plugin.model.base import FormatBase


# TODO Split out into more specific subclasses if necessary.
class TestPluginBase(unittest.TestCase):
    """Test harness for simplifying testing QIIME 2 plugins.

    ``TestPluginBase`` extends ``unittest.TestCase``, with a few extra helpers
    and assertions.

    Attributes
    ----------
    package : str
        The name of the plugin package to be tested.
    test_dir_prefix : str
        The prefix for the temporary testing dir created by the harness.

    """

    package = None
    test_dir_prefix = 'qiime2-plugin'

    def setUp(self):
        """Test runner setup hook.

        If overriding this hook in a test, call ``__super__`` to invoke this
        method in the overridden hook, otherwise the harness might not work
        as expected.

        """

        try:
            package = self.package.split('.')[0]
        except AttributeError:
            self.fail('Test class must have a package property.')

        # plugins are keyed by their names, so a search inside the plugin
        # object is required to match to the correct plugin
        plugin = None
        for name, plugin_ in qiime2.sdk.PluginManager().plugins.items():
            if plugin_.package == package:
                plugin = plugin_

        if plugin is not None:
            self.plugin = plugin
        else:
            self.fail('%s is not a registered QIIME 2 plugin.' % package)

        # TODO use qiime2 temp dir when ported to framework, and when the
        # configurable temp dir exists
        self.temp_dir = tempfile.TemporaryDirectory(
            prefix='%s-test-temp-' % self.test_dir_prefix)

    def tearDown(self):
        """Test runner teardown hook.

        If overriding this hook in a test, call ``__super__`` to invoke this
        method in the overridden hook, otherwise the harness might not work
        as expected.

        """

        self.temp_dir.cleanup()

    def get_data_path(self, filename):
        """Convenience method for getting a data asset while testing.

        Test data stored in the ``data/`` dir local to the running test
        can be accessed via this method.

        Parameters
        ----------
        filename : str
            The name of the file to look up.

        Returns
        -------
        filepath : str
            The materialized filepath to the requested test data.

        """

        return pkg_resources.resource_filename(self.package,
                                               'data/%s' % filename)

    def get_transformer(self, from_type, to_type):
        """Convenience method for getting a registered transformer.

        Parameters
        ----------
        from_type : A View Type
            The :term:`View` type of the source data.
        to_type : A View Type
            The :term:`View` type to transform to.

        Returns
        -------
        transformer : A Transformer Function
            The registered tranformer from ``from_type`` to ``to_type``.

        """

        try:
            transformer_record = self.plugin.transformers[from_type, to_type]
        except KeyError:
            self.fail(
                "Could not find registered transformer from %r to %r." %
                (from_type, to_type))

        return transformer_record.transformer

    def assertRegisteredSemanticType(self, semantic_type):
        """Test assertion for ensuring a plugin's semantic type is registered.

        Fails if the semantic type requested is not found in the Plugin
        Manager.

        Parameters
        ----------
        semantic_type : A Semantic Type
            The :term:`Semantic Type` to test the presence of.

        """

        try:
            semantic_type_record = self.plugin.types[semantic_type.name]
        except KeyError:
            self.fail(
                "Semantic type %r is not registered on the plugin." %
                semantic_type)

        obs_semantic_type = semantic_type_record.semantic_type

        self.assertEqual(obs_semantic_type, semantic_type)

    def assertSemanticTypeRegisteredToFormat(self, semantic_type, exp_format):
        """Test assertion for ensuring a semantic type is registered to a
           format.

        Fails if the semantic type requested is not registered to the format
        specified with ``exp_format``. Also fails if the semantic type isn't
        registered to **any** format.

        Parameters
        ----------
        semantic_type : A Semantic Type
            The :term:`Semantic Type` to check for.
        exp_format : A Format
            The :term:`Format` to check that the Semantic Type is registed on.

        """

        obs_format = None
        for type_format_record in self.plugin.type_formats:
            if type_format_record.type_expression == semantic_type:
                obs_format = type_format_record.format
                break

        self.assertIsNotNone(
            obs_format,
            "Semantic type %r is not registered to a format." % semantic_type)

        self.assertEqual(
            obs_format, exp_format,
            "Expected semantic type %r to be registered to format %r, not %r."
            % (semantic_type, exp_format, obs_format))

    def transform_format(self, source_format, target, filename=None,
                         filenames=None):
        """Helper utility for loading data and transforming it.

        Combines several other utilities in this class, will load files from
        ``data/``, as ``source_format``, then transform to the ``target`` view.

        Parameters
        ----------
        source_format : A Format
            The :term:`Format` to load the data as.
        target : A View Type
            The :term:`View Type <View>` to transform the data to.
        filename : str
            The name of the file to load from ``data``. Use this for formats
            that use a single file in their format definition. Mutually
            exclusive with the ``filenames`` parameter.
        filenames : list[str]
            The names of the files to load from ``data``. Use this for formats
            that use multiple files in their format definition. Mutually
            exclusive with the ``filename`` parameter.

        Returns
        -------
        input : A Format
            The data loaded from ``data`` as the specified ``source_format``.
        obs : A View Type
            The loaded data, transformed to the specified ``target`` view type.

        """

        # Guard any non-QIIME2 Format sources from being tested
        if not issubclass(source_format, FormatBase):
            raise ValueError("`source_format` must be a subclass of "
                             "FormatBase.")

        # Guard against invalid filename(s) usage
        if filename is not None and filenames is not None:
            raise ValueError("Cannot use both `filename` and `filenames` at "
                             "the same time.")

        # Handle format initialization
        source_path = None
        if filename:
            source_path = self.get_data_path(filename)
        elif filenames:
            source_path = self.temp_dir.name
            for filename in filenames:
                filepath = self.get_data_path(filename)
                shutil.copy(filepath, source_path)
        input = source_format(source_path, mode='r')

        transformer = self.get_transformer(source_format, target)
        obs = transformer(input)

        if issubclass(target, FormatBase):
            self.assertIsInstance(obs, (type(pathlib.Path()), str, target))
        else:
            self.assertIsInstance(obs, target)

        return input, obs

    def execute_examples(self):
        if self.plugin is not None:
            for _, action in self.plugin.actions.items():
                for name, example_f in action.examples.items():
                    with self.subTest(example=name):
                        use = usage.ExecutionUsage()
                        example_f(use)
