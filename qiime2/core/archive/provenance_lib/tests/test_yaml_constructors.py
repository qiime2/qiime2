import unittest
import warnings
import yaml

from .._yaml_constructors import MetadataInfo


class YamlConstructorTests(unittest.TestCase):
    """
    YAML Constructors are used to handle the custom YAML tags defined by the
    framework.

    """
    def test_unknown_tag(self):
        """
        Makes explicit the current handling of unimplemented custom tags. In
        future, we may want to deal with these more graciously (e.g. warn), but
        for now we're going to fail fast
        """
        tag = r"!foo 'this is not an implemented tag'"
        with self.assertRaisesRegex(yaml.constructor.ConstructorError,
                                    'could not determine a constructor.*!foo'):
            yaml.safe_load(tag)

    def test_citation_key_constructor(self):
        tag = r"!cite 'framework|qiime2:2020.6.0.dev0|0'"
        actual = yaml.safe_load(tag)
        self.assertEqual(actual, 'framework|qiime2:2020.6.0.dev0|0')

    def test_color_primitive_constructor(self):
        tag = r"!color '#57f289'"
        actual = yaml.safe_load(tag)
        self.assertEqual(actual, '#57f289')

    def test_forward_ref_action_plugin_ref(self):
        tag = r"plugin: !ref 'environment:plugins:diversity'"
        actual = yaml.safe_load(tag)
        self.assertEqual(actual, {'plugin': 'diversity'})

    def test_forward_ref_generic_ref(self):
        tag = r"plugin: !ref 'environment:framework:version'"
        actual = yaml.safe_load(tag)
        exp = {'plugin': ['environment', 'framework', 'version']}
        self.assertEqual(exp, actual)

    def test_metadata_path_constructor(self):
        tag = r"!metadata 'metadata.tsv'"
        actual = yaml.safe_load(tag)
        self.assertEqual(actual, MetadataInfo([], 'metadata.tsv'))

    def test_metadata_path_constructor_one_Artifact_as_md(self):
        tag = r"!metadata '415409a4-stuff-e3eaba5301b4:feature_metadata.tsv'"
        actual = yaml.safe_load(tag)
        self.assertEqual(
            actual,
            MetadataInfo(['415409a4-stuff-e3eaba5301b4'],
                         'feature_metadata.tsv'))

    def test_metadata_path_constructor_many_Artifacts_as_md(self):
        tag = (r"!metadata '415409a4-stuff-e3eaba5301b4,"
               r"12345-other-stuff-67890"
               r":feature_metadata.tsv'")
        actual = yaml.safe_load(tag)
        self.assertEqual(
            actual,
            MetadataInfo(['415409a4-stuff-e3eaba5301b4',
                          '12345-other-stuff-67890'],
                         'feature_metadata.tsv'))

    def test_no_provenance_constructor(self):
        tag = "!no-provenance '34b07e56-27a5-4f03-ae57-ff427b50aaa1'"
        with self.assertWarnsRegex(UserWarning,
                                   'Artifact 34b07e.*prior to provenance'):
            actual = yaml.safe_load(tag)
            self.assertEqual(actual, '34b07e56-27a5-4f03-ae57-ff427b50aaa1')

    def test_no_provenance_multiple_warnings_fire(self):
        tag_list = """
        - !no-provenance '34b07e56-27a5-4f03-ae57-ff427b50aaa1'
        - !no-provenance 'gerbil'
        """
        with warnings.catch_warnings(record=True) as w:
            # Just in case something else has modified the filter state
            warnings.simplefilter("default")
            yaml.safe_load(tag_list)
            # There should be exactly two warnings
            self.assertEqual(len(w), 2)

            # The first should be a Userwarning containing these strings
            self.assertEqual(UserWarning, w[0].category)
            self.assertIn('Artifact 34b07e', str(w[0].message))
            self.assertIn('prior to provenance', str(w[0].message))

            # And the second should look similar
            self.assertEqual(UserWarning, w[1].category)
            self.assertIn('gerbil', str(w[1].message))
            self.assertIn('prior to provenance', str(w[0].message))

    def test_set_ref(self):
        flow_tag = r"!set ['foo', 'bar', 'baz']"
        flow = yaml.safe_load(flow_tag)
        self.assertEqual(flow, {'foo', 'bar', 'baz'})

        # NOTE: we don't expect duplicate values here (because dumped values
        # were a set), but it doesn't hurt to test the behavior
        block_tag = '!set\n- spam\n- egg\n- spam\n'
        block = yaml.safe_load(block_tag)
        self.assertEqual(block, {'spam', 'egg'})
