import unittest
import unittest

from qiime2 import Artifact

class TestExample(unittest.TestCase):
  def test_example(self):
    # future = dummy.methods['params_only_method'].parsl('a', 1)
    # future.result()
    import qiime2.sdk

    pm = qiime2.sdk.PluginManager()
    self.dummy = pm.get_plugin(id='dummy_plugin')
    a = Artifact.import_data('Foo', "element 1", view_type=str)
    b = Artifact.import_data('Foo', "element 2", view_type=str)

    print(self.dummy)
    future = self.dummy.actions['constrained_input_visualization'].parsl(a=a, b=b)
    future.result()

  test_example2 = test_example
