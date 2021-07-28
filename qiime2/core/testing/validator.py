from .type import Kennel, Dog, Cat
from .plugin import dummy_plugin

@dummy_plugin.register_validator(Kennel[Dog])
def validator_test_null():
    print('Doesn\'t do anything')

@dummy_plugin.register_validator(Kennel[Dog | Cat])
def test_subset_or():
    pass

@dummy_plugin.register_validator(Kennel[Dog])
def validator_test_null2():
    print('does know everything')
