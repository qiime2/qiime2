from .type import Kennel, Dog
from .plugin import dummy_plugin

@dummy_plugin.register_validator(Kennel[Dog])
def validator_test_null():
    print('Doesn\'t do anything')
