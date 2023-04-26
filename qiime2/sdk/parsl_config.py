# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import psutil
import appdirs
import tomlkit
import threading
import importlib

from parsl.config import Config


PARSL_CONFIG = threading.local()
PARSL_CONFIG.parsl_config = None
PARSL_CONFIG.action_executor_mapping = {}

# Directs keys in the config whose values need to be objects to the module that
# contains the class they need to instantiate
module_paths = {
    'channel': 'parsl.channels',
    'channels': 'parsl.channels',
    'data_provider': 'parsl.data_provider',
    'data_providers': 'parsl.data_provider',
    'executor': 'parsl.executors',
    'executors': 'parsl.executors',
    'launcher': 'parsl.launcher',
    'launchers': 'parsl.launcher',
    'monitoring': 'parsl.monitoring',
    'provider': 'parsl.providers',
    'providers': 'parsl.providers'
}

VENDORED_FP = os.path.join(
    os.environ.get('CONDA_PREFIX'), 'etc', 'qiime2_config.toml')
VENDORED_CONFIG = {
    'parsl': {
        'strategy': 'None',
        'executors': [
            {'class': 'ThreadPoolExecutor', 'label': 'default',
             'max_threads': max(psutil.cpu_count() - 1, 1)},
            {'class': 'HighThroughputExecutor', 'label': 'htex',
             'max_workers': max(psutil.cpu_count() - 1, 1),
             'provider': {'class': 'LocalProvider'}}
            ],
        'executor_mapping': {'list_of_ints': 'htex'}
        }
    }


def setup_parsl(config_fp=None):
    """Sets the parsl config and action executor mapping from a file at a given
    path or looks through several default paths if no path is provided and
    loads a vendored config as a last resort
    """
    if config_fp is None:
        # If a config was already explicitly set and this function was called
        # again without an explicit filepath, do not replace the manually set
        # config with an internal one
        if PARSL_CONFIG.parsl_config is not None:
            return

        # Try to load custom config (exact order maybe subject to change)

        # 1. Check envvar
        config_fp = os.environ.get('QIIME2_CONFIG')

        if config_fp is None:
            # 2. Check in user writable location
            # appdirs.user_config_dir(appname='qiime2', author='...')
            if os.path.exists(fp_ := os.path.join(
                    appdirs.user_config_dir('qiime2'), 'qiime2_config.toml')):
                config_fp = fp_
            # 3. Check in admin writable location
            # /etc/
            # site_config_dir
            # appdirs.site_config_dir(appname='qiime2, author='...')
            elif config_fp is None and os.path.exists(fp_ := os.path.join(
                    appdirs.site_config_dir('qiime2'), 'qiime2_config.toml')):
                config_fp = fp_
            # 4. Check in conda env
            # ~/miniconda3/env/{env_name}/conf
            elif os.path.exists(fp_ := VENDORED_FP):
                config_fp = fp_
            # 5. Write the vendored config to the vendored location and use
            # that
            else:
                with open(VENDORED_FP, 'w') as fh:
                    tomlkit.dump(VENDORED_CONFIG, fh)
                config_fp = VENDORED_FP

    config_dict = _get_config(config_fp)
    mapping = _get_mapping(config_dict)

    processed_config = _process_config(config_dict)
    config = Config(**processed_config)

    PARSL_CONFIG.parsl_config = config
    PARSL_CONFIG.action_executor_mapping = mapping


def _get_config(fp):
    """Takes a config filepath and determines if the file exists and if so if
    it contains parsl config info.
    """
    with open(fp, 'r') as fh:
        config_dict = tomlkit.load(fh)

    return config_dict.get('parsl')


def _get_mapping(config_dict):
    """Takes a config filepath and determines if the file exists and if so if
    it contains action executor mapping info.
    """
    return config_dict.pop('executor_mapping', {})


def _process_config(config_dict):
    """Takes a path to a toml file describing a parsl.Config object and parses
    it into a dictionary of kwargs that can be used to instantiate a
    parsl.Config object.
    """
    config_kwargs = {}

    for key, value in config_dict.items():
        # Our value is a None
        if isinstance(value, str) and value.lower() == 'none':
            config_kwargs[key] = None
        # We have a list of values
        elif isinstance(value, list):
            config_kwargs[key] = []
            for item in value:
                config_kwargs[key].append(_process_key(key, item))
        # We have a single value
        else:
            config_kwargs[key] = _process_key(key, config_dict[key])

    return config_kwargs


def _process_key(key, value):
    """Takes a key given in the parsl config file and turns its value into the
    correct data type or class instance to be used in instantiating a
    parsl.Config object.
    """
    # Our key needs to point to some object.
    if key in module_paths:
        module = importlib.import_module(module_paths[key])
        class_ = getattr(module, value.pop('class'))
        kwargs = {}
        for k, v in value.items():
            kwargs[k] = _process_key(k, v)
        return class_(**kwargs)
    # Our key points to primitive data
    else:
        return value


class ParallelConfig():
    def __init__(self, parsl_config=None, action_executor_mapping={}):
        """Tell QIIME 2 how to parsl from the Python API

        action_executor_mapping: maps actions to executors. All unmapped
        actions will be run on the default executor. We check which executor a
        given action is supposed to use when we get ready to run the action, so
        errors will only occur if an action that is being run in a given
        QIIME 2 invocation has been mapped to an executor that does not exist

        parsl_config: Specifies which executors should be created and how they
        should be created
        """
        self.parsl_config = parsl_config
        self.action_executor_mapping = action_executor_mapping

    def __enter__(self):
        """Set this to be our Parsl config on the current thread local
        """
        self.backup_config = PARSL_CONFIG.parsl_config
        PARSL_CONFIG.parsl_config = self.parsl_config

        self.backup_map = PARSL_CONFIG.action_executor_mapping
        PARSL_CONFIG.action_executor_mapping = self.action_executor_mapping

    def __exit__(self, *args):
        """Set our Parsl config back to whatever it was before this one
        """
        PARSL_CONFIG.parsl_config = self.backup_config
        PARSL_CONFIG.action_executor_mapping = self.backup_map
