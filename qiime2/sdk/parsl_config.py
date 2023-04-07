# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import toml
import psutil
import appdirs
import threading
import importlib

# Do a parsl?
from parsl.config import Config
from parsl.providers import LocalProvider
from parsl.executors.threads import ThreadPoolExecutor
from parsl.executors import HighThroughputExecutor


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


def get_config(fp):
    """Takes a config filepath and determines if the file exists and if so if
    it contains parsl config info.
    """
    if fp is None:
        return None

    config_dict = toml.load(fp)
    return config_dict.pop('parsl')


def process_key(key, value):
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
            kwargs[k] = process_key(k, v)
        return class_(**kwargs)
    # Our key points to primitive data
    else:
        return value


def process_config(config_dict):
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
                config_kwargs[key].append(process_key(key, item))
        # We have a single value
        else:
            config_kwargs[key] = process_key(key, config_dict[key])

    return config_kwargs


def get_parsl_config():
    # If a config was already withed in by the user do nothing
    if PARSL_CONFIG.parsl_config is not None:
        return PARSL_CONFIG.parsl_config

    conda_env = os.environ.get('CONDA_DEFAULT_ENV')

    # Try to load custom config (exact order maybe subject to change)
    # Check envvar
    config_fp = os.environ.get('QIIME_CONF')

    if config_fp is None:
        # Check in user writable location
        #  appdirs.user_config_dir(appname='qiime2', author='...')
        if os.path.exists(fp_ := os.path.join(
                appdirs.user_config_dir('qiime2'), 'qiime_conf.toml')):
            config_fp = fp_
        # Check for conf in conda env
        # ~/miniconda3/env/{env_name}/conf
        # TODO: We don't enforce use of miniconda3 do we? Makes this
        # potentially problematic
        elif os.path.exists(fp_ := os.path.join(
                os.path.expanduser('~'), 'miniconda3', 'env', conda_env,
                'conf', 'qiime_conf.toml')):
            config_fp = fp_
        # If nothing check in admin writable location
        # /etc/
        # site_config_dir
        # appdirs.site_config_dir(appname='qiime2, author='...')
        elif config_fp is None and os.path.exists(fp_ := os.path.join(
                appdirs.site_config_dir('qiime2'), 'qiime_conf.toml')):
            config_fp = fp_

    # Check if we have a config file containing parsl config info
    config_dict = get_config(config_fp)

    # Load vendored default
    # If no custom config do this. This will probably end up in a vendored file
    if config_dict is None:
        config = Config(
            executors=[
                HighThroughputExecutor(
                    label='default',
                    max_workers=6,
                    provider=LocalProvider()
                ),
                ThreadPoolExecutor(
                    max_threads=max(psutil.cpu_count() - 1, 1),
                    label='tpool'
                )
            ],
            #  AdHoc Clusters should not be setup with scaling strategy.
            strategy='none',
        )
    else:
        config = Config(**process_config(config_dict))

    PARSL_CONFIG.parsl_config = config
    return PARSL_CONFIG.parsl_config


class ParallelConfig():
    def __init__(self, parsl_config=None, action_executor_mapping={}):
        """Tell QIIME 2 how to parsl

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
