# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import parsl
import psutil
import appdirs
import tomlkit
import threading
import importlib
import pkg_resources

PARALLEL_CONFIG = threading.local()
PARALLEL_CONFIG.dfk = None
PARALLEL_CONFIG.parallel_config = None
PARALLEL_CONFIG.action_executor_mapping = {}

# We write a default config to a location in the conda env if there is an
# active conda env. If there is not an active conda env (most likely because we
# are using Docker) then the path we want to write the default to will not
# exist, so we will not write a default, we will just load it from memory
CONDA_PREFIX = os.environ.get('CONDA_PREFIX', '')
VENDORED_FP = os.path.join(CONDA_PREFIX, 'etc', 'qiime2_config.toml')

VENDORED_CONFIG = {
    'parsl': {
        'strategy': 'None',
        'executors': [
            {'class': 'ThreadPoolExecutor', 'label': 'default',
                'max_threads': max(psutil.cpu_count() - 1, 1)},
            {'class': 'HighThroughputExecutor', 'label': 'htex',
                'max_workers': max(psutil.cpu_count() - 1, 1),
                'provider': {'class': 'LocalProvider'}}
            ]
        }
    }

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


def _setup_parallel(config_fp=None):
    """Sets the parsl config and action executor mapping from a file at a given
    path or looks through several default paths if no path is provided and
    loads a vendored config as a last resort
    """
    parallel_config = PARALLEL_CONFIG.parallel_config
    mapping = PARALLEL_CONFIG.action_executor_mapping

    # If we don't have a filepath or a currently existing config then get the
    # path to the vendored one. We do not want to get the vendored path if they
    # have a pre-existing config because we do not want to overwrite an
    # existing config with the vendored one
    if config_fp is None and PARALLEL_CONFIG.parallel_config is None:
        config_fp = _get_vendored_config()

    if config_fp is not None:
        parallel_config, mapping = get_config(config_fp)

        # If we don't have a config now, they gave a file that only contained a
        # mapping, so we want to load a default config assuming they do not
        # already have a loaded config
        if parallel_config is None and PARALLEL_CONFIG.parallel_config is None:
            config_fp = _get_vendored_config()
            parallel_config, = get_config(config_fp)

    _finalize_setup(parallel_config, mapping)


def _finalize_setup(parallel_config, mapping):
    """Finish loading the config and setting up our threadlocal
    """
    PARALLEL_CONFIG.dfk = parsl.load(parallel_config)
    PARALLEL_CONFIG.parallel_config = parallel_config
    if mapping != {}:
        PARALLEL_CONFIG.action_executor_mapping = mapping


def _cleanup_parallel():
    """Ask parsl to cleanup and then remove the currently active dfk
    """
    PARALLEL_CONFIG.dfk.cleanup()
    parsl.clear()


def get_config(fp):
    """Takes a config filepath and determines if the file exists and if so if
    it contains parsl config info.
    """
    with open(fp, 'r') as fh:
        config_dict = tomlkit.load(fh)

    return get_config_from_dict(config_dict)


def get_config_from_dict(config_dict):
    parallel_config_dict = config_dict.get('parsl')
    mapping = parallel_config_dict.pop('executor_mapping', {})

    processed_parallel_config_dict = _process_config(parallel_config_dict)

    # They could have given us a file that contained a mapping but no config.
    # This is technically valid.
    if processed_parallel_config_dict != {}:
        parallel_config = parsl.Config(**processed_parallel_config_dict)
    else:
        parallel_config = None

    return parallel_config, mapping


def _get_vendored_config():
    if 'QIIMETEST' in os.environ:
        return pkg_resources.resource_filename(
            'qiime2.sdk.tests', 'data/test_config.toml')

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
        elif os.path.exists(fp_ := os.path.join(
                appdirs.site_config_dir('qiime2'), 'qiime2_config.toml')):
            config_fp = fp_
        # NOTE: These next two are dependent on us being in a conda environment
        # 4. Check in conda env
        # ~/miniconda3/env/{env_name}/conf
        elif CONDA_PREFIX != '' and os.path.exists(fp_ := VENDORED_FP):
            config_fp = fp_
        # 5. Write the vendored config to the vendored location and use
        # that
        elif CONDA_PREFIX != '':
            with open(VENDORED_FP, 'w') as fh:
                tomlkit.dump(VENDORED_CONFIG, fh)
            config_fp = VENDORED_FP

    return config_fp


def _process_config(config_dict):
    """Takes a path to a toml file describing a parsl.Config object and parses
    it into a dictionary of kwargs that can be used to instantiate a
    parsl.Config object.
    """
    config_kwargs = {}

    for key, value in config_dict.items():
        # Parsl likes the string 'none' as opposed to None or 'None'
        if isinstance(value, str) and value.lower() == 'none':
            config_kwargs[key] = value.lower()
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
        # Get the module our class is from
        module = importlib.import_module(module_paths[key])

        _type = value['class']
        if _type == '_TEST_EXECUTOR_':
            # Only used for tests
            cls = _TEST_EXECUTOR_
        else:
            # Get the class we need to instantiate
            cls = getattr(module, value['class'])

        # Get the kwargs we need to pass to the class constructor
        kwargs = {}
        for k, v in value.items():
            # We already handled this key
            if k != 'class':
                kwargs[k] = _process_key(k, v)

        # Instantiate the class
        return cls(**kwargs)
    # Our key points to primitive data
    else:
        return value


class ParallelConfig():
    def __init__(self, parallel_config=None, action_executor_mapping={}):
        """Tell QIIME 2 how to parsl from the Python API

        action_executor_mapping: maps actions to executors. All unmapped
        actions will be run on the default executor. We check which executor a
        given action is supposed to use when we get ready to run the action, so
        errors will only occur if an action that is being run in a given
        QIIME 2 invocation has been mapped to an executor that does not exist

        parallel_config: Specifies which executors should be created and how
        they should be created. If this is None, it will use the default
        config.
        """
        self.parallel_config = parallel_config
        self.action_executor_mapping = action_executor_mapping

    def __enter__(self):
        """Set this to be our Parsl config on the current thread local
        """
        if PARALLEL_CONFIG.parallel_config is not None:
            raise ValueError('ParallelConfig already loaded, cannot nest '
                             'ParallelConfigs')

        PARALLEL_CONFIG.parallel_config = self.parallel_config
        PARALLEL_CONFIG.action_executor_mapping = self.action_executor_mapping

        _setup_parallel()

    def __exit__(self, *args):
        """Set our Parsl config back to whatever it was before this one
        """
        _cleanup_parallel()

        PARALLEL_CONFIG.dfk = None
        PARALLEL_CONFIG.parallel_config = None
        PARALLEL_CONFIG.action_executor_mapping = {}


def _check_env(cls):
    if 'QIIMETEST' not in os.environ:
        raise ValueError(
            f"Do not instantiate the class '{cls}' when not testing")


class _TEST_EXECUTOR_(parsl.executors.threads.ThreadPoolExecutor):
    """We needed multiple kinds of executor to ensure we were mapping things
    correctly, but the HighThroughputExecutor was leaking sockets, so we avoid
    creating those during the tests because so many sockets were being opened
    that we were getting "Too many open files" errors, so this gets used as the
    second executor type."""

    def __init__(self, *args, **kwargs):
        _check_env(self.__class__)
        super(_TEST_EXECUTOR_, self).__init__(*args, **kwargs)


# Used to test config loading behavior when outside of a conda environment
class _MASK_CONDA_ENV_():
    def __init__(self):
        _check_env(self.__class__)

    def __enter__(self):
        global CONDA_PREFIX, VENDORED_FP

        self.old_prefix = CONDA_PREFIX
        self.old_fp = VENDORED_FP

        CONDA_PREFIX = ''
        VENDORED_FP = None

    def __exit__(self, *args):
        global CONDA_PREFIX, VENDORED_FP

        CONDA_PREFIX = self.old_prefix
        VENDORED_FP = self.old_fp
