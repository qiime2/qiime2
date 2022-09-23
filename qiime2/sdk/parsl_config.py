# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import psutil
import threading
import os

# Do a parsl?
from parsl.config import Config
from parsl.providers import LocalProvider
from parsl.executors.threads import ThreadPoolExecutor
from parsl.executors import HighThroughputExecutor


PARSL_CONFIG = threading.local()
PARSL_CONFIG.parsl_config = None
PARSL_CONFIG.action_executor_mapping = {}


def get_parsl_config():
    # If a config was already withed in by the user do nothing
    if PARSL_CONFIG.parsl_config is not None:
        pass
    # Try to load custom config (exact order maybe subject to change)
    # Check envvar
    elif os.environ.get('QIIME_PARSL_CONF'):
        # Load file pointed at
        pass
    # Check in user writable location
    # NOTE: These are not actually kwargs I made them kwargs in the comment to
    # remind me of the order
    # appdirs.user_config_dir(appname='qiime2', author='...')
    elif 0:  # some check
        pass
    # Check for conf in conda env
    # ~/miniconda3/env/{env_name}/conf
    elif 0:
        pass
    # If nothing check in admin writable location on hpc/user writable local
    # /etc/
    # site_config_dir
    # appdirs.site_config_dir(appname='qiime2, author='...')
    elif 0:  # another check
        pass
    # Load vendored default
    # If no custom config do this. This will probably end up in a vendored file
    else:
        PARSL_CONFIG.parsl_config = Config(
            executors=[
                ThreadPoolExecutor(
                    max_threads=max(psutil.cpu_count() - 1, 1),
                    label='htex'
                ),
                HighThroughputExecutor(
                    label='default',
                    max_workers=6,

                    provider=LocalProvider()
                )
            ],
            #  AdHoc Clusters should not be setup with scaling strategy.
            strategy=None,
        )

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
