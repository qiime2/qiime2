# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import psutil
import threading
import os

# Do a parsl?
import parsl
# join_app or @python_app(join=True)
from parsl.app.app import python_app
from parsl.config import Config
from parsl.providers import AdHocProvider
from parsl.executors.threads import ThreadPoolExecutor
from parsl.executors import HighThroughputExecutor

import toml

import appdirs


LOCAL_CONFIG = threading.local()
LOCAL_CONFIG.config = None

def get_config():
    # Try to load custom config (exact order maybe subject to change)
    # Check envvar
    if os.environ.get('QIIME_PARSL_CONF'):
        # Load file pointed at
        pass
    # Check in user writable location
    # NOTE: These are not actually kwargs I made them kwargs in the comment to
    # remind me of the order
    # appdirs.user_config_dir(appname='qiime2', author='...')
    elif 0: # some check
        pass
    # Check for conf in conda env
    # ~/miniconda3/env/{env_name}/conf
    elif 0:
        pass
    # If nothing check in admin writable location on hpc/user writable local
    # /etc/
    # site_config_dir
    # appdirs.site_config_dir(appname='qiime2, authoe='...')
    elif 0: # another check
        pass
    # Load vendored default
    else:
        pass


    # If not custom config do this
    # DO NOT USE SAME CONFIG MORE THAN ONCE
    if True:#LOCAL_CONFIG.config is None:
        # LOCAL_CONFIG.config = Config(
        #     executors=[
        #         ThreadPoolExecutor(
        #             max_threads=max(psutil.cpu_count() - 1, 1),
        #             label='default_config'
        #         )
        #     ]
        # )
        # How to add executors?
        # How to specify whiche executor to use (perhaps the labels?)
        LOCAL_CONFIG.config = Config(
        executors=[
            HighThroughputExecutor(
                label='htex',
                max_workers=6,
                worker_logdir_root=os.getcwd(),
                provider=AdHocProvider(
                    # Command to be run before starting a worker, such as:
                    # 'module load Anaconda; source activate parsl_env'.
                    worker_init='',
                    channels=[parsl.channels.LocalChannel(script_dir=os.getcwd())]
                )
            )
        ],
        #  AdHoc Clusters should not be setup with scaling strategy.
        strategy=None,
    )

    return LOCAL_CONFIG.config
