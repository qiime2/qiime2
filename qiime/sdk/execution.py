# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor


class SubprocessExecutor:

    _python_executable = "python"

    def __call__(self, workflow, input_artifact_filepaths,
                 parameter_references, output_artifact_filepaths):
        input_artifact_abs_filepaths = \
            {k: os.path.abspath(v)
             for k, v in input_artifact_filepaths.items()}
        output_artifact_abs_filepaths = \
            {k: os.path.abspath(v)
             for k, v in output_artifact_filepaths.items()}
        job = workflow.to_script(input_artifact_abs_filepaths,
                                 parameter_references,
                                 output_artifact_abs_filepaths)
        temp_dir = tempfile.mkdtemp()
        pool = ProcessPoolExecutor(max_workers=1)
        py_filename = os.path.join(temp_dir, 'job.py')
        with open(py_filename, 'w') as py_file:
            py_file.write(job.code)
        # TODO: handle subproccess exceptions
        future = pool.submit(subprocess.run,
                             [self._python_executable, py_filename],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        # TODO: handle callback exceptions
        # TODO: make sure that tempdir is cleaned up even if there is an
        # exception in pool.submit or the callback
        future.add_done_callback(lambda _: shutil.rmtree(temp_dir))

        return future
