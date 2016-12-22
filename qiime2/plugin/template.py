# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import cookiecutter.main
import cookiecutter.exceptions


def plugin_init(output_dir='.'):
    """Initialize a plugin package from a template.

    Parameters
    ----------
    output_dir : str, optional
        Output directory in which to initialize plugin package subdirectory.

    Returns
    -------
    str
        Absolute path to plugin package directory created in `output_dir`.

    Raises
    ------
    OutputDirExistsException
        If the plugin package directory name already exists under `output_dir`.

    """
    try:
        return cookiecutter.main.cookiecutter(
            "https://github.com/qiime2/cookiecutter-plugin-template",
            output_dir=output_dir)
    # Raise an equivalent built-in exception so that packages wishing to catch
    # this error type do not have to explicitly depend on cookiecutter.
    except cookiecutter.exceptions.OutputDirExistsException as e:
        raise FileExistsError(str(e))
