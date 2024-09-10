# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import sys
import errno
import shutil
import importlib

import threading
import contextlib

_REDIRECTED_STDIO_LOCK = threading.Lock()


@contextlib.contextmanager
def redirected_stdio(stdout=None, stderr=None):
    """A context manager to redirect stdio to a new file (if provided).

    Parameters
    ----------
    stdout : file-like object or file-descriptor (int)
      The file to redirect stdout to. Does nothing to the process-table if not
      provided.
    stderr : file-like object or file-descriptor (int)
      The file to redirect stderr to. Does nothing to the process-table if not
      provided.

    Notes
    -----
    The current sys.stdout/sys.stderr must be backed by a file descriptor
    (i.e they cannot be replaced with BytesIO).
    """
    with _REDIRECTED_STDIO_LOCK:
        if stdout is not None:
            with _redirected_fd(to=stdout, stdio=sys.stdout):
                if stderr is not None:
                    with _redirected_fd(to=stderr, stdio=sys.stderr):
                        yield
                else:
                    yield
        elif stderr is not None:
            with _redirected_fd(to=stderr, stdio=sys.stderr):
                yield
        else:
            yield


# Taken whole-sale from: http://stackoverflow.com/a/22434262/579416
@contextlib.contextmanager
def _redirected_fd(to=os.devnull, stdio=None):
    if stdio is None:
        stdio = sys.stdout

    stdio_fd = _get_fileno(stdio)
    # copy stdio_fd before it is overwritten
    # NOTE: `copied` is inheritable on Windows when duplicating a standard
    # stream
    with os.fdopen(os.dup(stdio_fd), 'wb') as copied:
        stdio.flush()  # flush library buffers that dup2 knows nothing about
        try:
            os.dup2(_get_fileno(to), stdio_fd)  # $ exec >&to
        except ValueError:  # filename
            with open(to, 'wb') as to_file:
                os.dup2(to_file.fileno(), stdio_fd)  # $ exec > to
        try:
            yield stdio  # allow code to be run with the redirected stdio
        finally:
            # restore stdio to its previous value
            # NOTE: dup2 makes stdio_fd inheritable unconditionally
            stdio.flush()
            os.dup2(copied.fileno(), stdio_fd)  # $ exec >&copied


def _get_fileno(file_or_fd):
    fd = getattr(file_or_fd, 'fileno', lambda: file_or_fd)()
    if not isinstance(fd, int):
        raise ValueError("Expected a file (`.fileno()`) or a file descriptor")
    return fd


def duplicate(src, dst):
    """Alternative to :py:func:`shutil.copyfile`
    this will use :py:func:`os.link` when possible.

    See :py:func:`shutil.copyfile` for documention.
    Only ``src`` and ``dst`` are supported.
    Unlike copyfile, this will not overwrite the destination if it exists.

    """
    if os.path.isdir(src):
        # os.link will give a permission error
        raise OSError(errno.EISDIR, "Is a directory", src)
    if os.path.isdir(dst):
        # os.link will give a FileExists error
        raise OSError(errno.EISDIR, "Is a directory", dst)

    if os.path.exists(dst):
        # shutil.copyfile will overwrite the existing file
        raise OSError(errno.EEXIST, "File exists", src, "File exists", dst)

    try:
        os.link(src, dst)
    except OSError as e:
        if e.errno == errno.EXDEV:  # Invalid cross-device link
            shutil.copyfile(src, dst)
        elif e.errno == errno.EPERM:  # Permissions/ownership error
            shutil.copyfile(src, dst)
        elif e.errno == errno.ENOTSUP:  # Operation not supported
            shutil.copyfile(src, dst)
        else:
            raise


def get_filepath_from_package(package, relative_filepath):
    """Utility function for getting the absolute path to a packaged data asset.

    Parameters
    ----------
    package : str
        The name of the package or sub-package where the file can be found.
    relative_filepath : str
        The filepath to look up, relative to `package`.

    Returns
    -------
    fp : pathlib.Path
        The absolute filepath to the requested data asset.

    Raises
    ------
    FileNotFoundError
        If the filepath does not exist.

    """
    fp = importlib.resources.files(package) / relative_filepath
    if not fp.exists():
        raise FileNotFoundError(
            f'The requested data asset {fp} does not exist.')
    else:
        return fp