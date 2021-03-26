# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import sys
import errno
import shutil

import threading
import contextlib

_REDIRECTED_STDIO_LOCK = threading.Lock()


@contextlib.contextmanager
def redirected_stdio(stdout=None, stderr=None):
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
    """Alternative to shutil.copyfile, this will use os.link when possible.

    See shutil.copyfile for documention. Only `src` and `dst` are supported.
    Unlike copyfile, this will not overwrite the destination if it exists.

    """
    graft_file(src, dst)


def graft(src, dst, *, merge=False, overwrite=False, graft_within=False,
          remove_src=False, mkdirs=False, ignore=None):
    """Graft inodes onto another part of the filesystem (copying as needed).

    Your one-stop shop for migrating files from point A to point B (unless you
    really mean to duplicate the data on disk, in which case, don't use this).

    Cheat sheet:
        TODO

    Parameters
    ----------
    src : path (file or directory)
        The source to move, may be deleted by `remove_src`
    dst : path (file or directory)
        The destination recieving new files, may exist if `merge` or
        `graft_within` are True, otherwise, should not yet exist.
    merge : bool
        `src` must be a directory. If so, then the contents
        of `src` will be merged into `dst` (if it exists).
    overwrite : bool
        Whether files will be overwritten in the engraftment. This will not
        cause files to be replaced with directories or vice versa. Directories
        will never be overwritten (use `merge` or remove them before calling
        `graft`).
    graft_within : bool
        Places `src` *within* the `dst` directory. If `dst` does not exist,
        it will be assumed to be a directory.
    remove_src : bool
        Whether to remove the source file/directory after finishing (replicates
        `os.rename` behavior)
    mkdirs : bool
        Whether to create parent directories of `dst` as needed.
    ignore : filter
        ....

    Returns
    -------
    A list of (action, src, dst).
    """
    system_actions = []
    if dst == src:
        raise ValueError("destination cannot be the same as the source")

    if os.path.commonpath([src, dst]) == os.path.normpath(src):
        raise ValueError("destination cannot be nested within source")

    if not os.path.exists(src):
        raise OSError(errno.ENOENT, 'No such file or directory', src)

    if mkdirs:
        system_actions.append(
            lambda: os.makedirs(os.path.dirname(dst), exist_ok=True))
    elif os.path.dirname(dst) and not os.path.exists(os.path.dirname(dst)):
        raise OSError(errno.ENOENT, 'No such directory', os.path.dirname(dst))

    if graft_within:
        if os.path.isfile(dst):
            raise OSError(errno.ENOTDIR, "Not a directory", dst)
        elif not os.path.exists(dst):
            system_actions.append(lambda: os.mkdir(dst))

        dst = os.path.join(dst, os.path.basename(src))

    # dst may have changed
    if os.path.isdir(dst) and os.path.isdir(src) and not merge:
        raise ValueError("Cannot combine directory into an existing directory"
                         " without explicit `merge`. %s -> %s" % (src, dst))
    if merge and os.path.isfile(src):
        raise ValueError("Cannot use `merge` with a single file, use"
                         " `overwrite` instead.")

    if os.path.isfile(dst) and os.path.isdir(src):
        raise OSError(errno.EISDIR, "Is a directory", src)

    # most validation finished
    for action in system_actions:
        action()

    res = []
    if os.path.isdir(src):
        if remove_src and not merge:
            res.append(_checked_call(os.rename, src, dst))
        if not res:
            if not os.path.isdir(dst):
                os.mkdir(dst)
            for dir_entry in os.scandir(src):
                if ignore is not None and ignore(dir_entry):
                    continue
                merge = merge and dir_entry.is_dir()
                res.extend(graft(dir_entry.path, dst, graft_within=True,
                                 overwrite=overwrite, merge=merge,
                                 ignore=ignore))

        if remove_src and os.path.exists(src):
            shutil.rmtree(src)

    else:
        res.append(graft_file(src, dst, overwrite=overwrite,
                              remove_src=remove_src))
    return res


def _checked_call(func, src, dst, reraise=False):
    try:
        func(src, dst)
    except OSError:
        if reraise:
            raise
        return False
    return (func, src, dst)


def graft_file(src, dst, overwrite=False, remove_src=False):
    # Validate source
    if not os.path.exists(src):
        raise OSError(errno.ENOENT, 'No such file', src)
    elif os.path.isdir(src):
        raise OSError(errno.EISDIR, "Is a directory", src)

    # Validate destination
    if os.path.isfile(dst):
        if overwrite:
            os.remove(dst)
        else:
            raise OSError(errno.EEXIST, "File exists", dst)
    elif os.path.isdir(dst):
        raise OSError(errno.EISDIR, "Is a directory", dst)

    # Strategy: try, try, try, and try again
    succeeded = False
    # If we can destroy the source, then a rename is the fastest
    if remove_src:
        succeeded = _checked_call(os.rename, src, dst)
    # Try a hard-link
    if not succeeded:
        succeeded = _checked_call(os.link, src, dst)
    # Try a file copy
    if not succeeded:
        succeeded = _checked_call(shutil.copyfile, src, dst,
                                  reraise=True)
    # The engraftment is finished, cleanup source if needed
    if remove_src and os.path.exists(src):
        os.remove(src)

    return succeeded


def graft_iterable(iterable, dst, overwrite=False, remove_src=False,
                   mkdirs=False, ignore=None):
    for src in iterable:
        merge = True
        if os.path.isfile(src):
            merge = False
        graft(src, dst, merge=merge, graft_within=True,
              overwrite=overwrite,
              remove_src=remove_src,
              mkdirs=mkdirs,
              ignore=ignore)
