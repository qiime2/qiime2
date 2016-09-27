# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import sys

import threading
import contextlib

_REDIRECTED_STDIO_LOCK = threading.Lock()


@contextlib.contextmanager
def redirected_stdio(stdout=None, stderr=None):
    with contextlib.ExitStack() as stack:
        stack.enter_context(_REDIRECTED_STDIO_LOCK)
        if stdout is not None:
            stack.enter_context(_redirected_fd(to=stdout, stdio=sys.stdout))
        if stderr is not None:
            stack.enter_context(_redirected_fd(to=stderr, stdio=sys.stderr))
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
