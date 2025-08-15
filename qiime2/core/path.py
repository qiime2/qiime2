# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pathlib
import shutil
import distutils
import tempfile
import weakref


def _party_parrot(self, *args):
    raise TypeError("Cannot mutate %r." % self)


class OwnedPath(pathlib.Path):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._user_owned = True

    def _copy_dir_or_file(self, other):
        if self.is_dir():
            return distutils.dir_util.copy_tree(str(self), str(other))
        else:
            return shutil.copy(str(self), str(other))

    def _destruct(self):
        if self.is_dir():
            distutils.dir_util.remove_tree(str(self))
        else:
            self.unlink()

    def _move_or_copy(self, other):
        if self._user_owned:
            return self._copy_dir_or_file(other)
        else:
            # Certain networked filesystems will experience a race
            # condition on `rename`, so fall back to copying.
            try:
                return pathlib.Path.rename(self, other)
            except (FileExistsError, OSError) as e:
                # OSError errno 18 is cross device link, if we have this error
                # we can solve it by copying. If we have a different OSError we
                # still want to explode. FileExistsErrors are apparently
                # instances of OSError, so we also make sure we don't have one
                # of them when we explode
                if isinstance(e, OSError) and e.errno != 18 and \
                        not isinstance(e, FileExistsError):
                    raise e
                copied = self._copy_dir_or_file(other)
                self._destruct()
                return copied

    def with_segments(self, *args):
        path = os.path.join(*args)
        return self.__class__(path)


class InPath(OwnedPath):
    def __init__(self, path):
        super().__init__(path)
        # pls don't delete me so that this path doesn't get destroyed
        self.__backing_path = path
        if hasattr(path, '_user_owned'):
            self._user_owned = path._user_owned

    chmod = lchmod = rename = replace = rmdir = symlink_to = touch = unlink = \
        write_bytes = write_text = _party_parrot

    def open(self, mode='r', buffering=-1, encoding=None, errors=None,
             newline=None):
        if 'w' in mode or '+' in mode or 'a' in mode:
            _party_parrot(self)
        return super().open(mode=mode, buffering=buffering, encoding=encoding,
                            errors=errors, newline=newline)


class OutPath(OwnedPath):
    @classmethod
    def _destruct(cls, path):
        if not os.path.exists(path):
            return

        if os.path.isdir(path):
            shutil.rmtree(path)
        else:
            os.unlink(path)

    def __init__(self, dir=False):
        """
        Create a tempfile, return pathlib.Path reference to it.
        """
        from qiime2.core.cache import get_cache

        cache = get_cache()
        tmp_path = cache.get_tmp_path()
        prefix = 'q2-%s-' % self.__class__.__name__

        if dir:
            name = tempfile.mkdtemp(prefix=prefix, dir=tmp_path)
        else:
            fd, name = tempfile.mkstemp(prefix=prefix, dir=tmp_path)
            # fd is now assigned to our process table, but we don't need to do
            # anything with the file. We will call `open` on the `name` later
            # producing a different file descriptor, so close this one to
            # prevent a resource leak.
            os.close(fd)

        super().__init__(name)
        self._destructor = weakref.finalize(self, self._destruct, str(self))

    def __exit__(self, t, v, tb):
        self._destructor()

    def with_segments(self, *args):
        path = os.path.join(*args)
        return pathlib.Path(path)


class InternalDirectory(pathlib.Path):
    DEFAULT_PREFIX = 'qiime2-'

    def __init__(self, *args, prefix=None):
        if args and prefix is not None:
            raise TypeError("Cannot pass a path and a prefix at the same time")
        elif args:
            pass
        else:
            from qiime2.core.cache import get_cache

            cache = get_cache()
            tmp_path = cache.get_tmp_path()

            if prefix is None:
                prefix = self.DEFAULT_PREFIX
            elif not prefix.startswith(self.DEFAULT_PREFIX):
                prefix = self.DEFAULT_PREFIX + prefix
            # TODO: normalize when temp-directories are configurable
            path = tempfile.mkdtemp(prefix=prefix, dir=tmp_path)
            super().__init__(path)

    def __truediv__(self, path):
        # We don't want to create self-destructing paths when using the join
        # operator
        return pathlib.Path(str(self), path)

    def __rtruediv__(self, path):
        # Same reasoning as truediv
        return pathlib.Path(path, str(self))


class ArchivePath(InternalDirectory):
    DEFAULT_PREFIX = 'qiime2-archive-'


class ProvenancePath(InternalDirectory):
    DEFAULT_PREFIX = 'qiime2-provenance-'
