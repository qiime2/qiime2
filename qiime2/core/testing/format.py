# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import TextFileFormat
import qiime2.plugin.model as model


class IntSequenceFormat(TextFileFormat):
    """
    A sequence of integers stored on new lines in a file. Since this is a
    sequence, the integers have an order and repetition of elements is allowed.

    """
    def sniff(self):
        with self.open() as fh:
            for line, _ in zip(fh, range(5)):
                try:
                    int(line.rstrip('\n'))
                except (TypeError, ValueError):
                    return False
            return True


class MappingFormat(TextFileFormat):
    """
    A mapping of keys to values stored in a TSV file. Since this is a
    mapping, key-value pairs do not have an order and duplicate keys are
    disallowed.

    """
    def sniff(self):
        with self.open() as fh:
            for line, _ in zip(fh, range(5)):
                cells = line.rstrip('\n').split('\t')
                if len(cells) != 2:
                    return False
            return True


class SingleIntFormat(TextFileFormat):
    """
    Exactly one int on a single line in the file.

    """
    def sniff(self):
        with self.open() as fh:
            try:
                int(fh.readline().rstrip('\n'))
            except (TypeError, ValueError):
                return False
        return True


IntSequenceDirectoryFormat = model.SingleFileDirectoryFormat(
    'IntSequenceDirectoryFormat', 'ints.txt', IntSequenceFormat)


# This could have been a `SingleFileDirectoryFormat`, but isn't for testing
# purposes
class MappingDirectoryFormat(model.DirectoryFormat):
    mapping = model.File('mapping.tsv', format=MappingFormat)


class FourIntsDirectoryFormat(model.DirectoryFormat):
    """
    A sequence of exactly four integers stored across multiple files, some of
    which are in a nested directory. Each file contains a single integer.
    Since this is a sequence, the integers have an order (corresponding to
    filename) and repetition of elements is allowed.

    """
    single_ints = model.FileCollection(
        r'file[1-2]\.txt|nested/file[3-4]\.txt', format=SingleIntFormat)

    @single_ints.set_path_maker
    def single_ints_path_maker(self, num):
        if not 0 < num < 5:
            raise ValueError("`num` must be 1-4, not %r." % num)
        if num > 2:
            return 'nested/file%d.txt' % num
        else:
            return 'file%d.txt' % num
