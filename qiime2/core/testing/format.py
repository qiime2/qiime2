# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import TextFileFormat, ValidationError
import qiime2.plugin.model as model


class IntSequenceFormat(TextFileFormat):
    """
    A sequence of integers stored on new lines in a file. Since this is a
    sequence, the integers have an order and repetition of elements is allowed.

    """
    def validate(self):
        with self.open() as fh:
            for line, idx in zip(fh, range(1, 6)):
                try:
                    int(line.rstrip('\n'))
                except (TypeError, ValueError):
                    raise ValidationError("Line %d is not an integer." % idx)


class IntSequenceFormatV2(IntSequenceFormat):
    """
    Same as IntSequenceFormat, but has a header "VERSION 2"

    """
    def validate(self):
        with self.open() as fh:
            if fh.readline() != 'VERSION 2\n':
                raise ValidationError("Missing header: VERSION 2")


class MappingFormat(TextFileFormat):
    """
    A mapping of keys to values stored in a TSV file. Since this is a
    mapping, key-value pairs do not have an order and duplicate keys are
    disallowed.

    """
    def validate(self):
        with self.open() as fh:
            for line, idx in zip(fh, range(1, 6)):
                cells = line.rstrip('\n').split('\t')
                if len(cells) != 2:
                    raise ValidationError("Line %d does not have exactly 2 "
                                          "elements seperated by a tab." % idx)


class SingleIntFormat(TextFileFormat):
    """
    Exactly one int on a single line in the file.

    """
    def validate(self):
        with self.open() as fh:
            try:
                int(fh.readline().rstrip('\n'))
            except (TypeError, ValueError):
                raise ValidationError("File does not contain an integer")
            if fh.readline():
                raise ValidationError("Too many lines in file.")


IntSequenceDirectoryFormat = model.SingleFileDirectoryFormat(
    'IntSequenceDirectoryFormat', 'ints.txt', IntSequenceFormat)

IntSequenceV2DirectoryFormat = model.SingleFileDirectoryFormat(
    'IntSequenceV2DirectoryFormat', 'integers.txt', IntSequenceFormatV2)


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


class RedundantSingleIntDirectoryFormat(model.DirectoryFormat):
    """
    Two files of SingleIntFormat which are exactly the same.

    """
    int1 = model.File('file1.txt', format=SingleIntFormat)
    int2 = model.File('file2.txt', format=SingleIntFormat)

    def validate(self):
        if self.int1.view(int) != self.int2.view(int):
            raise ValidationError("file1.txt does not match file2.txt")


class UnimportableFormat(TextFileFormat):
    """
    Unimportable format used for testing.
    """


UnimportableDirectoryFormat = model.SingleFileDirectoryFormat(
    'UnimportableDirectoryFormat', 'ints.txt', UnimportableFormat)
