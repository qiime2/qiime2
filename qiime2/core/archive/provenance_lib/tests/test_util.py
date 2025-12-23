# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from .testing_utilities import CustomAssertions


class CustomAssertionsTests(CustomAssertions):
    def test_assert_re_appears_only_once(self):
        t = ("Lick an orange. It tastes like an orange.\n"
             "The strawberries taste like strawberries!\n"
             "The snozzberries taste like snozzberries!")
        self.assertREAppearsOnlyOnce(t, 'Lick an orange')
        self.assertREAppearsOnlyOnce(t, 'tastes like')
        with self.assertRaisesRegex(AssertionError, 'Regex.*match.*orange'):
            self.assertREAppearsOnlyOnce(t, 'orange')
        with self.assertRaisesRegex(AssertionError, 'Regex.*taste like'):
            self.assertREAppearsOnlyOnce(t, 'taste like')
        with self.assertRaisesRegex(AssertionError, 'Regex.*snozzberries'):
            self.assertREAppearsOnlyOnce(t, 'snozzberries')
        with self.assertRaisesRegex(AssertionError, 'Regex.*!'):
            self.assertREAppearsOnlyOnce(t, '!')
