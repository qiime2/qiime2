# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.plugin


IntSequence1 = qiime2.plugin.SemanticType('IntSequence1')
IntSequence2 = qiime2.plugin.SemanticType('IntSequence2')
Mapping = qiime2.plugin.SemanticType('Mapping')
FourInts = qiime2.plugin.SemanticType('FourInts')
SingleInt = qiime2.plugin.SemanticType('SingleInt')

Kennel = qiime2.plugin.SemanticType('Kennel', field_names='pet')
Dog = qiime2.plugin.SemanticType('Dog', variant_of=Kennel.field['pet'])
Cat = qiime2.plugin.SemanticType('Cat', variant_of=Kennel.field['pet'])
