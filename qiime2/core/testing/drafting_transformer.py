# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

# import pandas as pd

# from q2_types.feature_data import TSVTaxonomyFormat

# # Transformer
# def dataframe_to_tsv_format(df):
#     ff = TSVTaxonomyFormat()
#     df.to_csv(str(ff), sep='\t', header=True, index=True)
#     return ff

# # Testing transformer
# def test_dataframe_to_tsv_format(self):
#     index = pd.Index(['seq1', 'seq2'], name='Feature ID', dtype=object)
#     columns = ['Taxon', 'Foo', 'Bar']
#     df = pd.DataFrame([['taxon1', '42', 'foo'], ['taxon2', '43', 'bar']],
#                         index=index, columns=columns, dtype=object)
#     exp = (
#         'Feature ID\tTaxon\tFoo\tBar\n'
#         'seq1\ttaxon1\t42\tfoo\n'
#         'seq2\ttaxon2\t43\tbar\n'
#     )

#     transformer = self.get_transformer(pd.DataFrame, TSVTaxonomyFormat)
#     obs = transformer(df)

#     with obs.open() as fh:
#         self.assertEqual(fh.read(), exp)