# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import os
import os.path

import pandas as pd


# Multiple types of visualizations (index.html, index.tsv).
def most_common_viz(output_dir: str, ints: collections.Counter) -> None:
    df = pd.DataFrame(ints.most_common(), columns=["Integer", "Frequency"])

    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write('<html><body>\n')
        fh.write('<h3>Most common integers:</h3>\n')
        fh.write(df.to_html(index=False))
        fh.write('</body></html>')

    with open(os.path.join(output_dir, 'index.tsv'), 'w') as fh:
        fh.write(df.to_csv(sep='\t', index=False))


# Multiple html files (index.1.html, index.2.html)
def multi_html_viz(output_dir: str, ints: list) -> None:
    ints = [str(i) for i in ints]
    with open(os.path.join(output_dir, 'index.1.html'), 'w') as fh:
        fh.write('<html><body>\n')
        fh.write(' '.join(ints))
        fh.write('</body></html>')
    with open(os.path.join(output_dir, 'index.2.html'), 'w') as fh:
        fh.write('<html><body>\n')
        fh.write(' '.join(reversed(ints)))
        fh.write('</body></html>')


# No input artifacts, only parameters.
def params_only_viz(output_dir: str, name: str = 'Foo Bar', age: int = 42):
    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write('<html><body>\n')
        fh.write('Name: %s\n' % name)
        fh.write('Age: %s\n' % age)
        fh.write('</body></html>')


# No input artifacts or parameters.
def no_input_viz(output_dir: str):
    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write('<html><body>\n')
        fh.write('Hello, World!\n')
        fh.write('</body></html>')


# Multiple input artifacts and parameters, and a nested directory with required
# resources for rendering.
def mapping_viz(output_dir: str, mapping1: dict, mapping2: dict,
                key_label: str, value_label: str) -> None:
    df1 = _dict_to_dataframe(mapping1, key_label, value_label)
    df2 = _dict_to_dataframe(mapping2, key_label, value_label)

    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write('<html><head>')
        fh.write('<link rel="stylesheet" type="text/css" '
                 'href="css/style.css" />')
        fh.write('</head>\n')
        fh.write('<body>\n')
        fh.write('<h3>mapping1:</h3>\n')
        fh.write(df1.to_html(index=False, classes='dummy-class'))
        fh.write('<h3>mapping2:</h3>\n')
        fh.write(df2.to_html(index=False, classes='dummy-class'))
        fh.write('</body></html>')

    css_dir = os.path.join(output_dir, 'css')
    os.mkdir(css_dir)
    with open(os.path.join(css_dir, 'style.css'), 'w') as fh:
        fh.write(_css)


def _dict_to_dataframe(dict_, key_label, value_label):
    return pd.DataFrame(sorted(dict_.items()),
                        columns=[key_label, value_label])


# Example table styling taken from http://www.w3schools.com/css/css_table.asp
_css = """
.dummy-class {
    border-collapse: collapse;
    width: 100%;
}

.dummy-class th, td {
    padding: 8px;
    text-align: left;
    border-bottom: 1px solid #ddd;
}
"""
