# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import urllib.parse


def make_html(location):
    url = "/qiime2/redirect?location={location}".format(
        location=urllib.parse.quote(location))
    # This is dark magic. An image has an onload handler, which let's me
    # grab the parent dom in an anonymous way without needing to scope the
    # output cells of Jupyter with some kind of random ID.
    # Using transparent pixel from: https://stackoverflow.com/a/14115340/579416
    return ('<div><img onload="({anon_func})(this.parentElement, \'{url}\')"'
            ' src="data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAAL'
            'AAAAAABAAEAAAICRAEAOw==" /></div>'.format(
                anon_func=_anonymous_function, url=url))


# 404 - the extension isn't installed
# 428 - the result went out of scope, re-run cell
# 302->200 - set up the iframe for that location
_anonymous_function = '''\
function(div, url){
if (typeof require !== 'undefined') {
    var baseURL = require.toUrl('').split('/').slice(0, -2).join('/');
} else {
    var baseURL = JSON.parse(
        document.getElementById('jupyter-config-data').innerHTML
    ).baseUrl.slice(0, -1);
}
url = baseURL + url;
fetch(url).then(function(res) {
    if (res.status === 404) {
        div.innerHTML = 'Install QIIME 2 Jupyter extension with:<br />' +
                        '<code>jupyter serverextension enable --py qiime2' +
                        ' --sys-prefix</code><br />then restart your server.' +
                        '<br /><br />(Interactive output not available on ' +
                        'static notebook viewer services like nbviewer.)';
    } else if (res.status === 409) {
        div.innerHTML = 'Visualization no longer in scope. Re-run this cell' +
                        ' to see the visualization.';
    } else if (res.ok) {
        url = res.url;
        div.innerHTML = '<iframe src=\\'' + url + '\\' style=\\'' +
                        'width: 100%; height: 700px; border: 0;\\'>' +
                        '</iframe><hr />Open in a: <a href=\\'' + url + '\\'' +
                        ' target=\\'_blank\\'>new window</a>'
    } else {
        div.innerHTML = 'Something has gone wrong. Check notebook server for' +
                        ' errors.';
    }
});
}'''
