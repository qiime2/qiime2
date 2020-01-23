# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


def load_jupyter_server_extension(nb_server):
    from .handlers import QIIME2RedirectHandler, QIIME2ResultHandler
    from notebook.utils import url_path_join

    result_store = {}
    app = nb_server.web_app

    def route(path):
        return url_path_join(app.settings['base_url'], 'qiime2', path)

    app.add_handlers(r'.*', [
        (route(r'redirect'), QIIME2RedirectHandler,
         {'result_store': result_store}),
        (route(r'view/(.*)'), QIIME2ResultHandler,
         # This *is* odd, but it's because we are tricking StaticFileHandler
         {'path': result_store,
          'default_filename': 'index.html'})
    ])
