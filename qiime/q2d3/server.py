import tornado
import tornado.ioloop
import tornado.httpclient
import tornado.web
import webbrowser
import os
import io
import json
import time
import subprocess
import atexit
import uuid

from qiime.q2d3.context import Q2D3Context
from qiime.q2d3 import get_plugin_workflow


__UPLOADS__ = "uploads/"
RELATIVE = os.path.split(__file__)[0]
GB = 1024 * 1024 * 1024
class Index(tornado.web.RequestHandler):
    def get(self):
        self.render("static/index.html")

    def delete(self):
        context_manager = Q2D3Context(os.getcwd())
        uuid_str = self.get_query_argument('uuid')
        fp = context_manager.data[uuid.UUID(uuid_str)]
        os.remove(fp)

class Workflows(tornado.web.RequestHandler):
    def get(self):
        self.render("static/workflows.html")


class Artifacts(tornado.web.RequestHandler):
    def get(self):
        self.render("static/artifacts.html")

def instantiator(port):
    class Job(tornado.web.RequestHandler):
        def get(self, plugin, workflow):
            plugin, workflow = get_plugin_workflow(plugin, workflow)
            self.render("static/job.html", plugin=plugin, workflow=workflow)

        def post(self, plugin, workflow):
            plugin, workflow = get_plugin_workflow(plugin, workflow)
            input_artifacts = {}
            input_params = {}
            output_names = {}
            for argument in self.request.arguments:
                value = self.get_argument(argument)
                if argument.startswith('in-'):
                    name = argument[3:]
                    input_artifacts[name] = uuid.UUID(value)
                elif argument.startswith('param-'):
                    name = argument[6:]
                    input_params[name] = value
                elif argument.startswith('out-'):
                    name = argument[4:]
                    if not value:
                        raise TypeError("Must enter an output name for the artifact.")
                    output_names[name] = value

            context_manager = Q2D3Context(os.getcwd(), output_names=output_names)
            job = context_manager(workflow, input_artifacts, input_params)
            path = str(uuid.uuid4()) + '.md'
            with open(path, mode='w') as fh:
                fh.write(job)

            self.redirect('http://localhost:%d/notebooks/%s' % (port, path))

    return Job
#
# @tornado.web.stream_request_body
# class Upload(tornado.web.RequestHandler):
#
#     def post(self):
#         self._file.close()
#         self.finish(json.dumps({}))
#
#     def prepare(self):
#         self._boundary = b'--' + ''.join(self.request.headers.get('Content-Type').split('=')[1:]).encode('ascii')
#         self._end = self._boundary + b'--'
#
#         self._file = open(self.get_query_argument('type'), mode='wb')
#
#     def data_received(self, data):
#         start = data[:len(self._boundary)]
#         is_boundary = self._boundary == start
#         if(is_boundary):
#             index = data.find(b'\r\n\r\n')
#             data = data[index+4:]
#         else:
#             index = data.find(b'\r\n' + self._end)
#             if index > 0:
#                 data = data[:index]
#
#         self._file.write(data)




def start_server(port=4444):
    jport = port + 1
    application = tornado.web.Application([
            (r"/", Index),
            (r"/workflows", Workflows),
            (r"/job/(.*)/(.*)", instantiator(jport)),
            (r"/artifacts", Artifacts),
            (r'/static/(.*)', tornado.web.StaticFileHandler,
             {'path': RELATIVE + '/static/'}),
            # (r"/upload", Upload),
            ], debug=True)

    application.listen(port, max_buffer_size=10 * GB)
    process = subprocess.Popen("jupyter notebook --no-browser --port %d --port-retries 0" % jport, shell=True)
    @atexit.register
    def shutdown():
        subprocess.call("ps -ef | awk '$3 == \"%d\" {print $2}' | xargs kill -15" % process.pid, shell=True)
        process.kill()

    webbrowser.open('http://localhost:%d' % port, new=1, autoraise=True)
    tornado.ioloop.IOLoop.instance().start()
