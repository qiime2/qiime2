# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.sdk


class ProxyAction:
    """This class basically exists because .async/multiprocess.map is a thing
    """
    def __init__(self, action, locals):
        self._action = action
        self._locals = locals

    def __call__(self, *args, **kwargs):
        results = self._action(*args, **kwargs)
        self._locals += results  # Our own garbage collection scheme
        return results

    def __repr__(self):
        return "<ProxyAction of %r>" % self._action


class Context:
    def __init__(self):
        self._locals = []

    def import_action(self, plugin: str, action: str):
        # This singleton is shared across multiprocess on OS's with
        # copy-on-write semantics for `fork()` which means it isn't
        # reinitalized. This results in actions with the PID of the parent.
        # So instead, copy the action and "adopt" the right PID.
        pm = qiime2.sdk.PluginManager()
        action = pm.plugins[plugin].actions[action].copy()
        return ProxyAction(action, self._locals)

    def make_artifact(self, type, view, view_type=None):
        a = qiime2.sdk.Artifact.import_data(type, view, view_type=view_type)
        # This is needed so that in a subprocess situation, we can find
        # "internal" artifacts which will not be cleaned up correctly due to
        # `sys._exit`
        self._locals.append(a)
        return a
