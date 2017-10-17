# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.sdk


class Context:
    def __init__(self, parent=None):
        self._parent = parent
        self._scope = None

    def import_action(self, plugin: str, action: str):
        pm = qiime2.sdk.PluginManager()
        action = pm.plugins[plugin].actions[action]
        return action.bind(lambda: Context(parent=self))

    def make_artifact(self, type, view, view_type=None):
        artifact = qiime2.sdk.Artifact.import_data(type, view, view_type)
        self._scope.add_reference(artifact)
        return artifact

    def __enter__(self):
        if self._scope is not None:
            # Prevent odd things from happening to lifecycle cleanup
            raise Exception('Cannot enter a context twice.')
        self._scope = Scope(self)
        return self._scope

    def __exit__(self, exc_type, exc_value, exc_tb):
        if exc_type is not None:
            # Something went wrong, teardown everything
            self._scope.destroy()
        else:
            # Everything is fine, just cleanup internal references and hoist
            self._scope.destroy(internal_only=True)
            if self._parent is not None:
                for hoisted in self._scope.iter_hoisted():
                    self._parent._scope.add_reference(hoisted)


class Scope:
    def __init__(self, ctx):
        self.ctx = ctx
        self._locals = []
        self._hoists = []

    def add_reference(self, ref, hoist=False):
        if hoist:
            self._hoists.append(ref)
        else:
            self._locals.append(ref)

    def iter_hoisted(self):
        yield from self._hoists

    def destroy(self, internal_only=False):
        destructors = self._locals
        if not internal_only:
            destructors = destructors + self._hoists
        for ref in destructors:
            ref._destructor()

        self.ctx = None
