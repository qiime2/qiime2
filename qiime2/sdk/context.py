# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.core.cache import Cache
import qiime2.sdk
from qiime2.sdk.cache_config import CACHE_CONFIG
from qiime2.sdk.result import Artifact, Visualization


class Context:
    def __init__(self, parent=None):
        self._parent = parent
        self._scope = None

    def get_action(self, plugin: str, action: str):
        """Return a function matching the callable API of an action.

        This function is aware of the pipeline context and manages its own
        cleanup as appropriate.
        """
        pm = qiime2.sdk.PluginManager()
        plugin = plugin.replace('_', '-')
        try:
            plugin_obj = pm.plugins[plugin]
        except KeyError:
            raise ValueError("A plugin named %r could not be found." % plugin)

        try:
            action_obj = plugin_obj.actions[action]
        except KeyError:
            raise ValueError("An action named %r was not found for plugin %r"
                             % (action, plugin))
        # This factory will create new Contexts with this context as their
        # parent. This allows scope cleanup to happen recursively.
        # A factory is necessary so that independent applications of the
        # returned callable recieve their own Context objects.
        return action_obj._bind(lambda: Context(parent=self))

    def make_artifact(self, type, view, view_type=None):
        """Return a new artifact from a given view.

        This artifact is automatically tracked and cleaned by the pipeline
        context.
        """
        artifact = qiime2.sdk.Artifact.import_data(type, view, view_type)
        # self._scope WILL be defined at this point, as pipelines always enter
        # a scope before deferring to plugin code. (Otherwise cleanup wouldn't
        # happen)
        self._scope.add_reference(artifact)
        return artifact

    def __enter__(self):
        """For internal use only.

        Creates a scope API that can track references that need to be
        destroyed.

        """
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
            # Everything is fine, just cleanup internal references and pass
            # ownership off to the parent context.
            parent_refs = self._scope.destroy(local_references_only=True)
            if self._parent is not None:
                for ref in parent_refs:
                    self._parent._scope.add_reference(ref)


class Scope:
    def __init__(self, ctx):
        self.ctx = ctx
        self._locals = []
        self._parent_locals = []


    def add_reference(self, ref):
        """Add a reference to something destructable that is owned by this
           scope.
        """
        self._locals.append(ref)

    # NOTE: We end up with both the artifact and the pipeline alias of artifact
    # in the named cache in the end. We only have the pipeline alias in the
    # process cache
    def add_parent_reference(self, ref):
        """Add a reference to something destructable that will be owned by the
           parent scope. The reason it needs to be tracked is so that on
           failure, a context can still identify what will (no longer) be
           returned.
        """
        cache = Cache.get_cache()
        process_pool = cache.create_pool(process_pool=True, reuse=True)
        process_pool.save(ref)

        if CACHE_CONFIG.named_pool is not None:
            CACHE_CONFIG.named_pool.save(ref)

        self._parent_locals.append(ref)

        # Return an artifact backed by the data in the cache
        return process_pool.load(ref)

    def destroy(self, local_references_only=False):
        """Destroy all references and clear state.

        Parameters
        ----------
        local_references_only : bool
            Whether to destroy references that will belong to the parent scope.
        Returns
        -------
        list
            The list of references that were not destroyed.

        """
        local_refs = self._locals
        parent_refs = self._parent_locals

        # Unset instance state, handy to prevent cycles in GC, and also causes
        # catastrophic failure if some invariant is violated.
        del self._locals
        del self._parent_locals
        del self.ctx

        for ref in local_refs:
            # NOTE: This is getting a little weird. We're creating an instance
            # of a pool object, but we are not creating the pool itself, the
            # pool is on disk. We create a pool object referring to the
            # existing pool on disk then remove an item from it. If there is no
            # applicable pool to remove it from at this point, we want to
            # explode. There should always be one if everything worked out
            if isinstance(ref, Artifact) or isinstance(ref, Visualization):
                CACHE_CONFIG.cache.create_pool(process_pool=True, reuse=True).remove(ref)
            ref._destructor()

        if local_references_only:
            return parent_refs

        for ref in parent_refs:
            ref._destructor()

        return []
