# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.core.cache import get_cache
import qiime2.sdk
from qiime2.sdk.parsl_config import PARSL_CONFIG
from qiime2.sdk.result import Artifact, Visualization


class Context:
    def __init__(self, parent=None, parsl=False):
        if parent is not None:
            self.action_executor_mapping = parent.action_executor_mapping
            self.parsl = parent.parsl
            self.cache = parent.cache
        else:
            self.action_executor_mapping = PARSL_CONFIG.action_executor_mapping
            self.parsl = parsl
            self.cache = get_cache()

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
        def _bind_parsl_context(ctx):
            def _bind_parsl_args(*args, **kwargs):
                return action_obj._bind_parsl(ctx, *args, **kwargs)
            return _bind_parsl_args

        if self.parsl:
            return _bind_parsl_context(self)

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
            if self._parent is not None and self._parent._scope is not None:
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

    def add_parent_reference(self, ref):
        """Add a reference to something destructable that will be owned by the
           parent scope. The reason it needs to be tracked is so that on
           failure, a context can still identify what will (no longer) be
           returned.
        """
        # Add the ref to the process pool and get back a new ref backed by the
        # cache and the pool
        # NOTE: There are some problems with how we create this reference right
        # now. With a normal archive the archiver path looks something like
        # /qiime2-archive-<random_stuff>/uuid/the data and all that
        # With these pool refs it's just
        # /uuid/the data and all that
        # This discrepency is likely causing issues. For example, if we pass
        # pool_ref to named_pool.save we end up with the guts of the artifact
        # in the data directory (VERSION and metadata.yaml etc. in the top
        # level of the data directory)
        pool_ref = self.ctx.cache.process_pool.save(ref)

        # Add the ref to the named pool if one exists
        if self.ctx.cache.named_pool is not None:
            self.ctx.cache.named_pool.save(ref)

        self._parent_locals.append(ref)

        # Return an artifact backed by the data in the cache
        return pool_ref

    # TODO: Demote refs when they are aliased and remove those demoted refs
    # from the pool
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
        ctx = self.ctx

        # Unset instance state, handy to prevent cycles in GC, and also causes
        # catastrophic failure if some invariant is violated.
        del self._locals
        del self._parent_locals
        del self.ctx

        for ref in local_refs:
            ref._destructor()

            if isinstance(ref, (Artifact, Visualization)):
                ctx.cache.process_pool.remove(ref)

        if local_references_only:
            return parent_refs

        for ref in parent_refs:
            ref._destructor()

            if isinstance(ref, (Artifact, Visualization)):
                ctx.cache.process_pool.remove(ref)

        return []
