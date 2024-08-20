# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.core.type.util import is_collection_type
from qiime2.core.type import HashableInvocation
from qiime2.core.cache import get_cache
import qiime2.sdk
from qiime2.sdk.parallel_config import PARALLEL_CONFIG


class Context:
    def __init__(self, parent=None, parallel=False):
        if parent is not None:
            self.action_executor_mapping = parent.action_executor_mapping
            self.executor_name_type_mapping = parent.executor_name_type_mapping
            self.parallel = parent.parallel
            self.cache = parent.cache
        else:
            self.action_executor_mapping = \
                PARALLEL_CONFIG.action_executor_mapping
            # Cast type to str so yaml doesn't think it needs tpo instantiate
            # an executor object when we write this to then read this from
            # provenance
            self.executor_name_type_mapping = \
                None if PARALLEL_CONFIG.parallel_config is None \
                else {v.label: v.__class__.__name__
                      for v in PARALLEL_CONFIG.parallel_config.executors}
            self.parallel = parallel
            self.cache = get_cache()
            # Only ever do this on the root context. We only want to index the
            # pool once before we start adding our own stuff to it.
            if self.cache.named_pool is not None:
                self.cache.named_pool.create_index()

        self._parent = parent
        self._scope = None

    def get_action(self, plugin: str, action: str):
        """Return a function matching the callable API of an action.
        This function is aware of the pipeline context and manages its own
        cleanup as appropriate.
        """
        plugin = plugin.replace('_', '-')
        plugin_action = plugin + ':' + action

        pm = qiime2.sdk.PluginManager()
        try:
            plugin_obj = pm.plugins[plugin]
        except KeyError:
            raise ValueError("A plugin named %r could not be found." % plugin)

        try:
            action_obj = plugin_obj.actions[action]
        except KeyError:
            raise ValueError(
                "An action named %r was not found for plugin %r"
                % (action, plugin))

        # We return this callable which determines whether to return cached
        # results or to run the action requested.
        def deferred_action(*args, **kwargs):
            # The function is the first arg, we ditch that
            args = args[1:]

            # If we have a named_pool, we need to check for cached results that
            # we can reuse.
            #
            # We can short circuit our index checking if any of our arguments
            # are proxies because if we got a proxy as an argument, we know it
            # is a new thing we are computing from a prior step in the pipeline
            # and thus will not be cached. We can only have proxies if we are
            # executing with parsl
            if self.cache.named_pool is not None and (not self.parallel or (
                    self.parallel and not self._contains_proxies(
                        *args, **kwargs))):
                collated_inputs = action_obj.signature.collate_inputs(
                    *args, **kwargs)
                callable_args = action_obj.signature.coerce_user_input(
                    **collated_inputs)

                # Make args and kwargs look how they do when we read them out
                # of a .yaml file (list of single value dicts of
                # input_name: value)
                arguments = []
                for k, v in callable_args.items():
                    arguments.append({k: v})

                invocation = HashableInvocation(plugin_action, arguments)
                if invocation in self.cache.named_pool.index:
                    cached_outputs = self.cache.named_pool.index[invocation]
                    loaded_outputs = {}

                    for name, _type in action_obj.signature.outputs.items():
                        if is_collection_type(_type.qiime_type):
                            loaded_collection = qiime2.sdk.ResultCollection()
                            cached_collection = cached_outputs[name]

                            # Get the order we should load collection items in
                            collection_order = list(cached_collection.keys())
                            self._validate_collection(collection_order)
                            collection_order.sort(key=lambda x: x.idx)

                            for elem_info in collection_order:
                                elem = cached_collection[elem_info]
                                loaded_elem = self.cache.named_pool.load(elem)
                                loaded_collection[
                                    elem_info.item_name] = loaded_elem

                            loaded_outputs[name] = loaded_collection
                        else:
                            output = cached_outputs[name]
                            loaded_outputs[name] = \
                                self.cache.named_pool.load(output)

                    return qiime2.sdk.Results(
                        loaded_outputs.keys(), loaded_outputs.values())

            # If we didn't have cached results to reuse, we need to execute the
            # action.
            #
            # These factories will create new Contexts with this context as
            # their parent. This allows scope cleanup to happen recursively. A
            # factory is necessary so that independent applications of the
            # returned callable receive their own Context objects.
            #
            # The parsl factory is a bit more complicated because we need to
            # pass this exact Context along for a while longer until we run a
            # normal _bind in action/_run_parsl_action. Then we create a new
            # Context with this one as its parent inside of the parsl app
            def _bind_parsl_context(ctx):
                def _bind_parsl_args(*args, **kwargs):
                    return action_obj._bind_parsl(ctx, *args, **kwargs)
                return _bind_parsl_args

            if self.parallel:
                return _bind_parsl_context(self)(*args, **kwargs)

            return action_obj._bind(
                lambda: Context(parent=self))(*args, **kwargs)

        deferred_action = action_obj._rewrite_wrapper_signature(
            deferred_action)
        action_obj._set_wrapper_properties(deferred_action)
        return deferred_action

    def _contains_proxies(self, *args, **kwargs):
        """Returns True if any of the args or kwargs are proxies
        """
        return any(isinstance(arg, qiime2.sdk.proxy.Proxy) for arg in args) \
            or any(isinstance(value, qiime2.sdk.proxy.Proxy) for
                   value in kwargs.values())

    def _validate_collection(self, collection_order):
        """Validate that all indexed items in the collection agree on how
        large the collection should be and that we have that many elements.
        """
        assert all([elem.total == collection_order[0].total
                    for elem in collection_order])
        assert len(collection_order) == collection_order[0].total

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
        self.entries = 1
        self._locals = []
        self._parent_locals = []

    def add_reference(self, ref):
        """Add a reference to something destructable that is owned by this
           scope.
        """
        # if self.ctx._parent:
        #     self._parent_locals.append(ref)
        # else:
        self._locals.append(ref)

    # NOTE: We end up with both the artifact and the pipeline alias of artifact
    # in the named cache in the end. We only have the pipeline alias in the
    # process pool
    def add_parent_reference(self, ref):
        """Add a reference to something destructable that will be owned by the
           parent scope. The reason it needs to be tracked is so that on
           failure, a context can still identify what will (no longer) be
           returned.
        """
        new_ref = self.ctx.cache.process_pool.save(ref)

        if self.ctx.cache.named_pool is not None:
            self.ctx.cache.named_pool.save(new_ref)

        self._parent_locals.append(new_ref)
        self._parent_locals.append(ref)

        # Return an artifact backed by the data in the cache
        return new_ref

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
        self.entries -= 1

        if self.entries != 0:
            return []

        ctx = self.ctx
        local_refs = self._locals
        parent_refs = self._parent_locals

        # Unset instance state, handy to prevent cycles in GC, and also causes
        # catastrophic failure if some invariant is violated.
        del self.ctx
        del self._locals
        del self._parent_locals

        for ref in local_refs:
            ref._destructor()

        if local_references_only:
            return parent_refs

        for ref in parent_refs:
            ref._destructor()

        ctx.cache.garbage_collection()
        return []
