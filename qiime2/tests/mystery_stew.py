# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

# Note, the naming of this file doesn't include a `test_` prefix, this is the
# easiest way to ensure nose/pytest doesn't pick this up during its automatic
# test discovery

from qiime2.plugins import ArtifactAPIUsage

import pytest


def _labeler(val):
    if hasattr(val, 'id'):
        return val.id
    return val


def get_tests():
    import qiime2.sdk
    tests = []

    try:
        pm = qiime2.sdk.PluginManager.reuse_existing()
    except qiime2.sdk.UninitializedPluginManagerError:
        import os

        if 'MYSTERY_STEW' in os.environ:
            from q2_mystery_stew.plugin_setup import create_plugin

            the_stew = create_plugin()
            pm = qiime2.sdk.PluginManager(add_plugins=False)
            pm.add_plugin(the_stew)

        pm = qiime2.sdk.PluginManager()

    try:
        plugin = pm.plugins['mystery-stew']
    except KeyError:
        return tests
    for action in plugin.actions.values():
        for name in action.examples:
            tests.append((action, name))
    return tests


@pytest.mark.parametrize('action,example', get_tests(), ids=_labeler)
def test_mystery_stew(action, example):
    example_f = action.examples[example]
    use = ArtifactAPIUsage(enable_assertions=True)
    example_f(use)
    rendered = use.render()
    ctx = use.get_example_data()
    exec(rendered, ctx)
