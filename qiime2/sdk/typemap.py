import random
from itertools import combinations
import networkx as nx
import matplotlib.pyplot as plt

import qiime 2

q2methods = qiime2.sdk.util.actions_by_input_type('FeatureTable[Frequency]')

#############testing##################
q2p = qiime2.sdk.PluginManager()

#for plugin, value in q2p.plugins.items():
#    print(value.methods.keys())
##############testing#################



class TypeMap:
    """
    This class is is used to generate a graph of all methods in the Qiime2
    plugin and to create random testing paths from the methods.

    """


    def getNextParam(method, input=True):
        """Get a tuple of required/nonrequired inputs or outputs for each method

        Parameters
        ----------
        method : Qiime2.method (check the right signature for this)
            The method to get the semantic types from
        input : {True,False}
            Delineates if getting the method input or output types

        Returns
        -------
        Tuple - list of required and optional semantic types

        """

        req = []
        non_req = []

        if input:
            for k,v in method.signature.inputs.items():
                if not v.has_default():
