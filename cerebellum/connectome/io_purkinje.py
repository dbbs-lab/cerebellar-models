"""
    Module for the configuration node of the IO to purkinje cells ConnectionStrategy
"""

import numpy as np
from bsb import FixedIndegree, config, types


@config.node
class ConnectomeIoPurkinje(FixedIndegree):
    """
    BSB Connection strategy to connect IO to Purkinje cells with multiple synapses per pair.
    """

    contacts = config.attr(type=types.distribution(), default=1)
    """Number or distribution determining the amount of synaptic contacts one cell will form on another"""

    def connect(self, pre, post):
        super().connect(pre, post)
        for cs_name in self.get_output_names():
            # Draw the number of connection to create
            n = int(self.contacts.draw(1)[0])
            cs = self.scaffold.get_connectivity_set(cs_name)
            pre_ps = cs.pre_type.get_placement_set()
            post_ps = cs.pre_type.get_placement_set()
            if n > 1:
                pre_locs, post_locs = cs.load_connections().all()
                pre_locs = np.repeat(pre_locs, n - 1, axis=0)
                post_locs = np.repeat(post_locs, n - 1, axis=0)
                cs.connect(pre_ps, post_ps, pre_locs, post_locs)
