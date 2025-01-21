import numpy as np
from bsb import AfterConnectivityHook, config, refs, types


@config.node
class MultiplyConnContacts(AfterConnectivityHook):
    """
    BSB postprocessing node to create multiple synapse per connection pair.
    """

    conn_strat = config.ref(refs.connectivity_ref, required=True)
    """Connection Strategies for which to multiply the synapses."""

    contacts = config.attr(type=types.distribution(), default=1)
    """Number or distribution determining the amount of synaptic contacts one cell will form on another"""

    def postprocess(self):
        for pre_ct in self.conn_strat.presynaptic.cell_types:
            pre_ps = pre_ct.get_placement_set()
            for post_ct in self.conn_strat.postsynaptic.cell_types:
                post_ps = post_ct.get_placement_set()
                for cs_name in self.conn_strat.get_output_names(pre_ct, post_ct):
                    # Draw the number of connection to create
                    n = int(self.contacts.draw(1)[0])
                    cs = self.scaffold.get_connectivity_set(cs_name)
                    if n > 1:
                        pre_locs, post_locs = cs.load_connections().all()
                        pre_locs = np.repeat(pre_locs, n - 1, axis=0)
                        post_locs = np.repeat(post_locs, n - 1, axis=0)
                        cs.connect(pre_ps, post_ps, pre_locs, post_locs)
