import numpy as np
from bsb import ConnectionStrategy, HemitypeCollection, config


@config.node
class ConnectomePC_DCN(ConnectionStrategy):
    divergence = config.attr(type=int, required=True)

    # We need to connect a PC to #divergence DCN from the whole region,
    # so we just return all the chunks
    def get_region_of_interest(self, chunk):
        ct = self.postsynaptic.cell_types[0]
        chunks = ct.get_placement_set().get_all_chunks()
        return chunks

    # We need to override the default _get_connect_args_from_job function because
    # We need to get single-post-chunk, multi-pre-chunk ROI instead of the opposite.
    def _get_connect_args_from_job(self, chunk, roi):
        pre = HemitypeCollection(self.presynaptic, chunk)
        post = HemitypeCollection(self.postsynaptic, roi)
        return pre, post

    def connect(self, pre, post):
        for pre_ps in pre.placement:
            for post_ps in post.placement:
                self._connect_type(pre_ps, post_ps)

    def _connect_type(self, pre_ps, post_ps):

        pc_pos = pre_ps.load_positions()
        dcn_pos = post_ps.load_positions()
        n_pc = len(pc_pos)
        n_dcn = len(dcn_pos)
        max_synapses = n_pc * self.divergence
        pre_locs = np.full((max_synapses, 3), -1, dtype=int)
        post_locs = np.full((max_synapses, 3), -1, dtype=int)

        ptr = 0

        # We connect each PC to #divegence DCN
        for j, _ in enumerate(pc_pos):

            # Select randomly #divegence DCNs from all the DCNs
            selected_dcn_ids = np.random.choice(n_dcn, self.divergence, replace=False)
            pre_locs[ptr : ptr + self.divergence, 0] = j
            post_locs[ptr : ptr + self.divergence, 0] = selected_dcn_ids
            ptr = ptr + self.divergence

        self.connect_cells(pre_ps, post_ps, pre_locs, post_locs)
        # print("Connected", n_pc, "PCs to", n_dcn, "dcns fibers.")
