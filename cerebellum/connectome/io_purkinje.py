import numpy as np
from bsb import ConnectionStrategy, config
from bsb.mixins import NotParallel


@config.node
class ConnectomeIOPurkinje(NotParallel, ConnectionStrategy):

    def connect(self, pre, post):
        for pre_ps in pre.placement:
            for post_ps in post.placement:
                self._connect_type(pre_ps, post_ps)

    def _connect_type(self, pre_ps, post_ps):

        io_pos = pre_ps.load_positions()
        pc_pos = post_ps.load_positions()
        n_io = len(io_pos)
        n_pc = len(pc_pos)

        print("N io cells:", n_io)
        print("N pcs:", n_pc)

        # Each PC is connect to a single IO cell
        pre_locs = np.full((n_pc, 3), -1, dtype=int)
        post_locs = np.full((n_pc, 3), -1, dtype=int)

        # We connect each pc to an IO cell

        # TODO Select 1-3 IO to be connected using a negative exponential distribution
        # according to https://www.biorxiv.org/content/10.1101/2023.03.27.534425v2.full
        selected_io_ids = np.random.randint(0, high=n_io, size=n_pc, dtype=int)

        for j, _ in enumerate(pc_pos):

            pre_locs[j, 0] = selected_io_ids[j]
            post_locs[j, 0] = j

        self.connect_cells(pre_ps, post_ps, pre_locs, post_locs)
        print("Connected", n_io, "IO cells to", n_pc, "PCs.")
