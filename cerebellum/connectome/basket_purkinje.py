from itertools import chain
from types import CellType

import numpy as np
from bsb import config
from bsb.cell_types import CellType
from bsb.connectivity.strategy import ConnectionStrategy, HemitypeCollection
from bsb.reporting import report
from bsb.storage._chunks import Chunk
from scipy.stats.distributions import truncexpon


@config.node
class ConnectomeBasketPurkinje(ConnectionStrategy):
    axon_length = config.attr(type=float, required=True)
    basket_radius = config.attr(type=float, required=True)
    affinity = config.attr(type=float, required=True)

    def get_region_of_interest(self, chunk):

        ct = self.presynaptic.cell_types[0]
        chunks = ct.get_placement_set().get_all_chunks()
        selected_chunks = []

        # We consider only chunks whose y is less than or equal to the y of the current one

        for c in chunks:
            if (c[2]) * chunk.dimensions[
                2
            ] <= self.scaffold.partitions.granular_layer.thickness + self.scaffold.partitions.purkinje_layer.thickness + chunk.dimensions[
                2
            ] and (
                c[2] + 1
            ) * chunk.dimensions[
                2
            ] < self.scaffold.partitions.purkinje_layer.thickness + self.scaffold.partitions.granular_layer.thickness + self.scaffold.partitions.t_molecular_layer.thickness:
                selected_chunks.append(Chunk([c[0], c[1], c[2]], chunk.dimensions))

        return selected_chunks

    def connect(self, pre, post):
        for pre_ps in pre.placement:
            for post_ps in post.placement:
                self._connect_type(pre_ps, post_ps)

    def _connect_type(self, pre_ps, post_ps):

        basket_pos = pre_ps.load_positions()
        purkinje_pos = post_ps.load_positions()
        n_basket = len(basket_pos)
        n_purkinje = len(purkinje_pos)

        # We consider only one synapse per pair.
        pre_locs = np.full((n_basket * n_purkinje, 3), -1, dtype=int)
        post_locs = np.full((n_purkinje * n_purkinje, 3), -1, dtype=int)
        ptr = 0

        for i, basket in enumerate(basket_pos):
            x_dist = np.fabs(purkinje_pos[:, 0] - basket[0])
            p_dist = np.fabs(purkinje_pos[:, 2] - basket[2])

            # We connect the purkinje cells satifying the gemetric condition;
            # We fix the direction of the axon using the index (even/odd) of the basket cells.
            # In this way we have the same number of basket cells with axons in direction +x or -x
            to_connect = (
                (p_dist < self.basket_radius)
                & (x_dist < self.axon_length)
                & (purkinje_pos[:, 0] > np.power(-1, i))
            )

            pc_ids = np.nonzero(to_connect)[0]
            if len(pc_ids) > 0:
                if self.affinity < 1:
                    pc_ids = pc_ids[
                        np.random.choice(
                            len(pc_ids),
                            np.max(
                                [
                                    1,
                                    int(np.floor(self.affinity * len(pc_ids))),
                                ]
                            ),
                            replace=False,
                        )
                    ]

                post_locs[ptr : ptr + len(pc_ids), 0] = pc_ids
                pre_locs[ptr : ptr + len(pc_ids), 0] = i
                ptr += len(pc_ids)

        self.connect_cells(pre_ps, post_ps, pre_locs[:ptr], post_locs[:ptr])
