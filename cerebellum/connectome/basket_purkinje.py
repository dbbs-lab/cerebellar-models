from types import CellType
import numpy as np
from bsb.connectivity.strategy import ConnectionStrategy
from bsb.connectivity.strategy import HemitypeCollection
from bsb.storage import Chunk
from bsb import config
from bsb.cell_types import CellType
from itertools import chain
from bsb.reporting import report
from scipy.stats.distributions import truncexpon

@config.node
class ConnectomeBasketPurkinje(ConnectionStrategy):
    axon_length = config.attr(type=int, required=True)
    basket_radius = config.attr(type=int, required=True)
    affinity = config.attr(type=int, required=True)

    def get_region_of_interest(self, chunk):

        ct = self.presynaptic.cell_types[0]
        chunks = ct.get_placement_set().get_all_chunks()
        selected_chunks = []

        #We consider only chunks whose y is less than or equal to the y of the current one
        for c in chunks:
            #if ((c[1]+1)*chunk.dimensions[1] - self.scaffold.partitions.granular_layer.thickness >=0 and (c[1])*chunk.dimensions[1] < self.scaffold.partitions.purkinje_layer.thickness + self.scaffold.partitions.granular_layer.thickness):
            if ((c[1])*chunk.dimensions[1] <= self.scaffold.partitions.granular_layer.thickness + self.scaffold.partitions.purkinje_layer.thickness  and (c[1]+1)*chunk.dimensions[1] < self.scaffold.partitions.purkinje_layer.thickness + self.scaffold.partitions.granular_layer.thickness + self.scaffold.partitions.b_molecular_layer.thickness):
                selected_chunks.append(Chunk([c[0], c[1], c[2]], chunk.dimensions))

        return chunks

    def connect(self, pre, post):
        for pre_ct, pre_ps in pre.placement.items():
            for post_ct, post_ps in post.placement.items():
                self._connect_type(pre_ct, pre_ps, post_ct, post_ps)

    def _connect_type(self, pre_ct, pre_ps, post_ct, post_ps):
              
        basket_pos = pre_ps.load_positions()
        purkinje_pos = post_ps.load_positions()
        n_basket = len(basket_pos)
        n_purkinje = len(purkinje_pos)

        print("N basket:", n_basket)
        print("N purkinje:", n_purkinje)

        #We consider only one synapse per pair.
        pre_locs = np.full((n_basket*n_purkinje, 3), -1, dtype=int)
        post_locs = np.full((n_purkinje*n_purkinje, 3), -1, dtype=int)
        ptr = 0

        for i,basket in enumerate(basket_pos):
            x_dist = np.fabs(purkinje_pos[:,0]-basket[0])
            p_dist = np.fabs(purkinje_pos[:,2]-basket[2])
            to_connect = (p_dist < self.basket_radius) & (x_dist < self.axon_length)
            #print(to_connect)

            pc_ids = np.nonzero(to_connect)[0]
            #print(pc_ids)

            if (len(pc_ids)>0):
                if (self.affinity < 1):
                    pc_ids = pc_ids[
                                    np.random.choice(
                                        np.max(
                                            [
                                                1,
                                                int(
                                                    np.floor(
                                                        self.affinity * len(pc_ids)
                                                    )
                                                ),
                                            ]
                                        ),
                                    )
                                ]
                
                #print(len(pc_ids))
                post_locs[ptr:ptr+len(pc_ids),0] = pc_ids
                pre_locs[ptr:ptr+len(pc_ids),0] = i
                ptr += len(pc_ids)
        print(pre_locs[:ptr])
        print(post_locs[:ptr])
        self.connect_cells(pre_ps, post_ps, pre_locs[:ptr], post_locs[:ptr])
        print("Connected", n_basket, "basket to", n_purkinje, "purkinje.")
        #print(pre_locs)
        #print(post_locs)


