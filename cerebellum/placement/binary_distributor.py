from typing import List

import numpy as np
from bsb import config
from bsb.placement.distributor import MorphologyDistributor
from scipy.spatial.transform import Rotation


@config.node
class DepthDistributorBinary(MorphologyDistributor, classmap_entry="depth_binary"):
    # Distribute different morphologies according to the depth.
    separating_h = config.attr(type=float)
    below_morphologies = config.attr(type=list)
    above_morphologies = config.attr(type=list)

    def distribute(self, positions, loaders, context):

        # We assume there is only one partition
        p = context.partitions[0]
        # Relative positions of the cells wrt the ldc of the partition
        depth = positions[:, 2] - p.ldc[2]
        metas = [l.get_meta() for l in loaders]
        # We compute the heights of the cells, wrt the bottom of the partition
        below_ids = []
        above_ids = []
        # Find below morphos id
        for i, m in enumerate(metas):
            if m["name"] in self.below_morphologies:
                below_ids.append(i)
            else:
                above_ids.append(i)
        below_rng = np.random.choice(below_ids, size=len(below_ids))
        above_rng = np.random.choice(above_ids, size=len(above_ids))

        # We compute the heights of the cells, wrt the bottom of the partition
        morpho_heights = [m["mdc"][2] - m["ldc"][2] for m in metas]
        morphos_id = np.zeros((len(positions),), dtype=int)
        above_separating = depth > self.separating_h
        morphos_id[above_separating] = 1
        morphos_id[above_separating] = above_rng
        morphos_id[~above_separating] = below_rng

        return morphos_id
