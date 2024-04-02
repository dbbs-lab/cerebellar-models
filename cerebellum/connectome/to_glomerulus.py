import itertools
import numpy as np
from bsb import config, Chunk, InvertedRoI, ConnectionStrategy
from scipy.stats.distributions import truncexpon


@config.node
class ConnectomeGlomerulus(InvertedRoI, ConnectionStrategy):

    def connect(self, pre, post):
        # We use a truncated exponential distribution to favour the presynaptic fibers closer the
        # postsynaptic glomerulus.
        exp_dist = truncexpon(b=0.99)
        for pre_ps in pre.placement:
            for post_ps in post.placement:
                presyn_pos = pre_ps.load_positions()
                glomeruli_pos = post_ps.load_positions()
                n_glom = len(glomeruli_pos)

                # We work assuming that there is at least 1 presynaptic fiber in the ROI.
                # Otherwise, there is something wrong in the placement phase.
                pre_locs = np.full((n_glom, 3), -1, dtype=int)
                post_locs = np.full((n_glom, 3), -1, dtype=int)

                # Here truncexpon(b=0.99).rvs(size=n) provides n random numbers between [0, 1[ with
                # more likelihood for the lowest values. These n numbers are converted into ids.
                rolls = np.floor(len(presyn_pos) * exp_dist.rvs(size=n_glom)).astype(int)
                # We connect each glomerulus to a mossy fiber.
                for j, glomerulus in enumerate(glomeruli_pos):
                    dist = np.linalg.norm(glomerulus - presyn_pos, axis=1)
                    id_sorted_dist = np.argsort(dist)
                    # Sorting MF ids by distances so low ids are closer to glom.
                    # MF closer to the granule have a higher chance to be picked
                    pre_locs[j, 0] = id_sorted_dist[rolls[j]]
                    post_locs[j, 0] = j

                self.connect_cells(pre_ps, post_ps, pre_locs, post_locs)

@config.node
class ConnectomeMossyGlomerulus(ConnectomeGlomerulus):
    x_length = config.attr(type=float, required=True)
    z_length = config.attr(type=float, required=True)

    def get_region_of_interest(self, chunk):
        # Chunk here is a postsynaptic chunk because of InvertedRoI
        # We look for chunks containing mossy fibers that are within a rectangle of size
        # x_length * z_length centered on the postsynaptic chunk containing the glomerulus.
        chunks = set(
            itertools.chain.from_iterable(
                ct.get_placement_set().get_all_chunks() for ct in self.presynaptic.cell_types
            )
        )
        selected_chunks = []
        for c in chunks:
            x_dist = np.fabs(chunk[0] - c[0])
            z_dist = np.fabs(chunk[1] - c[1])
            x_dist = x_dist * chunk.dimensions[0]
            z_dist = z_dist * chunk.dimensions[2]

            if (x_dist < self.x_length / 2) and (z_dist < self.z_length / 2):
                selected_chunks.append(Chunk([c[0], c[1], c[2]], chunk.dimensions))
        return selected_chunks
