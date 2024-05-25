"""
    Module for the configuration node of the Glomerulus to UBC ConnectionStrategy
"""

import numpy as np
from bsb import ConnectionStrategy, config, ConfigurationError

from cerebellum.connectome.presyn_dist_strat import PresynDistStrat


@config.node
class ConnectomeGlomerulusUBC(PresynDistStrat, ConnectionStrategy):
    """
    BSB Connection strategy to connect any type of Glomerulus to UBC cells.
    """

    ratios_ubc: dict[str, float] = config.dict(type=float, required=True)
    """Dictionary that links a postsynaptic celltype name to the ratios of the presynaptic UBC 
    population that connects to it."""

    def boot(self):
        parsed_ratios = {k.name: (0. if k.name not in self.ratios_ubc else self.ratios_ubc[k.name]) for k in self.presynaptic.cell_types}
        sum_ = np.nansum(list(parsed_ratios.values()))
        if np.any(np.array(list(parsed_ratios.values())) < 0):
            raise ConfigurationError(
                "Presynaptic cell type ratios should be greater than 0"
            )
        if sum_ == 0:
            raise ConfigurationError(
                "At least one presynaptic ratio should be greater than 0"
            )
        for k, v in parsed_ratios.items():
            parsed_ratios[k] = v / sum_
        self.ratios_ubc = parsed_ratios.copy()

    def connect(self, pre, post):
        for post_ps in post.placement:
            ubc_pos = post_ps.load_positions()
            ubc_ids = np.random.permutation(len(ubc_pos))
            cum_sum = 0
            for pre_ps in pre.placement:
                # select the ratio of random ubc ids to connect
                loc_ratio = self.ratios_ubc[pre_ps.cell_type.name]
                ubc_ids = ubc_ids[cum_sum : cum_sum + int(len(ubc_pos) // loc_ratio)]
                cum_sum = int(len(ubc_pos) // loc_ratio)

                glom_pos = pre_ps.load_positions()
                selected_ubc = ubc_pos[ubc_ids]

                n_conn = len(selected_ubc)  # one ubc_glom per ubc
                pre_locs = np.full((n_conn, 3), -1, dtype=int)
                post_locs = np.full((n_conn, 3), -1, dtype=int)
                connected_gloms = np.full(len(glom_pos), False, dtype=bool)
                for i, ubc in enumerate(selected_ubc):
                    distances = np.linalg.norm(ubc - glom_pos[~connected_gloms], axis=1)
                    avaiable_gloms_ids = np.where(~connected_gloms)[0]
                    id_sorted_distances = np.argsort(distances)
                    sorted_distances = distances[id_sorted_distances]
                    mask = sorted_distances <= self.radius
                    if np.any(mask):
                        pre_id = avaiable_gloms_ids[id_sorted_distances[mask][0]]
                        post_locs[i] = ubc_ids[i]
                        pre_locs[i] = pre_id
                        connected_gloms[pre_id] = True

                self.connect_cells(pre_ps, post_ps, pre_locs, post_locs)
