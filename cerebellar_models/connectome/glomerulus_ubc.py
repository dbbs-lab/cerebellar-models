"""
Module for the configuration node of the Glomerulus to UBC ConnectionStrategy
"""

import numpy as np
from bsb import ConfigurationError, ConnectionStrategy, config

from cerebellar_models.connectome.presyn_dist_strat import PresynDistStrat


@config.node
class ConnectomeGlomerulusUBC(PresynDistStrat, ConnectionStrategy):
    """
    BSB Connection strategy to connect any type of Glomerulus to UBC cells.
    """

    ratios_ubc: dict[str, float] = config.dict(type=float, required=True)
    """Dictionary that links a postsynaptic celltype name to the ratios of the presynaptic UBC 
        population that connects to it."""

    def boot(self):
        parsed_ratios = {
            k.name: (0.0 if k.name not in self.ratios_ubc else self.ratios_ubc[k.name])
            for k in self.presynaptic.cell_types
        }
        sum_ = np.nansum(list(parsed_ratios.values()))
        if np.any(np.array(list(parsed_ratios.values())) < 0):
            raise ConfigurationError("Presynaptic cell type ratios should be greater than 0")
        if sum_ == 0:
            raise ConfigurationError("At least one presynaptic ratio should be greater than 0")
        for k, v in parsed_ratios.items():
            parsed_ratios[k] = v / sum_
        self.ratios_ubc = parsed_ratios.copy()

    def connect(self, pre, post):
        connected_gloms = {
            pre_ps.cell_type.name: np.full(len(pre_ps), False, dtype=bool)
            for pre_ps in pre.placement
        }

        for post_ps in post.placement:
            ubc_pos = post_ps.load_positions()
            ubc_ids = np.random.permutation(len(ubc_pos))
            ubc_pos = ubc_pos[ubc_ids]
            cum_sum = 0
            loc_ratio = 0
            for pre_ps in pre.placement:
                # select the ratio of random ubc ids to connect
                pre_ct = pre_ps.cell_type.name
                loc_ratio += self.ratios_ubc[pre_ct]
                new_ptr = min(len(ubc_ids), int(np.round(len(ubc_pos) * loc_ratio)))
                loc_ubc_ids = ubc_ids[cum_sum:new_ptr]
                cum_sum = new_ptr

                glom_pos = pre_ps.load_positions()
                selected_ubc = ubc_pos[loc_ubc_ids]

                n_conn = len(selected_ubc)  # one glom per ubc
                pre_locs = np.full((n_conn, 3), -1, dtype=int)
                post_locs = np.full((n_conn, 3), -1, dtype=int)
                for i, ubc in enumerate(selected_ubc):
                    filter_connected = ~connected_gloms[pre_ct]
                    if not np.any(filter_connected):
                        # if no glom is left un-connected
                        # we reset the connected list
                        connected_gloms[pre_ct][:] = False
                        filter_connected = ~connected_gloms[pre_ct]
                    avaiable_gloms_ids = np.where(filter_connected)[0]
                    distances = np.linalg.norm(ubc - glom_pos[filter_connected], axis=1)
                    # Take closest glom
                    pre_id = avaiable_gloms_ids[np.argsort(distances)[0]]
                    post_locs[i, 0] = loc_ubc_ids[i]
                    pre_locs[i, 0] = pre_id
                    connected_gloms[pre_ct][pre_id] = True

                self.connect_cells(pre_ps, post_ps, pre_locs, post_locs)
