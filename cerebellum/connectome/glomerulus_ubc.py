"""
    Module for the configuration node of the Glomerulus to UBC ConnectionStrategy
"""

import numpy as np
from bsb import ConfigurationError, ConnectionStrategy, ConnectivityError, config

from cerebellum.connectome.presyn_dist_strat import PresynDistStrat


@config.node
class ConnectomeGlomerulusUBC(PresynDistStrat, ConnectionStrategy):
    """
    BSB Connection strategy to connect any type of Glomerulus to UBC cells.
    """

    ratios_ubc: dict[str, float] = config.dict(type=float, required=True)
    """Dictionary that links a postsynaptic celltype name to the ratios of the presynaptic UBC 
        population that connects to it."""
    force_ratio: bool = config.attr(type=bool, required=False, default=False)
    """Flag to force indegree ratios to be maintained. Otherwise (default behavior, priority is given to the pre-post 
        distance constraint."""

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
        tot_pre = {}
        pre_cts = []
        connected_gloms = {}
        for pre_ps in pre.placement:
            pre_ct = pre_ps.cell_type.name
            connected_gloms[pre_ct] = np.full(len(pre_ps), False, dtype=bool)
            tot_pre[pre_ct] = len(pre_ps)
            pre_cts.append(pre_ct)
        tot_post = {ps.cell_type.name: len(ps) for ps in post.placement}
        if sum(tot_pre.values()) < sum(tot_post.values()):
            if self.force_ratio:
                raise ConnectivityError(
                    "Not enough presynpatic cells to match the provided ratios."
                )
            # If too many postsynaptic cells with respect to pre,
            # we reduce the number of connectable postsynaptic cells
            tot_post[post.cell_type.name] *= sum(tot_pre.values()) / sum(tot_post.values())

        if self.force_ratio:
            for pre_ct in pre_cts:
                if sum(tot_post.values()) * self.ratios_ubc[pre_ct] > tot_pre[pre_ct]:
                    raise ConnectivityError(
                        f"Not enough presynpatic cells of type {pre_ct} to match the provided ratios."
                    )
        for post_ps in post.placement:
            ubc_pos = post_ps.load_positions()
            ubc_ids = np.random.permutation(len(ubc_pos))[: tot_post[post_ps.cell_type.name]]
            ubc_pos = ubc_pos[ubc_ids]
            cum_sum = 0
            for pre_ct, pre_ps in zip(pre_cts, pre.placement):
                # select the ratio of random ubc ids to connect
                loc_ratio = self.ratios_ubc[pre_ct]
                loc_ubc_ids = ubc_ids[cum_sum : cum_sum + int(np.round(len(ubc_pos) * loc_ratio))]
                cum_sum = int(np.round(len(ubc_pos) * loc_ratio))

                glom_pos = pre_ps.load_positions()
                selected_ubc = ubc_pos[loc_ubc_ids]

                n_conn = len(selected_ubc)  # one glom per ubc
                pre_locs = np.full((n_conn, 3), -1, dtype=int)
                post_locs = np.full((n_conn, 3), -1, dtype=int)
                i = 0
                for ubc in selected_ubc:
                    distances = np.linalg.norm(ubc - glom_pos[~connected_gloms[pre_ct]], axis=1)
                    avaiable_gloms_ids = np.where(~connected_gloms[pre_ct])[0]
                    id_sorted_distances = np.argsort(distances)
                    sorted_distances = distances[id_sorted_distances]
                    if np.any(sorted_distances <= self.radius) or self.force_ratio:
                        pre_id = avaiable_gloms_ids[id_sorted_distances[0]]
                        post_locs[i] = loc_ubc_ids[i]
                        pre_locs[i] = pre_id
                        connected_gloms[pre_ct][pre_id] = True
                        i += 1

                self.connect_cells(pre_ps, post_ps, pre_locs[:i], post_locs[:i])
