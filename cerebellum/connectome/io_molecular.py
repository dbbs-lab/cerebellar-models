"""
    Module for the configuration node of the IO to molecular layer interneurons (MLI) ConnectionStrategy
"""

import numpy as np
from bsb import ConnectionStrategy, config, refs
from bsb.mixins import NotParallel


@config.node
class ConnectomeIO_MLI(NotParallel, ConnectionStrategy):
    """
    BSB Connection strategy to connect IO cells to molecular layer interneurons (MLI) cells through PC.
    IO cells which are connected to a PC should also connect to all the MLIs connected to this PC.
    """

    io_pc_connectivity: ConnectionStrategy = config.ref(refs.connectivity_ref, required=True)
    """Connection Strategy that links IO to PC."""
    mli_pc_connectivity = config.reflist(refs.connectivity_ref, required=True)
    """List of Connection Strategies that links MLI to PC."""
    pre_cell_pc = config.ref(refs.cell_type_ref, required=True)
    """Celltype used for to represent PC."""

    @config.property
    def depends_on(self):
        # Get the possibly missing `_depends_on` list.
        deps = getattr(self, "_depends_on", None) or []
        # Strat is required, but depends on a reference that isn't available when the config loads.
        strat_io = getattr(self, "io_pc_connectivity", None)
        strat_mli = getattr(self, "mli_pc_connectivity", None)
        return [*{*deps, strat_io, *strat_mli}]

    @depends_on.setter
    def depends_on(self, value):
        self._depends_on = value

    def connect(self, pre, post):
        for pre_ps in pre.placement:
            self._connect_type(pre_ps, post)

    def load_connectivity_set(self, connection_strat, cell_type):
        """
        Load the connection locations from a connection strategy that connects the provided cell type to PC.

        :param bsb.connectivity.strategy.ConnectionStrategy connection_strat: Connection strategy to load.
        :param str cell_type: Presynaptic cell type name.
        :return: A tuple containing:
                - an array of the presynaptic cell_type connection locations,
                - an array of the postsynaptic pc connection locations
        """
        cs = connection_strat.get_output_names(cell_type, self.pre_cell_pc)
        assert (
            len(cs) == 1
        ), f"Only one connection set should be given from {connection_strat.name}."
        cs = self.scaffold.get_connectivity_set(cs[0])
        return cs.load_connections().as_globals().all()

    def load_connections_mli_pc(self, post):
        """
        Load the connection locations for all the MLI to PC strategies.
        Will only keep one connection location information for each unique pair of MLI-PC.

        :param bsb.connectivity.strategy.HemitypeCollection post: Postsynaptic hemitype
        :return: A tuple containing:
                - an array of the presynaptic mli connection locations
                - an array of the postsynaptic pc connection locations
                - an array of the postsynaptic placement set indexes
        """
        loc_all_mli = []
        loc_all_pc_mli = []
        ct_ps_ids = []
        # For each mli type
        for i, post_ps in enumerate(post.placement):
            for strat in self.mli_pc_connectivity:
                try:
                    # Fetch the corresponding connectivity set info if it exists
                    loc_mli, loc_pc_mli = self.load_connectivity_set(strat, post_ps.cell_type)
                except ValueError:
                    continue
                # keep one unique pair of mli to pc
                to_keep = np.unique(
                    np.concatenate([[loc_mli[:, 0]], [loc_pc_mli[:, 0]]]).T,
                    axis=0,
                    return_index=True,
                )[1]
                loc_all_mli.extend(loc_mli[to_keep])
                loc_all_pc_mli.extend(loc_pc_mli[to_keep])
                ct_ps_ids.extend(np.repeat(i, len(to_keep)))

        return np.asarray(loc_all_mli), np.asarray(loc_all_pc_mli), np.asarray(ct_ps_ids)

    def _connect_type(self, pre_ps, post):
        # We retrieve the connectivity data for the mli-pc connectivity and io-pc connectivity
        loc_io, loc_pc = self.load_connectivity_set(self.io_pc_connectivity, pre_ps.cell_type)
        loc_mli, loc_mli_pc, ct_ps_ids = self.load_connections_mli_pc(post)

        u_purkinje = np.unique(loc_pc[:, 0])
        io_pc_list = []
        mli_per_pc_list = []
        grouped_ps_ids = []
        for current in u_purkinje:
            io_pc_list.append(loc_io[loc_pc[:, 0] == current][:, 0])
            mli_ids = loc_mli_pc[:, 0] == current
            mli_per_pc_list.append(loc_mli[mli_ids][:, 0])
            grouped_ps_ids.append(ct_ps_ids[mli_ids])

        max_len = np.max(
            [[len(mli_ids), len(io_ids)] for mli_ids, io_ids in zip(mli_per_pc_list, io_pc_list)],
            axis=0,
        )
        max_connections = np.prod(max_len) * u_purkinje.size
        pre_locs = np.full((max_connections, 3), -1, dtype=int)
        post_locs = np.full((max_connections, 3), -1, dtype=int)
        ps_locs = np.full(max_connections, -1, dtype=int)
        ptr = 0
        for i, io_ids in enumerate(io_pc_list):
            ln = len(mli_per_pc_list[i])
            for current in io_ids:
                pre_locs[ptr : ptr + ln, 0] = current
                post_locs[ptr : ptr + ln, 0] = mli_per_pc_list[i]
                ps_locs[ptr : ptr + ln] = grouped_ps_ids[i]
                ptr = ptr + ln

        # because we are using global indices we need to extract the global ps
        pre_ps = self.scaffold.get_placement_set(pre_ps.cell_type)
        for i, post_ps in enumerate(post.placement):
            post_ps = self.scaffold.get_placement_set(post_ps.cell_type)
            current = (ps_locs == i)[:ptr]
            self.connect_cells(pre_ps, post_ps, pre_locs[:ptr][current], post_locs[:ptr][current])
