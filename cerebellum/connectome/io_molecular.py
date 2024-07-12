import numpy as np
from bsb import ConnectionStrategy, config, refs
from bsb.mixins import NotParallel


@config.node
class ConnectomeIO_MLI(NotParallel, ConnectionStrategy):
    mli_pc_connectivity: ConnectionStrategy = config.ref(refs.connectivity_ref, required=True)
    io_pc_connectivity: ConnectionStrategy = config.ref(refs.connectivity_ref, required=True)
    pre_cell_pc = config.ref(refs.cell_type_ref, required=True)
    # pre_cell_io_pc = config.reflist(refs.cell_type_ref, required=True)

    @config.property
    def depends_on(self):
        # Get the possibly missing `_depends_on` list.
        deps = getattr(self, "_depends_on", None) or []
        # Strat is required, but depends on a reference that isn't available when the config loads.
        strat_io = getattr(self, "io_pc_connectivity", None)
        strat_mli = getattr(self, "mli_pc_connectivity", None)
        return [*{*deps, strat_io, strat_mli}]

    @depends_on.setter
    def depends_on(self, value):
        self._depends_on = value

    def connect(self, pre, post):
        for pre_ps in pre.placement:
            for post_ps in post.placement:
                self._connect_type(pre_ps, post_ps)

    def _connect_type(self, pre_ps, post_ps):
        mli_pos = post_ps.load_positions()

        # All the MLIs connected to the same PC must be connected to the IO cell which
        # is connected to that PC

        # We retrieve the connectivity data for the mli-pc connectivity
        cs = self.mli_pc_connectivity.get_output_names(post_ps.cell_type, self.pre_cell_pc)
        assert (
            len(cs) == 1
        ), f"Only one connection set should be given from {self.mli_pc_connectivity.name}."
        cs = self.scaffold.get_connectivity_set(cs[0])
        mli_per_pc_list = []

        loc_mli, loc_pc = cs.load_connections().as_globals().all()
        len_pc = len(self.scaffold.get_placement_set(self.pre_cell_pc))
        unique_pc = np.arange(0, len_pc, 1)

        for current in unique_pc:
            mli_ids = np.where(loc_pc[:, 0] == current)[0]
            mli_per_pc_list.append(loc_mli[mli_ids][:, 0])

        # We retrieve the connectivity data for the io-pc connectivity
        cs = self.io_pc_connectivity.get_output_names(pre_ps.cell_type, self.pre_cell_pc)
        assert (
            len(cs) == 1
        ), f"Only one connection set should be given from {self.io_pc_connectivity.name}."
        cs = self.scaffold.get_connectivity_set(cs[0])
        io_pc_list = []

        loc_io, loc_pc = cs.load_connections().as_globals().all()

        for current in unique_pc:
            io_ids = np.where(loc_pc[:, 0] == current)[0]
            io_pc_list.append(loc_io[io_ids][:, 0])

        n_purkinje = len(io_pc_list)
        max_connections = len(mli_pos) * n_purkinje
        if max_connections > 0:
            pre_locs = np.full((max_connections, 3), -1, dtype=int)
            post_locs = np.full((max_connections, 3), -1, dtype=int)

            ptr = 0

            for j, id_io in enumerate(io_pc_list):
                ln = len(mli_per_pc_list[j])
                pre_locs[ptr : ptr + ln, 0] = id_io
                post_locs[ptr : ptr + ln, 0] = mli_per_pc_list[j]
                ptr = ptr + ln

            self.connect_cells(pre_ps, post_ps, pre_locs, post_locs)
