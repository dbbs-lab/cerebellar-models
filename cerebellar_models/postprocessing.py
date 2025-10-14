from bsb import AfterConnectivityHook,config,DatasetNotFoundError,ConnectivityError
import numpy as np

@config.node
class CastSynapticLocation(AfterConnectivityHook):
    """
    Replaces pre/post connection locations on a connectivity set with the relayed location
    """

    # connections: list["ConnectivityReference"] = config.reflist(refs.connectivity_ref, required=True)
    connections: list[str] = config.list(required=True)
    suffix_name: str = config.attr(required=True)
    # new_locations:

    def postprocess(self):

        for connection in set(self.connections):
            try:
                cs = self.scaffold.get_connectivity_set(connection)
            except DatasetNotFoundError:
                raise ConnectivityError(
                    f"AfterConnectivityHook {self.suffix_name} do not find {connection} ConnectivitySet."
                )
            except ValueError as e:
                raise e

            pre, post = cs.load_connections().all()
            pre_ps = cs.pre_type.get_placement_set()
            post_ps = cs.post_type.get_placement_set()

            # if change pre
            morph = pre_ps.load_morphologies().get(0)
            idx_b = []
            idx_loc = []
            for i_b, b_i in enumerate(morph.branches):
                if b_i.contains_labels(["axon_initial_segment", "AIS"]):
                    print("branch with ais ", i_b)
                    idx_b.append(i_b)
                    idx_loc.append(
                        np.where(b_i.get_label_mask(["axon_initial_segment", "AIS"]))[0]
                    )

            if len(idx_b) > 0:
                idx_0 = idx_b[0]  # first branch with ais
                idx_loc_0 = idx_loc[0][len(idx_loc[0]) // 2]  # middle loc of first branch
            else:
                idx_0 = 0  # soma
                idx_loc_0 = 0

            pre[:, 1] = idx_0
            pre[:, 2] = idx_loc_0

            new_name_cs = connection + "_" + self.suffix_name
            self.scaffold.connect_cells(pre_ps, post_ps, pre, post, new_name_cs)