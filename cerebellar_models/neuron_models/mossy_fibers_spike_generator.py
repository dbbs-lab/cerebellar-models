from bsb import config
from bsb.config import types
from bsb_neuron.device import NeuronDevice
from bsb.simulation.results import SimulationRecorder
from bsb.simulation.targetting import LocationTargetting
from bsb.exceptions import *
from bsb.reporting import report, warn
import numpy as np
from patch import p
import random


@config.node
class MossyFibersSpikeGenerator(
    NeuronDevice, classmap_entry="mossy_fibers_spike_generator"
):
    locations = config.attr(type=LocationTargetting, default={"strategy": "soma"})
    synapses = config.attr(type=types.list())
    seed = config.attr(type=int, default=123)
    start = config.attr(type=float, default=0.0)
    end = config.attr(type=float, default=1e10)
    interval = config.attr(type=float, default=1.0)
    number = config.attr(type=int, default=1)
    noise = config.attr(type=float, default=0.0)
    weight = config.attr(type=float, default=0.1)
    delay = config.attr(type=float, default=0.0)
    # parameters = config.catch_all(type=types.any_())

    def implement(self, adapter, simulation, simdata):

        cs = self.scaffold.get_connectivity_set("glomerulus_to_granule")
        pre_glGr, post_gr = cs.load_connections().incoming().to(simdata.chunks).all()
        print("glom ", pre_glGr.shape, " to gr ", post_gr.shape)
        uniquePre_glGr, uniquePre_count_glGr = np.unique(
            pre_glGr[:, 0], axis=0, return_counts=True
        )
        print("[glom_i] to gr", np.mean(uniquePre_count_glGr))

        cs = self.scaffold.get_connectivity_set("glomerulus_to_golgi")
        pre_glGo, post_go = cs.load_connections().incoming().to(simdata.chunks).all()
        print("glom ", pre_glGo.shape, " to go ", post_go.shape)
        uniquePre_glGo, uniquePre_count_glGo = np.unique(
            pre_glGo[:, 0], axis=0, return_counts=True
        )
        print("[glom_i] to go", np.mean(uniquePre_count_glGo))

        cs = self.scaffold.get_connectivity_set("mossy_fibers_to_glomerulus")
        pre_mf, post_mfGl = cs.load_connections().all()
        print("mf ", pre_mf.shape, " to glom ", post_mfGl.shape)
        uniquePost, uniquePost_count = np.unique(
            post_mfGl[:, 0], axis=0, return_counts=True
        )
        print("[glom_i] from ", np.mean(uniquePost_count))

        ########################################
        ### interesect unique mfGl and glGr
        dtype = ", ".join([str(uniquePre_glGr.dtype)] * 1)
        _, idx_post, idx_glGr = np.intersect1d(
            uniquePost.view(dtype),
            uniquePre_glGr.view(dtype),
            assume_unique=True,
            return_indices=True,
        )

        # new pre/post matrix
        pre_mfGr = np.zeros(
            [np.dot(uniquePost_count[idx_post], uniquePre_count_glGr[idx_glGr]), 3],
            dtype=dtype,
        )
        post_mfGr = np.zeros(
            [np.dot(uniquePost_count[idx_post], uniquePre_count_glGr[idx_glGr]), 3],
            dtype=dtype,
        )

        # fill pre/post matrices
        j = 0
        for i, idx_post_i in enumerate(
            idx_post
        ):  # each index (from mf side) of common glom with Gr, same size if indeces in idx_glGr
            idx_glGr_i = idx_glGr[i]
            for idx_post_mfGl_i in np.where(post_mfGl[:, 0] == uniquePost[idx_post_i])[
                0
            ]:  # retreive mf connected to selected glom
                for idx_pre_glGr_i in np.where(
                    pre_glGr[:, 0] == uniquePre_glGr[idx_glGr_i]
                )[
                    0
                ]:  # retreive gr connected to selected glom
                    pre_mfGr[j] = pre_mf[idx_post_mfGl_i]
                    post_mfGr[j] = post_gr[idx_pre_glGr_i]
                    j += 1
        
        print('size pre/post cs matrices', pre_mfGr.shape, ' ', post_mfGr.shape, ' j:', j) 
        ########################################
        ### interesect unique mfGl and glGo
        dtype = ", ".join([str(uniquePre_glGo.dtype)] * 1)
        _, idx_post, idx_glGo = np.intersect1d(
            uniquePost.view(dtype),
            uniquePre_glGo.view(dtype),
            assume_unique=True,
            return_indices=True,
        )

        # new pre/post matrix
        pre_mfGo = np.zeros(
            [np.dot(uniquePost_count[idx_post], uniquePre_count_glGo[idx_glGo]), 3],
            dtype=dtype,
        )
        post_mfGo = np.zeros(
            [np.dot(uniquePost_count[idx_post], uniquePre_count_glGo[idx_glGo]), 3],
            dtype=dtype,
        )

        # fill matrices
        j = 0
        for i, idx_post_i in enumerate(idx_post):
            idx_glGo_i = idx_glGo[i]
            for idx_post_mfGl_i in np.where(post_mfGl[:, 0] == uniquePost[idx_post_i])[
                0
            ]:
                for idx_pre_glGo_i in np.where(
                    pre_glGo[:, 0] == uniquePre_glGo[idx_glGo_i]
                )[0]:
                    pre_mfGo[j] = pre_mf[idx_post_mfGl_i]
                    post_mfGo[j] = post_go[idx_pre_glGo_i]
                    j += 1
        
        print('size pre/post cs matrices', pre_mfGo.shape, ' ', post_mfGo.shape, ' j:', j)
        ########################################
        for post_cm, post_pop in simdata.populations.items():
            if post_cm.cell_type.name == "granule_cell":
                pop_gr = post_pop
            if post_cm.cell_type.name == "golgi_cell":
                pop_go = post_pop

        # spike time patter
        mf_ids = (
            simulation.scaffold.cell_types.get("mossy_fibers")
            .get_placement_set()
            .load_ids()
        )
        rng = np.random.default_rng(seed=self.seed)
        pattern = (
            np.cumsum(
                (1 - self.noise) * self.interval
                + rng.exponential(
                    self.noise * self.interval, [len(mf_ids), self.number]
                ),
                axis=1,
            )
            - (1 - self.noise) * self.interval
            + self.start
        )
        print('size pattern:', pattern.shape, ' first: ', pattern[0])
        
        syn_type_go = ["AMPA_MF", "NMDA"]
        for model, mf_ids in self.targetting.get_targets(
            adapter, simulation, simdata
        ).items():  # targetting done on mf id: pop is the list of id
            print("stim targetting ", model, " id mf ", mf_ids)
            orig_gr = (
                simulation.cell_models.get("granule_cell")
                .get_placement_set()
                .load_positions()
            )  # all cell
            orig_go = (
                simulation.cell_models.get("golgi_cell")
                .get_placement_set()
                .load_positions()
            )  # all cell
            # Insert and stimulate synapses in granule cells # important adding syn even if not stimulated in KO model, [glu] is a syn param
            for i, gr_i in enumerate(post_mfGr):
                for syn_type in self.synapses: #check if already there the syn (if 2 or more stimulator each one add a syn)
                    for syn in pop_gr[gr_i[0]].get_location(gr_i[1:]).section.synapses:
                        if syn.synapse_name == syn_type:
                            break
                    else:
                        syn = pop_gr[gr_i[0]].insert_synapse(syn_type, gr_i[1:])
                    if pre_mfGr[i, 0] in mf_ids:
                        print(
                            "pre id ",
                            pre_mfGr[i],
                            " id ",
                            pop_gr[gr_i[0]].id,
                            " ",
                            gr_i,
                            " pos ",
                            orig_gr[pop_gr[gr_i[0]].id],
                        )
                        syn.stimulate(
                            pattern=pattern[pre_mfGr[i, 0], pattern[pre_mfGr[i, 0]] < self.end],
                            weight=self.weight,
                            delay=self.delay,
                        )  # select pattern of pre_mfGr[i,0]
                        print("Inserted synapses in granule cells")

            # Insert and stimulate synapses in golgi cells
            for i, go_i in enumerate(post_mfGo):
                if pre_mfGo[i, 0] in mf_ids:
                    print(
                        "pre id ",
                        pre_mfGo[i],
                        " go id ",
                        pop_go[go_i[0]].id,
                        " ",
                        go_i,
                        " pos ",
                        orig_go[pop_go[go_i[0]].id],
                    )
                    for syn_type in syn_type_go: 
                        for syn in pop_go[go_i[0]].get_location(go_i[1:]).section.synapses: #check if already there
                            if syn.synapse_name == syn_type:
                                break
                        else:
                            syn = pop_go[go_i[0]].insert_synapse(syn_type, go_i[1:])
                        syn.stimulate(
                            pattern=pattern[pre_mfGo[i, 0], pattern[pre_mfGo[i, 0]] < self.end],
                            weight=self.weight,
                            delay=self.delay,
                        )  # select pattern
                        print("Inserted synapses in golgi cells")
