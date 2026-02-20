import random

import numpy as np
from bsb import config
from bsb.config import types
from bsb.exceptions import *
from bsb.reporting import report, warn
from bsb.simulation.results import SimulationRecorder
from bsb.simulation.targetting import LocationTargetting
from bsb_neuron.device import NeuronDevice
from patch import p


@config.node
class MossyFibersSpikeGenerator(NeuronDevice, classmap_entry="mossy_fibers_spike_generator"):
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


    def implement(self, adapter, simulation, simdata):

        # Retrieve populations of Granule and Golgi cells
        for post_cm, post_pop in simdata.populations.items():
            if post_cm.cell_type.name == "granule_cell":
                pop_gr = post_pop
            if post_cm.cell_type.name == "golgi_cell":
                pop_go = post_pop

        # generate a random spike time pattern for every mossy fiber
        mf_ids = simulation.scaffold.cell_types.get("mossy_fibers").get_placement_set().load_ids()
        rng = np.random.default_rng(seed=self.seed)
        pattern = (
            np.cumsum(
                (1 - self.noise) * self.interval
                + rng.exponential(self.noise * self.interval, [len(mf_ids), self.number]),
                axis=1,
            )
            - (1 - self.noise) * self.interval
            + self.start
        )
        print("size pattern:", pattern.shape, " first: ", pattern[0])
        orig_gr = (
            simulation.cell_models.get("granule_cell").get_placement_set().load_positions()
        )  # all cell
        orig_go = (
            simulation.cell_models.get("golgi_cell").get_placement_set().load_positions()
        )  # all cell

        syn_type_go = ["AMPA_MF", "NMDA"]
        # Here we assume that connectivity mossy to GoC or GrC are already fused
        cs_mgr = self.scaffold.get_connectivity_set("mossy_fibers_to_granule_cell")
        mossy_list, granule_list = cs_mgr.load_connections().incoming().to(simdata.chunks).all()
        cs_mgc = self.scaffold.get_connectivity_set("mossy_fibers_to_golgi_cell")
        mossy_list_go, golgi_list = cs_mgc.load_connections().incoming().to(simdata.chunks).all()
        for model, mf_ids in self.targetting.get_targets(
            adapter, simulation, simdata
        ).items():  # targetting done on mf id: pop is the list of id
            print("stim targetting ", model, " id mf ", mf_ids)

            find_mossy_gr = [ mossy[0] in mf_ids for mossy in mossy_list]
            find_mossy_go = [ mossy[0] in mf_ids for mossy in mossy_list_go]
            post_mfGr = granule_list[find_mossy_gr]
            post_mfGo = golgi_list[find_mossy_go]

            # Insert and stimulate synapses in granule cells # important adding syn even if not stimulated in KO model, [glu] is a syn param
            for i, gr_i in enumerate(post_mfGr):
                for (
                    syn_type
                ) in (
                    self.synapses
                ):  # check if already there is the syn (if 2 or more stimulator each one add a syn)
                    for syn in pop_gr[gr_i[0]].get_location(gr_i[1:]).section.synapses:
                        if syn.synapse_name == syn_type:
                            break
                    else:
                        syn = pop_gr[gr_i[0]].insert_synapse(syn_type, gr_i[1:])
                    if mossy_list[i, 0] in mf_ids:
                        print(
                            "pre id ",
                            mossy_list[i],
                            " id ",
                            pop_gr[gr_i[0]].id,
                            " ",
                            gr_i,
                            " pos ",
                            orig_gr[pop_gr[gr_i[0]].id],
                        )
                        syn.stimulate(
                            pattern=pattern[mossy_list[i, 0], pattern[mossy_list[i, 0]] < self.end],
                            weight=self.weight,
                            delay=self.delay,
                        )  # select pattern of pre_mfGr[i,0]
                        print("Inserted synapses in granule cells")

            # Insert and stimulate synapses in golgi cells
            for i, go_i in enumerate(post_mfGo):
                if mossy_list_go[i, 0] in mf_ids:
                    print(
                        "pre id ",
                        mossy_list_go[i],
                        " go id ",
                        pop_go[go_i[0]].id,
                        " ",
                        go_i,
                        " pos ",
                        orig_go[pop_go[go_i[0]].id],
                    )
                    for syn_type in syn_type_go:
                        for syn in (
                            pop_go[go_i[0]].get_location(go_i[1:]).section.synapses
                        ):  # check if already there
                            if syn.synapse_name == syn_type:
                                break
                        else:
                            syn = pop_go[go_i[0]].insert_synapse(syn_type, go_i[1:])
                        syn.stimulate(
                            pattern=pattern[mossy_list_go[i, 0], pattern[mossy_list_go[i, 0]] < self.end],
                            weight=self.weight,
                            delay=self.delay,
                        )  # select pattern
                        print("Inserted synapses in golgi cells")
