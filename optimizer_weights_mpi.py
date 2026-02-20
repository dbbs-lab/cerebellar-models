import logging
import os
import random
import warnings
import numpy as np
from bsb.core import from_storage
import bsb.options
from cerebellar_models.optimization.features import *
from cerebellar_models.optimization.fitness import *
from cerebellar_models.optimization.optimizer_weights import (
    ArchiveND,
    Constraint,
    WeightsOptimizer,
    _apply_constraints,
)
from mpi4py import MPI
from quantities import ms
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
bsb.options.verbosity = 1

warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
logging.getLogger().setLevel(logging.ERROR)
for name in ("deap", "matplotlib", "numexpr"):
    logging.getLogger(name).setLevel(logging.ERROR)
np.seterr(all="ignore")
os.environ.setdefault("PYTHONWARNINGS", "ignore")


def get_filt_spikes(all_spikes, time_from, time_to):
    return [st.time_slice(time_from * ms, time_to * ms) for st in all_spikes]

def extract_ct_device_name(scaffold ,device_name: str, simulation_name = 'basal_activity'):
    if "_record" in device_name:
        targetting = (
            scaffold.simulations[simulation_name].devices[device_name].targetting
        )
        ct = targetting.cell_models[0].name
        labels = targetting["labels"] if "labels" in targetting else set()
        return ct, labels
    else:
        return device_name, set()

def get_labelled_ct_name(ct_name, labels):
    extra = "_".join(labels)
    return ct_name + (extra if len(extra) == 0 else "_" + extra)

def extract_spikes(scaffold, results, target_features, simulation_name='basal_activity'):
    spikes_res = []
    cell_dict = {}
    current_id = 0

    for block in list(results):
        spiketrains = block.block.segments[0].spiketrains  # assume only one segment

        for st in spiketrains:
            st.segment = None  # remove spiketrain segment to allow merging
            cell_type, labels = extract_ct_device_name(scaffold,st.annotations["device"], simulation_name)
            if cell_type in list(target_features.keys()):
                cell_type_label = get_labelled_ct_name(cell_type, labels)
                if cell_type_label not in cell_dict:
                    cell_dict[cell_type_label] = current_id
                    current_id += 1
                    spikes_res.append([])
                spikes_res[cell_dict[cell_type_label]].append(st)
    all_spikes = []
    nb_neurons = np.zeros(len(cell_dict), dtype=int)
    for i, cell_type in enumerate(cell_dict):
        sts = spikes_res[cell_dict[cell_type]]
        merged = sts[0]
        for st in sts[1:]:
            merged = merged.merge(st)
        all_spikes.append(merged)
        nb_neurons[i] = all_spikes[i].annotations["pop_size"]
    return all_spikes, nb_neurons, list(cell_dict.keys())

def compute_metrics(spiketrains, nb_neurons, target_features, ct_names, order_ct):
    errors = np.zeros(len(ct_names) * 2)
    for nb , spikes, ct in zip(nb_neurons , spiketrains, ct_names):
        unique_counts = np.unique(spikes.array_annotations["senders"], return_counts=True)[1]
        unique_counts = np.concatenate(
            [unique_counts, np.zeros(nb - len(unique_counts))]
        )

        all_fr = unique_counts / ((spikes.t_stop - spikes.t_start) / 1000.0)
        all_fr = all_fr.magnitude
        mean_fr = np.absolute(np.mean(all_fr) - target_features[ct]["mean_fr"]) / target_features[ct]["mean_fr"]
        std_fr = np.absolute(np.std(all_fr) - target_features[ct]["std"]) / target_features[ct]["std"]
        errors[order_ct.index(ct)*2] = mean_fr
        errors[order_ct.index(ct)*2 +1] = std_fr
    return errors

def constrain_bounds(ind, ctx):
    bounds = ctx["bounds"]
    for pos, (a, b) in bounds.items():
        x = float(ind[pos]) if np.isfinite(ind[pos]) else (a + b) / 2.0
        if a > b:
            a, b = b, a
        ind[pos] = min(max(x, a), b)
    return ind


def random_init(optimizer):
    ind = [random.uniform(*optimizer.bounds[p]) for p in optimizer.opt_params]
    ind = _apply_constraints(ind, optimizer.constraints)
    return ind


def evaluate_network(opt, ind, simulation_name='basal_activity', transient = 300.):
    opt.set_weights(ind, simulation_name)
    results = opt.run_simulation(simulation_name)
    spiketrains,  nb_neurons, ct_names = extract_spikes(scaffold, [results], target_features)
    duration = scaffold.simulations['basal_activity'].duration
    spiketrains = get_filt_spikes(spiketrains, transient, duration)
    errors = compute_metrics(spiketrains, nb_neurons,target_features, ct_names, opt.order_ct)
    return errors


if __name__ == "__main__":
    print(f"Loading circuit for rank {rank}")
    scaffold = from_storage(f'granular_layer{rank}.hdf5', comm=comm.Split(rank))
    print("Done loading circuit")
    transient = 300.
    order_ct = ['granule_cell', 'golgi_cell']
    bounds = {
        "ascending_axon_to_golgi": (0, 1.3),
        "glomerulus_to_golgi": (0, 1.3),
        "golgi_to_glomerulus": (0, 1.3),
        "golgi_to_golgi": (0, 1.3),
        "parallel_fiber_to_golgi": (0, 1.3),
        "gap_goc":  (0, 1.3)
    }
    target_features = {
        'golgi_cell': {'mean_fr': 19, 'std': 15},
        'granule_cell': {'mean_fr': 0.81, 'std': 1.3},
    }
    if comm.Get_rank() == 0:
        fitness = {
            "Mean Firing Rate Granule": -1,
            "Standard Deviation Granule": -1,
            "Mean Firing Rate Golgi": -1,
            "Standard Deviation Golgi": -1,
        }
        archive = ArchiveND(eps_vec=[0.05] * 4, tau_vec=[0.15] * 4, cap=1000)

        model_params = {
            conn: scaffold.simulations["basal_activity"]
            .connection_models[conn]
            .synapses[0].weight
            for conn in bounds
        }

        optimizer = WeightsOptimizer(
            scaffold=scaffold,
            opt_params=list(bounds.keys()),
            model_params=model_params,
            fitness=fitness,
            bounds=bounds,
            comm = comm,
            archive=archive,
            knee_weights=[1.0, 1.0, 1.0, 1.0],
        )

        optimizer.order_ct = order_ct
        param_positions = {name: i for i, name in enumerate(optimizer.opt_params)}
        constraint = Constraint(
            func=constrain_bounds,
            name="bounds",
            ctx={"bounds": {param_positions[k]: v for k, v in bounds.items()}},
        )

        optimizer.add_constraint(constraint)
        optimizer.POP_SIZE = 28
        optimizer.N_GEN = 10
        optimizer.init_params_fn = random_init
        optimizer.evaluate_fn = evaluate_network
        optimizer.print_evolution = True

        output_folder = "network_opt"
        os.makedirs(output_folder, exist_ok=True)
        optimizer.output_path = output_folder

        toolbox = optimizer.set_optimizer()
        best, fit = optimizer.optimize()
        print("Optimization complete!")
        print("Best individual:", best)
        print("Fitness:", fit)
        wait_sim = comm.bcast(False, root=0)
    else:
        wait_sim = comm.recv(source=0)
        opt_params = list(bounds.keys())
        while wait_sim:
            weights = comm.recv(source=0)
            for i, conn_name in enumerate(opt_params):
                w = float(weights[i])
                synapses = (
                    scaffold
                    .simulations["basal_activity"]
                    .connection_models[conn_name]
                    .synapses
                )
                for s in synapses:
                    s.weight = w
            results = scaffold.run_simulation("basal_activity")
            print(f'Simulation for {rank} complete!')
            spiketrains,  nb_neurons, ct_names = extract_spikes(scaffold, [results], target_features)
            duration = scaffold.simulations['basal_activity'].duration
            spiketrains = get_filt_spikes(spiketrains, transient, duration)
            errors = compute_metrics(spiketrains, nb_neurons,target_features, ct_names, order_ct)
            comm.gather(errors, root=0)
            wait_sim = comm.recv(source=0)
