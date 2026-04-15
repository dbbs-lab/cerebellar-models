import logging
import os
import random
import time
import warnings

import bsb.options
import numpy as np
from bsb.core import from_storage
from mpi4py import MPI
from quantities import ms

from cerebellar_models.optimization.optimizer_weights import (
    ArchiveND,
    Constraint,
    WeightsOptimizer,
    _apply_constraints,
)

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


def extract_ct_device_name(scaffold, device_name: str, simulation_name="basal_activity"):
    if "_record" in device_name:
        targetting = scaffold.simulations[simulation_name].devices[device_name].targetting
        ct = targetting.cell_models[0].name
        labels = targetting["labels"] if "labels" in targetting else set()
        return ct, labels
    else:
        return device_name, set()


def get_labelled_ct_name(ct_name, labels):
    extra = "_".join(labels)
    return ct_name + (extra if len(extra) == 0 else "_" + extra)


def extract_spikes(scaffold, results, target_features, simulation_name="basal_activity"):
    spikes_res = []
    cell_dict = {}
    current_id = 0

    for block in list(results):
        spiketrains = block.block.segments[0].spiketrains  # assume only one segment

        for st in spiketrains:
            st.segment = None  # remove spiketrain segment to allow merging
            cell_type, labels = extract_ct_device_name(
                scaffold, st.annotations["device"], simulation_name
            )
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
    errors = np.zeros(len(order_ct) * 2, dtype=float)

    for nb, spikes, ct in zip(nb_neurons, spiketrains, ct_names):
        # --- Firing rates per neuron ---
        unique_counts = np.unique(spikes.array_annotations["senders"], return_counts=True)[1]
        unique_counts = np.concatenate([unique_counts, np.zeros(nb - len(unique_counts))])
        duration_s = float((spikes.t_stop - spikes.t_start).magnitude) / 1000.0
        all_fr = (unique_counts / duration_s) if duration_s > 0.0 else np.zeros(nb)

        mean_fr = np.mean(all_fr) if all_fr.size > 0 else 0.0
        std_fr = np.std(all_fr) if all_fr.size > 0 else 0.0

        tgt = target_features[ct]
        target_mean = tgt["mean_fr"]
        tol_mean = tgt["tol_mean_fr"]
        target_std = tgt["std"]
        tol_std = tgt["std_tol"]

        z = (mean_fr - target_mean) / tol_mean
        if z >= 0.0:
            mean_fr_error = 0.0 if z <= 1.0 else z
        else:
            mean_fr_error = abs(z) + 1.0

        # --- Std FR error ---
        z_std = (std_fr - target_std) / tol_std if tol_std != 0 else 0.0
        std_error = max(0.0, abs(z_std) - 1.0)

        idx = order_ct.index(ct)
        errors[idx * 2] = mean_fr_error
        errors[idx * 2 + 1] = std_error

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


if __name__ == "__main__":
    print(f"Loading circuit for rank {rank}")
    scaffold = from_storage(f"mouse_cerebellum{rank}.hdf5", comm=comm.Split(rank))
    print("Done loading circuit")
    transient = 300.0
    order_ct = ["purkinje_cell"]
    bounds = {
        "ascending_axon_to_purkinje": (0.1, 1.5),
        "stellate_to_purkinje": (0.1, 1.5),
        "basket_to_purkinje": (0.1, 1.5),
        "parallel_fiber_to_purkinje": (0.2, 1.5),
    }
    target_features = {
        "purkinje_cell": {"mean_fr": 31.0, "std": 1.0, "tol_mean_fr": 10.0, "std_tol": 2.0},
    }
    if comm.Get_rank() == 0:
        fitness = {
            "Mean Firing Rate Purkinje": -1,
            "Standard Deviation Purkinje": -1,
        }
        archive = ArchiveND(eps_vec=[0.05] * 2, tau_vec=[0.15] * 2, cap=1000)

        model_params = {
            conn: scaffold.simulations["basal_activity"].connection_models[conn].synapses[0].weight
            for conn in bounds
        }

        optimizer = WeightsOptimizer(
            scaffold=scaffold,
            opt_params=list(bounds.keys()),
            model_params=model_params,
            fitness=fitness,
            bounds=bounds,
            comm=comm,
            archive=archive,
            knee_weights=[1.0, 1.0],
        )

        optimizer.order_ct = order_ct
        param_positions = {name: i for i, name in enumerate(optimizer.opt_params)}
        constraint = Constraint(
            func=constrain_bounds,
            name="bounds",
            ctx={"bounds": {param_positions[k]: v for k, v in bounds.items()}},
        )

        optimizer.add_constraint(constraint)
        optimizer.POP_SIZE = 92
        optimizer.N_GEN = 150
        optimizer.init_params_fn = random_init
        optimizer.evaluate_fn = None
        optimizer.print_evolution = True
        optimizer.return_pareto = True

        output_folder = "network_opt_BC"
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
                    scaffold.simulations["basal_activity"].connection_models[conn_name].synapses
                )
                for s in synapses:
                    s.weight = w
            results = scaffold.run_simulation("basal_activity")
            print(f"Simulation for {rank} complete!")
            spiketrains, nb_neurons, ct_names = extract_spikes(scaffold, [results], target_features)
            duration = scaffold.simulations["basal_activity"].duration
            spiketrains = get_filt_spikes(spiketrains, transient, duration)
            errors = compute_metrics(spiketrains, nb_neurons, target_features, ct_names, order_ct)
            comm.gather(errors, root=0)
            wait_sim = comm.recv(source=0)
