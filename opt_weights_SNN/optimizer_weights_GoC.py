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

from cerebellar_models.optimization.features import *
from cerebellar_models.optimization.fitness import *
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
    bin_ms = 5.0
    fs = 1000.0 / bin_ms  # Hz
    theta_lo, theta_hi = 6.0, 12.0
    total_lo, total_hi = 1.0, 50.0

    for nb, spikes, ct_raw in zip(nb_neurons, spiketrains, ct_names):
        ct = ct_raw.split("_")[0] + "_" + ct_raw.split("_")[1] if ct_raw.count("_") >= 1 else ct_raw
        unique_counts = np.unique(spikes.array_annotations["senders"], return_counts=True)[1]
        unique_counts = np.concatenate([unique_counts, np.zeros(int(nb) - len(unique_counts))])

        duration_s = (spikes.t_stop - spikes.t_start) / 1000.0
        duration_s_val = (
            float(duration_s.magnitude) if hasattr(duration_s, "magnitude") else float(duration_s)
        )
        if duration_s_val <= 0.0:
            all_fr = np.zeros(int(nb), dtype=float)
        else:
            all_fr = (unique_counts / duration_s).magnitude

        mean_fr = np.mean(all_fr) if all_fr.size > 0 else 0.0
        std_fr = np.std(all_fr) if all_fr.size > 0 else 0.0
        if ct not in target_features:
            idx = order_ct.index(ct) if ct in order_ct else 0
            errors[idx * 2] = 1e3
            errors[idx * 2 + 1] = 1e3
            continue

        tgt = target_features[ct]
        target_mean = tgt["mean_fr"]
        tol_mean = tgt["tol_mean_fr"]
        target_std = tgt["std"]
        tol_std = tgt["std_tol"]

        z_tol = (mean_fr - target_mean) / tol_mean if tol_mean != 0 else 0.0
        mean_fr_error = 0.0
        penalty = 1.0  # extra penalty when violating the "undesired side" of target

        if ct == "golgi_cell":
            # Penalize above-target more strongly (avoid too high GoC firing).
            if z_tol > 0.0:
                mean_fr_error = z_tol + penalty
            else:
                z = abs(z_tol)
                mean_fr_error = 0.0 if z <= 1.0 else z

        elif ct == "granule_cell":
            # Penalize below-target more strongly (avoid too silent GrC).
            if z_tol >= 0.0:
                mean_fr_error = 0.0 if z_tol <= 1.0 else z_tol
            else:
                mean_fr_error = abs(z_tol) + penalty

        z_std = (std_fr - target_std) / tol_std if tol_std != 0 else 0.0
        std_error = max(0.0, abs(z_std) - 1.0)

        if ct == "golgi_cell":
            min_fr = 7.0
            max_silent_frac = 0.10  # at most 10% below min_fr

            silent_frac = np.mean(all_fr < min_fr) if all_fr.size else 1.0
            if silent_frac > max_silent_frac:
                silence_error = (silent_frac - max_silent_frac) * 2.0
                mean_fr_error += silence_error

            # Guardrail for a high-firing tail
            upper_fr = 20.0
            high_fr = all_fr[all_fr > upper_fr]
            if len(high_fr) > 0:
                high_error = np.mean((high_fr - upper_fr) / upper_fr)
                mean_fr_error += high_error

        theta_sync_error = 0.0
        if ct == "golgi_cell":
            # Build population spike-count time series by binning all spike times of merged train.
            t0 = float(spikes.t_start.magnitude)  # ms
            t1 = float(spikes.t_stop.magnitude)  # ms
            times_ms = spikes.times.magnitude

            # Need enough bins for a meaningful spectrum estimate
            # (rule of thumb: >= ~64 points)
            if (t1 > t0 + 64.0 * bin_ms) and (times_ms.size > 0):
                edges = np.arange(t0, t1 + bin_ms, bin_ms)
                counts, _ = np.histogram(times_ms, bins=edges)
                x = counts.astype(float)
                x = x - np.mean(x)

                n = len(x)
                if n >= 64:
                    # One-sided FFT power
                    X = np.fft.rfft(x)
                    P = (np.abs(X) ** 2) / n
                    freqs = np.fft.rfftfreq(n, d=1.0 / fs)

                    # Restrict to low frequency band where we expect GoC network rhythms
                    valid = (freqs >= total_lo) & (freqs <= total_hi)
                    freqs_v = freqs[valid]
                    P_v = P[valid]

                    band = (freqs_v >= theta_lo) & (freqs_v <= theta_hi)

                    total_power = np.trapz(P_v, freqs_v) if freqs_v.size > 1 else 0.0
                    theta_power = (
                        np.trapz(P_v[band], freqs_v[band])
                        if np.any(band) and freqs_v.size > 1
                        else 0.0
                    )

                    if total_power > 0.0:
                        theta_ratio = float(theta_power / total_power)

                        # Make it a "hinge" error: we want theta_ratio >= target - tol.
                        # Start soft. Tune later after you inspect baseline values.
                        target_theta_ratio = tgt.get("theta_ratio_target", 0.20)
                        tol_theta_ratio = tgt.get("theta_ratio_tol", 0.10)
                        theta_sync_error = max(
                            0.0, (target_theta_ratio - theta_ratio) / max(tol_theta_ratio, 1e-12)
                        )
                    else:
                        # No spectral power estimate -> treat as not meeting target
                        theta_sync_error = tgt.get("theta_ratio_target", 0.20) / max(
                            tgt.get("theta_ratio_tol", 0.10), 1e-12
                        )
                else:
                    theta_sync_error = tgt.get("theta_ratio_target", 0.20) / max(
                        tgt.get("theta_ratio_tol", 0.10), 1e-12
                    )
            else:
                theta_sync_error = tgt.get("theta_ratio_target", 0.20) / max(
                    tgt.get("theta_ratio_tol", 0.10), 1e-12
                )

        idx = order_ct.index(ct)

        if ct == "granule_cell":
            errors[idx * 2] = mean_fr_error
            errors[idx * 2 + 1] = std_error

        elif ct == "golgi_cell":
            errors[idx * 2] = mean_fr_error
            errors[idx * 2 + 1] = theta_sync_error

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


def evaluate_network(opt, ind, simulation_name="basal_activity", transient=300.0):
    opt.set_weights(ind, simulation_name)
    results = opt.run_simulation(simulation_name)
    spiketrains, nb_neurons, ct_names = extract_spikes(scaffold, [results], target_features)
    duration = scaffold.simulations["basal_activity"].duration
    spiketrains = get_filt_spikes(spiketrains, transient, duration)
    errors = compute_metrics(spiketrains, nb_neurons, target_features, ct_names, opt.order_ct)
    return errors


if __name__ == "__main__":
    print(f"Loading circuit for rank {rank}")
    scaffold = from_storage(f"granular_layer{rank}.hdf5", comm=comm.Split(rank))
    print("Done loading circuit")
    transient = 300.0
    order_ct = ["granule_cell", "golgi_cell"]
    bounds = {
        "ascending_axon_to_golgi": (0.01, 1.5),
        "glomerulus_to_golgi": (0.01, 1.5),
        # "golgi_to_glomerulus": (0.2, 1.3),
        "golgi_to_golgi": (0.01, 1.3),
        "parallel_fiber_to_golgi": (0.01, 1.5),
        "gap_goc": (0.01, 30.0),
        # "glomerulus_to_granule": (0.7, 1.2)
    }
    target_features = {
        "granule_cell": {"mean_fr": 0.81, "std": 1.3, "tol_mean_fr": 1.5, "std_tol": 1.0},
        "golgi_cell": {"mean_fr": 19, "std": 10, "tol_mean_fr": 10, "std_tol": 5.0},
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
        optimizer.POP_SIZE = 92
        optimizer.N_GEN = 150
        optimizer.init_params_fn = random_init
        optimizer.evaluate_fn = evaluate_network
        optimizer.print_evolution = True
        optimizer.return_pareto = True

        output_folder = "network_opt_GrC"
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
