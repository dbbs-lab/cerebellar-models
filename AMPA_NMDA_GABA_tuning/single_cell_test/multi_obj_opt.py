import matplotlib.pyplot as plt
import nest
import numpy as np
import pandas as pd
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.problem import Problem
from pymoo.optimize import minimize
from pymoo.termination.default import DefaultMultiObjectiveTermination


def single_sim(cell_params, sim_params, E_rev, duration):
    syn_spec = {"weight": 1.0, "delay": 1, "receptor_type": sim_params["receptor_type"]}
    nest.ResetKernel()
    nest.Install("cerebmodule")
    nest.resolution = 0.025

    cell = nest.Create("eglif_multisyn", 1, params=cell_params)
    input_spikes = nest.Create(
        "spike_generator", params={"spike_times": [sim_params["spike_time"]]}
    )
    nest.Connect(input_spikes, cell, syn_spec=dict(syn_spec, synapse_model="static_synapse"))

    multimeter = nest.Create(
        "multimeter",
        params={
            "interval": 0.025,
            "record_from": sim_params["record_from"],
            "record_to": "memory",
            "label": sim_params["label_multimeter"],
        },
    )
    spike_recorder = nest.Create(
        "spike_recorder",
        params={"record_to": "memory", "label": sim_params["label_spike_recorder"]},
    )
    nest.Connect(multimeter, cell)
    nest.Connect(cell, spike_recorder)
    nest.Simulate(duration)

    multimeter_data = multimeter.get()["events"]
    times = multimeter_data["times"]
    I_syn = multimeter_data[sim_params["current_name"]]
    V_m = multimeter_data["V_m"]
    g_syn = I_syn / (E_rev - V_m)
    if "Mg_block" in sim_params["record_from"]:
        g_syn = g_syn / multimeter_data["Mg_block"]
    if "I_nmda1" in sim_params["record_from"] and "I_nmda2" in sim_params["record_from"]:
        kernel1 = multimeter_data["I_nmda1"] / ((E_rev - V_m) * multimeter_data["Mg_block"])
        kernel2 = multimeter_data["I_nmda2"] / ((E_rev - V_m) * multimeter_data["Mg_block"])
        return times, g_syn, kernel1, kernel2
    elif "I_ampa1" in sim_params["record_from"] and "I_ampa2" in sim_params["record_from"]:
        kernel1 = multimeter_data["I_ampa1"] / (E_rev - V_m)
        kernel2 = multimeter_data["I_ampa2"] / (E_rev - V_m)
        return times, g_syn, kernel1, kernel2
    elif "I_gaba1" in sim_params["record_from"] and "I_gaba2" in sim_params["record_from"]:
        kernel1 = multimeter_data["I_gaba1"] / (E_rev - V_m)
        kernel2 = multimeter_data["I_gaba2"] / (E_rev - V_m)
        return times, g_syn, kernel1, kernel2
    else:
        return times, g_syn


def extract_trace(file_path, E_rev, v_clamp=-40.0, mg_block=None):
    df = pd.read_csv(file_path, delim_whitespace=True)
    df.columns = ["Time", "Current"]
    df["g"] = -df["Current"].values / (E_rev - v_clamp) * 1000
    if mg_block is not None:
        df["g"] = df["g"] / mg_block
    return df["Time"], df["g"]


def optimization(
    problem, xtol=1e-8, cvtol=1e-8, ftol=25e-4, period=30, n_gen=100, n_max_evals=10000
):
    algorithm = NSGA2(pop_size=32)

    termination = DefaultMultiObjectiveTermination(
        xtol=xtol, cvtol=cvtol, ftol=ftol, period=period, n_max_gen=n_gen, n_max_evals=n_max_evals
    )
    results = minimize(problem, algorithm, termination, seed=1, verbose=True)
    return results


def compute_metrics(times, g_syn, trace_g, t_min=250, t_max=1000):
    mask_t = np.logical_and(times <= t_max, times > t_min)
    trace_g_masked = trace_g[mask_t]
    g_masked = g_syn[mask_t]
    mse_val = np.mean((trace_g_masked - g_masked) ** 2)
    peak_penalty = abs(np.max(trace_g_masked) - np.max(g_masked)) / np.max(trace_g_masked)
    return mse_val, peak_penalty


def plot_pareto_front(res, sol_idx, threshold, title, save_path=None):
    mse_vals = res.F[:, 0]
    peak_penalty_vals = res.F[:, 1]
    plt.figure(figsize=(6, 4))
    plt.title(title)
    plt.scatter(mse_vals, peak_penalty_vals * 100, c="blue", alpha=0.5, edgecolors="none")
    plt.scatter(
        mse_vals[sol_idx],
        peak_penalty_vals[sol_idx] * 100,
        c="red",
        label="Selected solution",
        zorder=5,
        edgecolors="none",
    )
    plt.axhline(threshold, c="red", ls="--", alpha=0.5, label="Peak penalty threshold")
    plt.xlabel("MSE")
    plt.ylabel("Peak penalty (%)")
    plt.legend()
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path, dpi=300)
    plt.show()


def select_solution(F, threshold=0.05):
    mask = F[:, 1] <= threshold
    filtered = F[mask]
    idx_in_filtered = np.argmin(filtered[:, 0])
    idx = np.where(mask)[0][idx_in_filtered]
    return idx, F[idx]


def save_params(opt_params, output_filename):
    with open(output_filename, "w") as f:
        for k, v in zip(opt_params.keys(), opt_params.values()):
            f.write(f"{k}: {v}\n")


def plot_conductance(
    params,
    sim_params,
    trace_file,
    E_rev,
    duration,
    optimized_params=None,
    title="Synaptic conductance",
    save_path=None,
    v_clamp=-40.0,
    mg_block=None,
    xlim=None,
):
    for key, val in optimized_params.items():
        params[key] = val

    times_sim, g_sim, kernel1, kernel2 = single_sim(params, sim_params, E_rev, duration)
    trace_t, trace_g = extract_trace(trace_file, E_rev, v_clamp, mg_block)
    plt.figure(figsize=(8, 5))
    plt.title(title)
    plt.plot(times_sim, g_sim, label="NEST trace", color="black")
    plt.plot(trace_t, trace_g, label="NEURON trace", color="gray", linestyle="dashed")
    plt.plot(times_sim, kernel1, label="NEST kernel (fast)", color="red", alpha=0.3)
    plt.plot(times_sim, kernel2, label="NEST kernel (slow)", color="green", alpha=0.3)
    plt.xlabel("Time [ms]")
    plt.ylabel(r"$g_{syn}$ [ns]")
    plt.legend()
    if xlim is not None:
        plt.xlim(xlim)
    if save_path is not None:
        plt.savefig(save_path, dpi=300)
    plt.show()


class SynapseOptimization(Problem):
    def __init__(
        self,
        cell_params,
        sim_params,
        trace_file,
        param_names,
        E_rev,
        tmin,
        tmax,
        v_clamp=-40.0,
        mg_block=None,
        duration=1001,
    ):
        self.cell_params = cell_params
        self.sim_params = sim_params
        self.trace_file = trace_file
        self.param_names = param_names
        self.E_rev = E_rev
        self.v_clamp = v_clamp
        self.mg_block = mg_block
        self.duration = duration
        self.tmin = tmin
        self.tmax = tmax
        super().__init__(
            n_var=len(param_names),
            n_obj=2,
            n_constr=1,
            xl=np.array([0, 0, 0, 0, 0, 0, 0]),
            xu=np.array([5, 5.0, 5.0, 10.0, 5, 5, 5]),
        )
        self.trace_t, self.trace_g = extract_trace(trace_file, E_rev, v_clamp, mg_block)

    def _evaluate(self, X, out, *args, **kwargs):
        f1 = []
        f2 = []
        for params in X:
            for name, value in zip(self.param_names, params):
                self.cell_params[name] = value
            times, g_syn = single_sim(
                self.cell_params, self.sim_params, self.E_rev, duration=self.duration
            )
            mse_val, peak_penalty = compute_metrics(
                times, g_syn, self.trace_g, self.tmin, self.tmax
            )
            f1.append(mse_val)
            f2.append(peak_penalty)
        out["F"] = np.column_stack([f1, f2])
        g = X[:, 2] - X[:, 3]
        out["G"] = g.reshape(-1, 1)
