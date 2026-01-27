import os
import json
import matplotlib.pyplot as plt
import pickle
from cerebellar_models.optimization.features import (
    multicomp_features,
    point_neuron_features,
)
from cerebellar_models.optimization.fitness import *
from cerebellar_models.optimization.fitness import _scalar
from cerebellar_models.optimization.utils import _suppress_output


dict_colors = {
    'BC': [1, 0.647, 0, 1.0],
    'SC': [1, 0.84, 0, 1.0],
    'PC': [0.275, 0.800, 0.275, 1.0],
    'GoC': [0, 0.45, 0.7, 1.0],
    'GrC': [0.7, 0.15, 0.15, 1.]
}

def save_outputs(cell_name, cell_params_best, error_dict, out_dir="./results_opt"):
    os.makedirs(out_dir, exist_ok=True)
    with open(f"{out_dir}/{cell_name}_opt.json", "w") as f:
        json.dump(cell_params_best, f, indent=4)
    with open(f"{out_dir}/{cell_name}_opt_err.json", "w") as f:
        json.dump(error_dict, f, indent=4)

def load_best_record(pkl_path: str) -> dict:
    records = []
    with open(pkl_path, "rb") as f:
        while True:
            try:
                records.append(pickle.load(f))
            except EOFError:
                break
    if len(records) == 0:
        raise RuntimeError(f"Empty best_history file: {pkl_path}")
    return records[-1]

def load_cell_params_from_json(cell_name: str, in_dir: str = "./results_opt") -> dict:
    path = os.path.join(in_dir, f"{cell_name}_opt.json")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing params json: {path}")
    with open(path, "r") as f:
        params = json.load(f)
    return params

def extract_multicomp_features(multicomp_data, protocol):
    return multicomp_features(
        multicomp_data,
        threshold=protocol["threshold"],
        start_stim=protocol["start_stim"],
        end_stim=protocol["end_stim"],
    )

def nest_single_sim(
    current, cell_params, protocol, nest_model="eglif_multirec_opt", multimeter=False
):
    with _suppress_output():
        import nest

        try:
            nest.set_verbosity("M_FATAL")
        except:
            try:
                nest.SetKernelStatus({"print_time": False})
            except:
                pass
        nest.ResetKernel()
        nest.SetKernelStatus({"local_num_threads": 1})
        try:
            nest.Install("cerebmodule")
        except:
            pass
        nest.resolution = 0.1

    import nest

    cell = nest.Create(nest_model, 1, params=cell_params)
    cell.V_m = cell_params["E_L"]
    gen = nest.Create("dc_generator")
    nest.SetStatus(
        gen, {"amplitude": current, "start": protocol["start_stim"], "stop": protocol["end_stim"]}
    )
    sr = nest.Create("spike_recorder")

    mm = None
    if multimeter:
        mm = nest.Create(
            "multimeter", params={"interval": 0.1, "record_from": ["V_m"], "record_to": "memory"}
        )
        nest.Connect(mm, cell)

    nest.Connect(gen, cell)
    nest.Connect(cell, sr)
    nest.Simulate(protocol["duration"])

    spikes = nest.GetStatus(sr, "events")[0]["times"]
    if multimeter:
        vm_data = nest.GetStatus(mm, "events")[0]
        t = np.array(vm_data["times"])
        Vm = np.array(vm_data["V_m"])
        return {"spikes": spikes, "vm": (t, Vm)}
    return spikes

def nest_protocol(multicomp_features, cell_params, protocol, multimeter=False):
    return {
        I: nest_single_sim(I, cell_params, protocol, multimeter=multimeter)
        for I in multicomp_features["current"].values
    }

def extract_nest_features(multicomp_features, cell_params, protocol, multimeter=False):
    nest_results = nest_protocol(multicomp_features, cell_params, protocol, multimeter)
    spike_trains = [np.array(nest_results[I]) for I in multicomp_features["current"].values]
    nest_features = point_neuron_features(
        spike_trains=spike_trains,
        currents=multicomp_features["current"].values,
        start_stim=protocol["start_stim"],
        end_stim=protocol["end_stim"],
    )
    return nest_results, nest_features

def plot_results(
    multicomp_feat, results, cell_params, start_stim, stop_stim, cell_name, output_dir="."
):
    import matplotlib.ticker as mticker
    from matplotlib.gridspec import GridSpec

    os.makedirs(output_dir, exist_ok=True)

    def add_fake_spikes_to_vm(t, v, spikes, v_peak=20.0, sigma_ms=0.01):
        t = np.asarray(t, dtype=float)
        v = np.asarray(v, dtype=float)
        spikes = np.asarray(spikes if spikes is not None else [], dtype=float)
        v_out = v.copy()
        if spikes.size == 0:
            return v_out
        for ts in spikes:
            idx = int(np.argmin(np.abs(t - ts)))
            amp = max(0.0, v_peak - v[idx])
            if amp == 0.0:
                continue
            v_out += amp * np.exp(-0.5 * ((t - ts) / sigma_ms) ** 2)
        return np.minimum(v_out, v_peak)

    def style_vm_axis(ax):
        ax.tick_params(axis="both", which="major", labelsize=10, length=4, width=1)
        ax.tick_params(axis="both", which="minor", labelsize=10, length=2, width=0.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(1)
        ax.spines["bottom"].set_linewidth(1)

    def style_I_axis(axI, show_xlabel=False):
        axI.spines["top"].set_visible(False)
        axI.spines["left"].set_visible(False)
        axI.spines["bottom"].set_visible(True)
        axI.spines["bottom"].set_linewidth(1)
        axI.spines["right"].set_visible(True)
        axI.spines["right"].set_linewidth(1)
        axI.spines["right"].set_color("#7a1f2b")

        axI.yaxis.set_ticks_position("right")
        axI.yaxis.set_label_position("right")
        axI.tick_params(
            axis="y", which="major", labelsize=7, length=3, width=1, colors="#7a1f2b", pad=3
        )

        axI.tick_params(axis="x", which="major", labelsize=9, length=3, width=1, pad=10)
        axI.tick_params(axis="x", which="minor", length=2, width=0.8, pad=10)

        if show_xlabel:
            axI.set_xlabel("Time [ms]", fontsize=10, labelpad=6)
        else:
            axI.set_xlabel("")

    any_data = next(iter(results.values()))
    t_any = np.asarray(any_data["vm"][0], dtype=float)
    t_min, t_max = float(np.min(t_any)), float(np.max(t_any))

    V_th = float(cell_params["V_th"])
    v_peak = 20.0

    ordered_currents = list(multicomp_feat["current"].values.astype(float))
    Imax = max(1.0, max(abs(float(c)) for c in ordered_currents))
    Iylim = 1.10 * Imax
    Iticks = [-Imax, 0.0, Imax]

    major_locator = mticker.MaxNLocator(nbins=5)
    minor_locator = mticker.AutoMinorLocator(2)

    n = len(ordered_currents)
    fig_h = max(3.0, 2.55 * n)
    fig = plt.figure(figsize=(7.6, fig_h))

    gs = GridSpec(
        nrows=2 * n,
        ncols=1,
        height_ratios=sum(([4.0, 1.6] for _ in range(n)), []),
        hspace=0.42,
    )

    for i, Iamp in enumerate(ordered_currents):
        data = results[Iamp]
        axV = fig.add_subplot(gs[2 * i, 0])
        axI = fig.add_subplot(gs[2 * i + 1, 0], sharex=axV)

        axV.set_zorder(0)
        axV.patch.set_alpha(0.0)
        axI.set_zorder(2)

        t, v = data["vm"]
        spikes = np.asarray(data["spikes"], dtype=float)

        if Iamp!=0:
            axV.axvspan(
                start_stim, stop_stim, facecolor="#ffe680", alpha=0.18, edgecolor="none", zorder=-10
            )

        v_plot = add_fake_spikes_to_vm(t, v, spikes, v_peak=v_peak, sigma_ms=0.01)
        axV.plot(t, v_plot, linewidth=0.8, color=dict_colors[cell_name], zorder=2)
        axV.axhline(V_th, color="gray", linestyle="--", linewidth=0.8, zorder=1)

        axV.set_ylabel(r"$V_m$ [mV]", fontsize=11)
        axV.set_ylim(-150, 25)
        axV.set_xlim(t_min, t_max)
        style_vm_axis(axV)
        axV.tick_params(axis="x", labelbottom=False)

        # current step
        t = np.asarray(t, dtype=float)
        I = np.zeros_like(t, dtype=float)
        I[(t >= start_stim) & (t <= stop_stim)] = float(Iamp)

        axI.plot(t, I, drawstyle="steps-post", linewidth=1., color="#7a1f2b")
        axI.set_ylim(-Iylim, Iylim)
        axI.set_yticks(Iticks)
        axI.set_ylabel("I [pA]", color="#7a1f2b", fontsize=6, labelpad=3)

        axI.xaxis.set_major_locator(major_locator)
        axI.xaxis.set_minor_locator(minor_locator)

        show_xlabel = i == n - 1
        style_I_axis(axI, show_xlabel=show_xlabel)

        for lab in axI.get_xticklabels(which="both"):
            lab.set_clip_on(False)
            lab.set_zorder(10)

        txt = axI.text(
            0.98,
            0.30 if float(Iamp) >= 0 else 0.70,
            f"I = {float(Iamp):g} pA",
            transform=axI.transAxes,
            ha="right",
            va="center",
            fontsize=8,
            color="#7a1f2b",
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.80),
            zorder=10,
        )
        txt.set_clip_on(True)
        txt.set_clip_path(axI.patch)

    fig.subplots_adjust(top=0.98, bottom=0.06)
    fig.savefig(f"{output_dir}/best_exp.png", dpi=700, bbox_inches="tight")
    fig.savefig(f"{output_dir}/best_exp.pdf", dpi=700, bbox_inches="tight")
    plt.close(fig)

    # -------- f-I curve --------
    nest_freqs = []
    for Iamp in multicomp_feat["current"].values:
        spikes = np.array(results[Iamp]["spikes"], dtype=float)
        spikes_stim = spikes[(spikes >= start_stim) & (spikes <= stop_stim)]
        freq = len(spikes_stim) / ((stop_stim - start_stim) / 1000.0)
        nest_freqs.append(freq)

    multicomp_freqs = multicomp_feat["mean_frequency"].apply(_scalar).fillna(0).values
    currents = multicomp_feat["current"].values.astype(float)

    fig, ax = plt.subplots(figsize=(6.3, 4.9))
    ax.plot(
        currents,
        multicomp_freqs,
        "o",
        label="NEURON (target)",
        color="gray",
        markersize=5,
        linestyle="--",
        linewidth=1.0,
    )
    ax.plot(
        currents,
        nest_freqs,
        "s-",
        label="NEST (best)",
        color=dict_colors[cell_name],
        markersize=5,
        linewidth=1.3,
    )

    ax.set_xlabel("Injected current [pA]", fontsize=11)
    ax.set_ylabel("Firing rate [Hz]", fontsize=11)
    ax.set_title(f"{cell_name} f-I curve", fontsize=12)
    ax.legend(fontsize=10, frameon=False)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_linewidth(1)

    ax.tick_params(axis="x", which="major", labelsize=10, length=4, width=1)
    ax.set_ylim(-2, 200)
    # ax.tick_params(axis="y", which="major", labelsize=0, length=0)
    # ax.set_yticks([])

    fig.tight_layout()
    fig.savefig(f"{output_dir}/best_fi.png", dpi=700, bbox_inches="tight")
    fig.savefig(f"{output_dir}/best_fi.pdf", dpi=700, bbox_inches="tight")
    plt.close(fig)

def plot_error_bar(error_dict, cell_name, out_fig_dir):
    os.makedirs(out_fig_dir, exist_ok=True)

    names = list(error_dict.keys())
    vals = np.array([float(error_dict[k]) for k in names], dtype=float)
    vals = np.clip(vals, 0.0, 1.0)

    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    x = np.arange(len(names))
    bars = ax.bar(x, vals, color=dict_colors[cell_name], alpha=0.6)

    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("Error [%]")
    ax.set_title(f"{cell_name} – Optimization results")

    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=30, ha="right")

    # annotate as percentage
    for i, b in enumerate(bars):
        v = float(vals[i])
        ax.text(
            b.get_x() + b.get_width() / 2.0,
            min(1.0, b.get_height() + 0.03),
            f"{100.0*v:.1f}%",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(os.path.join(out_fig_dir, f"{cell_name}_errors.png"), dpi=700, bbox_inches="tight")
    fig.savefig(os.path.join(out_fig_dir, f"{cell_name}_errors.pdf"), dpi=700, bbox_inches="tight")
    plt.close(fig)
