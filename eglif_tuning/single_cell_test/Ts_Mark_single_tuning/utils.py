import json
import os

import matplotlib.pyplot as plt
import nest
import numpy as np
from matplotlib.ticker import MaxNLocator

dict_colors = {
    "BC": [1, 0.647, 0, 1.0],
    "SC": [1, 0.84, 0, 1.0],
    "PC": [0.275, 0.800, 0.275, 1.0],
    "GoC": [0, 0.45, 0.7, 1.0],
    "GrC": [0.7, 0.15, 0.15, 1.0],
}

dict_names = {
    "BC": "Basket Cell",
    "SC": "Stellate Cell",
    "PC": "Purkinje Cell",
    "GoC": "Golgi Cell",
    "GrC": "Granule Cell",
}


def simulate(
    cell_name,
    cell_params,
    receptor,
    syn_spec,
    receptor_dict,
    syn_name,
    freq,
    syn_model="static_synapse",
    sim_duration=5000.0,
    spike_times=None,
    PLOT=False,
    SAVE=False,
):

    nest.ResetKernel()
    nest.Install("cerebmodule")
    dt = 0.1
    nest.SetKernelStatus({"resolution": dt})
    if spike_times is None:
        spike_times = [200.0]
    input_spikes = nest.Create("spike_generator", params={"spike_times": spike_times})
    parrot = nest.Create("parrot_neuron", 1)
    neuron = nest.Create("eglif_multirec_opt", 1, params=cell_params)
    neuron.V_m = cell_params["E_L"]
    record_vars = ["V_m", receptor_dict[receptor][1]]
    if receptor == "2":
        record_vars.append("Mg_block")
    multimeter = nest.Create(
        "multimeter", params={"interval": dt, "record_from": record_vars, "record_to": "memory"}
    )

    nest.Connect(multimeter, neuron)
    nest.Connect(input_spikes, parrot)
    nest.Connect(parrot, neuron, syn_spec=dict(syn_spec, synapse_model=syn_model))

    nest.Simulate(sim_duration)

    events = multimeter.get()["events"]
    V_m = events["V_m"]
    I_syn = events[receptor_dict[receptor][1]]
    time = events["times"]

    if receptor == "2":
        Mg_block = events["Mg_block"]
        g = I_syn / ((cell_params[receptor_dict[receptor][2]] - V_m) * Mg_block)
    else:
        g = I_syn / (cell_params[receptor_dict[receptor][2]] - V_m)

    if PLOT:
        fig, ax = plt.subplots(figsize=(5.5, 3.5))
        color = dict_colors.get(cell_name, [0, 0, 0, 1.0])
        ax.plot(time - spike_times[0], g, color=color, lw=1.5)
        model_label = "Static" if syn_model == "static_synapse" else "with STP"
        freq_label = ""
        if len(spike_times) > 1:
            freq_label = f"- $f_{{input}} = {freq:.1f}\\,\\mathrm{{Hz}}$"

        title = f"{dict_names[cell_name]} {receptor_dict[receptor][0]} {model_label} {freq_label}"
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("Time [ms]", fontsize=10)
        ax.set_ylabel(r"$g_{syn}$ [nS]", fontsize=10)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))

        ax.tick_params(axis="both", labelsize=9)
        ax.set_xlim(0, spike_times[10] - spike_times[0])

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.tight_layout()
        if SAVE:
            os.makedirs("figs", exist_ok=True)
            plt.savefig(f"./figs/{syn_name}_{freq}.png")
        plt.show()

    return g, time


def mean_conductance(g_trace, t_trace, t_start):
    mask = t_trace >= t_start
    return np.trapz(g_trace[mask], t_trace[mask]) / (t_trace[mask][-1] - t_trace[mask][0])


def simulate_multirec(cell_params, spike_times, syn_model, list_syn_spec, rec_vars, sim_duration):
    nest.ResetKernel()
    nest.Install("cerebmodule")
    dt = 0.1
    nest.SetKernelStatus({"resolution": dt})
    input_spikes = nest.Create("spike_generator", params={"spike_times": spike_times})
    parrot = nest.Create("parrot_neuron", 1)
    neuron = nest.Create("eglif_multirec_opt", 1, params=cell_params)
    neuron.V_m = cell_params["E_L"]
    multimeter = nest.Create(
        "multimeter",
        params={
            "interval": 0.1,
            "record_from": rec_vars,
            "record_to": "memory",
            "label": "multimeter",
        },
    )
    spikes_rec = nest.Create("spike_recorder")
    nest.Connect(multimeter, neuron)
    nest.Connect(input_spikes, parrot)
    nest.Connect(neuron, spikes_rec)
    for syn_spec in list_syn_spec:
        nest.Connect(parrot, neuron, syn_spec=dict(syn_spec, synapse_model=syn_model))

    nest.Simulate(sim_duration)
    multimeter = multimeter.get()["events"]
    spikes = spikes_rec.get("events")

    return multimeter, spikes


def plot_ff_curve(freqs, fr_stat, fr_tm):
    plt.figure()
    plt.plot(freqs, fr_stat, label="fr_stat")
    plt.plot(freqs, fr_tm, label="fr_tm")
    plt.legend()
    plt.show()


def run_static_vs_tm_test(
    cell_name,
    cell_params,
    receptor_type,
    receptor_dict,
    syn_static,
    syn_tm,
    syn_base_name,
    freq_ref=5.0,
    test_freqs=(5.0, 50.0, 100.0, 150.0, 250.0, 500.0),
    n_spikes_ref=50,
    n_spikes_test=30,
    sim_duration=5000.0,
):

    isi = 1000.0 / freq_ref
    spike_times = np.round([20.0 + i * isi for i in range(n_spikes_ref)], 1)

    g_static, t_static = simulate(
        cell_name,
        cell_params=cell_params,
        receptor=receptor_type,
        receptor_dict=receptor_dict,
        syn_spec=syn_static,
        syn_name=f"{syn_base_name}_static",
        freq=freq_ref,
        spike_times=list(spike_times),
        sim_duration=sim_duration,
        PLOT=True,
    )

    g_mean_static = mean_conductance(g_static, t_static, t_start=spike_times[10])

    g_tm, t_tm = simulate(
        cell_name,
        cell_params=cell_params,
        receptor=receptor_type,
        receptor_dict=receptor_dict,
        syn_spec=syn_tm,
        syn_model="tsodyks_synapse",
        syn_name=f"{syn_base_name}_tm",
        freq=freq_ref,
        spike_times=list(spike_times),
        sim_duration=sim_duration,
        PLOT=True,
    )

    g_mean_tm = mean_conductance(g_tm, t_tm, t_start=spike_times[10])

    w_star = g_mean_static / g_mean_tm
    syn_tm["weight"] *= w_star

    print(f"[TUNING] Reference frequency: {freq_ref} Hz")
    print(f"[TUNING] Mean g static = {g_mean_static:.4e}")
    print(f"[TUNING] Mean g TM (w=1) = {g_mean_tm:.4e}")
    print(f"[TUNING] Scaled TM weight = {syn_tm['weight']:.4f}")

    os.makedirs("results_TM", exist_ok=True)
    out_path = f"results_TM/TM_{syn_base_name}.json"
    with open(out_path, "w") as f:
        json.dump(syn_tm, f, indent=4)

    print(f"TM parameters saved to {out_path}")

    for freq in test_freqs:
        isi = 1000.0 / freq
        spike_times = np.round([20.0 + i * isi for i in range(n_spikes_test)], 1)

        print(f"\n[TEST] Frequency = {freq} Hz")

        simulate(
            cell_name,
            cell_params=cell_params,
            receptor=receptor_type,
            receptor_dict=receptor_dict,
            syn_spec=syn_static,
            syn_name=f"{syn_base_name}_static",
            freq=freq,
            spike_times=list(spike_times),
            PLOT=True,
            SAVE=True,
        )

        simulate(
            cell_name,
            cell_params=cell_params,
            receptor=receptor_type,
            receptor_dict=receptor_dict,
            syn_spec=syn_tm,
            syn_model="tsodyks_synapse",
            syn_name=f"{syn_base_name}_tm",
            freq=freq,
            spike_times=list(spike_times),
            PLOT=True,
            SAVE=True,
        )
