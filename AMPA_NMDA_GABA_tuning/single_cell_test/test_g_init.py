"""
Test g_init on AMPA synapse
"""

import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import nest
import numpy as np
from scipy.stats import linregress


def plot_conductance_g(times_sim, g_sim, gi, g_init):
    norm = colors.Normalize(vmin=min(g_init), vmax=max(g_init))
    cmap = cm.Reds
    color = cmap(norm(gi))
    plt.plot(times_sim, g_sim, color=color, alpha=0.7, label=f"{gi}")


def plot_conductance_A(
    times_sim,
    g_sim,
    ar,
    Ar,
):
    norm = colors.Normalize(vmin=min(Ar), vmax=max(Ar))
    cmap = cm.Reds
    color = cmap(norm(ar))
    plt.plot(times_sim, g_sim, color=color, alpha=0.7, label=f"Ar=A1=A2={ar}")


def single_sim_AMPA(cell_params, syn_spec):
    nest.ResetKernel()
    nest.Install("cerebmodule")
    nest.resolution = 0.025
    grc = nest.Create("eglif_multisyn", 1, params=cell_params)
    input_spikes = nest.Create("spike_generator", params={"spike_times": [250]})
    nest.Connect(input_spikes, grc, syn_spec=dict(syn_spec, synapse_model="static_synapse"))

    multimeter = nest.Create(
        "multimeter",
        params={
            "interval": 0.025,
            "record_from": ["V_m", "I_syn_ampa"],
            "record_to": "memory",
            "label": "grc_multimeter",
        },
    )

    spike_recorder = nest.Create(
        "spike_recorder", params={"record_to": "memory", "label": "grc_recorder"}
    )

    nest.Connect(multimeter, grc)
    nest.Connect(grc, spike_recorder)
    duration = 1000.0
    nest.Simulate(duration)

    multimeter_data = multimeter.get()["events"]
    times = multimeter_data["times"]
    I_syn_ampa = multimeter_data["I_syn_ampa"]
    V_m = multimeter_data["V_m"]
    g_syn_ampa = I_syn_ampa / (cell_params["AMPA_E_rev"] - V_m)

    return times, g_syn_ampa


def single_sim_NMDA(cell_params, syn_spec):
    nest.ResetKernel()
    nest.Install("cerebmodule")
    nest.resolution = 0.025
    grc = nest.Create("eglif_multisyn", 1, params=cell_params)
    input_spikes = nest.Create("spike_generator", params={"spike_times": [250]})
    nest.Connect(input_spikes, grc, syn_spec=dict(syn_spec, synapse_model="static_synapse"))

    multimeter = nest.Create(
        "multimeter",
        params={
            "interval": 0.025,
            "record_from": ["V_m", "I_syn_nmda", "Mg_block"],
            "record_to": "memory",
            "label": "grc_multimeter",
        },
    )

    spike_recorder = nest.Create(
        "spike_recorder", params={"record_to": "memory", "label": "grc_recorder"}
    )

    nest.Connect(multimeter, grc)
    nest.Connect(grc, spike_recorder)
    duration = 1000.0
    nest.Simulate(duration)

    multimeter_data = multimeter.get()["events"]
    times = multimeter_data["times"]
    I_syn_ampa = multimeter_data["I_syn_nmda"]
    V_m = multimeter_data["V_m"]
    Mg_Block = multimeter_data["Mg_block"]
    g_syn_nmda = I_syn_ampa / ((cell_params["NMDA_E_rev"] - V_m) * Mg_Block)

    return times, g_syn_nmda


cell_params = {
    "t_ref": 1.5,
    "V_min": -150,
    "C_m": 7,
    "V_th": -41,
    "V_reset": -70,
    "E_L": -62,
    "I_e": -0.888,
    "V_m": -62,
    "lambda_0": 1,
    "tau_V": 0.3,
    "tau_m": 24.15,
    "k_adap": 0.022,
    "k_1": 0.311,
    "k_2": 0.041407868,
    "A1": 0.01,
    "A2": -0.94,
    "AMPA_E_rev": 0.0,
    "NMDA_E_rev": -3.7,
}


syn_spec_AMPA = {"weight": 1.0, "delay": 1, "receptor_type": 1}

syn_spec_NMDA = {"weight": 1.0, "delay": 1, "receptor_type": 2}

## AMPA
g_init = np.arange(0, 5.5, 0.5)
plt.title(f"mf_GrC AMPA Conductance for different g_init")
peaks = []
for gi in g_init:
    cell_params["AMPA_g_init"] = gi
    times, g_syn_ampa = single_sim_AMPA(cell_params, syn_spec_AMPA)
    plot_conductance_g(times, g_syn_ampa, gi, g_init)
    peaks.append(g_syn_ampa.max())


plt.xlabel("Time [ms]")
plt.ylabel(r"$g_{syn}$ [ns]")
plt.legend(title="g_init value", bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.0)
plt.xlim(240, 275)
plt.ylim(0, 1.6)
plt.tight_layout()
plt.savefig("./figs/g_init_ampa.png", dpi=300)
plt.show()

m, q, _, _, _ = linregress(g_init, peaks)
fit = m * g_init + q

plt.figure()
plt.title("mf_Grc AMPA conducatance Peak vs g_init relationship")
plt.scatter(g_init, peaks, c="black")
plt.plot(
    g_init, fit, "--", label=f"Linear regression: y = {m:.2f}x + {q:.2f}", alpha=0.5, color="orange"
)
plt.xlabel("g_init value")
plt.ylabel("Peak value")
plt.legend()
plt.savefig("./figs/g_init-peaks_AMPA.png", dpi=300)
plt.show()

## NMDA
g_init = np.arange(0, 0.2, 0.03)
plt.title(f"mf_GrC NMDA Conductance for different g_init")
peaks = []
for gi in g_init:
    cell_params["NMDA_g_init"] = gi
    times, g_syn_nmda = single_sim_NMDA(cell_params, syn_spec_NMDA)
    plot_conductance_g(times, g_syn_nmda, gi, g_init)
    peaks.append(g_syn_nmda.max())


plt.xlabel("Time [ms]")
plt.ylabel(r"$g_{syn}$ [ns]")
plt.legend(title="g_init value", bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.0)
plt.tight_layout()
plt.savefig("./figs/g_init_nmda.png", dpi=300)
plt.show()

m, q, _, _, _ = linregress(g_init, peaks)
fit = m * g_init + q

plt.figure()
plt.title("mf_Grc NMDA conducatance Peak vs g_init relationship")
plt.scatter(g_init, peaks, c="black")
plt.plot(
    g_init, fit, "--", label=f"Linear regression: y = {m:.2f}x + {q:.2f}", alpha=0.5, color="orange"
)
plt.xlabel("g_init value")
plt.ylabel("Peak value")
plt.legend()
plt.savefig("./figs/g_init-peaks_NMDA.png", dpi=300)
plt.show()


# TEST HYPOTHESIS:  g_init is a multiplicative common factor
"""cell_params1 = copy.deepcopy(cell_params)
cell_params1['AMPA_A_r'] = cell_params1['AMPA_A_r'] *cell_params['AMPA_g_init']
cell_params1['AMPA_A1'] = cell_params1['AMPA_A1'] *cell_params['AMPA_g_init']
cell_params1['AMPA_A2'] = cell_params1['AMPA_A2'] *cell_params['AMPA_g_init']
cell_params1['AMPA_g_init'] = 1.

times, g_syn_ampa = single_sim(cell_params, syn_spec)
times1, g_syn_ampa1 = single_sim(cell_params1, syn_spec)

plt.figure()
plt.title('Hypotesis: g_init common scaling factor (REJECTED)')
plt.plot(times, g_syn_ampa, label='g_syn_ampa', color = 'r', alpha=0.6)
plt.plot(times1, g_syn_ampa1, label='g_syn_ampa1', color ='black', ls='--',alpha = 0.6)
plt.legend()
plt.xlabel('Time [ms]')
plt.ylabel(r'$g_{syn}$ [ns]')
plt.xlim(240,280)
plt.tight_layout()
plt.show()"""
