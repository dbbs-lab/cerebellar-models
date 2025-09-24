import matplotlib.pyplot as plt
import nest

params_grc = {
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
    "AMPA_E_rev": 0,
    "NMDA_E_rev": -3.7,
    "GABA_E_rev": -70.0,
}


# activate AMPA and NMDA synapse
w = 2
syn_spec1 = {"weight": w, "delay": 1, "receptor_type": 1}
syn_spec2 = {"weight": w, "delay": 1, "receptor_type": 2}
syn_spec3 = {"weight": 1, "delay": 1, "receptor_type": 3}

# nest settings
nest.ResetKernel()
nest.Install("cerebmodule")
dt = 0.025
nest.SetKernelStatus({"resolution": dt})
grc = nest.Create("eglif_multisyn", 1, params=params_grc)

input_spikes = nest.Create("spike_generator", params={"spike_times": [250]})

nest.Connect(input_spikes, grc, syn_spec=dict(syn_spec1, synapse_model="static_synapse"))
nest.Connect(input_spikes, grc, syn_spec=dict(syn_spec2, synapse_model="static_synapse"))
nest.Connect(input_spikes, grc, syn_spec=dict(syn_spec3, synapse_model="static_synapse"))

multimeter = nest.Create(
    "multimeter",
    params={
        "interval": 0.025,
        "record_from": [
            "V_m",
            "I_syn_ampa",
            "I_syn_nmda",
            "I_ampa1",
            "I_ampa2",
            "I_nmda1",
            "I_nmda2",
            "I_syn_gaba",
            "I_gaba1",
            "I_gaba2",
            "Mg_block",
        ],
        "record_to": "memory",
        "label": "grc_multimeter",
    },
)
spike_recorder = nest.Create(
    "spike_recorder", params={"record_to": "memory", "label": "grc_spike_recorder"}
)
nest.Connect(multimeter, grc)
nest.Connect(grc, spike_recorder)

nest.Simulate(1000)

# Collect results
multimeter_grc = multimeter.get()["events"]
spikes_grc = spike_recorder.get("events")

spike_times_grc = spikes_grc["times"]
times_grc = multimeter_grc["times"]

V_m_grc = multimeter_grc["V_m"]
I_syn_NMDA = multimeter_grc["I_syn_nmda"]
I_syn_AMPA = multimeter_grc["I_syn_ampa"]
I_ampa1 = multimeter_grc["I_ampa1"]
I_ampa2 = multimeter_grc["I_ampa2"]
I_nmda1 = multimeter_grc["I_nmda1"]
I_nmda2 = multimeter_grc["I_nmda2"]
I_syn_GABA = multimeter_grc["I_syn_gaba"]
I_gaba1 = multimeter_grc["I_gaba1"]
I_gaba2 = multimeter_grc["I_gaba2"]
Mg_block = multimeter_grc["Mg_block"]


NMDA_E_rev = -3.7
g_syn_NMDA = I_syn_NMDA / ((NMDA_E_rev - V_m_grc) * Mg_block)
g_syn_nmda1 = I_nmda1 / ((NMDA_E_rev - V_m_grc) * Mg_block)
g_syn_nmda2 = I_nmda2 / ((NMDA_E_rev - V_m_grc) * Mg_block)

AMPA_E_rev = 0.0
g_syn_AMPA = I_syn_AMPA / (AMPA_E_rev - V_m_grc)
g_syn_ampa1 = I_ampa1 / (AMPA_E_rev - V_m_grc)
g_syn_ampa2 = I_ampa2 / (AMPA_E_rev - V_m_grc)

GABA_E_rev = -70.0
g_syn_GABA = I_syn_GABA / (GABA_E_rev - V_m_grc)
g_syn_gaba1 = I_gaba1 / (GABA_E_rev - V_m_grc)
g_syn_gaba2 = I_gaba2 / (GABA_E_rev - V_m_grc)


fig, ax = plt.subplots(4, 1, figsize=(10, 10), sharex=True, sharey=False)

ax[0].set_title(f"mf_GrC AMPA Conductance with kernels - w={w}")
ax[0].plot(times_grc, g_syn_AMPA, label="G_syn_AMPA", color="r")
ax[0].plot(
    times_grc, g_syn_ampa1, label="G_syn_ampa1 (fast)", color="gray", alpha=0.5, linestyle="--"
)
ax[0].plot(
    times_grc, g_syn_ampa2, label="G_syn_ampa2 (slow)", color="blue", alpha=0.5, linestyle="--"
)
ax[0].set_xlim(0, 1000)
ax[0].set_ylabel("Conductance [nS]")
ax[0].legend(loc="upper left")

ax[1].set_title(f"mf_GrC NMDA Conductance with kernels (\Mg_block)- w={w}")
ax[1].plot(times_grc, g_syn_NMDA, label="G_syn_NMDA", color="r")
ax[1].plot(
    times_grc, g_syn_nmda1, label="G_syn_nmda1 (fast)", color="gray", alpha=0.5, linestyle="--"
)
ax[1].plot(
    times_grc, g_syn_nmda2, label="G_syn_nmda2 (slow)", color="blue", alpha=0.5, linestyle="--"
)
ax[1].set_xlim(0, 1000)
ax[1].set_ylabel("Conductance [nS]")
ax[1].legend(loc="upper left")

ax[2].set_title(f"mf_GrC NMDA Conductance with kernels - w={w}")
ax[2].plot(times_grc, g_syn_NMDA * Mg_block, label="G_syn_NMDA", color="r")
ax[2].plot(
    times_grc,
    g_syn_nmda1 * Mg_block,
    label="G_syn_nmda1 (fast)",
    color="gray",
    alpha=0.5,
    linestyle="--",
)
ax[2].plot(
    times_grc,
    g_syn_nmda2 * Mg_block,
    label="G_syn_nmda2 (slow)",
    color="blue",
    alpha=0.5,
    linestyle="--",
)
ax[2].set_xlim(0, 1000)
ax[2].set_ylabel("Conductance [nS]")
ax[2].legend(loc="upper left")

ax[3].set_title(f"GoC_Grc GABA conductance with kernels - w=1")
ax[3].plot(times_grc, g_syn_GABA, label="G_syn_GABA", color="r")
ax[3].plot(
    times_grc, g_syn_gaba1, label="G_syn_gaba1 (fast)", color="gray", alpha=0.5, linestyle="--"
)
ax[3].plot(
    times_grc, g_syn_gaba2, label="G_syn_gaba2 (slow)", color="blue", alpha=0.5, linestyle="--"
)
ax[3].set_xlim(0, 1000)
ax[3].set_xlabel("Time [ms]")
ax[3].set_ylabel("Conductance [nS]")
ax[3].legend(loc="upper left")

plt.tight_layout()
plt.savefig(f"./figs/GrC_conducatances_w{w}.png", dpi=300)
plt.show()


fig, ax = plt.subplots(3, 1, figsize=(10, 10), sharex=True, sharey=False)

ax[0].set_title(f"mf_Grc AMPA Current (Excitatory)- w={w}")
ax[0].plot(times_grc, I_syn_AMPA, label="I_syn_AMPA", color="black")
ax[0].set_xlim(0, 1000)
ax[0].set_ylabel("Current [pA]")
ax[0].legend(loc="upper left")

ax[1].set_title(f"mf_GrC NMDA Current (Excitatory) - w={w}")
ax[1].plot(times_grc, I_syn_NMDA, label="I_syn_NMDA", color="black")
ax[1].set_xlim(0, 1000)
ax[1].legend(loc="upper left")
ax[1].set_ylabel("Current [pA]")

ax[2].set_title(f"GoC_GrC GABA current (Inhibitory)- w=1")
ax[2].plot(times_grc, I_syn_GABA, label="I_syn_GABA", color="black")
ax[2].set_xlim(0, 1000)
ax[2].set_xlabel("Time [ms]")
ax[2].set_ylabel("Current [pA]")
ax[2].legend(loc="upper left")

plt.tight_layout()
plt.savefig(f"./figs/GrC_currents_w{w}.png", dpi=300)
plt.show()
