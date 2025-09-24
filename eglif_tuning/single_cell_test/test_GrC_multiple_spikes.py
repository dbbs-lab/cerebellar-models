import matplotlib.pyplot as plt
import nest
import numpy as np

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
    "GABA_E_rev": -70.0000001,  # to ensure no division for 0
}


# activate AMPA NMDA and GABA synapses
k = 1.0
syn_spec1 = {"weight": k * 3.0, "delay": 1, "receptor_type": 1}
syn_spec2 = {"weight": k * 3.0, "delay": 1, "receptor_type": 2}
syn_spec3 = {"weight": k * 1.0, "delay": 1, "receptor_type": 3}

# nest settings
nest.ResetKernel()
nest.Install("cerebmodule")
dt = 0.025
nest.SetKernelStatus({"resolution": dt})

# input settings
freq = 25  # Hz
spike_interval = 1000 / freq
spike_times = [float(format(freq + i * 25, ".1f")) for i in range(10)]


grc = nest.Create("eglif_multisyn", 1, params=params_grc)
input_spikes = nest.Create("spike_generator", params={"spike_times": spike_times})

nest.Connect(input_spikes, grc, syn_spec=dict(syn_spec1, synapse_model="static_synapse"))
nest.Connect(input_spikes, grc, syn_spec=dict(syn_spec2, synapse_model="static_synapse"))
nest.Connect(input_spikes, grc, syn_spec=dict(syn_spec3, synapse_model="static_synapse"))

multimeter = nest.Create(
    "multimeter",
    params={
        "interval": 0.025,
        "record_from": [
            "V_m",
            "V_m_2",
            "I_dep",
            "I_adap",
            "I_syn",
            "I_syn_ampa",
            "I_syn_nmda",
            "I_syn_gaba",
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


duration = np.max(spike_times) + 200
print("The simulation duration is ", duration)
nest.Simulate(duration)

# Collect results
multimeter_grc = multimeter.get()["events"]
spikes_grc = spike_recorder.get("events")

spike_times_grc = spikes_grc["times"]
times_grc = multimeter_grc["times"]

V_m_grc = multimeter_grc["V_m"]
Mg_block_grc = multimeter_grc["Mg_block"]
I_syn_NMDA = multimeter_grc["I_syn_nmda"]
I_syn_AMPA = multimeter_grc["I_syn_ampa"]
I_syn_GABA = multimeter_grc["I_syn_gaba"]
I_syn = multimeter_grc["I_syn"]


g_syn_NMDA = I_syn_NMDA / (params_grc["NMDA_E_rev"] - V_m_grc)
g_syn_AMPA = I_syn_AMPA / (params_grc["AMPA_E_rev"] - V_m_grc)
g_syn_GABA = I_syn_GABA / (params_grc["GABA_E_rev"] - V_m_grc)

print("Number of emitted spikes: ", len(spike_times_grc))
if len(spikes_grc) > 0:
    print("Spikes times ", spike_times_grc)
    print("Number of emitted spikes: ", len(spike_times_grc))


fig = plt.figure(figsize=(10, 15))
gs = fig.add_gridspec(5, 2)

ax_i_syn = fig.add_subplot(gs[4, :])
ax_i_syn.set_title("Total Synaptic Current ($I_{syn}$)")
ax_i_syn.plot(times_grc, I_syn, "purple", label="Total $I_{syn}$")
ax_i_syn.set_ylabel(r"$I_{syn}$ [pA]")
ax_i_syn.set_xlabel("Time [ms]")
ax_i_syn.legend(loc="upper right", fontsize=8)
ax_i_syn.set_xlim(0, duration)
axs = np.empty((4, 2), dtype=object)

for i in range(4):
    for j in range(2):
        axs[i, j] = fig.add_subplot(gs[i, j], sharex=ax_i_syn)
        axs[i, j].set_xlim(0, duration)
        if i < 3:
            axs[i, j].tick_params(labelbottom=False)

# Membrane Voltage and Mg Block
axs[0][0].set_title("Membrane Voltage")
axs[0][0].plot(times_grc, V_m_grc, "black", label="NEST trace")
axs[0][0].axhline(params_grc["V_th"], color="r", linestyle="--", label=r"$V_{th}$")
axs[0][0].set_ylabel(r"$V_m$ [mV]")
axs[0][0].legend(loc="upper right", fontsize=8)
axs[0][0].tick_params(labelbottom=False)

axs[0][1].set_title("Mg Block")
axs[0][1].plot(times_grc, Mg_block_grc, "green")
axs[0][1].set_ylabel("Mg_Block")
axs[0][1].tick_params(labelbottom=False)

# NMDA Conductance and Current
axs[1][0].set_title("NMDA Conductance")
axs[1][0].plot(times_grc, g_syn_NMDA, "b")
axs[1][0].set_ylabel(r"$g_{syn_{NMDA}}$ [nS]")
axs[1][0].tick_params(labelbottom=False)

axs[1][1].set_title("NMDA Current")
axs[1][1].plot(times_grc, I_syn_NMDA, "red")
axs[1][1].set_ylabel(r"$I_{syn_{NMDA}}$ [pA]")
axs[1][1].tick_params(labelbottom=False)

# AMPA Conductance and Current
axs[2][0].set_title("AMPA Conductance")
axs[2][0].plot(times_grc, g_syn_AMPA, "b")
axs[2][0].set_ylabel(r"$g_{syn_{AMPA}}$ [nS]")
axs[2][0].tick_params(labelbottom=False)

axs[2][1].set_title("AMPA Current")
axs[2][1].plot(times_grc, I_syn_AMPA, "red")
axs[2][1].set_ylabel(r"$I_{syn_{AMPA}}$ [pA]")
axs[2][1].tick_params(labelbottom=False)

# GABA Conductance and Current
axs[3][0].set_title("GABA Conductance")
axs[3][0].plot(times_grc, g_syn_GABA, "b")
axs[3][0].set_ylabel(r"$g_{syn_{GABA}}$ [nS]")

axs[3][1].set_title("GABA Current")
axs[3][1].plot(times_grc, I_syn_GABA, "red")
axs[3][1].set_ylabel(r"$I_{syn_{GABA}}$ [pA]")


for k, spike in enumerate(spike_times):
    label = "input spikes" if k == 0 else None
    axs[0][0].axvline(spike, color="gray", linestyle="--", linewidth=0.5, label=label)


plt.tight_layout()
plt.savefig(f"./figs/AMPA_NMDA_w3_GABA_w1_10spikes.png", dpi=300)
plt.show()
