import matplotlib.pyplot as plt
import numpy as np
import nest
import pandas as pd

AMPA_g_scaling = 1.2 * 1 * 1.4
k =2.
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
    "AMPA_g_init": 49.99770425187628,
    "AMPA_Tau_r": 0.03860269857016296,
    "AMPA_Tau_d1": 0.32894309944849703,
    "AMPA_Tau_d2":3.71405056519806,
    "AMPA_A_r": 0.8708806368411421,
    "AMPA_A1": 0.07264536636611472,
    "AMPA_A2": 0.14681591004769137,
}

# activate AMPA synapse
syn_spec = {"weight":  1., "delay": 1, "receptor_type": 1}

# nest settings
nest.ResetKernel()
nest.Install('cerebmodule')
dt =0.025
nest.SetKernelStatus({"resolution": dt})

# input settings
freq = 25 # Hz
spike_interval = 1000 / freq
spike_times =  [float(format(freq + i * 25, ".1f")) for i in range(10)]
print('Input spike times: ', spike_times)


grc = nest.Create('eglif_cond_alpha_multisyn3', 1 , params=params_grc)
input_spikes = nest.Create("spike_generator", params={"spike_times": spike_times})

nest.Connect(input_spikes,grc, syn_spec = dict(syn_spec, synapse_model="static_synapse"))

multimeter = nest.Create("multimeter", params={"interval": 0.025,
                          "record_from": ["V_m","V_m_2","I_dep","I_adap","I_syn","I_syn_ampa","I_syn_nmda",
                                          "I_syn_gaba_a", "Mg_block"], "record_to": "memory", "label": "grc_multimeter"})
spike_recorder  = nest.Create("spike_recorder",params={"record_to": "memory", "label": "grc_spike_recorder"})
nest.Connect(multimeter, grc)
nest.Connect(grc, spike_recorder)

duration = np.max(spike_times)+200
print('The simulation duration is ', duration)

nest.Simulate(1000)

# Collect results
multimeter_grc= multimeter.get()['events']
spikes_grc = spike_recorder.get('events')

spike_times_grc = spikes_grc["times"]
times_grc = multimeter_grc['times']

V_m_grc = multimeter_grc['V_m']
I_dep = multimeter_grc['I_dep']
I_adap = multimeter_grc['I_adap']
I_syn = multimeter_grc['I_syn']
I_syn_AMPA = multimeter_grc['I_syn_ampa']

I_syn_NMDA = multimeter_grc['I_syn_nmda']
I_syn_GABA_A = multimeter_grc['I_syn_gaba_a']

AMPA_E_rev = 0.0
g_syn_AMPA = I_syn_AMPA / (AMPA_E_rev - V_m_grc)

print('Number of emitted spikes: ', len(spikes_grc))
if len(spikes_grc) > 0:
    print('Spikes times ', spike_times_grc)

# print in comparison with neuron trace
file_path = "../NEURON_traces/AMPA_GrC.csv"
df = pd.read_csv(file_path, delim_whitespace=True)
df.columns = ["Time", "Current"]
df["g_AMPA"] = df["Current"].values / (-40.) * 1000   # multiplying 10^3 to convert in nS


fig, axs = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

axs[0].set_title('AMPA conductance')
# axs[0].plot(df['Time'], df['g_AMPA'], color='gray', linestyle='dashed', label='Neuron trace')
axs[0].plot(times_grc, g_syn_AMPA, 'b', label='NEST trace')
axs[0].set_ylabel(r'$g_{syn_{AMPA}}$ [ns]')
axs[0].legend()
axs[0].set_xlim(0, 400)


axs[1].set_title('Membrane voltage')
axs[1].plot(times_grc, V_m_grc, 'black', label='NEST trace')
axs[1].axhline(params_grc['V_th'], color='r', linestyle='--', label=r'$V_{th}$')
axs[1].set_ylabel(r'$V_m$ [mV]')
axs[1].legend()
axs[1].set_xlim(0, 400)


axs[2].set_title('AMPA current')
axs[2].plot(times_grc, I_syn_AMPA, 'r', label='NEST trace')
axs[2].set_xlabel('Time [ms]')
axs[2].set_ylabel(r'$I_{syn_{AMPA}}$ [pA]')
axs[2].set_xlim(0, 400)


for i, ax in enumerate(axs):
    for j, spike in enumerate(spike_times):
        label = 'input spikes' if i == 0 and j == 0 else None
        ax.axvline(spike, color='gray', linestyle='--', linewidth=0.5, label=label)

axs[0].legend()
axs[1].legend()
axs[2].legend()

plt.tight_layout()
plt.savefig('AMPA_10inputspikes.png', dpi=300)
plt.show()