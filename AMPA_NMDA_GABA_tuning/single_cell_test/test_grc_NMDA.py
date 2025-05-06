"""Before executing this code fix MgBlock = 0.04 in the .nestml module to simulate voltage clamp protocol to compare synaptic shape with NEURON trace """

import matplotlib.pyplot as plt
import numpy as np
import nest



NMDA_g_scaling = 18.8 * 1 * 1.4     # https://github.com/dbbs-lab/catalogue/tree/main/dbbs_catalogue/mods
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
    "NMDA_block": 1,
    "NMDA_g_init": 7.618718942871291,
    "NMDA_Tau_r": 16.38748791597963,
    "NMDA_Tau_d1":9.363893772753693,
    "NMDA_Tau_d2": 64.42657189216386,
    "NMDA_A_r":  1.7501718153210402,
    "NMDA_A1": 0.007104204842947981,
    "NMDA_A2": 0.000848123147461344,
}

# activate AMPA synapse
k=2
syn_spec = {"weight": 1.*2 , "delay": 1, "receptor_type": 2}

# nest settings
nest.ResetKernel()
nest.Install('cerebmodule')
dt = 0.025
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
                                          "I_syn_gaba_a","Mg_block"], "record_to": "memory", "label": "grc_multimeter"})
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
Mg_block_grc = multimeter_grc['Mg_block']
I_syn_NMDA = multimeter_grc['I_syn_nmda']

NMDA_E_rev = -3.7
g_syn_NMDA = I_syn_NMDA / ((NMDA_E_rev - V_m_grc))#*Mg_block_grc)
print('Number of emitted spikes: ', len(spike_times_grc))
if len(spikes_grc) > 0:
    print('Spikes times ', spike_times_grc)


fig, axs = plt.subplots(4, 1, figsize=(8, 12), sharex=True)

axs[0].set_title('NMDA conductance')
axs[0].plot(times_grc, g_syn_NMDA, 'b', label='NEST trace')
axs[0].set_ylabel(r'$g_{syn_{NMDA}}$ [ns]')
axs[0].set_xlim(0, 400)


axs[1].set_title('Membrane voltage')
axs[1].plot(times_grc, V_m_grc, 'black',label='NEST trace')
axs[1].axhline(params_grc['V_th'], color='r', linestyle='--', label=r'$V_{th}$')
axs[1].set_ylabel(r'$V_m$ [mV]')
axs[1].set_xlim(0, 400)


axs[2].set_title('NMDA current')
axs[2].plot(times_grc, I_syn_NMDA, 'red', label='NEST trace')
axs[2].set_ylabel(r'$I_{syn_{NMDA}}$ [pA]')
axs[2].set_xlim(0, 400)

axs[3].set_title('Mg block')
axs[3].plot(times_grc, Mg_block_grc, 'green',label='NEST trace')
axs[3].set_xlabel('Time [ms]')
axs[3].set_ylabel(r'Mg_Block')
axs[3].set_xlim(0, 400)


for i, ax in enumerate(axs):
    for j, spike in enumerate(spike_times):
        label = 'input spikes' if i == 0 and j == 0 else None
        ax.axvline(spike, color='gray', linestyle='--', linewidth=0.5, label=label)


axs[0].legend()
axs[1].legend()
axs[2].legend()
axs[3].legend()


plt.tight_layout()
plt.savefig('NMDA_test_10inputspikes.png', dpi=300)
plt.show()

