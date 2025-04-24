import matplotlib.pyplot as plt
import numpy as np
import nest
import pandas as pd




k =10000.
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
    "NMDA_g_scaling" : 18.8 * 1 * 35,
    "NMDA_g_peak": 0.052*k,
    #"NMDA_Tau_decay1": 50.

}

# activate AMPA synapse
syn_spec = {"weight":  0.001*1 , "delay": 1, "receptor_type": 2}

# nest settings
nest.ResetKernel()
nest.Install('cerebmodule')
dt =0.025
nest.SetKernelStatus({"resolution": dt})

# input settings
freq = 25 # Hz
spike_interval = 1000 / freq
spike_times =  [float(format(freq + i * 25, ".1f")) for i in range(1)]
print('Input spike times: ', spike_times)


grc = nest.Create('eglif_cond_alpha_multisyn3', 1 , params=params_grc)
input_spikes = nest.Create("spike_generator", params={"spike_times": spike_times})

nest.Connect(input_spikes,grc, syn_spec = dict(syn_spec, synapse_model="static_synapse"))

multimeter = nest.Create("multimeter", params={"interval": 0.1,
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
I_dep = multimeter_grc['I_dep']
I_adap = multimeter_grc['I_adap']
I_syn = multimeter_grc['I_syn']
I_syn_NMDA = multimeter_grc['I_syn_nmda']

NMDA_E_rev = -3.7
g_syn_NMDA = I_syn_NMDA / (NMDA_E_rev - V_m_grc) * Mg_block_grc

print('Number of emitted spikes: ', len(spike_times_grc))
if len(spikes_grc) > 0:
    print('Spikes times ', spike_times_grc)

# print in comparison with neuron trace
file_path = "../NEURON_traces/NMDA_GrC.csv"
df = pd.read_csv(file_path, delim_whitespace=True)
df.columns = ["Time", "Current"]
df["g_NMDA"] = df["Current"].values / (-40.) * 1000   # multiplying 10^3 to convert in nS

plt.figure()
plt.title('NMDA conductance')
plt.plot(df['Time']-104,df['g_NMDA'], color='gray', linestyle='dashed', label='Neuron trace')
plt.plot(times_grc, g_syn_NMDA, 'b', label='NEST trace')
plt.xlabel('Time [ms]')
plt.ylabel(r'$g_{syn_{AMPA}}$ [ns]')
#plt.xlim(0, 450)
plt.legend()
plt.show()

plt.figure()
plt.title('Membrane voltage')
plt.plot(times_grc, V_m_grc, 'b')
plt.xlabel('Time [ms]')
plt.ylabel(r'$V_m$ [mV]')
plt.axhline(params_grc['V_th'], color='r', linestyle='--')
plt.show()

plt.figure()
plt.title('NMDA current')
plt.plot(times_grc, I_syn_NMDA, 'b')
plt.xlabel('Time [ms]')
plt.ylabel(r'$I_{syn_{NMDA}}$ [pA]')
plt.show()

plt.figure()
plt.title('Mg block')
plt.plot(times_grc, Mg_block_grc, 'b')
plt.xlabel('Time [ms]')
plt.ylabel(r'Mg_Block')
plt.show()

print('Max conductance NMDA: ', max(g_syn_NMDA))