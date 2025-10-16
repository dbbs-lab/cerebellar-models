import numpy as np
import nest
import matplotlib.pyplot as plt


def run_simulation_static(params_grc, tot_spikes):
    nest.ResetKernel()
    nest.Install("cerebmodule")

    glom = nest.Create("parrot_neuron", 1)
    grc = nest.Create("eglif_multirec_opt", 1, params=params_grc)

    sp_times = np.arange(1, 2001, tot_spikes, dtype=float)
    glom_generator_1 = nest.Create("spike_generator", params={"spike_times": sp_times})
    # glom_generator_2 = nest.Create("poisson_generator", params={"rate": 4, "start": 0, "stop": 1000})

    nest.Connect(glom, grc, syn_spec={"synapse_model": "static_synapse", "weight": 1, "delay": 1, "receptor_type" : 1})
    nest.Connect(glom, grc, syn_spec={"synapse_model": "static_synapse", "weight": 1, "delay": 1, "receptor_type": 2})
    #nest.Connect(glom, grc, syn_spec={"synapse_model": "static_synapse", "weight": 1, "delay": 1, "receptor_type": 3})
    mult_grc = nest.Create( "multimeter", params={"interval": 0.1, "record_from": ["V_m", "I_syn_ampa", "I_syn_nmda", "I_syn_gaba", "I_syn"], "record_to": "memory"})
    rec_grc = nest.Create("spike_recorder")

    nest.Connect(mult_grc, grc)
    nest.Connect(grc, rec_grc)
    nest.Connect(glom_generator_1, glom)

    nest.Simulate(2200)
    mult_events = mult_grc.get()["events"]
    rec_spikes = rec_grc.get()["events"]
    times = mult_events["times"]
    I_syn_ampa = mult_events["I_syn_ampa"]
    I_syn_nmda = mult_events["I_syn_nmda"]
    I_syn_gaba = mult_events["I_syn_gaba"]
    I_syn = mult_events["I_syn"]
    V_m = mult_events["V_m"]
    g_syn_ampa = I_syn_ampa / (params_grc["AMPA_E_rev"] - V_m)
    g_syn_nmda = I_syn_nmda / (params_grc["NMDA_E_rev"] - V_m)
    g_syn_gaba = I_syn_gaba / (params_grc["GABA_E_rev"] - V_m)

    return g_syn_ampa, g_syn_nmda, g_syn_gaba, times, I_syn_ampa, I_syn_nmda, I_syn_gaba, I_syn, V_m, rec_spikes

def run_simulation_stp(params_grc, tot_spikes):
    nest.ResetKernel()
    nest.Install("cerebmodule")

    glom = nest.Create("parrot_neuron", 1)
    grc = nest.Create("eglif_multirec_opt", 1, params=params_grc)

    sp_times = np.arange(1, 2001, tot_spikes, dtype=float)
    glom_generator_1 = nest.Create("spike_generator", params={"spike_times": sp_times})
    # glom_generator_2 = nest.Create("poisson_generator", params={"rate": 4, "start": 0, "stop": 1000})

    nest.Connect(glom, grc, syn_spec={"synapse_model" : "tsodyks2_synapse" ,"weight": 2, "U": 0.43, "x": 1, "tau_rec": 8, "tau_fac": 5, "delay": 1, "receptor_type" : 1})
    nest.Connect(glom, grc, syn_spec={"synapse_model": "tsodyks2_synapse", "weight": 2, "U": 0.43, "x": 1, "tau_rec": 8,"tau_fac": 5, "delay": 1, "receptor_type": 2})

    mult_grc = nest.Create( "multimeter", params={"interval": 0.1, "record_from": ["V_m", "I_syn_ampa", "I_syn_nmda", "I_syn_gaba", "I_syn"], "record_to": "memory"})
    rec_grc = nest.Create("spike_recorder")

    nest.Connect(mult_grc, grc)
    nest.Connect(grc, rec_grc)
    nest.Connect(glom_generator_1, glom)

    nest.Simulate(2200)
    mult_events = mult_grc.get()["events"]
    rec_spikes = rec_grc.get()["events"]
    times = mult_events["times"]
    I_syn_ampa = mult_events["I_syn_ampa"]
    I_syn_nmda = mult_events["I_syn_nmda"]
    I_syn_gaba = mult_events["I_syn_gaba"]
    I_syn = mult_events["I_syn"]
    V_m = mult_events["V_m"]
    g_syn_ampa = I_syn_ampa / (params_grc["AMPA_E_rev"] - V_m)
    g_syn_nmda = I_syn_nmda / (params_grc["NMDA_E_rev"] - V_m)
    g_syn_gaba = I_syn_gaba / (params_grc["GABA_E_rev"] - V_m)

    return g_syn_ampa, g_syn_nmda, g_syn_gaba, times, I_syn_ampa, I_syn_nmda, I_syn_gaba, I_syn, V_m, rec_spikes



#PLOT g BEHAVIOUR
params_grc = {'A1': 0.01,
                  'A2': -0.94,
                  'C_m': 7,
                  'E_L': -62,
                  'I_e': -0.888,
                  'V_m': -62.0,
                  'V_min': -150,
                  'V_reset': -70,
                  'V_th': -41,
                  'k_1': 0.311,
                  'k_2': 0.041407868,
                  'k_adap': 0.022,
                  'lambda_0': 1.0,
                  't_ref': 1.5,
                  'tau_V': 0.3,
                  'tau_m': 24.15,
                  # mf-GrC AMPA
                  'AMPA_E_rev': 0.0,
                  'AMPA_g_init': 1.3254051705870171,
                  'AMPA_Tau_r': 0.9966611702091053,
                  'AMPA_Tau_d1': 0.38598738306536895,
                  'AMPA_Tau_d2': 3.7021239405389736,
                  'AMPA_A_r': 6.973773179799838,
                  'AMPA_A1': 1.4551614318130006,
                  'AMPA_A2': 1.5956590886325035,
                  # mf-GrC NMDA
                  'NMDA_E_rev': -3.7,
                  'NMDA_g_init': 0.023436539864118,
                  'NMDA_Tau_r': 21.192329383729344,
                  'NMDA_Tau_d1': 14.860706054687927,
                  'NMDA_Tau_d2': 121.08581080624421,
                  'NMDA_A_r': 3.2972228905396754,
                  'NMDA_A1': 2.951290940616272,
                  'NMDA_A2': 0.21348858259846124,
                  # GoC-GrC GABA
                  'GABA_E_rev': -70,
                  'GABA_g_init': 0.16167187460716928,
                  'GABA_Tau_r': 0.37595979577132965,
                  'GABA_Tau_d1': 27.292807790043216,
                  'GABA_Tau_d2': 269.0589627801647,
                  'GABA_A_r': 4.540622636273912,
                  'GABA_A1': 10.603535445257966,
                  'GABA_A2': 3.144942430956719}
tot_spikes=[2000, 1000, 500, 200, 100, 20, 10, 5, 2]
freq=[0.5, 1, 2, 5, 10, 50, 100, 200, 500]

for i in range(len(tot_spikes)):
    g_syn_ampa, g_syn_nmda, g_syn_gaba, times, I_syn_ampa, I_syn_nmda, I_syn_gaba, I_syn, V_m, spikes = run_simulation_static(params_grc, tot_spikes[i])
    g_syn_ampa_p, g_syn_nmda_p, g_syn_gaba_p, times_p, I_syn_ampa_p, I_syn_nmda_p, I_syn_gaba_p, I_syn_p, V_m_p, spikes_p = run_simulation_stp(params_grc, tot_spikes[i])
    fig, axs = plt.subplots(4,2)
    fig.set_size_inches(20,10)
    axs[0,0].plot(times, g_syn_ampa)
    axs[0,0].set_title(f"g AMPA (static)")
    axs[0,0].set(xlabel="Time [ms]", ylabel="g [nS]")

    axs[1,0].plot(times, g_syn_nmda)
    axs[1,0].set_title(f"g NMDA (static)")
    axs[1,0].set(xlabel="Time [ms]", ylabel="g [nS]")

    axs[2, 0].plot(times, V_m)
    axs[2, 0].set_title(f"Voltage trace (static)")
    axs[2, 0].set(xlabel="Time [ms]", ylabel="V_m [mV]")

    axs[3, 0].scatter(spikes["times"], [1] * len(spikes["times"]), marker='|')
    axs[3, 0].set_title(f"Spike train (static)")
    axs[3, 0].set(xlabel="Time [ms]", ylabel="Spike")

    axs[0,1].plot(times_p, g_syn_ampa_p)
    axs[0,1].set_title(f"g AMPA (STP)")
    axs[0,1].set(xlabel="Time [ms]", ylabel="g [nS]")

    axs[1,1].plot(times_p , g_syn_nmda_p)
    axs[1,1].set_title(f"g NMDA (STP)")
    axs[1,1].set(xlabel="Time [ms]", ylabel="g [nS]")

    axs[2, 1].plot(times_p, V_m_p)
    axs[2, 1].set_title(f"Voltage trace (STP)")
    axs[2, 1].set(xlabel="Time [ms]", ylabel="V_m [mV]")

    axs[3, 1].scatter(spikes_p["times"], [1] * len(spikes_p["times"]), marker='|')
    axs[3, 1].set_title(f"Spike train (STP)")
    axs[3, 1].set(xlabel="Time [ms]", ylabel="Spike")

    fig.tight_layout()
    fig.savefig(f"static_vs_stp_{freq[i]}.png", dpi=300)
    fig.show()








############### MULTIPLE PLOTS (STATIC AND STP SEPARATED) ###############
# for i in range(len(tot_spikes)):
#     g_syn_ampa, g_syn_nmda, g_syn_gaba, times, I_syn_ampa, I_syn_nmda, I_syn_gaba, I_syn = run_simulation_static(params_grc, tot_spikes[i])
#     fig, axs = plt.subplots(3,1)
#     axs[0].plot(times, g_syn_ampa)
#     axs[0].set_title(f"AMPA conductance w/ f= {freq[i]} Hz (static)")
#     axs[0].set(xlabel="Time [ms]", ylabel="g [nS]")
#
#     axs[1].plot(times, g_syn_nmda)
#     axs[1].set_title(f"NMDA conductance w/ f= {freq[i]} Hz (static)")
#     axs[1].set(xlabel="Time [ms]", ylabel="g [nS]")
#
#     axs[2].plot(times, g_syn_gaba)
#     axs[2].set_title(f"GABA conductance w/ f= {freq[i]} Hz (static)")
#     axs[2].set(xlabel="Time [ms]", ylabel="g [nS]")
#
#     fig.tight_layout()
#     fig.show()
#
# for i in range(len(tot_spikes)):
#     g_syn_ampa, g_syn_nmda, g_syn_gaba, times, I_syn_ampa, I_syn_nmda, I_syn_gaba, I_syn = run_simulation_stp(params_grc, tot_spikes[i])
#
#     fig, axs = plt.subplots(3,1)
#     axs[0].plot(times, g_syn_ampa)
#     axs[0].set_title(f"AMPA conductance w/ f= {freq[i]} Hz (STP)")
#     axs[0].set(xlabel="Time [ms]", ylabel="g [nS]")
#
#     axs[1].plot(times, g_syn_nmda)
#     axs[1].set_title(f"NMDA conductance w/ f= {freq[i]} Hz (STP)")
#     axs[1].set(xlabel="Time [ms]", ylabel="g [nS]")
#
#     axs[2].plot(times, g_syn_gaba)
#     axs[2].set_title(f"GABA conductance w/ f= {freq[i]} Hz (STP)")
#     axs[2].set(xlabel="Time [ms]", ylabel="g [nS]")
#
#     fig.tight_layout()
#     fig.show()






############### SINGLE PLOTS ###############

# plt.figure()
# plt.plot(times, g_syn_ampa)
# plt.title("GrC AMPA conductance behaviour")
# plt.xlabel("Time [ms]")
# plt.ylabel("g [nS]")
# #plt.xlim([100,200])
# plt.show()
#
# plt.figure()
# plt.plot(times, g_syn_nmda)
# plt.title("GrC NMDA conductance behaviour")
# plt.xlabel("Time [ms]")
# plt.ylabel("g [nS]")
# plt.show()
#
# plt.figure()
# plt.plot(times, g_syn_gaba)
# plt.title("GrC GABA conductance behaviour")
# plt.xlabel("Time [ms]")
# plt.ylabel("g [nS]")
# plt.show()

# plt.figure()
# plt.plot(times, I_syn)
# plt.title("Total current behaviour")
# plt.xlabel("Time [ms]")
# plt.ylabel("I_syn_nmda [pA]")
# #plt.xlim([100,200])
# plt.show()