import numpy as np
import nest
import matplotlib.pyplot as plt
import quantities as pq
from neo.core import SpikeTrain
import elephant.statistics as es
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

# Specific for Glomerulus to Granule Cell synapse
def sim_static_grc(params_grc, tot_spikes):
    nest.ResetKernel()
    nest.Install("cerebmodule")

    glom = nest.Create("parrot_neuron", 1)
    grc = nest.Create("eglif_multirec_opt", 1, params=params_grc)

    sp_times = np.arange(1, 2001, tot_spikes, dtype=float)
    glom_generator_1 = nest.Create("spike_generator", params={"spike_times": sp_times})


    nest.Connect(glom, grc, syn_spec={"synapse_model": "static_synapse", "weight": 1, "delay": 0.1, "receptor_type": 1})
    nest.Connect(glom, grc, syn_spec={"synapse_model": "static_synapse", "weight": 1, "delay": 0.1, "receptor_type": 2})
    # nest.Connect(glom, grc, syn_spec={"synapse_model": "static_synapse", "weight": 1, "delay": 1, "receptor_type": 3})
    mult_grc = nest.Create("multimeter", params={"interval": 0.1,
                                                 "record_from": ["V_m", "I_syn_ampa", "I_syn_nmda", "I_syn_gaba",
                                                                 "I_syn"], "record_to": "memory"})
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
    if "NMDA_E_rev" in params_grc:
        g_syn_nmda = I_syn_nmda / (params_grc["NMDA_E_rev"] - V_m)
    else:
        g_syn_nmda = 0
    g_syn_gaba = I_syn_gaba / (params_grc["GABA_E_rev"] - V_m)

    sp_times = rec_grc.get()["events"]["times"]
    sp_times = sp_times[sp_times <= 2001]
    if sp_times is None:
        firing_rate = 0
    else:
        st = SpikeTrain(sp_times * pq.ms, t_start=1 * pq.ms, t_stop=2001 * pq.ms)
        firing_rate = es.mean_firing_rate(st).rescale(pq.Hz)

    return g_syn_ampa, g_syn_nmda, g_syn_gaba, times, I_syn_ampa, I_syn_nmda, I_syn_gaba, I_syn, V_m, rec_spikes, firing_rate


# Specific for Glomerulus to Granule Cell synapse
def sim_stp_grc(params_grc, tot_spikes):
    nest.ResetKernel()
    nest.Install("cerebmodule")

    glom = nest.Create("parrot_neuron", 1)
    grc = nest.Create("eglif_multirec_opt", 1, params=params_grc)

    sp_times = np.arange(1, 2001, tot_spikes, dtype=float)
    glom_generator_1 = nest.Create("spike_generator", params={"spike_times": sp_times})
    # glom_generator_2 = nest.Create("poisson_generator", params={"rate": 4, "start": 0, "stop": 1000})

    nest.Connect(glom, grc, syn_spec={"synapse_model": "tsodyks2_synapse", "weight": 2, "U": 0.43, "x": 1, "tau_rec": 8,
                                      "tau_fac": 5, "delay": 1, "receptor_type": 1})
    nest.Connect(glom, grc, syn_spec={"synapse_model": "tsodyks2_synapse", "weight": 2, "U": 0.43, "x": 1, "tau_rec": 8,
                                      "tau_fac": 5, "delay": 1, "receptor_type": 2})

    mult_grc = nest.Create("multimeter", params={"interval": 0.1,
                                                 "record_from": ["V_m", "I_syn_ampa", "I_syn_nmda", "I_syn_gaba",
                                                                 "I_syn"], "record_to": "memory"})
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
    if "NMDA_E_rev" in params_grc:
        g_syn_nmda = I_syn_nmda / (params_grc["NMDA_E_rev"] - V_m)
    else:
        g_syn_nmda = 0
    g_syn_gaba = I_syn_gaba / (params_grc["GABA_E_rev"] - V_m)

    sp_times = rec_grc.get()["events"]["times"]
    sp_times = sp_times[sp_times <= 2001]
    if sp_times is None:
        firing_rate = 0
    else:
        st = SpikeTrain(sp_times * pq.ms, t_start=1 * pq.ms, t_stop=2001 * pq.ms)
        firing_rate = es.mean_firing_rate(st).rescale(pq.Hz)

    return g_syn_ampa, g_syn_nmda, g_syn_gaba, times, I_syn_ampa, I_syn_nmda, I_syn_gaba, I_syn, V_m, rec_spikes, firing_rate



def sim_static(params_pre, params_post, weight, tot_spikes, aa_to_goc: bool, pf_to_goc: bool):
    nest.ResetKernel()
    nest.Install("cerebmodule")

    pre = nest.Create("eglif_multirec_opt", 1, params=params_pre)
    post = nest.Create("eglif_multirec_opt", 1, params=params_post)

    sp_times = np.arange(1, 2001, tot_spikes, dtype=float)
    sp_generator = nest.Create("spike_generator", params={"spike_times": sp_times})

    if "AMPA_E_rev" in params_post:
        nest.Connect(pre, post,
                     syn_spec={"synapse_model": "static_synapse", "weight": weight, "delay": 0.1, "receptor_type": 1})
    if "NMDA_E_rev" in params_post:
        nest.Connect(pre, post,
                     syn_spec={"synapse_model": "static_synapse", "weight": weight, "delay": 0.1, "receptor_type": 2})
    if "GABA_E_rev" in params_post:
        nest.Connect(pre, post,
                     syn_spec={"synapse_model": "static_synapse", "weight": weight, "delay": 0.1, "receptor_type": 3})
    if aa_to_goc:
        nest.Connect(pre, post,
                     syn_spec={"synapse_model": "static_synapse", "weight": weight, "delay": 0.1, "receptor_type": 4})
    if pf_to_goc:
        nest.Connect(pre, post,
                     syn_spec={"synapse_model": "static_synapse", "weight": weight, "delay": 0.1, "receptor_type": 5})

    mult = nest.Create("multimeter", params={"interval": 0.1,
                                             "record_from": ["V_m", "I_syn_ampa", "I_syn_ampa2", "I_syn_ampa3",
                                                             "I_syn_nmda", "I_syn_gaba", "I_syn"],
                                             "record_to": "memory"})

    rec = nest.Create("spike_recorder")

    nest.Connect(mult, post)
    nest.Connect(post, rec)
    nest.Connect(sp_generator, pre)

    nest.Simulate(2200)
    mult_events = mult.get()["events"]
    rec_spikes = rec.get()["events"]
    times = mult_events["times"]
    I_syn_ampa = mult_events["I_syn_ampa"]
    I_syn_ampa2 = mult_events["I_syn_ampa2"]
    I_sym_ampa3 = mult_events["I_sym_ampa3"]
    I_syn_nmda = mult_events["I_syn_nmda"]
    I_syn_gaba = mult_events["I_syn_gaba"]
    I_syn = mult_events["I_syn"]
    V_m = mult_events["V_m"]

    if "AMPA_E_rev" in params_post:
        g_syn_ampa = I_syn_ampa / (params_post["AMPA_E_rev"] - V_m)
    else:
        g_syn_ampa = 0

    if "AMPA2_E_rev" in params_post:
        g_syn_ampa2 = I_syn_ampa2 / (params_post["AMPA2_E_rev"] - V_m)
    else:
        g_syn_ampa2 = 0

    if "AMPA3_E_rev" in params_post:
        g_syn_ampa3 = I_sym_ampa3 / (params_post["AMPA3_E_rev"] - V_m)
    else:
        g_syn_ampa3 = 0

    if "NMDA_E_rev" in params_post:
        g_syn_nmda = I_syn_nmda / (params_post["NMDA_E_rev"] - V_m)
    else:
        g_syn_nmda = 0

    if "GABA_E_rev" in params_post:
        g_syn_gaba = I_syn_gaba / (params_post["GABA_E_rev"] - V_m)
    else:
        g_syn_gaba = 0

    sp_times = rec.get()["events"]["times"]
    sp_times = sp_times[sp_times <= 2001]
    if sp_times is None:
        firing_rate = 0
    else:
        st = SpikeTrain(sp_times * pq.ms, t_start=1 * pq.ms, t_stop=2001 * pq.ms)
        firing_rate = es.mean_firing_rate(st).rescale(pq.Hz)

    return g_syn_ampa, g_syn_ampa2, g_syn_ampa3, g_syn_nmda, g_syn_gaba, times, I_syn_ampa, I_syn_nmda, I_syn_gaba, I_syn, V_m, rec_spikes, firing_rate



def sim_stp(params_pre, params_post, params_syn, tot_spikes, aa_to_goc: bool, pf_to_goc: bool):
    nest.ResetKernel()
    nest.Install("cerebmodule")

    pre = nest.Create("eglif_multirec_opt", 1, params=params_pre)
    post = nest.Create("eglif_multirec_opt", 1, params=params_post)

    sp_times = np.arange(1, 2001, tot_spikes, dtype=float)
    sp_generator = nest.Create("spike_generator", params={"spike_times": sp_times})

    if "AMPA_E_rev" in params_post:
        nest.Connect(pre, post, syn_spec={"synapse_model": "tsodyks2_synapse", "weight": params_syn["weight"],
                                          "U": params_syn["U"], "x": params_syn["x"], "tau_rec": params_syn["tau_rec"],
                                          "tau_fac": params_syn["tau_fac"], "delay": 0.1, "receptor_type": 1})
    if "NMDA_E_rev" in params_post:
        nest.Connect(pre, post, syn_spec={"synapse_model": "tsodyks2_synapse", "weight": params_syn["weight"],
                                          "U": params_syn["U"], "x": params_syn["x"], "tau_rec": params_syn["tau_rec"],
                                          "tau_fac": params_syn["tau_fac"], "delay": 0.1, "receptor_type": 2})
    if "GABA_E_rev" in params_post:
        nest.Connect(pre, post, syn_spec={"synapse_model": "tsodyks2_synapse", "weight": params_syn["weight"],
                                          "U": params_syn["U"], "x": params_syn["x"], "tau_rec": params_syn["tau_rec"],
                                          "tau_fac": params_syn["tau_fac"], "delay": 0.1, "receptor_type": 3})
    if aa_to_goc:
        nest.Connect(pre, post, syn_spec={"synapse_model": "tsodyks2_synapse", "weight": params_syn["weight"],
                                          "U": params_syn["U"], "x": params_syn["x"], "tau_rec": params_syn["tau_rec"],
                                          "tau_fac": params_syn["tau_fac"], "delay": 0.1, "receptor_type": 4})
    if pf_to_goc:
        nest.Connect(pre, post, syn_spec={"synapse_model": "tsodyks2_synapse", "weight": params_syn["weight"],
                                          "U": params_syn["U"], "x": params_syn["x"], "tau_rec": params_syn["tau_rec"],
                                          "tau_fac": params_syn["tau_fac"], "delay": 0.1, "receptor_type": 5})

    mult = nest.Create("multimeter", params={"interval": 0.1,
                                             "record_from": ["V_m", "I_syn_ampa", "I_syn_ampa2", "I_syn_ampa3",
                                                             "I_syn_nmda", "I_syn_gaba", "I_syn"],
                                             "record_to": "memory"})

    rec = nest.Create("spike_recorder")

    nest.Connect(mult, post)
    nest.Connect(post, rec)
    nest.Connect(sp_generator, pre)

    nest.Simulate(2200)
    mult_events = mult.get()["events"]
    rec_spikes = rec.get()["events"]
    times = mult_events["times"]
    I_syn_ampa = mult_events["I_syn_ampa"]
    I_syn_ampa2 = mult_events["I_syn_ampa2"]
    I_sym_ampa3 = mult_events["I_sym_ampa3"]
    I_syn_nmda = mult_events["I_syn_nmda"]
    I_syn_gaba = mult_events["I_syn_gaba"]
    I_syn = mult_events["I_syn"]
    V_m = mult_events["V_m"]

    if "AMPA_E_rev" in params_post:
        g_syn_ampa = I_syn_ampa / (params_post["AMPA_E_rev"] - V_m)
    else:
        g_syn_ampa = 0

    if "AMPA2_E_rev" in params_post:
        g_syn_ampa2 = I_syn_ampa2 / (params_post["AMPA2_E_rev"] - V_m)
    else:
        g_syn_ampa2 = 0

    if "AMPA3_E_rev" in params_post:
        g_syn_ampa3 = I_sym_ampa3 / (params_post["AMPA3_E_rev"] - V_m)
    else:
        g_syn_ampa3 = 0

    if "NMDA_E_rev" in params_post:
        g_syn_nmda = I_syn_nmda / (params_post["NMDA_E_rev"] - V_m)
    else:
        g_syn_nmda = 0

    if "GABA_E_rev" in params_post:
        g_syn_gaba = I_syn_gaba / (params_post["GABA_E_rev"] - V_m)
    else:
        g_syn_gaba = 0

    sp_times = rec.get()["events"]["times"]
    sp_times = sp_times[sp_times <= 2001]
    if sp_times is None:
        firing_rate = 0
    else:
        st = SpikeTrain(sp_times * pq.ms, t_start=1 * pq.ms, t_stop=2001 * pq.ms)
        firing_rate = es.mean_firing_rate(st).rescale(pq.Hz)

    return g_syn_ampa, g_syn_ampa2, g_syn_ampa3, g_syn_nmda, g_syn_gaba, times, I_syn_ampa, I_syn_nmda, I_syn_gaba, I_syn, V_m, rec_spikes, firing_rate


# Makes the specified simulations, plots and saves the plots into png images
def run_and_plot(params_pre, params_post, params_syn, weight, tot_spikes, freq, aa_to_goc: bool, pf_to_goc: bool, f_rates_static, f_rates_stp):
    for i in range(len(tot_spikes)):
        if 'parrot_neuron' in params_pre:
            g_syn_ampa, g_syn_nmda, g_syn_gaba, times, I_syn_ampa, I_syn_nmda, I_syn_gaba, I_syn, V_m, spikes, fr_static = sim_static_grc(params_post, tot_spikes[i])
            g_syn_ampa_p, g_syn_nmda_p, g_syn_gaba_p, times_p, I_syn_ampa_p, I_syn_nmda_p, I_syn_gaba_p, I_syn_p, V_m_p, spikes_p, fr_p = sim_stp_grc(params_post, tot_spikes[i])
        else:
            g_syn_ampa, g_syn_ampa2, g_syn_ampa3, g_syn_nmda, g_syn_gaba, times, I_syn_ampa, I_syn_nmda, I_syn_gaba, I_syn, V_m, spikes, fr_static = sim_static(
                params_pre, params_post, weight, tot_spikes[i], aa_to_goc, pf_to_goc)
            g_syn_ampa_p, g_syn_ampa2_p, g_syn_ampa3_p, g_syn_nmda_p, g_syn_gaba_p, times_p, I_syn_ampa_p, I_syn_nmda_p, I_syn_gaba_p, I_syn_p, V_m_p, spikes_p, fr_p = sim_stp(
                params_pre, params_post, params_syn, tot_spikes[i], aa_to_goc, pf_to_goc)

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
        fig.savefig(f"imgs/static_vs_stp_{freq[i]}.png", dpi=300)
        fig.show()

        f_rates_static.append(fr_static)
        f_rates_stp.append(fr_p)



    fig, axs = plt.subplots()
    axs.plot(freq, f_rates_static, marker='o')
    axs.plot(freq, f_rates_stp, marker='o')
    axs.legend(["f-static", "f-stp"], loc="lower right")
    # plt.xlim([0,10])
    # plt.ylim([0,25])
    axs.set_title("f-f Curve Glom_to_Grc")
    axs.set_xlabel("Frequency Input [Hz]")
    axs.set_ylabel("Frequency Output [Hz]")
    axins = inset_axes(axs, width="30%", height="40%", loc="upper left")
    t_zoom_center = 5
    t_zoom_width = 25

    axins.plot(freq, f_rates_static, marker='o')
    axins.plot(freq, f_rates_stp, marker='o')
    axins.yaxis.tick_right()
    axins.set_xlim(-3, 20)
    axins.set_ylim(-0.5, 8)
    mark_inset(axs, axins, loc1=3, loc2=4, fc="none", ec="0.1", alpha=0.5)
    plt.show()
    fig.savefig(f"imgs/static_vs_stp_f.png", dpi=300)



# params_grc = {'A1': 0.01,
#                   'A2': -0.94,
#                   'C_m': 7,
#                   'E_L': -62,
#                   'I_e': -0.888,
#                   'V_m': -62.0,
#                   'V_min': -150,
#                   'V_reset': -70,
#                   'V_th': -41,
#                   'k_1': 0.311,
#                   'k_2': 0.041407868,
#                   'k_adap': 0.022,
#                   'lambda_0': 1.0,
#                   't_ref': 1.5,
#                   'tau_V': 0.3,
#                   'tau_m': 24.15,
#                   # mf-GrC AMPA
#                   'AMPA_E_rev': 0.0,
#                   'AMPA_g_init': 1.3254051705870171,
#                   'AMPA_Tau_r': 0.9966611702091053,
#                   'AMPA_Tau_d1': 0.38598738306536895,
#                   'AMPA_Tau_d2': 3.7021239405389736,
#                   'AMPA_A_r': 6.973773179799838,
#                   'AMPA_A1': 1.4551614318130006,
#                   'AMPA_A2': 1.5956590886325035,
#                   # mf-GrC NMDA
#                   'NMDA_E_rev': -3.7,
#                   'NMDA_g_init': 0.023436539864118,
#                   'NMDA_Tau_r': 21.192329383729344,
#                   'NMDA_Tau_d1': 14.860706054687927,
#                   'NMDA_Tau_d2': 121.08581080624421,
#                   'NMDA_A_r': 3.2972228905396754,
#                   'NMDA_A1': 2.951290940616272,
#                   'NMDA_A2': 0.21348858259846124,
#                   # GoC-GrC GABA
#                   'GABA_E_rev': -70,
#                   'GABA_g_init': 0.16167187460716928,
#                   'GABA_Tau_r': 0.37595979577132965,
#                   'GABA_Tau_d1': 27.292807790043216,
#                   'GABA_Tau_d2': 269.0589627801647,
#                   'GABA_A_r': 4.540622636273912,
#                   'GABA_A1': 10.603535445257966,
#                   'GABA_A2': 3.144942430956719}