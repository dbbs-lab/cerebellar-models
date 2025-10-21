import numpy as np
import nest
import matplotlib.pyplot as plt
import quantities as pq
from neo.core import SpikeTrain
import elephant.statistics as es
import yaml

def run_sim_static_grc(params_grc, tot_spikes):
    nest.ResetKernel()
    nest.Install("cerebmodule")

    glom = nest.Create("parrot_neuron", 1)
    grc = nest.Create("eglif_multirec_opt", 1, params=params_grc)

    sp_times = np.arange(1, 2001, tot_spikes, dtype=float)
    glom_generator_1 = nest.Create("spike_generator", params={"spike_times": sp_times})
    # glom_generator_2 = nest.Create("poisson_generator", params={"rate": 4, "start": 0, "stop": 1000})

    nest.Connect(glom, grc, syn_spec={"synapse_model": "static_synapse", "weight": 1, "delay": 0.1, "receptor_type" : 1})
    nest.Connect(glom, grc, syn_spec={"synapse_model": "static_synapse", "weight": 1, "delay": 0.1, "receptor_type": 2})
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

def run_sim_stp_grc(params_grc, tot_spikes):
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

def run_sim_static(params_pre, params_post, tot_spikes):
    nest.ResetKernel()
    nest.Install("cerebmodule")

    pre = nest.Create("eglif_multirec_opt", 1, params=params_pre)
    post = nest.Create("eglif_multirec_opt", 1, params=params_post)

    sp_times = np.arange(1, 2001, tot_spikes, dtype=float)
    sp_generator = nest.Create("spike_generator", params={"spike_times": sp_times})

    nest.Connect(pre, post, syn_spec={"synapse_model": "static_synapse", "weight": 1, "delay": 1, "receptor_type": 1})
    nest.Connect(pre, post, syn_spec={"synapse_model": "static_synapse", "weight": 1, "delay": 1, "receptor_type": 2})
    nest.Connect(pre, post, syn_spec={"synapse_model": "static_synapse", "weight": 1, "delay": 1, "receptor_type": 3})

    mult = nest.Create("multimeter", params={"interval": 0.1,"record_from": ["V_m", "I_syn_ampa", "I_syn_nmda", "I_syn_gaba","I_syn"], "record_to": "memory"})
    rec = nest.Create("spike_recorder")

    nest.Connect(mult, post)
    nest.Connect(post, rec)
    nest.Connect(sp_generator, pre)

    nest.Simulate(2200)
    mult_events = mult.get()["events"]
    rec_spikes = rec.get()["events"]
    times = mult_events["times"]
    I_syn_ampa = mult_events["I_syn_ampa"]
    I_syn_nmda = mult_events["I_syn_nmda"]
    I_syn_gaba = mult_events["I_syn_gaba"]
    I_syn = mult_events["I_syn"]
    V_m = mult_events["V_m"]
    g_syn_ampa = I_syn_ampa / (params_post["AMPA_E_rev"] - V_m)
    if "NMDA_E_rev" in params_post:
        g_syn_nmda = I_syn_nmda / (params_post["NMDA_E_rev"] - V_m)
    else:
        g_syn_nmda = 0
    g_syn_gaba = I_syn_gaba / (params_post["GABA_E_rev"] - V_m)

    sp_times = rec.get()["events"]["times"]
    sp_times = sp_times[sp_times <= 2001]
    if sp_times is None:
        firing_rate = 0
    else:
        st = SpikeTrain(sp_times * pq.ms, t_start=1 * pq.ms, t_stop=2001 * pq.ms)
        firing_rate = es.mean_firing_rate(st).rescale(pq.Hz)

    return g_syn_ampa, g_syn_nmda, g_syn_gaba, times, I_syn_ampa, I_syn_nmda, I_syn_gaba, I_syn, V_m, rec_spikes, firing_rate

def run_sum_stp(params_pre, params_post, tot_spikes):
    nest.ResetKernel()
    nest.Install("cerebmodule")

    pre = nest.Create("eglif_multirec_opt", 1, params=params_pre)
    post = nest.Create("eglif_multirec_opt", 1, params=params_post)

    sp_times = np.arange(1, 2001, tot_spikes, dtype=float)
    sp_generator = nest.Create("spike_generator", params={"spike_times": sp_times})

    nest.Connect(pre, post, syn_spec={"synapse_model": "tsodyks2_synapse", "weight": 2, "U": 0.43, "x": 1, "tau_rec": 8,"tau_fac": 5, "delay": 1, "receptor_type": 1})
    nest.Connect(pre, post, syn_spec={"synapse_model": "tsodyks2_synapse", "weight": 2, "U": 0.43, "x": 1, "tau_rec": 8,"tau_fac": 5, "delay": 1, "receptor_type": 2})
    nest.Connect(pre, post, syn_spec={"synapse_model": "tsodyks2_synapse", "weight": 2, "U": 0.43, "x": 1, "tau_rec": 8,"tau_fac": 5, "delay": 1, "receptor_type": 3})

    mult = nest.Create("multimeter", params={"interval": 0.1,"record_from": ["V_m", "I_syn_ampa", "I_syn_nmda", "I_syn_gaba", "I_syn"],"record_to": "memory"})
    rec = nest.Create("spike_recorder")

    nest.Connect(mult, post)
    nest.Connect(post, rec)
    nest.Connect(sp_generator, pre)

    nest.Simulate(2200)
    mult_events = mult.get()["events"]
    rec_spikes = rec.get()["events"]
    times = mult_events["times"]
    I_syn_ampa = mult_events["I_syn_ampa"]
    I_syn_nmda = mult_events["I_syn_nmda"]
    I_syn_gaba = mult_events["I_syn_gaba"]
    I_syn = mult_events["I_syn"]
    V_m = mult_events["V_m"]
    g_syn_ampa = I_syn_ampa / (params_post["AMPA_E_rev"] - V_m)
    if "NMDA_E_rev" in params_post:
        g_syn_nmda = I_syn_nmda / (params_post["NMDA_E_rev"] - V_m)
    else:
        g_syn_nmda = 0
    g_syn_gaba = I_syn_gaba / (params_post["GABA_E_rev"] - V_m)

    sp_times = rec.get()["events"]["times"]
    sp_times = sp_times[sp_times <= 2001]
    if sp_times is None:
        firing_rate = 0
    else:
        st = SpikeTrain(sp_times * pq.ms, t_start=1 * pq.ms, t_stop=2001 * pq.ms)
        firing_rate = es.mean_firing_rate(st).rescale(pq.Hz)

    return g_syn_ampa, g_syn_nmda, g_syn_gaba, times, I_syn_ampa, I_syn_nmda, I_syn_gaba, I_syn, V_m, rec_spikes, firing_rate


if __name__ == "__main__":

    static_path = "../../../configurations/mouse/nest_examples/basal_vitro_rec.yaml"
    stp_path = "../../../configurations/mouse/in-vitro/nest/connection_models/tsodyks2_synapse.yaml"

    with open(static_path) as stream_s:
        configuration = yaml.safe_load(stream_s)
        configuration = configuration["simulations"]["nest_basal_activity"]["cell_models"]

    with open(stp_path) as stream:
        config_stp = yaml.safe_load(stream)
        config_stp = config_stp["simulations"]["basal_activity"]["connection_models"]

    cell_params = {}
    stp_params = {}
    static_params = {}
    for name,parameters in configuration.items():
        if "glomerulus" in name:
            cell_params[name] = parameters["model"]
        elif "mossy_fibers" in name:
            cell_params[name] = parameters["model"]
        else:
            cell_params[name] = parameters["constants"]

    # for name, parameters in configuration.simulations.nest_basal_activity.cell_models.items():

    for name,parameters in config_stp.items():
        stp_params[name] = parameters["synapse"]
        del stp_params[name]["model"]

    tot_spikes=[2000, 1000, 500, 200, 100, 20, 10, 5, 2]
    freq=[0.5, 1, 2, 5, 10, 50, 100, 200, 500]
    f_rates_static = []
    f_rates_stp = []

    for i in range(len(tot_spikes)):
        g_syn_ampa, g_syn_nmda, g_syn_gaba, times, I_syn_ampa, I_syn_nmda, I_syn_gaba, I_syn, V_m, spikes, fr_static = run_sim_static_grc(cell_params["granule_cell"], tot_spikes[i])
        g_syn_ampa_p, g_syn_nmda_p, g_syn_gaba_p, times_p, I_syn_ampa_p, I_syn_nmda_p, I_syn_gaba_p, I_syn_p, V_m_p, spikes_p, fr_p = run_sim_stp_grc(cell_params["granule_cell"], tot_spikes[i])
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


    fig = plt.figure()
    plt.plot(freq, f_rates_static, marker='o')
    plt.plot(freq, f_rates_stp, marker='o')
    plt.legend(["f-static", "f-stp"])
    # plt.xlim([0,10])
    # plt.ylim([0,25])
    plt.title("f-f Curve Glom_to_Grc")
    plt.xlabel("Frequency Input [Hz]")
    plt.ylabel("Frequency Output [Hz]")
    plt.show()
    fig.savefig(f"imgs/static_vs_stp_f.png", dpi=300)



# axins = inset_axes(axs[4], width="30%", height="40%", loc="center right", borderpad=2)
# t_zoom_center = t_start + (t_end - t_start) / 2
# t_zoom_width = 2.5
# x1, x2 = t_zoom_center - t_zoom_width, t_zoom_center + t_zoom_width
# axins.plot(times, Mg_block, "green")
# axins.set_xlim(x1, x2)
# axins.set_ylim(0.9 * mg_max, 1.02 * mg_max)
# axins.set_xticks([])
# axins.set_yticks([])
# max_idx = np.argmax(Mg_block)
# axins.text(
#     times[max_idx],
#     Mg_block[max_idx] + 0.002,
#     f"t_Mgunblock = {t_Mgunblock:.2f} ms",
#     fontsize=7,
#     ha="center",
#     va="bottom",
#     color="black",
# )
#
# mark_inset(axs[4], axins, loc1=3, loc2=1, fc="none", ec="0.1", alpha=0.5)



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