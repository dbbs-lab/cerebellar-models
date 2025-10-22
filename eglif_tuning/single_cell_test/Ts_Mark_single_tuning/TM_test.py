import yaml
import rec_functions as rec

if __name__ == "__main__":

    # Load yaml files necessary for neuron and synapse models/parameters
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

    # Structure all models into a clean dict variable
    for name,parameters in configuration.items():
        if "glomerulus" in name:
            cell_params[name] = parameters["model"]
        elif "mossy_fibers" in name:
            cell_params[name] = parameters["model"]
        else:
            cell_params[name] = parameters["constants"]

    for name,parameters in config_stp.items():
        stp_params[name] = parameters["synapse"]
        del stp_params[name]["model"]

    # Input firing rates
    tot_spikes=[2000, 1000, 500, 200, 100, 20, 10, 5, 2]
    freq=[0.5, 1, 2, 5, 10, 50, 100, 200, 500]
    f_rates_static = []
    f_rates_stp = []

    # run the simulations for synapse testing (static vs. Tsodyks-Markram Short-Term Plasticity)
    rec.run_and_plot(cell_params["glomerulus"], cell_params["granule_cell"], stp_params["glomerulus_to_granule"],
                     1, tot_spikes , freq, False, False, f_rates_static, f_rates_stp)


    # for i in range(len(tot_spikes)):
    #     g_syn_ampa, g_syn_nmda, g_syn_gaba, times, I_syn_ampa, I_syn_nmda, I_syn_gaba, I_syn, V_m, spikes, fr_static = rec.sim_static_grc(cell_params["granule_cell"], tot_spikes[i])
    #     g_syn_ampa_p, g_syn_nmda_p, g_syn_gaba_p, times_p, I_syn_ampa_p, I_syn_nmda_p, I_syn_gaba_p, I_syn_p, V_m_p, spikes_p, fr_p = rec.sim_stp_grc(cell_params["granule_cell"], tot_spikes[i])
    #     fig, axs = plt.subplots(4,2)
    #     fig.set_size_inches(20,10)
    #     axs[0,0].plot(times, g_syn_ampa)
    #     axs[0,0].set_title(f"g AMPA (static)")
    #     axs[0,0].set(xlabel="Time [ms]", ylabel="g [nS]")
    #
    #     axs[1,0].plot(times, g_syn_nmda)
    #     axs[1,0].set_title(f"g NMDA (static)")
    #     axs[1,0].set(xlabel="Time [ms]", ylabel="g [nS]")
    #
    #     axs[2, 0].plot(times, V_m)
    #     axs[2, 0].set_title(f"Voltage trace (static)")
    #     axs[2, 0].set(xlabel="Time [ms]", ylabel="V_m [mV]")
    #
    #     axs[3, 0].scatter(spikes["times"], [1] * len(spikes["times"]), marker='|')
    #     axs[3, 0].set_title(f"Spike train (static)")
    #     axs[3, 0].set(xlabel="Time [ms]", ylabel="Spike")
    #
    #     axs[0,1].plot(times_p, g_syn_ampa_p)
    #     axs[0,1].set_title(f"g AMPA (STP)")
    #     axs[0,1].set(xlabel="Time [ms]", ylabel="g [nS]")
    #
    #     axs[1,1].plot(times_p , g_syn_nmda_p)
    #     axs[1,1].set_title(f"g NMDA (STP)")
    #     axs[1,1].set(xlabel="Time [ms]", ylabel="g [nS]")
    #
    #     axs[2, 1].plot(times_p, V_m_p)
    #     axs[2, 1].set_title(f"Voltage trace (STP)")
    #     axs[2, 1].set(xlabel="Time [ms]", ylabel="V_m [mV]")
    #
    #     axs[3, 1].scatter(spikes_p["times"], [1] * len(spikes_p["times"]), marker='|')
    #     axs[3, 1].set_title(f"Spike train (STP)")
    #     axs[3, 1].set(xlabel="Time [ms]", ylabel="Spike")
    #
    #     fig.tight_layout()
    #     fig.savefig(f"imgs/static_vs_stp_{freq[i]}.png", dpi=300)
    #     fig.show()
    #
    #     f_rates_static.append(fr_static)
    #     f_rates_stp.append(fr_p)
    #
    #
    # fig, axs = plt.subplots()
    # axs.plot(freq, f_rates_static, marker='o')
    # axs.plot(freq, f_rates_stp, marker='o')
    # axs.legend(["f-static", "f-stp"], loc="lower right")
    # # plt.xlim([0,10])
    # # plt.ylim([0,25])
    # axs.set_title("f-f Curve Glom_to_Grc")
    # axs.set_xlabel("Frequency Input [Hz]")
    # axs.set_ylabel("Frequency Output [Hz]")
    # axins = inset_axes(axs, width="30%", height="40%", loc="upper left")
    # t_zoom_center = 5
    # t_zoom_width = 25
    #
    # axins.plot(freq, f_rates_static, marker='o')
    # axins.plot(freq, f_rates_stp, marker='o')
    # axins.yaxis.tick_right()
    # axins.set_xlim(-3, 20)
    # axins.set_ylim(-0.5, 8)
    # mark_inset(axs, axins, loc1=3, loc2=4, fc="none", ec="0.1", alpha=0.5)
    # plt.show()
    # fig.savefig(f"imgs/static_vs_stp_f.png", dpi=300)