import nest
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

def single_sim(cell_params, syn_spec, spike_times):
    nest.ResetKernel()
    nest.Install('cerebmodule')
    nest.resolution = 0.025
    grc = nest.Create('eglif_multisyn', 1, params=cell_params)
    input_spikes = nest.Create("spike_generator", params={"spike_times": spike_times})
    nest.Connect(input_spikes, grc, syn_spec=dict(syn_spec, synapse_model="static_synapse"))

    multimeter = nest.Create("multimeter", params={
        "interval": 0.025,
        "record_from": ['V_m', 'I_syn_nmda', 'Mg_block'],
        "record_to": "memory",
        "label":'grc_multimeter'
    })

    spike_recorder = nest.Create("spike_recorder", params={
        "record_to": "memory",
        "label": 'grc_recorder'
    })

    nest.Connect(multimeter, grc)
    nest.Connect(grc, spike_recorder)
    duration = 1000.
    nest.Simulate(duration)

    multimeter_data = multimeter.get()['events']
    times = multimeter_data['times']
    I_syn_nmda = multimeter_data['I_syn_nmda']
    V_m = multimeter_data["V_m"]
    Mg_block = multimeter_data["Mg_block"]
    g_syn_nmda1 = I_syn_nmda / (cell_params['NMDA_E_rev'] - V_m)
    g_syn_nmda2 = g_syn_nmda1/Mg_block

    return times, g_syn_nmda1, g_syn_nmda2, V_m, Mg_block, I_syn_nmda


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
    "NMDA_E_rev": -3.7
}

k = 10. # to ensure at least a spike will occur
syn_spec = {"weight": k * 1. , "delay": 1, "receptor_type": 2}

# spike_times = list(np.arange(1,400,50))
spike_times = [50.]
times, g_syn1, g_syn2, V_m, Mg_block, I_syn = single_sim(cell_params, syn_spec, spike_times)

# verify Mg unblock duration
mg_max = np.max(Mg_block)
is_max = (Mg_block == mg_max)
t_start = times[np.argmax(is_max)]
t_end = times[len(is_max) - np.argmax(is_max[::-1]) - 1]
t_Mgunblock = (t_end - t_start)


fig, axs = plt.subplots(5, 1, figsize=(8, 12), sharex=True)

axs[0].set_title('NMDA conductance (/ MgBlock)')
axs[0].plot(times, g_syn2, 'b', label='NEST trace')
axs[0].set_ylabel(r'$g_{syn_{NMDA}}$ [ns]')
axs[0].set_xlim(0, 200)

axs[1].set_title('NMDA conductance (with MgBlock)')
axs[1].plot(times, g_syn1, 'b', label='NEST trace')
axs[1].set_ylabel(r'$g_{syn_{NMDA}}$ [ns]')
axs[1].set_xlim(0, 200)

axs[2].set_title('Membrane voltage')
axs[2].plot(times, V_m, 'black',label='NEST trace')
axs[2].axhline(cell_params['V_th'], color='r', linestyle='--', label=r'$V_{th}$')
axs[2].set_ylabel(r'$V_m$ [mV]')
axs[2].set_xlim(0, 200)


axs[3].set_title('NMDA current')
axs[3].plot(times, I_syn, 'red', label='NEST trace')
axs[3].set_ylabel(r'$I_{syn_{NMDA}}$ [pA]')
axs[3].set_xlim(0, 200)

axs[4].set_title('Mg block')
axs[4].plot(times, Mg_block, 'green',label='NEST trace')
axs[4].set_ylabel(r'Mg_Block')
axs[4].set_xlim(0, 200)
axs[4].set_xlabel('Time [ms]')

axins = inset_axes(axs[4], width="30%", height="40%", loc='center right', borderpad=2)
t_zoom_center = t_start + (t_end - t_start) / 2
t_zoom_width = 2.5
x1, x2 = t_zoom_center - t_zoom_width, t_zoom_center + t_zoom_width
axins.plot(times, Mg_block, 'green')
axins.set_xlim(x1, x2)
axins.set_ylim(0.9 * mg_max, 1.02 * mg_max)
axins.set_xticks([])
axins.set_yticks([])
max_idx = np.argmax(Mg_block)
axins.text(times[max_idx], Mg_block[max_idx] + 0.002,
           f't_Mgunblock = {t_Mgunblock:.2f} ms',
           fontsize=7, ha='center', va='bottom', color='black')

mark_inset(axs[4], axins, loc1=3, loc2=1, fc="none", ec="0.1", alpha=0.5)

for i, ax in enumerate(axs):
    for j, spike in enumerate(spike_times):
        label = 'input spikes' if i == 0 and j == 0 else None
        ax.axvline(spike, color='gray', linestyle='--', linewidth=0.5, label=label)

for ax in axs[:-1]:
    ax.tick_params(labelbottom=False, bottom=False)

axs[0].legend()
axs[1].legend()
axs[2].legend()
axs[3].legend()
axs[4].legend()

plt.setp(axs[-1].get_xticklabels(), visible=True)
plt.tight_layout()
plt.savefig('./figs/Mg_block_test.png', dpi=300)
plt.show()



