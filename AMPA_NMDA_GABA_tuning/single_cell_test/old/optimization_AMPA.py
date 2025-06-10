from scipy.optimize import minimize
import matplotlib.pyplot as plt
import numpy as np
import nest
import pandas as pd

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
    "AMPA_g_init": 2.5,
    "AMPA_Tau_r": 1.,
    "AMPA_Tau_d1": 1.,
    "AMPA_Tau_d2": 5.,
    "AMPA_A_r": 2.5,
    "AMPA_A1": 2.5,
    "AMPA_A2": 2.5,
}

def mse(y_true, y_pred):
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    return np.mean((y_true - y_pred) ** 2)

def single_sim(params_grc):
    syn_spec = {"weight": 1., "delay": 1, "receptor_type": 1}
    nest.ResetKernel()
    nest.Install('cerebmodule')
    nest.resolution = 0.025

    grc = nest.Create('eglif_multisyn', 1, params=params_grc)
    input_spikes = nest.Create("spike_generator", params={"spike_times": [(250.-31*0.025)]})
    nest.Connect(input_spikes, grc, syn_spec=dict(syn_spec, synapse_model="static_synapse"))
    multimeter = nest.Create("multimeter", params={"interval": 0.1,
                                                   "record_from": ["V_m","I_syn_ampa"], "record_to": "memory",
                                                   "label": "grc_multimeter",
                                                   "interval": 0.025})
    spike_recorder = nest.Create("spike_recorder", params={"record_to": "memory", "label": "grc_spike_recorder"})
    nest.Connect(multimeter, grc)
    nest.Connect(grc, spike_recorder)
    nest.Simulate(1001)

    multimeter_grc = multimeter.get()['events']
    times_grc = multimeter_grc['times']
    I_syn_AMPA = multimeter_grc['I_syn_ampa']
    V_m_grc = multimeter_grc['V_m']
    AMPA_E_rev = 0.0
    g_syn_AMPA = I_syn_AMPA / (AMPA_E_rev-V_m_grc)

    return times_grc, g_syn_AMPA

def extract_trace():
    file_path = "../../NEURON_traces/Masoli_synapses/AMPA_GrC.csv"
    df = pd.read_csv(file_path, delim_whitespace=True)
    df.columns = ["Time", "Current"]
    df["g_AMPA"] = df["Current"].values / (-40.) * 1000

    return df["Time"], df["g_AMPA"]

def objective_function(params, peak_penalty = True, w=3e-5):
    AMPA_g_init, AMPA_Tau_r, AMPA_Tau_d1, AMPA_Tau_d2, AMPA_A_r, AMPA_A1, AMPA_A2 = params

    params_grc["AMPA_g_init"] = AMPA_g_init
    params_grc["AMPA_Tau_r"] = AMPA_Tau_r
    params_grc["AMPA_Tau_d1"] = AMPA_Tau_d1
    params_grc["AMPA_Tau_d2"] = AMPA_Tau_d2
    params_grc["AMPA_A_r"] = AMPA_A_r
    params_grc["AMPA_A1"] = AMPA_A1
    params_grc["AMPA_A2"] = AMPA_A2

    times_grc, g_syn_AMPA = single_sim(params_grc)
    trace_t, trace_g = extract_trace()
    mask_t = np.logical_and(times_grc <= 1000, times_grc > 250)
    error = mse(trace_g[mask_t], g_syn_AMPA[mask_t])

    if peak_penalty:
        peak_penalty = abs(np.max(trace_g[mask_t]) - np.max(g_syn_AMPA[mask_t])) / np.max(trace_g[mask_t])
        #print('Mse:', error)
        error = error + w * peak_penalty
        #print('Peak penalty:', w*peak_penalty)



    return error


initial_params = [
    2.5,    # AMPA_g_init
    1.,     # AMPA_Tau_r
    1.,     # AMPA_Tau_d1
    5.,     # AMPA_Tau_d2
    2.5 ,   # AMPA_A_r
    2.5,    # AMPA_A1
    2.5     # AMPA_A2
]


bounds = [
    (0, 5),         # AMPA_g_init
    (0, 5),         # AMPA_Tau_r
    (0, 5),         # AMPA_Tau_d1
    (0, 10),        # AMPA_Tau_d2
    (0, 5),         # AMPA_A_r
    (0, 5),         # AMPA_A1
    (0, 5)          # AMPA_A2
]


result = minimize(objective_function, initial_params, bounds=bounds, method='L-BFGS-B', tol = 1e-24)

print("Initial parameters: ", initial_params)
print(f"Optimization completed for {result.message}. Optimal parameters:")
print(f"AMPA_g_init:  {result.x[0]}")
print(f"AMPA_Tau_r: {result.x[1]}")
print(f"AMPA_Tau_d1: {result.x[2]}")
print(f"AMPA_Tau_d2: {result.x[3]}")
print(f"AMPA_A_r: {result.x[4]}")
print(f"AMPA_A1: {result.x[5]}")
print(f"AMPA_A2: {result.x[6]}")
print(f"MAE ottimizzato: {result.fun}")


params_grc["AMPA_g_init"] = result.x[0]
params_grc["AMPA_Tau_r"] = result.x[1]
params_grc["AMPA_Tau_d1"] = result.x[2]
params_grc["AMPA_Tau_d2"] = result.x[3]
params_grc["AMPA_A_r"] = result.x[4]
params_grc["AMPA_A1"] = result.x[5]
params_grc["AMPA_A2"] = result.x[6]

times_grc, g_syn_AMPA = single_sim(params_grc)
trace_t, trace_g = extract_trace()

plt.title('AMPA conductance')
plt.plot(times_grc, g_syn_AMPA, label=f"NEST trace")
plt.plot(times_grc, trace_g, label=f"NEURON trace", linestyle="dashed")
plt.xlabel('Time [ms]')
plt.ylabel(r'$g_{syn_{AMPA}}$ [ns]')
plt.legend()
plt.xlim(240, 280)
#plt.savefig('opt_AMPA.png', dpi=300)
plt.show()