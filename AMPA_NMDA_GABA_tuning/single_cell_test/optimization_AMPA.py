from scipy.optimize import minimize
import matplotlib.pyplot as plt
import numpy as np
import nest
import pandas as pd



AMPA_g_scaling = 1.2 * 1 * 1.4
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
    "AMPA_g_peak": 0.4,
    "AMPA_Tau_r": 0.35863884770941956,
    "AMPA_Tau_d1": 0.3397690354352023,
    "AMPA_Tau_d2": 3.743538818882957,
    "AMPA_A_r": 0.4220193749982735*AMPA_g_scaling,
    "AMPA_A_d1": 0.06369577418128455*AMPA_g_scaling,
    "AMPA_A_d2": 0.17478084529487067*AMPA_g_scaling,
}

def mae(y_true, y_pred):
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    return np.mean(np.abs(y_true - y_pred))

def mse(y_true, y_pred):
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    return np.mean((y_true - y_pred) ** 2)

def single_sim(params_grc):
    # dynamic change of params


    syn_spec = {"weight": 1., "delay": 1, "receptor_type": 1}
    nest.ResetKernel()
    nest.Install('cerebmodule')
    nest.resolution = 0.025

    grc = nest.Create('eglif_cond_alpha_multisyn3', 1, params=params_grc)
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
    g_syn_AMPA = I_syn_AMPA / (AMPA_E_rev - V_m_grc)

    return times_grc, g_syn_AMPA

def extract_trace():
    file_path = "../NEURON_traces/AMPA_GrC.csv"
    df = pd.read_csv(file_path, delim_whitespace=True)
    df.columns = ["Time", "Current"]
    df["g_AMPA"] = df["Current"].values / (-40.) * 1000

    return df["Time"], df["g_AMPA"]

def objective_function(params):
    # Estrai i parametri ottimizzati
    AMPA_g_peak, AMPA_Tau_r, AMPA_Tau_d1, AMPA_Tau_d2, AMPA_A_r, AMPA_A_d1, AMPA_A_d2 = params

    # Aggiorna i parametri del modello
    params_grc["AMPA_g_peak"] = AMPA_g_peak
    params_grc["AMPA_Tau_r"] = AMPA_Tau_r
    params_grc["AMPA_Tau_d1"] = AMPA_Tau_d1
    params_grc["AMPA_Tau_d2"] = AMPA_Tau_d2
    params_grc["AMPA_A_r"] = AMPA_A_r
    params_grc["AMPA_A_d1"] = AMPA_A_d1
    params_grc["AMPA_A_d2"] = AMPA_A_d2

    times_grc, g_syn_AMPA = single_sim(params_grc)
    trace_t, trace_g = extract_trace()
    mask_t = np.logical_and(times_grc <= 500, times_grc > 250)

    return mae(trace_g[mask_t], g_syn_AMPA[mask_t])


# Definisci i valori iniziali dei parametri
initial_params = [
    0.4,  # AMPA_g_peak
    0.35863884770941956,  # AMPA_Tau_r
    0.3397690354352023,  # AMPA_Tau_d1
    3.743538818882957,  # AMPA_Tau_d2
    0.4220193749982735 * AMPA_g_scaling,  # AMPA_A_r
    0.06369577418128455 * AMPA_g_scaling,  # AMPA_A_d1
    0.17478084529487067 * AMPA_g_scaling  # AMPA_A_d2
]

# Definisci i limiti (se necessario)
k = 100
bounds = [
    (0, initial_params[0]*k),  # AMPA_g_peak
    (0, initial_params[1]*k),  # AMPA_Tau_r
    (0, initial_params[2]*k),  # AMPA_Tau_d1
    (0, initial_params[3]*k),  # AMPA_Tau_d2
    (0, initial_params[4]*k),  # AMPA_A_r
    (0, initial_params[5]*k),  # AMPA_A_d1
    (0, initial_params[6]*k)  # AMPA_A_d2
]

# Ottimizzazione usando il metodo 'Nelder-Mead' o 'L-BFGS-B'
result = minimize(objective_function, initial_params, bounds=bounds, method='L-BFGS-B', tol = 1e-24)

print(result.message)
# Mostra il risultato ottimizzato
print("Ottimizzazione completata. Parametri ottimali:")
print(f"AMPA_g_peak: {result.x[0]}")
print(f"AMPA_Tau_r: {result.x[1]}")
print(f"AMPA_Tau_d1: {result.x[2]}")
print(f"AMPA_Tau_d2: {result.x[3]}")
print(f"AMPA_A_r: {result.x[4]}")
print(f"AMPA_A_d1: {result.x[5]}")
print(f"AMPA_A_d2: {result.x[6]}")
print(f"MAE ottimizzato: {result.fun}")


params_grc["AMPA_g_peak"] = result.x[0]
params_grc["AMPA_Tau_r"] = result.x[1]
params_grc["AMPA_Tau_d1"] = result.x[2]
params_grc["AMPA_Tau_d2"] = result.x[3]
params_grc["AMPA_A_r"] = result.x[4]
params_grc["AMPA_A_d1"] = result.x[5]
params_grc["AMPA_A_d2"] = result.x[6]

print("Parametri iniziali: ", initial_params)
times_grc, g_syn_AMPA = single_sim(params_grc)
trace_t, trace_g = extract_trace()

plt.title('AMPA conductance')
plt.plot(times_grc, g_syn_AMPA, label=f"AMPA")
plt.plot(times_grc, trace_g, label=f"Recordings", linestyle="dashed")
plt.xlabel('Time [ms]')
plt.ylabel(r'$g_{syn_{AMPA}}$ [ns]')
plt.legend()
plt.xlim(240, 280)
plt.savefig('opt_AMPA.png')
plt.show()