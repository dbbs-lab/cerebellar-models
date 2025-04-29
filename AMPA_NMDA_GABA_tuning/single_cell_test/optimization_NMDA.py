from scipy.optimize import minimize
import matplotlib.pyplot as plt
import numpy as np
import nest
import pandas as pd


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
    "NMDA_g_peak": 0.0983349791,
    "NMDA_Tau_rise": 5.1937904737528555,
    "NMDA_Tau_decay1": 0.1538075620377694,
    "NMDA_Tau_decay2": 27.190522726451153,
    "NMDA_A_rise":  0.004458015923894151*NMDA_g_scaling,
    "NMDA_A1_decay1": 9.999999999752668e-07*NMDA_g_scaling,
    "NMDA_A2_decay2": 0.0037741253463865533*NMDA_g_scaling,
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


    syn_spec = {"weight": 1., "delay": 1, "receptor_type": 2}
    nest.ResetKernel()
    nest.Install('cerebmodule')
    nest.resolution = 0.025

    grc = nest.Create('eglif_cond_alpha_multisyn3', 1, params=params_grc)
    input_spikes = nest.Create("spike_generator", params={"spike_times": [(250.-31*0.025)]})
    nest.Connect(input_spikes, grc, syn_spec=dict(syn_spec, synapse_model="static_synapse"))
    multimeter = nest.Create("multimeter", params={"interval": 0.1,
                                                   "record_from": ["V_m", "I_syn_nmda"], "record_to": "memory",
                                                   "label": "grc_multimeter",
                                                    "interval": 0.025})
    spike_recorder = nest.Create("spike_recorder", params={"record_to": "memory", "label": "grc_spike_recorder"})
    nest.Connect(multimeter, grc)
    nest.Connect(grc, spike_recorder)
    nest.Simulate(1001)

    multimeter_grc = multimeter.get()['events']
    times_grc = multimeter_grc['times']
    I_syn_NMDA = multimeter_grc['I_syn_nmda']
    V_m_grc = multimeter_grc['V_m']
    NMDA_E_rev = -3.7
    g_syn_NMDA = I_syn_NMDA / ((NMDA_E_rev - V_m_grc) * 0.04)

    return times_grc, g_syn_NMDA

def extract_trace():
    file_path = "../NEURON_traces/NMDA_GrC.csv"
    df = pd.read_csv(file_path, delim_whitespace=True)
    df.columns = ["Time", "Current"]
    df["g_NMDA"] = df["Current"].values / ((-40.) * 1000 * 0.47)

    return df["Time"], df["g_NMDA"]

def objective_function(params):
    # Estrai i parametri ottimizzati
    NMDA_g_peak, NMDA_Tau_rise, NMDA_Tau_decay1, NMDA_Tau_decay2, NMDA_A_rise, NMDA_A1_decay1, NMDA_A2_decay2 = params

    # Aggiorna i parametri del modello
    params_grc["NMDA_g_peak"] = NMDA_g_peak
    params_grc["NMDA_Tau_rise"] = NMDA_Tau_rise
    params_grc["NMDA_Tau_decay1"] = NMDA_Tau_decay1
    params_grc["NMDA_Tau_decay2"] = NMDA_Tau_decay2
    params_grc["NMDA_A_rise"] =NMDA_A_rise
    params_grc["NMDA_A1_decay1"] = NMDA_A1_decay1
    params_grc["NMDA_A2_decay2"] = NMDA_A2_decay2

    times_grc, g_syn_NMDA= single_sim(params_grc)
    trace_t, trace_g = extract_trace()
    mask_t = np.logical_and(times_grc <= 300, times_grc > 250)

    return mse(trace_g[mask_t], g_syn_NMDA[mask_t])


# Definisci i valori iniziali dei parametri
initial_params = [
    0.0510351113,
    5.1937904737528555,
    0.1538075620377694,
    27.190522726451153,
    0.004458015923894151 * NMDA_g_scaling,
    9.999999999752668e-07 * NMDA_g_scaling,
    0.0037741253463865533 * NMDA_g_scaling,
]

# Definisci i limiti (se necessario)
k = 5000.
bounds = [
    (0, initial_params[0]*k),  # AMPA_g_peak
    (0, initial_params[1]*k),  # AMPA_Tau_r
    (0, initial_params[2]*k),  # AMPA_Tau_d1
    (0, initial_params[3]*k),  # AMPA_Tau_d2
    (0, initial_params[4]*k**2),  # AMPA_A_r
    (0, initial_params[5]*k**2),  # AMPA_A_d1
    (0, initial_params[6]*k**2)  # AMPA_A_d2
]

# Ottimizzazione usando il metodo 'Nelder-Mead' o 'L-BFGS-B'
result = minimize(objective_function, initial_params, bounds=bounds, method='L-BFGS-B', tol = 1e-24)

print(result.message)
# Mostra il risultato ottimizzato
print("Ottimizzazione completata. Parametri ottimali:")
print(f"NMDA_g_peak: {result.x[0]}")
print(f"NMDA_Tau_r: {result.x[1]}")
print(f"NMDA_Tau_d1: {result.x[2]}")
print(f"NMDA_Tau_d2: {result.x[3]}")
print(f"NMDA_A_r: {result.x[4]}")
print(f"NMDA_A_d1: {result.x[5]}")
print(f"NMDA_A_d2: {result.x[6]}")
print(f"MAE ottimizzato: {result.fun}")


params_grc["NMDA_g_peak"] = result.x[0]
params_grc["NMDA_Tau_rise"] = result.x[1]
params_grc["NMDA_Tau_decay1"] = result.x[2]
params_grc["NMDA_Tau_decay2"] = result.x[3]
params_grc["NMDA_A_rise"] = result.x[4]
params_grc["NMDA_A1_decay1"] = result.x[5]
params_grc["NMDA_A2_decay2"] = result.x[6]

print("Parametri iniziali: ", initial_params)
times_grc, g_syn_NMDA = single_sim(params_grc)
trace_t, trace_g = extract_trace()

plt.plot(times_grc, g_syn_NMDA, label=f"NMDA")
plt.plot(times_grc, trace_g, label=f"Recordings", linestyle="dashed")
plt.legend()
plt.xlim(240, 600)
plt.show()