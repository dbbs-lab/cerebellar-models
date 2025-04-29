import matplotlib.pyplot as plt
import numpy as np
import nest
import pandas as pd
from scipy.interpolate import interp1d

AMPA_g_scaling = 1.2 * 1 * 35
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
    "AMPA_Tau_r": 0.36,
    "AMPA_Tau_d1": 0.34,
    "AMPA_Tau_d2": 3.7,
    "AMPA_A_r": 0.017*AMPA_g_scaling,
    "AMPA_A_d1": 0.0025*AMPA_g_scaling,
    "AMPA_A_d2": 0.007*AMPA_g_scaling,
}

def mae(y_true, y_pred):
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    return np.mean(np.abs(y_true - y_pred))

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

if __name__ == "__main__":
    ks = np.arange(0, 2, 0.1)
    actual_maes = []
    for k in ks:
        params_grc["AMPA_g_peak"] =  0.4 * 3 *k,
        params_grc["AMPA_Tau_r"] = 0.36*k,
        params_grc["AMPA_Tau_d1"]=  0.34*k,
        params_grc["AMPA_Tau_d2"]= 3.7*k
        params_grc["AMPA_A_r"]=  0.017 * AMPA_g_scaling*k,
        params_grc["AMPA_A_d1"]= 0.0025 * AMPA_g_scaling*k,
        params_grc["AMPA_A_d2"]=  0.007 * AMPA_g_scaling*k,
        times_grc, g_syn_AMPA = single_sim(params_grc)
        trace_t, trace_g = extract_trace()


        plt.plot(times_grc, g_syn_AMPA, label=f"AMPA {k}")
        plt.plot(times_grc, trace_g, label=f"Recordings  {k}")
        plt.legend()
        plt.xlim(240, 260)
        plt.show()
        actual_mae = mae(trace_g, g_syn_AMPA)
        actual_maes.append([k, actual_mae])
        print(f"Actual MAE for k={k}:", actual_mae)




plt.plot(np.array(actual_maes).T[0],np.array(actual_maes).T[1], '+')
plt.xlabel("k scaling")
plt.ylabel("MAE")
plt.show()