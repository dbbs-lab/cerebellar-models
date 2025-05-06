import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

params = {
    "V_m": -40,  # Membrane potential (mV)
    "E_rev": 0,  # Reversal potential (mV)
    "g_max": 1.2, #  nS
    "gmax_factor": 1,
    "gbar_Q10": 1.4,
}


def double_exponential_decay(t, A1,A2, t_d1,t_d2, B_d):
    params['g'] = params['g_max'] * params['gmax_factor'] * params['gbar_Q10']
    return params['g'] * (A1*np.exp(-t/t_d1) + A2*np.exp(-t / t_d2)) + B_d

def single_exponential_rise(t, A_r, t_r, B_r):
    params['g'] = params['g_max'] * params['gmax_factor'] * params['gbar_Q10']
    return  A_r*params['g'] * (1-np.exp(-t / t_r)) + B_r

if __name__ == '__main__':
    file_path = "../NEURON_traces/AMPA_GrC.csv"
    df = pd.read_csv(file_path, delim_whitespace=True)
    df.columns = ["Time", "Current"]
    df["g_AMPA"] = -df["Current"].values * 1000 / (params['V_m']-params['E_rev'])   # multiplying 10^3 to convert conductance in nS

    max_g_ampa = df["g_AMPA"].min()
    time_max =  df.loc[df["g_AMPA"].idxmin(), "Time"]
    min_g_ampa = df ["g_AMPA"].max()
    time_min = df.loc[df["g_AMPA"].idxmax(), "Time"]

    df_rise = df[df[(df["g_AMPA"] < 0)].index[0]-1:][["Time", "g_AMPA"]]
    init_time_rise = df_rise["Time"].iloc[0]
    df_rise["Time"] = df["Time"] - df_rise["Time"].min()
    time_min_g_ampa = df_rise.loc[df_rise["g_AMPA"].idxmin(), "Time"]
    df_rise = df_rise[df_rise["Time"] <= time_min_g_ampa][["Time", "g_AMPA"]]
    B_initial_rise = df_rise["g_AMPA"].min()
    opt_pars_rise = curve_fit(single_exponential_rise, df_rise["Time"], df_rise["g_AMPA"],method = 'trf',p0=[0.356161, 0.6785, B_initial_rise],bounds=(-np.inf, np.inf))[0]
    print(f'Fitted parameters: \n')
    print(f'A_r: {abs(opt_pars_rise[0])}\nTau_r: {abs(opt_pars_rise[1])}\nB_r: {abs(opt_pars_rise[2])}\n')

    init_time_decay = df.loc[df["g_AMPA"].idxmin(), "Time"]
    df_decay = df[df["Time"] >= init_time_decay][["Time", "g_AMPA"]]
    df_decay["Time"] = df_decay["Time"] - df_decay["Time"].min()  # Normalize time to start from 0
    B_initial_decay = df_decay["g_AMPA"].median()
    opt_pars_decay = curve_fit(double_exponential_decay, df_decay["Time"], df_decay["g_AMPA"],method = 'trf',p0=[-0.00040, -0.00040, 0.0251, 1.042, B_initial_decay],bounds=(-np.inf, np.inf))[0]

    print(f'Fitted parameters: \n')
    print(f'A1: {abs(opt_pars_decay[0])}\nA2: {abs(opt_pars_decay[1])}\nTau_d1: {abs(opt_pars_decay[2])})\n', f'Tau_d2: {abs(opt_pars_decay[3])})\nB_d:{abs(opt_pars_decay[4])}')
    print(f"Initial Baseline Guess for Decay: {B_initial_decay}")

    max_idx = df["g_AMPA"].idxmin()
    max_time = df["Time"].iloc[max_idx]
    max_value = df["g_AMPA"].iloc[max_idx]

    plt.figure()
    plt.title('Fit on AMPA conductance')
    plt.plot(df_rise["Time"]+ init_time_rise, np.abs(single_exponential_rise(df_rise["Time"], *opt_pars_rise)), label='Fit rise', color='r')
    plt.plot(df["Time"], np.abs(df["g_AMPA"]), label='NEURON trace', color='gray', linestyle='dashed', alpha=0.5)
    plt.plot(df_decay["Time"]+init_time_decay, np.abs(double_exponential_decay(df_decay["Time"], *opt_pars_decay)), label='Fit decay', color='blue')
    #plt.plot(max_time, np.abs(max_value), marker='*', color='black', markersize=12, label='Max NEURON trace')
    plt.xlabel('Time [ms]')
    plt.ylabel(r'$g_{syn_{AMPA}}$ [nS]')
    #plt.xlim(200,350)
    plt.legend()
    plt.savefig('./figs/Fit1_ampa.png', dpi=300)
    plt.show()