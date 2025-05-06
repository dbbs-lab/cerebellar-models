import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd


# function to estimate the intial time constant for fitting
def estimate_tau(t, g):
    """
    Function to estimate the initial time constant to fit.
    """
    g = np.abs(g)
    g = g[g > 0]  # Avoid log(0)
    t = t[:len(g)]
    log_g = np.log(g)  # convert the Converts the exponential decay into a linear form:
    slope, _ = np.polyfit(t, log_g, 1)  # Fits a straight line to the log-transformed data: log_g = slope * t + intercept.
    return -1 / slope if slope < 0 else 100  # fallback If the slope is positive (which would imply growth, not decay), it returns a fallback value of 100 ms.
    #used linear regression to find the slope and derive τ.

# Single exponential decay function
def single_exponential_fixed_start(t, A1, A2, t_d1, t_d2):
    g0 = df_decay["g_NMDA"].iloc[0]  # actual starting value
    exp_sum = A1 * np.exp(-t / t_d1) + A2 * np.exp(-t / t_d2)
    exp_sum_at_0 = A1 + A2  # since exp(0)=1
    B = g0 - exp_sum_at_0 # force fit to pass through first decay point
    return exp_sum + B

def single_exponential_rise(t, A, t_d, B):
    params['g'] = params['g_max'] * params['gmax_factor'] * params['gbar_Q10']
    return  A* params['g'] * (1-np.exp(-t / t_d)) + B

# Parameters for the system
params = {

    "E_rev": -3.7, # mV
    "V_m": -40,  # mV
    "v0_block": -20 , # Half-max voltage for Mg²⁺ block (mV)
    "k_block": 13 , # Steepness of Mg²⁺ block
    "g_max": 18.8, #  nS
    "gmax_factor": 1,
    "gbar_Q10":  1.4,
}

if __name__ == '__main__':
    file_path = "../NEURON_traces/NMDA_GrC.csv"
    df = pd.read_csv(file_path, delim_whitespace=True)
    df.columns = ["Time", "Current"]
    df["g_NMDA"] = -df["Current"].values / ((params['V_m']-params['E_rev']) * 0.177) * 1000

    max_g_nmda = df["g_NMDA"].min()
    time_max =  df.loc[df["g_NMDA"].idxmin(), "Time"]
    min_g_nmda = df ["g_NMDA"].max()
    time_min = df.loc[df["g_NMDA"].idxmax(), "Time"]

    df_rise = df[df[(df["g_NMDA"] < 0)].index[0]-1:][["Time", "g_NMDA"]]
    init_time_rise = df_rise["Time"].iloc[0]
    df_rise["Time"] = df["Time"] - df_rise["Time"].min()
    time_min_g_nmda = df_rise.loc[df_rise["g_NMDA"].idxmin(), "Time"]
    df_rise = df_rise[df_rise["Time"] <= time_min_g_nmda][["Time", "g_NMDA"]]


    gnmda_max = df_rise["g_NMDA"].min()  # most negative point = max magnitude
    g_scale = params['g_max'] * params['gmax_factor'] * params['gbar_Q10']
    B_initial_rise = df_rise["g_NMDA"].iloc[1]
    #B_initial_rise = df_rise["g_NMDA"].iloc[:2].mean()
    A_lower = (gnmda_max - 0.1) / g_scale

    # Ensure A_upper is negative and conservative
    # A_upper = min(A_upper, 0)
    opt_pars_rise = curve_fit(single_exponential_rise, df_rise["Time"], df_rise["g_NMDA"],method = 'trf',p0=[A_lower, 1, B_initial_rise],bounds=([A_lower,1,B_initial_rise],[0,np.inf,0]))[0]
    gnmda_max = df_rise["g_NMDA"].min()

    print(f"Maximum g_NMDA: {max_g_nmda:.10f} nS")
    print(f'Fitted parameters: \n')
    print(f'A1_Tau_rise: {abs(opt_pars_rise[0])}\nTau_rise: {abs(opt_pars_rise[1])}\nBaseline_rise: {abs(opt_pars_rise[2])}')
    print(f"Initial Baseline Guess for rise: {B_initial_rise}")

    init_time_decay = df.loc[df["g_NMDA"].idxmin(), "Time"]
    df_decay = df[df["Time"] >= init_time_decay][["Time", "g_NMDA"]]
    df_decay["Time"] = df_decay["Time"] - df_decay["Time"].min()  # Normalize time to start from 0


    # Make sure df_decay starts exactly at the max point (no extra points before or after)
    max_idx = df["g_NMDA"].idxmin()  # Index of max (minimum value, most negative)

    df_decay = df.iloc[max_idx:][["Time", "g_NMDA"]]  # Start exactly from the peak
    df_decay["Time"] = df_decay["Time"] - df_decay["Time"].iloc[0]  # Normalize to 0
    B_initial_decay = df_decay["g_NMDA"].iloc[1]
    # Decay fitting with fixed start constraint
    g_values = df_decay["g_NMDA"].values
    t_values = df_decay["Time"].values
    # Get two amplitude estimates from early and late segments
    early_window = g_values[:int(len(g_values)*0.2)]  # First 20% — fast decay region.
    late_window = g_values[int(len(g_values)*0.7):]   # last 30% — slow decay region.
    A1_init = early_window[0] - early_window[-1] # intial amplitude of early window
    A2_init = late_window[0] - late_window[-1] # inital amplitude of slow window
    tau1_init = estimate_tau(t_values[:int(len(t_values)*0.2)], early_window)
    tau2_init = estimate_tau(t_values[int(len(t_values)*0.7):], late_window)
    # Initial parameter guess and bounds
    p0 = [A1_init, A2_init, max(1, tau1_init), max(10, tau2_init)]
    lower_bounds = [-np.inf, -np.inf, 1, 10]
    upper_bounds = [0, 0, 500, 2000]
    # Curve fitting with constrained start
    opt_pars_decay, _ = curve_fit(
        single_exponential_fixed_start,
        df_decay["Time"],
        df_decay["g_NMDA"],
        method='trf',
        p0=p0,
        bounds=(lower_bounds, upper_bounds)
    )

    # # Show decay parameters
    print(f'Fitted decay parameters (constrained start):\n')
    print(f'A1_decay1: {abs(opt_pars_decay[0])}')
    print(f'A2_decay2: {abs(opt_pars_decay[1])}')
    print(f'Tau_decay1: {abs(opt_pars_decay[2])}')
    print(f'Tau_decay2: {abs(opt_pars_decay[3])}')

    max_idx = df["g_NMDA"].idxmin()
    max_time = df["Time"].iloc[max_idx]
    max_value = df["g_NMDA"].iloc[max_idx]


    plt.figure()
    plt.plot(df_rise["Time"] + init_time_rise, np.abs(single_exponential_rise(df_rise["Time"], *opt_pars_rise)), label='Fit rise', color='r')
    plt.plot(df["Time"], np.abs(df["g_NMDA"]), label='NEURON trace', color='gray', linestyle='dashed', alpha=0.5)
    plt.plot(df_decay["Time"] + init_time_decay, np.abs(single_exponential_fixed_start(df_decay["Time"], *opt_pars_decay)), label='Fit decay (constrained)', color='blue')
    #plt.plot(max_time, max_value, marker='*', color='black', markersize=12, label='Max NEURON trace')
    plt.xlim(0, 1000)
    plt.xlabel("Time [ms]")
    plt.ylabel(r"$g_{syn_{NMDA}}$ [nS]")
    plt.legend()
    plt.title("Fit on NMDA conductance")
    plt.savefig('./figs/Fit1_nmda.png', dpi=300)
    plt.show()






