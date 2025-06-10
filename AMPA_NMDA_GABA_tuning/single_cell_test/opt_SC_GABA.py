import os

from multi_obj_opt import *

cell_params = {
    "t_ref": 1.59,
    "C_m": 14.6,
    "V_th": -53,
    "V_reset": -78,
    "E_L": -68,
    "V_m": -60.0,
    "lambda_0": 1.8,
    "tau_V": 1,
    "tau_m": 9.125,
    "I_e": 3.711,
    "k_adap": 2.025,
    "k_1": 1.887,
    "k_2": 1.096,
    "A1": 5.953,
    "A2": 5.863,
    "GABA_E_rev": -70.0,
}
sim_params = {
    "receptor_type": 3,
    "record_from": ["V_m", "I_syn_gaba"],
    "label_multimeter": "SC_multimeter",
    "label_spike_recorder": "SC_spike_recorder",
    "current_name": "I_syn_gaba",
    "spike_time": 200.0 - 36 * 0.025,
}
trace_file = "../NEURON_traces/Danilo_synapses/SCGABA.txt"
opt_param_names = [
    "GABA_g_init",
    "GABA_Tau_r",
    "GABA_Tau_d1",
    "GABA_Tau_d2",
    "GABA_A_r",
    "GABA_A1",
    "GABA_A2",
]
problem = SynapseOptimization(
    cell_params,
    sim_params,
    trace_file,
    opt_param_names,
    E_rev=cell_params["GABA_E_rev"],
    duration=1201,
    tmin=200,
    tmax=1200,
)
problem.xl = np.array([0, 0.0, 0, 10, 0, 0, 0])
problem.xu = np.array([5, 5.0, 10.0, 100.0, 5, 10, 5])
res = optimization(problem, n_gen=600)
sol_idx, solution = select_solution(res.F, threshold=0.01)
opt_params_dict = dict(zip(opt_param_names, res.X[sol_idx]))
os.makedirs("./Optimized_params", exist_ok=True)
output_filename = "./Optimized_params/params_GABA_SC.txt"
save_params(opt_params_dict, output_filename)
plot_pareto_front(
    res,
    sol_idx,
    threshold=1,
    title="Pareto Front for GABA SC",
    save_path="./figs/Pareto_GABA_SC.png",
)
sim_params["record_from"] = ["V_m", "I_syn_gaba", "I_gaba1", "I_gaba2"]
plot_conductance(
    params=cell_params,
    sim_params=sim_params,
    trace_file=trace_file,
    E_rev=cell_params["GABA_E_rev"],
    optimized_params=opt_params_dict,
    title="GABA conductance SC",
    xlim=(190, 600),
    duration=1201,
    save_path="./figs/opt_GABA_SC.png",
)
