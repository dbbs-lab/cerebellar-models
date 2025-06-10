import os

from multi_obj_opt import *

cell_params = {
    "t_ref": 2,
    "C_m": 145,
    "V_th": -55,
    "V_reset": -75,
    "E_L": -62,
    "V_min": -150,
    "V_m": -60.0,
    "lambda_0": 0.15,
    "tau_V": 0.09,
    "tau_m": 44,
    "I_e": 16.214,
    "k_adap": 0.217,
    "k_1": 0.031,
    "k_2": 0.023,
    "A1": 259.988,
    "A2": 178.01,
    "GABA_E_rev": -70.0,
}
sim_params = {
    "receptor_type": 3,
    "record_from": ["V_m", "I_syn_gaba"],
    "label_multimeter": "GoC_multimeter",
    "label_spike_recorder": "GoC_spike_recorder",
    "current_name": "I_syn_gaba",
    "spike_time": 200.0 - 36 * 0.025,
}
trace_file = "../NEURON_traces/Danilo_synapses/GoCGABA.txt"
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
problem.xu = np.array([5, 5.0, 10.0, 100.0, 15, 10, 10])
res = optimization(problem, n_gen=600)
sol_idx, solution = select_solution(res.F, threshold=0.01)
opt_params_dict = dict(zip(opt_param_names, res.X[sol_idx]))
os.makedirs("./Optimized_params", exist_ok=True)
output_filename = "./Optimized_params/params_GABA_GoC.txt"
save_params(opt_params_dict, output_filename)
plot_pareto_front(
    res,
    sol_idx,
    threshold=1,
    title="Pareto Front for GABA GoC",
    save_path="./figs/Pareto_GABA_GoC.png",
)
sim_params["record_from"] = ["V_m", "I_syn_gaba", "I_gaba1", "I_gaba2"]
plot_conductance(
    params=cell_params,
    sim_params=sim_params,
    trace_file=trace_file,
    E_rev=cell_params["GABA_E_rev"],
    optimized_params=opt_params_dict,
    title="GABA conductance GoC",
    xlim=(190, 400),
    duration=1201,
    save_path="./figs/opt_GABA_GoC.png",
)
