import os

from multi_obj_opt import *

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
    "AMPA_E_rev": 0.0,
}
sim_params = {
    "receptor_type": 1,
    "record_from": ["V_m", "I_syn_ampa"],
    "label_multimeter": "grc_multimeter",
    "label_spike_recorder": "grc_spike_recorder",
    "current_name": "I_syn_ampa",
    "spike_time": 250.0 - 31 * 0.025,
}
trace_file = "../NEURON_traces/Masoli_synapses/AMPA_GrC.csv"
opt_param_names = [
    "AMPA_g_init",
    "AMPA_Tau_r",
    "AMPA_Tau_d1",
    "AMPA_Tau_d2",
    "AMPA_A_r",
    "AMPA_A1",
    "AMPA_A2",
]
problem = SynapseOptimization(
    cell_params,
    sim_params,
    trace_file,
    opt_param_names,
    E_rev=cell_params["AMPA_E_rev"],
    tmin=250,
    tmax=1000,
    duration=1001,
)
problem.xl = np.array([0, 0, 0, 0, 0, 0, 0])
problem.xu = np.array([5, 1.0, 1.0, 10.0, 10, 10, 10])
res = optimization(problem, n_gen=500)
sol_idx, solution = select_solution(res.F, threshold=0.01)
opt_params_dict = dict(zip(opt_param_names, res.X[sol_idx]))
os.makedirs("Optimized_params", exist_ok=True)
output_filename = "Optimized_params/params_AMPA_mf_Grc.txt"
save_params(opt_params_dict, output_filename)
plot_pareto_front(
    res,
    sol_idx,
    threshold=1,
    title="Pareto Front for AMPA mf_Grc",
    save_path="./figs/Pareto_AMPA_mf_GrC.png",
)
sim_params["record_from"] = ["V_m", "I_syn_ampa", "I_ampa1", "I_ampa2"]
plot_conductance(
    params=cell_params,
    sim_params=sim_params,
    trace_file=trace_file,
    E_rev=cell_params["AMPA_E_rev"],
    optimized_params=opt_params_dict,
    title="AMPA conductance mf_Grc",
    xlim=(240, 280),
    save_path="./figs/opt_AMPA_mf_GrC.png",
    duration=1001,
)
