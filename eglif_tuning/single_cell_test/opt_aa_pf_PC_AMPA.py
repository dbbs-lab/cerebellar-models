import os

from multi_obj_opt import *

cell_params = {
    "t_ref": 0.5,
    "V_min": -350,
    "C_m": 334,
    "V_th": -43,
    "V_reset": -69,
    "E_L": -59,
    "V_m": -59.0,
    "lambda_0": 0.5,
    "tau_V": 2.0,
    "tau_m": 47,
    "I_e": 590.0,
    "k_adap": 1.491,
    "k_1": 0.195,
    "k_2": 0.041,
    "A1": 157.622,
    "A2": 172.622,
    "AMPA_E_rev": 0.0,
}

sim_params = {
    "receptor_type": 1,
    "record_from": ["V_m", "I_syn_ampa"],
    "label_multimeter": "pc_multimeter",
    "label_spike_recorder": "pc_spike_recorder",
    "current_name": "I_syn_ampa",
    "spike_time": 200.0 - 36 * 0.025,
}
trace_file = "../NEURON_traces/Danilo_synapses/PCaa_pfAMPA.txt"
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
    tmin=200,
    tmax=1200,
    duration=1201,
)
problem.xl = np.array([0, 0, 0, 0, 0, 0, 0])
problem.xu = np.array([5, 1.0, 5.0, 5.0, 5, 15, 5])
res = optimization(problem, n_gen=500)
sol_idx, solution = select_solution(res.F, threshold=0.01)
opt_params_dict = dict(zip(opt_param_names, res.X[sol_idx]))
os.makedirs("./Optimized_params", exist_ok=True)
output_filename = "./Optimized_params/params_AMPA_aa_pf_PC.txt"
save_params(opt_params_dict, output_filename)
plot_pareto_front(
    res,
    sol_idx,
    threshold=1,
    title="Pareto Front for AMPA aa_pf_PC",
    save_path="./figs/Pareto_AMPA_aa_pf_PC.png",
)
sim_params["record_from"] = ["V_m", "I_syn_ampa", "I_ampa1", "I_ampa2"]
plot_conductance(
    params=cell_params,
    sim_params=sim_params,
    trace_file=trace_file,
    E_rev=cell_params["AMPA_E_rev"],
    optimized_params=opt_params_dict,
    title="AMPA conductance aa_pf_PC",
    xlim=(190, 215),
    save_path="./figs/opt_AMPA_aa_pf_PC.png",
    duration=1201,
)
