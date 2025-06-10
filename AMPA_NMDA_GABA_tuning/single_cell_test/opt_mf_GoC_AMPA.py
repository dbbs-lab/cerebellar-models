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
    "AMPA_E_rev": 0.
}

sim_params = {
    "receptor_type": 1,
    "record_from": ["V_m", "I_syn_ampa"],
    "label_multimeter": "goc_multimeter",
    "label_spike_recorder": "goc_spike_recorder",
    "current_name": "I_syn_ampa",
    "spike_time": 200.-36*0.025
}
trace_file = "../NEURON_traces/Danilo_synapses/GoCmfAMPA.txt"
opt_param_names = ["AMPA_g_init", "AMPA_Tau_r", "AMPA_Tau_d1", "AMPA_Tau_d2", "AMPA_A_r", "AMPA_A1", "AMPA_A2"]
problem = SynapseOptimization(cell_params, sim_params, trace_file, opt_param_names, E_rev=cell_params['AMPA_E_rev'],
                              tmin=200, tmax=1200, duration=1201)
problem.xl = np.array([0, 0, 0, 0, 0, 0, 0])
problem.xu = np.array([5, 1., 5., 10., 5, 10, 5])
res = optimization(problem, n_gen=600)
sol_idx, solution = select_solution(res.F, threshold=0.01)
opt_params_dict = dict(zip(opt_param_names, res.X[sol_idx]))
os.makedirs('./Optimized_params', exist_ok=True)
output_filename = "./Optimized_params/params_AMPA_mf_GoC.txt"
save_params(opt_params_dict, output_filename)
plot_pareto_front(res, sol_idx, threshold=1, title='Pareto Front for AMPA mf_GoC', save_path='./figs/Pareto_AMPA_mf_GoC.png')
sim_params['record_from']=["V_m", "I_syn_ampa", 'I_ampa1', 'I_ampa2']
plot_conductance(params=cell_params, sim_params=sim_params, trace_file=trace_file, E_rev=cell_params['AMPA_E_rev'],
                 optimized_params=opt_params_dict, title="AMPA conductance mf_Goc", xlim = (195,250),
                 save_path="./figs/opt_AMPA_mf_GoC.png", duration=1201)