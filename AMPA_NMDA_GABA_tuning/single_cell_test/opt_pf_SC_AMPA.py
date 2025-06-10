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
    "AMPA_E_rev": 0.
}
sim_params = {
    "receptor_type": 1,
    "record_from": ["V_m", "I_syn_ampa"],
    "label_multimeter": "SC_multimeter",
    "label_spike_recorder": "SC_spike_recorder",
    "current_name": "I_syn_ampa",
    "spike_time": 200.-36*0.025
}
trace_file = "../NEURON_traces/Danilo_synapses/SCpfAMPA.txt"
opt_param_names = ["AMPA_g_init", "AMPA_Tau_r", "AMPA_Tau_d1", "AMPA_Tau_d2", "AMPA_A_r", "AMPA_A1", "AMPA_A2"]
problem = SynapseOptimization(cell_params, sim_params, trace_file, opt_param_names, E_rev=cell_params['AMPA_E_rev'],
                              tmin=200, tmax=1200, duration=1201)
problem.xl = np.array([0, 0, 0, 0, 0, 0, 0])
problem.xu = np.array([5, 1., 1., 5., 10, 20, 15])
res = optimization(problem, n_gen=500)
sol_idx, solution = select_solution(res.F, threshold=0.01)
opt_params_dict = dict(zip(opt_param_names, res.X[sol_idx]))
os.makedirs('./Optimized_params', exist_ok=True)
output_filename = "./Optimized_params/params_AMPA_pf_SC.txt"
save_params(opt_params_dict, output_filename)
plot_pareto_front(res, sol_idx, threshold=1, title='Pareto Front for AMPA pf_SC', save_path='./figs/Pareto_AMPA_pf_SC.png')
sim_params['record_from']=["V_m", "I_syn_ampa", 'I_ampa1', 'I_ampa2']
plot_conductance(params=cell_params, sim_params=sim_params, trace_file=trace_file, E_rev=cell_params['AMPA_E_rev'],
                 optimized_params=opt_params_dict, title="AMPA conductance pf_SC", xlim = (190,215),
                 save_path="./figs/opt_AMPA_pf_SC.png", duration=1201)