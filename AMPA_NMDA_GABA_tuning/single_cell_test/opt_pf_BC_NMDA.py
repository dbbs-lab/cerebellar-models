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
    "tau_V": 1.0,
    "tau_m": 9.125,
    "I_e": 3.711,
    "k_adap": 2.025,
    "k_1": 1.887,
    "k_2": 1.096,
    "A1": 5.953,
    "A2": 5.863,
    "NMDA_E_rev": -3.7
}
sim_params = {
    "receptor_type": 2,
    "record_from": ["V_m", "I_syn_nmda", 'Mg_block'],
    "label_multimeter": "bc_multimeter",
    "label_spike_recorder": "bc_spike_recorder",
    "current_name": "I_syn_nmda",
    "spike_time": 200.-36*0.025
}
trace_file = "../NEURON_traces/Danilo_synapses/BCpfNMDA.txt"
opt_param_names = ["NMDA_g_init", "NMDA_Tau_r", "NMDA_Tau_d1", "NMDA_Tau_d2", "NMDA_A_r", "NMDA_A1", "NMDA_A2"]
problem = SynapseOptimization(cell_params, sim_params, trace_file, opt_param_names, E_rev=cell_params['NMDA_E_rev'],
                              mg_block=0.177, tmin=200, tmax=1200, duration=1201)   # For -40 mV patch clamp
problem.xl = np.array([0., 17.,  0, 100, 0, 0, 0])
problem.xu = np.array([5, 50.0, 100, 1000, 10, 10, 10])
res = optimization(problem, n_gen=700)
sol_idx, solution = select_solution(res.F, threshold=0.01)
opt_params_dict = dict(zip(opt_param_names, res.X[sol_idx]))
os.makedirs('./Optimized_params', exist_ok=True)
output_filename = "./Optimized_params/params_NMDA_pf_BC.txt"
save_params(opt_params_dict, output_filename)
plot_pareto_front(res, sol_idx, threshold=1, title='Pareto Front for NMDA pf_BC', save_path='./figs/Pareto_NMDA_pf_BC.png')
sim_params['record_from']=["V_m", "I_syn_nmda", 'I_nmda1', 'I_nmda2', 'Mg_block']
plot_conductance(params=cell_params, sim_params=sim_params, trace_file=trace_file, E_rev=cell_params['NMDA_E_rev'],
                 optimized_params=opt_params_dict, title="NMDA conductance pf_BC", xlim = (190,700), mg_block=0.177, duration=1201,
                 save_path="./figs/opt_NMDA_pf_BC.png")