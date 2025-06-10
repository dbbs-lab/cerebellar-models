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
    "NMDA_E_rev": -3.7
}
sim_params = {
    "receptor_type": 2,
    "record_from": ["V_m", "I_syn_nmda", 'Mg_block'],
    "label_multimeter": "grc_multimeter",
    "label_spike_recorder": "grc_spike_recorder",
    "current_name": "I_syn_nmda",
    "spike_time": 250.-31*0.025
}
trace_file = "../NEURON_traces/Masoli_synapses/NMDA_GrC.csv"
opt_param_names = ["NMDA_g_init", "NMDA_Tau_r", "NMDA_Tau_d1", "NMDA_Tau_d2", "NMDA_A_r", "NMDA_A1", "NMDA_A2"]
problem = SynapseOptimization(cell_params, sim_params, trace_file, opt_param_names, E_rev=cell_params['NMDA_E_rev'], mg_block=0.177, tmin=250, tmax=1000)   # For -40 mV patch clamp
problem.xl = np.array([0, 3.0,  10, 100, 0, 0, 0])
problem.xu = np.array([5, 25.0, 60, 200, 5, 5, 5])
res = optimization(problem, n_gen=500)
sol_idx, solution = select_solution(res.F, threshold=0.01)
opt_params_dict = dict(zip(opt_param_names, res.X[sol_idx]))
os.makedirs('Optimized_params', exist_ok=True)
output_filename = "Optimized_params/params_NMDA_mf_GrC.txt"
save_params(opt_params_dict, output_filename)
plot_pareto_front(res, sol_idx, threshold=1, title='Pareto Front for NMDA mf_Grc', save_path='./figs/Pareto_NMDA_mf_GrC.png')
sim_params['record_from']=["V_m", "I_syn_nmda", 'I_nmda1', 'I_nmda2', 'Mg_block']
plot_conductance(params=cell_params, sim_params=sim_params, trace_file=trace_file, E_rev=cell_params['NMDA_E_rev'],
                 optimized_params=opt_params_dict, title="NMDA conductance mf_Grc", xlim = (240,600), mg_block=0.177,
                 save_path="./figs/opt_NMDA_mf_GrC.png", duration=1000)