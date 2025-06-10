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
    "GABA_E_rev": -70.
}
sim_params = {
    "receptor_type": 3,
    "record_from": ["V_m", "I_syn_gaba"],
    "label_multimeter": "grc_multimeter",
    "label_spike_recorder": "grc_spike_recorder",
    "current_name": "I_syn_gaba",
    "spike_time": 200.-36*0.025
}
trace_file = "../NEURON_traces/Danilo_synapses/GrCGABA.txt"
opt_param_names = ["GABA_g_init", "GABA_Tau_r", "GABA_Tau_d1", "GABA_Tau_d2", "GABA_A_r", "GABA_A1", "GABA_A2"]
problem = SynapseOptimization(cell_params, sim_params, trace_file, opt_param_names, E_rev=cell_params['GABA_E_rev'], duration=1201,
                              tmin=200, tmax=1200)
problem.xl = np.array([0, 0., 0, 50, 0, 0, 0])
problem.xu = np.array([5, 1., 40., 350., 5, 15, 5])
res = optimization(problem, n_gen=600)
sol_idx, solution = select_solution(res.F, threshold=0.01)
opt_params_dict = dict(zip(opt_param_names, res.X[sol_idx]))
os.makedirs('./Optimized_params', exist_ok=True)
output_filename = "./Optimized_params/params_GABA_GoC_GrC.txt"
save_params(opt_params_dict, output_filename)
plot_pareto_front(res, sol_idx, threshold=1, title='Pareto Front for GABA GoC_Grc', save_path='./figs/Pareto_GABA_GoC_GrC.png')
sim_params['record_from']=["V_m", "I_syn_gaba", 'I_gaba1', 'I_gaba2']
plot_conductance(params=cell_params, sim_params=sim_params, trace_file=trace_file, E_rev=cell_params['GABA_E_rev'],
                 optimized_params=opt_params_dict, title="GABA conductance GoC_GrC", xlim = (190,600), duration=1201,
                 save_path="./figs/opt_GABA_GoC_GrC.png")