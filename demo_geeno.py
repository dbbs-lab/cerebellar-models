import numpy as np
import nest
import matplotlib.pyplot as plt

#from test_general import offset


# from test_general import syn_spec_dict


def run_simulation():
    nest.ResetKernel()
    nest.Install("custom_stdp_module")
    #nest.Install("cerebmodule")

    params_pc = {"t_ref": 0.5,
          "Vmin": -350,
          "C_m": 334,
          "V_th": -43,
          "V_reset": -69,
          "E_L": -59,
          "V_m": -59.0,
          "lambda_0": 0.001,
          "tau_V": 0.5,
          "tau_m": 47,
          "I_e": 742.543,
          "kadap": 1.491,
          "k1": 0.195,
          "k2": 0.041,
          "A1": 157.622,
          "A2": 172.622,
          "tau_syn1": 1.1,
          "tau_syn2": 2.8,
          "tau_syn3": 0.4,
          "E_rev1": 0,
          "E_rev2": -80,
          "E_rev3": 0}

    params_io= {"t_ref": 1,
          "C_m": 189,
          "V_th": -35,
          "V_reset": -45,
          "E_L": -45,
          "I_e": 18.101,
          "V_m": -45,
          "lambda_0": 1.2,
          "tau_V": 0.8,
          "tau_m": 11,
          "kadap": 1.928,
          "k1": 0.191,
          "k2": 0.091,
          "A1": 1810.923,
          "A2": 1358.197,
          "tau_syn1": 1,
          "tau_syn2": 60,
          "E_rev1": 0,
          "E_rev2": -80}

    params_gr = {"t_ref": 1.5,
          "Vmin": -150,
          "C_m": 7,
          "V_th": -41,
          "V_reset": -70,
          "E_L": -62,
          "I_e": -0.888,
          "V_m": -62.0,
          "lambda_0": 1.0,
          "tau_V": 0.3,
          "tau_m": 24.15,
          "kadap": 0.022,
          "k1": 0.311,
          "k2": 0.041407868,
          "A1": 0.01,
          "A2": -0.94,
          "tau_syn1": 5.8,
          "tau_syn2": 13.61,
          "E_rev1": 0,
          "E_rev2": -80}


    io = nest.Create("eglif_io_nestml", 1, params = params_io)
    pc = nest.Create("eglif_pc_nestml", 1, params = params_pc)
    gr = nest.Create("eglif_gr_nestml", 1, params=params_gr)
    wr = nest.Create("weight_recorder")
    nest.CopyModel("stdp_synapse_sinexp", "stdp_rec", {"weight_recorder": wr})
    #nest.Connect(pre, post, syn_spec={"synapse_model" : "stdp_synapse_sinexp" , "receptor_type" : 5})
    nest.Connect(io, pc, syn_spec={"synapse_model" : "static_synapse" ,"receptor_type" : 5})
    nest.Connect(gr, pc, syn_spec={"synapse_model": "stdp_rec", "receptor_type" : 1})
    #nest.SetStatus(pre, {"offset": 1})
    #nest.SetStatus(post, {"offset": 1})

    io_generator = nest.Create("spike_generator", params = {"spike_times": [10]})
    gr_generator = nest.Create("spike_generator", params={"spike_times": [10, 100, 200]})
    #gr_2_gen = nest.Create("poisson_generator", params={"rate": 500, "start": 800, "stop": 810})

    sr_io, sr_pc, sr_gr = nest.Create("spike_recorder", 3)
    nest.Connect(io, sr_io)
    nest.Connect(pc, sr_pc)
    nest.Connect(gr, sr_gr)
    nest.Connect(io_generator, io, syn_spec={"receptor_type" : 1, "weight" : 5})
    nest.Connect(gr_generator, gr, syn_spec={"receptor_type" : 1,"weight" : 1})
    #nest.Connect(gr_2_gen, gr, syn_spec={"receptor_type": 5})
    #nest.Connect(post, pre, syn_spec={"synapse_model": "stdp_synapse_sinexp", "receptor_type" : 5})

    #pre.I_e = 1E4    # [pA]
    #post.I_e = .8E3    # [pA]

    pf_pc_conns = nest.GetConnections(synapse_model="stdp_rec")

    weights = np.zeros((1000, 1))
    for i in range(1000):
        nest.Simulate(1)
        weights[i] = nest.GetStatus(pf_pc_conns, "weight")
        # if weights[i] > 0.14:
        #     print("timestamp complex spike: ", i)
        # if i > 500:
        #     print("roba succede: ", weights[i])
        # if weights[i] < 0.14 and weights[i] > 0.12:
        #     print("SIUUUUUUUUUUM: ", i)
        #     print("PAZZESCOOOOOO: ", weights[i])


    return sr_io.events["times"], sr_gr.events["times"], sr_pc.events["times"], weights


io, gr, pc, w = run_simulation()
print(np.diff(w, 0))
print("IO SPIKES: ", io)
print("GR SPIKES: ", gr)
print("PC SPIKES: ", pc)
i =  np.zeros((999, 1))
counter = 0
for z in range(999):
    i[z] = counter
    counter = counter +1
plt.figure()
plt.plot(i, np.diff(w, axis =0))
plt.xlim(0,600)
plt.show()

plt.figure()
plt.plot(np.append(i, 1000), w)
plt.xlim(0,600)
plt.show()
# counter = 0
# for i in range(1000):
#     if w[i] < 0.14:
#         counter = counter +1
# print(counter)
# pre_t_sp_with_offset_trick, post_t_sp_with_offset_trick, w_with_offset_trick = run_simulation(use_offset_trick=True)
# pre_t_sp_without_offset_trick, post_t_sp_without_offset_trick, w_without_offset_trick = run_simulation(use_offset_trick=False)
#
# np.testing.assert_allclose(pre_t_sp_with_offset_trick, pre_t_sp_without_offset_trick)
# np.testing.assert_allclose(post_t_sp_with_offset_trick, post_t_sp_without_offset_trick)
# np.testing.assert_allclose(w_with_offset_trick, w_without_offset_trick)
