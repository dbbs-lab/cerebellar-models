import numpy as np
import nest
import matplotlib.pyplot as plt

#from test_general import offset


# from test_general import syn_spec_dict


def run_simulation():
    nest.ResetKernel()
    nest.Install("custom_stdp_module")
    nest.Install("cerebmodule")

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
          "I_e": -18.101,
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

    params_dcni = {"t_ref": 3,
          "C_m": 56,
          "V_th": -39,
          "V_reset": -55,
          "E_L": -40,
          "I_e": 10., # 2.384
          "V_m": -40,
          "lambda_0": 1., # 0.001
          "tau_V": 0.3,  # 0.5
          "tau_m": 56,
          "k_adap": 0.079,
          "k_1": 0.041,
          "k_2": 0.044,
          "A1": 176.358,
          "A2": 176.358,
          "tau_syn1": 3.64,
          "tau_syn2": 1.14,
          "E_rev1": 0,
          "E_rev2": -80}

    io = nest.Create("eglif_io_nestml", 1, params = params_io)
    parrot = nest.Create("parrot_neuron", 1)
    pc = nest.Create("eglif_pc_nestml", 1, params = params_pc)
    gr = nest.Create("eglif_gr_nestml", 1, params=params_gr)
    dcn_i = nest.Create("eglif_cond_alpha_multisyn", 1, params=params_dcni)
    wr = nest.Create("weight_recorder")
    nest.CopyModel("stdp_synapse_sinexp", "stdp_rec", {"weight_recorder": wr})
    #nest.Connect(pre, post, syn_spec={"synapse_model" : "stdp_synapse_sinexp" , "receptor_type" : 5})
    # nest.Connect(io, pc, syn_spec={"synapse_model" : "static_synapse" ,"receptor_type" : 5})
    nest.Connect(parrot, pc, syn_spec={"synapse_model" : "static_synapse" ,"receptor_type" : 5})
    nest.Connect(gr, pc, syn_spec={"synapse_model": "stdp_rec", "receptor_type" : 1})
    nest.Connect(dcn_i, io, syn_spec={"synapse_model" : "static_synapse" ,"weight": 1.25, "receptor_type" : 2})
    #nest.SetStatus(pre, {"offset": 1})
    #nest.SetStatus(post, {"offset": 1})

    io_generator = nest.Create("spike_generator", params = {"spike_times": [298, 315]})
    gr_generator = nest.Create("spike_generator", params={"spike_times": [ 225, 230, 235, 400 ]})
    # io_pois_1 = nest.Create("poisson_generator", params={"rate": 500, "start": 300, "stop": 310})
    # io_pois_2 = nest.Create("poisson_generator", params={"rate": 500, "start": 700, "stop": 710})

    sr_io, sr_pc, sr_gr = nest.Create("spike_recorder", 3)
    nest.Connect(parrot, sr_io)
    nest.Connect(pc, sr_pc)
    nest.Connect(gr, sr_gr)
    # nest.Connect(io_pois_1, io, syn_spec={"receptor_type" : 1, "weight" : 20})
    # nest.Connect(io_pois_2, io, syn_spec={"receptor_type": 1, "weight": 20})
    nest.Connect(gr_generator, gr, syn_spec={"receptor_type" : 1,"weight" : 1})
    nest.Connect(io_generator, parrot, syn_spec={"weight" : 1})
    #nest.Connect(gr_2_gen, gr, syn_spec={"receptor_type": 5})
    #nest.Connect(post, pre, syn_spec={"synapse_model": "stdp_synapse_sinexp", "receptor_type" : 5})

    #pre.I_e = 1E4    # [pA]
    #post.I_e = .8E3    # [pA]

    pf_pc_conns = nest.GetConnections(synapse_model="stdp_rec")

    io_pot = nest.Create("voltmeter")
    nest.Connect(io_pot, pc)

    # nest.Simulate(1000)
    # nest.voltage_trace.from_device(io_pot)
    # plt.axvline(x= 50, linestyle='--', color="r")
    # plt.axvline(x= 400, linestyle='--', color="r")
    # plt.axvline(x= 500, linestyle='--', color="r")
    # # plt.axhline(y=params_pc["V_th"], linestyle="--", color="r")
    # plt.show()

    # nest.Simulate(1000)
    # weights = nest.GetStatus(pf_pc_conns, "weight")
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


    return sr_io.events, sr_gr.events, sr_pc.events, weights


io, gr, pc, w = run_simulation()
#print(np.diff(w, 0))
print("IO SPIKES: ", io["times"])
print("GR SPIKES: ", gr["times"])
print("PC SPIKES: ", pc["times"])
plt.scatter(io["times"], io["senders"], c="b", s=7)
plt.scatter(gr["times"], gr["senders"], c="r", s=7)
plt.scatter(pc["times"], pc["senders"], c="g", s=7)
plt.title("Scatter plot for plasticity test")
plt.axvline(x= 100, linestyle='--', color="k")
plt.axvline(x= 300, linestyle='--', color="k")
plt.axvline(x= 600, linestyle='--', color="c")
plt.axvline(x= 800, linestyle='--', color="c")
plt.legend({"IO", "GR", "PC"}, loc="best")
i =  np.zeros((999, 1))
counter = 0
for z in range(999):
    i[z] = counter
    counter = counter +1
plt.figure()
plt.plot(i, np.diff(w, axis =0))
plt.title("Weight diff")
plt.axvline(x= 100, linestyle='--', color="k")
plt.axvline(x= 300, linestyle='--', color="k")
plt.axvline(x= 600, linestyle='--', color="c")
plt.axvline(x= 800, linestyle='--', color="c")
#plt.xlim(0,600)
plt.show()

plt.figure()
plt.plot(np.append(i, 1000), w)
plt.title("Weight evolution")
plt.axvline(x= 100, linestyle='--', color="k")
plt.axvline(x= 300, linestyle='--', color="k")
plt.axvline(x= 600, linestyle='--', color="c")
plt.axvline(x= 800, linestyle='--', color="c")
#plt.xlim(0,600)
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
LTD = np.array([0.00203423,0.00225064,0.00247454, 0.00270388,0.00293637, 0.00316948, 0.00340048, 0.00363473 ,0.00385241,0.00405902,0.00424474,0.00442108,0.00457795,0.00471289,0.00482011,0.0049061,0.00496489,0.00499551,0.00499747,0.00497072,0.00491574,0.00483695,0.00472961,0.00459798,0.00444419,0.00427063,0.00407997,0.00387504,0.00365878,0.0034342,0.00320425,0.00297182,0.00274821,0.00251869,0.00230242,0.00208474,0.00187588,0.00167731,0.00149692,0.0013217,0.00115938,0.00101028,0.000879261,0.000756088,0.000645696,0.000550989,0.000464057,0.00039065,0.000324329,0.000269213,0.000220209,0.000180139,0.000143907,0.000112998,8.40677e-05,5.72521e-05,3.28335e-05,1.0543e-05])
gr_spikes = np.array([102.3, 105,  107.7, 110.4, 113.1, 115.8, 118.5, 121.3, 124,  126.7, 129.3, 132,  134.7, 137.4, 140,  142.7, 145.4, 148.1, 150.8, 153.5, 156.2, 158.8, 161.5, 164.2, 166.9, 169.6, 172.3, 175,  177.7,180.4, 183.1, 185.8, 188.4, 191.1, 193.7, 196.4, 199.1, 201.8, 204.4, 207.1, 209.8, 212.5, 215.1, 217.8, 220.5, 223.1, 225.8, 228.4, 231.1, 233.7, 236.4, 239,  241.8, 244.7, 248.1, 252.3, 258,  268.5])
a=0
for i in range(np.size(LTD)):
        if LTD[i] == max(LTD):
           a = i
           print(gr_spikes[a])
           print(max(LTD))
           break
plt.plot(gr_spikes, LTD)
plt.plot(gr_spikes[a], max(LTD), marker='.', markersize=10, color='green')
plt.annotate('[t_Gr, LTD_value] = [150.8, 0.00499747]', xy=(gr_spikes[a], max(LTD)))
plt.axvline(x= 298, linestyle='--', color="r")
plt.axvline(x= 98, linestyle='--', color="r")
plt.title('LTD curve')
plt.legend({ "Time window", "Max value","LTD"}, loc='best')
plt.show()
