import numpy as np
import nest
import matplotlib.pyplot as plt

#from test_general import offset


# from test_general import syn_spec_dict


def run_simulation():
    nest.ResetKernel()
    nest.Install("custom_stdp_module")
    nest.Install("cerebmodule")

    params_mli = {"t_ref": 1.59,
                 "C_m": 14.6,
                 "V_th": -53,
                 "V_reset": -78,
                 "E_L": -68,
                 "V_m": -60.0,
                 "lambda_0": 1.8,
                 "tau_V": 1.,
                 "tau_m": 9.125,
                 "I_e": 3.711,
                 "kadap": 2.025,
                 "k1": 1.887,
                 "k2": 1.096,
                 "A1": 5.953,
                 "A2": 5.863,
                 "tau_syn1": 0.64,
                 "tau_syn2": 2,
                 "tau_syn3": 1.2,
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
          "k_adap": 1.928,
          "k_1": 0.191,
          "k_2": 0.091,
          "A1": 1810.923,
          "A2": 1358.197,
          "tau_syn1": 1,
          "tau_syn2": 60,
          "E_rev1": 0,
          "E_rev2": -80,
          "E_rev3": 0}

    params_gr = {"t_ref": 1.5,
          "V_min": -150,
          "C_m": 7,
          "V_th": -41,
          "V_reset": -70,
          "E_L": -62,
          "I_e": -0.888,
          "V_m": -62.0,
          "lambda_0": 1.0,
          "tau_V": 0.3,
          "tau_m": 24.15,
          "k_adap": 0.022,
          "k_1": 0.311,
          "k_2": 0.041407868,
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

    io = nest.Create("eglif_cond_alpha_multisyn", 1, params = params_io)
    mli = nest.Create("eglif_mli", 1, params = params_mli)
    gr = nest.Create("eglif_cond_alpha_multisyn", 1, params=params_gr)
    dcn_i = nest.Create("eglif_cond_alpha_multisyn", 1, params=params_dcni)
    wr = nest.Create("weight_recorder")
    nest.CopyModel("stdp_synapse_alpha", "stdp_rec", {"weight_recorder": wr})
    nest.Connect(io, mli, syn_spec={"synapse_model": "static_synapse" , "receptor_type": 3})
    nest.Connect(dcn_i, io, syn_spec={"synapse_model" : "static_synapse" ,"weight": 1.25, "receptor_type" : 2})
    nest.Connect(gr, mli, syn_spec={"synapse_model" : "stdp_rec", "receptor_type" : 1})

    #io_generator = nest.Create("spike_generator", params = {"spike_times": [300, 800]})
    gr_generator = nest.Create("spike_generator", params={"spike_times": [ 50, 60, 70, 80, 120, 500, 600, 730, 900]})
    io_pois_1 = nest.Create("poisson_generator", params={"rate": 500, "start": 298, "stop": 310})
    io_pois_2 = nest.Create("poisson_generator", params={"rate": 500, "start": 800, "stop": 810})

    sr_io, sr_mli, sr_gr = nest.Create("spike_recorder", 3)
    nest.Connect(io, sr_io)
    nest.Connect(mli, sr_mli)
    nest.Connect(gr, sr_gr)
    nest.Connect(io_pois_1, io, syn_spec={"receptor_type" : 1, "weight" : 20})
    nest.Connect(io_pois_2, io, syn_spec={"receptor_type": 1, "weight": 20})
    nest.Connect(gr_generator, gr, syn_spec={"receptor_type" : 1,"weight" : 1})
    #nest.Connect(io_generator, io, syn_spec={"weight" : 1})
    #nest.Connect(gr_2_gen, gr, syn_spec={"receptor_type": 5})


    pf_mli_conns = nest.GetConnections(synapse_model="stdp_rec")

    io_pot = nest.Create("voltmeter")
    nest.Connect(io_pot, mli)

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
        weights[i] = nest.GetStatus(pf_mli_conns, "weight")


    return sr_io.events, sr_gr.events, sr_mli.events, weights


io, gr, mli, w = run_simulation()
print(np.diff(w, 0))
print("IO SPIKES: ", io["times"])
print("GR SPIKES: ", gr["times"])
print("MLI SPIKES: ", mli["times"])
plt.scatter(io["times"], io["senders"], c="b", s=7, label = "IO")
plt.scatter(gr["times"], gr["senders"], c="r", s=7, label = "GR")
plt.scatter(mli["times"], mli["senders"], c="g", s=7, label = "MLI")
plt.title("Scatter plot for plasticity test")
plt.axvline(x= 100, linestyle='--', color="k")
plt.axvline(x= 300, linestyle='--', color="k")
plt.axvline(x= 600, linestyle='--', color="c")
plt.axvline(x= 800, linestyle='--', color="c")
plt.legend(loc="best")
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



# LTP = [0.205817,0.21491,0.224355,0.233817,0.243624,0.253417,0.263546,0.27402,0.284844,0.296026,0.307571,0.319935,0.33224,0.345403,
#         0.357996,
#         0.371455,
#         0.385305,
#         0.399549,
#         0.414187,
#         0.429219,
#         0.444645,
#         0.460461,
#         0.476663,
#         0.492625,
#         0.509566,
#         0.526224,
#         0.543868,
#         0.56185,
#         0.580153,
#         0.598757,
#         0.61764,
#         0.636775,
#         0.656131,
#         0.675676,
#         0.695369,
#         0.715168,
#         0.735023,
#         0.75488,
#         0.774677,
#         0.794347,
#         0.813099,
#         0.832294,
#         0.851117,
#         0.869465,
#         0.88723,
#         0.904289,
#         0.92051,
#         0.935747,
#         0.949841,
#         0.962618,
#         0.973889,
#         0.983448,
#         0.991068,
#         0.996507,0.999499,0.999797,0.997123,0.991074,0.981284,0.967357,0.949635,0.926306,0.897465,0.862572,0.821043,0.772251,0.715517,0.64751,0.572276,0.483319,0.366719,0.205351]
#
# gr_sp = [105.7, 109.4, 112.7, 115.8, 118.7, 121.6, 124.4, 127.2, 129.9, 132.6, 135.3, 138.,
#         140.7, 143.4, 146.2, 148.9, 151.7, 154.3, 157.,  159.7, 162.4, 165.1, 167.8, 170.5,
#         173.2, 175.9, 178.5, 181.2, 183.8, 186.5, 189.2, 191.9, 194.6, 197.3, 200.,  202.7,
#         205.4, 208.1, 210.8, 213.5, 216.2, 218.9, 221.6, 224.2, 226.9, 229.6, 232.3, 235.,
#         237.7, 240.4, 243.1, 245.8, 248.5, 251.2, 253.9, 256.6, 259.3, 262.,  264.6, 267.3,
#         270.,  272.7, 275.4, 278.,  280.7, 283.4, 286.1, 288.8, 291.5, 294.2, 297.,  299.7]
#
# idx_max = 0
# for i in range(len(LTP)):
#     if LTP[i] == max(LTP):
#         idx_max = i
#
# plt.figure()
# plt.plot(gr_sp, LTP, label="LTP")
# plt.title("LTP curve")
# plt.axvline(x= 100, linestyle='--', color="r", label="Time window")
# plt.axvline(x= 300, linestyle='--', color="r")
# plt.plot(gr_sp[idx_max], max(LTP), marker='o', label="Max value", color="g", markersize = 7)
# plt.annotate("[t_Gr, LTP_value] = [256.6, 0.999797]", (gr_sp[idx_max], max(LTP)), (135, max(LTP)))
# plt.legend()
# plt.show()