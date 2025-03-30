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

    params_dcnp = {"t_ref": 1.5,
          "C_m": 142,
          "V_th": -36,
          "V_reset": -55,
          "E_L": -45,
          "I_e": 75.385,
          "V_m": -45,
          "lambda_0": 3.0,
          "tau_V": 3.0,
          "tau_m": 33,
          "k_adap": 0.408,
          "k_1": 0.697,
          "k_2": 0.047,
          "A1": 13.857,
          "A2": 3.477,
          "tau_syn1": 1.0,
          "tau_syn2": 0.7,
          "E_rev1": 0,
          "E_rev2": -80}

    dcn_p = nest.Create("eglif_dcnp", 1, params=params_dcnp)
    pc = nest.Create("eglif_pc_nestml", 1, params=params_pc)
    mf = nest.create("parrot_neuron", 1)
    wr = nest.Create("weight_recorder")
    nest.CopyModel("stdp_synapse_cosexp", "stdp_rec", {"weight_recorder": wr})
    nest.Connect(pc, dcn_p, syn_spec={"synapse_model" : "static_synapse" ,"weight": 1.25, "receptor_type" : 2})
    nest.Connect(mf, dcn_p, syn_spec={"synapse_model" : "stdp_rec", "receptor_type" : 1})

    mf_gen = nest.Create("spike_generator", params={"spike_times": [50, 100, 150, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 250]})
    pc_gen = nest.Create("poisson_generator", params={"rate": 500, "start": 200, "stop": 210})

    sr_mf, sr_pc, sr_dcn = nest.Create("spike_recorder", 3)
    nest.Connect(mf, sr_mf)
    nest.Connect(pc, sr_pc)
    nest.Connect(dcn_p, sr_dcn)
    nest.Connect(mf_gen, mf, syn_spec={"receptor_type" : 1, "weight" : 20})
    nest.Connect(pc_gen, pc, syn_spec={"receptor_type" : 1,"weight" : 1})
    #nest.Connect(io_generator, io, syn_spec={"weight" : 1})
    #nest.Connect(gr_2_gen, gr, syn_spec={"receptor_type": 5})


    mf_dcn_conns = nest.GetConnections(synapse_model="stdp_rec")

    # pc_pot = nest.Create("voltmeter")
    # nest.Connect(pc_pot, pc)

    # nest.Simulate(1000)
    # nest.voltage_trace.from_device(pc_pot)
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
        weights[i] = nest.GetStatus(mf_dcn_conns, "weight")


    return sr_mf.events, sr_pc.events, sr_dcn.events, weights


mf, pc, dcn, w = run_simulation()
print(np.diff(w, 0))
print("IO SPIKES: ", mf["times"])
print("GR SPIKES: ", pc["times"])
print("MLI SPIKES: ", dcn["times"])
plt.scatter(mf["times"], mf["senders"], c="b", s=7, label = "MF")
plt.scatter(pc["times"], pc["senders"], c="r", s=7, label = "PC")
plt.scatter(dcn["times"], dcn["senders"], c="g", s=7, label = "DCN")
plt.title("Scatter plot for plasticity test")
plt.axvline(x= 200, linestyle='--', color="k")
plt.axvline(x= 210, linestyle='--', color="k")
plt.legend(loc="best")
i =  np.zeros((999, 1))
counter = 0
for z in range(999):
    i[z] = counter
    counter = counter +1
plt.figure()
plt.plot(i, np.diff(w, axis =0))
plt.title("Weight diff")
plt.axvline(x= 200, linestyle='--', color="k")
plt.axvline(x= 210, linestyle='--', color="k")

#plt.xlim(0,600)
plt.show()

plt.figure()
plt.plot(np.append(i, 1000), w)
plt.title("Weight evolution")
plt.axvline(x= 200, linestyle='--', color="k")
plt.axvline(x= 210, linestyle='--', color="k")
#plt.xlim(0,600)
plt.show()
