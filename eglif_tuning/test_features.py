import matplotlib.pyplot as plt
import nest

from cerebellar_models.optimization.features import *


def nest_protocol(cell_params, currents, start_stim, stop_stim, duration):
    def nest_simulation(cp, curr, t_i, t_f, T):
        nest.ResetKernel()
        nest.SetKernelStatus({"local_num_threads": 1})
        nest.Install("cerebmodule")
        nest.resolution = 0.1

        cell = nest.Create("eglif_multirec", 1, params=cp)

        gen = nest.Create("dc_generator")
        nest.SetStatus(gen, {"amplitude": curr, "start": t_i, "stop": t_f})

        sr = nest.Create("spike_recorder")

        nest.Connect(gen, cell)
        nest.Connect(cell, sr)

        nest.Simulate(T)

        spike_times = nest.GetStatus(sr, "events")[0]["times"]

        return spike_times

    return {I: nest_simulation(cell_params, I, start_stim, stop_stim, duration) for I in currents}


if __name__ == "__main__":
    PLOT = False

    ### 1. TEST FEATURES EXTRACTION FROM MULTICOMPARTIMENTAL NEURON EXPERIMENTS ###
    feature_names = ["peak_time", "time_to_first_spike", "time_to_second_spike", "mean_frequency"]
    threshold_GrC = -41.0  # Example for GrC
    data_folder = os.path.join(os.path.dirname(__file__), "tofit_eglif/results_tofitEglif/GrC/")
    neuron_features = multicomp_features(
        data_folder=data_folder,
        threshold=threshold_GrC,
        features=feature_names,
        start_stim=100,
        end_stim=600,
    )
    print("Neuron experiments current steps: ", neuron_features["current"].values)

    current_values = 10.0
    selected_features = neuron_features[neuron_features["current"] == current_values]
    print("Test one feature: ", selected_features["time_to_first_spike"].values)

    if PLOT:
        currents = neuron_features["current"].values
        freqs = neuron_features["mean_frequency"].fillna(0).values

        x = np.linspace(0, max(currents), 100)

        def line(x, m, q):
            return m * x + q

        plt.figure()
        plt.scatter(currents, freqs)
        plt.show()

    ### 2.  TEST FEATURES EXTRACTION FROM NEST EXPERIMENTS ###
    currents = neuron_features["current"].values
    cell_params = {
        "t_ref": 1.5,
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
    }
    # Example for GrC
    nest_results = nest_protocol(cell_params, currents, start_stim=100, stop_stim=600, duration=700)
    spike_trains = []
    for I in currents:
        spike_trains.append(np.array(nest_results[I]))

    nest_features = point_neuron_features(
        spike_trains=spike_trains,
        currents=currents,
        features=feature_names,
        start_stim=100,
        end_stim=600,
    )

    print("Nest experiments current steps: ", nest_features["current"].values)
    selected_features = nest_features[nest_features["current"] == current_values]
    print("Test one feature ", selected_features["time_to_first_spike"].values)
