from bsb.services.mpi import MPIService
from bsb_nest import NestAdapter

from bsb import from_storage
from cerebellum.analysis.spiking_results import SimResultsTable

import numpy as np
from elephant.statistics import isi
# visualization
import matplotlib.pyplot as plt
import torch
import matplotlib
import corner
# sbi
from sbi import analysis as analysis
from sbi import utils as utils
from sbi.inference import NPE, simulate_for_sbi
from sbi.utils.user_input_checks import (
    check_sbi_inputs,
    process_prior,
    process_simulator,
)
# from bsb.services import MPI
from mpi4py import MPI

def sim_core(weights):
    weights = weights.numpy()
    weights = np.append(1, weights)


    # Clear NEST
    adapter.reset_kernel()
    #print(f"ADAPTER RANK: {mpi.Get_rank()}")
    # Set the weights
    for i, conn in enumerate(simulation.connection_models):
        print(simulation.connection_models[conn].name)
        simulation.connection_models[conn].synapse.weight = weights[i]
    # Let the adapter translate the simulation config into
    # simulator specific instructions
    adapter.prepare(simulation)

    # Let the adapter run the simulation and collect the output.
    results = adapter.run(simulation)
    results = adapter.collect(results)[0]
    # Extract spiketrains
    #if mpi.Get_rank() == 0:
    nb_neurons = np.zeros(len(populations), dtype=int)
    order = np.zeros(len(populations), dtype=int)
    for j, st in enumerate(results.spiketrains):
        cell_type = st.annotations["device"].split("_rec")[0]
        if cell_type not in populations:
            cell_type += "_cell"
        if cell_type in populations:
            i = populations.index(cell_type)
            nb_neurons[i] = st.annotations["pop_size"]
            order[i] = j
    all_spikes = [results.spiketrains[i] for i in order]
    print(f"Populations: {populations} | Order: {order}")

    # Compute firing rates
    sim_res = SimResultsTable(
        (10, 10),
        scaffold,
        simulation_name,
        0,
        1000,
        all_spikes,
        nb_neurons,
        populations,
    )
    sim_res.update()
    sim_results = [(np.mean(fr) if len(fr) > 0 else 0.) for fr in sim_res.get_firing_rates().values()]
    # mean_isi = [np.mean(result_isi) if len(result_isi) > 0 else 0. for result_isi in sim_res.get_isis_values().values()]
    # sim_results = np.append(fr_results, mean_isi)
    #results = torch.as_tensor(sim_results)

    return sim_results, sim_res

def sim_wrapper(weights):
    weights = weights.numpy()
    results, _ = sim_core(weights)

    return results

if __name__ == "__main__":
    simulation_name = "basal_activity"
    scaffold = from_storage("mouse_cereb_nest_vitro.hdf5")
    populations = sorted(list(set(scaffold.cell_types.keys()) - {"glomerulus", "mossy_fibers"}))
    simulation = scaffold.get_simulation(simulation_name)

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()


    # mpi = MPI.COMM_WORLD
    adapter = NestAdapter(comm = comm)
    # adapter = NestAdapter()

    # defining range of priors
    #prior_min = np.zeros((1, 14))
    # prior_max = np.ones((1,14))
    #prior_max = [0.6, 0.6, 0.6, 0.01, 1, 0.1, 0.1, 0.6, 0.1, 0.01, 0.6, 0.1, 0.01, 0.6]
    n_sims = 1000 #25
    prior_min = np.zeros((1,14))
    #prior_max = (5 * np.ones((1,14)))
    prior_max = np.ones((1,14))
    priors = utils.torchutils.BoxUniform(low=torch.as_tensor(prior_min), high=torch.as_tensor(prior_max))

    if rank == 0:
        samples = priors.sample((n_sims,))
    else:
        samples = None

    samples = comm.bcast(samples, root=0)

    local_indices = list(range(rank, n_sims, size))
    local_results = []

    for i in local_indices:
        result, _ = sim_core(samples[i])
        local_results.append(result)

    local_results = torch.tensor(local_results, dtype=torch.float32)
    all_results = comm.gather(local_results, root=0)

    if rank == 0:
        x = torch.cat(all_results, dim=0)
        samples = torch.squeeze(samples)
        print(f"Completed {n_sims} simulations. Final result shapes: x = {x.shape} || theta = {samples.shape}")
        priors, num_parameters, prior_returns_numpy = process_prior(priors)
        simulator = process_simulator(sim_wrapper, priors, prior_returns_numpy)
        inference = NPE(prior=priors)

        density_estimator = inference.append_simulations(samples, x).train(show_train_summary=True)
        posterior = inference.build_posterior(density_estimator)

        print("POSTERIOR: ", posterior)

        # Saving all data

        torch.save({
            "theta": samples,
            "x": x,
            "density_estimator": density_estimator,
            "posterior": posterior,
            "n_sims": n_sims,
            "timestamp": "2025-06-30",
        }, "sbi_results_1000.pt")

    # ----- SEACH FOR THE BEST PARAMETER SET W/ OBSERVATIONAL DATA -----
    # obs_mean_fr = torch.tensor([38.846666666666664, 25.577464788732396, 13.724188371053595, 50.88405797101449, 27.93],dtype=torch.float32)
    # obs = torch.tensor([16, 11, 3.4, 47, 13, 77, 120, 330, 22, 180], dtype=torch.float32)
    obs = torch.tensor([16, 11, 3.4, 47, 13], dtype=torch.float32)
    # range = [(8.3, 23.7), (5.6, 16.4), (0, 6.8),
    # (48.0, 52.0),
    # (26.0, 30.0),
    # (32, 122),
    # (26, 214),
    #     (),
    #     (),
    #     ()
    # ]
    # std = np.array([7.7, 5.4, 3.4, 1.1, 11, 45, 94, 270, 0.48, 240])
    std = np.array([7.7, 5.4, 3.4, 1.1, 11])
    # Cell name	Mean Firing rate (Hz)	Mean ISI (ms) BASAL ACTIVITY IN VITRO
    # Mossy cell	4.0±1.4	250±140
    # Granule cell	3.4±3.4	330±270
    # Golgi cell	11±5.4	120±94
    # Purkinje cell	47±1.1	22±0.48
    # Basket cell	16±7.7	77±45
    # Stellate cell	13±11	180±240
    if rank == 0:
        # posterior_samples = posterior.sample((n_sims,), x=obs_mean_fr)
        posterior_samples = posterior.sample((n_sims,), x=obs)
    else:
        posterior_samples = None

    posterior_samples = comm.bcast(posterior_samples, root=0)

    distances = []
    all_sim_results = []
    mses = []
    weighted_mses = []
    pens = []

    for i in local_indices:
        # Simula ogni campione e calcola la distanza dai dati osservati
        #for theta_i in posterior_samples:
        #sim_result, _ = sim_core(theta_i)  # Output: torch.tensor([...])
        sim_result, _ = sim_core(posterior_samples[i])
        #sim_result = torch.tensor(sim_result, dtype=torch.float32)
        all_sim_results.append(sim_result)

        # Distanza Euclidea tra simulazione e dato osservato
        #distance = torch.norm(sim_result - obs_mean_fr)
        error = sim_result - obs.numpy()
        distance = np.linalg.norm(error)
        mse = np.mean(error**2)
        # penalties = np.where(np.abs(error) > std, 2 * np.abs(error), 0)
        # penalties = []
        # for k in range(error):
        #     if error[k] > std[k]:
        #         penalties.append(2 * np.abs(error[k]))
        penalties = np.where(np.abs(error) > std, 2 * np.abs(error), 0)
        weighted_mse = mse + np.sum(penalties)

        distances.append(distance)
        mses.append(mse)
        weighted_mses.append(weighted_mse)
        pens.append(penalties)



    distances = torch.tensor(distances, dtype=torch.float32)
    #distances = torch.squeeze(distances)
    all_distances = comm.gather(distances, root=0)
    all_sim_results = torch.tensor(all_sim_results, dtype=torch.float32)
    #all_sim_results = torch.squeeze(all_sim_results)
    total_sims = comm.gather(all_sim_results, root=0)
    mses = torch.tensor(mses, dtype=torch.float32)
    all_mses = comm.gather(mses, root=0)
    weighted_mses = torch.tensor(weighted_mses, dtype=torch.float32)
    all_weighted_mses = comm.gather(weighted_mses, root=0)
    pens = torch.tensor(pens, dtype=torch.float32)
    all_penalties = comm.gather(pens, root=0)

    if rank == 0:
        distances = torch.cat(all_distances, dim=0)
        posterior_samples = torch.squeeze(posterior_samples)
        total = torch.cat(total_sims, dim=0)
        mses = torch.cat(all_mses, dim=0)
        weighted_mses = torch.cat(all_weighted_mses, dim=0)
        pens = torch.cat(all_penalties, dim=0)
        # Trova il theta con distanza minima
        print(f"Completed {n_sims} samples. Final result shapes: x_obs = {total.shape} || theta = {posterior_samples.shape}")

        best_distance_idx = torch.argmin(distances)
        best_idx = torch.argmin(weighted_mses)
        print(f"Best index: {best_idx} and posterior samples: {posterior_samples.shape}")
        print(f"Best distance: {best_distance_idx} and relative theta: {posterior_samples[best_distance_idx]}")
        best_theta = posterior_samples[best_idx]
        best_sim_output = total[best_idx]
        # Stampa i risultati
        print("Miglior theta trovato:", best_theta)
        print("Output simulato con miglior theta:", best_sim_output)

        # Salva il miglior risultato trovato
        torch.save({
            "best_theta": best_theta,
            "obs": obs,
            "best_sim_output": best_sim_output,
            "posterior_samples": posterior_samples,
            "total_sims": total_sims,
            "distances": distances,
            "mses": mses,
            "weighted_mses": weighted_mses,
            "best_distance": posterior_samples[best_distance_idx],
            "total": total,
            "penalties": pens,
        }, "sbi_best_fit_1000.pt")

        # Visualizzazione dei risultati simulati vs osservati (boxplot)
        sim_matrix = np.vstack(total)
        labels = [f"Pop {i + 1}" for i in range(sim_matrix.shape[1])]

        matplotlib.use("Agg")

        plt.figure(figsize=(10, 6))
        plt.boxplot(total, labels=labels)
        plt.plot(range(1, len(obs) + 1), obs.numpy(), 'ro', label='Dati osservati')
        plt.plot(range(1, len(obs) + 1), best_sim_output.numpy(), 'bo', label='Migliori dati test')
        plt.ylabel("Firing rate (Hz)")
        plt.title("Posterior predictive check")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        # plt.show()
        plt.savefig("posterior_predictive_check_1000.png")

        # Corner plot dei parametri theta
        theta_matrix = posterior_samples.numpy()
        corner_labels = [f"θ{i + 1}" for i in range(theta_matrix.shape[1])]

        fig = corner.corner(
            theta_matrix,
            labels=corner_labels,
            truths=best_theta.numpy(),
            show_titles=True,
            title_fmt=".2f",
            title_kwargs={"fontsize": 10},
            label_kwargs={"fontsize": 10}
        )
        plt.suptitle("Posterior distribution of θ (corner plot)", fontsize=14)
        plt.tight_layout()
        # plt.show()
        fig.savefig("posterior_corner_plot_1000.png")





    # samples = priors.sample((n_sims,))
    # samples = torch.squeeze(samples)
    # #samples = torch.transpose(samples, 0, 1)
    # x = []
    # y = []
    #
    # for k in range(n_sims):
    #     y,_ = sim_core(samples[k])
    #     x.append(y)
    #
    # x = torch.tensor(x, dtype=torch.float32)

    # print("END!!!")
