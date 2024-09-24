import report
from cerebellum.analysis.structure_analysis import StructureReport
from cerebellum.analysis.spiking_results import BasicSimulationReport
from bsb import from_storage
from cerebellum.analysis.spiking_results import SimResultsTable, FiringRatesPlot

reco_file = "mouse_cereb_dcn_nest_copy.hdf5"
nio_folder = "res_nio/res_stim_dcn_nest_copy"
reco_pdf = "report_dcn_nest_copy.pdf"
sim_basal_pdf = "report_dcn_sim_basal_copy.pdf"
sim_stim_pdf = "report_dcn_sim_stim_copy.pdf"
scaffold = from_storage(reco_file)
# RECONSTRUCTION REPORT
# report_struct = StructureReport(reco_file)
# report_struct.print_report(reco_pdf)

# SIMULATION REPORT - basal activity
# report_sim_bas = BasicSimulationReport(scaffold,simulation_name='basal_activity',folder_nio=nio_folder)
# report_sim_bas.print_report(sim_basal_pdf)

# SIMULATION REPORT - mossy fiber stimulus
report_sim_stim = BasicSimulationReport(scaffold, simulation_name='mf_stimulus', folder_nio=nio_folder,
                                        time_from=1200, time_to=1250)
table = SimResultsTable(
    fig_size=(10,5),
    scaffold = report_sim_stim.scaffold,
    simulation_name= report_sim_stim.simulation_name,
    time_from= report_sim_stim.time_from,
    time_to= report_sim_stim.time_to,
    all_spikes=report_sim_stim.all_spikes,
    nb_neurons = report_sim_stim.nb_neurons,
    populations= report_sim_stim.populations
)
table.set_axis_off()
table.save_figure('table_sim_stim_copy.png', dpi=200)

report_sim_stim = BasicSimulationReport(scaffold, simulation_name='mf_stimulus', folder_nio=nio_folder)
firing = FiringRatesPlot(
    fig_size=(10,5),
    scaffold = report_sim_stim.scaffold,
    simulation_name= report_sim_stim.simulation_name,
    time_from= 0,
    time_to= 5000,
    all_spikes=report_sim_stim.all_spikes,
    nb_neurons = report_sim_stim.nb_neurons,
    populations= report_sim_stim.populations,
    dict_colors=report_sim_stim.colors
)

firing.save_figure('firing_sim_stim_copy.png', dpi=200)
