from bsb import from_storage

from cerebellum.analysis.spiking_results import (
    BasicSimulationReport,
    FiringRatesPlot,
    RasterPSTHPlot,
    SimResultsTable,
)
from cerebellum.analysis.structure_analysis import StructureReport
from cerebellum.analysis.report import PlotTypeInfo

reco_file = "mouse_cereb_microzones_nest.hdf5"
nio_folder = "nio_files/micro_stim"
reco_pdf = "report_microzones.pdf"
sim_basal_pdf = "report_micro_sim_basal.pdf"
sim_stim_pdf = "report_micro_sim_stim.pdf"
scaffold = from_storage(reco_file)


LIST_CT_INFO = [
    PlotTypeInfo("mossy_fibers", [0.847, 0, 0.451, 1.0], "mf"),
    PlotTypeInfo("glomerulus", [0.847, 0, 0.451, 1.0], "glom"),
    PlotTypeInfo("granule_cell", [0.7, 0.15, 0.15, 0.5], "GrC"),
    PlotTypeInfo("ascending_axon", [0.7, 0.15, 0.15, 0.5], "aa"),
    PlotTypeInfo("parallel_fiber", [0.7, 0.15, 0.15, 0.5], "pf"),
    PlotTypeInfo("unipolar_brush_cell", [0.196, 0.808, 0.988, 1.0], "UBC"),
    PlotTypeInfo("ubc_glomerulus", [0.196, 0.808, 0.988, 1.0], "ubc_glom"),
    PlotTypeInfo("golgi_cell", [0, 0.45, 0.7, 1.0], "GoC"),
    PlotTypeInfo("purkinje_cell", [0.275, 0.800, 0.275, 1.0], "PC+"),
    PlotTypeInfo("purkinje_cell_minus", [0.275, 0.550, 0.275, 1.0], "PC-"),
    PlotTypeInfo("basket_cell", [1, 0.647, 0, 1.0], "BC"),
    PlotTypeInfo("stellate_cell", [1, 0.84, 0, 1.0], "SC"),
    PlotTypeInfo("dcn_p_plus", [0.3, 0.3, 0.3, 1.0], "DCN_P+"),
    PlotTypeInfo("dcn_p_minus", [0.1, 0.1, 0.1, 1.0], "DCN_P-"),
    PlotTypeInfo("dcn_i_plus", [0.635, 0, 0.145, 1.0], "DCN_I+"),
    PlotTypeInfo("dcn_i_minus", [0.435, 0, 0.145, 1.0], "DCN_I-"),
    PlotTypeInfo("io_plus", [0.46, 0.376, 0.54, 1.0], "IO+"),
    PlotTypeInfo("io_minus", [0.76, 0.276, 0.74, 1.0], "IO-"),
]


# RECONSTRUCTION REPORT
report_struct = StructureReport(reco_file, cell_type_info=LIST_CT_INFO)
report_struct.print_report(reco_pdf)

# SIMULATION REPORT
report_sim_stim = BasicSimulationReport(scaffold,simulation_name='mf_cf_stimulus',folder_nio=nio_folder, cell_types_info=LIST_CT_INFO)
report_sim_stim.print_report(sim_stim_pdf)

# SIMULATION REPORT - mossy fiber stimulus
# report_sim_stim = BasicSimulationReport(
#     scaffold, simulation_name="mf_stimulus", folder_nio=nio_folder, time_from=1200, time_to=1250
# )
# table = SimResultsTable(
#     fig_size=(10, 5),
#     scaffold=report_sim_stim.scaffold,
#     simulation_name=report_sim_stim.simulation_name,
#     time_from=report_sim_stim.time_from,
#     time_to=report_sim_stim.time_to,
#     all_spikes=report_sim_stim.all_spikes,
#     nb_neurons=report_sim_stim.nb_neurons,
#     populations=report_sim_stim.populations,
# )
# table.set_axis_off()
# table.save_figure("table_sim_stim_dcn_io.png", dpi=200)

# report_sim_stim = BasicSimulationReport(
#     scaffold,
#     simulation_name="mf_stimulus",
#     folder_nio=nio_folder
# )
# firing = RasterPSTHPlot(
#     fig_size=(10, 5),
#     scaffold=report_sim_stim.scaffold,
#     simulation_name=report_sim_stim.simulation_name,
#     time_from=report_sim_stim.time_from,
#     time_to=report_sim_stim.time_to,
#     all_spikes=report_sim_stim.all_spikes,
#     nb_neurons=report_sim_stim.nb_neurons,
#     populations=report_sim_stim.populations,
#     dict_colors=report_sim_stim.colors,
# )
#
# firing.save_figure("raster_sim_stim.png", dpi=200)
