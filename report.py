from os.path import join

from cerebellum.analysis import BasicSimulationReport, StructureReport

root_path = "./"
cereb_file = join(root_path, "mouse_cereb_alice_pc_to_dcnp_1_5_40redPC.hdf5")
nio_folder = join(root_path, "nio_files", "pc_to_dcnp_1_5_40redPC")

report = StructureReport(cereb_file)
# report.plots["placement_3d"].show()
report.print_report("report_struc_dcn_pc_to_dcnp_1_5_40redPC.pdf")
report = BasicSimulationReport(cereb_file, "basal_activity", nio_folder)
report.print_report("report_sim_dcn_pc_to_dcnp_1_5_40redPC.pdf")
