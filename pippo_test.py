from bsb import from_storage
import numpy as np


scaffold = from_storage("mouse_cerebellum_microzones.hdf5")
print(scaffold.get_connectivity_sets())

# cs = scaffold.connectivity['io_to_purkinje_plus']
# a = cs.get_output_names()
# b = scaffold.get_connectivity_set(a[0])
# dio = b.load_connections().as_globals().all()
# print("aaaaaaaa: ", dio)

#PER VEDERE GLI ID DELLE CONNESSIONI, TI PRINTI LA TABELLA DEL CONNECTIVITY SET CON TUTTE LE COPPIE --> cs.load_connections().as_globals().all()
#POI TI PRINTI CON GET_LABELLED, SPECIFICANDO LA LABEL, I NEURONI CON QUELLA LABEL E VEDI SE COMBACIANO

# for name,ct in scaffold.cell_types.items():
#     if name == "io":
#         ct1= ct
#     elif name == "purkinje_cell":
#         ct2=ct
#
# print(cs.get_output_names(ct1,ct2))



# cs_m = scaffold.connectivity['io_to_purkinje_minus']
# a_m = cs_m.get_output_names()
# b_m = scaffold.get_connectivity_set(a_m[0])
# dio_m = b_m.load_connections().as_globals().all()
# print("aaaaaaaa: ", dio_m)
# dio = np.asarray(dio)
# dio_m = np.asarray(dio_m)
# ids = np.hstack([dio[:,:,0], dio_m[:,:,0]])
#
# ids_new = np.where(np.isin(ids[0,:], ids[1,:]))[0]

for type_name, cell_type in scaffold.cell_types.items():
    ps = cell_type.get_placement_set()
    if type_name == "purkinje_cell":
        pc_plus = ps.get_labelled(["plus"])
        pc_minus = ps.get_labelled(["minus"])
        pos = ps.load_positions()
        print("PC+:", pc_plus)
        print("PC-:", pc_minus)

for type_name, cell_type in scaffold.cell_types.items():
    ps = cell_type.get_placement_set()
    if type_name == "io":
        io_plus = ps.get_labelled(["plus"])
        io_minus = ps.get_labelled(["minus"])
        pos = ps.load_positions()
        print("IO+: ", io_plus)
        print("IO-: ", io_minus)

cs = scaffold.connectivity['io_to_mli_plus']
a = cs.get_output_names()
b = scaffold.get_connectivity_set(a[0])
dio = b.load_connections().as_globals().all()
print("aaaaaaaa: ", dio)
cs_m = scaffold.connectivity['io_to_mli_minus']
a_m = cs_m.get_output_names()
b_m = scaffold.get_connectivity_set(a_m[0])
dio_m = b_m.load_connections().as_globals().all()
print("aaaaaaaa: ", dio_m)
dio = np.asarray(dio)
dio_m = np.asarray(dio_m)
ids = np.hstack([dio[:,:,0], dio_m[:,:,0]])

ids_new = np.where(np.isin(ids[0,:], ids[1,:]))[0]




























# for name,ct in scaffold.cell_types.items():
#     if name == "io":
#         ct1= ct
#     elif name == "purkinje_cell":
#         ct2=ct
#
# print(cs_m.get_output_names(ct1,ct2))
#
#
# cs_bp = scaffold.connectivity['basket_to_purkinje_plus']
#
# for name,ct in scaffold.cell_types.items():
#     if name == "basket_cell":
#         ct1= ct
#     elif name == "purkinje_cell":
#         ct2=ct
#
# print(cs_bp.get_output_names(ct1,ct2))
#
#
# cs_bm = scaffold.connectivity['basket_to_purkinje_minus']
#
# for name,ct in scaffold.cell_types.items():
#     if name == "basket_cell":
#         ct1= ct
#     elif name == "purkinje_cell":
#         ct2=ct
#
# print(cs_bm.get_output_names(ct1,ct2))
#
#
# cs_sp = scaffold.connectivity['stellate_to_purkinje_plus']
#
# for name,ct in scaffold.cell_types.items():
#     if name == "stellate_cell":
#         ct1= ct
#     elif name == "purkinje_cell":
#         ct2=ct
#
# print(cs_sp.get_output_names(ct1,ct2))
#
#
# cs_sm = scaffold.connectivity['stellate_to_purkinje_minus']
#
# for name,ct in scaffold.cell_types.items():
#     if name == "stellate_cell":
#         ct1= ct
#     elif name == "purkinje_cell":
#         ct2=ct
#
# print(cs_sm.get_output_names(ct1,ct2))
#
#
# cs_mlip = scaffold.connectivity['io_to_mli_plus']
#
# for name,ct in scaffold.cell_types.items():
#     if name == "io":
#         ct1= ct
#     elif name == "basket_cell":
#         ct2=ct
#
# print(cs_mlip.get_output_names(ct1,ct2))
#
#
# cs_mlim = scaffold.connectivity['io_to_mli_minus']
#
# for name,ct in scaffold.cell_types.items():
#     if name == "io":
#         ct1= ct
#     elif name == "basket_cell":
#         ct2=ct
#
# print(cs_mlim.get_output_names(ct1,ct2))