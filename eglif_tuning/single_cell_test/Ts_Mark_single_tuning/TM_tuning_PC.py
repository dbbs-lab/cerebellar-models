import json
import os

from utils import *

receptor_dict = {
    "1": ["AMPA", "I_syn_ampa", "AMPA_E_rev"],
    "2": ["NMDA", "I_syn_nmda", "NMDA_E_rev"],
    "3": ["GABA", "I_syn_gaba", "GABA_E_rev"],
}


if __name__ == "__main__":
    cell_name = "PC"
    params_PC = {
        "C_m": 334.0,
        "E_L": -59.0,
        "V_m": -59.0,
        "V_min": -350.0,
        "V_reset": -69.0,
        "V_th": -43.0,
        "t_ref": 0.5,
        "tau_m": 47.0,
        "A1": 1945.2104835367854,
        "A2": 4.962674236410129,
        "I_e": 158.59305781086863,
        "k_1": 0.2882491796788571,
        "k_2": 0.03495569814544013,
        "k_adap": 0.4807699301937877,
        # aa-PC & pf-PC AMPA
        "AMPA_E_rev": 0.0,
        "AMPA_g_init": 0.19294907221534321,
        "AMPA_Tau_r": 0.686334499921023,
        "AMPA_Tau_d1": 0.08637409507467056,
        "AMPA_Tau_d2": 0.10086052627970157,
        "AMPA_A_r": 0.7097862389833014,
        "AMPA_A1": 4.011865516012454,
        "AMPA_A2": 3.988541673612559,
        # SC-PC & BC-PC GABA
        "GABA_E_rev": -70.0,
        "GABA_g_init": 2.0879285694406864,
        "GABA_Tau_r": 0.7622223660610657,
        "GABA_Tau_d1": 2.814505586125079,
        "GABA_Tau_d2": 28.744462070299914,
        "GABA_A_r": 7.382885340300089,
        "GABA_A1": 1.7476237681931543,
        "GABA_A2": 1.2018061113172258,
    }
    syn_aa_PC_ampa = {"weight": 1.0, "delay": 0.1, "receptor_type": 1}
    syn_aa_PC_ampa_tm = {
        "U": 0.13,
        "u": 0,
        "y": 0,
        "x": 1.0,
        "tau_rec": 35.1,
        "tau_fac": 54.0,
        "tau_psc": 1.0,
        "delay": 0.1,
        "receptor_type": 1,
        "weight": 1.0,
    }
    syn_BC_PC_gaba = {"weight": 1.0, "delay": 0.1, "receptor_type": 3}
    syn_BC_PC_gaba_tm = {
        "U": 0.35,
        "u": 0,
        "y": 0,
        "x": 1.0,
        "tau_rec": 15.0,
        "tau_fac": 4.0,
        "tau_psc": 1.0,
        "receptor_type": 3,
        "weight": 1.0,
        "delay": 0.1,
    }
    syn_pf_PC_ampa = {"weight": 1.0, "delay": 0.1, "receptor_type": 1}
    syn_pf_PC_ampa_tm = {
        "U": 0.13,
        "u": 0,
        "y": 0,
        "x": 1.0,
        "tau_rec": 35.1,
        "tau_fac": 54.0,
        "tau_psc": 1.0,
        "receptor_type": 1,
        "weight": 1.0,
        "delay": 0.1,
    }
    syn_SC_PC_gaba = {"weight": 1.0, "delay": 0.1, "receptor_type": 3}
    syn_SC_PC_gaba_tm = {
        "U": 0.35,
        "u": 0,
        "y": 0,
        "x": 1.0,
        "tau_rec": 15.0,
        "tau_fac": 4.0,
        "tau_psc": 1.0,
        "receptor_type": 3,
        "weight": 1.0,
        "delay": 0.1,
    }

    run_static_vs_tm_test(
        cell_name=cell_name,
        cell_params=params_PC,
        receptor_type=str(1),
        receptor_dict=receptor_dict,
        syn_static=syn_aa_PC_ampa,
        syn_tm=syn_aa_PC_ampa_tm,
        syn_base_name="aa_PC_ampa",
    )

    run_static_vs_tm_test(
        cell_name=cell_name,
        cell_params=params_PC,
        receptor_type=str(1),
        receptor_dict=receptor_dict,
        syn_static=syn_pf_PC_ampa,
        syn_tm=syn_pf_PC_ampa_tm,
        syn_base_name="pf_PC_ampa",
    )

    run_static_vs_tm_test(
        cell_name=cell_name,
        cell_params=params_PC,
        receptor_type=str(3),
        receptor_dict=receptor_dict,
        syn_static=syn_BC_PC_gaba,
        syn_tm=syn_BC_PC_gaba_tm,
        syn_base_name="BC_PC_gaba",
    )

    run_static_vs_tm_test(
        cell_name=cell_name,
        cell_params=params_PC,
        receptor_type=str(3),
        receptor_dict=receptor_dict,
        syn_static=syn_SC_PC_gaba,
        syn_tm=syn_SC_PC_gaba_tm,
        syn_base_name="SC_PC_gaba",
    )
