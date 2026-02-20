import json
import os

from utils import *

receptor_dict = {
    "1": ["AMPA", "I_syn_ampa", "AMPA_E_rev"],
    "2": ["NMDA", "I_syn_nmda", "NMDA_E_rev"],
    "3": ["GABA", "I_syn_gaba", "GABA_E_rev"],
    "4": ["AMPA2", "I_syn_ampa2", "AMPA2_E_rev"],
    "5": ["AMPA3", "I_syn_ampa3", "AMPA3_E_rev"],
}

if __name__ == "__main__":
    cell_name = "SC"
    params_SC = {
        "t_ref": 1.59,
        "C_m": 14.6,
        "V_th": -53,
        "V_reset": -78,
        "E_L": -68,
        "V_m": -68,
        "tau_m": 9.125,
        "I_e": 11.148339491906636,
        "k_adap": 0.17534246575342463,
        "k_2": 0.11229818783162836,
        "A1": 209.4867066866642,
        "A2": 163.99782404683813,
        "k_1": 0.15727774426530935,
        # pf-SC AMPA
        "AMPA_E_rev": 0.0,
        "AMPA_g_init": 0.09324813802356961,
        "AMPA_Tau_r": 0.832811417468158,
        "AMPA_Tau_d1": 0.16858899559030588,
        "AMPA_Tau_d2": 1.1220894106498769,
        "AMPA_A_r": 4.881100222623705,
        "AMPA_A1": 17.07322080857238,
        "AMPA_A2": 14.764149216255884,
        # pf-SC NMDA
        "NMDA_E_rev": -3.7,
        "NMDA_g_init": 0.01514802993681726,
        "NMDA_Tau_r": 37.622573319469055,
        "NMDA_Tau_d1": 53.15105534517695,
        "NMDA_Tau_d2": 892.4760819720726,
        "NMDA_A_r": 4.743977499613861,
        "NMDA_A1": 2.854964381660505,
        "NMDA_A2": 0.35317007422511043,
        # SC-SC GABA
        "GABA_E_rev": -65.0,
        "GABA_g_init": 1.0230134014130856,
        "GABA_Tau_r": 0.31230551452717137,
        "GABA_Tau_d1": 2.9081700675470903,
        "GABA_Tau_d2": 29.3943338444682,
        "GABA_A_r": 3.8423880947398508,
        "GABA_A1": 6.293522738150232,
        "GABA_A2": 4.524929470710231,
    }

    # Static synapses for SC
    pf_SC_ampa = {"weight": 1, "delay": 0.1, "receptor_type": 1}
    pf_SC_nmda = {"weight": 1, "delay": 0.1, "receptor_type": 2}
    SC_SC_gaba = {"weight": 1, "delay": 0.1, "receptor_type": 3}

    # TM params for SC synapses
    pf_SC_ampa_tm = {
        "u": 0.0,
        "y": 0.0,
        "x": 1.0,
        "U": 0.15,
        "tau_rec": 35.1,
        "tau_fac": 10.8,
        "tau_psc": 1,
        "weight": 1,
        "delay": 0.1,
        "receptor_type": 1,
    }
    pf_SC_nmda_tm = {
        "u": 0.0,
        "y": 0.0,
        "x": 1.0,
        "U": 0.15,
        "tau_rec": 35.1,
        "tau_fac": 10.8,
        "tau_psc": 1,
        "weight": 1,
        "delay": 0.1,
        "receptor_type": 2,
    }
    SC_SC_gaba_tm = {
        "u": 0.0,
        "x": 1.0,
        "y": 0.0,
        "U": 0.42,
        "tau_rec": 38.7,
        "tau_fac": 4.0,
        "tau_psc": 1,
        "weight": 1,
        "delay": 0.1,
        "receptor_type": 3,
    }

    run_static_vs_tm_test(
        cell_name=cell_name,
        cell_params=params_SC,
        receptor_type=str(1),
        receptor_dict=receptor_dict,
        syn_static=pf_SC_ampa,
        syn_tm=pf_SC_ampa_tm,
        syn_base_name="pf_SC_ampa",
    )

    run_static_vs_tm_test(
        cell_name=cell_name,
        cell_params=params_SC,
        receptor_type=str(2),
        receptor_dict=receptor_dict,
        syn_static=pf_SC_nmda,
        syn_tm=pf_SC_nmda_tm,
        syn_base_name="pf_SC_nmda",
    )

    run_static_vs_tm_test(
        cell_name=cell_name,
        cell_params=params_SC,
        receptor_type=str(3),
        receptor_dict=receptor_dict,
        syn_static=SC_SC_gaba,
        syn_tm=SC_SC_gaba_tm,
        syn_base_name="SC_SC_gaba",
    )
