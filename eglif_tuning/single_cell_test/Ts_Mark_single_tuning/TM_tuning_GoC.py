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
    cell_name = "GoC"
    params_GoC = {
        "t_ref": 2,
        "C_m": 145,
        "V_th": -55,
        "V_reset": -75,
        "E_L": -62,
        "V_min": -150,
        "V_m": -62,
        "tau_m": 44,
        "k_2": 0.022727272727272728,
        "I_e": 24.485214227491436,
        "k_adap": 0.10973532752297317,
        "k_1": 0.0683128131324206,
        "A1": 345.4705973413434,
        "A2": 94.19948918571492,
        # glom-GoC AMPA
        "AMPA_E_rev": 0.0,
        "AMPA_g_init": 1.115965180410803,
        "AMPA_Tau_r": 0.06828289728594421,
        "AMPA_Tau_d1": 1.7041177058103865,
        "AMPA_Tau_d2": 7.020524126553463,
        "AMPA_A_r": 1.886292737184908,
        "AMPA_A1": 7.63072301638014,
        "AMPA_A2": 2.5824424428033272,
        # aa-GoC AMPA
        "AMPA2_E_rev": 0.0,
        "AMPA2_g_init": 0.3545867855750201,
        "AMPA2_Tau_r": 0.07998662076204388,
        "AMPA2_Tau_d1": 0.8822965880874376,
        "AMPA2_Tau_d2": 1.3586665277793426,
        "AMPA2_A_r": 0.9852197655470364,
        "AMPA2_A1": 12.993176835484471,
        "AMPA2_A2": 2.1378395883863957,
        # aa-GoC & mf-GoC NMDA
        "NMDA_E_rev": -3.7,
        "NMDA_g_init": 0.06507491002133782,
        "NMDA_Tau_r": 7.492933042162157,
        "NMDA_Tau_d1": 49.78301745858933,
        "NMDA_Tau_d2": 784.4672149211721,
        "NMDA_A_r": 1.3691425668447825,
        "NMDA_A1": 2.511748085213512,
        "NMDA_A2": 0.3665220416477684,
        # pf-GoC AMPA
        "AMPA3_E_rev": 0.0,
        "AMPA3_g_init": 0.7261241439376229,
        "AMPA3_Tau_r": 0.06698111728388399,
        "AMPA3_Tau_d1": 0.8482199889588548,
        "AMPA3_Tau_d2": 1.9498307579186067,
        "AMPA3_A_r": 0.7566796251124199,
        "AMPA3_A1": 5.690348485137838,
        "AMPA3_A2": 0.6780028154381605,
        # GoC-GoC GABA
        "GABA_E_rev": -65.0,
        "GABA_g_init": 1.8938290878206254,
        "GABA_Tau_r": 0.1846418908724778,
        "GABA_Tau_d1": 2.639679024851566,
        "GABA_Tau_d2": 30.94312372047139,
        "GABA_A_r": 4.657518450694935,
        "GABA_A1": 8.885629800234906,
        "GABA_A2": 5.752371937918636,
    }

    # Static synapses to GoC
    aa_to_GoC_ampa = {"weight": 1, "delay": 0.1, "receptor_type": 4}
    aa_to_GoC_nmda = {"weight": 1, "delay": 0.1, "receptor_type": 2}
    mf_to_GoC_ampa = {"weight": 1, "delay": 0.1, "receptor_type": 1}
    mf_to_GoC_nmda = {"weight": 1, "delay": 0.1, "receptor_type": 2}
    pf_to_GoC_ampa = {"weight": 1, "delay": 0.1, "receptor_type": 5}

    # TM params for synapses to GoC
    aa_to_GoC_ampa_tm = {
        "u": 0,
        "y": 0,
        "x": 1,
        "U": 0.4,
        "tau_fac": 54,
        "tau_rec": 35.1,
        "tau_psc": 1.0,
        "weight": 1,
        "delay": 0.1,
        "receptor_type": 4,
    }
    aa_to_GoC_nmda_tm = {
        "u": 0,
        "y": 0,
        "x": 1,
        "U": 0.4,
        "tau_fac": 54,
        "tau_rec": 35.1,
        "tau_psc": 1.0,
        "weight": 1,
        "delay": 0.1,
        "receptor_type": 2,
    }
    mf_to_GoC_ampa_tm = {
        "u": 0,
        "y": 0,
        "x": 1,
        "U": 0.43,
        "tau_rec": 8.0,
        "tau_psc": 1.0,
        "tau_fac": 5,
        "weight": 1,
        "delay": 0.1,
        "receptor_type": 1,
    }
    mf_to_GoC_nmda_tm = {
        "u": 0,
        "y": 0,
        "x": 1,
        "U": 0.43,
        "tau_rec": 8.0,
        "tau_psc": 1.0,
        "tau_fac": 5,
        "weight": 1,
        "delay": 0.1,
        "receptor_type": 2,
    }
    pf_to_GoC_ampa_tm = {
        "u": 0.0,
        "y": 0.0,
        "x": 1,
        "U": 0.4,
        "tau_fac": 54,
        "tau_rec": 35.1,
        "tau_psc": 1.0,
        "weight": 1,
        "delay": 0.1,
        "receptor_type": 5,
    }

    run_static_vs_tm_test(
        cell_name=cell_name,
        cell_params=params_GoC,
        receptor_type=str(4),
        receptor_dict=receptor_dict,
        syn_static=aa_to_GoC_ampa,
        syn_tm=aa_to_GoC_ampa_tm,
        syn_base_name="aa_GoC_ampa",
    )

    run_static_vs_tm_test(
        cell_name=cell_name,
        cell_params=params_GoC,
        receptor_type=str(2),
        receptor_dict=receptor_dict,
        syn_static=aa_to_GoC_nmda,
        syn_tm=aa_to_GoC_nmda_tm,
        syn_base_name="aa_GoC_nmda",
    )

    run_static_vs_tm_test(
        cell_name=cell_name,
        cell_params=params_GoC,
        receptor_type=str(1),
        receptor_dict=receptor_dict,
        syn_static=mf_to_GoC_ampa,
        syn_tm=mf_to_GoC_ampa_tm,
        syn_base_name="mf_GoC_ampa",
    )

    run_static_vs_tm_test(
        cell_name=cell_name,
        cell_params=params_GoC,
        receptor_type=str(2),
        receptor_dict=receptor_dict,
        syn_static=mf_to_GoC_nmda,
        syn_tm=mf_to_GoC_nmda_tm,
        syn_base_name="mf_GoC_nmda",
    )

    run_static_vs_tm_test(
        cell_name=cell_name,
        cell_params=params_GoC,
        receptor_type=str(5),
        receptor_dict=receptor_dict,
        syn_static=pf_to_GoC_ampa,
        syn_tm=pf_to_GoC_ampa_tm,
        syn_base_name="pf_GoC_ampa",
    )
