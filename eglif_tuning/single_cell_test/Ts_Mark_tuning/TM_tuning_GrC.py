from utils import *

receptor_dict = {
    "1": ["AMPA", "I_syn_ampa", "AMPA_E_rev"],
    "2": ["NMDA", "I_syn_nmda", "NMDA_E_rev"],
    "3": ["GABA", "I_syn_gaba", "GABA_E_rev"],
    "4": ["AMPA2", "I_syn_ampa2", "AMPA2_E_rev"],
    "5": ["AMPA3", "I_syn_ampa3", "AMPA3_E_rev"],
}

if __name__ == "__main__":
    cell_name = "GrC"
    params_grc = {
        "t_ref": 1.5,
        "V_min": -150,
        "C_m": 7,
        "V_th": -41,
        "V_reset": -70,
        "E_L": -62,
        "I_e": 4.914404288084161,
        "tau_m": 24.15,
        "k_adap": 0.06746499109544667,
        "k_1": 0.16041441608311624,
        "k_2": 0.041407867494824016,
        "A1": 20.74311152075074,
        "A2": 4.1951490212728775,
        "AMPA_E_rev": 0.0,
        "AMPA_g_init": 1.3254051705870171,
        "AMPA_Tau_r": 0.9966611702091053,
        "AMPA_Tau_d1": 0.38598738306536895,
        "AMPA_Tau_d2": 3.7021239405389736,
        "AMPA_A_r": 6.973773179799838,
        "AMPA_A1": 1.4551614318130006,
        "AMPA_A2": 1.5956590886325035,
        "NMDA_E_rev": -3.7,
        "NMDA_g_init": 0.023436539864118,
        "NMDA_Tau_r": 21.192329383729344,
        "NMDA_Tau_d1": 14.860706054687927,
        "NMDA_Tau_d2": 121.08581080624421,
        "NMDA_A_r": 3.2972228905396754,
        "NMDA_A1": 2.951290940616272,
        "NMDA_A2": 0.21348858259846124,
        "GABA_E_rev": -65.0,
        "GABA_g_init": 0.16167187460716928,
        "GABA_Tau_r": 0.37595979577132965,
        "GABA_Tau_d1": 27.292807790043216,
        "GABA_Tau_d2": 269.0589627801647,
        "GABA_A_r": 4.540622636273912,
        "GABA_A1": 10.603535445257966,
        "GABA_A2": 3.144942430956719,
    }

    syn_mf_Grc_ampa = {"weight": 1., "delay": 0.1, "receptor_type": 1}
    syn_mf_Grc_nmda = {"weight": 1., "delay": 0.1, "receptor_type": 2}
    syn_GoC_GrC_gaba = {"weight": 1., "delay": 0.1, "receptor_type": 3}
    syn_mf_Grc_ampa_tm = {
        "U": 0.416,
        "u": 0.0,
        "y": 0.0,
        "x": 1,
        "tau_psc": 3.0,
        "tau_rec": 35.1,
        "tau_fac": 10.8,
        "weight": 1.0,
        "delay": 0.1,
        "receptor_type": 1,
    }
    syn_mf_Grc_nmda_tm = {
        "U": 0.416,
        "u": 0.0,
        "y": 0.0,
        "x": 1,
        "tau_psc": 3.0,
        "tau_rec": 35.1,
        "tau_fac": 10.8,
        "weight": 1.0,
        "delay": 0.1,
        "receptor_type": 2,
    }
    syn_GoC_GrC_gaba_tm = {
        "U": 0.35,
        "u": 0.0,
        "y": 0.0,
        "x": 1,
        "tau_rec": 36,
        "tau_fac": 58,
        "tau_psc": 1.0,
        "weight": 1.0,
        "delay": 0.1,
        "receptor_type": 3,
    }

    run_static_vs_tm_test(
        cell_name=cell_name,
        cell_params=params_grc,
        receptor_type=str(1),
        receptor_dict=receptor_dict,
        syn_static=syn_mf_Grc_ampa,
        syn_tm=syn_mf_Grc_ampa_tm,
        syn_base_name="mf_Grc_ampa",
    )

    run_static_vs_tm_test(
        cell_name=cell_name,
        cell_params=params_grc,
        receptor_type=str(2),
        receptor_dict=receptor_dict,
        syn_static=syn_mf_Grc_nmda,
        syn_tm=syn_mf_Grc_nmda_tm,
        syn_base_name="mf_Grc_nmda",
    )

    run_static_vs_tm_test(
        cell_name=cell_name,
        cell_params=params_grc,
        receptor_type=str(3),
        receptor_dict=receptor_dict,
        syn_static=syn_GoC_GrC_gaba,
        syn_tm=syn_GoC_GrC_gaba_tm,
        syn_base_name="GoC_GrC_gaba",
    )
