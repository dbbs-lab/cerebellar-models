from utils import *

receptor_dict = {
    "1": ["AMPA", "I_syn_ampa", "AMPA_E_rev"],
    "2": ["NMDA", "I_syn_nmda", "NMDA_E_rev"],
    "3": ["GABA", "I_syn_gaba", "GABA_E_rev"],
    "4": ["AMPA2", "I_syn_ampa2", "AMPA2_E_rev"],
    "5": ["AMPA3", "I_syn_ampa3", "AMPA3_E_rev"],
}

if __name__ == "__main__":
    cell_name = "BC"
    params_BC = {
        "t_ref": 1.59,
        "C_m": 14.6,
        "V_th": -53.0,
        "V_reset": -78.0,
        "E_L": -68.0,
        "V_m": -68.0,
        "tau_m": 9.125,
        "I_e": 5.043407259038052,
        "k_adap": 0.17585099665700937,
        "k_1": 0.1509723879307489,
        "A1": 213.73301121613932,
        "A2": 151.19435263704997,
        "k_2": 0.11428245921171565,
        # pf-BC AMPA
        "AMPA_E_rev": 0.0,
        "AMPA_g_init": 0.5986151464078818,
        "AMPA_Tau_r": 0.2885454465778898,
        "AMPA_Tau_d1": 0.361214360260991,
        "AMPA_Tau_d2": 1.1702221289816273,
        "AMPA_A_r": 2.537866943371753,
        "AMPA_A1": 0.7686340271346837,
        "AMPA_A2": 1.3031463264578722,
        # pf-BC NMDA
        "NMDA_E_rev": -3.7,
        "NMDA_g_init": 0.005421817538203816,
        "NMDA_Tau_r": 47.40777695927046,
        "NMDA_Tau_d1": 51.58825950692386,
        "NMDA_Tau_d2": 769.3667455068031,
        "NMDA_A_r": 6.155174515299148,
        "NMDA_A1": 3.7422751919354367,
        "NMDA_A2": 0.4980434420823697,
        # BC-BC GABA
        "GABA_E_rev": -65.0,
        "GABA_g_init": 1.1487687456539493,
        "GABA_Tau_r": 0.37383077569256745,
        "GABA_Tau_d1": 3.1073175264777047,
        "GABA_Tau_d2": 29.20898014832495,
        "GABA_A_r": 6.80395085901533,
        "GABA_A1": 9.254883332147838,
        "GABA_A2": 6.738700331790382,
    }

    # Static synapse for BC
    pf_BC_ampa = {"weight": 1.0, "delay": 0.1, "receptor_type": 1}
    pf_BC_nmda = {"weight": 1.0, "delay": 0.1, "receptor_type": 2}
    BC_BC_gaba = {"weight": 1.0, "delay": 0.1, "receptor_type": 3}

    # TM parameters for BC synapses
    pf_BC_ampa_tm = {
        "u": 0.0,
        "y": 0.0,
        "x": 1.0,
        "U": 0.15,
        "tau_rec": 35.1,
        "tau_fac": 10.8,
        "tau_psc": 1,
        "weight": 1.0,
        "delay": 0.1,
        "receptor_type": 1,
    }
    pf_BC_nmda_tm = {
        "u": 0.0,
        "y": 0.0,
        "x": 1.0,
        "U": 0.15,
        "tau_rec": 35.1,
        "tau_fac": 10.8,
        "tau_psc": 1,
        "weight": 1.0,
        "delay": 0.1,
        "receptor_type": 2,
    }
    BC_BC_gaba_tm = {
        "u": 0,
        "y": 0,
        "x": 1.0,
        "U": 0.42,
        "tau_rec": 38.7,
        "tau_fac": 4,
        "tau_psc": 1,
        "weight": 1.0,
        "delay": 0.1,
        "receptor_type": 3,
    }

    run_static_vs_tm_test(
        cell_name=cell_name,
        cell_params=params_BC,
        receptor_type=str(1),
        receptor_dict=receptor_dict,
        syn_static=pf_BC_ampa,
        syn_tm=pf_BC_ampa_tm,
        syn_base_name="pf_BC_ampa",
    )

    run_static_vs_tm_test(
        cell_name=cell_name,
        cell_params=params_BC,
        receptor_type=str(2),
        receptor_dict=receptor_dict,
        syn_static=pf_BC_nmda,
        syn_tm=pf_BC_nmda_tm,
        syn_base_name="pf_BC_nmda",
    )

    run_static_vs_tm_test(
        cell_name=cell_name,
        cell_params=params_BC,
        receptor_type=str(3),
        receptor_dict=receptor_dict,
        syn_static=BC_BC_gaba,
        syn_tm=BC_BC_gaba_tm,
        syn_base_name="BC_BC_gaba",
    )
