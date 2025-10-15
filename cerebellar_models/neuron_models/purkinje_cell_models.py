from arborize import define_model

PurkinjeCellModel = define_model(
    {
        "synapse_types": {
            "AMPA": {
                "mechanism": "AMPA",
                "parameters": {
                    "tau_facil": 54,
                    "tau_rec": 35.1,
                    "tau_1": 6,
                    "gmax": 1200,
                    "U": 0.13,
                },
            },
            "GABA": {
                "mechanism": "GABA",
                "parameters": {
                    "tau_facil": 4,
                    "tau_rec": 15,
                    "tau_1": 1,
                    "Erev": -70,
                    "gmaxA1": 2600,
                    "U": 0.35,
                },
            },
            "Iext": {"mechanism": "Iext", "parameters": {"amp": 0.00}},
        },
        "cable_types": {
            "soma": {
                "cable": {"Ra": 122, "cm": 1},
                "ions": {
                    "na": {"rev_pot": 60},
                    "k": {"rev_pot": -88},
                    "h": {"rev_pot": -34.4},
                    "ca": {"rev_pot": 137.52625},
                },
                "mechanisms": {
                    "Leak": {"e": -61, "gmax": 0.001},
                    "Nav1_6": {"gbar": 0.19835422387652},
                    "Kv1_1": {"gbar": 0.00185494695569},
                    "Kv3_4": {"gkbar": 0.05684366437809},
                    "Kir2_3": {"gkbar": 2.100544604e-05},
                    "Cav2_1": {"pcabar": 0.0001530123689},
                    "Cav3_1": {"pcabar": 4.19659246e-06},
                    "Cav3_2": {"gcabar": 0.00067268089296},
                    "Cav3_3": {"pcabar": 0.00015401124325},
                    "Kca1_1": {"gbar": 0.00877088561313},
                    "Kca2_2": {"gkbar": 0.00085748826176},
                    "Kca3_1": {"gkbar": 0.00978230722317},
                    "HCN1": {"gbar": 0.00192222696826},
                    ("cdp5", "CAM"): {"TotalPump": 2e-08},
                },
                "synapses": {
                    "GABA": {
                        "mechanism": "GABA",
                        "parameters": {
                            "tau_facil": 4,
                            "tau_rec": 15,
                            "tau_1": 1,
                            "Erev": -70,
                            "gmaxA1": 2600,
                            "U": 0.35,
                        },
                    },
                },
            },
            "dendrites": {
                "cable": {
                    # "cm": lambda d: (11.510294 * math.exp(-1.376463 * d) + 2.120503),
                    # Formula for cm varied from 2 to 6, aprroximated with avg
                    "cm": 5.577701084442752,
                    "Ra": 122,
                },
                "ions": {
                    "na": {"rev_pot": 60},
                    "k": {"rev_pot": -88},
                    "h": {"rev_pot": -34.4},
                    "ca": {"rev_pot": 137.52625},
                },
                "mechanisms": {
                    "Leak": {"e": -61, "gmax": 0.0003},
                    "Kv1_1": {"gbar": 0.00105296938319},
                    "Kv1_5": {"gKur": 0.00018237109331},
                    "Kv3_3": {"gbar": 0.01192448545734},
                    "Kv4_3": {"gkbar": 0.00082020749325},
                    "Cav2_1": {"pcabar": 0.001594103024736},
                    "Cav3_3": {"pcabar": 0.00013000262047},
                    "Kca1_1": {"gbar": 0.03803999223084},
                    "HCN1": {"gbar": 2.49811225e-06},
                    ("cdp5", "CAM"): {"TotalPump": 5e-08},
                },
                "synapses": {
                    "AMPA": {
                        "mechanism": "AMPA",
                        "parameters": {
                            "tau_facil": 54,
                            "tau_rec": 35.1,
                            "tau_1": 6,
                            "gmax": 1200,
                            "U": 0.13,
                        },
                    },
                },
            },
            "basal_dendrites": {
                "cable": {"Ra": 122, "cm": 1},
                "ions": {
                    "na": {"rev_pot": 60},
                    "k": {"rev_pot": -88},
                    "ca": {"rev_pot": 137.52625},
                },
                "mechanisms": {
                    "Leak": {"e": -61, "gmax": 0.0003},
                    "Kir2_3": {"gkbar": 1.135399111e-05},
                    "Cav3_1": {"pcabar": 4.0066819e-06},
                    "Cav3_2": {"gcabar": 0.00194175819441},
                    "Kca2_2": {"gkbar": 0.00078534639818},
                    "Kca3_1": {"gkbar": 0.00330481139341},
                },
            },
            "sodium_dendrites": {
                "cable": {"Ra": 122, "cm": 1},
                "ions": {
                    "na": {"rev_pot": 60},
                    "k": {"rev_pot": -88},
                    "ca": {"rev_pot": 137.52625},
                },
                "mechanisms": {
                    "Leak": {"e": -61, "gmax": 0.0003},
                    "Kir2_3": {"gkbar": 1.135399111e-05},
                    "Cav3_1": {"pcabar": 4.0066819e-06},
                    "Cav3_2": {"gcabar": 0.00194175819441},
                    "Kca2_2": {"gkbar": 0.00078534639818},
                    "Kca3_1": {"gkbar": 0.00330481139341},
                    "Nav1_6": {"gbar": 0.01568012827236},
                },
            },
            "aa_targets": {"cable": {}, "ions": {}, "mechanisms": {}},
            "pf_targets": {"cable": {}, "ions": {}, "mechanisms": {}},
            "sc_targets": {"cable": {}, "ions": {}, "mechanisms": {}},
            "AIS": {
                "cable": {"Ra": 122, "cm": 1},
                "ions": {
                    "na": {"rev_pot": 60},
                    "k": {"rev_pot": -88},
                    "ca": {"rev_pot": 137.52625},
                },
                "mechanisms": {
                    "Leak": {"e": -61, "gmax": 0.0003},
                    "Nav1_6": {"gbar": 0.76520432560779},
                    "Kv3_4": {"gkbar": 0.01470016164534},
                    "Cav2_1": {"pcabar": 0.00028786836482},
                    "Cav3_1": {"pcabar": 6.43429659e-06},
                    ("cdp5", "CAM"): {"TotalPump": 2e-08},
                },
            },
            "axon": {"cable": {}, "ions": {}, "mechanisms": {}},
            "AIS_K": {
                "cable": {"Ra": 122, "cm": 1},
                "ions": {"k": {"rev_pot": -88}},
                "mechanisms": {
                    "Leak": {"e": -61, "gmax": 0.0003},
                    "Kv1_1": {"gbar": 0.01011038066089},
                },
            },
            "axonmyelin": {
                "cable": {"Ra": 122, "cm": 1.87e-11},
                "ions": {},
                "mechanisms": {"pas": {"e": -61, "g": 5.6e-09}},
            },
            "nodes": {
                "cable": {"Ra": 122, "cm": 1},
                "ions": {
                    "na": {"rev_pot": 60},
                    "k": {"rev_pot": -88},
                    "ca": {"rev_pot": 137.52625},
                },
                "mechanisms": {
                    "Leak": {"e": -61, "gmax": 0.0003},
                    "Nav1_6": {"gbar": 0.02749309338198},
                    "Kv3_4": {"gkbar": 0.02227585026371},
                    "Cav2_1": {"pcabar": 0.00012353396215},
                    "Cav3_1": {"pcabar": 1.508814156e-05},
                    ("cdp5", "CAM"): {"TotalPump": 5e-07},
                },
            },
        },
    },
    use_defaults=True,
)

PurkinjeCellModel.swc_tags = {
    6: ["axon", "AIS"],
    7: ["axon", "AIS_K"],
    8: ["axon", "axonmyelin"],
    9: ["axon", "nodes"],
    10: ["dendrites", "basal_dendrites"],
    11: ["dendrites", "pf_targets", "sc_targets"],
    12: ["dendrites", "aa_targets", "sc_targets"],
    13: ["dendrites", "sodium_dendrites"],
}
