from arborize import define_model

GolgiCellModel = define_model(
    {
        "synapse_types": {
            "AMPA_PF": {
                "mechanism": "AMPA",
                "parameters": {
                    "tau_facil": 54,
                    "tau_rec": 35.1,
                    "tau_1": 30,
                    "gmax": 1200,
                    "U": 0.4,
                },
            },
            "AMPA_AA": {
                "mechanism": "AMPA",
                "parameters": {
                    "tau_facil": 54,
                    "tau_rec": 35.1,
                    "tau_1": 30,
                    "gmax": 1200,
                    "U": 0.4,
                },
            },
            "AMPA_MF": {
                "mechanism": ("AMPA", "granule"),
                "parameters": {
                    "tau_facil": 8,
                    "tau_rec": 5,
                    "tau_1": 1,
                    "gmax": 1200,
                    "U": 0.43,
                },
            },
            "NMDA": {
                "mechanism": ("NMDA", "stellate"),
                "parameters": {
                    "tau_facil": 5,
                    "tau_rec": 8,
                    "tau_1": 1,
                    "gmax": 10000,
                    "U": 0.43,
                },
            },
            "GABA": {
                "mechanism": "GABA",
                "parameters": {
                    "tau_facil": 0,
                    "tau_rec": 38.7,
                    "tau_1": 1,
                    "gmaxA1": 3200,
                    "U": 0.42,
                    "Erev": -65,
                },
            },
            "gap_s": {"mechanism": "gap_junction", "source": "vgap"},
            "gap": {"mechanism": ("gap_junction", "parallel"), "source": "vgap"},
        },
        "cable_types": {
            "soma": {
                "cable": {"Ra": 122, "cm": 1},
                "ions": {
                    "na": {"rev_pot": 60},
                    "k": {"rev_pot": -80},
                    "ca": {"rev_pot": 140.70872293},
                },
                "mechanisms": {
                    "Leak": {"e": -55, "gmax": 3e-05},
                    "Nav1_6": {"gbar": 0.14927733727426},
                    "Kv1_1": {"gbar": 0.00549507510519},
                    "Kv3_4": {"gkbar": 0.14910988921938},
                    "Kv4_3": {"gkbar": 0.00406420380423},
                    "Kca1_1": {"gbar": 0.01764345789036},
                    "Kca3_1": {"gkbar": 0.10177335775222},
                    "Ca": {"gcabar": 0.0087689418803},
                    "Cav3_1": {"pcabar": 3.407734319e-05},
                    ("cdp5", "CAM_GoC"): {"TotalPump": 1e-07},
                },
            },
            "dendrites": {"cable": {}, "ions": {}, "mechanisms": {}},
            "basal_dendrites": {
                "cable": {"Ra": 122, "cm": 2.5},
                "ions": {
                    "na": {"rev_pot": 60},
                    "k": {"rev_pot": -80},
                    "ca": {"rev_pot": 140.70872293},
                },
                "mechanisms": {
                    "Leak": {"e": -55, "gmax": 3e-05},
                    "Nav1_6": {"gbar": 0.0080938853146},
                    "Kca1_1": {"gbar": 0.01226052748146},
                    "Kca2_2": {"gkbar": 0.01650689958385},
                    "Ca": {"gcabar": 0.00139885617712},
                    ("cdp5", "CAM_GoC"): {"TotalPump": 2e-09},
                },
                "synapses": {
                    "AMPA_PF": {
                        "mechanism": "AMPA",
                        "parameters": {
                            "tau_facil": 54,
                            "tau_rec": 35.1,
                            "tau_1": 30,
                            "gmax": 1200,
                            "U": 0.4,
                        },
                    },
                    "AMPA_AA": {
                        "mechanism": "AMPA",
                        "parameters": {
                            "tau_facil": 54,
                            "tau_rec": 35.1,
                            "tau_1": 30,
                            "gmax": 1200,
                            "U": 0.4,
                        },
                    },
                    "AMPA_MF": {
                        "mechanism": ("AMPA", "granule"),
                    },
                    "NMDA": {
                        "mechanism": ("NMDA", "stellate"),
                        "parameters": {
                            "tau_facil": 5,
                            "tau_rec": 8,
                            "tau_1": 1,
                            "gmax": 10000,
                            "U": 0.43,
                        },
                    },
                    "GABA": {
                        "mechanism": "GABA",
                        "parameters": {
                            "tau_facil": 0,
                            "tau_rec": 38.7,
                            "tau_1": 1,
                            "gmaxA1": 3200,
                            "U": 0.42,
                            "Erev": -65,
                        },
                    },
                },
            },
            "apical_dendrites": {
                "cable": {"Ra": 122, "cm": 2.5},
                "ions": {
                    "na": {"rev_pot": 60},
                    "k": {"rev_pot": -80},
                    "ca": {"rev_pot": 140.70872293},
                },
                "mechanisms": {
                    "Leak": {"e": -55, "gmax": 3e-05},
                    "Nav1_6": {"gbar": 0.00499506303209},
                    "Kca1_1": {"gbar": 0.01016375552607},
                    "Kca2_2": {"gkbar": 0.00247172479141},
                    "Cav2_3": {"gcabar": 0.00128859564935},
                    "Cav3_1": {"pcabar": 3.690771983e-05},
                    ("cdp5", "CAM_GoC"): {"TotalPump": 5e-09},
                },
                "synapses": {
                    "AMPA_PF": {
                        "mechanism": "AMPA",
                        "parameters": {
                            "tau_facil": 54,
                            "tau_rec": 35.1,
                            "tau_1": 30,
                            "gmax": 1200,
                            "U": 0.4,
                        },
                    },
                },
            },
            "axon": {
                "cable": {"Ra": 122, "cm": 1},
                "ions": {
                    "na": {"rev_pot": 60},
                    "k": {"rev_pot": -80},
                    "ca": {"rev_pot": 140.70872293},
                },
                "mechanisms": {
                    "Leak": {"e": -55, "gmax": 1e-06},
                    "Nav1_6": {"gbar": 0.0115},
                    "Kv3_4": {"gkbar": 0.0091},
                    ("cdp5", "CAM_GoC"): {"TotalPump": 1e-08},
                },
            },
            "axon_initial_segment": {
                "cable": {"Ra": 122, "cm": 1},
                "ions": {
                    "na": {"rev_pot": 60},
                    "k": {"rev_pot": -80},
                    "ca": {"rev_pot": 140.70872293},
                },
                "mechanisms": {
                    "Leak": {"e": -55, "gmax": 3e-05},
                    ("HCN1", "golgi"): {"gbar": 0.0003371456442},
                    "HCN2": {"gbar": 0.00030643090764},
                    "Nav1_6": {"gbar": 0.17233663543619},
                    "Ca": {"gcabar": 0.00595046001148},
                    "Kca1_1": {"gbar": 0.10008178886943},
                    "Km": {"gkbar": 0.00024381226198},
                    ("cdp5", "CAM_GoC"): {"TotalPump": 1e-08},
                },
            },
        },
    },
    use_defaults=True,
)

GolgiCellModel.swc_tags = {
    6: ["dendrites", "basal_dendrites"],
    7: ["dendrites", "apical_dendrites"],
    8: ["axon", "axon_initial_segment"],
}
