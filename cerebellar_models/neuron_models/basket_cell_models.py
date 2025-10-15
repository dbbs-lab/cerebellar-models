from arborize import define_model

BasketCellModel = define_model(
    {
        "synapse_types": {
            "AMPA": {
                "mechanism": "AMPA",
                "parameters": {
                    "tau_facil": 10.8,
                    "tau_rec": 35.1,
                    "tau_1": 10,
                    "gmax": 1200,
                    "U": 0.13,
                },
            },
            "NMDA": {
                "mechanism": ("NMDA", "stellate"),
                "parameters": {
                    "tau_facil": 10.8,
                    "tau_rec": 35.1,
                    "tau_1": 10,
                    "gmax": 5000,
                    "U": 0.13,
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
        "cable_types": {
            "soma": {
                "cable": {"Ra": 122, "cm": 1},
                "ions": {
                    "na": {"rev_pot": 60},
                    "k": {"rev_pot": -80},
                    "ca": {"rev_pot": 137.5},
                    "h": {"rev_pot": -34},
                },
                "mechanisms": {
                    "Leak": {"e": -55, "gmax": 0.00004},
                    "Nav1_1": {"gbar": 0.2},
                    "Cav3_2": {"gcabar": 0.0001},
                    "Cav1_2": {"gbar": 0.0007},
                    "Cav1_3": {"gbar": 5e-06},
                    "Kir2_3": {"gkbar": 0.0001},
                    "Kv3_4": {"gkbar": 0.097},
                    "Kv4_3": {"gkbar": 0.01},
                    "Kca3_1": {"gkbar": 0.001},
                    "HCN1": {"gbar": 0.001},
                    "cdp5": {"TotalPump": 2e-09},
                },
            },
            "dendrites": {
                "cable": {"Ra": 122, "cm": 1},
                "ions": {"k": {"rev_pot": -80}},
                "mechanisms": {
                    "Leak": {"e": -55, "gmax": 1e-05},
                    "Cav3_2": {"gcabar": 0.00005},
                    "Cav1_2": {"gbar": 0.0002},
                    "Cav1_3": {"gbar": 5e-06},
                    "Kv4_3": {"gkbar": 0.00987201764943},
                    "Kca2_2": {"gkbar": 0.0065},
                    "cdp5": {"TotalPump": 2e-09},
                },
            },
            "axon": {
                "cable": {"Ra": 122, "cm": 1},
                "ions": {
                    "na": {"rev_pot": 60},
                    "k": {"rev_pot": -80},
                    "h": {"rev_pot": -34},
                },
                "mechanisms": {
                    "Leak": {"e": -55, "gmax": 1e-06},
                    "Kv1_1": {"gbar": 0.0005},
                    "Nav1_6": {"gbar": 0.001},
                    "Kv3_4": {"gkbar": 0.001},
                    "HCN1": {"gbar": 0.0001},
                    "Kca1_1": {"gbar": 0.001},
                    "Cav2_1": {"pcabar": 0.00008},
                    "cdp5": {"TotalPump": 1e-09},
                },
            },
            "axon_initial_segment": {
                "cable": {"Ra": 122, "cm": 1},
                "ions": {
                    "na": {"rev_pot": 60},
                    "k": {"rev_pot": -80},
                    "h": {"rev_pot": -34},
                },
                "mechanisms": {
                    "Leak": {"e": -55, "gmax": 1e-05},
                    "HCN1": {"gbar": 0.001},
                    "Nav1_6": {"gbar": 0.25},
                    "Kv1_1": {"gbar": 0.02333516534870},
                    "Kv3_4": {"gkbar": 0.002},
                    "Kca1_1": {"gbar": 0.01},
                    "Cav2_1": {"pcabar": 2.2e-4},
                    "cdp5": {"TotalPump": 1e-09},
                },
            },
        },
    },
    use_defaults=True,
)

BasketCellModel.swc_tags = {
    6: ["axon", "axon_initial_segment"],
}
