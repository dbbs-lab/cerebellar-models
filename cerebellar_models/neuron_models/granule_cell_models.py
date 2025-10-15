from arborize import define_model

GranuleCellModel = define_model(
    {
        "synapse_types": {
            "AMPA": {
                "mechanism": ("AMPA", "granule"),
                "parameters": {
                    "tau_facil": 5,
                    "tau_rec": 8,
                    "tau_1": 1,
                    "gmax": 1200,
                    "U": 0.43,
                },
            },
            "NMDA": {
                "mechanism": ("NMDA", "granule"),
                "parameters": {
                    "tau_facil": 5,
                    "tau_rec": 8,
                    "tau_1": 1,
                    "gmax": 18800,
                    "U": 0.43,
                },
            },
            "GABA": {"mechanism": ("GABA", "granule"), "parameters": {"U": 0.35}},
        },
        "cable_types": {
            "soma": {
                "cable": {"Ra": 100, "cm": 2},
                "ions": {"k": {"rev_pot": -88}, "ca": {"rev_pot": 137.5}},
                "mechanisms": {
                    "Leak": {"e": -60, "gmax": 0.00029038073716},
                    "Kv3_4": {"gkbar": 0.00076192450952},
                    "Kv4_3": {"gkbar": 0.00281496839061},
                    "Kir2_3": {"gkbar": 0.00074725514702},
                    "Ca": {"gcabar": 0.00060938071784},
                    "Kv1_1": {"gbar": 0.00569738264555},
                    "Kv1_5": {"gKur": 0.00083407556714},
                    "Kv2_2": {"gKv2_2bar": 1.203410852e-05},
                    ("cdp5", "CR"): {},
                },
            },
            "dendrites": {
                "cable": {"Ra": 100, "cm": 2.5},
                "ions": {"k": {"rev_pot": -88}, "ca": {"rev_pot": 137.5}},
                "mechanisms": {
                    "Leak": {"e": -60, "gmax": 0.00025029700737},
                    ("Leak", "GABA"): {},
                    "Ca": {"gcabar": 0.00500128008459},
                    "Kca1_1": {"gbar": 0.01001807454651},
                    "Kv1_1": {"gbar": 0.00381819207934},
                    ("cdp5", "CR"): {},
                },
            },
            "axon": {"cable": {}, "ions": {}, "mechanisms": {}},
            "ascending_axon": {
                "cable": {"Ra": 100, "cm": 1},
                "ions": {
                    "na": {"rev_pot": 87.39},
                    "k": {"rev_pot": -88},
                    "ca": {"rev_pot": 137.5},
                },
                "mechanisms": {
                    ("Na", "granule_cell"): {"gnabar": 0.02630163681502},
                    "Kv3_4": {"gkbar": 0.00237386061632},
                    "Leak": {"e": -60, "gmax": 9.364092125e-05},
                    "Ca": {"gcabar": 0.00068197420273},
                    ("cdp5", "CR"): {},
                },
            },
            "axon_initial_segment": {
                "cable": {"Ra": 100, "cm": 1},
                "ions": {
                    "na": {"rev_pot": 87.39},
                    "k": {"rev_pot": -88},
                    "ca": {"rev_pot": 137.5},
                },
                "mechanisms": {
                    ("Na", "granule_cell_FHF"): {"gnabar": 1.28725006737226},
                    "Kv3_4": {"gkbar": 0.00649595340654},
                    "Leak": {"e": -60, "gmax": 0.00029276697557},
                    "Ca": {"gcabar": 0.00031198539472},
                    "Km": {"gkbar": 0.00056671971737},
                    ("cdp5", "CR"): {},
                },
            },
            "axon_hillock": {
                "cable": {"Ra": 100, "cm": 2},
                "ions": {
                    "na": {"rev_pot": 87.39},
                    "k": {"rev_pot": -88},
                    "ca": {"rev_pot": 137.5},
                },
                "mechanisms": {
                    "Leak": {"e": -60, "gmax": 0.0003695818972},
                    ("Na", "granule_cell_FHF"): {"gnabar": 0.00928805851462},
                    "Kv3_4": {"gkbar": 0.02037346310915},
                    "Ca": {"gcabar": 0.00057726155447},
                    ("cdp5", "CR"): {},
                },
            },
            "parallel_fiber": {
                "cable": {"Ra": 100, "cm": 1},
                "ions": {
                    "na": {"rev_pot": 87.39},
                    "k": {"rev_pot": -88},
                    "ca": {"rev_pot": 137.5},
                },
                "mechanisms": {
                    ("Na", "granule_cell"): {"gnabar": 0.01771848449261},
                    "Kv3_4": {"gkbar": 0.00817568047037},
                    "Leak": {"e": -60, "gmax": 3.5301616e-07},
                    "Ca": {"gcabar": 0.0002085683353},
                    ("cdp5", "CR"): {},
                },
            },
        },
    },
    use_defaults=True,
)

GranuleCellModel.swc_tags = {
    6: ["axon", "axon_hillock"],
    7: ["axon", "axon_initial_segment"],
    8: ["axon", "ascending_axon"],
    9: ["axon", "parallel_fiber"],
}
