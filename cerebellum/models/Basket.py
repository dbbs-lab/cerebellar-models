from arborize import define_model

definitionBasket = define_model({
    "synapse_types" : {
        "AMPA": {
            "mechanism": "AMPA",
            "attributes": {
                "tau_facil": 54,
                "tau_rec": 35.1,
                "tau_1": 6,
                "gmax": 12000,
                "U": 0.4,
            },
        },
        "NMDA": {
            "mechanism": ("NMDA", "stellate"),
            "attributes": {
                "tau_facil": 5,
                "tau_rec": 8,
                "tau_1": 1,
                "gmax": 50000,
                "U": 0.15,
            },
        },
        "GABA": {
            "mechanism": "GABA",
            "attributes": {
                "tau_facil": 0,
                "tau_rec": 38.7,
                "tau_1": 1,
                "gmax": 32000,
                "U": 0.42,
                "Erev": -65,
                },
            },
        },
    "cable_types" : {
        "soma": {
            "cable": {"Ra": 122, "cm": 1},
            "ions": {
                "na": {"e": 60},
                "k": {"e": -80},
                "ca": {"e": 137.5},
                "h": {"e": -34},
            },
            "mechanisms": {
                "Leak": {"e": -60, "gmax": 3e-05},
                "Nav1_1": {"gbar": 0.10946415489712},
                "Cav3_2": {"gcabar": 0.0006295539},
                "Cav3_3": {"pcabar": 0.0007543986},
                # "CaL13": {"pcaLbar": 0.0002593795},
                "Kir2_3": {"gkbar": 0.0012770833},
                "Kv3_4": {"gkbar": 0.0277114781},
                "Kca1_1": {"gbar": 0.0044251081},
                "Cav2_1": {"pcabar": 0.000846789},
                "HCN1": {"gbar": 0.0006902196},
                "cdp5": {"TotalPump": 1e-09},
            },
        },
        "dendrites": {
            "cable": {"Ra": 122, "cm": 1},
            "ions": {"k": {"e": -80}},
            "mechanisms": {
                "Leak": {"e": -60, "gmax": 3e-05},
                "Cav2_1": {"pcabar": 0.0004965596},
                "Kca1_1": {"gbar": 0.0020575902},
                "Kv1_1": {"gbar": 0.0285137286},
                "cdp5": {"TotalPump": 1e-09},
            },
            "synapses": {
            "AMPA": {
            "mechanism": "AMPA",
            "attributes": {
                "tau_facil": 54,
                "tau_rec": 35.1,
                "tau_1": 6,
                "gmax": 1200,
                "U": 0.4,
                },
            },
            "NMDA": {
                "mechanism": ("NMDA", "stellate"),
                "attributes": {
                    "tau_facil": 5,
                    "tau_rec": 8,
                    "tau_1": 1,
                    "gmax": 5000,
                    "U": 0.15,
                },
            },
            "GABA": {
                "mechanism": "GABA",
                "attributes": {
                    "tau_facil": 0,
                    "tau_rec": 38.7,
                    "tau_1": 1,
                    "gmax": 3200,
                    "U": 0.42,
                    "Erev": -65,
                },
            },
            },
        },
        "axon": {
            "cable": {"Ra": 122, "cm": 1},
            "ions": {"na": {"e": 60}, "k": {"e": -80}, "h": {"e": -34}},
            "mechanisms": {
                "Leak": {"e": -60, "gmax": 3e-05},
                "Kv1_1": {"gbar": 0.0069654709},
                "Nav1_6": {"gbar": 0.0072983226},
                "Kv3_4": {"gkbar": 0.0151487764},
                "HCN1": {"gbar": 0.0034633208},
                "cdp5": {},
            },
        },
        "axon_initial_segment": {
            "cable": {"Ra": 122, "cm": 1},
            "ions": {"na": {"e": 60}, "k": {"e": -80}, "h": {"e": -34}},
            "mechanisms": {
                "Leak": {"e": -60, "gmax": 3e-05},
                "HCN1": {"gbar": 0.0048096086},
                "Nav1_6": {"gbar": 0.5724695612},
                "Kv1_1": {"gbar": 0.0827297077},
                "Kv3_4": {"gkbar": 0.0300388404},
                "cdp5": {},
            },
            },
    },
})

tagsBk = {
    16: ["axon", "axon_initial_segment"],
}
