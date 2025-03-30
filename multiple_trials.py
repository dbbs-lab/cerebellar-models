from bsb import parse_configuration_file,open_storage, Scaffold
import numpy as np
from bsb_nest.devices import PoissonGenerator
import os
from mpi4py import MPI 

if __name__ == "__main__":
    cfg = parse_configuration_file('./configurations/mouse/dcn-io/dcn_io_awake_nest.yaml', parser="yaml")
    storage = open_storage('mouse_cereb_dcn_io.hdf5')
    trials = np.arange(2,16)
    duration = 1000
    cfg.simulations['mf_cf_stimulus'].duration = duration*(len(trials)+1)
    for i, trial in enumerate(trials):
        cfg.simulations['mf_cf_stimulus'].devices[f'mf_stimulus_{i+1}'] = PoissonGenerator(**{
        "rate": 40,
        "start": 500 + (i+1) * duration,
        "stop": 760 + (i+1) * duration,
        "targetting": {
            "strategy": "cell_model",
            "cell_models": [
                "mossy_fibers"
            ]
        },
        "weight": 1,
        "delay": 0.1
        })
        cfg.simulations['mf_cf_stimulus'].devices[f'cf_stimulus_{i + 1}'] = PoissonGenerator(**{
            "rate": 500,
            "start": 750+ (i+1) * duration,
            "stop": 760+ (i+1) * duration,
            "targetting": {
                "strategy": "cell_model",
                "cell_models": [
                    "io"
                ]
            },
            "receptor_type": 1,
            "weight": 55,
            "delay": 0.1
        })

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    scaffold = Scaffold(cfg, storage)
    results = scaffold.run_simulation('mf_cf_stimulus')
    results.write(f"./simulations_plasticity/plasticity_healthy_{rank}.nio", "ow")


