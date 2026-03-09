import unittest

import numpy as np
from bsb import Configuration, ConfigurationError, Scaffold
from bsb_test import NumpyTestCase, RandomStorageFixture


class TestIoMolecular(
    RandomStorageFixture,
    NumpyTestCase,
    unittest.TestCase,
    engine_name="hdf5",
):
    def setUp(self):
        super().setUp()
        self.radius = 40
        self.convergence = 4
        self.chunk_size = np.array([30, 30, 30])
        self.cfg = Configuration.default(
            network=dict(
                chunk_size=self.chunk_size,
                x=self.chunk_size[0] * 2,
                y=self.chunk_size[1] * 2,
                z=self.chunk_size[2] * 2,
            ),
            cell_types=dict(
                io=dict(spatial=dict(radius=2, count=16)),
                pc=dict(spatial=dict(radius=2, count=24)),
                mli=dict(spatial=dict(radius=2, count=32)),
                mli2=dict(spatial=dict(radius=2, count=16)),
                mli3=dict(spatial=dict(radius=2, count=16)),
            ),
            partitions=dict(
                layer=dict(
                    type="rhomboid",
                    dimensions=[
                        self.chunk_size[0] * 2,
                        self.chunk_size[1] * 2,
                        self.chunk_size[2] * 2,
                    ],
                )
            ),
            placement=dict(
                random_placement=dict(
                    strategy="bsb.placement.RandomPlacement",
                    partitions=["layer"],
                    cell_types=["io", "pc", "mli", "mli2", "mli3"],
                ),
            ),
            connectivity=dict(
                mli_to_pc=dict(
                    strategy="bsb.connectivity.AllToAll",
                    presynaptic=dict(cell_types=["mli", "mli2"]),
                    postsynaptic=dict(cell_types=["pc"]),
                ),
                mli_to_mli2=dict(
                    strategy="bsb.connectivity.AllToAll",
                    presynaptic=dict(cell_types=["mli"]),
                    postsynaptic=dict(cell_types=["mli2"]),
                ),
                io_to_pc=dict(
                    strategy="bsb.connectivity.AllToAll",
                    presynaptic=dict(cell_types=["io"]),
                    postsynaptic=dict(cell_types=["pc"]),
                ),
            ),
        )
        self.network = Scaffold(self.cfg, self.storage)
        self.network.compile(skip_connectivity=True)

    def test_errors_io_mli_strat(self):
        with self.assertRaises(ConfigurationError):
            # wrong mli_to_pc strat
            self.cfg.connectivity.add(
                "io_to_mli",
                dict(
                    strategy="cerebellar_models.connectome.io_molecular.ConnectomeIO_MLI",
                    presynaptic=dict(cell_types=["io"]),
                    postsynaptic=dict(cell_types=["mli", "mli2"]),
                    io_pc_connectivity=["io_to_pc"],
                    mli_pc_connectivity=["mli_to_mli2"],
                    pre_cell_pc="pc",
                ),
            )
        with self.assertRaises(ConfigurationError):
            # not all mli in mli_to_pc strat
            self.cfg.connectivity.add(
                "io_to_mli2",
                dict(
                    strategy="cerebellar_models.connectome.io_molecular.ConnectomeIO_MLI",
                    presynaptic=dict(cell_types=["io"]),
                    postsynaptic=dict(cell_types=["mli", "mli2", "mli3"]),
                    io_pc_connectivity=["io_to_pc"],
                    mli_pc_connectivity=["mli_to_pc"],
                    pre_cell_pc="pc",
                ),
            )
        with self.assertRaises(ConfigurationError):
            # wrong io_to_pc strat
            self.cfg.connectivity.add(
                "io_to_mli3",
                dict(
                    strategy="cerebellar_models.connectome.io_molecular.ConnectomeIO_MLI",
                    presynaptic=dict(cell_types=["io"]),
                    postsynaptic=dict(cell_types=["mli", "mli2"]),
                    io_pc_connectivity=["mli_to_mli2"],
                    mli_pc_connectivity=["mli_to_pc"],
                    pre_cell_pc="pc",
                ),
            )
        with self.assertRaises(ConfigurationError):
            # not all io in io_to_pc strat
            self.cfg.connectivity.add(
                "io_to_mli4",
                dict(
                    strategy="cerebellar_models.connectome.io_molecular.ConnectomeIO_MLI",
                    presynaptic=dict(cell_types=["io", "mli3"]),
                    postsynaptic=dict(cell_types=["mli", "mli2"]),
                    io_pc_connectivity=["io_to_pc"],
                    mli_pc_connectivity=["mli_to_pc"],
                    pre_cell_pc="pc",
                ),
            )

    def test_io_mli(self):
        self.cfg.connectivity.add(
            "io_to_mli",
            dict(
                strategy="cerebellar_models.connectome.io_molecular.ConnectomeIO_MLI",
                presynaptic=dict(cell_types=["io"]),
                postsynaptic=dict(cell_types=["mli", "mli2"]),
                io_pc_connectivity=["io_to_pc"],
                mli_pc_connectivity=["mli_to_pc"],
                pre_cell_pc="pc",
            ),
        )
        self.network.compile(append=True, skip_placement=True)
        all_io = self.network.get_placement_set("io")
        all_pc = self.network.get_placement_set("pc")
        for post_mli in ["mli", "mli2"]:
            all_mli = self.network.get_placement_set(post_mli)
            connected = np.zeros((len(all_io), len(all_mli)), dtype=int)
            cs = self.network.get_connectivity_set(f"io_to_mli_io_to_{post_mli}")
            for from_, to_ in cs.load_connections().as_globals():
                connected[from_[0], to_[0]] += 1
            self.assertAll(connected == len(all_pc))

    def test_regression_68(self):
        self.cfg.connectivity = {}
        self.cfg.after_placement = dict(
            label_cell=dict(
                strategy="cerebellar_models.placement.microzones.LabelCells",
                cell_types=["io", "pc"],
                labels=["p", "m"],
            ),
        )

        self.cfg.connectivity.add(
            "mli_to_pc",
            dict(
                strategy="bsb.connectivity.AllToAll",
                presynaptic=dict(cell_types=["mli", "mli2"]),
                postsynaptic=dict(cell_types=["pc"]),
            ),
        )
        self.cfg.connectivity.add(
            "io_to_pc_p",
            dict(
                strategy="bsb.connectivity.AllToAll",
                presynaptic=dict(cell_types=["io"], labels=["p"]),
                postsynaptic=dict(cell_types=["pc"], labels=["p"]),
            ),
        )
        self.cfg.connectivity.add(
            "io_to_mli_p",
            dict(
                strategy="cerebellar_models.connectome.io_molecular.ConnectomeIO_MLI",
                presynaptic=dict(cell_types=["io"], labels=["p"]),
                postsynaptic=dict(cell_types=["mli", "mli2"]),
                io_pc_connectivity=["io_to_pc_p"],
                mli_pc_connectivity=["mli_to_pc"],
                pre_cell_pc="pc",
            ),
        )
        self.cfg.connectivity.add(
            "io_to_pc_m",
            dict(
                strategy="bsb.connectivity.AllToAll",
                presynaptic=dict(cell_types=["io"], labels=["m"]),
                postsynaptic=dict(cell_types=["pc"], labels=["m"]),
            ),
        )
        self.cfg.connectivity.add(
            "io_to_mli_m",
            dict(
                strategy="cerebellar_models.connectome.io_molecular.ConnectomeIO_MLI",
                presynaptic=dict(cell_types=["io"], labels=["m"]),
                postsynaptic=dict(cell_types=["mli", "mli2"]),
                io_pc_connectivity=["io_to_pc_m"],
                mli_pc_connectivity=["mli_to_pc"],
                pre_cell_pc="pc",
            ),
        )
        self.network.compile(append=True, skip_placement=True)
        ps_io = self.network.get_placement_set("io")
        io_plus = ps_io.get_labelled(["p"])
        io_minus = ps_io.get_labelled(["m"])
        io_not_plus = set(ps_io.load_ids()) - set(io_plus)
        self.assertEqual(io_not_plus, set(io_minus), "minus and plus set should not overlap")
        io_mli_m = self.network.get_connectivity_set("io_to_mli_m_io_to_mli")
        io_mli_p = self.network.get_connectivity_set("io_to_mli_p_io_to_mli")

        io_m = io_mli_m.load_connections().all()[0][:, 0]
        io_p = io_mli_p.load_connections().all()[0][:, 0]
        self.assertAll(np.unique(io_m) == io_minus, "All io minus cells should be connected")
        self.assertAll(np.unique(io_p) == io_plus, "All io minus cells should be connected")
