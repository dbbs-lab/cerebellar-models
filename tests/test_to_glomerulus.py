"""
Unit tests of the `ConnectomeMossyGlomerulus`, `ConnectomeGlomerulus` strategies
"""

import unittest

import numpy as np
from bsb import Scaffold
from bsb_test import (
    FixedPosConfigFixture,
    NetworkFixture,
    NumpyTestCase,
    RandomStorageFixture,
)


def _test_fallback_to_closest_presyn(self, conn_name, strategy, extra_params):
    """
    Shared helper: place one glomerulus (test_cell) and two pre cells whose positions are
    ALL outside the geometric ROI.  After compiling, every glomerulus must be connected and
    the connected pre cell must be the closest one (fallback branch).
    Each pre cell also has a glomerulus at its own position, which it reserves, so that the
    glomerulus under test is not reserved.

    The chunk size is pinned to 100 so all three positions sit in chunk (0,0,0) and are
    therefore loaded together during connect(), ensuring presyn_pos is never empty.
    """
    chunk_size = 100.0
    glom_pos = np.array([[50.0, 50.0, 50.0]])
    # Both pre cells are at x-distance 40 and 45 from the glomerulus.
    # With the small ROIs used by the callers (x_length=5 or radius=5) they are
    # intentionally outside the selection region, triggering the fallback.
    pre_close = np.array([90.0, 50.0, 50.0])
    pre_far = np.array([95.0, 50.0, 50.0])

    # Pin chunk size and override the single glomerulus position.
    self.cfg.network.chunk_size = chunk_size
    self.cfg.cell_types["test_cell"].spatial.count = 3
    self.cfg.placement.ch4_c25.positions = np.vstack([glom_pos, pre_close, pre_far])

    # Add a separate pre cell type with fixed positions.
    self.cfg.cell_types.add("pre_cell", dict(spatial=dict(radius=2.5, count=2)))
    self.cfg.placement.add(
        "place_pre",
        dict(
            strategy="bsb.placement.strategy.FixedPositions",
            partitions=[],
            cell_types=["pre_cell"],
        ),
    )
    self.cfg.placement["place_pre"].positions = np.vstack([pre_close, pre_far])

    # Add the connectivity rule under test.
    self.cfg.connectivity.add(
        conn_name,
        dict(
            strategy=strategy,
            presynaptic=dict(cell_types=["pre_cell"]),
            postsynaptic=dict(cell_types=["test_cell"]),
            **extra_params,
        ),
    )
    self.network = Scaffold(self.cfg, self.storage)
    self.network.compile(clear=True)

    cs = self.network.get_connectivity_set(conn_name)
    pre_positions = self.network.get_placement_set("pre_cell").load_positions()
    glom_positions = self.network.get_placement_set("test_cell").load_positions()

    # Every glomerulus must receive a connection via the fallback.
    self.assertEqual(len(cs), 3, "Glomerulus must be connected even when no pre cell is in the ROI")

    for from_, to_ in cs.load_connections().as_globals():
        if not np.all(glom_positions[to_[0]] == glom_pos[0]):
            continue
        connected_pre_pos = pre_positions[from_[0]]
        # The fallback must select the geometrically closest pre cell.
        expected_closest = pre_positions[
            np.argmin(np.linalg.norm(pre_positions - glom_pos[0], axis=1))
        ]
        self.assertClose(
            connected_pre_pos,
            expected_closest,
            "Closest presynaptic cell must be selected as the fallback target",
        )


def _add_pre_cells(self, positions, conn_name="x_to_glomerulus"):
    """
    Add a presynaptic cell type with fixed positions, used as presynaptic cells of the
    ``conn_name`` connectivity rule.
    """
    self.cfg.cell_types.add("pre_cell", dict(spatial=dict(radius=2.5, count=len(positions))))
    self.cfg.placement.add(
        "place_pre",
        dict(
            strategy="bsb.placement.strategy.FixedPositions",
            partitions=[],
            cell_types=["pre_cell"],
        ),
    )
    self.cfg.placement["place_pre"].positions = positions
    self.cfg.connectivity[conn_name].presynaptic.cell_types = ["pre_cell"]


def _test_distance_to_glomerulus(self, nb_trials=50):
    # Override positions
    self.chunk_size = 100.0
    pos_1 = np.array([0.5, 0.5, 0.5]) * self.chunk_size
    pos_2 = np.array([0.9, 0.5, 0.5]) * self.chunk_size  # close enough
    pos_3 = np.array([0.5, 1.0, 0.5]) * self.chunk_size  # too far away
    # Glomerulus under test, close to pos_1. Each presynaptic cell has a glomerulus at its
    # own position, which it reserves, so the glomerulus under test is drawn.
    glom_pos = pos_1 + np.array([0.0, 0.0, 1.0])
    self.cfg.network.chunk_size = self.chunk_size
    self.cfg.cell_types["test_cell"].spatial.count = 4
    self.cfg.placement.ch4_c25.positions = np.vstack((glom_pos, pos_1, pos_2, pos_3))
    _add_pre_cells(self, np.vstack((pos_1, pos_2, pos_3)))

    sources = np.full(nb_trials, -1)
    for i in range(nb_trials):
        # Connectivity now draws from the network's own configured randomness, keyed on
        # the chunk pair rather than the global generator, so repeating the same network
        # repeats the same draw. Give each trial its own root seed to sample independently.
        self.cfg.rng = {"seed": i}
        self.network = Scaffold(self.cfg, self.storage)
        self.network.compile(clear=True)
        cs = self.network.get_connectivity_set("x_to_glomerulus")
        glom_positions = self.network.get_placement_set("test_cell").load_positions()
        pre_positions = self.network.get_placement_set("pre_cell").load_positions()
        for from_, to_ in cs.load_connections().as_globals():
            if np.all(glom_positions[to_[0]] == glom_pos):
                if np.all(pre_positions[from_[0]] == pos_2):
                    sources[i] = 1
                elif np.all(pre_positions[from_[0]] == pos_1):
                    sources[i] = 0
                else:
                    sources[i] = 2
                break
    self.assertClose(
        np.unique(sources),
        np.array([0, 1]),
        "pos_3 should be unreachable.\npos_2 is less likely but should still happen",
    )
    self.assertTrue(
        np.count_nonzero(sources == 0) > np.count_nonzero(sources == 1),
        "Close targets should be more likely",
    )


def _test_every_presyn_connected(self, extra_params):
    """
    Shared helper: a presynaptic cell far from every glomerulus, which is never drawn by any
    glomerulus, must still be connected to its closest free glomerulus. Glomeruli are spread
    over two chunks, so that each connectivity job only sees part of them.
    """
    self.cfg.network.chunk_size = 100.0
    glom_pos = np.array(
        [[10.0, 50.0, 50.0], [20.0, 50.0, 50.0], [110.0, 50.0, 50.0], [120.0, 50.0, 50.0]]
    )
    self.cfg.cell_types["test_cell"].spatial.count = len(glom_pos)
    self.cfg.placement.ch4_c25.positions = glom_pos
    # The last presynaptic cell is outside of every glomerulus' region of interest.
    pre_pos = np.array([[15.0, 50.0, 50.0], [115.0, 50.0, 50.0], [190.0, 50.0, 50.0]])
    _add_pre_cells(self, pre_pos)
    for key, value in extra_params.items():
        setattr(self.cfg.connectivity["x_to_glomerulus"], key, value)
    self.network = Scaffold(self.cfg, self.storage)
    self.network.compile(clear=True)

    cs = self.network.get_connectivity_set("x_to_glomerulus")
    glom_positions = self.network.get_placement_set("test_cell").load_positions()
    pre_positions = self.network.get_placement_set("pre_cell").load_positions()
    pre_locs, post_locs = cs.load_connections().as_globals().all()
    self.assertAll(
        np.bincount(post_locs[:, 0], minlength=len(glom_pos)) == 1,
        "Each glomerulus has exactly one presynaptic cell",
    )
    self.assertAll(
        np.bincount(pre_locs[:, 0], minlength=len(pre_pos)) >= 1,
        "Each presynaptic cell has at least one glomerulus",
    )
    far_pre = np.nonzero(np.all(pre_positions == pre_pos[2], axis=1))[0][0]
    self.assertClose(
        glom_positions[post_locs[pre_locs[:, 0] == far_pre, 0]],
        glom_pos[[3]],
        "The far presynaptic cell is connected to its closest free glomerulus",
    )


class TestReserveClosestTargets(NumpyTestCase, unittest.TestCase):
    def test_reserve(self):
        from cerebellar_models.connectome.to_glomerulus import reserve_closest_targets

        post = np.array([[0.0, 0, 0], [1.0, 0, 0], [10.0, 0, 0]])
        # Both cells are closest to post 0: the closest one wins, the other one gets the
        # next closest free target.
        pre = np.array([[0.5, 0, 0], [0.1, 0, 0]])
        self.assertAll(reserve_closest_targets(pre, post) == [1, 0])
        # Ties are broken by presynaptic index.
        pre = np.array([[0.0, 0, 0], [0.0, 0, 0]])
        self.assertAll(reserve_closest_targets(pre, post) == [0, 1])
        # More presynaptic cells than targets: the farthest one is left out.
        pre = np.array([[0.0, 0, 0], [1.0, 0, 0], [10.0, 0, 0], [50.0, 0, 0]])
        self.assertAll(reserve_closest_targets(pre, post) == [0, 1, 2, -1])
        # The search extends beyond the first k=8 neighbours.
        post = np.arange(20.0)[:, None] * [1.0, 0, 0]
        pre = np.zeros((12, 3))
        self.assertAll(np.sort(reserve_closest_targets(pre, post)) == np.arange(12))
        self.assertAll(reserve_closest_targets(np.zeros((0, 3)), post) == [])
        self.assertAll(reserve_closest_targets(pre, np.zeros((0, 3))) == -1)


class TestConnectomeMossyGlomerulus(
    RandomStorageFixture,
    NetworkFixture,
    FixedPosConfigFixture,
    NumpyTestCase,
    unittest.TestCase,
    engine_name="hdf5",
):
    def setUp(self):
        super().setUp()
        self.x_length = 60
        self.y_length = 20
        self.cfg.connectivity.add(
            "x_to_glomerulus",
            dict(
                strategy="cerebellar_models.connectome.to_glomerulus.ConnectomeMossyGlomerulus",
                presynaptic=dict(cell_types=["test_cell"]),
                postsynaptic=dict(cell_types=["test_cell"]),
                x_length=self.x_length,
                y_length=self.y_length,
            ),
        )
        self.network = Scaffold(self.cfg, self.storage)

    def test_connect(self):
        self.network.compile()
        cs = self.network.get_connectivity_set("x_to_glomerulus")
        cell_positions = self.network.get_placement_set("test_cell").load_positions()

        cell_targets = np.array([-1, -1])
        distances = np.full(len(cell_positions), -1.0)
        self.assertEqual(len(cs), len(cell_positions), "As many connection as postsyn cell")
        for from_, to_ in cs.load_connections().as_globals():
            self.assertAll(from_[1:] == cell_targets)
            self.assertAll(to_[1:] == cell_targets)
            diff = np.absolute(cell_positions[from_[0]] - cell_positions[to_[0]])
            self.assertTrue(diff[0] <= self.x_length)
            self.assertTrue(diff[1] <= self.y_length)
            distances[to_[0]] = np.linalg.norm(diff)
        self.assertAll(distances >= 0, "Each postsyn cell has a connection")

    def test_distance(self):
        _test_distance_to_glomerulus(self)

    def test_fallback_to_closest_when_no_presyn_in_roi(self):
        """When no mossy fiber is within the x/y box, the closest one is connected."""
        _test_fallback_to_closest_presyn(
            self,
            conn_name="mf_to_glom",
            strategy="cerebellar_models.connectome.to_glomerulus.ConnectomeMossyGlomerulus",
            # x_length=5: both pre cells are at x-distance 40 and 45 → outside the box.
            extra_params={"x_length": 5, "y_length": 5},
        )

    def test_every_presyn_connected(self):
        _test_every_presyn_connected(self, {"x_length": 10, "y_length": 10})


class TestConnectomeUBCGlomerulus(
    RandomStorageFixture,
    NetworkFixture,
    FixedPosConfigFixture,
    NumpyTestCase,
    unittest.TestCase,
    engine_name="hdf5",
):
    def setUp(self):
        super().setUp()
        self.radius = 45
        self.cfg.connectivity.add(
            "x_to_glomerulus",
            dict(
                strategy="cerebellar_models.connectome.to_glomerulus.ConnectomeUBCGlomerulus",
                presynaptic=dict(cell_types=["test_cell"]),
                postsynaptic=dict(cell_types=["test_cell"]),
                radius=self.radius,
            ),
        )
        self.network = Scaffold(self.cfg, self.storage)

    def test_connect(self):
        self.network.compile()
        cs = self.network.get_connectivity_set("x_to_glomerulus")
        cell_positions = self.network.get_placement_set("test_cell").load_positions()

        cell_targets = np.array([-1, -1])
        distances = np.full(len(cell_positions), -1.0)
        self.assertEqual(len(cs), len(cell_positions), "As many connection as postsyn cell")
        for from_, to_ in cs.load_connections().as_globals():
            self.assertAll(from_[1:] == cell_targets)
            self.assertAll(to_[1:] == cell_targets)
            diff = np.linalg.norm(cell_positions[from_[0]] - cell_positions[to_[0]])
            self.assertTrue(diff <= self.radius)
            distances[to_[0]] = diff
        self.assertAll(distances >= 0, "Each postsyn cell has a connection")

    def test_distance(self):
        _test_distance_to_glomerulus(self)

    def test_fallback_to_closest_when_no_presyn_in_roi(self):
        """When no UBC is within the sphere, the closest one is connected."""
        _test_fallback_to_closest_presyn(
            self,
            conn_name="ubc_to_glom",
            strategy="cerebellar_models.connectome.to_glomerulus.ConnectomeUBCGlomerulus",
            # radius=5: both pre cells are at distance 40 and 45 → outside the sphere.
            extra_params={"radius": 5},
        )

    def test_every_presyn_connected(self):
        _test_every_presyn_connected(self, {"radius": 10})
