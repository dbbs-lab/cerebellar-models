import os
import unittest

import numpy as np
from bsb import Scaffold


def get_tiny_network():
    scaffold_filename = "_test_tiny_network.hdf5"
    scaffold_path = os.path.join(os.path.dirname(__file__), "..", scaffold_filename)
    if not os.path.exists(scaffold_path):
        config = JSONConfig(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                "configs",
                "3_9_mouse.json",
            )
        )
        config.output_formatter.file = scaffold_path
        scaffold = Scaffold(config)
        mf_ids = scaffold.create_entities(scaffold.configuration.cell_types["mossy_fibers"], 3)
        glom_ids = scaffold.place_cells(
            scaffold.configuration.cell_types["glomerulus"],
            None,
            np.array([[0, 0, 0]] * 12),
        )
        grc_ids = scaffold.place_cells(
            scaffold.configuration.cell_types["granule_cell"],
            None,
            np.array([[0, 0, 0]] * 12),
        )
        goc_ids = scaffold.place_cells(
            scaffold.configuration.cell_types["golgi_cell"],
            None,
            np.array([[0, 0, 0]] * 2),
        )
        pc_ids = scaffold.place_cells(
            scaffold.configuration.cell_types["purkinje_cell"],
            None,
            np.array([[0, 0, 0]]),
        )
        sc_ids = scaffold.place_cells(
            scaffold.configuration.cell_types["stellate_cell"],
            None,
            np.array([[0, 0, 0]] * 3),
        )
        bc_ids = scaffold.place_cells(
            scaffold.configuration.cell_types["basket_cell"],
            None,
            np.array([[0, 0, 0]] * 3),
        )

        scaffold.connect_cells(
            scaffold.configuration.connection_types["mossy_to_glomerulus"],
            np.array(
                [
                    [mf_ids[0], glom_ids[0 * 4]],
                    [mf_ids[0], glom_ids[0 * 4 + 1]],
                    [mf_ids[0], glom_ids[0 * 4 + 2]],
                    [mf_ids[0], glom_ids[0 * 4 + 3]],
                    [mf_ids[1], glom_ids[1 * 4]],
                    [mf_ids[1], glom_ids[1 * 4 + 1]],
                    [mf_ids[1], glom_ids[1 * 4 + 2]],
                    [mf_ids[1], glom_ids[1 * 4 + 3]],
                    [mf_ids[2], glom_ids[2 * 4]],
                    [mf_ids[2], glom_ids[2 * 4 + 1]],
                    [mf_ids[2], glom_ids[2 * 4 + 2]],
                    [mf_ids[2], glom_ids[2 * 4 + 3]],
                ]
            ),
        )
        scaffold.connect_cells(
            scaffold.configuration.connection_types["glomerulus_to_granule"],
            np.array(
                [
                    [glom_ids[0], grc_ids[0]],
                    [glom_ids[0], grc_ids[1]],
                    [glom_ids[0], grc_ids[2]],
                    [glom_ids[0], grc_ids[3]],
                    [glom_ids[0], grc_ids[4]],
                    [glom_ids[1], grc_ids[0]],
                    [glom_ids[1], grc_ids[1]],
                    [glom_ids[1], grc_ids[2]],
                    [glom_ids[1], grc_ids[3]],
                    [glom_ids[1], grc_ids[4]],
                    [glom_ids[2], grc_ids[0]],
                    [glom_ids[2], grc_ids[1]],
                    [glom_ids[2], grc_ids[2]],
                    [glom_ids[2], grc_ids[3]],
                    [glom_ids[2], grc_ids[4]],
                    [glom_ids[3], grc_ids[5]],
                    [glom_ids[3], grc_ids[6]],
                    [glom_ids[3], grc_ids[7]],
                    [glom_ids[3], grc_ids[8]],
                    [glom_ids[3], grc_ids[9]],
                    [glom_ids[4], grc_ids[5]],
                    [glom_ids[4], grc_ids[6]],
                    [glom_ids[4], grc_ids[7]],
                    [glom_ids[4], grc_ids[8]],
                    [glom_ids[4], grc_ids[9]],
                    [glom_ids[5], grc_ids[5]],
                    [glom_ids[5], grc_ids[6]],
                    [glom_ids[5], grc_ids[7]],
                    [glom_ids[5], grc_ids[8]],
                    [glom_ids[5], grc_ids[9]],
                    [glom_ids[6], grc_ids[5]],
                    [glom_ids[6], grc_ids[6]],
                    [glom_ids[6], grc_ids[7]],
                    [glom_ids[6], grc_ids[8]],
                    [glom_ids[6], grc_ids[9]],
                    [glom_ids[7], grc_ids[10]],
                    [glom_ids[8], grc_ids[11]],
                    [glom_ids[9], grc_ids[10]],
                    [glom_ids[10], grc_ids[0]],
                    [glom_ids[11], grc_ids[1]],
                ]
            ),
            compartments=(
                c := np.array(
                    [
                        *([[-1, 113]] * 5),
                        *([[-1, 122]] * 5),
                        *([[-1, 131]] * 5),
                        *([[-1, 113]] * 5),
                        *([[-1, 122]] * 5),
                        *([[-1, 131]] * 5),
                        *([[-1, 140]] * 5),
                        *([[-1, 113]] * 2),
                        [-1, 122],
                        [-1, 140],
                        [-1, 140],
                    ]
                )
            ),
            morphologies=np.array([["GranuleCell"] * 2] * c.shape[0]),
        )

        scaffold.connect_cells(
            scaffold.configuration.connection_types["parallel_fiber_to_golgi"],
            np.array(
                [
                    [grc_ids[0], goc_ids[0]],
                    [grc_ids[1], goc_ids[0]],
                    [grc_ids[1], goc_ids[0]],
                    [grc_ids[3], goc_ids[0]],
                    [grc_ids[4], goc_ids[0]],
                    [grc_ids[5], goc_ids[0]],
                    [grc_ids[6], goc_ids[0]],
                    [grc_ids[7], goc_ids[0]],
                    [grc_ids[8], goc_ids[0]],
                    [grc_ids[9], goc_ids[0]],
                    [grc_ids[10], goc_ids[0]],
                    [grc_ids[11], goc_ids[0]],
                    [grc_ids[0], goc_ids[1]],
                    [grc_ids[1], goc_ids[1]],
                    [grc_ids[1], goc_ids[1]],
                    [grc_ids[3], goc_ids[1]],
                    [grc_ids[4], goc_ids[1]],
                    [grc_ids[5], goc_ids[1]],
                    [grc_ids[6], goc_ids[1]],
                    [grc_ids[7], goc_ids[1]],
                    [grc_ids[8], goc_ids[1]],
                    [grc_ids[9], goc_ids[1]],
                    [grc_ids[10], goc_ids[1]],
                    [grc_ids[11], goc_ids[1]],
                ]
            ),
            compartments=np.array(
                [
                    *([[58, 2930]] * 6),
                    *([[59, 2929]] * 6),
                    *([[60, 2928]] * 6),
                    *([[61, 2927]] * 6),
                ]
            ),
            morphologies=np.array([["GranuleCell", "GolgiCell"]] * 24),
        )
        scaffold.connect_cells(
            scaffold.configuration.connection_types["ascending_axon_to_golgi"],
            np.array(
                [
                    [grc_ids[0], goc_ids[1]],
                    [grc_ids[1], goc_ids[1]],
                    [grc_ids[1], goc_ids[1]],
                    [grc_ids[3], goc_ids[1]],
                    [grc_ids[4], goc_ids[1]],
                    [grc_ids[5], goc_ids[1]],
                    [grc_ids[6], goc_ids[1]],
                    [grc_ids[7], goc_ids[1]],
                    [grc_ids[8], goc_ids[1]],
                    [grc_ids[9], goc_ids[1]],
                    [grc_ids[10], goc_ids[1]],
                    [grc_ids[11], goc_ids[1]],
                    [grc_ids[0], goc_ids[0]],
                    [grc_ids[1], goc_ids[0]],
                    [grc_ids[1], goc_ids[0]],
                    [grc_ids[3], goc_ids[0]],
                    [grc_ids[4], goc_ids[0]],
                    [grc_ids[5], goc_ids[0]],
                    [grc_ids[6], goc_ids[0]],
                    [grc_ids[7], goc_ids[0]],
                    [grc_ids[8], goc_ids[0]],
                    [grc_ids[9], goc_ids[0]],
                    [grc_ids[10], goc_ids[0]],
                    [grc_ids[11], goc_ids[0]],
                ]
            ),
            compartments=np.array(
                [
                    *([[58, 2931]] * 6),
                    *([[59, 2932]] * 6),
                    *([[60, 2933]] * 6),
                    *([[61, 2802]] * 6),
                ]
            ),
            morphologies=np.array([["GranuleCell", "GolgiCell"]] * 24),
        )
        scaffold.connect_cells(
            scaffold.configuration.connection_types["golgi_to_golgi"],
            np.array([*([[goc_ids[0], goc_ids[1]]] * 100), *([[goc_ids[1], goc_ids[0]]] * 100)]),
            compartments=np.tile(
                np.column_stack((np.arange(396, 496), np.arange(495, 395, -1))), (2, 1)
            ),
            morphologies=np.tile(np.array(["GolgiCell"]), (200, 2)),
        )
        scaffold.connect_cells(
            scaffold.configuration.connection_types["gap_goc"],
            np.array([*([[goc_ids[0], goc_ids[1]]] * 5), *([[goc_ids[1], goc_ids[0]]] * 5)]),
            compartments=np.ones((10, 2)) * 395,
            morphologies=np.tile(["GolgiCell"], (10, 2)),
        )
        scaffold.connect_cells(
            scaffold.configuration.connection_types["parallel_fiber_to_stellate"],
            np.column_stack(
                (
                    np.tile(grc_ids, len(sc_ids)),
                    np.tile(sc_ids, (len(grc_ids), 1)).T.flatten(),
                )
            ),
            compartments=np.tile([80, 2580], (len(grc_ids) * len(sc_ids), 1)),
            morphologies=np.tile(["GranuleCell", "StellateCell"], (len(grc_ids) * len(sc_ids), 1)),
        )
        scaffold.connect_cells(
            scaffold.configuration.connection_types["parallel_fiber_to_basket"],
            np.column_stack(
                (
                    np.tile(grc_ids, len(bc_ids)),
                    np.tile(bc_ids, (len(grc_ids), 1)).T.flatten(),
                )
            ),
            compartments=np.tile([79, 2740], (len(grc_ids) * len(bc_ids), 1)),
            morphologies=np.tile(["GranuleCell", "BasketCell"], (len(grc_ids) * len(bc_ids), 1)),
        )
        scaffold.connect_cells(
            scaffold.configuration.connection_types["parallel_fiber_to_purkinje"],
            np.column_stack(
                (
                    np.tile(grc_ids, len(pc_ids)),
                    np.tile(pc_ids, (len(grc_ids), 1)).T.flatten(),
                )
            ),
            compartments=np.tile([81, 2371], (len(grc_ids) * len(pc_ids), 1)),
            morphologies=np.tile(["GranuleCell", "PurkinjeCell"], (len(grc_ids) * len(pc_ids), 1)),
        )
        scaffold.connect_cells(
            scaffold.configuration.connection_types["stellate_to_purkinje"],
            np.column_stack(
                (
                    np.tile(sc_ids, len(pc_ids)),
                    np.tile(pc_ids, (len(sc_ids), 1)).T.flatten(),
                )
            ),
            compartments=np.tile([2870, 2875], (len(sc_ids) * len(pc_ids), 1)),
            morphologies=np.tile(["StellateCell", "PurkinjeCell"], (len(sc_ids) * len(pc_ids), 1)),
        )
        scaffold.connect_cells(
            scaffold.configuration.connection_types["basket_to_purkinje"],
            np.column_stack(
                (
                    np.tile(bc_ids, len(pc_ids)),
                    np.tile(pc_ids, (len(bc_ids), 1)).T.flatten(),
                )
            ),
            compartments=np.tile([6129, 10], (len(bc_ids) * len(pc_ids), 1)),
            morphologies=np.tile(["BasketCell", "PurkinjeCell"], (len(bc_ids) * len(pc_ids), 1)),
        )
        scaffold.compile_output()
    else:
        scaffold = from_hdf5(scaffold_path)


@unittest.skip("Needs to be updated to bsb v4.0.1")
class TestCerebellumPlacement(unittest.TestCase):
    """
    Check if the placement of all cell types is correct
    """

    @classmethod
    def setUpClass(self):
        super().setUpClass()
        self.scaffold = get_test_network(200, 200)

    def test_bounds(self):
        # Create different sized test networks
        for dimensions in [[100, 100], [200, 200]]:
            scaffold = get_test_network(*dimensions)
            # Get all non-entity cell types.
            for cell_type in scaffold.get_cell_types(entities=False):
                # Start an out-of-bounds check.
                with self.subTest(
                    cell_type=cell_type.name,
                    dimensions=dimensions,
                    placement=cell_type.placement.__class__.__name__,
                ):
                    # Treat satellites differently because they can belong to
                    # multiple layers of their planet cell types.
                    if isinstance(cell_type.placement, Satellite):
                        self._test_satellite_bounds(cell_type)
                        continue
                    # Get layer bounds and cell positions.
                    layer = cell_type.placement.layer_instance
                    positions = scaffold.get_cells_by_type(cell_type.name)[:, 2:5]
                    min = layer.origin
                    max = layer.origin + layer.dimensions
                    # If any element in the positional data is smaller than its
                    # corresponding element in `min` or larger than that in `max`,
                    # np.where will add its index to an array which needs to be empty
                    # for the test to pass.
                    i = len(np.where(((positions < min) | (positions > max)))[0])
                    self.assertEqual(i, 0, "Cells placed out of bounds.")

    def _test_satellite_bounds(self, cell_type):
        # Use the after array to get all planet cell type layers.
        after = cell_type.placement.after
        planet_cell_types = [self.scaffold.get_cell_type(n) for n in after]
        layers = [
            planet_cell_type.placement.layer_instance for planet_cell_type in planet_cell_types
        ]
        # Get satellite positions
        positions = self.scaffold.get_cells_by_type(cell_type.name)[:, 2:5]
        in_bounds = set()
        # Check layer per layer whether the satellites are in a layer Update a set per
        # layer to keep track of all satellites that belong to at least one planet layer
        for layer in layers:
            min = layer.origin
            max = layer.origin + layer.dimensions
            in_layer_bounds = np.where(np.all((positions > min) & (positions < max), axis=1))[0]
            in_bounds.update(in_layer_bounds)
        # Verify all satellite cells belong to at least one planet layer.
        out_of_bounds = len(positions) - len(in_bounds)
        message = str(out_of_bounds) + " cells placed out of bounds"
        self.assertEqual(out_of_bounds, 0, message)

    def test_purkinje(self):
        import scipy.spatial.distance as dist

        config = self.scaffold.configuration
        layer = config.layers["purkinje_layer"]
        pc = config.cell_types["purkinje_cell"]
        self.scaffold.reset_network_cache()
        pc.placement.place()
        pcCount = self.scaffold.cells_by_type["purkinje_cell"].shape[0]
        density = pcCount / layer.width / layer.depth
        pc_pos = self.scaffold.cells_by_type["purkinje_cell"][:, [2, 3, 4]]
        Dist2D = dist.pdist(np.column_stack((pc_pos[:, 0], pc_pos[:, 2])), "euclidean")
        overlapSomata = np.where(Dist2D < 80 / 100 * pc.spatial.radius)[0]  #
        Dist1Dsqr = np.zeros((pc_pos.shape[0], pc_pos.shape[0], 2))
        Dist1Dsqr[:, :, 0] = dist.squareform(dist.pdist(pc_pos[:, [0]], "euclidean"))
        Dist1Dsqr[:, :, 1] = dist.squareform(dist.pdist(pc_pos[:, [2]], "euclidean"))
        overlapDend_whichPairs = np.where(
            np.logical_and(
                Dist1Dsqr[:, :, 0] < 80 / 100 * pc.placement.extension_x,
                Dist1Dsqr[:, :, 1] < 80 / 100 * pc.placement.extension_z,
            )
        )
        overlapDend_whichPairs = np.column_stack(
            (overlapDend_whichPairs[0], overlapDend_whichPairs[1])
        )
        overlapDend_whichPairs = overlapDend_whichPairs[
            np.where(overlapDend_whichPairs[:, 0] != overlapDend_whichPairs[:, 1])[0], :
        ]

        # Asserts
        self.assertAlmostEqual(overlapSomata.shape[0], 0, delta=pcCount * 1 / 100)
        self.assertAlmostEqual(overlapDend_whichPairs.shape[0] / 2, 0, delta=pcCount * 4 / 100)
