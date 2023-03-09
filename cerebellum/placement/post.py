from bsb import config
from bsb.config import types
from bsb.postprocessing import PostProcessingHook

import numpy as np
from scipy.stats import truncnorm


@config.node
class LabelMicrozones(PostProcessingHook):
    targets = config.attr(type=types.list(), required=True)

    def after_placement(self):
        # Divide the volume into two sub-parts (one positive and one negative)
        for neurons_2b_labeled in self.targets:
            ps = self.scaffold.get_placement_set(neurons_2b_labeled)
            ids = ps.load_identifiers()
            zeds = ps.load_positions()[:, 2]
            z_sep = np.median(zeds)
            index_pos = np.where(zeds >= z_sep)[0]
            index_neg = np.where(zeds < z_sep)[0]

            labels = {
                "microzone-positive": ids[index_pos],
                "microzone-negative": ids[index_neg],
            }

            self.scaffold.label_cells(
                ids[index_pos],
                label="microzone-positive",
            )
            self.scaffold.label_cells(
                ids[index_neg],
                label="microzone-negative",
            )

            self.label_satellites(neurons_2b_labeled, labels)

    def label_satellites(self, planet_type, labels):
        for possible_satellites in self.scaffold.get_cell_types():
            # Find all cell types that specify this type as their planet type
            if (
                    hasattr(possible_satellites.placement, "planet_types")
                    and planet_type in possible_satellites.placement.planet_types
            ):
                # Get the IDs of this sattelite cell type.
                ps = self.scaffold.get_placement_set(possible_satellites.name)
                satellites = ps.identifiers
                # Retrieve the planet map for this satellite type. A planet map is an
                # array that lists the planet for each satellite. `sattelite_map[n]`` will
                # hold the planet ID for sattelite `n`, where `n` the index of the
                # satellites in their cell type, not their scaffold ID.
                satellite_map = self.scaffold._planets[possible_satellites.name].copy()
                # Create counters for each label for the report below
                satellite_label_count = {l: 0 for l in labels.keys()}
                # Iterate each label to check for any planets with that label, and label
                # their sattelite with the same label. After iterating all labels, each
                # satellite should have the same labels as their planet.
                for label, labelled_cells in labels.items():
                    for i, satellite in enumerate(satellites):
                        planet = satellite_map[i]
                        if planet in labelled_cells:
                            self.scaffold.label_cells([satellite], label=label)
                            # Increase the counter of this label
                            satellite_label_count[label] += 1


class AscendingAxonLengths(PostProcessingHook):
    def after_placement(self):
        granule_type = self.scaffold.cell_types.granule_cell
        granules = self.scaffold.get_placement_set(granule_type.name)
        granule_geometry = granule_type.spatial.geometry
        parallel_fibers = np.zeros((len(granules)))
        pf_height = granule_geometry.pf_height
        pf_height_sd = granule_geometry.pf_height_sd
        molecular_layer = self.scaffold.partitions.molecular_layer
        floor_ml = molecular_layer.boundaries.y
        roof_ml = floor_ml + molecular_layer.boundaries.height

        for idx, granule in enumerate(granules.load_cells()):
            granule_y = granule.position[1]
            # Determine min and max height so that the parallel fiber stays inside the
            # molecular layer
            pf_height_min = floor_ml - granule_y
            pf_height_max = roof_ml - granule_y
            # Determine the shape parameters a and b of the truncated normal
            # distribution. See https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.truncnorm.html
            a, b = (
                (pf_height_min - pf_height) / pf_height_sd,
                (pf_height_max - pf_height) / pf_height_sd,
            )
            # Draw a sample for the parallel fiber height from a truncated normal
            # distribution with sd `pf_height_sd` and mean `pf_height`, truncated by
            # the molecular layer bounds.
            parallel_fibers[idx] = (
                    truncnorm.rvs(a, b, size=1) * pf_height_sd + pf_height
            )
        granules.append_additional(
            chunk, "ascending_axon_lengths", data=parallel_fibers
        )


class DCNRotations(PostProcessingHook):
    """
    Create a matrix of planes tilted between -45° and 45°,
    storing id and the planar coefficients a, b, c and d for each DCN cell
    """

    def after_placement(self):
        ps = self.scaffold.get_placement_set("dcn_cell")
        dend_tree_coeff = np.zeros((len(ps), 4))
        positions = ps.positions
        for i in range(len(ps)):
            # Make the planar coefficients a, b and c.
            dend_tree_coeff[i] = np.random.rand(4) * 2.0 - 1.0
            # Calculate the last planar coefficient d from ax + by + cz - d = 0
            # => d = - (ax + by + cz)
            dend_tree_coeff[i, 3] = -np.sum(dend_tree_coeff[i, 0:3] * positions[i, 0:3])
        ps.create_additional("dcn_orientations", data=dend_tree_coeff)
