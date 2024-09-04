import numpy as np
from bsb import AllenStructure, MorphologyGenerator, TopologyError, config
from bsb._util import rotation_matrix_from_vectors
from scipy.spatial.transform import Rotation

from cerebellum.placement.morphology_bender import MorphologyBender, signed_modulo
from cerebellum.placement.utils import get_diff_angle


@config.node
class BenderGenerator(MorphologyBender, MorphologyGenerator, classmap_entry="morpho_bender"):
    """BSB MorphologyGenerator based on morphology bending"""

    may_be_empty = config.provide(True)

    def generate(self, positions, morphologies, context):
        if len(context.partitions) > 1 or not isinstance(context.partitions[0], AllenStructure):
            raise TopologyError("BenderGenerator should be used with a single AllenAtlas partition")
        self.cell_type = context.indicator.cell_type
        self.partition = context.partitions[0]
        return self.process(
            positions,
            np.array([m.load() for m in morphologies]),
        )


@config.node
class GranuleBender(MorphologyBender):
    """Specific bender for granule cell."""

    SIZE_ASCENDING_AXON = 140.0  # um
    NB_ASCENDING_AXON_SEGMENTS = 12  # +1 to include the soma
    nb_section_left = NB_ASCENDING_AXON_SEGMENTS

    rescale = False
    ratio_gr = None

    def process_scaling(self, point):
        # We calculate the scaling ratio so that the tip of the ascending axon's depth ratio within
        # the molecular layer is equal to the depth ratio of the soma within the granular layer.
        distances = self.partition.voxel_data_of(point, self.thicknesses)
        top_dist = distances[1]  # dist to gr/mol boundary
        if self.ratio_gr is None:
            self.ratio_gr = top_dist / (top_dist + distances[2])  # ratio depth within gr.
        if "mo" in self.get_lay_abv(point):
            ascending_axon_length = (1 - self.ratio_gr) * (distances[0] + top_dist) - top_dist
        else:
            ascending_axon_length = top_dist + (1 - self.ratio_gr) * (distances[0] - top_dist)

        self.scaling = np.maximum(
            ascending_axon_length
            / self.nb_section_left
            / self.SIZE_ASCENDING_AXON
            * self.NB_ASCENDING_AXON_SEGMENTS,
            0.1,
        )
        self.nb_section_left -= 1
        return self.scaling  # ratio according to default axon length

    def scale_morpho(self, branch, i, scaling):
        old_deriv = branch.points[i] - branch.points[i - 1]
        rescale, new_scaling = super().scale_morpho(branch, i, scaling)
        self.has_rescale = False
        if not rescale:
            if self.nb_section_left <= 0 and self.get_lay_abv(branch.points[i]) != "mo":
                # If the last point does not land on molecular layer, add an ascending axon point.
                self.nb_section_left = 1
                branch.introduce_point(i, branch.points[i])
                branch.points[i + 1] = branch.points[i] + old_deriv
                # translate all the children of the branch
                for child in branch.children:
                    child.translate(old_deriv)

            # If the difference of orientation between the current point of the ascending axon
            # according to its first point yields an angle greater than 90 degrees, the ascending
            # axon has crossed a frontier of the region, and its last two points can be removed.
            rescale = np.any(
                np.absolute(
                    get_diff_angle(
                        self.orientation_field,
                        self.partition.voxel_of(branch.points[i]),
                        self.partition.voxel_of(branch.points[0]),
                    )
                )
                > np.pi / 2
            )
            if rescale:
                self.has_rescale = True
                old_coord = np.copy(branch.points[i])
                branch.points[i] = np.copy(branch.points[i - 2])
                # translate all the points of the branch starting at i
                branch.points[i + 1 :] += branch.points[i] - old_coord
                # translate all the children of the branch
                for child in branch.children:
                    child.points += branch.points[i] - old_coord
        return rescale, new_scaling

    def delete_point(self, branch, i):
        # Delete 2 points instead of 1
        branch.delete_point(i - 1)
        if self.has_rescale:
            branch.delete_point(i - 1)
            self.has_rescale = False

    def rotate_point(self, source, branch, i, old_rots):
        rotation = super().rotate_point(source, branch, i, old_rots)
        if np.isin(list(branch.labelsets[branch.labels[i]]), ["parallel_fiber"]).any():
            target = branch.points[i]
            target_ = Rotation.from_euler("xyz", rotation).apply(target - source) + source
            distances = self.partition.voxel_data_of(target_, self.thicknesses)
            if np.sum(distances[:1]) > 1e-6 and np.linalg.norm(source - target) > 1e-3:
                ratio_mol = distances[0] / (distances[1] + distances[0])
                if ratio_mol > self.ratio_gr + (1 - self.ratio_gr) / 3:  # go too deep
                    angle = signed_modulo(
                        np.asarray(
                            Rotation.from_matrix(
                                rotation_matrix_from_vectors(
                                    self.partition.voxel_orient(
                                        self.orientation_field,
                                        target_,
                                    ),
                                    target_ - source,
                                )
                            ).as_euler("xyz")
                        ),
                        2 * np.pi,
                    )
                    # bring back the fiber at 90 degree with respect to the orientation field.
                    # works only on parallel fibers !
                    correction_angle = np.sign(angle) * np.pi / 2 - angle
                    next_loc = (
                        Rotation.from_euler("xyz", rotation + correction_angle).apply(
                            target - source
                        )
                        + source
                    )
                    distances = self.partition.voxel_data_of(next_loc, self.thicknesses)
                    if (
                        not self.is_target_wrong(source, next_loc)
                        and distances[0] / (distances[1] + distances[0]) < ratio_mol
                    ):
                        return rotation + correction_angle
        return rotation

    def deform_morphology(self, morphology):
        self.nb_section_left = self.NB_ASCENDING_AXON_SEGMENTS
        self.ratio_gr = None
        return super().deform_morphology(morphology)


@config.node
class GolgiGenerator(BenderGenerator, classmap_entry="golgi_bender"):
    """Specific bender for golgi cell."""

    def deform_morphology(self, morphology):
        morphology = super().deform_morphology(morphology)
        # Re-label the dendrites based on their layer location.
        id_den_gr = id_den_mol = -1
        for k, v in morphology.labelsets.items():
            if "dendrites" in v:
                if "basal" in v:
                    id_den_gr = k
                elif "apical" in v:
                    id_den_mol = k
        for branch in morphology.branches:
            is_dendrite = np.array(
                [
                    np.isin(list(branch.labelsets[branch.labels[i]]), ["dendrites"]).any()
                    for i in range(len(branch.points))
                ]
            )
            if np.any(is_dendrite):
                branch.labels[is_dendrite] = np.array(
                    [
                        (
                            id_den_mol
                            if (
                                self.partition.voxel_data_of(point, self.annotations) != 0
                                and "mo" in self.get_lay_abv(point)
                            )
                            else id_den_gr
                        )
                        for point in branch.points[is_dendrite]
                    ]
                )
        return morphology


@config.node
class BasketGenerator(BenderGenerator, classmap_entry="basket_bender"):
    """Specific bender for basket cell."""

    def is_target_wrong(self, source, new_target, branch_labels=None):
        if super().is_target_wrong(source, new_target, branch_labels):
            return True

        # dendrites should remain in molecular layer
        current_abv = self.get_lay_abv(new_target)
        if branch_labels is not None and np.all(["axon" not in label for label in branch_labels]):
            return "mo" not in current_abv

        # axon should get closer to Purkinje layer.
        current_dist = self.partition.voxel_data_of(new_target, self.thicknesses)
        old_dist = self.partition.voxel_data_of(source, self.thicknesses)
        if "mo" in current_abv:
            return current_dist[1] > 62.5 and old_dist[1] < current_dist[1]
        elif "gr" in current_abv:
            return 37.5 < current_dist[1] < old_dist[1]


@config.node
class PurkinjeGenerator(BenderGenerator, classmap_entry="purkinje_bender"):
    """Specific bender for Purkinje cell."""

    def is_target_wrong(self, source, new_target, branch_labels=None):
        if super().is_target_wrong(source, new_target, branch_labels):
            return True

        # dendrites should remain in molecular and purkinje layer
        current_abv = self.get_lay_abv(new_target)
        if branch_labels is not None and np.all(["axon" not in label for label in branch_labels]):
            return "mo" not in current_abv and "pu" not in current_abv
        return False
