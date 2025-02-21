import functools
from collections import deque

import numpy as np
from bsb import AllenStructure, NrrdDependencyNode, config
from bsb.config._attrs import cfgdict
from scipy.spatial.transform import Rotation

from cerebellum.placement.utils import (
    boundaries_index_of,
    bresenham_line,
    signed_modulo,
)


class RotationReminder:
    """Utility class to keep information on last applied rotations."""

    def __init__(self, last_rotation, old_diff_rotation, rotation_to_correct=None):
        """

        :param last_rotation: Last rotation applied to the morphology's points.
        """
        self.last_rotation = last_rotation
        if rotation_to_correct is None:
            self.rotation_to_correct = np.zeros(3)
        else:
            self.rotation_to_correct = np.copy(rotation_to_correct)
        self.old_diff_rotation = np.copy(old_diff_rotation)

    def copy(self):
        return RotationReminder(
            self.last_rotation, self.old_diff_rotation, self.rotation_to_correct
        )


class MorphologyBender:
    default_depth: cfgdict[str, float] = config.dict(
        required=False, type=float, default={"mo": 150.0, "gr": 150.0, "pu": 150}
    )

    rescale: list[str] = config.list(required=False, type=str, default=["axon", "dendrites"])

    deform: list[str] = config.list(required=False, type=str, default=["axon", "dendrites"])

    fixed_dimension: int = config.attr(required=False, type=int, default=-1)
    """axis on which the orientation field will not be considered."""

    no_turn_back: bool = config.attr(required=False, type=bool, default=True)

    partition: AllenStructure = None

    @property
    def region_map(self):
        """
        Return RegionMap instance to manipulate the Allen mouse brain region hierarchy.

        :rtype: voxcell.region_mdeformed_morphoap.RegionMap
        """
        return self.partition.region_map

    @property
    def annotations(self):
        """
        Return the mouse brain annotation voxel array

        :rtype: numpy.ndarray
        """
        return self.partition.annotations.raw

    @functools.cached_property
    def orientation_field(self):
        """
        Return the brain orientation field pointing towards the outer shell of molecular layer

        :rtype: numpy.ndarray
        """
        loc_orient = self.partition.datasets["orientations"].raw
        loc_orient /= np.linalg.norm(loc_orient, axis=3)[..., np.newaxis]
        if 0 <= self.fixed_dimension <= 2:
            loc_orient[..., self.fixed_dimension] = 0.0
            loc_orient /= np.linalg.norm(loc_orient, axis=3)[..., np.newaxis]
        return loc_orient

    @functools.cached_property
    def thicknesses(self):
        """
        Return the brain depth field, i.e. the distance of each voxel from its layer boundaries

        :rtype: numpy.ndarray
        """
        return self.partition.datasets["thicknesses"].raw * self.partition.voxel_size

    @functools.cached_property
    def boundaries(self):
        """
        Return a boolean array which tells for each voxel transition (3,3,3), if it remains in the
        current orientation space.

        :rtype: numpy.ndarray
        """
        return self.partition.datasets["boundaries"].raw.reshape(self.annotations.shape + (3, 3, 3))

    def get_lay_abv(self, point):
        """
        Return the annotation layer abbreviation at the point location

        :param numpy.ndarray point: point location
        :return: layer abbreviation
        :rtype: str
        """
        return self.region_map.get(
            self.partition.mask_source.voxel_data_of(point, self.annotations), "acronym"
        )[-2:]

    def _ann_to_abv(self, id_reg):
        """
        Return layer index for thickness estimation and abbreviation.
        Return Nones if region is not part of Cerebellar cortex

        :param int id_reg: region id
        :return: layer index in thickness and its abbreviation.
        :rtype: Tuple(int, str) | Tuple(None, None)
        """
        expected = {"mo": 0, "pu": 1, "gr": 1}
        if id_reg == 0 or id_reg is None:
            return None, None
        lay = self.region_map.get(id_reg, "acronym")[-2:]
        for exp, i in expected.items():
            if exp in lay:
                return i, exp
        return None, None

    def test_voxels_between(self, old_vox, new_vox):
        """
        Check that every voxel between a source and target remain in the source rotation
        space.

        :param numpy.ndarray old_vox: starting voxel
        :param numpy.ndarray new_vox: target voxel
        :return: True if a border has been crossed, False otherwise.
        :rtype: bool
        """
        if np.any(new_vox != old_vox):
            last_vox = old_vox
            for voxel in bresenham_line(old_vox, new_vox)[1:]:
                voxel = np.array(voxel)
                if not self.boundaries[last_vox[0], last_vox[1], last_vox[2]][
                    boundaries_index_of(last_vox, voxel)
                ]:
                    # we hit the border of the region
                    return True
                last_vox = np.copy(voxel)
        return False

    def is_target_wrong(self, source, new_target, branch_labels=None):
        """
        Check if the segment between source and new_target remains in the correct rotation space.

        :param numpy.ndarray source: source point
        :param numpy.ndarray new_target: target point
        :param list[str] branch_labels: list of labels attached to the current segment.
        :return: True if the target is incorrect in the source rotation space.
        :rtype: bool
        """
        return self.test_voxels_between(
            self.partition.mask_source.voxel_of(source),
            self.partition.mask_source.voxel_of(new_target),
        )

    def delete_point(self, branch, i):
        """
        Delete a point in a morphology branch and translate the subsequent points and children
        branches to its parent location.

        :param bsb.morphologies.Branch branch: morphology branch to modify
        :param int i: index of point to delete
        :return: number of points deleted
        :rtype: int
        """
        delta = (branch.points[i - 1] - branch.points[i]) if i > 0 else np.zeros(3)
        # translate all the points of the branch starting at i
        branch.points[i + 1 :] += delta
        # translate all the children of the branch
        for child in branch.children:
            child.translate(delta)
        branch.delete_point(i)

    def rotate_point(self, source, branch, i, old_rots):
        """
        Compute the rotation to apply at a source point so that it follows the change in the
        orientation field while making sure the target remains within the frontiers of the region.

        :param numpy.ndarray source: center of the rotation to apply
        :param bsb.morphologies.Branch branch: Branch on which the rotation should be applied
        :param int i: current index in branch
        :param RotationReminder old_rots: Previous rotations applied to the previous points.
        :return: Euler angle of the rotation to apply at the source point
        :rtype: numpy.ndarray
        """
        max_angle = np.pi / 2 if self.no_turn_back else np.pi
        target = branch.points[i]
        to_rotate = True
        new_rotation = self.partition.mask_source.voxel_rotation_of(self.orientation_field, source)
        diff_rotation = signed_modulo(
            new_rotation.as_euler("xyz") - old_rots.last_rotation.as_euler("xyz"), 2 * np.pi
        )
        scaled_diff_rotation = None
        inc = 1.0
        branch_labels = list(branch.labelsets[branch.labels[i]])
        while to_rotate:
            scaled_diff_rotation = diff_rotation * inc
            if (np.absolute(scaled_diff_rotation) > max_angle).any():
                if inc >= 1:
                    diff_rotation -= np.sign(diff_rotation) * 1e-3
                    inc = -1.0
                    continue
                else:
                    raise ValueError("Hit a wall. Stopping")
            if to_rotate := self.is_target_wrong(
                source,
                Rotation.from_euler("xyz", scaled_diff_rotation).apply(target - source) + source,
                branch_labels,
            ):
                if np.linalg.norm(diff_rotation) == 0:
                    diff_rotation = np.copy(old_rots.old_diff_rotation)
                else:
                    inc += inc / 4

        if (
            np.linalg.norm(old_rots.rotation_to_correct) > 1e-5
            and inc == 1.0
            and not self.is_target_wrong(
                source,
                Rotation.from_euler(
                    "xyz", scaled_diff_rotation + old_rots.rotation_to_correct
                ).apply(target - source)
                + source,
                branch_labels,
            )
        ):
            scaled_diff_rotation += old_rots.rotation_to_correct
        # Update previous rotations.
        old_rots.last_rotation = new_rotation
        if np.linalg.norm(diff_rotation) > 0:
            old_rots.old_diff_rotation = np.copy(diff_rotation)
        old_rots.rotation_to_correct = signed_modulo(
            old_rots.rotation_to_correct - scaled_diff_rotation + diff_rotation, 2 * np.pi
        )
        return scaled_diff_rotation

    def process_scaling(self, point):
        """
        Calculate the local scaling factor to apply to a morphology segment based on the local
        thickness at its location.

        :param numpy.ndarray point: segment location
        :return: scaling factor at the current location
        :rtype: float
        """
        curr_ann = self.partition.mask_source.voxel_data_of(point, self.annotations)
        thick, lay = self._ann_to_abv(curr_ann)
        if lay is not None:
            return np.maximum(
                np.sum(
                    self.partition.mask_source.voxel_data_of(point, self.thicknesses)[
                        thick : thick + 2
                    ]
                )
                / self.default_depth[lay],
                0.1,
            )
        else:  # out of the annotations / depth / orientations fields.
            return 0.1

    def scale_morpho(self, branch, i, scaling):
        """
        Scale a morphology's branch at a position i, if the new position remains within
        the frontiers of the region.

        :param bsb.morphologies.Branch branch: morphology's branch to scale
        :param int i: position in the branch
        :param float scaling: factor of scaling.
        :return: False if the scaling could be performed, True else.
            Returns also the updated scaling
        :rtype: Tuple(bool, float)
        """

        # Update the scaling
        new_scaling = self.process_scaling(branch.points[i])
        if not np.isnan(new_scaling) and not np.isinf(new_scaling):
            scaling = new_scaling
        old_coord = np.copy(branch.points[i])
        new_coord = branch.points[i - 1] + (old_coord - branch.points[i - 1]) * scaling
        if 0 <= self.fixed_dimension <= 2:
            new_coord[self.fixed_dimension] = old_coord[self.fixed_dimension]
        # check that every voxel between the points remain within the region boundary
        rescale = self.test_voxels_between(
            self.partition.mask_source.voxel_of(branch.points[i - 1]),
            self.partition.mask_source.voxel_of(new_coord),
        )
        # if scaling resulted in an overshoot, set to old point coordinate
        branch.points[i] = np.copy(new_coord) if not rescale else np.copy(branch.points[i - 1])
        # translate all the points of the branch starting at i
        branch.points[i + 1 :] += branch.points[i] - old_coord
        # translate all the children of the branch
        for child in branch.children:
            child.translate(branch.points[i] - old_coord)

        return rescale, scaling

    def _init_stack(self, morphology):
        """
        Initialize stack and perform the roots' rotation.
        :param bsb.morphologies.Morphology morphology: Morphology to deform
        :return: stack of branch to deform
        """
        stack_data = []
        for branch in morphology.roots:
            try:
                rotation = self.partition.mask_source.voxel_rotation_of(
                    self.orientation_field,
                    branch.points[0],
                )
                branch.root_rotate(rotation)
                curr_scaling = self.process_scaling(branch.points[0])
                axis_max = max(
                    enumerate(
                        np.absolute(
                            self.partition.mask_source.voxel_data_of(
                                branch.points[0], self.orientation_field
                            )
                        )
                    ),
                    key=lambda x: x[1],
                )[0]
                old_diff_rotation = np.zeros(3)
                old_diff_rotation[axis_max] = 1e-2
                stack_data.append(
                    (branch, RotationReminder(rotation, old_diff_rotation), curr_scaling)
                )
            except ValueError as _:
                continue
            return stack_data
        return stack_data

    def deform_morphology(self, morphology):
        """
        Visit each point of the morphology and rescale and or deform them according to the
        orientation field and the space available.
        :param bsb.morphologies.Morphology morphology: Morphology to deform
        :return: deformed morphology
        :rtype: bsb.morphologies.Morphology
        """
        stack_data = self._init_stack(morphology)
        stack = deque(stack_data)

        while True:
            try:
                branch, old_rots, curr_scaling = stack.pop()
            except IndexError:
                break
            else:
                last_point = branch.points[0]
                last_index = 0
                i = 0
                while i < len(branch.points):
                    branch_labels = branch.labelsets[branch.labels[i]]
                    if np.isin(list(branch_labels), self.rescale).any() and i > 0:
                        rescale, curr_scaling = self.scale_morpho(branch, i, curr_scaling)
                        if rescale:
                            self.delete_point(branch, i)
                            continue
                    if np.isin(list(branch_labels), self.deform).any():
                        try:
                            rotation = Rotation.from_euler(
                                "xyz", self.rotate_point(last_point, branch, i, old_rots)
                            )
                            branch.root_rotate(rotation, downstream_of=last_index)
                        except ValueError as _:
                            self.delete_point(branch, i)
                            continue
                    last_point = branch.points[i]
                    last_index = i
                    i += 1

                if len(branch.points) == 0:  # removed all the points, so we cut the branch
                    branch.detach()
                else:
                    if len(branch.children):
                        stack.extend(
                            [
                                (child, old_rots.copy(), curr_scaling)
                                for child in branch.children
                                if NrrdDependencyNode.is_within(
                                    self.partition.mask_source.voxel_of(child.points[0]),
                                    self.annotations,
                                )
                            ]
                        )
        morphology.close_gaps()
        return morphology

    def process(self, positions, morphologies):
        """
        Generate a list of deformed morphologies for each corresponding voxel at the provided
        locations.

        :param List[numpy.ndarray positions: List of cell positions
        :param List[bsb.morphologies.Morphology] morphologies: list of morphologies to associate at
            each location. Each location will randomly choose one of the available morphologies.
        :return: A list of deformed morphologies for each location provided.
        :rtype: List[bsb.morphologies.Morphology]
        """
        if len(positions) == 0:
            return []
        morpho_ids = np.random.default_rng().integers(
            np.unique(morphologies).size, size=len(positions)
        )
        morphology_list = morphologies[morpho_ids]
        deformed_list = np.zeros_like(morphology_list)
        voxel_pos = self.partition.mask_source.voxel_of(positions)
        uniques, indexes = np.unique(voxel_pos, axis=0, return_inverse=True)

        for i, uniq_vox in enumerate(uniques):
            filter_pos = indexes == i
            u_morpho, u_index = np.unique(morpho_ids[filter_pos], return_inverse=True)
            u_morpho = morphology_list[u_morpho]
            # filter for positions inside the orientation and depth field.
            if NrrdDependencyNode.is_within(uniq_vox, self.orientation_field):
                translation_vec = (
                    uniq_vox + 0.5
                ) * self.partition.voxel_size + self.partition.mask_source.space_origin
                for j, morpho in enumerate(u_morpho):
                    deformed_morpho = morpho.copy()
                    deformed_morpho.translate(translation_vec)
                    # Bend the morphology according to orientation field.
                    try:
                        deformed_morpho = self.deform_morphology(deformed_morpho)
                    except Exception as e:
                        print(f"Error with morphology: {morpho._meta['name']} at {uniq_vox}.")
                        raise e
                    deformed_morpho.center()
                    u_morpho[j] = deformed_morpho
            deformed_list[filter_pos] = u_morpho[u_index]

        return deformed_list
