from collections import deque
from copy import deepcopy

import numpy as np
from bsb import AllenStructure, NrrdDependencyNode, config, pool_cache, types
from bsb.config._attrs import cfgdict
from bsb.voxels import voxel_rotation_of
from scipy.spatial.transform import Rotation

from cerebellar_models.placement.utils import boundaries_index_of, bresenham_line


class RotationReminder:
    """Utility class to keep information on last applied rotations."""

    def __init__(self, last_rotation, old_diff_rotation, rotation_to_correct=None):
        """

        :param last_rotation: Last rotation applied to the morphology's points.
        """
        self.last_rotation = deepcopy(last_rotation)
        self.original_rotation = deepcopy(last_rotation)
        if rotation_to_correct is None:
            self.rotation_to_correct = Rotation.from_euler("xyz", np.zeros(3))
        else:
            self.rotation_to_correct = deepcopy(rotation_to_correct)
        self.old_diff_rotation = deepcopy(old_diff_rotation)

    def copy(self):
        return RotationReminder(
            self.last_rotation, self.old_diff_rotation, self.rotation_to_correct
        )


def get_branch_labels(branch, index_branch):
    return list(branch.labelsets[branch.labels[index_branch]])


def has_label(branch_labels, label):
    return np.any([label in l for l in branch_labels])


class MorphologyBender:

    orientations_field: NrrdDependencyNode = config.ref(config.refs.vox_dset_ref)

    depths_field: NrrdDependencyNode = config.ref(config.refs.vox_dset_ref)

    boundaries_field: NrrdDependencyNode = config.ref(config.refs.vox_dset_ref)

    default_depth: cfgdict[str, float] = config.dict(
        required=False, type=float, default={"mo": 150.0, "pu": 165.0, "gr": 150.0}
    )
    """reference layers' thickness to use during rescale"""

    rescale: list[str] = config.list(required=False, type=str, default=["axon", "dendrites"])
    """list of labels to filter morphologies' branches to rescale"""

    deform: list[str] = config.list(required=False, type=str, default=["axon", "dendrites"])
    """list of labels to filter morphologies' branches to deform"""

    fixed_dimensions = config.attr(required=False, type=types.or_(int, dict[str:int]), default=-1)
    """axis `x` on which the orientation field will not be considered or dictionary 
    linking morphology label to their corresponding `x` axis"""

    no_turn_back: bool = config.attr(required=False, type=bool, default=True)
    """allow for branches section to rotate with respect to their parent with an angle greater than 90 degrees"""

    partition: AllenStructure = None

    @property
    def region_map(self):
        """
        Return RegionMap instance to manipulate the Allen mouse brain region hierarchy.

        :rtype: voxcell.region_map.RegionMap
        """
        return self.partition.region_map

    @property
    def annotations(self):
        """
        Return the mouse brain annotation voxel array

        :rtype: numpy.ndarray
        """
        return self.partition.annotations.raw

    @property
    def default_vector(self):
        return self.partition.mask_source.default_vector

    def orientation_field(self):
        """
        Return the brain orientation field pointing towards the outer shell of molecular layer

        :rtype: numpy.ndarray
        """
        loc_orient = np.asarray(self.orientations_field.load_object().raw, dtype=np.float32)
        loc_orient /= np.linalg.norm(loc_orient, axis=3)[..., np.newaxis]
        return loc_orient

    @pool_cache
    def _fixed_orientation(self, fixed_dim=None):
        loc_orient = np.copy(self.orientation_field())
        if fixed_dim is not None:
            loc_orient[..., fixed_dim] = 0.0
            loc_orient /= np.linalg.norm(loc_orient, axis=3)[..., np.newaxis]
        return loc_orient

    def fixed_dimension(self, branch_labels):
        if type(self.fixed_dimensions) == int:
            return self.fixed_dimensions
        for l in branch_labels:
            if l in self.fixed_dimensions:
                return self.fixed_dimensions[l]
        return -1

    def fix_orientation(self, branch_labels):
        """
        Get orientation field with fixed dimension.

        :param list[str] branch_labels: list of labels attached to the current segment.
        :rtype: numpy.ndarray
        """
        return self._fixed_orientation(self.fixed_dimension(branch_labels))

    @pool_cache
    def thicknesses(self):
        """
        Return the brain depth field, i.e. the distance of each voxel from its layer boundaries

        :rtype: numpy.ndarray
        """
        return np.asarray(
            self.depths_field.load_object().raw * self.depths_field.voxel_size, dtype=np.float32
        )

    @pool_cache
    def boundaries(self):
        """
        Return a boolean array which tells for each voxel transition (3,3,3), if it remains in the
        current orientation space.

        :rtype: numpy.ndarray
        """
        return self.boundaries_field.load_object().raw.reshape(self.annotations.shape + (3, 3, 3))

    def get_lay_abv(self, point):
        """
        Return the annotation layer abbreviation at the point location

        :param numpy.ndarray point: point location
        :return: layer abbreviation
        :rtype: str
        """
        return self.region_map.get(self.voxel_data_of(point, self.annotations), "acronym")[-2:]

    @staticmethod
    def is_within(vox, dataset):
        """
        Check if a voxel location is within a dataset's dimension, based on its shape.
        :param numpy.ndarray vox: 3D position of the voxel
        :param numpy.ndarray dataset: array to test
        :return: True if vox is within the dataset.
        :rtype: bool
        """
        return (
            len(dataset.shape) >= 3
            and len(vox) >= 3
            and np.all(vox >= 0)
            and (vox[0] < dataset.shape[0])
            and (vox[1] < dataset.shape[1])
            and (vox[2] < dataset.shape[2])
        )

    def voxel_data_of(self, point, dataset):
        """
        Retrieve voxel information from a dataset.
        :param numpy.ndarray point: floating point
        :param numpy.ndarray dataset: 3D numpy dataset
        :return: data stored at the point position.
        """
        loc_dataset = np.asarray(dataset)
        vox = self.partition.mask_source.voxel_of(point)
        if self.is_within(vox, loc_dataset):
            return loc_dataset[vox[0], vox[1], vox[2]]
        else:
            raise ValueError(
                f"Position is outside of the dataset.\n"
                f"Shape: {loc_dataset.shape}, Resolution: {self.partition.mask_source.voxel_size}."
            )

    def _ann_to_abv(self, id_reg):
        """
        Return layer abbreviation.
        Return None if region is not part of Cerebellar cortex

        :param int id_reg: region id
        :return: layer index in thickness and its abbreviation.
        :rtype: str | None
        """
        if id_reg is None or id_reg <= 0:
            return None
        return self.region_map.get(id_reg, "acronym")[-2:]

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
                if not self.boundaries()[last_vox[0], last_vox[1], last_vox[2]][
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

    def _test_new_rotation(self, rotation, source, target, branch_labels, last_voxel=None):
        """
        Check effect of rotation on target

        :param scipy.spatial.transform.Rotation rotation: rotation to apply
        :param numpy.ndarray source: source point
        :param numpy.ndarray target: target point
        :param numpy.ndarray last_voxel: last voxel
        """
        new_target = Rotation.from_euler("xyz", rotation).apply(target - source) + source
        new_voxel = self.partition.mask_source.voxel_of(new_target)
        return (
            np.all(new_voxel == last_voxel)
            or self.is_target_wrong(source, new_target, branch_labels),
            new_voxel,
        )

    def rotate_point(self, source, branch, i, branch_labels, old_rots):
        """
        Compute the rotation to apply at a source point so that it follows the change in the
        orientation field while making sure the target remains within the frontiers of the region.

        :param numpy.ndarray source: center of the rotation to apply
        :param bsb.morphologies.Branch branch: Branch on which the rotation should be applied
        :param int i: current index in branch
        :param list[str] branch_labels: list of labels attached to the current segment.
        :param RotationReminder old_rots: Previous rotations applied to the previous points.
        :return: Euler angle of the rotation to apply at the source point
        :rtype: scipy.spatial.transform.Rotation
        """
        target = branch.points[i]
        new_rotation = voxel_rotation_of(
            self.fix_orientation(branch_labels),
            self.partition.mask_source.voxel_of(source),
            self.default_vector,
        )
        diff_rotation = (new_rotation * old_rots.last_rotation.inv()).as_euler("xyz")
        diff_rotation[np.absolute(diff_rotation) < 1e-5] = 0
        inc = 1.0
        max_angle = np.pi / 2 if self.no_turn_back else np.pi
        if np.all(
            (
                self.partition.mask_source.voxel_of(target)
                == self.partition.mask_source.voxel_of(source)
            )
            * (diff_rotation == 0)
        ):
            # no changes of rotation and target lands in the same voxel
            # we assume the source was correct so unless the branch_labels changed
            scaled_diff_rotation = diff_rotation
        else:
            to_rotate = True
            old_voxel_pos = None
            old_voxel_neg = None
            skip_positive = False
            skip_negative = False
            while to_rotate:
                scaled_diff_rotation = diff_rotation * inc
                if skip_positive or (np.absolute(scaled_diff_rotation) > max_angle).any():
                    skip_positive = True
                else:
                    # test positive rotation
                    to_rotate, old_voxel_pos = self._test_new_rotation(
                        scaled_diff_rotation, source, target, branch_labels, old_voxel_pos
                    )
                if to_rotate and inc > 1:
                    scaled_diff_rotation = diff_rotation - scaled_diff_rotation
                    if skip_negative or (np.absolute(scaled_diff_rotation) > max_angle).any():
                        skip_negative = True
                    else:
                        # test negative rotation
                        to_rotate, old_voxel_neg = self._test_new_rotation(
                            scaled_diff_rotation, source, target, branch_labels, old_voxel_neg
                        )
                if to_rotate:
                    if np.linalg.norm(diff_rotation) == 0:
                        diff_rotation = old_rots.old_diff_rotation.as_euler("xyz")
                    elif skip_positive and skip_negative:
                        raise ValueError("Hit a wall. Stopping")
                    else:
                        inc += inc / 4
        scaled_diff_rotation = Rotation.from_euler("xyz", scaled_diff_rotation)
        if (
            np.linalg.norm(old_rots.rotation_to_correct.as_euler("xyz")) > 1e-5
            and inc == 1.0
            and not self.is_target_wrong(
                source,
                (scaled_diff_rotation * old_rots.rotation_to_correct).apply(target - source)
                + source,
                branch_labels,
            )
        ):
            scaled_diff_rotation = scaled_diff_rotation * old_rots.rotation_to_correct
        # Update previous rotations.
        old_rots.last_rotation = new_rotation
        if np.linalg.norm(diff_rotation) > 0:
            old_rots.old_diff_rotation = Rotation.from_euler("xyz", diff_rotation)
        new_to_correct = (
            old_rots.rotation_to_correct
            * scaled_diff_rotation.inv()
            * Rotation.from_euler("xyz", diff_rotation)
        )
        if (np.absolute(new_to_correct.as_euler("xyz")) > max_angle).any():
            # corrections pilled up to a point that we are going in the wrong direction aborting
            raise ValueError("Hit a wall. Stopping")
        old_rots.rotation_to_correct = new_to_correct
        return scaled_diff_rotation

    def process_scaling(self, point):
        """
        Calculate the local scaling factor to apply to a morphology segment based on the local
        thickness at its location.

        :param numpy.ndarray point: segment location
        :return: scaling factor at the current location
        :rtype: float
        """
        lay = self._ann_to_abv(int(self.voxel_data_of(point, self.annotations)))
        if lay is not None and lay in ["mo", "pu", "gr"]:
            if lay == "mo":  # take only mo thickness
                thick = np.sum(self.voxel_data_of(point, self.thicknesses())[0:2])
            elif lay == "pu":  # pu is too thin add mo
                thick = np.sum(self.voxel_data_of(point, self.thicknesses())[0:3])
            else:  # parts of pu is in gr
                thick = np.sum(self.voxel_data_of(point, self.thicknesses())[1:4])
            return np.maximum(thick / self.default_depth[lay], 0.1)
        else:  # out of the annotations / depth / orientations fields.
            return 0.1

    def _find_valid_scale(
        self, source, old_coord, max_scale, branch_labels, fixed_dimension, n_iter=8
    ):
        """Binary search for the maximum scale in [0, max_scale] where the segment
        source -> source + (old_coord - source) * scale doesn't cross a boundary."""
        direction = old_coord - source
        lo, hi = 0.0, max_scale
        best = 0.0
        best_coord = np.copy(source)
        for _ in range(n_iter):
            mid = (lo + hi) / 2
            new_coord = source + direction * mid
            if 0 <= fixed_dimension <= 2:
                new_coord[fixed_dimension] = old_coord[fixed_dimension]
            if self.is_target_wrong(source, new_coord, branch_labels):
                hi = mid
            else:
                best = mid
                best_coord = np.copy(new_coord)
                lo = mid
        return best, best_coord

    def scale_morpho(self, branch, i, scaling, branch_labels):
        """
        Scale a morphology's branch at a position i, if the new position remains within
        the frontiers of the region.

        :param bsb.morphologies.Branch branch: morphology's branch to scale
        :param int i: position in the branch
        :param float scaling: factor of scaling.
        :param list[str] branch_labels: list of labels attached to the current segment.
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
        fixed_dimension = self.fixed_dimension(branch_labels)
        if 0 <= fixed_dimension <= 2:
            new_coord[fixed_dimension] = old_coord[fixed_dimension]
        # check that every voxel between the points remain within the region boundary
        rescale = self.is_target_wrong(branch.points[i - 1], new_coord, branch_labels)
        # if scaling resulted in an overshoot
        if rescale:
            # Binary search: keep the point alive at the furthest valid position
            best_scale, new_coord = self._find_valid_scale(
                branch.points[i - 1], old_coord, scaling, branch_labels, fixed_dimension
            )
            if best_scale > 1e-3:  # threshold: avoid degenerate near-zero segments
                scaling = best_scale
                rescale = False
        # set to old point coordinate if rescale failed
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
                branch_labels = get_branch_labels(branch, 0)
                rotation = voxel_rotation_of(
                    self.fix_orientation(branch_labels),
                    self.partition.mask_source.voxel_of(branch.points[0]),
                    self.default_vector,
                )
                branch.root_rotate(rotation)
                curr_scaling = self.process_scaling(branch.points[0])
                axis_max = max(
                    enumerate(
                        np.absolute(
                            self.voxel_data_of(
                                branch.points[0], self.fix_orientation(branch_labels)
                            )
                        )
                    ),
                    key=lambda x: x[1],
                )[0]
                old_diff_rotation = np.zeros(3)
                old_diff_rotation[axis_max] = 1e-2
                stack_data.append(
                    (
                        branch,
                        RotationReminder(rotation, Rotation.from_euler("xyz", old_diff_rotation)),
                        curr_scaling,
                    )
                )
            except ValueError as _:
                continue
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
                    branch_labels = get_branch_labels(branch, i)
                    has_scale = False
                    if np.isin(list(branch_labels), self.rescale).any() and i > 0:
                        has_scale = True
                        fail_rescale, curr_scaling = self.scale_morpho(
                            branch, i, curr_scaling, branch_labels
                        )
                        if fail_rescale:
                            self.delete_point(branch, i)
                            continue
                    if np.isin(list(branch_labels), self.deform).any():
                        try:
                            rotation = self.rotate_point(
                                last_point, branch, i, branch_labels, old_rots
                            )
                            branch.root_rotate(rotation, downstream_of=last_index)
                            old_rots.original_rotation = old_rots.original_rotation * rotation
                        except ValueError as _:
                            self.delete_point(branch, i)
                            continue
                    elif (
                        i > 0
                        and not has_scale
                        and self.is_target_wrong(
                            last_point,
                            branch.points[i],
                            branch_labels,
                        )
                    ):
                        # case when the next branch point was neither rotated nor scaled
                        # hence it could land outside the region.
                        # we delete the current point so that the last one land inside
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
                                if self.is_within(
                                    self.partition.mask_source.voxel_of(child.points[0]),
                                    self.annotations,
                                )
                            ]
                        )
        morphology.close_gaps()
        morphology.optimize(force=True)
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
            if self.is_within(uniq_vox, self._fixed_orientation()):
                translation_vec = (uniq_vox + 0.5) * self.partition.voxel_size
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
