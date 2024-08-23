import numpy as np
from bsb._util import rotation_matrix_from_vectors
from scipy.spatial.transform import Rotation


def signed_modulo(a, mod):
    """
    Returns the signed modulo of an array of values a modulo mod, ie between -mod/2 and mod/2

    :param numpy.ndarray a: array of values
    :param int | float | double mod: modulo value
    :return: a modulo mod
    :rtype: numpy.ndarray
    """
    n = np.absolute(a) % mod
    s = np.sign(a)
    filt = n > mod / 2
    b = n * s
    b[filt] = -s[filt] * (mod - n[filt])
    return b


def _line_point(src, tgt, deriv, di):
    source = np.copy(src)
    target = np.copy(tgt)
    derivate = np.copy(deriv)
    diff = np.copy(di)
    result = [np.copy(source)]
    p1 = 2 * derivate[1] - derivate[0]
    p2 = 2 * derivate[2] - derivate[0]
    while source[0] != target[0]:
        source[0] += diff[0]
        if p1 >= 0:
            source[1] += diff[1]
            p1 -= 2 * derivate[0]
        if p2 >= 0:
            source[2] += diff[2]
            p2 -= 2 * derivate[0]
        p1 += 2 * derivate[1]
        p2 += 2 * derivate[2]
        result.append(np.copy(source))
    return np.array(result)


def bresenham_line(source, target):
    """
    Draw a lign between a source and a target position, and returns all the crossing
    3D integer position. This follows the Bresenham algorithm.

    :param numpy.ndarray source: 3D position of the starting point
    :param numpy.ndarray target: 3D position of the ending point
    :return: list of 3D integer position between source and target
    :rtype: np.ndarray
    """

    derivate = np.abs(target - source)
    diff = np.ones(3, dtype=int)
    diff[target <= source] = -1

    # Driving axis is X-axis
    if derivate[0] >= derivate[1] and derivate[0] >= derivate[2]:
        return _line_point(source, target, derivate, diff)
    # Driving axis is Y-axis
    elif derivate[1] >= derivate[0] and derivate[1] >= derivate[2]:
        return np.roll(
            _line_point(
                np.roll(source, -1), np.roll(target, -1), np.roll(derivate, -1), np.roll(diff, -1)
            ),
            1,
            axis=1,
        )
    # Driving axis is Z-axis
    else:
        return np.roll(
            _line_point(
                np.roll(source, -2), np.roll(target, -2), np.roll(derivate, -2), np.roll(diff, -2)
            ),
            2,
            axis=1,
        )


def boundaries_index_of(vox, new_vox):
    """
    Convert the difference between two voxels into a transition in a (3,3) matrix.

    :param numpy.ndarray vox: source voxel
    :param numpy.ndarray new_vox: target voxel
    :return: transition indexes from source to target
    :rtype: Tuple(int, int, int)
    """
    return tuple(np.minimum(np.maximum((new_vox - vox), -1), 1) + 1)


def get_diff_angle(orientation_field, vox, old_vox):
    """
    Get the Euler angle between the orientation vectors of two different position.

    :param numpy.ndarray orientation_field: brain orientation field
    :param numpy.ndarray vox: first voxel position
    :param numpy.ndarray old_vox: second voxel position
    :return: Euler angle between the two orientation field vectors.
    :rtype: numpy.ndarray
    """
    loc_rot = Rotation.from_matrix(
        rotation_matrix_from_vectors(
            orientation_field[old_vox[0], old_vox[1], old_vox[2]],
            orientation_field[vox[0], vox[1], vox[2]],
        )
    )
    return loc_rot.as_euler("xyz")
