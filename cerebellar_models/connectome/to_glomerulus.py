"""
Module for the configuration node of every presynaptic cell to Glomerulus ConnectionStrategy
"""

import abc
import itertools

import numpy as np
from bsb import Chunk, ConnectionStrategy, InvertedRoI, config, warn
from bsb.connectivity.strategy import HemitypeCollection, roi_key
from scipy.spatial import cKDTree

from cerebellar_models.connectome.presyn_dist_strat import PresynDistStrat


def norm_exp_dist(rng: np.random.Generator, size: int = 1, b: float = 2.0):
    """
    Normalized exponential random generator for distance based selection

    :param numpy.random.Generator rng: generator to draw from.
    :param int size: number of random to sample
    :param float b: strength of the exponential decay.
    :return: random numbers sampled
    :rtype: numpy.ndarray
    """
    return (np.exp(-b * rng.random(size)) - np.exp(-b)) / (1 - np.exp(-b))


def reserve_closest_targets(presyn_pos: np.ndarray, postsyn_pos: np.ndarray):
    """
    Pair each presynaptic cell with a distinct postsynaptic cell, as close as possible.

    Candidate pairs are greedily assigned by increasing distance, ties being broken by the
    presynaptic and then the postsynaptic index, so the result only depends on the positions.
    If there are more presynaptic than postsynaptic cells, the remaining presynaptic cells
    are left without a reserved target.

    :param numpy.ndarray presyn_pos: positions of the presynaptic cells.
    :param numpy.ndarray postsyn_pos: positions of the postsynaptic cells.
    :return: for each presynaptic cell, the index of its reserved postsynaptic cell, or -1.
    :rtype: numpy.ndarray
    """
    n_pre, n_post = len(presyn_pos), len(postsyn_pos)
    reserved = np.full(n_pre, -1, dtype=int)
    if n_pre == 0 or n_post == 0:
        return reserved
    tree = cKDTree(postsyn_pos)
    taken = np.zeros(n_post, dtype=bool)
    pending = np.arange(n_pre)
    k = min(8, n_post)
    while len(pending) > 0 and not np.all(taken):
        # Look for the k closest targets of each pending presynaptic cell.
        dists, targets = tree.query(presyn_pos[pending], k=k)
        dists = np.reshape(dists, (len(pending), k)).ravel()
        targets = np.reshape(targets, (len(pending), k)).ravel()
        pre_ids = np.repeat(pending, k)
        for i in np.lexsort((targets, pre_ids, dists)):
            if reserved[pre_ids[i]] < 0 and not taken[targets[i]]:
                reserved[pre_ids[i]] = targets[i]
                taken[targets[i]] = True
        pending = pending[reserved[pending] < 0]
        if k == n_post:
            break
        k = min(2 * k, n_post)
    return reserved


class ConnectomeGlomerulus(InvertedRoI, ConnectionStrategy):
    """
    BSB Connection strategy to connect a presynaptic cell to Glomeruli.

    Each glomerulus is connected to exactly one presynaptic cell. Each presynaptic cell first
    gets its closest glomerulus reserved (see :func:`reserve_closest_targets`), which guarantees
    that every presynaptic cell connects to at least one glomerulus, as long as there are enough
    glomeruli. The remaining glomeruli draw their presynaptic cell among the ones selected by
    :meth:`pre_selection`.
    """

    def _reserved_presyn(self, presyn_pos, all_post_ps, post_ps):
        """
        Compute, for the glomeruli of ``post_ps``, the presynaptic cell they are reserved for.

        The reservation is computed over the whole presynaptic and postsynaptic populations, so
        every job, whichever its chunk, draws the same reservation.

        :param numpy.ndarray presyn_pos: positions of the whole presynaptic population.
        :param bsb.storage.interfaces.PlacementSet all_post_ps: whole postsynaptic population.
        :param bsb.storage.interfaces.PlacementSet post_ps: postsynaptic cells of this job.
        :return: for each glomerulus of ``post_ps``, the id of its reserved presynaptic cell
            in ``presyn_pos``, or -1.
        :rtype: numpy.ndarray
        """
        all_post_pos = all_post_ps.load_positions()
        if len(presyn_pos) > len(all_post_pos):
            warn(
                f"'{self.name}': {len(presyn_pos)} presynaptic cells for only "
                f"{len(all_post_pos)} glomeruli, some presynaptic cells will not be connected."
            )
        reserved = reserve_closest_targets(presyn_pos, all_post_pos)
        owners = np.full(len(all_post_pos), -1, dtype=int)
        has_target = reserved >= 0
        owners[reserved[has_target]] = np.nonzero(has_target)[0]
        # Global ids are sorted, so locate the glomeruli of this job in the whole population.
        return owners[np.searchsorted(all_post_ps.load_ids(), post_ps.load_ids())]

    def connect(self, pre, post):
        # The whole presynaptic population is loaded, whatever the region of interest, so that
        # each job computes the same reservation of glomeruli.
        all_pre_sets = HemitypeCollection(pre.hemitype, None).placement
        all_post_sets = HemitypeCollection(post.hemitype, None).placement
        # We use a truncated exponential distribution to favour the presynaptic fibers closer the
        # postsynaptic glomerulus.
        for pre_ps in all_pre_sets:
            for post_ps, all_post_ps in zip(post.placement, all_post_sets):
                presyn_pos = pre_ps.load_positions()
                glomeruli_pos = post_ps.load_positions()
                n_glom = len(glomeruli_pos)
                reserved = self._reserved_presyn(presyn_pos, all_post_ps, post_ps)

                # We work assuming that there is at least 1 presynaptic fiber in the ROI.
                # Otherwise, there is something wrong in the placement phase.
                pre_locs = np.full((n_glom, 3), -1, dtype=int)
                post_locs = np.full((n_glom, 3), -1, dtype=int)

                # Keyed on the strategy, the cell type pair and the chunks involved, so
                # every rank draws the same connections for the same chunk pair.
                rng = self.get_rng(
                    key=(
                        "connectivity",
                        self.name,
                        pre_ps.cell_type.name,
                        post_ps.cell_type.name,
                        roi_key(pre),
                        roi_key(post),
                    ),
                )
                # We connect each glomerulus to a presynaptic cell.
                for j, glomerulus in enumerate(glomeruli_pos):
                    if reserved[j] >= 0:
                        pre_locs[j, 0] = reserved[j]
                        post_locs[j, 0] = j
                        continue
                    pre_ids = self.pre_selection(presyn_pos, glomerulus)
                    if len(pre_ids) > 0:
                        pre_locs[j, 0] = pre_ids[
                            int(np.floor(len(pre_ids) * norm_exp_dist(rng)[0]))
                        ]
                    else:
                        # if there is no presyn within the box use the closest one
                        pre_locs[j, 0] = np.argmin(np.linalg.norm(presyn_pos - glomerulus, axis=1))
                    post_locs[j, 0] = j
                self.connect_cells(pre_ps, post_ps, pre_locs, post_locs)

    @abc.abstractmethod
    def pre_selection(
        self,
        presyn_pos,
        glom_pos,
    ):  # pragma: nocover
        """
        Order presynaptic cell ids based on their respective distance to glomerulus

        :param numpy.ndarray presyn_pos: list of presynaptic cell positions
        :param numpy.ndarray glom_pos: single glomerulus position
        :return: presynaptic cell ids sorted by distance to glomerulus
        :rtype: numpy.ndarray
        """
        pass


@config.node
class ConnectomeMossyGlomerulus(ConnectomeGlomerulus):
    """
    BSB Connection strategy to connect Mossy fibers to Glomeruli.
    """

    x_length: float = config.attr(type=float, required=True)
    """Length of the box along the x axis surrounding the glomerulus cell soma in which the 
        presynaptic cell can be connected."""
    y_length: float = config.attr(type=float, required=True)
    """Length of the box along the y axis surrounding the glomerulus cell soma in which the 
        presynaptic cell can be connected."""

    def pre_selection(
        self,
        presyn_pos,
        glom_pos,
    ):
        diff = np.absolute(glom_pos - presyn_pos)
        ids_to_keep = np.where((diff[:, 0] <= self.x_length) * (diff[:, 1] <= self.y_length))[0]
        dist = np.linalg.norm(diff[ids_to_keep], axis=1)
        return ids_to_keep[np.argsort(dist)]


@config.node
class ConnectomeUBCGlomerulus(ConnectomeGlomerulus, PresynDistStrat):
    """
    BSB Connection strategy to connect UBC to Glomeruli.
    """

    def pre_selection(
        self,
        presyn_pos,
        glom_pos,
    ):
        dist = np.linalg.norm(presyn_pos - glom_pos, axis=1)
        ids_to_keep = np.where(dist <= self.radius)[0]
        return ids_to_keep[np.argsort(dist[ids_to_keep])]
