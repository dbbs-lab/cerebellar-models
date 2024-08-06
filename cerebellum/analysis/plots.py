"""
    Module for abstract classes interfacing matplotlib plots.
"""

import abc
from typing import List, Tuple

import numpy as np
from bsb import Scaffold
from matplotlib import patches
from matplotlib import pyplot as plt


class Plot(abc.ABC):
    """
    Matplotlib plot class interface.
    """

    def __init__(
        self,
        fig_size: Tuple[float, float] | np.ndarray,
        nb_rows: int = 1,
        nb_cols: int = 1,
        dict_colors=None,
        **kwargs,
    ):
        if nb_rows <= 0 or nb_cols <= 0:
            raise ValueError("Number of columns and rows must be strictly positive.")
        self.nb_rows = nb_rows
        """Number of sub-panel rows in the plot."""
        self.nb_cols = nb_cols
        """Number of sub-panel columns in the plot."""
        self.fig_size = fig_size
        """Tuple size of the figure in inches."""
        self.is_plotted = False
        """Flag to indicate if this plot has been plotted."""
        self.is_updated = False
        """Flag to indicate if this plot has been updated."""
        self.figure = None
        """Matplotlib Figure of the plot."""
        self.axes = None
        """Matplotlib Axes of the plot"""
        self.dict_colors = {}
        """Dictionary of element name to their color"""
        if dict_colors is not None:
            for key, value in dict_colors.items():
                self.set_color(key, value)
        self.init_plot(**kwargs)

    def clear(self):
        """
        Clear the figure axes
        """
        for ax in self.get_axes():
            ax.clear()
        self.is_plotted = False

    def init_plot(self, **kwargs):
        """
        Initialize the plot and return figure and axes.
        """
        self.clear()
        self.figure, self.axes = plt.subplots(
            nrows=self.nb_rows, ncols=self.nb_cols, figsize=self.fig_size, **kwargs
        )

    def get_axes(self):
        """
        Return figure axes as a flat list.
        """
        if self.axes is None:
            return []
        elif self.nb_cols == 1 and self.nb_rows == 1:
            return [self.axes]
        elif self.nb_cols == 1 or self.nb_rows == 1:
            return self.axes
        else:
            return [ax for col_ax in self.axes for ax in col_ax]

    def get_ax(self, idx: int = 0):
        """
        Return the axis at the given index.
        """
        if idx >= self.nb_cols * self.nb_rows or idx < 0:
            raise IndexError(
                f"Index of matplotlib ax out of range. Max {self.nb_cols * self.nb_rows}"
            )
        if self.nb_cols == 1 and self.nb_rows == 1:
            return self.axes
        elif self.nb_cols == 1 or self.nb_rows == 1:
            return self.axes[idx]
        else:
            return self.axes[idx // self.nb_cols][idx % self.nb_cols]

    def set_color(self, key: str, color: np.ndarray[float]):
        """
        Set the color dictionary for a given key.
        Colors must be an array of type RGB or RGBA.
        """
        if len(color) != 3 and len(color) != 4:
            raise ValueError("Color must be an array of size 3 or 4.")
        self.dict_colors[key] = np.array(color, dtype=float)
        # colors are updated so the figure should be updated too.
        if self.is_plotted:
            self.clear()

    def save_figure(self, output_name: str, dpi: int, pad=0, **kwargs):
        """
        Save the figure as a file.
        """
        if not self.is_plotted:
            self.plot()
        self.figure.tight_layout(pad=pad)
        self.figure.savefig(output_name, dpi=dpi, facecolor="white", **kwargs)

    def update(self):
        """
        Update function to prepare the data before plotting.
        """
        self.is_updated = True

    def plot(self, *args, **kwargs):
        """
        Plot or replot the figure.
        Calls the update function if needed.
        """
        if self.is_plotted:
            self.clear()
        if not self.is_updated:
            self.update()
        self.is_plotted = True

    def show(self):  # pragma:nocover
        """
        Show the figure.
        The figure will be plotted if needed.
        """
        if not self.is_plotted:
            self.plot()
        self.figure.show()

    def set_axis_off(self, axes=None):
        """
        Removes the borders of the provided axes.
        If none are provided, all axes' borders will be removed.

        :param axes: List of matplotlib axes
        """
        if axes is None:
            axes = self.get_axes()
        for ax in axes:
            ax.axis("off")


class ScaffoldPlot(Plot):
    """
    Matplotlib plot interface for BSB Scaffold.
    """

    def __init__(
        self,
        fig_size: Tuple[float, float],
        scaffold: Scaffold,
        nb_rows: int = 1,
        nb_cols: int = 1,
        dict_colors: dict = None,
        **kwargs,
    ):
        super().__init__(fig_size, nb_rows, nb_cols, dict_colors, **kwargs)
        self.scaffold = scaffold

    def set_scaffold(self, scaffold):
        """
        Set the scaffold of the plot.
        """
        is_different = scaffold != self.scaffold
        if is_different:
            self.scaffold = scaffold
            self.is_updated = False
            if self.is_plotted:
                self.clear()
        return is_different


class Legend(Plot):
    """
    Patch Legend plot
    """

    def __init__(
        self,
        fig_size: Tuple[float, float],
        nb_cols: int = 1,
        dict_colors: dict = None,
        dict_legend=None,
        **kwargs,
    ):
        super().__init__(fig_size, 1, 1, dict_colors, **kwargs)
        self.dict_legend = dict_legend or dict()
        self.cols_legend = nb_cols

    def plot(self, **kwargs):
        super().plot()
        patchs = [
            patches.Patch(color=color[:3], label=label) for label, color in self.dict_colors.items()
        ]
        dict_plot = self.dict_legend.copy()
        dict_plot.update(kwargs)
        self.get_ax().legend(patchs, self.dict_colors.keys(), ncol=self.cols_legend, **dict_plot)

    def remove_ct(self, to_keep: List[str], to_ignore: List[str] = None):
        """
        Remove cell types that are not in the list to keep or in the list to ignore

        :params to_keep List[str]: list of cell types to keep
        :params to_ignore List[str]: list of cell types to remove

        """
        to_ignore = to_ignore or []
        self.dict_colors = {
            ct: color
            for ct, color in self.dict_colors.items()
            if ct not in to_ignore and (ct in to_keep or ct + "_cell" in to_keep)
        }
