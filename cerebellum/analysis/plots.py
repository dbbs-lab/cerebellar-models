import abc
from typing import Tuple

import numpy as np
from bsb import Scaffold
from matplotlib import patches
from matplotlib import pyplot as plt
from matplotlib.axes import Axes


class Plot(abc.ABC):
    def __init__(
        self,
        fig_size: Tuple[float, float],
        nb_rows: int = 1,
        nb_cols: int = 1,
        dict_colors=None,
        **kwargs,
    ):
        self.nb_rows = nb_rows
        self.nb_cols = nb_cols
        self.fig_size = fig_size
        self.is_plotted = False
        self.figure, self.axes = self.init_plot(**kwargs)

        self.dict_colors = dict_colors if dict_colors is not None else {}

    def init_plot(self, **kwargs):
        self.is_plotted = False
        return plt.subplots(nrows=self.nb_rows, ncols=self.nb_cols, figsize=self.fig_size, **kwargs)

    @classmethod
    def as_subfig(cls, figure: plt.Figure, axes: plt.Axes, dict_colors=None):
        result = cls((1, 1), dict_colors=dict_colors)
        plt.close(result.figure)
        result.figure = figure
        result.axes = axes

    def get_axes(self):
        if self.nb_cols == 1 and self.nb_rows == 1:
            return [self.axes]
        elif self.nb_cols == 1 or self.nb_rows == 1:
            return self.axes
        else:
            return [ax for col_ax in self.axes for ax in col_ax]

    def get_ax(self, idx: int = 0) -> Axes:
        if self.nb_cols == 1 and self.nb_rows == 1:
            return self.axes
        elif self.nb_cols == 1 or self.nb_rows == 1:
            return self.axes[idx]
        else:
            return self.axes[idx // self.nb_cols][idx % self.nb_cols]

    def clear(self):
        for ax in self.get_axes():
            ax.clear()
        self.is_plotted = False

    def set_color(self, key: str, color: np.ndarray[float]):
        assert len(color) == 3 or len(color) == 4
        self.dict_colors[key] = color
        # colors are updated so the figure should be updated too.
        if self.is_plotted:
            self.clear()

    def save_figure(self, output_name: str, dpi: int, pad=0, **kwargs):
        if not self.is_plotted:
            self.plot()
        self.figure.tight_layout(pad=pad)
        self.figure.savefig(output_name, dpi=dpi, facecolor="white", **kwargs)

    def plot(self, *args, **kwargs):
        self.is_plotted = True

    def show(self):
        if not self.is_plotted:
            self.plot()
        self.figure.show()

    def set_axis_off(self, axes=None):
        if axes is None:
            axes = self.get_axes()
        for ax in axes:
            ax.axis("off")

    def add_legend(self, id_ax: int, elem: dict, **kwargs):
        ax = self.get_ax(id_ax)
        ax.legend(elem.values(), elem.keys(), **kwargs)


class ScaffoldPlot(Plot):
    def __init__(
        self,
        fig_size: Tuple[float, float],
        scaffold: Scaffold,
        nb_rows: int = 1,
        nb_cols: int = 1,
        dict_colors=None,
        **kwargs,
    ):
        super().__init__(fig_size, nb_rows, nb_cols, dict_colors, **kwargs)
        self.scaffold = scaffold

    def set_scaffold(self, scaffold):
        self.scaffold = scaffold
        if self.is_plotted:
            self.clear()


class Legend(Plot):
    def __init__(self, fig_size: Tuple[float, float], nb_cols: int = 1, dict_legend=None, **kwargs):
        super().__init__(fig_size, 1, 1, **kwargs)
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
