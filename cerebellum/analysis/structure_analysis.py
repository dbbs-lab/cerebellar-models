from typing import List, Tuple

import numpy as np
from bsb import CellType, ConnectivitySet, Scaffold, warn
from matplotlib import pyplot as plt

from cerebellum.analysis.plots import Legend, ScaffoldPlot
from cerebellum.analysis.report import LIST_CT_INFO, CellTypeInfo, Report


class TablePlot:
    values = []
    rows = []
    columns = []

    def update_values(self):
        self.values = []
        self.rows = []

    def plot_table(self, **kwargs):
        if len(self.values) == 0:
            warn("No values to plot")
            return
        dict_plot = dict(
            rowColours=np.full((len(self.rows), 3), [0.8, 0.8, 0.8]),
            rowLoc="right",
            colColours=np.full((len(self.columns), 3), [0.8, 0.8, 0.8]),
            loc="center",
        )
        dict_plot.update(kwargs)
        self.get_ax().table(
            cellText=self.values, rowLabels=self.rows, colLabels=self.columns, **dict_plot
        )


class PlacementTable(TablePlot, ScaffoldPlot):
    def __init__(
        self,
        fig_size: Tuple[float, float],
        scaffold: Scaffold = None,
        dict_colors=None,
        dict_abv=None,
        **kwargs,
    ):
        super().__init__(
            fig_size,
            scaffold=scaffold,
            dict_colors=dict_colors,
            **kwargs,
        )
        self.columns = ["Cell counts", "Cell densities [$\mu m^{-3}$]"]
        self.dict_abv = dict_abv or {ct.name: ct.abbreviation for ct in LIST_CT_INFO}

    def extract_ct_name(self, ct: CellType):
        key = ct.name.split("_cell")[0]
        return self.dict_abv[key] if key in self.dict_abv else ct.name

    def update(self):
        super().update()
        self.update_values()
        for i, ps in enumerate(self.scaffold.get_placement_sets()):
            ct = ps.cell_type
            self.rows.append(self.extract_ct_name(ct))
            counts = ps.load_positions().shape[0]
            volume = [p.volume() for place in ct.get_placement() for p in place.partitions]
            self.values.append(["{:.2E}".format(counts), "{:.2E}".format(counts / np.sum(volume))])

    def plot(self, **kwargs):
        super().plot()
        self.plot_table(**kwargs)


class ConnectivityTable(TablePlot, ScaffoldPlot):
    def __init__(
        self,
        fig_size: Tuple[float, float],
        scaffold: Scaffold = None,
        dict_colors=None,
        dict_abv=None,
        **kwargs,
    ):
        super().__init__(
            fig_size,
            scaffold=scaffold,
            dict_colors=dict_colors,
            **kwargs,
        )
        self.columns = ["Nb. Synapses", "Synapses per pair", "Convergence", "Divergence"]
        self.dict_abv = dict_abv or {ct.name: ct.abbreviation for ct in LIST_CT_INFO}
        self.update_values()

    def extract_strat_name(self, ps: ConnectivitySet):
        splits = ps.tag.split("_")
        tag = []
        found = False
        to_convert = []
        for text in splits:
            if "to" == text:
                tag.extend(to_convert)
                to_convert = []
                if not found:
                    tag.append(text)
                    found = True
                else:
                    tag = [tag[-1], "to"]
            else:
                if "cell" == text:
                    continue
                to_convert.append(text)
                to_convert = to_convert[-2:]
                if "_".join(to_convert) in self.dict_abv:
                    tag.append(self.dict_abv["_".join(to_convert)])
                    to_convert = []
        if len(to_convert):
            tag.extend(to_convert)
        return " ".join(tag)

    def update(self):
        super().update()
        self.update_values()
        for ps in self.scaffold.get_connectivity_sets():
            # Get the ConnectivityIterator for the current connectivity strategy
            cs = self.scaffold.get_connectivity_set(ps.tag).load_connections().as_globals()
            pre_locs, post_locs = cs.all()
            # Find the pairs of pre-post neurons (combos)
            # and count how many synapses there are between each pair (combo_counts)
            combos, combo_counts = np.unique(
                np.column_stack((pre_locs[:, 0], post_locs[:, 0])), axis=0, return_counts=True
            )

            # Find the unique post and pre neurons
            uniquePre, uniquePre_count = np.unique(combos[:, 0], axis=0, return_counts=True)
            niquePost, uniquePost_count = np.unique(combos[:, 1], axis=0, return_counts=True)
            self.rows.append(self.extract_strat_name(ps))
            self.values.append(
                [
                    len(pre_locs),
                    "{:.2} $\pm$ {:.2}".format(np.mean(combo_counts), np.std(combo_counts)),
                    "{:.2} $\pm$ {:.2}".format(np.mean(uniquePost_count), np.std(uniquePost_count)),
                    "{:.2} $\pm$ {:.2}".format(np.mean(uniquePre_count), np.std(uniquePre_count)),
                ]
            )

    def plot(self, **kwargs):
        super().plot()
        self.plot_table(**kwargs)


class CellPlacement3D(ScaffoldPlot):
    def __init__(
        self,
        fig_size: Tuple[float, float],
        scaffold: Scaffold = None,
        dict_colors=None,
        ignored_ct=None,
        **kwargs,
    ):
        super().__init__(fig_size, scaffold, dict_colors=dict_colors, **kwargs)
        self.ignored_ct = ignored_ct or ["mossy_fibers", "glomerulus", "ubc_glomerulus"]

    def init_plot(self, **kwargs):
        self.is_plotted = False
        fig = plt.figure(figsize=self.fig_size, **kwargs)
        ax = fig.add_subplot(111, projection="3d")
        return fig, ax

    @staticmethod
    def set_axes_equal(ax):
        """
        Make axes of 3D plot have equal scale so that spheres appear as spheres,
        cubes as cubes, etc.

        Input
          ax: a matplotlib axis, e.g., as output from plt.gca().
        """

        x_limits = ax.get_xlim3d()
        y_limits = ax.get_ylim3d()
        z_limits = ax.get_zlim3d()

        x_range = abs(x_limits[1] - x_limits[0])
        x_middle = np.mean(x_limits)
        y_range = abs(y_limits[1] - y_limits[0])
        y_middle = np.mean(y_limits)
        z_range = abs(z_limits[1] - z_limits[0])
        z_middle = np.mean(z_limits)

        # The plot bounding box is a sphere in the sense of the infinity
        # norm, hence I call half the max range the plot radius.
        plot_radius = 0.5 * max([x_range, y_range, z_range])

        ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
        ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
        ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

    def plot(self, **kwargs):
        super().plot()
        ax = self.get_ax()
        for i, ps in enumerate(self.scaffold.get_placement_sets()):
            ct = ps.cell_type
            ct_name = ct.name.split("_cell")[0]
            if not ct.name in self.ignored_ct and ct_name in self.dict_colors:
                *color, alpha = self.dict_colors[ct_name]
                scale = np.power(ct.spatial.radius, 2)
                positions = ps.load_positions()
                ax.scatter(
                    positions[:, 0],
                    positions[:, 1],
                    positions[:, 2],
                    c=np.repeat([color], len(positions), axis=0),
                    alpha=np.repeat([alpha], len(positions), axis=0),
                    s=scale,
                )
        ax.set_xlabel("x in $\mu m$")
        ax.set_ylabel("y in $\mu m$")
        ax.set_zlabel("z in $\mu m$")
        self.set_axes_equal(ax)
        ax.set_title("Placement results", fontsize=40)


class StructureReport(Report):
    def __init__(self, pathname: str, cell_type_info: List[CellTypeInfo] = None):
        super().__init__(pathname, cell_type_info)
        legend = Legend(
            (10, 2.5),
            3,
            dict_legend=dict(columnspacing=2.0, handletextpad=0.1, fontsize=20, loc="lower center"),
        )
        density_table = PlacementTable((5, 2.5), scaffold=self.network, dict_abv=self.abbreviations)
        connectivity_table = ConnectivityTable(
            (10, 5), scaffold=self.network, dict_abv=self.abbreviations
        )
        plot3d = CellPlacement3D((10, 10), scaffold=self.network)
        self.add_plot("density_table", density_table)
        self.add_plot("connectivity_table", connectivity_table)
        self.add_plot("placement_3d", plot3d)
        self.add_plot("legend", legend)
        legend.set_axis_off()
        density_table.set_axis_off()
        connectivity_table.set_axis_off()
