"""
Module for the plots and reports related to the structural analysis of BSB scaffold.
"""

from typing import List, Tuple, Union

import numpy as np
from bsb import AfterConnectivityHook, CellType, ConnectivitySet, Scaffold, config, warn
from matplotlib import pyplot as plt

from cerebellum.analysis.plots import Legend, ScaffoldPlot
from cerebellum.analysis.report import LIST_CT_INFO, BSBReport, PlotTypeInfo


class TablePlot:
    """
    Mixin for plotting tables with matplotlib
    """

    _values = []
    """List of list of values from the table"""
    table_values = []
    """List of list of string values to put in the table"""
    rows = []
    """Names of the table's rows"""
    columns = []
    """Names of the table's columns"""

    def reset_table(self):
        """
        Update the values of the table.
        """
        self._values = []
        self.table_values = []
        self.rows = []

    def plot_table(self, **kwargs):
        """
        Plot the table in the Figure.
        """
        if len(self.table_values) == 0:
            warn("No values to plot", UserWarning)
            return
        dict_plot = dict(
            rowColours=np.full((len(self.rows), 3), [0.8, 0.8, 0.8]),
            rowLoc="right",
            colColours=np.full((len(self.columns), 3), [0.8, 0.8, 0.8]),
            loc="center",
        )
        dict_plot.update(kwargs)
        self.get_ax().table(
            cellText=self.table_values, rowLabels=self.rows, colLabels=self.columns, **dict_plot
        )


class PlacementTable(TablePlot, ScaffoldPlot):
    """
    Table plot of the results of the placement for BSB Scaffold.
    This includes the counts and density of each cell type.
    """

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
        self.dict_abv = dict_abv or {}
        """Dictionary of abbreviations for cell types"""

    def extract_ct_name(self, ct: CellType):
        """
        Convert the name of a cell type to its abbreviation.
        """
        return self.dict_abv[ct.name] if ct.name in self.dict_abv else ct.name

    def update(self):
        super().update()
        self.reset_table()
        for i, ps in enumerate(self.scaffold.get_placement_sets()):
            ct = ps.cell_type
            self.rows.append(self.extract_ct_name(ct))
            counts = ps.load_positions().shape[0]
            volume = np.sum([p.volume() for place in ct.get_placement() for p in place.partitions])
            self._values.append([counts, volume])
            self.table_values.append(
                [
                    "{:.2E}".format(counts),
                    "{:.2E}".format((counts / volume) if volume > 0.0 else 0.0),
                ]
            )

    def plot(self, **kwargs):
        super().plot()
        self.plot_table(**kwargs)

    def get_volumes(self):
        """
        Return a dictionary which gives for each cell type
        the volume occupied by its cells in the Scaffold.
        The plot needs to be updated.

        :rtype: Dict[str, int]
        """
        return {ct: line[1] for ct, line in zip(self.rows, self._values)}

    def get_counts(self):
        """
        Return a dictionary which gives for each cell type
        the number placed in the Scaffold.
        The plot needs to be updated.

        :rtype: Dict[str, int]
        """
        return {ct: line[0] for ct, line in zip(self.rows, self._values)}

    def get_densities(self):
        return {ct: line[0] / line[1] for ct, line in zip(self.rows, self._values)}


class ConnectivityTable(TablePlot, ScaffoldPlot):
    """
    Table plot of the results of the connectivity for BSB Scaffold.
    This includes for each pair of connected cell types:

    - the number of synapses formed
    - the number of synapses per unique pair of cell.
    - the convergence ratio defined as the mean number of afferent
      connections created with a single postsynaptic cell
    - the divergence ratio defined as the mean number of efferent
      connections created with a single presynaptic cell
    """

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
        self.dict_abv = dict_abv or {}
        """Dictionary of abbreviations for cell types"""

    def extract_strat_name(self, ps: ConnectivitySet):
        """
        Convert the name of a connection into its abbreviation.
        Connection strategy name is assumed to be in the form A_B, where
        A is the presynaptic type name and B is the postsynaptic type name.
        """
        for cs in self.scaffold.configuration.connectivity:
            if cs in ps.tag:
                if ps.tag == cs:
                    splits = ps.tag
                else:
                    splits = ps.tag.split(cs + "_", 1)[-1]
                break
        splits = splits.split("_")
        tag = []
        to_convert = []

        def parse_remains():
            if len(to_convert) > 0:
                to_test = "_".join(to_convert) + "_cell"
                if to_test in self.dict_abv:
                    return [self.dict_abv[to_test]]
                else:
                    return to_convert
            return []

        for text in splits:
            if "to" == text:
                tag.extend(parse_remains())
                to_convert = []
                tag.append("to")
            else:
                to_convert.append(text)
                to_test = "_".join(to_convert)
                if to_test in self.dict_abv:
                    tag.append(self.dict_abv[to_test])
                    to_convert = []
        tag.extend(parse_remains())

        return " ".join(tag)

    def update(self):
        super().update()
        self.reset_table()
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
            _, uniquePre_count = np.unique(combos[:, 0], axis=0, return_counts=True)
            _, uniquePost_count = np.unique(combos[:, 1], axis=0, return_counts=True)
            self.rows.append(self.extract_strat_name(ps))
            self._values.append([len(pre_locs), combo_counts, uniquePost_count, uniquePre_count])
            self.table_values.append(
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

    def get_synapse_counts(self):
        """
        Return a dictionary which gives for each connection name
        the number of synapses formed.
        The plot needs to be updated.

        :rtype: Dict[str, int]
        """
        return {ct: line[0] for ct, line in zip(self.rows, self._values)}

    def get_nb_synapse_per_pair(self):
        """
        Return a dictionary which gives for each connection name
        the number of synapses per unique pair of cell.
        The plot needs to be updated.

        :rtype: Dict[str, numpy.ndarray[int]]
        """
        return {ct: line[1] for ct, line in zip(self.rows, self._values)}

    def get_convergences(self):
        """
        Return a dictionary which gives for each connection name
        the divergence number of each unique postsynaptic cell.
        The plot needs to be updated.

        :rtype: Dict[str, numpy.ndarray[int]]
        """
        return {ct: line[2] for ct, line in zip(self.rows, self._values)}

    def get_divergences(self):
        """
        Return a dictionary which gives for each connection name
        the divergence number of each unique presynaptic cell.
        The plot needs to be updated.

        :rtype: Dict[str, numpy.ndarray[int]]
        """
        return {ct: line[3] for ct, line in zip(self.rows, self._values)}


class CellPlacement3D(ScaffoldPlot):
    """
    Plot the position of the cells in the scaffold in 3D space.
    """

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
        """List of cell type names to ignore in the plot."""

    def init_plot(self, **kwargs):
        if self.is_initialized:
            plt.close(self.figure)
        self.is_initialized = True
        self.is_plotted = False
        self.figure = plt.figure(figsize=self.fig_size, **kwargs)
        self.axes = self.figure.add_subplot(111, projection="3d")

    @staticmethod
    def set_axes_equal(ax):
        """
        Make axes of 3D plot have equal scale so that spheres appear as spheres,
        cubes as cubes, etc.

        :param ax: a matplotlib axis, e.g., as output from plt.gca().
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
            if not ct.name in self.ignored_ct:
                positions = ps.load_positions()
                if len(positions) > 0:
                    if ct.name not in self.dict_colors:
                        # default color is grey
                        color = [0.6, 0.6, 0.6]
                        alpha = 1.0
                    elif len(self.dict_colors[ct.name]) == 3:
                        color = self.dict_colors[ct.name]
                        alpha = 1.0
                    else:
                        *color, alpha = self.dict_colors[ct.name]
                    scale = np.power(ct.spatial.radius, 2)
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


class StructureReport(BSBReport):
    """
    Report of the scaffold neural network structure containing:

    - a table listing the number and density of each cell type placed.
    - a table listing the number, convergence and divergence of connections
      for each connected cell type pair
    - a 3D plot showing the soma location of each cell in the circuit.
    - a legend plot
    """

    def __init__(self, scaffold: Union[str, Scaffold], cell_type_info: List[PlotTypeInfo] = None):
        super().__init__(scaffold, cell_type_info)
        legend = Legend(
            (10, 2.5),
            3,
            dict_legend=dict(columnspacing=2.0, handletextpad=0.1, fontsize=20, loc="lower center"),
        )
        density_table = PlacementTable(
            (5, 0.25 * (len(self.scaffold.get_placement_sets()) + 1)),
            scaffold=self.scaffold,
            dict_abv=self.abbreviations,
        )
        connectivity_table = ConnectivityTable(
            (10, 0.25 * (len(self.scaffold.get_connectivity_sets()) + 1)),
            scaffold=self.scaffold,
            dict_abv=self.abbreviations,
        )
        plot3d = CellPlacement3D((10, 10), scaffold=self.scaffold)
        self.add_plot("density_table", density_table)
        self.add_plot("connectivity_table", connectivity_table)
        self.add_plot("placement_3d", plot3d)
        self.add_plot("legend", legend)
        legend.remove_ct(self.cell_names, ["glomerulus", "ubc_glomerulus"])

    def postprocessing(self):
        self.plots["legend"].set_axis_off()
        self.plots["density_table"].set_axis_off()
        self.plots["connectivity_table"].set_axis_off()


@config.node
class RunStructureReport(AfterConnectivityHook):
    """
    BSB postprocessing node to generate a scaffold structural report after running the connectivity jobs.
    """

    output_filename: str = config.attr(required=True)
    """Name of the pdf file to save the report."""

    def postprocess(self):
        report = StructureReport(self.scaffold, LIST_CT_INFO)
        report.print_report(self.output_filename)
