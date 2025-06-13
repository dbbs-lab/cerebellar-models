"""
Module for the Report class.
"""

from typing import List, Union
from warnings import warn

import matplotlib.backends.backend_pdf
import numpy as np
from bsb import Scaffold, from_storage

from cerebellar_models.analysis.plots import Plot, ScaffoldPlot


class PlotTypeInfo:
    """
    Class storing information about element plotted, usually cell types or fibers.
    """

    def __init__(self, name: str, color, abbreviation: str = None):
        self.name = name
        """Name of the element."""
        self.abbreviation = abbreviation or name
        """Abbreviation of the element to print in the plot."""
        self.color = color
        """Color of the element."""


LIST_CT_INFO = [
    PlotTypeInfo("mossy_fibers", [0.847, 0, 0.451, 1.0], "mf"),
    PlotTypeInfo("glomerulus", [0.847, 0, 0.451, 1.0], "glom"),
    PlotTypeInfo("granule_cell", [0.7, 0.15, 0.15, 0.5], "GrC"),
    PlotTypeInfo("ascending_axon", [0.7, 0.15, 0.15, 0.5], "aa"),
    PlotTypeInfo("parallel_fiber", [0.7, 0.15, 0.15, 0.5], "pf"),
    PlotTypeInfo("unipolar_brush_cell", [0.196, 0.808, 0.988, 1.0], "UBC"),
    PlotTypeInfo("ubc_glomerulus", [0.196, 0.808, 0.988, 1.0], "ubc_glom"),
    PlotTypeInfo("golgi_cell", [0, 0.45, 0.7, 1.0], "GoC"),
    PlotTypeInfo("purkinje_cell", [0.275, 0.800, 0.275, 1.0], "PC"),
    PlotTypeInfo("basket_cell", [1, 0.647, 0, 1.0], "BC"),
    PlotTypeInfo("stellate_cell", [1, 0.84, 0, 1.0], "SC"),
    PlotTypeInfo("dcn_p", [0.3, 0.3, 0.3, 1.0], "DCN_P"),
    PlotTypeInfo("dcn_i", [0.635, 0, 0.145, 1.0], "DCN_I"),
    PlotTypeInfo("io", [0.46, 0.376, 0.54, 1.0], "IO"),
]
"""Cell and fiber information for plotting cerebellum BSB circuit."""


class Report:
    """
    Class interfacing a list of matplotlib plots to save in an external file.
    """

    def __init__(self, cell_types_info: List[PlotTypeInfo] = None):
        self.plots = {}
        """Dictionary mapping the report's plots' name to their instance"""
        self.cell_types_info = cell_types_info or []
        """List of PlotTypeInfo for each element to plot."""

    @property
    def colors(self):
        """
        Dictionary from the name of the elements to plot to its color.
        """
        return {ct.name: ct.color for ct in self.cell_types_info}

    @property
    def abbreviations(self):
        """
        Dictionary from the name of the elements to plot to its abbreviation.
        """
        return {ct.name: ct.abbreviation for ct in self.cell_types_info}

    def set_plot_colors(self, plot: Plot):
        """
        Set the plot color dictionary according to the report's
        """
        for k, v in self.colors.items():
            plot.set_color(k, v)

    def set_color(self, key: str, color: np.ndarray[float]):
        """
        Set a color for an element to plot.
        Colors must be an array of type RGB or RGBA.
        """
        if len(color) != 3 and len(color) != 4:
            raise ValueError("Color must be an array of size 3 or 4.")
        info_names = [i.name for i in self.cell_types_info]
        new_color = np.array(color, dtype=float)
        if key in info_names:
            self.cell_types_info[info_names.index(key)].color = new_color
        else:
            self.cell_types_info.append(PlotTypeInfo(key, new_color, key))
        for plot in self.plots.values():
            self.set_plot_colors(plot)

    def add_plot(self, name: str, plot: Plot):
        """
        Add a plot to the list of the report's plots.
        """
        if name in list(self.plots.keys()):
            warn("A plot named '{}' already exists in the report. Skipping it".format(name))
            return
        self.set_plot_colors(plot)
        self.plots[name] = plot

    def preprocessing(self):
        """
        Function to apply modifications to the report's plots
        after the plotting is done.
        """
        pass

    def print_report(self, output_name: str, dpi: int = 200, pad: int = 0, **kwargs):
        """
        Print the report and export it in a pdf file.
        """
        pdf = matplotlib.backends.backend_pdf.PdfPages(output_name)
        self.preprocessing()
        for name, plot in self.plots.items():
            if not plot.is_plotted:
                plot.plot()
            plot.figure.tight_layout(pad=pad)
            plot.figure.savefig(pdf, format="pdf", dpi=dpi, **kwargs)
        pdf.close()

    def save_plot(self, plot_name: str, output_name: str, dpi: int = 200, pad: int = 0, **kwargs):
        """
        Save one of the report's plots as a separate figure.
        """
        if plot_name in self.plots:
            self.plots[plot_name].save_figure(output_name, dpi, pad, **kwargs)

    def show_plot(self, plot_name):  # pragma:nocover
        """
        Show one of the report's plots.
        """
        if plot_name in self.plots:
            self.plots[plot_name].show()

    def show(self):  # pragma:nocover
        """
        Show all report's plots one by one.
        """
        for plot in self.plots.values():
            plot.show()


class BSBReport(Report):
    """
    Class interfacing a list of matplotlib plots for analysis of BSB reconstructions.
    """

    def __init__(self, scaffold: Union[str, Scaffold], cell_types_info: List[PlotTypeInfo] = None):
        super().__init__(cell_types_info or LIST_CT_INFO)
        if isinstance(scaffold, Scaffold):
            self.scaffold = scaffold
        else:
            self.scaffold = from_storage(scaffold)

    @property
    def cell_names(self):
        """
        List of the name of the elements to plot.
        """
        return list(self.scaffold.cell_types.keys())

    @property
    def labelled_cell_names(self):
        ct_names = []
        for ps in self.scaffold.get_placement_sets():
            ct_name = ps.cell_type.name
            for labels in ps.get_unique_labels():
                ct_names.append(ScaffoldPlot.get_labelled_ct_name(ct_name, labels))
        return ct_names

    @property
    def labelled_abbreviations(self):
        dict_abv = self.abbreviations.copy()
        for ct_name in self.labelled_cell_names:
            if ct_name not in self.cell_names:
                ct, label = ct_name.rsplit("_", 1)
                dict_abv[ct_name] = self.abbreviations[ct] + "_" + label
        return dict_abv

    def add_plot(self, name: str, plot: Plot):
        super().add_plot(name, plot)
        if isinstance(plot, ScaffoldPlot):
            plot.set_scaffold(self.scaffold)
