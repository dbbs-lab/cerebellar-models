#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from bsb import config
from bsb_neuron.cell import ArborizedModel

# from library_function import adjust_morpho_diam


@config.node
class AddIext(ArborizedModel):
    amp = config.attr(type=float, default=0.0)

    def create(self, id, pos, morpho, rot, additional):
        cell = super().create(id, pos, morpho, rot, additional)
        if self.amp > 0:
            cell.soma[0].synapse_types["Iext"].parameters["amp"] = self.amp
            locs = cell.soma[0].locations
            cell.insert_synapse("Iext", locs[0])
        return cell


@config.node
class MorphoAdjModel(AddIext):

    def create(self, id, pos, morpho, rot, additional):
        cell = super().create(id, pos, morpho, rot, additional)
        return self.adjust_morpho_diam(cell)

    def adjust_morpho_diam(self, cell):

        for i_ses, ses_i in enumerate(cell.sections):
            if len(ses_i.children()) > 1:
                for i_chil, chil_i in enumerate(ses_i.children()):
                    if not (chil_i in cell.soma):
                        chil_i.pt3dchange(0, chil_i.diam3d(1))
        return cell


@config.node
class MorphoAdjGrCModel(MorphoAdjModel):

    def create(self, id, pos, morpho, rot, additional):
        cell = super().create(id, pos, morpho, rot, additional)
        params = {"ascending_axon": 18, "parallel_fiber": 142}
        return self.set_nseg(cell, params)

    def set_nseg(self, cell, params):
        for label_i, val_i in params.items():
            secs = cell.get_sections_with_label(label_i)
            for i_secs, secs_i in enumerate(secs):
                secs_i.set_segments(val_i)
        return cell


@config.node
class MakeAutisticModel(MorphoAdjGrCModel):

    def create(self, id, pos, morpho, rot, additional):
        cell = super().create(id, pos, morpho, rot, additional)
        cell.soma[1].diam = 5.3
        return cell
