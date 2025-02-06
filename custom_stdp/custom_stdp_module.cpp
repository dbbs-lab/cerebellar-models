
/*
*  nestml_allmodels_module.cpp
*
*  This file is part of NEST.
*
*  Copyright (C) 2004 The NEST Initiative
*
*  NEST is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 2 of the License, or
*  (at your option) any later version.
*
*  NEST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
*
*  2024-11-05 14:48:37.050316
*/

// Include from NEST
#include "nest_extension_interface.h"

// include headers with your own stuff
#include "iaf_cond_alpha_neuron_nestml.h"
#include "eglif_io_nestml.h"
#include "eglif_pc_nestml.h"
#include "eglif_gr_nestml.h"
#include "eglif_mli.h"
#include "stdp_synapse_sinexp.h"
#include "stdp_synapse_cosexp.h"
#include "stdp_synapse_alpha.h"



class custom_stdp_module : public nest::NESTExtensionInterface
{
  public:
    custom_stdp_module() {}
    ~custom_stdp_module() {}

    void initialize() override;
};

// Define module instance outside of namespace to avoid name-mangling problems
custom_stdp_module custom_stdp_module_LTX_module;

void custom_stdp_module::initialize()
{
    // register neurons, synapse or device models
    register_iaf_cond_alpha_neuron_nestml("iaf_cond_alpha_neuron_nestml");
    nest::register_eglif_io_nestml("eglif_io_nestml");
    nest::register_eglif_pc_nestml("eglif_pc_nestml");
    nest::register_eglif_gr_nestml("eglif_gr_nestml");
    nest::register_eglif_mli("eglif_mli");
    nest::register_stdp_synapse_sinexp( "stdp_synapse_sinexp" );
    nest::register_stdp_synapse_sinexp( "stdp_synapse_cosexp" );
    nest::register_stdp_synapse_alpha( "stdp_synapse_alpha" );

}
