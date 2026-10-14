// #define DEBUG 1
/*
 *  eglif_cond_alpha_multisyn.cpp
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
 *  Generated from NESTML at time: 2025-05-22 16:02:48.475281
**/

// C++ includes:
#include <limits>

// Includes from libnestutil:
#include "numerics.h"

// Includes from nestkernel:
#include "exceptions.h"
#include "kernel_manager.h"
#include "nest_impl.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"
#include "lockptrdatum.h"

#include "eglif_cond_alpha_multisyn.h"
void
register_eglif_cond_alpha_multisyn( const std::string& name )
{
  nest::register_node_model< eglif_cond_alpha_multisyn >( name );
}

// ---------------------------------------------------------------------------
//   Recordables map
// ---------------------------------------------------------------------------
nest::RecordablesMap<eglif_cond_alpha_multisyn> eglif_cond_alpha_multisyn::recordablesMap_;
namespace nest
{

  // Override the create() method with one call to RecordablesMap::insert_()
  // for each quantity to be recorded.
template <> void RecordablesMap<eglif_cond_alpha_multisyn>::create()
  {
    // add state variables to recordables map
   insert_(eglif_cond_alpha_multisyn_names::_V_m, &eglif_cond_alpha_multisyn::get_V_m);
   insert_(eglif_cond_alpha_multisyn_names::_I_dep, &eglif_cond_alpha_multisyn::get_I_dep);
   insert_(eglif_cond_alpha_multisyn_names::_I_adap, &eglif_cond_alpha_multisyn::get_I_adap);
   insert_(eglif_cond_alpha_multisyn_names::_syn_kernel2__X__syn2_spike, &eglif_cond_alpha_multisyn::get_syn_kernel2__X__syn2_spike);
   insert_(eglif_cond_alpha_multisyn_names::_syn_kernel2__X__syn2_spike__d, &eglif_cond_alpha_multisyn::get_syn_kernel2__X__syn2_spike__d);
   insert_(eglif_cond_alpha_multisyn_names::_syn_kernel4__X__syn4_spike, &eglif_cond_alpha_multisyn::get_syn_kernel4__X__syn4_spike);
   insert_(eglif_cond_alpha_multisyn_names::_syn_kernel4__X__syn4_spike__d, &eglif_cond_alpha_multisyn::get_syn_kernel4__X__syn4_spike__d);
   insert_(eglif_cond_alpha_multisyn_names::_syn_kernel1__X__syn1_spike, &eglif_cond_alpha_multisyn::get_syn_kernel1__X__syn1_spike);
   insert_(eglif_cond_alpha_multisyn_names::_syn_kernel1__X__syn1_spike__d, &eglif_cond_alpha_multisyn::get_syn_kernel1__X__syn1_spike__d);
   insert_(eglif_cond_alpha_multisyn_names::_syn_kernel3__X__syn3_spike, &eglif_cond_alpha_multisyn::get_syn_kernel3__X__syn3_spike);
   insert_(eglif_cond_alpha_multisyn_names::_syn_kernel3__X__syn3_spike__d, &eglif_cond_alpha_multisyn::get_syn_kernel3__X__syn3_spike__d);
    // add recordable inline expressions to recordables map
	insert_(eglif_cond_alpha_multisyn_names::_I_syn, &eglif_cond_alpha_multisyn::get_I_syn);

    // Add vector variables  
  }
}
std::vector< std::tuple< int, int > > eglif_cond_alpha_multisyn::rport_to_nestml_buffer_idx =
{
  { eglif_cond_alpha_multisyn::SYN1_SPIKE, eglif_cond_alpha_multisyn::PORT_NOT_AVAILABLE },
  { eglif_cond_alpha_multisyn::SYN2_SPIKE, eglif_cond_alpha_multisyn::PORT_NOT_AVAILABLE },
  { eglif_cond_alpha_multisyn::SYN3_SPIKE, eglif_cond_alpha_multisyn::PORT_NOT_AVAILABLE },
  { eglif_cond_alpha_multisyn::SYN4_SPIKE, eglif_cond_alpha_multisyn::PORT_NOT_AVAILABLE },
};

// ---------------------------------------------------------------------------
//   Default constructors defining default parameters and state
//   Note: the implementation is empty. The initialization is of variables
//   is a part of eglif_cond_alpha_multisyn's constructor.
// ---------------------------------------------------------------------------

eglif_cond_alpha_multisyn::Parameters_::Parameters_()
{
}

eglif_cond_alpha_multisyn::State_::State_()
{
}

// ---------------------------------------------------------------------------
//   Parameter and state extractions and manipulation functions
// ---------------------------------------------------------------------------

eglif_cond_alpha_multisyn::Buffers_::Buffers_(eglif_cond_alpha_multisyn &n):
  logger_(n)
  , spike_inputs_( std::vector< nest::RingBuffer >( NUM_SPIKE_RECEPTORS ) )
  , spike_inputs_grid_sum_( std::vector< double >( NUM_SPIKE_RECEPTORS ) )
  , __s( nullptr ), __c( nullptr ), __e( nullptr )
{
  // Initialization of the remaining members is deferred to init_buffers_().
}

eglif_cond_alpha_multisyn::Buffers_::Buffers_(const Buffers_ &, eglif_cond_alpha_multisyn &n):
  logger_(n)
  , spike_inputs_( std::vector< nest::RingBuffer >( NUM_SPIKE_RECEPTORS ) )
  , spike_inputs_grid_sum_( std::vector< double >( NUM_SPIKE_RECEPTORS ) )
  , __s( nullptr ), __c( nullptr ), __e( nullptr )
{
  // Initialization of the remaining members is deferred to init_buffers_().
}

// ---------------------------------------------------------------------------
//   Default constructor for node
// ---------------------------------------------------------------------------

eglif_cond_alpha_multisyn::eglif_cond_alpha_multisyn():ArchivingNode(), P_(), S_(), B_(*this)
{
  init_state_internal_();
  recordablesMap_.create();
  pre_run_hook();
}

// ---------------------------------------------------------------------------
//   Copy constructor for node
// ---------------------------------------------------------------------------

eglif_cond_alpha_multisyn::eglif_cond_alpha_multisyn(const eglif_cond_alpha_multisyn& __n):
  ArchivingNode(), P_(__n.P_), S_(__n.S_), B_(__n.B_, *this) {

  // copy parameter struct P_
  P_.C_m = __n.P_.C_m;
  P_.tau_m = __n.P_.tau_m;
  P_.E_L = __n.P_.E_L;
  P_.t_ref = __n.P_.t_ref;
  P_.I_e = __n.P_.I_e;
  P_.V_min = __n.P_.V_min;
  P_.V_th = __n.P_.V_th;
  P_.lambda_0 = __n.P_.lambda_0;
  P_.tau_V = __n.P_.tau_V;
  P_.V_reset = __n.P_.V_reset;
  P_.k_1 = __n.P_.k_1;
  P_.k_2 = __n.P_.k_2;
  P_.k_adap = __n.P_.k_adap;
  P_.A1 = __n.P_.A1;
  P_.A2 = __n.P_.A2;
  P_.E_rev1 = __n.P_.E_rev1;
  P_.tau_syn1 = __n.P_.tau_syn1;
  P_.E_rev2 = __n.P_.E_rev2;
  P_.tau_syn2 = __n.P_.tau_syn2;
  P_.E_rev3 = __n.P_.E_rev3;
  P_.tau_syn3 = __n.P_.tau_syn3;
  P_.E_rev4 = __n.P_.E_rev4;
  P_.tau_syn4 = __n.P_.tau_syn4;

  // copy state struct S_
  S_.ode_state[State_::V_m] = __n.S_.ode_state[State_::V_m];
  S_.ode_state[State_::I_dep] = __n.S_.ode_state[State_::I_dep];
  S_.ode_state[State_::I_adap] = __n.S_.ode_state[State_::I_adap];
  S_.r = __n.S_.r;
  S_.ode_state[State_::syn_kernel2__X__syn2_spike] = __n.S_.ode_state[State_::syn_kernel2__X__syn2_spike];
  S_.ode_state[State_::syn_kernel2__X__syn2_spike__d] = __n.S_.ode_state[State_::syn_kernel2__X__syn2_spike__d];
  S_.ode_state[State_::syn_kernel4__X__syn4_spike] = __n.S_.ode_state[State_::syn_kernel4__X__syn4_spike];
  S_.ode_state[State_::syn_kernel4__X__syn4_spike__d] = __n.S_.ode_state[State_::syn_kernel4__X__syn4_spike__d];
  S_.ode_state[State_::syn_kernel1__X__syn1_spike] = __n.S_.ode_state[State_::syn_kernel1__X__syn1_spike];
  S_.ode_state[State_::syn_kernel1__X__syn1_spike__d] = __n.S_.ode_state[State_::syn_kernel1__X__syn1_spike__d];
  S_.ode_state[State_::syn_kernel3__X__syn3_spike] = __n.S_.ode_state[State_::syn_kernel3__X__syn3_spike];
  S_.ode_state[State_::syn_kernel3__X__syn3_spike__d] = __n.S_.ode_state[State_::syn_kernel3__X__syn3_spike__d];

  // copy internals V_
  V_.RefractoryCounts = __n.V_.RefractoryCounts;
  V_.__h = __n.V_.__h;
  V_.__P__I_dep__I_dep = __n.V_.__P__I_dep__I_dep;
  V_.__P__syn_kernel2__X__syn2_spike__syn_kernel2__X__syn2_spike = __n.V_.__P__syn_kernel2__X__syn2_spike__syn_kernel2__X__syn2_spike;
  V_.__P__syn_kernel2__X__syn2_spike__syn_kernel2__X__syn2_spike__d = __n.V_.__P__syn_kernel2__X__syn2_spike__syn_kernel2__X__syn2_spike__d;
  V_.__P__syn_kernel2__X__syn2_spike__d__syn_kernel2__X__syn2_spike = __n.V_.__P__syn_kernel2__X__syn2_spike__d__syn_kernel2__X__syn2_spike;
  V_.__P__syn_kernel2__X__syn2_spike__d__syn_kernel2__X__syn2_spike__d = __n.V_.__P__syn_kernel2__X__syn2_spike__d__syn_kernel2__X__syn2_spike__d;
  V_.__P__syn_kernel4__X__syn4_spike__syn_kernel4__X__syn4_spike = __n.V_.__P__syn_kernel4__X__syn4_spike__syn_kernel4__X__syn4_spike;
  V_.__P__syn_kernel4__X__syn4_spike__syn_kernel4__X__syn4_spike__d = __n.V_.__P__syn_kernel4__X__syn4_spike__syn_kernel4__X__syn4_spike__d;
  V_.__P__syn_kernel4__X__syn4_spike__d__syn_kernel4__X__syn4_spike = __n.V_.__P__syn_kernel4__X__syn4_spike__d__syn_kernel4__X__syn4_spike;
  V_.__P__syn_kernel4__X__syn4_spike__d__syn_kernel4__X__syn4_spike__d = __n.V_.__P__syn_kernel4__X__syn4_spike__d__syn_kernel4__X__syn4_spike__d;
  V_.__P__syn_kernel1__X__syn1_spike__syn_kernel1__X__syn1_spike = __n.V_.__P__syn_kernel1__X__syn1_spike__syn_kernel1__X__syn1_spike;
  V_.__P__syn_kernel1__X__syn1_spike__syn_kernel1__X__syn1_spike__d = __n.V_.__P__syn_kernel1__X__syn1_spike__syn_kernel1__X__syn1_spike__d;
  V_.__P__syn_kernel1__X__syn1_spike__d__syn_kernel1__X__syn1_spike = __n.V_.__P__syn_kernel1__X__syn1_spike__d__syn_kernel1__X__syn1_spike;
  V_.__P__syn_kernel1__X__syn1_spike__d__syn_kernel1__X__syn1_spike__d = __n.V_.__P__syn_kernel1__X__syn1_spike__d__syn_kernel1__X__syn1_spike__d;
  V_.__P__syn_kernel3__X__syn3_spike__syn_kernel3__X__syn3_spike = __n.V_.__P__syn_kernel3__X__syn3_spike__syn_kernel3__X__syn3_spike;
  V_.__P__syn_kernel3__X__syn3_spike__syn_kernel3__X__syn3_spike__d = __n.V_.__P__syn_kernel3__X__syn3_spike__syn_kernel3__X__syn3_spike__d;
  V_.__P__syn_kernel3__X__syn3_spike__d__syn_kernel3__X__syn3_spike = __n.V_.__P__syn_kernel3__X__syn3_spike__d__syn_kernel3__X__syn3_spike;
  V_.__P__syn_kernel3__X__syn3_spike__d__syn_kernel3__X__syn3_spike__d = __n.V_.__P__syn_kernel3__X__syn3_spike__d__syn_kernel3__X__syn3_spike__d;
}

// ---------------------------------------------------------------------------
//   Destructor for node
// ---------------------------------------------------------------------------

eglif_cond_alpha_multisyn::~eglif_cond_alpha_multisyn()
{
  // GSL structs may not have been allocated, so we need to protect destruction

  if (B_.__s)
  {
    gsl_odeiv_step_free( B_.__s );
  }

  if (B_.__c)
  {
    gsl_odeiv_control_free( B_.__c );
  }

  if (B_.__e)
  {
    gsl_odeiv_evolve_free( B_.__e );
  }
}

// ---------------------------------------------------------------------------
//   Node initialization functions
// ---------------------------------------------------------------------------
void eglif_cond_alpha_multisyn::calibrate_time( const nest::TimeConverter& tc )
{
  LOG( nest::M_WARNING,
    "eglif_cond_alpha_multisyn",
    "Simulation resolution has changed. Internal state and parameters of the model have been reset!" );

  init_state_internal_();
}
void eglif_cond_alpha_multisyn::init_state_internal_()
{
#ifdef DEBUG
  std::cout << "eglif_cond_alpha_multisyn::init_state_internal_()" << std::endl;
#endif

  const double __resolution = nest::Time::get_resolution().get_ms();  // do not remove, this is necessary for the resolution() function

  // use a default "good enough" value for the absolute error. It can be adjusted via `node.set()`
  P_.__gsl_error_tol = 1e-3;
  // initial values for parameters
    
    
    
    
    P_.C_m = 14.6; // as pF
    P_.tau_m = 9.125; // as ms
    P_.E_L = (-68.0); // as mV
    P_.t_ref = 1.59; // as ms
    P_.I_e = 3.711; // as pA
    P_.V_min = (-150.0); // as mV
    P_.V_th = (-53.0); // as mV
    P_.lambda_0 = 1.8 / 1; // as real
    P_.tau_V = 1.1; // as mV
    P_.V_reset = (-78.0); // as mV
    P_.k_1 = 1.887 / 1; // as real
    P_.k_2 = 1.096 / 1; // as real
    P_.k_adap = 2.025 / (1.0 * 1.0); // as real
    P_.A1 = 5.953; // as pA
    P_.A2 = 5.863; // as pA
    P_.E_rev1 = 0; // as mV
    P_.tau_syn1 = 0.2; // as ms
    P_.E_rev2 = (-80.0); // as mV
    P_.tau_syn2 = 2.0; // as ms
    P_.E_rev3 = 0.0; // as mV
    P_.tau_syn3 = 2.0; // as ms
    P_.E_rev4 = (-80.0); // as mV
    P_.tau_syn4 = 2.0; // as ms

  recompute_internal_variables();
  // initial values for state variables
    

    S_.ode_state[State_::V_m] = (-60.0); // as mV
    

    S_.ode_state[State_::I_dep] = 0; // as pA
    

    S_.ode_state[State_::I_adap] = 0; // as pA
    

    S_.r = 0; // as integer
    

    S_.ode_state[State_::syn_kernel2__X__syn2_spike] = 0; // as real
    

    S_.ode_state[State_::syn_kernel2__X__syn2_spike__d] = 0 * pow(1000.0, (-1)); // as 1 / s
    

    S_.ode_state[State_::syn_kernel4__X__syn4_spike] = 0; // as real
    

    S_.ode_state[State_::syn_kernel4__X__syn4_spike__d] = 0 * pow(1000.0, (-1)); // as 1 / s
    

    S_.ode_state[State_::syn_kernel1__X__syn1_spike] = 0; // as real
    

    S_.ode_state[State_::syn_kernel1__X__syn1_spike__d] = 0 * pow(1000.0, (-1)); // as 1 / s
    

    S_.ode_state[State_::syn_kernel3__X__syn3_spike] = 0; // as real
    

    S_.ode_state[State_::syn_kernel3__X__syn3_spike__d] = 0 * pow(1000.0, (-1)); // as 1 / s
}

void eglif_cond_alpha_multisyn::init_buffers_()
{
#ifdef DEBUG
  std::cout << "eglif_cond_alpha_multisyn::init_buffers_()" << std::endl;
#endif
  // spike input buffers
  get_spike_inputs_().clear();
  get_spike_inputs_grid_sum_().clear();

  // continuous time input buffers  

  get_I_stim().clear();
  B_.I_stim_grid_sum_ = 0;

  B_.logger_.reset();



  if ( not B_.__s )
  {
    B_.__s = gsl_odeiv_step_alloc( gsl_odeiv_step_rkf45, State_::STATE_VEC_SIZE );
  }
  else
  {
    gsl_odeiv_step_reset( B_.__s );
  }

  if ( not B_.__c )
  {
    B_.__c = gsl_odeiv_control_y_new( P_.__gsl_error_tol, 0.0 );
  }
  else
  {
    gsl_odeiv_control_init( B_.__c, P_.__gsl_error_tol, 0.0, 1.0, 0.0 );
  }

  if ( not B_.__e )
  {
    B_.__e = gsl_odeiv_evolve_alloc( State_::STATE_VEC_SIZE );
  }
  else
  {
    gsl_odeiv_evolve_reset( B_.__e );
  }

  B_.__sys.function = eglif_cond_alpha_multisyn_dynamics;
  B_.__sys.jacobian = nullptr;
  B_.__sys.dimension = State_::STATE_VEC_SIZE;
  B_.__sys.params = reinterpret_cast< void* >( this );
  B_.__step = nest::Time::get_resolution().get_ms();
  B_.__integration_step = nest::Time::get_resolution().get_ms();
}

void eglif_cond_alpha_multisyn::recompute_internal_variables(bool exclude_timestep) {
  const double __resolution = nest::Time::get_resolution().get_ms();  // do not remove, this is necessary for the resolution() function

  if (exclude_timestep) {    
      

      V_.RefractoryCounts = nest::Time(nest::Time::ms((double) (P_.t_ref))).get_steps(); // as integer
      

      V_.__P__I_dep__I_dep = 1.0 * std::exp((-V_.__h) * P_.k_1); // as real
      

      V_.__P__syn_kernel2__X__syn2_spike__syn_kernel2__X__syn2_spike = 1.0 * (V_.__h + P_.tau_syn2) * std::exp((-V_.__h) / P_.tau_syn2) / P_.tau_syn2; // as real
      

      V_.__P__syn_kernel2__X__syn2_spike__syn_kernel2__X__syn2_spike__d = 1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn2); // as real
      

      V_.__P__syn_kernel2__X__syn2_spike__d__syn_kernel2__X__syn2_spike = (-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn2) / pow(P_.tau_syn2, 2); // as real
      

      V_.__P__syn_kernel2__X__syn2_spike__d__syn_kernel2__X__syn2_spike__d = 1.0 * ((-V_.__h) + P_.tau_syn2) * std::exp((-V_.__h) / P_.tau_syn2) / P_.tau_syn2; // as real
      

      V_.__P__syn_kernel4__X__syn4_spike__syn_kernel4__X__syn4_spike = 1.0 * (V_.__h + P_.tau_syn4) * std::exp((-V_.__h) / P_.tau_syn4) / P_.tau_syn4; // as real
      

      V_.__P__syn_kernel4__X__syn4_spike__syn_kernel4__X__syn4_spike__d = 1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn4); // as real
      

      V_.__P__syn_kernel4__X__syn4_spike__d__syn_kernel4__X__syn4_spike = (-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn4) / pow(P_.tau_syn4, 2); // as real
      

      V_.__P__syn_kernel4__X__syn4_spike__d__syn_kernel4__X__syn4_spike__d = 1.0 * ((-V_.__h) + P_.tau_syn4) * std::exp((-V_.__h) / P_.tau_syn4) / P_.tau_syn4; // as real
      

      V_.__P__syn_kernel1__X__syn1_spike__syn_kernel1__X__syn1_spike = 1.0 * (V_.__h + P_.tau_syn1) * std::exp((-V_.__h) / P_.tau_syn1) / P_.tau_syn1; // as real
      

      V_.__P__syn_kernel1__X__syn1_spike__syn_kernel1__X__syn1_spike__d = 1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn1); // as real
      

      V_.__P__syn_kernel1__X__syn1_spike__d__syn_kernel1__X__syn1_spike = (-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn1) / pow(P_.tau_syn1, 2); // as real
      

      V_.__P__syn_kernel1__X__syn1_spike__d__syn_kernel1__X__syn1_spike__d = 1.0 * ((-V_.__h) + P_.tau_syn1) * std::exp((-V_.__h) / P_.tau_syn1) / P_.tau_syn1; // as real
      

      V_.__P__syn_kernel3__X__syn3_spike__syn_kernel3__X__syn3_spike = 1.0 * (V_.__h + P_.tau_syn3) * std::exp((-V_.__h) / P_.tau_syn3) / P_.tau_syn3; // as real
      

      V_.__P__syn_kernel3__X__syn3_spike__syn_kernel3__X__syn3_spike__d = 1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn3); // as real
      

      V_.__P__syn_kernel3__X__syn3_spike__d__syn_kernel3__X__syn3_spike = (-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn3) / pow(P_.tau_syn3, 2); // as real
      

      V_.__P__syn_kernel3__X__syn3_spike__d__syn_kernel3__X__syn3_spike__d = 1.0 * ((-V_.__h) + P_.tau_syn3) * std::exp((-V_.__h) / P_.tau_syn3) / P_.tau_syn3; // as real
  }
  else {    
      

      V_.RefractoryCounts = nest::Time(nest::Time::ms((double) (P_.t_ref))).get_steps(); // as integer
      

      V_.__h = __resolution; // as ms
      

      V_.__P__I_dep__I_dep = 1.0 * std::exp((-V_.__h) * P_.k_1); // as real
      

      V_.__P__syn_kernel2__X__syn2_spike__syn_kernel2__X__syn2_spike = 1.0 * (V_.__h + P_.tau_syn2) * std::exp((-V_.__h) / P_.tau_syn2) / P_.tau_syn2; // as real
      

      V_.__P__syn_kernel2__X__syn2_spike__syn_kernel2__X__syn2_spike__d = 1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn2); // as real
      

      V_.__P__syn_kernel2__X__syn2_spike__d__syn_kernel2__X__syn2_spike = (-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn2) / pow(P_.tau_syn2, 2); // as real
      

      V_.__P__syn_kernel2__X__syn2_spike__d__syn_kernel2__X__syn2_spike__d = 1.0 * ((-V_.__h) + P_.tau_syn2) * std::exp((-V_.__h) / P_.tau_syn2) / P_.tau_syn2; // as real
      

      V_.__P__syn_kernel4__X__syn4_spike__syn_kernel4__X__syn4_spike = 1.0 * (V_.__h + P_.tau_syn4) * std::exp((-V_.__h) / P_.tau_syn4) / P_.tau_syn4; // as real
      

      V_.__P__syn_kernel4__X__syn4_spike__syn_kernel4__X__syn4_spike__d = 1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn4); // as real
      

      V_.__P__syn_kernel4__X__syn4_spike__d__syn_kernel4__X__syn4_spike = (-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn4) / pow(P_.tau_syn4, 2); // as real
      

      V_.__P__syn_kernel4__X__syn4_spike__d__syn_kernel4__X__syn4_spike__d = 1.0 * ((-V_.__h) + P_.tau_syn4) * std::exp((-V_.__h) / P_.tau_syn4) / P_.tau_syn4; // as real
      

      V_.__P__syn_kernel1__X__syn1_spike__syn_kernel1__X__syn1_spike = 1.0 * (V_.__h + P_.tau_syn1) * std::exp((-V_.__h) / P_.tau_syn1) / P_.tau_syn1; // as real
      

      V_.__P__syn_kernel1__X__syn1_spike__syn_kernel1__X__syn1_spike__d = 1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn1); // as real
      

      V_.__P__syn_kernel1__X__syn1_spike__d__syn_kernel1__X__syn1_spike = (-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn1) / pow(P_.tau_syn1, 2); // as real
      

      V_.__P__syn_kernel1__X__syn1_spike__d__syn_kernel1__X__syn1_spike__d = 1.0 * ((-V_.__h) + P_.tau_syn1) * std::exp((-V_.__h) / P_.tau_syn1) / P_.tau_syn1; // as real
      

      V_.__P__syn_kernel3__X__syn3_spike__syn_kernel3__X__syn3_spike = 1.0 * (V_.__h + P_.tau_syn3) * std::exp((-V_.__h) / P_.tau_syn3) / P_.tau_syn3; // as real
      

      V_.__P__syn_kernel3__X__syn3_spike__syn_kernel3__X__syn3_spike__d = 1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn3); // as real
      

      V_.__P__syn_kernel3__X__syn3_spike__d__syn_kernel3__X__syn3_spike = (-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn3) / pow(P_.tau_syn3, 2); // as real
      

      V_.__P__syn_kernel3__X__syn3_spike__d__syn_kernel3__X__syn3_spike__d = 1.0 * ((-V_.__h) + P_.tau_syn3) * std::exp((-V_.__h) / P_.tau_syn3) / P_.tau_syn3; // as real
  }
}
void eglif_cond_alpha_multisyn::pre_run_hook() {
  B_.logger_.init();

  // parameters might have changed -- recompute internals
  recompute_internal_variables();

  // buffers B_
  B_.spike_inputs_.resize(NUM_SPIKE_RECEPTORS);
  B_.spike_inputs_grid_sum_.resize(NUM_SPIKE_RECEPTORS);
}

// ---------------------------------------------------------------------------
//   Update and spike handling functions
// ---------------------------------------------------------------------------

extern "C" inline int eglif_cond_alpha_multisyn_dynamics(double __time, const double ode_state[], double f[], void* pnode)
{
  typedef eglif_cond_alpha_multisyn::State_ State_;
  // get access to node so we can almost work as in a member function
  assert( pnode );
  const eglif_cond_alpha_multisyn& node = *( reinterpret_cast< eglif_cond_alpha_multisyn* >( pnode ) );

  // ode_state[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.ode_state[].
  f[State_::V_m] = ((-node.P_.E_L) + ode_state[State_::V_m]) / node.P_.tau_m + 1.0 * node.P_.E_rev1 * ode_state[State_::syn_kernel1__X__syn1_spike] / node.P_.C_m + 1.0 * node.P_.E_rev2 * ode_state[State_::syn_kernel2__X__syn2_spike] / node.P_.C_m + node.P_.I_e / node.P_.C_m + (1.0 * node.P_.E_rev3 * ode_state[State_::syn_kernel3__X__syn3_spike] + 1.0 * node.P_.E_rev4 * ode_state[State_::syn_kernel4__X__syn4_spike] - ode_state[State_::I_adap] + ode_state[State_::I_dep] + node.B_.I_stim_grid_sum_ - 1.0 * ode_state[State_::V_m] * ode_state[State_::syn_kernel1__X__syn1_spike] - 1.0 * ode_state[State_::V_m] * ode_state[State_::syn_kernel2__X__syn2_spike] - 1.0 * ode_state[State_::V_m] * ode_state[State_::syn_kernel3__X__syn3_spike] - 1.0 * ode_state[State_::V_m] * ode_state[State_::syn_kernel4__X__syn4_spike]) / node.P_.C_m;
  f[State_::I_dep] = (-1.0) * ode_state[State_::I_dep] * node.P_.k_1;
  f[State_::I_adap] = (-1.0) * node.P_.E_L * node.P_.k_adap - 1.0 * ode_state[State_::I_adap] * node.P_.k_2 + 1.0 * ode_state[State_::V_m] * node.P_.k_adap;
  f[State_::syn_kernel2__X__syn2_spike] = 1.0 * ode_state[State_::syn_kernel2__X__syn2_spike__d];
  f[State_::syn_kernel2__X__syn2_spike__d] = (-ode_state[State_::syn_kernel2__X__syn2_spike]) / pow(node.P_.tau_syn2, 2) - 2 * ode_state[State_::syn_kernel2__X__syn2_spike__d] / node.P_.tau_syn2;
  f[State_::syn_kernel4__X__syn4_spike] = 1.0 * ode_state[State_::syn_kernel4__X__syn4_spike__d];
  f[State_::syn_kernel4__X__syn4_spike__d] = (-ode_state[State_::syn_kernel4__X__syn4_spike]) / pow(node.P_.tau_syn4, 2) - 2 * ode_state[State_::syn_kernel4__X__syn4_spike__d] / node.P_.tau_syn4;
  f[State_::syn_kernel1__X__syn1_spike] = 1.0 * ode_state[State_::syn_kernel1__X__syn1_spike__d];
  f[State_::syn_kernel1__X__syn1_spike__d] = (-ode_state[State_::syn_kernel1__X__syn1_spike]) / pow(node.P_.tau_syn1, 2) - 2 * ode_state[State_::syn_kernel1__X__syn1_spike__d] / node.P_.tau_syn1;
  f[State_::syn_kernel3__X__syn3_spike] = 1.0 * ode_state[State_::syn_kernel3__X__syn3_spike__d];
  f[State_::syn_kernel3__X__syn3_spike__d] = (-ode_state[State_::syn_kernel3__X__syn3_spike]) / pow(node.P_.tau_syn3, 2) - 2 * ode_state[State_::syn_kernel3__X__syn3_spike__d] / node.P_.tau_syn3;
  return GSL_SUCCESS;
}

void eglif_cond_alpha_multisyn::update(nest::Time const & origin,const long from, const long to)
{
  const double __resolution = nest::Time::get_resolution().get_ms();  // do not remove, this is necessary for the resolution() function



  for ( long lag = from ; lag < to ; ++lag )
  {
    auto get_t = [origin, lag](){ return nest::Time( nest::Time::step( origin.get_steps() + lag + 1) ).get_ms(); };

    for (long i = 0; i < NUM_SPIKE_RECEPTORS; ++i)
    {
        get_spike_inputs_grid_sum_()[i] = get_spike_inputs_()[i].get_value(lag);
    }
    B_.I_stim_grid_sum_ = get_I_stim().get_value(lag);

    // NESTML generated code for the update block
  S_.ode_state[State_::V_m] = std::max(S_.ode_state[State_::V_m], P_.V_min);
  double I_dep__tmp = S_.ode_state[State_::I_dep] * V_.__P__I_dep__I_dep;
  double syn_kernel2__X__syn2_spike__tmp = V_.__P__syn_kernel2__X__syn2_spike__syn_kernel2__X__syn2_spike * S_.ode_state[State_::syn_kernel2__X__syn2_spike] + V_.__P__syn_kernel2__X__syn2_spike__syn_kernel2__X__syn2_spike__d * S_.ode_state[State_::syn_kernel2__X__syn2_spike__d];
  double syn_kernel2__X__syn2_spike__d__tmp = V_.__P__syn_kernel2__X__syn2_spike__d__syn_kernel2__X__syn2_spike * S_.ode_state[State_::syn_kernel2__X__syn2_spike] + V_.__P__syn_kernel2__X__syn2_spike__d__syn_kernel2__X__syn2_spike__d * S_.ode_state[State_::syn_kernel2__X__syn2_spike__d];
  double syn_kernel4__X__syn4_spike__tmp = V_.__P__syn_kernel4__X__syn4_spike__syn_kernel4__X__syn4_spike * S_.ode_state[State_::syn_kernel4__X__syn4_spike] + V_.__P__syn_kernel4__X__syn4_spike__syn_kernel4__X__syn4_spike__d * S_.ode_state[State_::syn_kernel4__X__syn4_spike__d];
  double syn_kernel4__X__syn4_spike__d__tmp = V_.__P__syn_kernel4__X__syn4_spike__d__syn_kernel4__X__syn4_spike * S_.ode_state[State_::syn_kernel4__X__syn4_spike] + V_.__P__syn_kernel4__X__syn4_spike__d__syn_kernel4__X__syn4_spike__d * S_.ode_state[State_::syn_kernel4__X__syn4_spike__d];
  double syn_kernel1__X__syn1_spike__tmp = V_.__P__syn_kernel1__X__syn1_spike__syn_kernel1__X__syn1_spike * S_.ode_state[State_::syn_kernel1__X__syn1_spike] + V_.__P__syn_kernel1__X__syn1_spike__syn_kernel1__X__syn1_spike__d * S_.ode_state[State_::syn_kernel1__X__syn1_spike__d];
  double syn_kernel1__X__syn1_spike__d__tmp = V_.__P__syn_kernel1__X__syn1_spike__d__syn_kernel1__X__syn1_spike * S_.ode_state[State_::syn_kernel1__X__syn1_spike] + V_.__P__syn_kernel1__X__syn1_spike__d__syn_kernel1__X__syn1_spike__d * S_.ode_state[State_::syn_kernel1__X__syn1_spike__d];
  double syn_kernel3__X__syn3_spike__tmp = V_.__P__syn_kernel3__X__syn3_spike__syn_kernel3__X__syn3_spike * S_.ode_state[State_::syn_kernel3__X__syn3_spike] + V_.__P__syn_kernel3__X__syn3_spike__syn_kernel3__X__syn3_spike__d * S_.ode_state[State_::syn_kernel3__X__syn3_spike__d];
  double syn_kernel3__X__syn3_spike__d__tmp = V_.__P__syn_kernel3__X__syn3_spike__d__syn_kernel3__X__syn3_spike * S_.ode_state[State_::syn_kernel3__X__syn3_spike] + V_.__P__syn_kernel3__X__syn3_spike__d__syn_kernel3__X__syn3_spike__d * S_.ode_state[State_::syn_kernel3__X__syn3_spike__d];
  double __t = 0;
  // numerical integration with adaptive step size control:
  // ------------------------------------------------------
  // gsl_odeiv_evolve_apply performs only a single numerical
  // integration step, starting from t and bounded by step;
  // the while-loop ensures integration over the whole simulation
  // step (0, step] if more than one integration step is needed due
  // to a small integration step size;
  // note that (t+IntegrationStep > step) leads to integration over
  // (t, step] and afterwards setting t to step, but it does not
  // enforce setting IntegrationStep to step-t; this is of advantage
  // for a consistent and efficient integration across subsequent
  // simulation intervals
  while ( __t < B_.__step )
  {

    const int status = gsl_odeiv_evolve_apply(B_.__e,
                                              B_.__c,
                                              B_.__s,
                                              &B_.__sys,              // system of ODE
                                              &__t,                   // from t
                                              B_.__step,              // to t <= step
                                              &B_.__integration_step, // integration step size
                                              S_.ode_state);          // neuronal state

    if ( status != GSL_SUCCESS )
    {
      throw nest::GSLSolverFailure( get_name(), status );
    }
  }
  /* replace analytically solvable variables with precisely integrated values  */
  S_.ode_state[State_::I_dep] = I_dep__tmp;
  S_.ode_state[State_::syn_kernel2__X__syn2_spike] = syn_kernel2__X__syn2_spike__tmp;
  S_.ode_state[State_::syn_kernel2__X__syn2_spike__d] = syn_kernel2__X__syn2_spike__d__tmp;
  S_.ode_state[State_::syn_kernel4__X__syn4_spike] = syn_kernel4__X__syn4_spike__tmp;
  S_.ode_state[State_::syn_kernel4__X__syn4_spike__d] = syn_kernel4__X__syn4_spike__d__tmp;
  S_.ode_state[State_::syn_kernel1__X__syn1_spike] = syn_kernel1__X__syn1_spike__tmp;
  S_.ode_state[State_::syn_kernel1__X__syn1_spike__d] = syn_kernel1__X__syn1_spike__d__tmp;
  S_.ode_state[State_::syn_kernel3__X__syn3_spike] = syn_kernel3__X__syn3_spike__tmp;
  S_.ode_state[State_::syn_kernel3__X__syn3_spike__d] = syn_kernel3__X__syn3_spike__d__tmp;
  S_.ode_state[State_::syn_kernel2__X__syn2_spike__d] += ((0.001 * B_.spike_inputs_grid_sum_[SYN2_SPIKE - MIN_SPIKE_RECEPTOR])) * (numerics::e / P_.tau_syn2) / (1 / 1000.0);
  S_.ode_state[State_::syn_kernel4__X__syn4_spike__d] += ((0.001 * B_.spike_inputs_grid_sum_[SYN4_SPIKE - MIN_SPIKE_RECEPTOR])) * (numerics::e / P_.tau_syn4) / (1 / 1000.0);
  S_.ode_state[State_::syn_kernel1__X__syn1_spike__d] += ((0.001 * B_.spike_inputs_grid_sum_[SYN1_SPIKE - MIN_SPIKE_RECEPTOR])) * (numerics::e / P_.tau_syn1) / (1 / 1000.0);
  S_.ode_state[State_::syn_kernel3__X__syn3_spike__d] += ((0.001 * B_.spike_inputs_grid_sum_[SYN3_SPIKE - MIN_SPIKE_RECEPTOR])) * (numerics::e / P_.tau_syn3) / (1 / 1000.0);
  if (S_.r > 0)
  {  
    S_.r -= 1;
    S_.ode_state[State_::V_m] = P_.V_reset;
  }
  else
  {  
    double lambda = P_.lambda_0 * std::exp((S_.ode_state[State_::V_m] - P_.V_th) / P_.tau_V);
    if (((0) + (1) * nest::get_vp_specific_rng( get_thread() )->drand()) < (-numerics::expm1((-lambda) * __resolution)))
    {  
      S_.r = V_.RefractoryCounts;
      S_.ode_state[State_::I_dep] = P_.A1;
      S_.ode_state[State_::I_adap] += P_.A2;
      S_.ode_state[State_::V_m] = P_.V_reset;
      set_spiketime(nest::Time::step(origin.get_steps()+lag+1));
      nest::SpikeEvent se;
      nest::kernel().event_delivery_manager.send(*this, se, lag);
    }
  }
    // voltage logging
    B_.logger_.record_data(origin.get_steps() + lag);
  }
}

// Do not move this function as inline to h-file. It depends on
// universal_data_logger_impl.h being included here.
void eglif_cond_alpha_multisyn::handle(nest::DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}


void eglif_cond_alpha_multisyn::handle(nest::SpikeEvent &e)
{
  assert(e.get_delay_steps() > 0);
  assert( e.get_rport() < B_.spike_inputs_.size() );

  double weight = e.get_weight();
  size_t nestml_buffer_idx = 0;
  if ( weight >= 0.0 )
  {
    nestml_buffer_idx = std::get<0>(rport_to_nestml_buffer_idx[e.get_rport()]);
  }
  else
  {
    nestml_buffer_idx = std::get<1>(rport_to_nestml_buffer_idx[e.get_rport()]);
    if ( nestml_buffer_idx == eglif_cond_alpha_multisyn::PORT_NOT_AVAILABLE )
    {
      nestml_buffer_idx = std::get<0>(rport_to_nestml_buffer_idx[e.get_rport()]);
    }
    weight = -weight;
  }
  B_.spike_inputs_[ nestml_buffer_idx - MIN_SPIKE_RECEPTOR ].add_value(
    e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin() ),
    weight * e.get_multiplicity() );
}

void eglif_cond_alpha_multisyn::handle(nest::CurrentEvent& e)
{
  assert(e.get_delay_steps() > 0);

  const double current = e.get_current();     // we assume that in NEST, this returns a current in pA
  const double weight = e.get_weight();
  get_I_stim().add_value(
               e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
               weight * current );
}

