// #define DEBUG 1
/*
 *  eglif_dcnp.cpp
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
 *  Generated from NESTML at time: 2024-10-25 09:35:56.352139
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

#include "eglif_dcnp.h"

using namespace nest;
// ---------------------------------------------------------------------------
//   Recordables map
// ---------------------------------------------------------------------------
nest::RecordablesMap<eglif_dcnp> eglif_dcnp::recordablesMap_;

void
nest::register_eglif_dcnp( const std::string& name )
{
  register_node_model< eglif_dcnp >( name );
}


namespace nest
{

  // Override the create() method with one call to RecordablesMap::insert_()
  // for each quantity to be recorded.
template <> 
void 
RecordablesMap<eglif_dcnp>::create()
  {
    // add state variables to recordables map
   insert_(eglif_dcnp_names::_V_m, &eglif_dcnp::get_V_m);
   insert_(eglif_dcnp_names::_I_adap, &eglif_dcnp::get_I_adap);
   insert_(eglif_dcnp_names::_I_dep, &eglif_dcnp::get_I_dep);
   insert_(eglif_dcnp_names::_lambda, &eglif_dcnp::get_lambda);
   insert_(eglif_dcnp_names::_tick, &eglif_dcnp::get_tick);
   insert_(eglif_dcnp_names::_cf_buffer, &eglif_dcnp::get_cf_buffer);
   insert_(eglif_dcnp_names::_gr_buffer, &eglif_dcnp::get_gr_buffer);
   insert_(eglif_dcnp_names::_last_io, &eglif_dcnp::get_last_io);
   insert_(eglif_dcnp_names::_complex_flag, &eglif_dcnp::get_complex_flag);
   insert_(eglif_dcnp_names::_g1__X__rec1, &eglif_dcnp::get_g1__X__rec1);
   insert_(eglif_dcnp_names::_g1__X__rec1__d, &eglif_dcnp::get_g1__X__rec1__d);
   insert_(eglif_dcnp_names::_g4__X__rec4, &eglif_dcnp::get_g4__X__rec4);
   insert_(eglif_dcnp_names::_g4__X__rec4__d, &eglif_dcnp::get_g4__X__rec4__d);
   insert_(eglif_dcnp_names::_g2__X__rec2, &eglif_dcnp::get_g2__X__rec2);
   insert_(eglif_dcnp_names::_g2__X__rec2__d, &eglif_dcnp::get_g2__X__rec2__d);
   insert_(eglif_dcnp_names::_g3__X__rec3, &eglif_dcnp::get_g3__X__rec3);
   insert_(eglif_dcnp_names::_g3__X__rec3__d, &eglif_dcnp::get_g3__X__rec3__d);
    // add recordable inline expressions to recordables map
	insert_(eglif_dcnp_names::_I_syn, &eglif_dcnp::get_I_syn);
	insert_(eglif_dcnp_names::_I_tot, &eglif_dcnp::get_I_tot);

    // Add vector variables  
  }
}
std::vector< std::tuple< int, int > > eglif_dcnp::rport_to_nestml_buffer_idx =
{
  { eglif_dcnp::REC1, eglif_dcnp::PORT_NOT_AVAILABLE },
  { eglif_dcnp::REC2, eglif_dcnp::PORT_NOT_AVAILABLE },
  { eglif_dcnp::REC3, eglif_dcnp::PORT_NOT_AVAILABLE },
  { eglif_dcnp::REC4, eglif_dcnp::PORT_NOT_AVAILABLE },
};

// ---------------------------------------------------------------------------
//   Default constructors defining default parameters and state
//   Note: the implementation is empty. The initialization is of variables
//   is a part of eglif_dcnp's constructor.
// ---------------------------------------------------------------------------

eglif_dcnp::Parameters_::Parameters_()
{
}

eglif_dcnp::State_::State_()
{
}

// ---------------------------------------------------------------------------
//   Parameter and state extractions and manipulation functions
// ---------------------------------------------------------------------------

eglif_dcnp::Buffers_::Buffers_(eglif_dcnp &n):
  logger_(n)
  , spike_inputs_( std::vector< nest::RingBuffer >( NUM_SPIKE_RECEPTORS ) )
  , spike_inputs_grid_sum_( std::vector< double >( NUM_SPIKE_RECEPTORS ) )
  , spike_input_received_( std::vector< nest::RingBuffer >( NUM_SPIKE_RECEPTORS ) )
  , spike_input_received_grid_sum_( std::vector< double >( NUM_SPIKE_RECEPTORS ) )
  , __s( nullptr ), __c( nullptr ), __e( nullptr )
{
  // Initialization of the remaining members is deferred to init_buffers_().
}

eglif_dcnp::Buffers_::Buffers_(const Buffers_ &, eglif_dcnp &n):
  logger_(n)
  , spike_inputs_( std::vector< nest::RingBuffer >( NUM_SPIKE_RECEPTORS ) )
  , spike_inputs_grid_sum_( std::vector< double >( NUM_SPIKE_RECEPTORS ) )
  , spike_input_received_( std::vector< nest::RingBuffer >( NUM_SPIKE_RECEPTORS ) )
  , spike_input_received_grid_sum_( std::vector< double >( NUM_SPIKE_RECEPTORS ) )
  , __s( nullptr ), __c( nullptr ), __e( nullptr )
{
  // Initialization of the remaining members is deferred to init_buffers_().
}

// ---------------------------------------------------------------------------
//   Default constructor for node
// ---------------------------------------------------------------------------

eglif_dcnp::eglif_dcnp():ExtendedPostHistoryArchivingNode(), P_(), S_(), B_(*this)
{
  init_state_internal_();
  recordablesMap_.create();
  pre_run_hook();
}

// ---------------------------------------------------------------------------
//   Copy constructor for node
// ---------------------------------------------------------------------------

eglif_dcnp::eglif_dcnp(const eglif_dcnp& __n):
  ExtendedPostHistoryArchivingNode(), P_(__n.P_), S_(__n.S_), B_(__n.B_, *this)
{

  // copy parameter struct P_
  P_.C_m = __n.P_.C_m;
  P_.tau_m = __n.P_.tau_m;
  P_.E_L = __n.P_.E_L;
  P_.t_ref = __n.P_.t_ref;
  P_.V_reset = __n.P_.V_reset;
  P_.V_th = __n.P_.V_th;
  P_.Vmin = __n.P_.Vmin;
  P_.I_e = __n.P_.I_e;
  P_.Vinit = __n.P_.Vinit;
  P_.lambda_0 = __n.P_.lambda_0;
  P_.tau_V = __n.P_.tau_V;
  P_.kadap = __n.P_.kadap;
  P_.k2 = __n.P_.k2;
  P_.k1 = __n.P_.k1;
  P_.A1 = __n.P_.A1;
  P_.A2 = __n.P_.A2;
  P_.E_rev1 = __n.P_.E_rev1;
  P_.E_rev2 = __n.P_.E_rev2;
  P_.E_rev3 = __n.P_.E_rev3;
  P_.E_rev4 = __n.P_.E_rev4;
  P_.tau_syn1 = __n.P_.tau_syn1;
  P_.tau_syn2 = __n.P_.tau_syn2;
  P_.tau_syn3 = __n.P_.tau_syn3;
  P_.tau_syn4 = __n.P_.tau_syn4;
  P_.offset = __n.P_.offset;
  // copy state struct S_
  S_.ode_state[State_::V_m] = __n.S_.ode_state[State_::V_m];
  S_.ode_state[State_::I_adap] = __n.S_.ode_state[State_::I_adap];
  S_.ode_state[State_::I_dep] = __n.S_.ode_state[State_::I_dep];
  S_.r = __n.S_.r;
  S_.lambda = __n.S_.lambda;
  S_.complex_flag = __n.S_.complex_flag;
  S_.tick = __n.S_.tick;
  S_.cf_buffer = __n.S_.cf_buffer;
  S_.gr_buffer = __n.S_.gr_buffer;
  S_.last_io = __n.S_.last_io;
  S_.ode_state[State_::g1__X__rec1] = __n.S_.ode_state[State_::g1__X__rec1];
  S_.ode_state[State_::g1__X__rec1__d] = __n.S_.ode_state[State_::g1__X__rec1__d];
  S_.ode_state[State_::g4__X__rec4] = __n.S_.ode_state[State_::g4__X__rec4];
  S_.ode_state[State_::g4__X__rec4__d] = __n.S_.ode_state[State_::g4__X__rec4__d];
  S_.ode_state[State_::g2__X__rec2] = __n.S_.ode_state[State_::g2__X__rec2];
  S_.ode_state[State_::g2__X__rec2__d] = __n.S_.ode_state[State_::g2__X__rec2__d];
  S_.ode_state[State_::g3__X__rec3] = __n.S_.ode_state[State_::g3__X__rec3];
  S_.ode_state[State_::g3__X__rec3__d] = __n.S_.ode_state[State_::g3__X__rec3__d];

  // copy internals V_
  V_.RefractoryCounts = __n.V_.RefractoryCounts;
  V_.__h = __n.V_.__h;
  V_.__P__I_dep__I_dep = __n.V_.__P__I_dep__I_dep;
  V_.__P__g1__X__rec1__g1__X__rec1 = __n.V_.__P__g1__X__rec1__g1__X__rec1;
  V_.__P__g1__X__rec1__g1__X__rec1__d = __n.V_.__P__g1__X__rec1__g1__X__rec1__d;
  V_.__P__g1__X__rec1__d__g1__X__rec1 = __n.V_.__P__g1__X__rec1__d__g1__X__rec1;
  V_.__P__g1__X__rec1__d__g1__X__rec1__d = __n.V_.__P__g1__X__rec1__d__g1__X__rec1__d;
  V_.__P__g4__X__rec4__g4__X__rec4 = __n.V_.__P__g4__X__rec4__g4__X__rec4;
  V_.__P__g4__X__rec4__g4__X__rec4__d = __n.V_.__P__g4__X__rec4__g4__X__rec4__d;
  V_.__P__g4__X__rec4__d__g4__X__rec4 = __n.V_.__P__g4__X__rec4__d__g4__X__rec4;
  V_.__P__g4__X__rec4__d__g4__X__rec4__d = __n.V_.__P__g4__X__rec4__d__g4__X__rec4__d;
  V_.__P__g2__X__rec2__g2__X__rec2 = __n.V_.__P__g2__X__rec2__g2__X__rec2;
  V_.__P__g2__X__rec2__g2__X__rec2__d = __n.V_.__P__g2__X__rec2__g2__X__rec2__d;
  V_.__P__g2__X__rec2__d__g2__X__rec2 = __n.V_.__P__g2__X__rec2__d__g2__X__rec2;
  V_.__P__g2__X__rec2__d__g2__X__rec2__d = __n.V_.__P__g2__X__rec2__d__g2__X__rec2__d;
  V_.__P__g3__X__rec3__g3__X__rec3 = __n.V_.__P__g3__X__rec3__g3__X__rec3;
  V_.__P__g3__X__rec3__g3__X__rec3__d = __n.V_.__P__g3__X__rec3__g3__X__rec3__d;
  V_.__P__g3__X__rec3__d__g3__X__rec3 = __n.V_.__P__g3__X__rec3__d__g3__X__rec3;
  V_.__P__g3__X__rec3__d__g3__X__rec3__d = __n.V_.__P__g3__X__rec3__d__g3__X__rec3__d;
}

// ---------------------------------------------------------------------------
//   Destructor for node
// ---------------------------------------------------------------------------

eglif_dcnp::~eglif_dcnp()
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
void eglif_dcnp::calibrate_time( const nest::TimeConverter& tc )
{
  LOG( nest::M_WARNING,
    "eglif_dcnp",
    "Simulation resolution has changed. Internal state and parameters of the model have been reset!" );

  init_state_internal_();
}
void eglif_dcnp::init_state_internal_()
{

  const double __resolution = nest::Time::get_resolution().get_ms();  // do not remove, this is necessary for the resolution() function
  // by default, integrate all variables
  // use a default "good enough" value for the absolute error. It can be adjusted via `node.set()`
  P_.__gsl_error_tol = 1e-3;
  // initial values for parameters
  P_.C_m = 14.6; // as pF
  P_.tau_m = 9.125; // as ms
  P_.E_L = (-68.0); // as mV
  P_.t_ref = 1.59; // as ms
  P_.V_reset = (-78.0); // as mV
  P_.V_th = (-53.0); // as mV
  P_.Vmin = (-150.0); // as mV
  P_.I_e = 3.711; // as pA
  P_.Vinit = (-60.0); // as mV
  P_.lambda_0 = 1.8 / 1; // as 1 / ms
  P_.tau_V = 1.1; // as mV
  P_.kadap = 2.025 / (1.0 * 1.0); // as pA / (mV ms)
  P_.k2 = 1.096 / 1; // as 1 / ms
  P_.k1 = 1.887 / 1; // as 1 / ms
  P_.A1 = 5.953; // as pA
  P_.A2 = 5.863; // as pA
  P_.E_rev1 = 0.0; // as mV
  P_.E_rev2 = (-80.0); // as mV
  P_.E_rev3 = 0.0; // as mV
  P_.E_rev4 = (-80.0); // as mV
  P_.tau_syn1 = 0.2; // as ms
  P_.tau_syn2 = 2.0; // as ms
  P_.tau_syn3 = 2.0; // as ms
  P_.tau_syn4 = 2.0; // as ms
  P_.offset = 0; // as integer

  V_.__h = nest::Time::get_resolution().get_ms();
  recompute_internal_variables();
  // initial values for state variables
  S_.ode_state[State_::V_m] = P_.Vinit; // as mV
  S_.ode_state[State_::I_adap] = 0.0; // as pA
  S_.ode_state[State_::I_dep] = 0.0; // as pA
  S_.r = 0; // as integer
  S_.lambda = pow(0, (-1)); // as 1 / ms
  S_.complex_flag = 0.0;
  S_.tick = 0; // as ms
  S_.cf_buffer = 0; // as real
  S_.gr_buffer = 0; // as real
  S_.last_io = 1000; // as ms
  S_.ode_state[State_::g1__X__rec1] = 0; // as real
  S_.ode_state[State_::g1__X__rec1__d] = 0; // as 1 / s
  S_.ode_state[State_::g4__X__rec4] = 0; // as real
  S_.ode_state[State_::g4__X__rec4__d] = 0; // as 1 / s
  S_.ode_state[State_::g2__X__rec2] = 0; // as real
  S_.ode_state[State_::g2__X__rec2__d] = 0; // as 1 / s
  S_.ode_state[State_::g3__X__rec3] = 0; // as real
  S_.ode_state[State_::g3__X__rec3__d] = 0; // as 1 / s
}

void eglif_dcnp::init_buffers_()
{
#ifdef DEBUG
  std::cout << "eglif_dcnp::init_buffers_()" << std::endl;
#endif
  // spike input buffers
  get_spike_inputs_().clear();
  get_spike_inputs_grid_sum_().clear();
  get_spike_input_received_().clear();
  get_spike_input_received_grid_sum_().clear();


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

  // B_.__sys.function = eglif_dcnp_dynamics; // will be set just prior to the call to gsl_odeiv_evolve_apply()
  B_.__sys.jacobian = nullptr;
  B_.__sys.dimension = State_::STATE_VEC_SIZE;
  B_.__sys.params = reinterpret_cast< void* >( this );
  B_.__step = nest::Time::get_resolution().get_ms();
  B_.__integration_step = nest::Time::get_resolution().get_ms();
}

void eglif_dcnp::recompute_internal_variables(bool exclude_timestep)
{
  const double __resolution = nest::Time::get_resolution().get_ms();  // do not remove, this is necessary for the resolution() function

  if (exclude_timestep)
  {    
    V_.RefractoryCounts = nest::Time(nest::Time::ms((double) (P_.t_ref))).get_steps(); // as integer
    V_.__P__I_dep__I_dep = 1.0 * std::exp((-V_.__h) * P_.k1); // as real
    V_.__P__g1__X__rec1__g1__X__rec1 = 1.0 * (V_.__h + P_.tau_syn1) * std::exp((-V_.__h) / P_.tau_syn1) / P_.tau_syn1; // as real
    V_.__P__g1__X__rec1__g1__X__rec1__d = 1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn1); // as real
    V_.__P__g1__X__rec1__d__g1__X__rec1 = (-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn1) / pow(P_.tau_syn1, 2); // as real
    V_.__P__g1__X__rec1__d__g1__X__rec1__d = 1.0 * ((-V_.__h) + P_.tau_syn1) * std::exp((-V_.__h) / P_.tau_syn1) / P_.tau_syn1; // as real
    V_.__P__g4__X__rec4__g4__X__rec4 = 1.0 * (V_.__h + P_.tau_syn4) * std::exp((-V_.__h) / P_.tau_syn4) / P_.tau_syn4; // as real
    V_.__P__g4__X__rec4__g4__X__rec4__d = 1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn4); // as real
    V_.__P__g4__X__rec4__d__g4__X__rec4 = (-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn4) / pow(P_.tau_syn4, 2); // as real
    V_.__P__g4__X__rec4__d__g4__X__rec4__d = 1.0 * ((-V_.__h) + P_.tau_syn4) * std::exp((-V_.__h) / P_.tau_syn4) / P_.tau_syn4; // as real
    V_.__P__g2__X__rec2__g2__X__rec2 = 1.0 * (V_.__h + P_.tau_syn2) * std::exp((-V_.__h) / P_.tau_syn2) / P_.tau_syn2; // as real
    V_.__P__g2__X__rec2__g2__X__rec2__d = 1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn2); // as real
    V_.__P__g2__X__rec2__d__g2__X__rec2 = (-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn2) / pow(P_.tau_syn2, 2); // as real
    V_.__P__g2__X__rec2__d__g2__X__rec2__d = 1.0 * ((-V_.__h) + P_.tau_syn2) * std::exp((-V_.__h) / P_.tau_syn2) / P_.tau_syn2; // as real
    V_.__P__g3__X__rec3__g3__X__rec3 = 1.0 * (V_.__h + P_.tau_syn3) * std::exp((-V_.__h) / P_.tau_syn3) / P_.tau_syn3; // as real
    V_.__P__g3__X__rec3__g3__X__rec3__d = 1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn3); // as real
    V_.__P__g3__X__rec3__d__g3__X__rec3 = (-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn3) / pow(P_.tau_syn3, 2); // as real
    V_.__P__g3__X__rec3__d__g3__X__rec3__d = 1.0 * ((-V_.__h) + P_.tau_syn3) * std::exp((-V_.__h) / P_.tau_syn3) / P_.tau_syn3; // as real
  }
  else {    
    V_.RefractoryCounts = nest::Time(nest::Time::ms((double) (P_.t_ref))).get_steps(); // as integer
    V_.__h = __resolution; // as ms
    V_.__P__I_dep__I_dep = 1.0 * std::exp((-V_.__h) * P_.k1); // as real
    V_.__P__g1__X__rec1__g1__X__rec1 = 1.0 * (V_.__h + P_.tau_syn1) * std::exp((-V_.__h) / P_.tau_syn1) / P_.tau_syn1; // as real
    V_.__P__g1__X__rec1__g1__X__rec1__d = 1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn1); // as real
    V_.__P__g1__X__rec1__d__g1__X__rec1 = (-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn1) / pow(P_.tau_syn1, 2); // as real
    V_.__P__g1__X__rec1__d__g1__X__rec1__d = 1.0 * ((-V_.__h) + P_.tau_syn1) * std::exp((-V_.__h) / P_.tau_syn1) / P_.tau_syn1; // as real
    V_.__P__g4__X__rec4__g4__X__rec4 = 1.0 * (V_.__h + P_.tau_syn4) * std::exp((-V_.__h) / P_.tau_syn4) / P_.tau_syn4; // as real
    V_.__P__g4__X__rec4__g4__X__rec4__d = 1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn4); // as real
    V_.__P__g4__X__rec4__d__g4__X__rec4 = (-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn4) / pow(P_.tau_syn4, 2); // as real
    V_.__P__g4__X__rec4__d__g4__X__rec4__d = 1.0 * ((-V_.__h) + P_.tau_syn4) * std::exp((-V_.__h) / P_.tau_syn4) / P_.tau_syn4; // as real
    V_.__P__g2__X__rec2__g2__X__rec2 = 1.0 * (V_.__h + P_.tau_syn2) * std::exp((-V_.__h) / P_.tau_syn2) / P_.tau_syn2; // as real
    V_.__P__g2__X__rec2__g2__X__rec2__d = 1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn2); // as real
    V_.__P__g2__X__rec2__d__g2__X__rec2 = (-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn2) / pow(P_.tau_syn2, 2); // as real
    V_.__P__g2__X__rec2__d__g2__X__rec2__d = 1.0 * ((-V_.__h) + P_.tau_syn2) * std::exp((-V_.__h) / P_.tau_syn2) / P_.tau_syn2; // as real
    V_.__P__g3__X__rec3__g3__X__rec3 = 1.0 * (V_.__h + P_.tau_syn3) * std::exp((-V_.__h) / P_.tau_syn3) / P_.tau_syn3; // as real
    V_.__P__g3__X__rec3__g3__X__rec3__d = 1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn3); // as real
    V_.__P__g3__X__rec3__d__g3__X__rec3 = (-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn3) / pow(P_.tau_syn3, 2); // as real
    V_.__P__g3__X__rec3__d__g3__X__rec3__d = 1.0 * ((-V_.__h) + P_.tau_syn3) * std::exp((-V_.__h) / P_.tau_syn3) / P_.tau_syn3; // as real
  }
}
void eglif_dcnp::pre_run_hook()
{
  B_.logger_.init();

  // parameters might have changed -- recompute internals
  V_.__h = nest::Time::get_resolution().get_ms();
  recompute_internal_variables();

  // buffers B_
  B_.spike_inputs_.resize(NUM_SPIKE_RECEPTORS);
  B_.spike_inputs_grid_sum_.resize(NUM_SPIKE_RECEPTORS);
  B_.spike_input_received_.resize(NUM_SPIKE_RECEPTORS);
  B_.spike_input_received_grid_sum_.resize(NUM_SPIKE_RECEPTORS);
}

// ---------------------------------------------------------------------------
//   Update and spike handling functions
// ---------------------------------------------------------------------------

extern "C" inline int eglif_dcnp_dynamics(double __time, const double ode_state[], double f[], void* pnode)
{
  typedef eglif_dcnp::State_ State_;
  // get access to node so we can almost work as in a member function
  assert( pnode );
  const eglif_dcnp& node = *( reinterpret_cast< eglif_dcnp* >( pnode ) );

  // ode_state[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.ode_state[].
  f[State_::V_m] = ((-node.P_.E_L) + std::max(ode_state[State_::V_m], node.P_.Vmin)) / node.P_.tau_m + ((-ode_state[State_::I_adap]) + ode_state[State_::I_dep] + node.P_.I_e + node.B_.I_stim_grid_sum_ + 1.0 * ode_state[State_::g1__X__rec1] * (node.P_.E_rev1 - ode_state[State_::V_m]) + 1.0 * ode_state[State_::g2__X__rec2] * (node.P_.E_rev2 - ode_state[State_::V_m]) + 1.0 * ode_state[State_::g3__X__rec3] * (node.P_.E_rev3 - ode_state[State_::V_m]) + 1.0 * ode_state[State_::g4__X__rec4] * (node.P_.E_rev4 - ode_state[State_::V_m])) / node.P_.C_m;
  f[State_::I_dep] = (-ode_state[State_::I_dep]) * node.P_.k1;
  f[State_::I_adap] = (-ode_state[State_::I_adap]) * node.P_.k2 + node.P_.kadap * ((-node.P_.E_L) + std::max(ode_state[State_::V_m], node.P_.Vmin));
  f[State_::g1__X__rec1] = 1.0 * ode_state[State_::g1__X__rec1__d];
  f[State_::g1__X__rec1__d] = (-ode_state[State_::g1__X__rec1]) / pow(node.P_.tau_syn1, 2) - 2 * ode_state[State_::g1__X__rec1__d] / node.P_.tau_syn1;
  f[State_::g4__X__rec4] = 1.0 * ode_state[State_::g4__X__rec4__d];
  f[State_::g4__X__rec4__d] = (-ode_state[State_::g4__X__rec4]) / pow(node.P_.tau_syn4, 2) - 2 * ode_state[State_::g4__X__rec4__d] / node.P_.tau_syn4;
  f[State_::g2__X__rec2] = 1.0 * ode_state[State_::g2__X__rec2__d];
  f[State_::g2__X__rec2__d] = (-ode_state[State_::g2__X__rec2]) / pow(node.P_.tau_syn2, 2) - 2 * ode_state[State_::g2__X__rec2__d] / node.P_.tau_syn2;
  f[State_::g3__X__rec3] = 1.0 * ode_state[State_::g3__X__rec3__d];
  f[State_::g3__X__rec3__d] = (-ode_state[State_::g3__X__rec3]) / pow(node.P_.tau_syn3, 2) - 2 * ode_state[State_::g3__X__rec3__d] / node.P_.tau_syn3;
  return GSL_SUCCESS;
}
void eglif_dcnp::update(nest::Time const & origin,const long from, const long to)
{
  const double __resolution = nest::Time::get_resolution().get_ms();  // do not remove, this is necessary for the resolution() function

  for ( long lag = from ; lag < to ; ++lag )
  {
    auto get_t = [origin, lag](){ return nest::Time( nest::Time::step( origin.get_steps() + lag + 1) ).get_ms(); };
    /**
     * buffer spikes from spiking input ports
    **/

    for (long i = 0; i < NUM_SPIKE_RECEPTORS; ++i)
    {
      get_spike_inputs_grid_sum_()[i] = get_spike_inputs_()[i].get_value(lag);
      get_spike_input_received_grid_sum_()[i] = get_spike_input_received_()[i].get_value(lag);

    }

    /**
     * subthreshold updates of the convolution variables
     *
     * step 1: regardless of whether and how integrate_odes() will be called, update variables due to convolutions
    **/
const double g1__X__rec1__tmp_ = V_.__P__g1__X__rec1__g1__X__rec1 * S_.ode_state[State_::g1__X__rec1] + V_.__P__g1__X__rec1__g1__X__rec1__d * S_.ode_state[State_::g1__X__rec1__d];
const double g1__X__rec1__d__tmp_ = V_.__P__g1__X__rec1__d__g1__X__rec1 * S_.ode_state[State_::g1__X__rec1] + V_.__P__g1__X__rec1__d__g1__X__rec1__d * S_.ode_state[State_::g1__X__rec1__d];
const double g4__X__rec4__tmp_ = V_.__P__g4__X__rec4__g4__X__rec4 * S_.ode_state[State_::g4__X__rec4] + V_.__P__g4__X__rec4__g4__X__rec4__d * S_.ode_state[State_::g4__X__rec4__d];
const double g4__X__rec4__d__tmp_ = V_.__P__g4__X__rec4__d__g4__X__rec4 * S_.ode_state[State_::g4__X__rec4] + V_.__P__g4__X__rec4__d__g4__X__rec4__d * S_.ode_state[State_::g4__X__rec4__d];
const double g2__X__rec2__tmp_ = V_.__P__g2__X__rec2__g2__X__rec2 * S_.ode_state[State_::g2__X__rec2] + V_.__P__g2__X__rec2__g2__X__rec2__d * S_.ode_state[State_::g2__X__rec2__d];
const double g2__X__rec2__d__tmp_ = V_.__P__g2__X__rec2__d__g2__X__rec2 * S_.ode_state[State_::g2__X__rec2] + V_.__P__g2__X__rec2__d__g2__X__rec2__d * S_.ode_state[State_::g2__X__rec2__d];
const double g3__X__rec3__tmp_ = V_.__P__g3__X__rec3__g3__X__rec3 * S_.ode_state[State_::g3__X__rec3] + V_.__P__g3__X__rec3__g3__X__rec3__d * S_.ode_state[State_::g3__X__rec3__d];
const double g3__X__rec3__d__tmp_ = V_.__P__g3__X__rec3__d__g3__X__rec3 * S_.ode_state[State_::g3__X__rec3] + V_.__P__g3__X__rec3__d__g3__X__rec3__d * S_.ode_state[State_::g3__X__rec3__d];


    /**
     * Begin NESTML generated code for the update block(s)
    **/

    S_.tick = get_t();
    // in this condition, if there is an PC spike received from the receptor 2,
    // it will trigger a "complex flag" change, utilized for the depression in the MF-DCNp synapses
    S_.cf_buffer = (0.001 * B_.spike_inputs_grid_sum_[REC2 - MIN_SPIKE_RECEPTOR]);
    if (S_.cf_buffer != 0)
    {
        set_spiketime(nest::Time::step(origin.get_steps() + lag), 1);
        nest::SpikeEvent se;
        nest::kernel().event_delivery_manager.send(*this, se, lag);
    }
    if (S_.r == 0)
    {  

        // start rendered code for integrate_odes()

        // analytic solver: integrating state variables (first step): I_dep, g1__X__rec1, g1__X__rec1__d, g4__X__rec4, g4__X__rec4__d, g2__X__rec2, g2__X__rec2__d, g3__X__rec3, g3__X__rec3__d, 
        const double I_dep__tmp = S_.ode_state[State_::I_dep] * V_.__P__I_dep__I_dep;
        const double g1__X__rec1__tmp = V_.__P__g1__X__rec1__g1__X__rec1 * S_.ode_state[State_::g1__X__rec1] + V_.__P__g1__X__rec1__g1__X__rec1__d * S_.ode_state[State_::g1__X__rec1__d];
        const double g1__X__rec1__d__tmp = V_.__P__g1__X__rec1__d__g1__X__rec1 * S_.ode_state[State_::g1__X__rec1] + V_.__P__g1__X__rec1__d__g1__X__rec1__d * S_.ode_state[State_::g1__X__rec1__d];
        const double g4__X__rec4__tmp = V_.__P__g4__X__rec4__g4__X__rec4 * S_.ode_state[State_::g4__X__rec4] + V_.__P__g4__X__rec4__g4__X__rec4__d * S_.ode_state[State_::g4__X__rec4__d];
        const double g4__X__rec4__d__tmp = V_.__P__g4__X__rec4__d__g4__X__rec4 * S_.ode_state[State_::g4__X__rec4] + V_.__P__g4__X__rec4__d__g4__X__rec4__d * S_.ode_state[State_::g4__X__rec4__d];
        const double g2__X__rec2__tmp = V_.__P__g2__X__rec2__g2__X__rec2 * S_.ode_state[State_::g2__X__rec2] + V_.__P__g2__X__rec2__g2__X__rec2__d * S_.ode_state[State_::g2__X__rec2__d];
        const double g2__X__rec2__d__tmp = V_.__P__g2__X__rec2__d__g2__X__rec2 * S_.ode_state[State_::g2__X__rec2] + V_.__P__g2__X__rec2__d__g2__X__rec2__d * S_.ode_state[State_::g2__X__rec2__d];
        const double g3__X__rec3__tmp = V_.__P__g3__X__rec3__g3__X__rec3 * S_.ode_state[State_::g3__X__rec3] + V_.__P__g3__X__rec3__g3__X__rec3__d * S_.ode_state[State_::g3__X__rec3__d];
        const double g3__X__rec3__d__tmp = V_.__P__g3__X__rec3__d__g3__X__rec3 * S_.ode_state[State_::g3__X__rec3] + V_.__P__g3__X__rec3__d__g3__X__rec3__d * S_.ode_state[State_::g3__X__rec3__d];

         
        // numeric solver: integrating state variables: V_m, I_dep, I_adap, g1__X__rec1, g1__X__rec1__d, g4__X__rec4, g4__X__rec4__d, g2__X__rec2, g2__X__rec2__d, g3__X__rec3, g3__X__rec3__d, 
        double __t = 0;
        B_.__sys.function = eglif_dcnp_dynamics;
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
        // analytic solver: integrating state variables (second step): I_dep, g1__X__rec1, g1__X__rec1__d, g4__X__rec4, g4__X__rec4__d, g2__X__rec2, g2__X__rec2__d, g3__X__rec3, g3__X__rec3__d, 
        /* replace analytically solvable variables with precisely integrated values  */
        S_.ode_state[State_::I_dep] = I_dep__tmp;
        S_.ode_state[State_::g1__X__rec1] = g1__X__rec1__tmp;
        S_.ode_state[State_::g1__X__rec1__d] = g1__X__rec1__d__tmp;
        S_.ode_state[State_::g4__X__rec4] = g4__X__rec4__tmp;
        S_.ode_state[State_::g4__X__rec4__d] = g4__X__rec4__d__tmp;
        S_.ode_state[State_::g2__X__rec2] = g2__X__rec2__tmp;
        S_.ode_state[State_::g2__X__rec2__d] = g2__X__rec2__d__tmp;
        S_.ode_state[State_::g3__X__rec3] = g3__X__rec3__tmp;
        S_.ode_state[State_::g3__X__rec3__d] = g3__X__rec3__d__tmp;
    }
    else if (S_.r > 0)
    {  
        S_.r -= 1;
    }
    S_.lambda = P_.lambda_0 * std::exp((S_.ode_state[State_::V_m] - P_.V_th) / P_.tau_V);
    if (S_.lambda > 0 / 1.0)
    {  
        double rnd = ((0) + (1) * nest::get_vp_specific_rng( get_thread() )->drand());
        double thr = 1 - std::exp((-S_.lambda) * __resolution);
        if (rnd < thr)
        {  
            S_.r = V_.RefractoryCounts;
            S_.ode_state[State_::V_m] = P_.V_reset;
            S_.ode_state[State_::I_adap] += P_.A2;
            S_.ode_state[State_::I_dep] = P_.A1;

            /**
             * generated code for emit_spike() function
            **/

            set_spiketime(nest::Time::step(origin.get_steps() + lag + 1));
            nest::SpikeEvent se;
            //std::cout << "PC: Emitting simple spike.\n";
            nest::kernel().event_delivery_manager.send(*this, se, lag);
        }
    }

    /**
     * Begin NESTML generated code for the onReceive block(s)
    **/

    /**
     * subthreshold updates of the convolution variables
     *
     * step 2: regardless of whether and how integrate_odes() was called, update variables due to convolutions. Set to the updated values at the end of the timestep.
    **/

    S_.ode_state[State_::g1__X__rec1] = g1__X__rec1__tmp_;
    S_.ode_state[State_::g1__X__rec1__d] = g1__X__rec1__d__tmp_;
    S_.ode_state[State_::g4__X__rec4] = g4__X__rec4__tmp_;
    S_.ode_state[State_::g4__X__rec4__d] = g4__X__rec4__d__tmp_;
    S_.ode_state[State_::g2__X__rec2] = g2__X__rec2__tmp_;
    S_.ode_state[State_::g2__X__rec2__d] = g2__X__rec2__d__tmp_;
    S_.ode_state[State_::g3__X__rec3] = g3__X__rec3__tmp_;
    S_.ode_state[State_::g3__X__rec3__d] = g3__X__rec3__d__tmp_;


    /**
     * spike updates due to convolutions
    **/

    S_.ode_state[State_::g1__X__rec1__d] += ((0.001 * B_.spike_inputs_grid_sum_[REC1 - MIN_SPIKE_RECEPTOR])) * (numerics::e / P_.tau_syn1) / (1 / 1000.0);
    S_.ode_state[State_::g4__X__rec4__d] += ((0.001 * B_.spike_inputs_grid_sum_[REC4 - MIN_SPIKE_RECEPTOR])) * (numerics::e / P_.tau_syn4) / (1 / 1000.0);
    S_.ode_state[State_::g2__X__rec2__d] += ((0.001 * B_.spike_inputs_grid_sum_[REC2 - MIN_SPIKE_RECEPTOR])) * (numerics::e / P_.tau_syn2) / (1 / 1000.0);
    S_.ode_state[State_::g3__X__rec3__d] += ((0.001 * B_.spike_inputs_grid_sum_[REC3 - MIN_SPIKE_RECEPTOR])) * (numerics::e / P_.tau_syn3) / (1 / 1000.0);

    /**
     * Begin NESTML generated code for the onCondition block(s)
    **/


    /**
     * handle continuous input ports
    **/
    B_.I_stim_grid_sum_ = get_I_stim().get_value(lag);
    // voltage logging
    B_.logger_.record_data(origin.get_steps() + lag);
  }
}

// Do not move this function as inline to h-file. It depends on
// universal_data_logger_impl.h being included here.
void eglif_dcnp::handle(nest::DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}


void eglif_dcnp::handle(nest::SpikeEvent &e)
{
  assert(e.get_delay_steps() > 0);
  assert( e.get_rport() < B_.spike_inputs_.size() );

  double weight = e.get_weight();
  size_t nestml_buffer_idx = 0;
  if ( weight > 0.0 )
  {
    nestml_buffer_idx = std::get<0>(rport_to_nestml_buffer_idx[e.get_rport()]);
  }
  else if ( weight < 0.0 )
  {
    nestml_buffer_idx = std::get<1>(rport_to_nestml_buffer_idx[e.get_rport()]);
    if ( nestml_buffer_idx == eglif_dcnp::PORT_NOT_AVAILABLE )
    {
      nestml_buffer_idx = std::get<0>(rport_to_nestml_buffer_idx[e.get_rport()]);
    }
    weight = -weight;
  }
  if ( weight != 0.0 ){
      B_.spike_inputs_[ nestml_buffer_idx - MIN_SPIKE_RECEPTOR ].add_value(
        e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin() ),
        weight * e.get_multiplicity() );
      B_.spike_input_received_[ nestml_buffer_idx - MIN_SPIKE_RECEPTOR ].add_value(
        e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin() ),
        1. );
  }
}

void eglif_dcnp::handle(nest::CurrentEvent& e)
{
  assert(e.get_delay_steps() > 0);

  const double current = e.get_current();     // we assume that in NEST, this returns a current in pA
  const double weight = e.get_weight();
  get_I_stim().add_value(
               e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
               weight * current );
}



