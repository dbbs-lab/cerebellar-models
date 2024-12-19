
/*
 *  iaf_cond_alpha_neuron_nestml.cpp
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
 *  Generated from NESTML at time: 2024-11-05 14:48:36.475194
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

#include "iaf_cond_alpha_neuron_nestml.h"

// uncomment the next line to enable printing of detailed debug information
#define DEBUG


void
register_iaf_cond_alpha_neuron_nestml( const std::string& name )
{
  nest::register_node_model< iaf_cond_alpha_neuron_nestml >( name );
}

// ---------------------------------------------------------------------------
//   Recordables map
// ---------------------------------------------------------------------------
nest::RecordablesMap<iaf_cond_alpha_neuron_nestml> iaf_cond_alpha_neuron_nestml::recordablesMap_;
namespace nest
{

  // Override the create() method with one call to RecordablesMap::insert_() for each quantity to be recorded.
  template <> void RecordablesMap<iaf_cond_alpha_neuron_nestml>::create()
  {
    // add state variables to recordables map
    insert_(iaf_cond_alpha_neuron_nestml_names::_V_m, &iaf_cond_alpha_neuron_nestml::get_V_m);
    insert_(iaf_cond_alpha_neuron_nestml_names::_refr_t, &iaf_cond_alpha_neuron_nestml::get_refr_t);
    insert_(iaf_cond_alpha_neuron_nestml_names::_g_inh__X__inh_spikes__DOT__weight, &iaf_cond_alpha_neuron_nestml::get_g_inh__X__inh_spikes__DOT__weight);
    insert_(iaf_cond_alpha_neuron_nestml_names::_g_inh__X__inh_spikes__DOT__weight__d, &iaf_cond_alpha_neuron_nestml::get_g_inh__X__inh_spikes__DOT__weight__d);
    insert_(iaf_cond_alpha_neuron_nestml_names::_g_exc__X__exc_spikes__DOT__weight, &iaf_cond_alpha_neuron_nestml::get_g_exc__X__exc_spikes__DOT__weight);
    insert_(iaf_cond_alpha_neuron_nestml_names::_g_exc__X__exc_spikes__DOT__weight__d, &iaf_cond_alpha_neuron_nestml::get_g_exc__X__exc_spikes__DOT__weight__d);

    // Add vector variables
  }
}

// ---------------------------------------------------------------------------
//   Default constructors defining default parameters and state
//   Note: the implementation is empty. The initialization is of variables
//   is a part of iaf_cond_alpha_neuron_nestml's constructor.
// ---------------------------------------------------------------------------

iaf_cond_alpha_neuron_nestml::Parameters_::Parameters_()
{
}

iaf_cond_alpha_neuron_nestml::State_::State_()
{
}

// ---------------------------------------------------------------------------
//   Parameter and state extractions and manipulation functions
// ---------------------------------------------------------------------------

iaf_cond_alpha_neuron_nestml::Buffers_::Buffers_(iaf_cond_alpha_neuron_nestml &n):
  logger_(n)
  , spike_input_exc_spikes__DOT__weight_( nest::RingBuffer() )
  , spike_input_exc_spikes__DOT__weight_grid_sum_( 0. )
  , spike_input_exc_spikes_spike_input_received_( nest::RingBuffer() )
  , spike_input_exc_spikes_spike_input_received_grid_sum_( 0. )
  , spike_input_inh_spikes__DOT__weight_( nest::RingBuffer() )
  , spike_input_inh_spikes__DOT__weight_grid_sum_( 0. )
  , spike_input_inh_spikes_spike_input_received_( nest::RingBuffer() )
  , spike_input_inh_spikes_spike_input_received_grid_sum_( 0. )
  , __s( nullptr ), __c( nullptr ), __e( nullptr )
{
  // Initialization of the remaining members is deferred to init_buffers_().
}

iaf_cond_alpha_neuron_nestml::Buffers_::Buffers_(const Buffers_ &, iaf_cond_alpha_neuron_nestml &n):
  logger_(n)
  , spike_input_exc_spikes__DOT__weight_( nest::RingBuffer() )
  , spike_input_exc_spikes__DOT__weight_grid_sum_( 0. )
  , spike_input_exc_spikes_spike_input_received_( nest::RingBuffer() )
  , spike_input_exc_spikes_spike_input_received_grid_sum_( 0. )
  , spike_input_inh_spikes__DOT__weight_( nest::RingBuffer() )
  , spike_input_inh_spikes__DOT__weight_grid_sum_( 0. )
  , spike_input_inh_spikes_spike_input_received_( nest::RingBuffer() )
  , spike_input_inh_spikes_spike_input_received_grid_sum_( 0. )
  , __s( nullptr ), __c( nullptr ), __e( nullptr )
{
  // Initialization of the remaining members is deferred to init_buffers_().
}

// ---------------------------------------------------------------------------
//   Default constructor for node
// ---------------------------------------------------------------------------

iaf_cond_alpha_neuron_nestml::iaf_cond_alpha_neuron_nestml():ExtendedPostHistoryArchivingNode(), P_(), S_(), B_(*this)
{
  init_state_internal_();
  recordablesMap_.create();
  pre_run_hook();
}

// ---------------------------------------------------------------------------
//   Copy constructor for node
// ---------------------------------------------------------------------------

iaf_cond_alpha_neuron_nestml::iaf_cond_alpha_neuron_nestml(const iaf_cond_alpha_neuron_nestml& __n):
  ExtendedPostHistoryArchivingNode(), P_(__n.P_), S_(__n.S_), B_(__n.B_, *this)
{
  // copy parameter struct P_
  P_.C_m = __n.P_.C_m;
  P_.g_L = __n.P_.g_L;
  P_.E_L = __n.P_.E_L;
  P_.refr_T = __n.P_.refr_T;
  P_.V_th = __n.P_.V_th;
  P_.V_reset = __n.P_.V_reset;
  P_.use_offset_trick = __n.P_.use_offset_trick;
  P_.E_exc = __n.P_.E_exc;
  P_.E_inh = __n.P_.E_inh;
  P_.tau_syn_exc = __n.P_.tau_syn_exc;
  P_.tau_syn_inh = __n.P_.tau_syn_inh;
  P_.I_e = __n.P_.I_e;

  // copy state struct S_
  S_.ode_state[State_::V_m] = __n.S_.ode_state[State_::V_m];
  S_.ode_state[State_::refr_t] = __n.S_.ode_state[State_::refr_t];
  S_.ode_state[State_::g_inh__X__inh_spikes__DOT__weight] = __n.S_.ode_state[State_::g_inh__X__inh_spikes__DOT__weight];
  S_.ode_state[State_::g_inh__X__inh_spikes__DOT__weight__d] = __n.S_.ode_state[State_::g_inh__X__inh_spikes__DOT__weight__d];
  S_.ode_state[State_::g_exc__X__exc_spikes__DOT__weight] = __n.S_.ode_state[State_::g_exc__X__exc_spikes__DOT__weight];
  S_.ode_state[State_::g_exc__X__exc_spikes__DOT__weight__d] = __n.S_.ode_state[State_::g_exc__X__exc_spikes__DOT__weight__d];

  // copy internals V_
  V_.__h = __n.V_.__h;
  V_.__P__refr_t__refr_t = __n.V_.__P__refr_t__refr_t;
  V_.__P__g_inh__X__inh_spikes__DOT__weight__g_inh__X__inh_spikes__DOT__weight = __n.V_.__P__g_inh__X__inh_spikes__DOT__weight__g_inh__X__inh_spikes__DOT__weight;
  V_.__P__g_inh__X__inh_spikes__DOT__weight__g_inh__X__inh_spikes__DOT__weight__d = __n.V_.__P__g_inh__X__inh_spikes__DOT__weight__g_inh__X__inh_spikes__DOT__weight__d;
  V_.__P__g_inh__X__inh_spikes__DOT__weight__d__g_inh__X__inh_spikes__DOT__weight = __n.V_.__P__g_inh__X__inh_spikes__DOT__weight__d__g_inh__X__inh_spikes__DOT__weight;
  V_.__P__g_inh__X__inh_spikes__DOT__weight__d__g_inh__X__inh_spikes__DOT__weight__d = __n.V_.__P__g_inh__X__inh_spikes__DOT__weight__d__g_inh__X__inh_spikes__DOT__weight__d;
  V_.__P__g_exc__X__exc_spikes__DOT__weight__g_exc__X__exc_spikes__DOT__weight = __n.V_.__P__g_exc__X__exc_spikes__DOT__weight__g_exc__X__exc_spikes__DOT__weight;
  V_.__P__g_exc__X__exc_spikes__DOT__weight__g_exc__X__exc_spikes__DOT__weight__d = __n.V_.__P__g_exc__X__exc_spikes__DOT__weight__g_exc__X__exc_spikes__DOT__weight__d;
  V_.__P__g_exc__X__exc_spikes__DOT__weight__d__g_exc__X__exc_spikes__DOT__weight = __n.V_.__P__g_exc__X__exc_spikes__DOT__weight__d__g_exc__X__exc_spikes__DOT__weight;
  V_.__P__g_exc__X__exc_spikes__DOT__weight__d__g_exc__X__exc_spikes__DOT__weight__d = __n.V_.__P__g_exc__X__exc_spikes__DOT__weight__d__g_exc__X__exc_spikes__DOT__weight__d;
}

// ---------------------------------------------------------------------------
//   Destructor for node
// ---------------------------------------------------------------------------

iaf_cond_alpha_neuron_nestml::~iaf_cond_alpha_neuron_nestml()
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
void iaf_cond_alpha_neuron_nestml::calibrate_time( const nest::TimeConverter& tc )
{
  LOG( nest::M_WARNING,
    "iaf_cond_alpha_neuron_nestml",
    "Simulation resolution has changed. Internal state and parameters of the model have been reset!" );

  init_state_internal_();
}
void iaf_cond_alpha_neuron_nestml::init_state_internal_()
{
#ifdef DEBUG
  std::cout << "[neuron " << this << "] iaf_cond_alpha_neuron_nestml::init_state_internal_()" << std::endl;
#endif

  const double __timestep = nest::Time::get_resolution().get_ms();  // do not remove, this is necessary for the timestep() function
  // by default, integrate all variables with a conservative tolerance, in the sense that we err on the side of integrating very precisely at the expense of extra computation
  P_.__gsl_abs_error_tol = 1e-6;
  P_.__gsl_rel_error_tol = 1e-6;
  // initial values for parameters
  P_.C_m = 250; // as pF
  P_.g_L = 16.6667; // as nS
  P_.E_L = (-70); // as mV
  P_.refr_T = 2; // as ms
  P_.V_th = (-55); // as mV
  P_.V_reset = (-60); // as mV
  P_.E_exc = 0; // as mV
  P_.use_offset_trick = 0;
  P_.E_inh = (-85); // as mV
  P_.tau_syn_exc = 0.2; // as ms
  P_.tau_syn_inh = 2; // as ms
  P_.I_e = 0; // as pA

  V_.__h = nest::Time::get_resolution().get_ms();
  recompute_internal_variables();
  // initial values for state variables
  S_.ode_state[State_::V_m] = P_.E_L; // as mV
  S_.ode_state[State_::refr_t] = 0; // as ms
  S_.ode_state[State_::g_inh__X__inh_spikes__DOT__weight] = 0; // as real
  S_.ode_state[State_::g_inh__X__inh_spikes__DOT__weight__d] = 0; // as real
  S_.ode_state[State_::g_exc__X__exc_spikes__DOT__weight] = 0; // as real
  S_.ode_state[State_::g_exc__X__exc_spikes__DOT__weight__d] = 0; // as real
}

void iaf_cond_alpha_neuron_nestml::init_buffers_()
{
#ifdef DEBUG
  std::cout << "[neuron " << this << "] iaf_cond_alpha_neuron_nestml::init_buffers_()" << std::endl;
#endif
  // spike input buffers -- note that .clear() includes a resize
  B_.spike_input_exc_spikes__DOT__weight_.clear();
  B_.spike_input_exc_spikes_spike_input_received_.clear();
  B_.spike_input_inh_spikes__DOT__weight_.clear();
  B_.spike_input_inh_spikes_spike_input_received_.clear();


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
    B_.__c = gsl_odeiv_control_y_new( P_.__gsl_abs_error_tol, P_.__gsl_rel_error_tol );
  }
  else
  {
    gsl_odeiv_control_init( B_.__c, P_.__gsl_abs_error_tol, P_.__gsl_rel_error_tol, 1.0, 0.0 );

  }

  if ( not B_.__e )
  {
    B_.__e = gsl_odeiv_evolve_alloc( State_::STATE_VEC_SIZE );
  }
  else
  {
    gsl_odeiv_evolve_reset( B_.__e );
  }

  // B_.__sys.function = iaf_cond_alpha_neuron_nestml_dynamics; // will be set just prior to the call to gsl_odeiv_evolve_apply()
  B_.__sys.jacobian = nullptr;
  B_.__sys.dimension = State_::STATE_VEC_SIZE;
  B_.__sys.params = reinterpret_cast< void* >( this );
  B_.__step = nest::Time::get_resolution().get_ms();
  B_.__integration_step = nest::Time::get_resolution().get_ms();
}

void iaf_cond_alpha_neuron_nestml::recompute_internal_variables(bool exclude_timestep)
{
  const double __timestep = nest::Time::get_resolution().get_ms();  // do not remove, this is necessary for the timestep() function

  if (exclude_timestep)
  {
    V_.__P__refr_t__refr_t = 1.0; // as real
    V_.__P__g_inh__X__inh_spikes__DOT__weight__g_inh__X__inh_spikes__DOT__weight = 1.0 * (V_.__h + P_.tau_syn_inh) * std::exp((-V_.__h) / P_.tau_syn_inh) / P_.tau_syn_inh; // as real
    V_.__P__g_inh__X__inh_spikes__DOT__weight__g_inh__X__inh_spikes__DOT__weight__d = 1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn_inh); // as real
    V_.__P__g_inh__X__inh_spikes__DOT__weight__d__g_inh__X__inh_spikes__DOT__weight = (-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn_inh) / pow(P_.tau_syn_inh, 2); // as real
    V_.__P__g_inh__X__inh_spikes__DOT__weight__d__g_inh__X__inh_spikes__DOT__weight__d = 1.0 * ((-V_.__h) + P_.tau_syn_inh) * std::exp((-V_.__h) / P_.tau_syn_inh) / P_.tau_syn_inh; // as real
    V_.__P__g_exc__X__exc_spikes__DOT__weight__g_exc__X__exc_spikes__DOT__weight = 1.0 * (V_.__h + P_.tau_syn_exc) * std::exp((-V_.__h) / P_.tau_syn_exc) / P_.tau_syn_exc; // as real
    V_.__P__g_exc__X__exc_spikes__DOT__weight__g_exc__X__exc_spikes__DOT__weight__d = 1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn_exc); // as real
    V_.__P__g_exc__X__exc_spikes__DOT__weight__d__g_exc__X__exc_spikes__DOT__weight = (-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn_exc) / pow(P_.tau_syn_exc, 2); // as real
    V_.__P__g_exc__X__exc_spikes__DOT__weight__d__g_exc__X__exc_spikes__DOT__weight__d = 1.0 * ((-V_.__h) + P_.tau_syn_exc) * std::exp((-V_.__h) / P_.tau_syn_exc) / P_.tau_syn_exc; // as real
  }
  else {
    V_.__h = nest::Time::get_resolution().get_ms(); // as ms
    V_.__P__refr_t__refr_t = 1.0; // as real
    V_.__P__g_inh__X__inh_spikes__DOT__weight__g_inh__X__inh_spikes__DOT__weight = 1.0 * (V_.__h + P_.tau_syn_inh) * std::exp((-V_.__h) / P_.tau_syn_inh) / P_.tau_syn_inh; // as real
    V_.__P__g_inh__X__inh_spikes__DOT__weight__g_inh__X__inh_spikes__DOT__weight__d = 1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn_inh); // as real
    V_.__P__g_inh__X__inh_spikes__DOT__weight__d__g_inh__X__inh_spikes__DOT__weight = (-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn_inh) / pow(P_.tau_syn_inh, 2); // as real
    V_.__P__g_inh__X__inh_spikes__DOT__weight__d__g_inh__X__inh_spikes__DOT__weight__d = 1.0 * ((-V_.__h) + P_.tau_syn_inh) * std::exp((-V_.__h) / P_.tau_syn_inh) / P_.tau_syn_inh; // as real
    V_.__P__g_exc__X__exc_spikes__DOT__weight__g_exc__X__exc_spikes__DOT__weight = 1.0 * (V_.__h + P_.tau_syn_exc) * std::exp((-V_.__h) / P_.tau_syn_exc) / P_.tau_syn_exc; // as real
    V_.__P__g_exc__X__exc_spikes__DOT__weight__g_exc__X__exc_spikes__DOT__weight__d = 1.0 * V_.__h * std::exp((-V_.__h) / P_.tau_syn_exc); // as real
    V_.__P__g_exc__X__exc_spikes__DOT__weight__d__g_exc__X__exc_spikes__DOT__weight = (-1.0) * V_.__h * std::exp((-V_.__h) / P_.tau_syn_exc) / pow(P_.tau_syn_exc, 2); // as real
    V_.__P__g_exc__X__exc_spikes__DOT__weight__d__g_exc__X__exc_spikes__DOT__weight__d = 1.0 * ((-V_.__h) + P_.tau_syn_exc) * std::exp((-V_.__h) / P_.tau_syn_exc) / P_.tau_syn_exc; // as real
  }
}
void iaf_cond_alpha_neuron_nestml::pre_run_hook()
{
#ifdef DEBUG
  std::cout << "[neuron " << this << "] iaf_cond_alpha_neuron_nestml::pre_run_hook()" << std::endl;
#endif

  B_.logger_.init();

  // parameters might have changed -- recompute internals
  V_.__h = nest::Time::get_resolution().get_ms();
  recompute_internal_variables();
}

// ---------------------------------------------------------------------------
//   Update and spike handling functions
// ---------------------------------------------------------------------------

extern "C" inline int iaf_cond_alpha_neuron_nestml_dynamics_refr_t(double __time, const double ode_state[], double f[], void* pnode)
{
  typedef iaf_cond_alpha_neuron_nestml::State_ State_;
  // get access to node so we can almost work as in a member function
  assert( pnode );
  const iaf_cond_alpha_neuron_nestml& node = *( reinterpret_cast< iaf_cond_alpha_neuron_nestml* >( pnode ) );

  // ode_state[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.ode_state[].
  f[State_::V_m] = 0;
  f[State_::refr_t] = (-1.0);
  f[State_::g_inh__X__inh_spikes__DOT__weight] = 1.0 * ode_state[State_::g_inh__X__inh_spikes__DOT__weight__d];
  f[State_::g_inh__X__inh_spikes__DOT__weight__d] = (-ode_state[State_::g_inh__X__inh_spikes__DOT__weight]) / pow(node.P_.tau_syn_inh, 2) - 2 * ode_state[State_::g_inh__X__inh_spikes__DOT__weight__d] / node.P_.tau_syn_inh;
  f[State_::g_exc__X__exc_spikes__DOT__weight] = 1.0 * ode_state[State_::g_exc__X__exc_spikes__DOT__weight__d];
  f[State_::g_exc__X__exc_spikes__DOT__weight__d] = (-ode_state[State_::g_exc__X__exc_spikes__DOT__weight]) / pow(node.P_.tau_syn_exc, 2) - 2 * ode_state[State_::g_exc__X__exc_spikes__DOT__weight__d] / node.P_.tau_syn_exc;
  return GSL_SUCCESS;
}extern "C" inline int iaf_cond_alpha_neuron_nestml_dynamics_V_m(double __time, const double ode_state[], double f[], void* pnode)
{
  typedef iaf_cond_alpha_neuron_nestml::State_ State_;
  // get access to node so we can almost work as in a member function
  assert( pnode );
  const iaf_cond_alpha_neuron_nestml& node = *( reinterpret_cast< iaf_cond_alpha_neuron_nestml* >( pnode ) );

  // ode_state[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.ode_state[].
  f[State_::V_m] = (node.P_.I_e + node.B_.I_stim_grid_sum_ - node.P_.g_L * ((-node.P_.E_L) + ode_state[State_::V_m]) - ode_state[State_::g_exc__X__exc_spikes__DOT__weight] * ((-node.P_.E_exc) + ode_state[State_::V_m]) - ode_state[State_::g_inh__X__inh_spikes__DOT__weight] * ((-node.P_.E_inh) + ode_state[State_::V_m])) / node.P_.C_m;
  f[State_::refr_t] = 0;
  f[State_::g_inh__X__inh_spikes__DOT__weight] = 1.0 * ode_state[State_::g_inh__X__inh_spikes__DOT__weight__d];
  f[State_::g_inh__X__inh_spikes__DOT__weight__d] = (-ode_state[State_::g_inh__X__inh_spikes__DOT__weight]) / pow(node.P_.tau_syn_inh, 2) - 2 * ode_state[State_::g_inh__X__inh_spikes__DOT__weight__d] / node.P_.tau_syn_inh;
  f[State_::g_exc__X__exc_spikes__DOT__weight] = 1.0 * ode_state[State_::g_exc__X__exc_spikes__DOT__weight__d];
  f[State_::g_exc__X__exc_spikes__DOT__weight__d] = (-ode_state[State_::g_exc__X__exc_spikes__DOT__weight]) / pow(node.P_.tau_syn_exc, 2) - 2 * ode_state[State_::g_exc__X__exc_spikes__DOT__weight__d] / node.P_.tau_syn_exc;
  return GSL_SUCCESS;
}

void iaf_cond_alpha_neuron_nestml::update(nest::Time const & origin, const long from, const long to)
{
  const double __timestep = nest::Time::get_resolution().get_ms();  // do not remove, this is necessary for the timestep() function

  for ( long lag = from ; lag < to ; ++lag )
  {
    auto get_t = [origin, lag](){ return nest::Time( nest::Time::step( origin.get_steps() + lag + 1) ).get_ms(); };

    /**
     * buffer spikes from spiking input ports
    **/
    B_.spike_input_exc_spikes__DOT__weight_grid_sum_ = B_.spike_input_exc_spikes__DOT__weight_.get_value(lag);
    B_.spike_input_exc_spikes_spike_input_received_grid_sum_ = B_.spike_input_exc_spikes_spike_input_received_.get_value(lag);
    B_.spike_input_inh_spikes__DOT__weight_grid_sum_ = B_.spike_input_inh_spikes__DOT__weight_.get_value(lag);
    B_.spike_input_inh_spikes_spike_input_received_grid_sum_ = B_.spike_input_inh_spikes_spike_input_received_.get_value(lag);

    /**
     * subthreshold updates of the convolution variables
     *
     * step 1: regardless of whether and how integrate_odes() will be called, update variables due to convolutions
    **/

    const double g_inh__X__inh_spikes__DOT__weight__tmp_ = V_.__P__g_inh__X__inh_spikes__DOT__weight__g_inh__X__inh_spikes__DOT__weight * S_.ode_state[State_::g_inh__X__inh_spikes__DOT__weight] + V_.__P__g_inh__X__inh_spikes__DOT__weight__g_inh__X__inh_spikes__DOT__weight__d * S_.ode_state[State_::g_inh__X__inh_spikes__DOT__weight__d];
    const double g_inh__X__inh_spikes__DOT__weight__d__tmp_ = V_.__P__g_inh__X__inh_spikes__DOT__weight__d__g_inh__X__inh_spikes__DOT__weight * S_.ode_state[State_::g_inh__X__inh_spikes__DOT__weight] + V_.__P__g_inh__X__inh_spikes__DOT__weight__d__g_inh__X__inh_spikes__DOT__weight__d * S_.ode_state[State_::g_inh__X__inh_spikes__DOT__weight__d];
    const double g_exc__X__exc_spikes__DOT__weight__tmp_ = V_.__P__g_exc__X__exc_spikes__DOT__weight__g_exc__X__exc_spikes__DOT__weight * S_.ode_state[State_::g_exc__X__exc_spikes__DOT__weight] + V_.__P__g_exc__X__exc_spikes__DOT__weight__g_exc__X__exc_spikes__DOT__weight__d * S_.ode_state[State_::g_exc__X__exc_spikes__DOT__weight__d];
    const double g_exc__X__exc_spikes__DOT__weight__d__tmp_ = V_.__P__g_exc__X__exc_spikes__DOT__weight__d__g_exc__X__exc_spikes__DOT__weight * S_.ode_state[State_::g_exc__X__exc_spikes__DOT__weight] + V_.__P__g_exc__X__exc_spikes__DOT__weight__d__g_exc__X__exc_spikes__DOT__weight__d * S_.ode_state[State_::g_exc__X__exc_spikes__DOT__weight__d];


    /**
     * Begin NESTML generated code for the update block(s)
    **/

  if (S_.ode_state[State_::refr_t] > 0)
  {

    // start rendered code for integrate_odes(refr_t)

    // analytic solver: integrating state variables (first step): refr_t,
    const double refr_t__tmp = V_.__P__refr_t__refr_t * S_.ode_state[State_::refr_t] - 1.0 * V_.__h;


    // numeric solver: integrating state variables: refr_t,
    double __t = 0;
    B_.__sys.function = iaf_cond_alpha_neuron_nestml_dynamics_refr_t;
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
    // analytic solver: integrating state variables (second step): refr_t,
    /* replace analytically solvable variables with precisely integrated values  */
    S_.ode_state[State_::refr_t] = refr_t__tmp;
  }
  else
  {

    // start rendered code for integrate_odes(V_m)



    // numeric solver: integrating state variables: V_m,
    double __t = 0;
    B_.__sys.function = iaf_cond_alpha_neuron_nestml_dynamics_V_m;
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
  }

    /**
     * Begin NESTML generated code for the onReceive block(s)
    **/


    /**
     * subthreshold updates of the convolution variables
     *
     * step 2: regardless of whether and how integrate_odes() was called, update variables due to convolutions. Set to the updated values at the end of the timestep.
    **/

    S_.ode_state[State_::g_inh__X__inh_spikes__DOT__weight] = g_inh__X__inh_spikes__DOT__weight__tmp_;
    S_.ode_state[State_::g_inh__X__inh_spikes__DOT__weight__d] = g_inh__X__inh_spikes__DOT__weight__d__tmp_;
    S_.ode_state[State_::g_exc__X__exc_spikes__DOT__weight] = g_exc__X__exc_spikes__DOT__weight__tmp_;
    S_.ode_state[State_::g_exc__X__exc_spikes__DOT__weight__d] = g_exc__X__exc_spikes__DOT__weight__d__tmp_;

    /**
     * spike updates due to convolutions
    **/

    S_.ode_state[State_::g_inh__X__inh_spikes__DOT__weight__d] += (B_.spike_input_inh_spikes__DOT__weight_grid_sum_) * (numerics::e / P_.tau_syn_inh) / (1.0);
    S_.ode_state[State_::g_exc__X__exc_spikes__DOT__weight__d] += (B_.spike_input_exc_spikes__DOT__weight_grid_sum_) * (numerics::e / P_.tau_syn_exc) / (1.0);

    /**
     * Begin NESTML generated code for the onCondition block(s)
    **/

    if (S_.ode_state[State_::refr_t] <= 0 && S_.ode_state[State_::V_m] >= P_.V_th)
    {
      S_.ode_state[State_::refr_t] = P_.refr_T;
      S_.ode_state[State_::V_m] = P_.V_reset;

      // begin generated code for emit_spike() function

      #ifdef DEBUG
      std::cout << "Emitting a spike at t = " << nest::Time(nest::Time::step(origin.get_steps() + lag + 1)).get_ms() << "\n";
      #endif
      if (get_use_offset_trick() and get_t() < 10.)
      {
        set_spiketime( nest::Time::step( origin.get_steps() + lag + 1 ), 1. );
      }
      else
      {
        set_spiketime( nest::Time::step( origin.get_steps() + lag + 1 ), 0. );
      }

      nest::SpikeEvent se;
      nest::kernel().event_delivery_manager.send(*this, se, lag);
      // end generated code for emit_spike() function
    }

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
void iaf_cond_alpha_neuron_nestml::handle(nest::DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}


void iaf_cond_alpha_neuron_nestml::handle(nest::SpikeEvent &e)
{
#ifdef DEBUG
  std::cout << "[neuron " << this << "] iaf_cond_alpha_neuron_nestml::handle(SpikeEvent) on rport " << e.get_rport() << std::endl;
#endif

  assert(e.get_delay_steps() > 0);
  assert(e.get_rport() < 2);
    B_.spike_input_exc_spikes__DOT__weight_.add_value(
      e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin() ),
      e.get_weight() * e.get_multiplicity() );

    B_.spike_input_exc_spikes_spike_input_received_.add_value(
      e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin() ),
      1. );
}

void iaf_cond_alpha_neuron_nestml::handle(nest::CurrentEvent& e)
{
#ifdef DEBUG
  std::cout << "[neuron " << this << "] iaf_cond_alpha_neuron_nestml::handle(CurrentEvent)" << std::endl;
#endif
  assert(e.get_delay_steps() > 0);

  const double current = e.get_current();     // we assume that in NEST, this returns a current in pA
  const double weight = e.get_weight();
  get_I_stim().add_value(
               e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
               weight * current );
}

// -------------------------------------------------------------------------
//   Methods corresponding to onReceive blocks
// -------------------------------------------------------------------------

