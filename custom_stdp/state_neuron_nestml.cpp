// #define DEBUG 1
/*
 *  state_neuron_nestml.cpp
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
 *  Generated from NESTML at time: 2024-10-11 08:49:24.519807
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

#include "state_neuron_nestml.h"
void
register_state_neuron_nestml( const std::string& name )
{
  nest::register_node_model< state_neuron_nestml >( name );
}

// ---------------------------------------------------------------------------
//   Recordables map
// ---------------------------------------------------------------------------
namespace nest
{

  // Override the create() method with one call to RecordablesMap::insert_()
  // for each quantity to be recorded.
template <> void DynamicRecordablesMap<state_neuron_nestml>::create(state_neuron_nestml& host)
  {
    insert("in_rate", host.get_data_access_functor( state_neuron_nestml::State_::IN_RATE ));
    insert("out_rate", host.get_data_access_functor( state_neuron_nestml::State_::OUT_RATE ));
    insert("mean_fbk", host.get_data_access_functor( state_neuron_nestml::State_::MEAN_FBK ));
    insert("mean_pred", host.get_data_access_functor( state_neuron_nestml::State_::MEAN_PRED ));
    insert("var_fbk", host.get_data_access_functor( state_neuron_nestml::State_::VAR_FBK ));
    insert("var_pred", host.get_data_access_functor( state_neuron_nestml::State_::VAR_PRED ));
    insert("CV_fbk", host.get_data_access_functor( state_neuron_nestml::State_::CV_FBK ));
    insert("CV_pred", host.get_data_access_functor( state_neuron_nestml::State_::CV_PRED ));
    insert("total_CV", host.get_data_access_functor( state_neuron_nestml::State_::TOTAL_CV ));
    insert("lambda_poisson", host.get_data_access_functor( state_neuron_nestml::State_::LAMBDA_POISSON ));

    // Add vector variables  
      host.insert_recordables();
  }
}
std::vector< std::tuple< int, int > > state_neuron_nestml::rport_to_nestml_buffer_idx =
{
  { state_neuron_nestml::FBK_SPIKES_0, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_1, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_2, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_3, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_4, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_5, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_6, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_7, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_8, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_9, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_10, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_11, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_12, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_13, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_14, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_15, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_16, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_17, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_18, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_19, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_20, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_21, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_22, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_23, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_24, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_25, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_26, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_27, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_28, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_29, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_30, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_31, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_32, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_33, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_34, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_35, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_36, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_37, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_38, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_39, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_40, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_41, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_42, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_43, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_44, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_45, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_46, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_47, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_48, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::FBK_SPIKES_49, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_0, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_1, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_2, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_3, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_4, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_5, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_6, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_7, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_8, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_9, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_10, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_11, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_12, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_13, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_14, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_15, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_16, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_17, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_18, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_19, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_20, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_21, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_22, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_23, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_24, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_25, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_26, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_27, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_28, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_29, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_30, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_31, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_32, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_33, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_34, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_35, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_36, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_37, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_38, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_39, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_40, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_41, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_42, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_43, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_44, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_45, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_46, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_47, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_48, state_neuron_nestml::PORT_NOT_AVAILABLE },
  { state_neuron_nestml::PRED_SPIKES_49, state_neuron_nestml::PORT_NOT_AVAILABLE },
};
  std::string state_neuron_nestml::get_var_name(size_t elem, std::string var_name)
  {
    std::stringstream n;
    n << var_name << elem;
    return n.str();
  }

  void state_neuron_nestml::insert_recordables(size_t first)
  {
      for (size_t i = 0; i < 
P_.N_fbk; i++)
      {
        size_t elem = state_neuron_nestml::State_::CURRENT_FBK_INPUT + i;
        recordablesMap_.insert(get_var_name(i, "CURRENT_FBK_INPUT_"), this->get_data_access_functor(elem));
      }
      for (size_t i = 0; i < 
P_.N_pred; i++)
      {
        size_t elem = state_neuron_nestml::State_::CURRENT_PRED_INPUT + i;
        recordablesMap_.insert(get_var_name(i, "CURRENT_PRED_INPUT_"), this->get_data_access_functor(elem));
      }
      for (size_t i = 0; i < 
P_.fbk_bf_size; i++)
      {
        size_t elem = state_neuron_nestml::State_::FBK_BUFFER + i;
        recordablesMap_.insert(get_var_name(i, "FBK_BUFFER_"), this->get_data_access_functor(elem));
      }
      for (size_t i = 0; i < 
P_.pred_bf_size; i++)
      {
        size_t elem = state_neuron_nestml::State_::PRED_BUFFER + i;
        recordablesMap_.insert(get_var_name(i, "PRED_BUFFER_"), this->get_data_access_functor(elem));
      }
      for (size_t i = 0; i < 
P_.N_fbk; i++)
      {
        size_t elem = state_neuron_nestml::State_::FBK_COUNTS + i;
        recordablesMap_.insert(get_var_name(i, "FBK_COUNTS_"), this->get_data_access_functor(elem));
      }
      for (size_t i = 0; i < 
P_.N_pred; i++)
      {
        size_t elem = state_neuron_nestml::State_::PRED_COUNTS + i;
        recordablesMap_.insert(get_var_name(i, "PRED_COUNTS_"), this->get_data_access_functor(elem));
      }
  }

  nest::DataAccessFunctor< state_neuron_nestml >
  state_neuron_nestml::get_data_access_functor( size_t elem )
  {
    return nest::DataAccessFunctor< state_neuron_nestml >( *this, elem );
  }

// ---------------------------------------------------------------------------
//   Default constructors defining default parameters and state
//   Note: the implementation is empty. The initialization is of variables
//   is a part of state_neuron_nestml's constructor.
// ---------------------------------------------------------------------------

state_neuron_nestml::Parameters_::Parameters_()
{
}

state_neuron_nestml::State_::State_()
{
}

// ---------------------------------------------------------------------------
//   Parameter and state extractions and manipulation functions
// ---------------------------------------------------------------------------

state_neuron_nestml::Buffers_::Buffers_(state_neuron_nestml &n):
  logger_(n)
  , spike_inputs_( std::vector< nest::RingBuffer >( NUM_SPIKE_RECEPTORS ) )
  , spike_inputs_grid_sum_( std::vector< double >( NUM_SPIKE_RECEPTORS ) )
  , spike_input_received_( std::vector< nest::RingBuffer >( NUM_SPIKE_RECEPTORS ) )
  , spike_input_received_grid_sum_( std::vector< double >( NUM_SPIKE_RECEPTORS ) )
{
  // Initialization of the remaining members is deferred to init_buffers_().
}

state_neuron_nestml::Buffers_::Buffers_(const Buffers_ &, state_neuron_nestml &n):
  logger_(n)
  , spike_inputs_( std::vector< nest::RingBuffer >( NUM_SPIKE_RECEPTORS ) )
  , spike_inputs_grid_sum_( std::vector< double >( NUM_SPIKE_RECEPTORS ) )
  , spike_input_received_( std::vector< nest::RingBuffer >( NUM_SPIKE_RECEPTORS ) )
  , spike_input_received_grid_sum_( std::vector< double >( NUM_SPIKE_RECEPTORS ) )
{
  // Initialization of the remaining members is deferred to init_buffers_().
}

// ---------------------------------------------------------------------------
//   Default constructor for node
// ---------------------------------------------------------------------------

state_neuron_nestml::state_neuron_nestml():StructuralPlasticityNode(), P_(), S_(), B_(*this)
{
  init_state_internal_();
  recordablesMap_.create(*this);
  pre_run_hook();
}

// ---------------------------------------------------------------------------
//   Copy constructor for node
// ---------------------------------------------------------------------------

state_neuron_nestml::state_neuron_nestml(const state_neuron_nestml& __n):
  StructuralPlasticityNode(), P_(__n.P_), S_(__n.S_), B_(__n.B_, *this)
{

  // copy parameter struct P_
  P_.kp = __n.P_.kp;
  P_.pos = __n.P_.pos;
  P_.base_rate = __n.P_.base_rate;
  P_.buffer_size = __n.P_.buffer_size;
  P_.simulation_steps = __n.P_.simulation_steps;
  P_.N_fbk = __n.P_.N_fbk;
  P_.N_pred = __n.P_.N_pred;
  P_.fbk_bf_size = __n.P_.fbk_bf_size;
  P_.pred_bf_size = __n.P_.pred_bf_size;
  P_.time_wait = __n.P_.time_wait;
  P_.time_trial = __n.P_.time_trial;

  // copy state struct S_
  S_.in_rate = __n.S_.in_rate;
  S_.out_rate = __n.S_.out_rate;
  S_.spike_count_out = __n.S_.spike_count_out;
  S_.current_fbk_input = __n.S_.current_fbk_input;
  S_.current_pred_input = __n.S_.current_pred_input;
  S_.fbk_buffer = __n.S_.fbk_buffer;
  S_.pred_buffer = __n.S_.pred_buffer;
  S_.fbk_counts = __n.S_.fbk_counts;
  S_.pred_counts = __n.S_.pred_counts;
  S_.tick = __n.S_.tick;
  S_.position_count = __n.S_.position_count;
  S_.mean_fbk = __n.S_.mean_fbk;
  S_.mean_pred = __n.S_.mean_pred;
  S_.var_fbk = __n.S_.var_fbk;
  S_.var_pred = __n.S_.var_pred;
  S_.CV_fbk = __n.S_.CV_fbk;
  S_.CV_pred = __n.S_.CV_pred;
  S_.total_CV = __n.S_.total_CV;
  S_.lambda_poisson = __n.S_.lambda_poisson;

  // copy internals V_
  V_.res = __n.V_.res;
  V_.buffer_steps = __n.V_.buffer_steps;
  V_.trial_steps = __n.V_.trial_steps;
  V_.wait_steps = __n.V_.wait_steps;
  V_.__h = __n.V_.__h;
  recordablesMap_.create(*this);
}

// ---------------------------------------------------------------------------
//   Destructor for node
// ---------------------------------------------------------------------------

state_neuron_nestml::~state_neuron_nestml()
{
}

// ---------------------------------------------------------------------------
//   Node initialization functions
// ---------------------------------------------------------------------------
void state_neuron_nestml::calibrate_time( const nest::TimeConverter& tc )
{
  LOG( nest::M_WARNING,
    "state_neuron_nestml",
    "Simulation resolution has changed. Internal state and parameters of the model have been reset!" );

  init_state_internal_();
}
void state_neuron_nestml::init_state_internal_()
{
#ifdef DEBUG
  std::cout << "state_neuron_nestml::init_state_internal_()" << std::endl;
#endif

  const double __resolution = nest::Time::get_resolution().get_ms();  // do not remove, this is necessary for the resolution() function
  // initial values for parameters
  P_.kp = 1; // as real
  P_.pos = true; // as boolean
  P_.base_rate = 0; // as Hz
  P_.buffer_size = 100; // as ms
  P_.simulation_steps = 1000; // as integer
  P_.N_fbk = 50; // as integer
  P_.N_pred = 50; // as integer
  P_.fbk_bf_size = 10000; // as integer
  P_.pred_bf_size = 10000; // as integer
  P_.time_wait = 150.0; // as ms
  P_.time_trial = 650.0; // as ms

  V_.__h = nest::Time::get_resolution().get_ms();
  recompute_internal_variables();
  // initial values for state variables
  S_.in_rate = 0; // as Hz
  S_.out_rate = 0; // as Hz
  S_.spike_count_out = 0; // as integer
  S_.current_fbk_input.resize(
  P_.N_fbk, 0);
  S_.current_pred_input.resize(
  P_.N_pred, 0);
  S_.fbk_buffer.resize(
  P_.fbk_bf_size, 0);
  S_.pred_buffer.resize(
  P_.pred_bf_size, 0);
  S_.fbk_counts.resize(
  P_.N_fbk, 0);
  S_.pred_counts.resize(
  P_.N_pred, 0);
  S_.tick = 0; // as integer
  S_.position_count = 0; // as integer
  S_.mean_fbk = 0.0; // as real
  S_.mean_pred = 0.0; // as real
  S_.var_fbk = 0.0; // as real
  S_.var_pred = 0.0; // as real
  S_.CV_fbk = 0.0; // as real
  S_.CV_pred = 0.0; // as real
  S_.total_CV = 0.0; // as real
  S_.lambda_poisson = 0; // as real
}

void state_neuron_nestml::init_buffers_()
{
#ifdef DEBUG
  std::cout << "state_neuron_nestml::init_buffers_()" << std::endl;
#endif
  // spike input buffers
  get_spike_inputs_().clear();
  get_spike_inputs_grid_sum_().clear();
  get_spike_input_received_().clear();
  get_spike_input_received_grid_sum_().clear();


  B_.logger_.reset();


}

void state_neuron_nestml::recompute_internal_variables(bool exclude_timestep)
{
  const double __resolution = nest::Time::get_resolution().get_ms();  // do not remove, this is necessary for the resolution() function

  if (exclude_timestep)
  {    
    V_.res = __resolution; // as ms
    V_.buffer_steps = nest::Time(nest::Time::ms((double) (P_.buffer_size))).get_steps(); // as integer
    V_.trial_steps = nest::Time(nest::Time::ms((double) (P_.time_trial))).get_steps(); // as integer
    V_.wait_steps = nest::Time(nest::Time::ms((double) (P_.time_wait))).get_steps(); // as integer
  }
  else {    
    V_.res = __resolution; // as ms
    V_.__h = __resolution; // as ms
    V_.buffer_steps = nest::Time(nest::Time::ms((double) (P_.buffer_size))).get_steps(); // as integer
    V_.trial_steps = nest::Time(nest::Time::ms((double) (P_.time_trial))).get_steps(); // as integer
    V_.wait_steps = nest::Time(nest::Time::ms((double) (P_.time_wait))).get_steps(); // as integer
  }
}
void state_neuron_nestml::pre_run_hook()
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


void state_neuron_nestml::update(nest::Time const & origin,const long from, const long to)
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


    /**
     * Begin NESTML generated code for the update block(s)
    **/

    S_.tick = nest::Time(nest::Time::ms((double) (get_t()))).get_steps();
    long i = 0;
    for ( i = 0;
                     i<(P_.N_fbk - 1);
         i += 1 )
    {
      S_.current_fbk_input[i] = (0.001 * B_.spike_inputs_grid_sum_[FBK_SPIKES_0 + i - MIN_SPIKE_RECEPTOR]);
    }
    long j = 0;
    for ( j = 0;
                     j<(P_.N_pred - 1);
         j += 1 )
    {
      S_.current_pred_input[j] = (0.001 * B_.spike_inputs_grid_sum_[PRED_SPIKES_0 + j - MIN_SPIKE_RECEPTOR]);
    }
    long index = 0;
    for ( i = 0;
                     i<(P_.N_fbk - 1);
         i += 1 )
    {
      index = S_.position_count * P_.N_fbk + i;
      S_.fbk_buffer[index] = S_.current_fbk_input[i];
    }
    for ( j = 0;
                     j<(P_.N_pred - 1);
         j += 1 )
    {
      index = S_.position_count * P_.N_pred + j;
      S_.pred_buffer[index] = S_.current_pred_input[j];
    }
    S_.position_count += 1;
    if (S_.position_count > V_.buffer_steps - 1)
    {  
        S_.position_count = 0;
    }
    long k = 0;
    long jump = 0;
    for ( k = 0;
                     k<(P_.N_fbk - 1);
         k += 1 )
    {
      S_.fbk_counts[k] = 0;
      for ( jump = 0;
                       jump<(V_.buffer_steps - 1);
           jump += 1 )
      {
        index = P_.N_fbk * jump + k;
        if (S_.fbk_buffer[index] != 0)
        {  
            S_.fbk_counts[k] += 1;
        }
      }
    }
    long m = 0;
    for ( m = 0;
                     m<(P_.N_pred - 1);
         m += 1 )
    {
      S_.pred_counts[m] = 0;
      for ( jump = 0;
                       jump<(V_.buffer_steps - 1);
           jump += 1 )
      {
        index = (P_.N_pred * jump) + m;
        if (S_.pred_buffer[index] != 0)
        {  
            S_.pred_counts[m] += 1;
        }
      }
    }
    S_.mean_fbk = 0.0;
    if (P_.N_fbk == 0)
    {  
        S_.CV_fbk = pow(10, 6);
    }
    else
    {  
        for ( k = 0;
                         k<(P_.N_fbk - 1);
             k += 1 )
        {
          S_.mean_fbk += S_.fbk_counts[k];
        }
        S_.mean_fbk /= P_.N_fbk;
        if (S_.mean_fbk != 0)
        {  
            S_.var_fbk = 0.0;
            for ( k = 0;
                             k<(P_.N_fbk - 1);
                 k += 1 )
            {
              S_.var_fbk += pow((S_.fbk_counts[k] - S_.mean_fbk), 2);
            }
            S_.var_fbk /= P_.N_fbk;
            S_.CV_fbk = (S_.var_fbk / S_.mean_fbk);
        }
        else
        {  
            S_.CV_fbk = 3.0;
        }
    }
    S_.mean_pred = 0.0;
    if (P_.N_pred == 0)
    {  
        S_.CV_pred = pow(10, 6);
    }
    else
    {  
        for ( m = 0;
                         m<(P_.N_pred - 1);
             m += 1 )
        {
          S_.mean_pred += S_.pred_counts[m];
        }
        S_.mean_pred /= P_.N_pred;
        if (S_.mean_pred != 0)
        {  
            S_.var_pred = 0.0;
            for ( m = 0;
                             m<P_.N_pred;
                 m += 1 )
            {
              S_.var_pred += pow((S_.pred_counts[m] - S_.mean_pred), 2);
            }
            S_.var_pred /= P_.N_pred;
            S_.CV_pred = (S_.var_pred / S_.mean_pred);
        }
        else
        {  
            S_.CV_pred = 3.0;
        }
    }
    S_.total_CV = S_.CV_fbk + S_.CV_pred;
    S_.in_rate = (1000.0 * ((S_.mean_pred * S_.CV_fbk / S_.total_CV + S_.mean_fbk * S_.CV_pred / S_.total_CV) / P_.buffer_size));
    S_.out_rate = P_.base_rate + P_.kp * S_.in_rate;
    S_.lambda_poisson = S_.out_rate * __resolution * 0.001;
    S_.spike_count_out = ([&]() -> int { nest::poisson_distribution::param_type poisson_params(S_.lambda_poisson); int sample = poisson_dev_( nest::get_vp_specific_rng( get_thread() ), poisson_params); return sample; })();
    S_.spike_count_out = ([&]() -> int { nest::poisson_distribution::param_type poisson_params(S_.lambda_poisson); int sample = poisson_dev_( nest::get_vp_specific_rng( get_thread() ), poisson_params); return sample; })();
    if (S_.spike_count_out > 0 && (S_.tick % V_.trial_steps) > V_.wait_steps)
    {  

        /**
         * generated code for emit_spike() function
        **/

        set_spiketime(nest::Time::step(origin.get_steps() + lag + 1));
        nest::SpikeEvent se;
        nest::kernel().event_delivery_manager.send(*this, se, lag);


    }

    /**
     * Begin NESTML generated code for the onReceive block(s)
    **/


    /**
     * subthreshold updates of the convolution variables
     *
     * step 2: regardless of whether and how integrate_odes() was called, update variables due to convolutions. Set to the updated values at the end of the timestep.
    **/



    /**
     * spike updates due to convolutions
    **/


    /**
     * Begin NESTML generated code for the onCondition block(s)
    **/


    /**
     * handle continuous input ports
    **/
    // voltage logging
    B_.logger_.record_data(origin.get_steps() + lag);
  }
}

// Do not move this function as inline to h-file. It depends on
// universal_data_logger_impl.h being included here.
void state_neuron_nestml::handle(nest::DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}


void state_neuron_nestml::handle(nest::SpikeEvent &e)
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
    if ( nestml_buffer_idx == state_neuron_nestml::PORT_NOT_AVAILABLE )
    {
      nestml_buffer_idx = std::get<0>(rport_to_nestml_buffer_idx[e.get_rport()]);
    }
    weight = -weight;
  }
  B_.spike_inputs_[ nestml_buffer_idx - MIN_SPIKE_RECEPTOR ].add_value(
    e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin() ),
    weight * e.get_multiplicity() );
  B_.spike_input_received_[ nestml_buffer_idx - MIN_SPIKE_RECEPTOR ].add_value(
    e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin() ),
    1. );
}

// -------------------------------------------------------------------------
//   Methods corresponding to event handlers
// -------------------------------------------------------------------------

