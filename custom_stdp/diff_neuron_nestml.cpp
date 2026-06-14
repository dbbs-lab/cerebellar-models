// #define DEBUG 1
/*
 *  diff_neuron_nestml.cpp
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
 *  Generated from NESTML at time: 2024-10-11 08:49:24.053495
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

#include "diff_neuron_nestml.h"
void
register_diff_neuron_nestml( const std::string& name )
{
  nest::register_node_model< diff_neuron_nestml >( name );
}

// ---------------------------------------------------------------------------
//   Recordables map
// ---------------------------------------------------------------------------
namespace nest
{

  // Override the create() method with one call to RecordablesMap::insert_()
  // for each quantity to be recorded.
template <> void DynamicRecordablesMap<diff_neuron_nestml>::create(diff_neuron_nestml& host)
  {
    insert("in_rate", host.get_data_access_functor( diff_neuron_nestml::State_::IN_RATE ));
    insert("in_rate_pre", host.get_data_access_functor( diff_neuron_nestml::State_::IN_RATE_PRE ));
    insert("out_rate", host.get_data_access_functor( diff_neuron_nestml::State_::OUT_RATE ));
    insert("spike_count_in", host.get_data_access_functor( diff_neuron_nestml::State_::SPIKE_COUNT_IN ));
    insert("spike_count_in_pos", host.get_data_access_functor( diff_neuron_nestml::State_::SPIKE_COUNT_IN_POS ));
    insert("spike_count_in_neg", host.get_data_access_functor( diff_neuron_nestml::State_::SPIKE_COUNT_IN_NEG ));
    insert("lambda_poisson", host.get_data_access_functor( diff_neuron_nestml::State_::LAMBDA_POISSON ));

    // Add vector variables  
      host.insert_recordables();
  }
}
  std::string diff_neuron_nestml::get_var_name(size_t elem, std::string var_name)
  {
    std::stringstream n;
    n << var_name << elem;
    return n.str();
  }

  void diff_neuron_nestml::insert_recordables(size_t first)
  {
      for (size_t i = 0; i < 
P_.simulation_steps; i++)
      {
        size_t elem = diff_neuron_nestml::State_::SPIKES_BUFFER + i;
        recordablesMap_.insert(get_var_name(i, "SPIKES_BUFFER_"), this->get_data_access_functor(elem));
      }
  }

  nest::DataAccessFunctor< diff_neuron_nestml >
  diff_neuron_nestml::get_data_access_functor( size_t elem )
  {
    return nest::DataAccessFunctor< diff_neuron_nestml >( *this, elem );
  }

// ---------------------------------------------------------------------------
//   Default constructors defining default parameters and state
//   Note: the implementation is empty. The initialization is of variables
//   is a part of diff_neuron_nestml's constructor.
// ---------------------------------------------------------------------------

diff_neuron_nestml::Parameters_::Parameters_()
{
}

diff_neuron_nestml::State_::State_()
{
}

// ---------------------------------------------------------------------------
//   Parameter and state extractions and manipulation functions
// ---------------------------------------------------------------------------

diff_neuron_nestml::Buffers_::Buffers_(diff_neuron_nestml &n):
  logger_(n)
  , spike_inputs_( std::vector< nest::RingBuffer >( NUM_SPIKE_RECEPTORS ) )
  , spike_inputs_grid_sum_( std::vector< double >( NUM_SPIKE_RECEPTORS ) )
  , spike_input_received_( std::vector< nest::RingBuffer >( NUM_SPIKE_RECEPTORS ) )
  , spike_input_received_grid_sum_( std::vector< double >( NUM_SPIKE_RECEPTORS ) )
{
  // Initialization of the remaining members is deferred to init_buffers_().
}

diff_neuron_nestml::Buffers_::Buffers_(const Buffers_ &, diff_neuron_nestml &n):
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

diff_neuron_nestml::diff_neuron_nestml():StructuralPlasticityNode(), P_(), S_(), B_(*this)
{
  init_state_internal_();
  recordablesMap_.create(*this);
  pre_run_hook();
}

// ---------------------------------------------------------------------------
//   Copy constructor for node
// ---------------------------------------------------------------------------

diff_neuron_nestml::diff_neuron_nestml(const diff_neuron_nestml& __n):
  StructuralPlasticityNode(), P_(__n.P_), S_(__n.S_), B_(__n.B_, *this)
{

  // copy parameter struct P_
  P_.kp = __n.P_.kp;
  P_.pos = __n.P_.pos;
  P_.base_rate = __n.P_.base_rate;
  P_.buffer_size = __n.P_.buffer_size;
  P_.simulation_steps = __n.P_.simulation_steps;
  P_.time_wait = __n.P_.time_wait;

  // copy state struct S_
  S_.in_rate = __n.S_.in_rate;
  S_.in_rate_pre = __n.S_.in_rate_pre;
  S_.out_rate = __n.S_.out_rate;
  S_.spike_count_in = __n.S_.spike_count_in;
  S_.spike_count_in_pos = __n.S_.spike_count_in_pos;
  S_.spike_count_in_neg = __n.S_.spike_count_in_neg;
  S_.spike_count_out = __n.S_.spike_count_out;
  S_.tick = __n.S_.tick;
  S_.lambda_poisson = __n.S_.lambda_poisson;
  S_.spikes_buffer = __n.S_.spikes_buffer;

  // copy internals V_
  V_.res = __n.V_.res;
  V_.window_counts = __n.V_.window_counts;
  V_.__h = __n.V_.__h;
  recordablesMap_.create(*this);
}

// ---------------------------------------------------------------------------
//   Destructor for node
// ---------------------------------------------------------------------------

diff_neuron_nestml::~diff_neuron_nestml()
{
}

// ---------------------------------------------------------------------------
//   Node initialization functions
// ---------------------------------------------------------------------------
void diff_neuron_nestml::calibrate_time( const nest::TimeConverter& tc )
{
  LOG( nest::M_WARNING,
    "diff_neuron_nestml",
    "Simulation resolution has changed. Internal state and parameters of the model have been reset!" );

  init_state_internal_();
}
void diff_neuron_nestml::init_state_internal_()
{
#ifdef DEBUG
  std::cout << "diff_neuron_nestml::init_state_internal_()" << std::endl;
#endif

  const double __resolution = nest::Time::get_resolution().get_ms();  // do not remove, this is necessary for the resolution() function
  // initial values for parameters
  P_.kp = 1.0; // as real
  P_.pos = true; // as boolean
  P_.base_rate = 0; // as Hz
  P_.buffer_size = 100; // as ms
  P_.simulation_steps = 1000; // as integer
  P_.time_wait = 300.0; // as real

  V_.__h = nest::Time::get_resolution().get_ms();
  recompute_internal_variables();
  // initial values for state variables
  S_.in_rate = 0; // as Hz
  S_.in_rate_pre = 0.0; // as real
  S_.out_rate = 0; // as Hz
  S_.spike_count_in = 0.0; // as real
  S_.spike_count_in_pos = 0.0; // as real
  S_.spike_count_in_neg = 0.0; // as real
  S_.spike_count_out = 0; // as integer
  S_.tick = 0; // as integer
  S_.lambda_poisson = 0; // as real
  S_.spikes_buffer.resize(
  P_.simulation_steps, 0);
}

void diff_neuron_nestml::init_buffers_()
{
#ifdef DEBUG
  std::cout << "diff_neuron_nestml::init_buffers_()" << std::endl;
#endif
  // spike input buffers
  get_spike_inputs_().clear();
  get_spike_inputs_grid_sum_().clear();
  get_spike_input_received_().clear();
  get_spike_input_received_grid_sum_().clear();


  B_.logger_.reset();


}

void diff_neuron_nestml::recompute_internal_variables(bool exclude_timestep)
{
  const double __resolution = nest::Time::get_resolution().get_ms();  // do not remove, this is necessary for the resolution() function

  if (exclude_timestep)
  {    
    V_.res = __resolution; // as ms
    V_.window_counts = nest::Time(nest::Time::ms((double) (P_.buffer_size))).get_steps(); // as integer
  }
  else {    
    V_.res = __resolution; // as ms
    V_.__h = __resolution; // as ms
    V_.window_counts = nest::Time(nest::Time::ms((double) (P_.buffer_size))).get_steps(); // as integer
  }
}
void diff_neuron_nestml::pre_run_hook()
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


void diff_neuron_nestml::update(nest::Time const & origin,const long from, const long to)
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
    S_.spikes_buffer[S_.tick] = (0.001 * B_.spike_inputs_grid_sum_[SPIKES - MIN_SPIKE_RECEPTOR]);
    long i = 0;
    long index = 0;
    S_.spike_count_in = 0;
    S_.spike_count_in_pos = 0;
    S_.spike_count_in_neg = 0;
    for ( i = 0;
                     i<V_.window_counts;
         i += 1 )
    {
      index = S_.tick - i;
      if ((index >= 0 && S_.spikes_buffer[index] != 0))
      {  
          S_.spike_count_in += S_.spikes_buffer[index];
          if (S_.spikes_buffer[index] > 0)
          {  
              S_.spike_count_in_pos += S_.spikes_buffer[index];
          }
          else
          {  
              S_.spike_count_in_neg += S_.spikes_buffer[index];
          }
      }
    }
    S_.in_rate_pre = (1000.0 * std::abs(S_.spike_count_in)) / P_.buffer_size;
    long lambda_exp = 0;
    lambda_exp = std::max(std::abs(S_.spike_count_in_pos), std::abs(S_.spike_count_in_neg));
    double thresh = pow((2 * lambda_exp), 0.5);
    if (std::abs(S_.spike_count_in) < thresh)
    {  
        S_.spike_count_in = 0;
    }
    else if (S_.spike_count_in > 0)
    {  
        S_.spike_count_in = S_.spike_count_in - thresh;
    }
    else if (S_.spike_count_in < 0)
    {  
        S_.spike_count_in = S_.spike_count_in + thresh;
    }
    if (((S_.spike_count_in < 0 && P_.pos == true) || (S_.spike_count_in > 0 && P_.pos == false)))
    {  
        S_.spike_count_in = 0;
    }
    S_.in_rate = (1000.0 * ((1000.0 * std::abs(S_.spike_count_in)) / P_.buffer_size));
    S_.out_rate = P_.base_rate + P_.kp * S_.in_rate;
    S_.lambda_poisson = S_.out_rate * __resolution * 0.001;
    S_.spike_count_out = ([&]() -> int { nest::poisson_distribution::param_type poisson_params(S_.lambda_poisson); int sample = poisson_dev_( nest::get_vp_specific_rng( get_thread() ), poisson_params); return sample; })();
    if (S_.spike_count_out > 0)
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
void diff_neuron_nestml::handle(nest::DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}


void diff_neuron_nestml::handle(nest::SpikeEvent &e)
{
  assert(e.get_delay_steps() > 0);
  assert( e.get_rport() < B_.spike_inputs_.size() );

  double weight = e.get_weight();
  size_t nestml_buffer_idx = 0;
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

