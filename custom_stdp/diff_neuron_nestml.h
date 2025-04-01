
/**
 *  diff_neuron_nestml.h
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
#ifndef DIFF_NEURON_NESTML
#define DIFF_NEURON_NESTML

#ifndef HAVE_LIBLTDL
#error "NEST was compiled without support for dynamic loading. Please install libltdl and recompile NEST."
#endif

// C++ includes:
#include <cmath>

#include "config.h"

// Includes for random number generator
#include <random>

// Includes from nestkernel:
#include "structural_plasticity_node.h"
#include "connection.h"
#include "dict_util.h"
#include "event.h"
#include "nest_types.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"

// Includes from sli:
#include "dictdatum.h"

namespace nest
{
namespace diff_neuron_nestml_names
{
    const Name _in_rate( "in_rate" );
    const Name _in_rate_pre( "in_rate_pre" );
    const Name _out_rate( "out_rate" );
    const Name _spike_count_in( "spike_count_in" );
    const Name _spike_count_in_pos( "spike_count_in_pos" );
    const Name _spike_count_in_neg( "spike_count_in_neg" );
    const Name _spike_count_out( "spike_count_out" );
    const Name _tick( "tick" );
    const Name _lambda_poisson( "lambda_poisson" );
    const Name _spikes_buffer( "spikes_buffer" );
    const Name _kp( "kp" );
    const Name _pos( "pos" );
    const Name _base_rate( "base_rate" );
    const Name _buffer_size( "buffer_size" );
    const Name _simulation_steps( "simulation_steps" );
    const Name _time_wait( "time_wait" );
}
}




#include "nest_time.h"
  typedef size_t nest_port_t;
  typedef size_t nest_rport_t;

/* BeginDocumentation
  Name: diff_neuron_nestml

  Description:

    

  Parameters:
  The following parameters can be set in the status dictionary.
kp [real]  Gain
pos [boolean]  Sign sensitivity of the neuron
base_rate [Hz]  Base firing rate
buffer_size [ms]  Size of the sliding window
simulation_steps [integer]  Number of simulation steps (simulation_time/resolution())


  Dynamic state variables:
in_rate [Hz]  Input firing rate: to be computed from spikes
out_rate [Hz]  Output firing rate: defined accordingly to the input firing rate
spike_count_in [real]  Total incoming spikes (both excitatory and inhibitory)
spike_count_in_pos [real]  Incoming excitatory spikes
spike_count_in_neg [real]  Incoming inhibitory spikes
spike_count_out [integer]  Outgoing spikes
tick [integer]  Tick 
lambda_poisson [real]  Parameter of the Poisson distribution defining generator behavior
spikes_buffer [real]  Buffer for incoming spikes


  Sends: nest::SpikeEvent

  Receives: Spike,  DataLoggingRequest
*/

// Register the neuron model
void register_diff_neuron_nestml( const std::string& name );

class diff_neuron_nestml : public nest::StructuralPlasticityNode
{
public:
  /**
   * The constructor is only used to create the model prototype in the model manager.
  **/
  diff_neuron_nestml();

  /**
   * The copy constructor is used to create model copies and instances of the model.
   * @node The copy constructor needs to initialize the parameters and the state.
   *       Initialization of buffers and interal variables is deferred to
   *       @c init_buffers_() and @c pre_run_hook() (or calibrate() in NEST 3.3 and older).
  **/
  diff_neuron_nestml(const diff_neuron_nestml &);

  /**
   * Destructor.
  **/
  ~diff_neuron_nestml() override;

  // -------------------------------------------------------------------------
  //   Import sets of overloaded virtual functions.
  //   See: Technical Issues / Virtual Functions: Overriding, Overloading,
  //        and Hiding
  // -------------------------------------------------------------------------

  using nest::Node::handles_test_event;
  using nest::Node::handle;

  /**
   * Used to validate that we can send nest::SpikeEvent to desired target:port.
  **/
  nest_port_t send_test_event(nest::Node& target, nest_rport_t receptor_type, nest::synindex, bool) override;


  // -------------------------------------------------------------------------
  //   Functions handling incoming events.
  //   We tell nest that we can handle incoming events of various types by
  //   defining handle() for the given event.
  // -------------------------------------------------------------------------


  void handle(nest::SpikeEvent &) override;        //! accept spikes

  void handle(nest::DataLoggingRequest &) override;//! allow recording with multimeter
  nest_port_t handles_test_event(nest::SpikeEvent&, nest_port_t) override;
  nest_port_t handles_test_event(nest::DataLoggingRequest&, nest_port_t) override;

  // -------------------------------------------------------------------------
  //   Functions for getting/setting parameters and state values.
  // -------------------------------------------------------------------------

  void get_status(DictionaryDatum &) const override;
  void set_status(const DictionaryDatum &) override;


  // -------------------------------------------------------------------------
  //   Getters/setters for state block
  // -------------------------------------------------------------------------

  inline double get_in_rate() const
  {
    return S_.in_rate;
  }

  inline void set_in_rate(const double __v)
  {
    S_.in_rate = __v;
  }

  inline double get_in_rate_pre() const
  {
    return S_.in_rate_pre;
  }

  inline void set_in_rate_pre(const double __v)
  {
    S_.in_rate_pre = __v;
  }

  inline double get_out_rate() const
  {
    return S_.out_rate;
  }

  inline void set_out_rate(const double __v)
  {
    S_.out_rate = __v;
  }

  inline double get_spike_count_in() const
  {
    return S_.spike_count_in;
  }

  inline void set_spike_count_in(const double __v)
  {
    S_.spike_count_in = __v;
  }

  inline double get_spike_count_in_pos() const
  {
    return S_.spike_count_in_pos;
  }

  inline void set_spike_count_in_pos(const double __v)
  {
    S_.spike_count_in_pos = __v;
  }

  inline double get_spike_count_in_neg() const
  {
    return S_.spike_count_in_neg;
  }

  inline void set_spike_count_in_neg(const double __v)
  {
    S_.spike_count_in_neg = __v;
  }

  inline long get_spike_count_out() const
  {
    return S_.spike_count_out;
  }

  inline void set_spike_count_out(const long __v)
  {
    S_.spike_count_out = __v;
  }

  inline long get_tick() const
  {
    return S_.tick;
  }

  inline void set_tick(const long __v)
  {
    S_.tick = __v;
  }

  inline double get_lambda_poisson() const
  {
    return S_.lambda_poisson;
  }

  inline void set_lambda_poisson(const double __v)
  {
    S_.lambda_poisson = __v;
  }

  inline std::vector< double >  get_spikes_buffer() const
  {
    return S_.spikes_buffer;
  }

  inline void set_spikes_buffer(const std::vector< double >  __v)
  {
    S_.spikes_buffer = __v;
  }


  // -------------------------------------------------------------------------
  //   Getters/setters for parameters
  // -------------------------------------------------------------------------

  inline double get_kp() const
  {
    return P_.kp;
  }

  inline void set_kp(const double __v)
  {
    P_.kp = __v;
  }

  inline bool get_pos() const
  {
    return P_.pos;
  }

  inline void set_pos(const bool __v)
  {
    P_.pos = __v;
  }

  inline double get_base_rate() const
  {
    return P_.base_rate;
  }

  inline void set_base_rate(const double __v)
  {
    P_.base_rate = __v;
  }

  inline double get_buffer_size() const
  {
    return P_.buffer_size;
  }

  inline void set_buffer_size(const double __v)
  {
    P_.buffer_size = __v;
  }

  inline long get_simulation_steps() const
  {
    return P_.simulation_steps;
  }

  inline void set_simulation_steps(const long __v)
  {
    P_.simulation_steps = __v;
  }

  inline double get_time_wait() const
  {
    return P_.time_wait;
  }

  inline void set_time_wait(const double __v)
  {
    P_.time_wait = __v;
  }


  // -------------------------------------------------------------------------
  //   Getters/setters for internals
  // -------------------------------------------------------------------------

  inline double get_res() const
  {
    return V_.res;
  }

  inline void set_res(const double __v)
  {
    V_.res = __v;
  }
  inline long get_window_counts() const
  {
    return V_.window_counts;
  }

  inline void set_window_counts(const long __v)
  {
    V_.window_counts = __v;
  }
  inline double get___h() const
  {
    return V_.__h;
  }

  inline void set___h(const double __v)
  {
    V_.__h = __v;
  }


  // -------------------------------------------------------------------------
  //   Methods corresponding to event handlers
  // -------------------------------------------------------------------------

  

  // -------------------------------------------------------------------------
  //   Initialization functions
  // -------------------------------------------------------------------------
  void calibrate_time( const nest::TimeConverter& tc ) override;

protected:

private:
  void recompute_internal_variables(bool exclude_timestep=false);

private:

  static const nest_port_t MIN_SPIKE_RECEPTOR = 0;
  static const nest_port_t PORT_NOT_AVAILABLE = -1;

  enum SynapseTypes
  {
    SPIKES = 0,
    MAX_SPIKE_RECEPTOR = 1
  };

  static const size_t NUM_SPIKE_RECEPTORS = MAX_SPIKE_RECEPTOR - MIN_SPIKE_RECEPTOR;



  /**
   * Reset state of neuron.
  **/

  void init_state_internal_();

  /**
   * Reset internal buffers of neuron.
  **/
  void init_buffers_() override;

  /**
   * Initialize auxiliary quantities, leave parameters and state untouched.
  **/
  void pre_run_hook() override;

  /**
   * Take neuron through given time interval
  **/
  void update(nest::Time const &, const long, const long) override;

  // The next two classes need to be friends to access the State_ class/member
  friend class nest::DynamicRecordablesMap< diff_neuron_nestml >;
  friend class nest::DynamicUniversalDataLogger< diff_neuron_nestml >;
  friend class nest::DataAccessFunctor< diff_neuron_nestml >;

  /**
   * Free parameters of the neuron.
   *


   *
   * These are the parameters that can be set by the user through @c `node.set()`.
   * They are initialized from the model prototype when the node is created.
   * Parameters do not change during calls to @c update() and are not reset by
   * @c ResetNetwork.
   *
   * @note Parameters_ need neither copy constructor nor @c operator=(), since
   *       all its members are copied properly by the default copy constructor
   *       and assignment operator. Important:
   *       - If Parameters_ contained @c Time members, you need to define the
   *         assignment operator to recalibrate all members of type @c Time . You
   *         may also want to define the assignment operator.
   *       - If Parameters_ contained members that cannot copy themselves, such
   *         as C-style arrays, you need to define the copy constructor and
   *         assignment operator to copy those members.
  **/
  struct Parameters_
  {    
    //!  Gain
    double kp;
    //!  Sign sensitivity of the neuron
    bool pos;
    //!  Base firing rate
    double base_rate;
    //!  Size of the sliding window
    double buffer_size;
    //!  Number of simulation steps (simulation_time/resolution())
    long simulation_steps;
    double time_wait;

    /**
     * Initialize parameters to their default values.
    **/
    Parameters_();
  };

  /**
   * Dynamic state of the neuron.
   *
   *
   *
   * These are the state variables that are advanced in time by calls to
   * @c update(). In many models, some or all of them can be set by the user
   * through @c `node.set()`. The state variables are initialized from the model
   * prototype when the node is created. State variables are reset by @c ResetNetwork.
   *
   * @note State_ need neither copy constructor nor @c operator=(), since
   *       all its members are copied properly by the default copy constructor
   *       and assignment operator. Important:
   *       - If State_ contained @c Time members, you need to define the
   *         assignment operator to recalibrate all members of type @c Time . You
   *         may also want to define the assignment operator.
   *       - If State_ contained members that cannot copy themselves, such
   *         as C-style arrays, you need to define the copy constructor and
   *         assignment operator to copy those members.
  **/
  struct State_
  {
enum StateVecVars {
    IN_RATE = 0,
    IN_RATE_PRE = 1,
    OUT_RATE = 2,
    SPIKE_COUNT_IN = 3,
    SPIKE_COUNT_IN_POS = 4,
    SPIKE_COUNT_IN_NEG = 5,
    LAMBDA_POISSON = 6,
    SPIKES_BUFFER = 7,
};    
    //!  Input firing rate: to be computed from spikes
    double in_rate;
    double in_rate_pre;
    //!  Output firing rate: defined accordingly to the input firing rate
    double out_rate;
    //!  Total incoming spikes (both excitatory and inhibitory)
    double spike_count_in;
    //!  Incoming excitatory spikes
    double spike_count_in_pos;
    //!  Incoming inhibitory spikes
    double spike_count_in_neg;
    //!  Outgoing spikes
    long spike_count_out;
    //!  Tick 
    long tick;
    //!  Parameter of the Poisson distribution defining generator behavior
    double lambda_poisson;
    //!  Buffer for incoming spikes
    std::vector< double >  spikes_buffer;

    State_();
  };

  struct DelayedVariables_
  {
  };

  /**
   * Internal variables of the neuron.
   *
   *
   *
   * These variables must be initialized by @c pre_run_hook (or calibrate in NEST 3.3 and older), which is called before
   * the first call to @c update() upon each call to @c Simulate.
   * @node Variables_ needs neither constructor, copy constructor or assignment operator,
   *       since it is initialized by @c pre_run_hook() (or calibrate() in NEST 3.3 and older). If Variables_ has members that
   *       cannot destroy themselves, Variables_ will need a destructor.
  **/
  struct Variables_
  {
    double res;
    //!  Number of ticks corresponding to the window size
    long window_counts;
    double __h;
  };

  /**
   * Buffers of the neuron.
   * Usually buffers for incoming spikes and data logged for analog recorders.
   * Buffers must be initialized by @c init_buffers_(), which is called before
   * @c pre_run_hook() (or calibrate() in NEST 3.3 and older) on the first call to @c Simulate after the start of NEST,
   * ResetKernel or ResetNetwork.
   * @node Buffers_ needs neither constructor, copy constructor or assignment operator,
   *       since it is initialized by @c init_nodes_(). If Buffers_ has members that
   *       cannot destroy themselves, Buffers_ will need a destructor.
  **/
  struct Buffers_
  {
    Buffers_(diff_neuron_nestml &);
    Buffers_(const Buffers_ &, diff_neuron_nestml &);

    /**
     * Logger for all analog data
    **/
    nest::DynamicUniversalDataLogger<diff_neuron_nestml> logger_;

    // -----------------------------------------------------------------------
    //   Spike buffers and sums of incoming spikes/currents per timestep
    // -----------------------------------------------------------------------    



    /**
     * Buffer containing the incoming spikes
    **/
    inline std::vector< nest::RingBuffer >& get_spike_inputs_()
    {
        return spike_inputs_;
    }
    std::vector< nest::RingBuffer > spike_inputs_;

    /**
     * Buffer containing the sum of all the incoming spikes
    **/
    inline std::vector< double >& get_spike_inputs_grid_sum_()
    {
        return spike_inputs_grid_sum_;
    }
    std::vector< double > spike_inputs_grid_sum_;

    /**
     * Buffer containing a flag whether incoming spikes have been received on a given port
    **/
    inline std::vector< nest::RingBuffer >& get_spike_input_received_()
    {
        return spike_input_received_;
    }
    std::vector< nest::RingBuffer > spike_input_received_;

    /**
     * Buffer containing a flag whether incoming spikes have been received on a given port
    **/
    inline std::vector< double >& get_spike_input_received_grid_sum_()
    {
        return spike_input_received_grid_sum_;
    }
    std::vector< double > spike_input_received_grid_sum_;

    // -----------------------------------------------------------------------
    //   Continuous-input buffers
    // -----------------------------------------------------------------------

    
  };

  // -------------------------------------------------------------------------
  //   Getters/setters for inline expressions
  // -------------------------------------------------------------------------
  

  // -------------------------------------------------------------------------
  //   Getters/setters for input buffers
  // -------------------------------------------------------------------------  




  /**
   * Buffer containing the incoming spikes
  **/
  inline std::vector< nest::RingBuffer >& get_spike_inputs_()
  {
      return B_.get_spike_inputs_();
  }

  /**
   * Buffer containing the sum of all the incoming spikes
  **/
  inline std::vector< double >& get_spike_inputs_grid_sum_()
  {
      return B_.get_spike_inputs_grid_sum_();
  }

  /**
   * Buffer containing a flag whether incoming spikes have been received on a given port
  **/
  inline std::vector< nest::RingBuffer >& get_spike_input_received_()
  {
      return B_.get_spike_input_received_();
  }

  /**
   * Buffer containing a flag whether incoming spikes have been received on a given port
  **/
  inline std::vector< double >& get_spike_input_received_grid_sum_()
  {
      return B_.get_spike_input_received_grid_sum_();
  }

  // -------------------------------------------------------------------------
  //   Member variables of neuron model.
  //   Each model neuron should have precisely the following four data members,
  //   which are one instance each of the parameters, state, buffers and variables
  //   structures. Experience indicates that the state and variables member should
  //   be next to each other to achieve good efficiency (caching).
  //   Note: Devices require one additional data member, an instance of the
  //   ``Device`` child class they belong to.
  // -------------------------------------------------------------------------


  Parameters_       P_;        //!< Free parameters.
  State_            S_;        //!< Dynamic state.
  DelayedVariables_ DV_;       //!< Delayed state variables.
  Variables_        V_;        //!< Internal Variables
  Buffers_          B_;        //!< Buffers.

  //! Mapping of recordables names to access functions
  nest::DynamicRecordablesMap<diff_neuron_nestml> recordablesMap_;
  nest::DataAccessFunctor< diff_neuron_nestml > get_data_access_functor( size_t elem );
  std::string get_var_name(size_t elem, std::string var_name);
  void insert_recordables(size_t first=0);


inline double get_state_element(size_t elem)
  {
    if
    (elem == State_::IN_RATE)
    {
      return S_.in_rate;
    }
    else if
    (elem == State_::IN_RATE_PRE)
    {
      return S_.in_rate_pre;
    }
    else if
    (elem == State_::OUT_RATE)
    {
      return S_.out_rate;
    }
    else if
    (elem == State_::SPIKE_COUNT_IN)
    {
      return S_.spike_count_in;
    }
    else if
    (elem == State_::SPIKE_COUNT_IN_POS)
    {
      return S_.spike_count_in_pos;
    }
    else if
    (elem == State_::SPIKE_COUNT_IN_NEG)
    {
      return S_.spike_count_in_neg;
    }
    else if
    (elem == State_::LAMBDA_POISSON)
    {
      return S_.lambda_poisson;
    }
    else

    {
      return S_.spikes_buffer[ elem - State_::SPIKES_BUFFER ];
    }
  }
  nest::normal_distribution normal_dev_; //!< random deviate generator
  nest::poisson_distribution poisson_dev_; //!< random deviate generator

}; /* neuron diff_neuron_nestml */

inline nest_port_t diff_neuron_nestml::send_test_event(nest::Node& target, nest_rport_t receptor_type, nest::synindex, bool)
{
  // You should usually not change the code in this function.
  // It confirms that the target of connection @c c accepts @c nest::SpikeEvent on
  // the given @c receptor_type.
  nest::SpikeEvent e;
  e.set_sender(*this);
  return target.handles_test_event(e, receptor_type);
}

inline nest_port_t diff_neuron_nestml::handles_test_event(nest::SpikeEvent&, nest_port_t receptor_type)
{
    // You should usually not change the code in this function.
    // It confirms to the connection management system that we are able
    // to handle @c SpikeEvent on port 0. You need to extend the function
    // if you want to differentiate between input ports.
    if (receptor_type != 0)
    {
      throw nest::UnknownReceptorType(receptor_type, get_name());
    }
    return 0;
}

inline nest_port_t diff_neuron_nestml::handles_test_event(nest::DataLoggingRequest& dlr, nest_port_t receptor_type)
{
  // You should usually not change the code in this function.
  // It confirms to the connection management system that we are able
  // to handle @c DataLoggingRequest on port 0.
  // The function also tells the built-in UniversalDataLogger that this node
  // is recorded from and that it thus needs to collect data during simulation.
  if (receptor_type != 0)
  {
    throw nest::UnknownReceptorType(receptor_type, get_name());
  }

  return B_.logger_.connect_logging_device(dlr, recordablesMap_);
}

inline void diff_neuron_nestml::get_status(DictionaryDatum &__d) const
{
  // parameters
  def< double >(__d, nest::diff_neuron_nestml_names::_kp, get_kp());
  def< bool >(__d, nest::diff_neuron_nestml_names::_pos, get_pos());
  def< double >(__d, nest::diff_neuron_nestml_names::_base_rate, get_base_rate());
  def< double >(__d, nest::diff_neuron_nestml_names::_buffer_size, get_buffer_size());
  def< long >(__d, nest::diff_neuron_nestml_names::_simulation_steps, get_simulation_steps());
  def< double >(__d, nest::diff_neuron_nestml_names::_time_wait, get_time_wait());

  // initial values for state variables in ODE or kernel
  def< double >(__d, nest::diff_neuron_nestml_names::_in_rate, get_in_rate());
  def< double >(__d, nest::diff_neuron_nestml_names::_in_rate_pre, get_in_rate_pre());
  def< double >(__d, nest::diff_neuron_nestml_names::_out_rate, get_out_rate());
  def< double >(__d, nest::diff_neuron_nestml_names::_spike_count_in, get_spike_count_in());
  def< double >(__d, nest::diff_neuron_nestml_names::_spike_count_in_pos, get_spike_count_in_pos());
  def< double >(__d, nest::diff_neuron_nestml_names::_spike_count_in_neg, get_spike_count_in_neg());
  def< long >(__d, nest::diff_neuron_nestml_names::_spike_count_out, get_spike_count_out());
  def< long >(__d, nest::diff_neuron_nestml_names::_tick, get_tick());
  def< double >(__d, nest::diff_neuron_nestml_names::_lambda_poisson, get_lambda_poisson());
  def< std::vector< double >  >(__d, nest::diff_neuron_nestml_names::_spikes_buffer, get_spikes_buffer());

  StructuralPlasticityNode::get_status( __d );

  (*__d)[nest::names::recordables] = recordablesMap_.get_list();
}

inline void diff_neuron_nestml::set_status(const DictionaryDatum &__d)
{
  // parameters
  double tmp_kp = get_kp();
  nest::updateValueParam<double>(__d, nest::diff_neuron_nestml_names::_kp, tmp_kp, this);
  // Resize vectors
  if (tmp_kp != get_kp())
  {
  }
  bool tmp_pos = get_pos();
  nest::updateValueParam<bool>(__d, nest::diff_neuron_nestml_names::_pos, tmp_pos, this);
  // Resize vectors
  if (tmp_pos != get_pos())
  {
  }
  double tmp_base_rate = get_base_rate();
  nest::updateValueParam<double>(__d, nest::diff_neuron_nestml_names::_base_rate, tmp_base_rate, this);
  // Resize vectors
  if (tmp_base_rate != get_base_rate())
  {
  }
  double tmp_buffer_size = get_buffer_size();
  nest::updateValueParam<double>(__d, nest::diff_neuron_nestml_names::_buffer_size, tmp_buffer_size, this);
  // Resize vectors
  if (tmp_buffer_size != get_buffer_size())
  {
  }
  long tmp_simulation_steps = get_simulation_steps();
  nest::updateValueParam<long>(__d, nest::diff_neuron_nestml_names::_simulation_steps, tmp_simulation_steps, this);
  // Resize vectors
  if (tmp_simulation_steps != get_simulation_steps())
  {
    std::vector< double >  _tmp_spikes_buffer = get_spikes_buffer();
    _tmp_spikes_buffer.resize(tmp_simulation_steps, 0.);
    set_spikes_buffer(_tmp_spikes_buffer);
  }
  double tmp_time_wait = get_time_wait();
  nest::updateValueParam<double>(__d, nest::diff_neuron_nestml_names::_time_wait, tmp_time_wait, this);
  // Resize vectors
  if (tmp_time_wait != get_time_wait())
  {
  }

  // initial values for state variables in ODE or kernel
  double tmp_in_rate = get_in_rate();
  nest::updateValueParam<double>(__d, nest::diff_neuron_nestml_names::_in_rate, tmp_in_rate, this);
  double tmp_in_rate_pre = get_in_rate_pre();
  nest::updateValueParam<double>(__d, nest::diff_neuron_nestml_names::_in_rate_pre, tmp_in_rate_pre, this);
  double tmp_out_rate = get_out_rate();
  nest::updateValueParam<double>(__d, nest::diff_neuron_nestml_names::_out_rate, tmp_out_rate, this);
  double tmp_spike_count_in = get_spike_count_in();
  nest::updateValueParam<double>(__d, nest::diff_neuron_nestml_names::_spike_count_in, tmp_spike_count_in, this);
  double tmp_spike_count_in_pos = get_spike_count_in_pos();
  nest::updateValueParam<double>(__d, nest::diff_neuron_nestml_names::_spike_count_in_pos, tmp_spike_count_in_pos, this);
  double tmp_spike_count_in_neg = get_spike_count_in_neg();
  nest::updateValueParam<double>(__d, nest::diff_neuron_nestml_names::_spike_count_in_neg, tmp_spike_count_in_neg, this);
  long tmp_spike_count_out = get_spike_count_out();
  nest::updateValueParam<long>(__d, nest::diff_neuron_nestml_names::_spike_count_out, tmp_spike_count_out, this);
  long tmp_tick = get_tick();
  nest::updateValueParam<long>(__d, nest::diff_neuron_nestml_names::_tick, tmp_tick, this);
  double tmp_lambda_poisson = get_lambda_poisson();
  nest::updateValueParam<double>(__d, nest::diff_neuron_nestml_names::_lambda_poisson, tmp_lambda_poisson, this);
  std::vector< double >  tmp_spikes_buffer = get_spikes_buffer();
  updateValue<std::vector< double > >(__d, nest::diff_neuron_nestml_names::_spikes_buffer, tmp_spikes_buffer);
   
  // Check if the new vector size matches its original size
  if ( tmp_spikes_buffer.size() != tmp_simulation_steps )
  {
    std::stringstream msg;
    msg << "The vector \"spikes_buffer\" does not match its size: " << tmp_simulation_steps;
    throw nest::BadProperty(msg.str());
  }

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  StructuralPlasticityNode::set_status(__d);

  // if we get here, temporaries contain consistent set of properties
  set_kp(tmp_kp);
  set_pos(tmp_pos);
  set_base_rate(tmp_base_rate);
  set_buffer_size(tmp_buffer_size);
  set_simulation_steps(tmp_simulation_steps);
  set_time_wait(tmp_time_wait);
  set_in_rate(tmp_in_rate);
  set_in_rate_pre(tmp_in_rate_pre);
  set_out_rate(tmp_out_rate);
  set_spike_count_in(tmp_spike_count_in);
  set_spike_count_in_pos(tmp_spike_count_in_pos);
  set_spike_count_in_neg(tmp_spike_count_in_neg);
  set_spike_count_out(tmp_spike_count_out);
  set_tick(tmp_tick);
  set_lambda_poisson(tmp_lambda_poisson);
  set_spikes_buffer(tmp_spikes_buffer);





  // recompute internal variables in case they are dependent on parameters or state that might have been updated in this call to set_status()
  recompute_internal_variables();
};



#endif /* #ifndef DIFF_NEURON_NESTML */
