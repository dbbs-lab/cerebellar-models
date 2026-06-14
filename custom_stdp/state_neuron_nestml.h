
/**
 *  state_neuron_nestml.h
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
#ifndef STATE_NEURON_NESTML
#define STATE_NEURON_NESTML

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
namespace state_neuron_nestml_names
{
    const Name _in_rate( "in_rate" );
    const Name _out_rate( "out_rate" );
    const Name _spike_count_out( "spike_count_out" );
    const Name _current_fbk_input( "current_fbk_input" );
    const Name _current_pred_input( "current_pred_input" );
    const Name _fbk_buffer( "fbk_buffer" );
    const Name _pred_buffer( "pred_buffer" );
    const Name _fbk_counts( "fbk_counts" );
    const Name _pred_counts( "pred_counts" );
    const Name _tick( "tick" );
    const Name _position_count( "position_count" );
    const Name _mean_fbk( "mean_fbk" );
    const Name _mean_pred( "mean_pred" );
    const Name _var_fbk( "var_fbk" );
    const Name _var_pred( "var_pred" );
    const Name _CV_fbk( "CV_fbk" );
    const Name _CV_pred( "CV_pred" );
    const Name _total_CV( "total_CV" );
    const Name _lambda_poisson( "lambda_poisson" );
    const Name _kp( "kp" );
    const Name _pos( "pos" );
    const Name _base_rate( "base_rate" );
    const Name _buffer_size( "buffer_size" );
    const Name _simulation_steps( "simulation_steps" );
    const Name _N_fbk( "N_fbk" );
    const Name _N_pred( "N_pred" );
    const Name _fbk_bf_size( "fbk_bf_size" );
    const Name _pred_bf_size( "pred_bf_size" );
    const Name _time_wait( "time_wait" );
    const Name _time_trial( "time_trial" );
}
}




#include "nest_time.h"
  typedef size_t nest_port_t;
  typedef size_t nest_rport_t;

/* BeginDocumentation
  Name: state_neuron_nestml

  Description:

    

  Parameters:
  The following parameters can be set in the status dictionary.
kp [real]  Gain
pos [boolean]  Sign sensitivity of the neuron
base_rate [Hz]  Base firing rate
buffer_size [ms]  Size of the sliding window
simulation_steps [integer]  Number of simulation steps (simulation_time/resolution())
N_fbk [integer]  Population size for sensory feedback
N_pred [integer]  Population size for sensory prediction


  Dynamic state variables:
in_rate [Hz]  Input firing rate: to be computed from spikes
out_rate [Hz]  Output firing rate: defined accordingly to the input firing rate
spike_count_out [integer]  Outgoing spikes
fbk_buffer [real]  Buffer for sensory feedback spikes
pred_buffer [real]  Buffer for sensory prediction spikes
fbk_counts [real]  Counts of incoming feedback spikes
pred_counts [real]  Counts of incoming prediction spikes
tick [integer]  Tick 
mean_fbk [real]  Mean sensory feedback
mean_pred [real]  Mean sensory prediction
var_fbk [real]  Variance of sensory feedback
var_pred [real]  Variance of sensory prediction
CV_fbk [real]  Coefficient of variation of sensory feedback
CV_pred [real]  Coefficient of variation of sensory prediction
lambda_poisson [real]  Parameter of the Poisson distribution defining generator behavior


  Sends: nest::SpikeEvent

  Receives: Spike,  DataLoggingRequest
*/

// Register the neuron model
void register_state_neuron_nestml( const std::string& name );

class state_neuron_nestml : public nest::StructuralPlasticityNode
{
public:
  /**
   * The constructor is only used to create the model prototype in the model manager.
  **/
  state_neuron_nestml();

  /**
   * The copy constructor is used to create model copies and instances of the model.
   * @node The copy constructor needs to initialize the parameters and the state.
   *       Initialization of buffers and interal variables is deferred to
   *       @c init_buffers_() and @c pre_run_hook() (or calibrate() in NEST 3.3 and older).
  **/
  state_neuron_nestml(const state_neuron_nestml &);

  /**
   * Destructor.
  **/
  ~state_neuron_nestml() override;

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

  inline double get_out_rate() const
  {
    return S_.out_rate;
  }

  inline void set_out_rate(const double __v)
  {
    S_.out_rate = __v;
  }

  inline long get_spike_count_out() const
  {
    return S_.spike_count_out;
  }

  inline void set_spike_count_out(const long __v)
  {
    S_.spike_count_out = __v;
  }

  inline std::vector< double >  get_current_fbk_input() const
  {
    return S_.current_fbk_input;
  }

  inline void set_current_fbk_input(const std::vector< double >  __v)
  {
    S_.current_fbk_input = __v;
  }

  inline std::vector< double >  get_current_pred_input() const
  {
    return S_.current_pred_input;
  }

  inline void set_current_pred_input(const std::vector< double >  __v)
  {
    S_.current_pred_input = __v;
  }

  inline std::vector< double >  get_fbk_buffer() const
  {
    return S_.fbk_buffer;
  }

  inline void set_fbk_buffer(const std::vector< double >  __v)
  {
    S_.fbk_buffer = __v;
  }

  inline std::vector< double >  get_pred_buffer() const
  {
    return S_.pred_buffer;
  }

  inline void set_pred_buffer(const std::vector< double >  __v)
  {
    S_.pred_buffer = __v;
  }

  inline std::vector< double >  get_fbk_counts() const
  {
    return S_.fbk_counts;
  }

  inline void set_fbk_counts(const std::vector< double >  __v)
  {
    S_.fbk_counts = __v;
  }

  inline std::vector< double >  get_pred_counts() const
  {
    return S_.pred_counts;
  }

  inline void set_pred_counts(const std::vector< double >  __v)
  {
    S_.pred_counts = __v;
  }

  inline long get_tick() const
  {
    return S_.tick;
  }

  inline void set_tick(const long __v)
  {
    S_.tick = __v;
  }

  inline long get_position_count() const
  {
    return S_.position_count;
  }

  inline void set_position_count(const long __v)
  {
    S_.position_count = __v;
  }

  inline double get_mean_fbk() const
  {
    return S_.mean_fbk;
  }

  inline void set_mean_fbk(const double __v)
  {
    S_.mean_fbk = __v;
  }

  inline double get_mean_pred() const
  {
    return S_.mean_pred;
  }

  inline void set_mean_pred(const double __v)
  {
    S_.mean_pred = __v;
  }

  inline double get_var_fbk() const
  {
    return S_.var_fbk;
  }

  inline void set_var_fbk(const double __v)
  {
    S_.var_fbk = __v;
  }

  inline double get_var_pred() const
  {
    return S_.var_pred;
  }

  inline void set_var_pred(const double __v)
  {
    S_.var_pred = __v;
  }

  inline double get_CV_fbk() const
  {
    return S_.CV_fbk;
  }

  inline void set_CV_fbk(const double __v)
  {
    S_.CV_fbk = __v;
  }

  inline double get_CV_pred() const
  {
    return S_.CV_pred;
  }

  inline void set_CV_pred(const double __v)
  {
    S_.CV_pred = __v;
  }

  inline double get_total_CV() const
  {
    return S_.total_CV;
  }

  inline void set_total_CV(const double __v)
  {
    S_.total_CV = __v;
  }

  inline double get_lambda_poisson() const
  {
    return S_.lambda_poisson;
  }

  inline void set_lambda_poisson(const double __v)
  {
    S_.lambda_poisson = __v;
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

  inline long get_N_fbk() const
  {
    return P_.N_fbk;
  }

  inline void set_N_fbk(const long __v)
  {
    P_.N_fbk = __v;
  }

  inline long get_N_pred() const
  {
    return P_.N_pred;
  }

  inline void set_N_pred(const long __v)
  {
    P_.N_pred = __v;
  }

  inline long get_fbk_bf_size() const
  {
    return P_.fbk_bf_size;
  }

  inline void set_fbk_bf_size(const long __v)
  {
    P_.fbk_bf_size = __v;
  }

  inline long get_pred_bf_size() const
  {
    return P_.pred_bf_size;
  }

  inline void set_pred_bf_size(const long __v)
  {
    P_.pred_bf_size = __v;
  }

  inline double get_time_wait() const
  {
    return P_.time_wait;
  }

  inline void set_time_wait(const double __v)
  {
    P_.time_wait = __v;
  }

  inline double get_time_trial() const
  {
    return P_.time_trial;
  }

  inline void set_time_trial(const double __v)
  {
    P_.time_trial = __v;
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
  inline long get_buffer_steps() const
  {
    return V_.buffer_steps;
  }

  inline void set_buffer_steps(const long __v)
  {
    V_.buffer_steps = __v;
  }
  inline long get_trial_steps() const
  {
    return V_.trial_steps;
  }

  inline void set_trial_steps(const long __v)
  {
    V_.trial_steps = __v;
  }
  inline long get_wait_steps() const
  {
    return V_.wait_steps;
  }

  inline void set_wait_steps(const long __v)
  {
    V_.wait_steps = __v;
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
/**
   * Synapse types to connect to
   * @note Excluded lower and upper bounds are defined as MIN_, MAX_.
   *       Excluding port 0 avoids accidental connections.
  **/
  static const nest_port_t MIN_SPIKE_RECEPTOR = 1;
  static const nest_port_t PORT_NOT_AVAILABLE = -1;

  enum SynapseTypes
  {
    FBK_SPIKES_0 = 1,
    FBK_SPIKES_1 = 2,
    FBK_SPIKES_2 = 3,
    FBK_SPIKES_3 = 4,
    FBK_SPIKES_4 = 5,
    FBK_SPIKES_5 = 6,
    FBK_SPIKES_6 = 7,
    FBK_SPIKES_7 = 8,
    FBK_SPIKES_8 = 9,
    FBK_SPIKES_9 = 10,
    FBK_SPIKES_10 = 11,
    FBK_SPIKES_11 = 12,
    FBK_SPIKES_12 = 13,
    FBK_SPIKES_13 = 14,
    FBK_SPIKES_14 = 15,
    FBK_SPIKES_15 = 16,
    FBK_SPIKES_16 = 17,
    FBK_SPIKES_17 = 18,
    FBK_SPIKES_18 = 19,
    FBK_SPIKES_19 = 20,
    FBK_SPIKES_20 = 21,
    FBK_SPIKES_21 = 22,
    FBK_SPIKES_22 = 23,
    FBK_SPIKES_23 = 24,
    FBK_SPIKES_24 = 25,
    FBK_SPIKES_25 = 26,
    FBK_SPIKES_26 = 27,
    FBK_SPIKES_27 = 28,
    FBK_SPIKES_28 = 29,
    FBK_SPIKES_29 = 30,
    FBK_SPIKES_30 = 31,
    FBK_SPIKES_31 = 32,
    FBK_SPIKES_32 = 33,
    FBK_SPIKES_33 = 34,
    FBK_SPIKES_34 = 35,
    FBK_SPIKES_35 = 36,
    FBK_SPIKES_36 = 37,
    FBK_SPIKES_37 = 38,
    FBK_SPIKES_38 = 39,
    FBK_SPIKES_39 = 40,
    FBK_SPIKES_40 = 41,
    FBK_SPIKES_41 = 42,
    FBK_SPIKES_42 = 43,
    FBK_SPIKES_43 = 44,
    FBK_SPIKES_44 = 45,
    FBK_SPIKES_45 = 46,
    FBK_SPIKES_46 = 47,
    FBK_SPIKES_47 = 48,
    FBK_SPIKES_48 = 49,
    FBK_SPIKES_49 = 50,
    PRED_SPIKES_0 = 51,
    PRED_SPIKES_1 = 52,
    PRED_SPIKES_2 = 53,
    PRED_SPIKES_3 = 54,
    PRED_SPIKES_4 = 55,
    PRED_SPIKES_5 = 56,
    PRED_SPIKES_6 = 57,
    PRED_SPIKES_7 = 58,
    PRED_SPIKES_8 = 59,
    PRED_SPIKES_9 = 60,
    PRED_SPIKES_10 = 61,
    PRED_SPIKES_11 = 62,
    PRED_SPIKES_12 = 63,
    PRED_SPIKES_13 = 64,
    PRED_SPIKES_14 = 65,
    PRED_SPIKES_15 = 66,
    PRED_SPIKES_16 = 67,
    PRED_SPIKES_17 = 68,
    PRED_SPIKES_18 = 69,
    PRED_SPIKES_19 = 70,
    PRED_SPIKES_20 = 71,
    PRED_SPIKES_21 = 72,
    PRED_SPIKES_22 = 73,
    PRED_SPIKES_23 = 74,
    PRED_SPIKES_24 = 75,
    PRED_SPIKES_25 = 76,
    PRED_SPIKES_26 = 77,
    PRED_SPIKES_27 = 78,
    PRED_SPIKES_28 = 79,
    PRED_SPIKES_29 = 80,
    PRED_SPIKES_30 = 81,
    PRED_SPIKES_31 = 82,
    PRED_SPIKES_32 = 83,
    PRED_SPIKES_33 = 84,
    PRED_SPIKES_34 = 85,
    PRED_SPIKES_35 = 86,
    PRED_SPIKES_36 = 87,
    PRED_SPIKES_37 = 88,
    PRED_SPIKES_38 = 89,
    PRED_SPIKES_39 = 90,
    PRED_SPIKES_40 = 91,
    PRED_SPIKES_41 = 92,
    PRED_SPIKES_42 = 93,
    PRED_SPIKES_43 = 94,
    PRED_SPIKES_44 = 95,
    PRED_SPIKES_45 = 96,
    PRED_SPIKES_46 = 97,
    PRED_SPIKES_47 = 98,
    PRED_SPIKES_48 = 99,
    PRED_SPIKES_49 = 100,
    MAX_SPIKE_RECEPTOR = 101
  };

  static const size_t NUM_SPIKE_RECEPTORS = MAX_SPIKE_RECEPTOR - MIN_SPIKE_RECEPTOR;

static std::vector< std::tuple< int, int > > rport_to_nestml_buffer_idx;

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
  friend class nest::DynamicRecordablesMap< state_neuron_nestml >;
  friend class nest::DynamicUniversalDataLogger< state_neuron_nestml >;
  friend class nest::DataAccessFunctor< state_neuron_nestml >;

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
    //!  Population size for sensory feedback
    long N_fbk;
    //!  Population size for sensory prediction
    long N_pred;
    long fbk_bf_size;
    long pred_bf_size;
    double time_wait;
    double time_trial;

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
    OUT_RATE = 1,
    CURRENT_FBK_INPUT = 2,
    CURRENT_PRED_INPUT = 52,
    FBK_BUFFER = 102,
    PRED_BUFFER = 10102,
    FBK_COUNTS = 20102,
    PRED_COUNTS = 20152,
    MEAN_FBK = 20202,
    MEAN_PRED = 20203,
    VAR_FBK = 20204,
    VAR_PRED = 20205,
    CV_FBK = 20206,
    CV_PRED = 20207,
    TOTAL_CV = 20208,
    LAMBDA_POISSON = 20209,
};    
    //!  Input firing rate: to be computed from spikes
    double in_rate;
    //!  Output firing rate: defined accordingly to the input firing rate
    double out_rate;
    //!  Outgoing spikes
    long spike_count_out;
    std::vector< double >  current_fbk_input;
    std::vector< double >  current_pred_input;
    //!  Buffer for sensory feedback spikes
    std::vector< double >  fbk_buffer;
    //!  Buffer for sensory prediction spikes
    std::vector< double >  pred_buffer;
    //!  Counts of incoming feedback spikes
    std::vector< double >  fbk_counts;
    //!  Counts of incoming prediction spikes
    std::vector< double >  pred_counts;
    //!  Tick 
    long tick;
    long position_count;
    //!  Mean sensory feedback
    double mean_fbk;
    //!  Mean sensory prediction
    double mean_pred;
    //!  Variance of sensory feedback
    double var_fbk;
    //!  Variance of sensory prediction
    double var_pred;
    //!  Coefficient of variation of sensory feedback
    double CV_fbk;
    //!  Coefficient of variation of sensory prediction
    double CV_pred;
    double total_CV;
    //!  Parameter of the Poisson distribution defining generator behavior
    double lambda_poisson;

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
    long buffer_steps;
    long trial_steps;
    long wait_steps;
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
    Buffers_(state_neuron_nestml &);
    Buffers_(const Buffers_ &, state_neuron_nestml &);

    /**
     * Logger for all analog data
    **/
    nest::DynamicUniversalDataLogger<state_neuron_nestml> logger_;

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
  nest::DynamicRecordablesMap<state_neuron_nestml> recordablesMap_;
  nest::DataAccessFunctor< state_neuron_nestml > get_data_access_functor( size_t elem );
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
    (elem == State_::OUT_RATE)
    {
      return S_.out_rate;
    }
    else if(elem >= State_::CURRENT_FBK_INPUT and elem < State_::CURRENT_FBK_INPUT + 
P_.N_fbk)
    {
      return S_.current_fbk_input[ elem - State_::CURRENT_FBK_INPUT ];
    }
    else if(elem >= State_::CURRENT_PRED_INPUT and elem < State_::CURRENT_PRED_INPUT + 
P_.N_pred)
    {
      return S_.current_pred_input[ elem - State_::CURRENT_PRED_INPUT ];
    }
    else if(elem >= State_::FBK_BUFFER and elem < State_::FBK_BUFFER + 
P_.fbk_bf_size)
    {
      return S_.fbk_buffer[ elem - State_::FBK_BUFFER ];
    }
    else if(elem >= State_::PRED_BUFFER and elem < State_::PRED_BUFFER + 
P_.pred_bf_size)
    {
      return S_.pred_buffer[ elem - State_::PRED_BUFFER ];
    }
    else if(elem >= State_::FBK_COUNTS and elem < State_::FBK_COUNTS + 
P_.N_fbk)
    {
      return S_.fbk_counts[ elem - State_::FBK_COUNTS ];
    }
    else if(elem >= State_::PRED_COUNTS and elem < State_::PRED_COUNTS + 
P_.N_pred)
    {
      return S_.pred_counts[ elem - State_::PRED_COUNTS ];
    }
    else if
    (elem == State_::MEAN_FBK)
    {
      return S_.mean_fbk;
    }
    else if
    (elem == State_::MEAN_PRED)
    {
      return S_.mean_pred;
    }
    else if
    (elem == State_::VAR_FBK)
    {
      return S_.var_fbk;
    }
    else if
    (elem == State_::VAR_PRED)
    {
      return S_.var_pred;
    }
    else if
    (elem == State_::CV_FBK)
    {
      return S_.CV_fbk;
    }
    else if
    (elem == State_::CV_PRED)
    {
      return S_.CV_pred;
    }
    else if
    (elem == State_::TOTAL_CV)
    {
      return S_.total_CV;
    }
    else
    {
      return S_.lambda_poisson;
    }
  }
  nest::normal_distribution normal_dev_; //!< random deviate generator
  nest::poisson_distribution poisson_dev_; //!< random deviate generator

}; /* neuron state_neuron_nestml */

inline nest_port_t state_neuron_nestml::send_test_event(nest::Node& target, nest_rport_t receptor_type, nest::synindex, bool)
{
  // You should usually not change the code in this function.
  // It confirms that the target of connection @c c accepts @c nest::SpikeEvent on
  // the given @c receptor_type.
  nest::SpikeEvent e;
  e.set_sender(*this);
  return target.handles_test_event(e, receptor_type);
}

inline nest_port_t state_neuron_nestml::handles_test_event(nest::SpikeEvent&, nest_port_t receptor_type)
{
    assert( B_.spike_inputs_.size() == NUM_SPIKE_RECEPTORS );
    if ( receptor_type < MIN_SPIKE_RECEPTOR or receptor_type >= MAX_SPIKE_RECEPTOR )
    {
      throw nest::UnknownReceptorType( receptor_type, get_name() );
    }
    return receptor_type - MIN_SPIKE_RECEPTOR;
}

inline nest_port_t state_neuron_nestml::handles_test_event(nest::DataLoggingRequest& dlr, nest_port_t receptor_type)
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

inline void state_neuron_nestml::get_status(DictionaryDatum &__d) const
{
  // parameters
  def< double >(__d, nest::state_neuron_nestml_names::_kp, get_kp());
  def< bool >(__d, nest::state_neuron_nestml_names::_pos, get_pos());
  def< double >(__d, nest::state_neuron_nestml_names::_base_rate, get_base_rate());
  def< double >(__d, nest::state_neuron_nestml_names::_buffer_size, get_buffer_size());
  def< long >(__d, nest::state_neuron_nestml_names::_simulation_steps, get_simulation_steps());
  def< long >(__d, nest::state_neuron_nestml_names::_N_fbk, get_N_fbk());
  def< long >(__d, nest::state_neuron_nestml_names::_N_pred, get_N_pred());
  def< long >(__d, nest::state_neuron_nestml_names::_fbk_bf_size, get_fbk_bf_size());
  def< long >(__d, nest::state_neuron_nestml_names::_pred_bf_size, get_pred_bf_size());
  def< double >(__d, nest::state_neuron_nestml_names::_time_wait, get_time_wait());
  def< double >(__d, nest::state_neuron_nestml_names::_time_trial, get_time_trial());

  // initial values for state variables in ODE or kernel
  def< double >(__d, nest::state_neuron_nestml_names::_in_rate, get_in_rate());
  def< double >(__d, nest::state_neuron_nestml_names::_out_rate, get_out_rate());
  def< long >(__d, nest::state_neuron_nestml_names::_spike_count_out, get_spike_count_out());
  def< std::vector< double >  >(__d, nest::state_neuron_nestml_names::_current_fbk_input, get_current_fbk_input());
  def< std::vector< double >  >(__d, nest::state_neuron_nestml_names::_current_pred_input, get_current_pred_input());
  def< std::vector< double >  >(__d, nest::state_neuron_nestml_names::_fbk_buffer, get_fbk_buffer());
  def< std::vector< double >  >(__d, nest::state_neuron_nestml_names::_pred_buffer, get_pred_buffer());
  def< std::vector< double >  >(__d, nest::state_neuron_nestml_names::_fbk_counts, get_fbk_counts());
  def< std::vector< double >  >(__d, nest::state_neuron_nestml_names::_pred_counts, get_pred_counts());
  def< long >(__d, nest::state_neuron_nestml_names::_tick, get_tick());
  def< long >(__d, nest::state_neuron_nestml_names::_position_count, get_position_count());
  def< double >(__d, nest::state_neuron_nestml_names::_mean_fbk, get_mean_fbk());
  def< double >(__d, nest::state_neuron_nestml_names::_mean_pred, get_mean_pred());
  def< double >(__d, nest::state_neuron_nestml_names::_var_fbk, get_var_fbk());
  def< double >(__d, nest::state_neuron_nestml_names::_var_pred, get_var_pred());
  def< double >(__d, nest::state_neuron_nestml_names::_CV_fbk, get_CV_fbk());
  def< double >(__d, nest::state_neuron_nestml_names::_CV_pred, get_CV_pred());
  def< double >(__d, nest::state_neuron_nestml_names::_total_CV, get_total_CV());
  def< double >(__d, nest::state_neuron_nestml_names::_lambda_poisson, get_lambda_poisson());

  StructuralPlasticityNode::get_status( __d );
  DictionaryDatum __receptor_type = new Dictionary();
    ( *__receptor_type )[ "FBK_SPIKES_0" ] = 1,
    ( *__receptor_type )[ "FBK_SPIKES_1" ] = 2,
    ( *__receptor_type )[ "FBK_SPIKES_2" ] = 3,
    ( *__receptor_type )[ "FBK_SPIKES_3" ] = 4,
    ( *__receptor_type )[ "FBK_SPIKES_4" ] = 5,
    ( *__receptor_type )[ "FBK_SPIKES_5" ] = 6,
    ( *__receptor_type )[ "FBK_SPIKES_6" ] = 7,
    ( *__receptor_type )[ "FBK_SPIKES_7" ] = 8,
    ( *__receptor_type )[ "FBK_SPIKES_8" ] = 9,
    ( *__receptor_type )[ "FBK_SPIKES_9" ] = 10,
    ( *__receptor_type )[ "FBK_SPIKES_10" ] = 11,
    ( *__receptor_type )[ "FBK_SPIKES_11" ] = 12,
    ( *__receptor_type )[ "FBK_SPIKES_12" ] = 13,
    ( *__receptor_type )[ "FBK_SPIKES_13" ] = 14,
    ( *__receptor_type )[ "FBK_SPIKES_14" ] = 15,
    ( *__receptor_type )[ "FBK_SPIKES_15" ] = 16,
    ( *__receptor_type )[ "FBK_SPIKES_16" ] = 17,
    ( *__receptor_type )[ "FBK_SPIKES_17" ] = 18,
    ( *__receptor_type )[ "FBK_SPIKES_18" ] = 19,
    ( *__receptor_type )[ "FBK_SPIKES_19" ] = 20,
    ( *__receptor_type )[ "FBK_SPIKES_20" ] = 21,
    ( *__receptor_type )[ "FBK_SPIKES_21" ] = 22,
    ( *__receptor_type )[ "FBK_SPIKES_22" ] = 23,
    ( *__receptor_type )[ "FBK_SPIKES_23" ] = 24,
    ( *__receptor_type )[ "FBK_SPIKES_24" ] = 25,
    ( *__receptor_type )[ "FBK_SPIKES_25" ] = 26,
    ( *__receptor_type )[ "FBK_SPIKES_26" ] = 27,
    ( *__receptor_type )[ "FBK_SPIKES_27" ] = 28,
    ( *__receptor_type )[ "FBK_SPIKES_28" ] = 29,
    ( *__receptor_type )[ "FBK_SPIKES_29" ] = 30,
    ( *__receptor_type )[ "FBK_SPIKES_30" ] = 31,
    ( *__receptor_type )[ "FBK_SPIKES_31" ] = 32,
    ( *__receptor_type )[ "FBK_SPIKES_32" ] = 33,
    ( *__receptor_type )[ "FBK_SPIKES_33" ] = 34,
    ( *__receptor_type )[ "FBK_SPIKES_34" ] = 35,
    ( *__receptor_type )[ "FBK_SPIKES_35" ] = 36,
    ( *__receptor_type )[ "FBK_SPIKES_36" ] = 37,
    ( *__receptor_type )[ "FBK_SPIKES_37" ] = 38,
    ( *__receptor_type )[ "FBK_SPIKES_38" ] = 39,
    ( *__receptor_type )[ "FBK_SPIKES_39" ] = 40,
    ( *__receptor_type )[ "FBK_SPIKES_40" ] = 41,
    ( *__receptor_type )[ "FBK_SPIKES_41" ] = 42,
    ( *__receptor_type )[ "FBK_SPIKES_42" ] = 43,
    ( *__receptor_type )[ "FBK_SPIKES_43" ] = 44,
    ( *__receptor_type )[ "FBK_SPIKES_44" ] = 45,
    ( *__receptor_type )[ "FBK_SPIKES_45" ] = 46,
    ( *__receptor_type )[ "FBK_SPIKES_46" ] = 47,
    ( *__receptor_type )[ "FBK_SPIKES_47" ] = 48,
    ( *__receptor_type )[ "FBK_SPIKES_48" ] = 49,
    ( *__receptor_type )[ "FBK_SPIKES_49" ] = 50,
    ( *__receptor_type )[ "PRED_SPIKES_0" ] = 51,
    ( *__receptor_type )[ "PRED_SPIKES_1" ] = 52,
    ( *__receptor_type )[ "PRED_SPIKES_2" ] = 53,
    ( *__receptor_type )[ "PRED_SPIKES_3" ] = 54,
    ( *__receptor_type )[ "PRED_SPIKES_4" ] = 55,
    ( *__receptor_type )[ "PRED_SPIKES_5" ] = 56,
    ( *__receptor_type )[ "PRED_SPIKES_6" ] = 57,
    ( *__receptor_type )[ "PRED_SPIKES_7" ] = 58,
    ( *__receptor_type )[ "PRED_SPIKES_8" ] = 59,
    ( *__receptor_type )[ "PRED_SPIKES_9" ] = 60,
    ( *__receptor_type )[ "PRED_SPIKES_10" ] = 61,
    ( *__receptor_type )[ "PRED_SPIKES_11" ] = 62,
    ( *__receptor_type )[ "PRED_SPIKES_12" ] = 63,
    ( *__receptor_type )[ "PRED_SPIKES_13" ] = 64,
    ( *__receptor_type )[ "PRED_SPIKES_14" ] = 65,
    ( *__receptor_type )[ "PRED_SPIKES_15" ] = 66,
    ( *__receptor_type )[ "PRED_SPIKES_16" ] = 67,
    ( *__receptor_type )[ "PRED_SPIKES_17" ] = 68,
    ( *__receptor_type )[ "PRED_SPIKES_18" ] = 69,
    ( *__receptor_type )[ "PRED_SPIKES_19" ] = 70,
    ( *__receptor_type )[ "PRED_SPIKES_20" ] = 71,
    ( *__receptor_type )[ "PRED_SPIKES_21" ] = 72,
    ( *__receptor_type )[ "PRED_SPIKES_22" ] = 73,
    ( *__receptor_type )[ "PRED_SPIKES_23" ] = 74,
    ( *__receptor_type )[ "PRED_SPIKES_24" ] = 75,
    ( *__receptor_type )[ "PRED_SPIKES_25" ] = 76,
    ( *__receptor_type )[ "PRED_SPIKES_26" ] = 77,
    ( *__receptor_type )[ "PRED_SPIKES_27" ] = 78,
    ( *__receptor_type )[ "PRED_SPIKES_28" ] = 79,
    ( *__receptor_type )[ "PRED_SPIKES_29" ] = 80,
    ( *__receptor_type )[ "PRED_SPIKES_30" ] = 81,
    ( *__receptor_type )[ "PRED_SPIKES_31" ] = 82,
    ( *__receptor_type )[ "PRED_SPIKES_32" ] = 83,
    ( *__receptor_type )[ "PRED_SPIKES_33" ] = 84,
    ( *__receptor_type )[ "PRED_SPIKES_34" ] = 85,
    ( *__receptor_type )[ "PRED_SPIKES_35" ] = 86,
    ( *__receptor_type )[ "PRED_SPIKES_36" ] = 87,
    ( *__receptor_type )[ "PRED_SPIKES_37" ] = 88,
    ( *__receptor_type )[ "PRED_SPIKES_38" ] = 89,
    ( *__receptor_type )[ "PRED_SPIKES_39" ] = 90,
    ( *__receptor_type )[ "PRED_SPIKES_40" ] = 91,
    ( *__receptor_type )[ "PRED_SPIKES_41" ] = 92,
    ( *__receptor_type )[ "PRED_SPIKES_42" ] = 93,
    ( *__receptor_type )[ "PRED_SPIKES_43" ] = 94,
    ( *__receptor_type )[ "PRED_SPIKES_44" ] = 95,
    ( *__receptor_type )[ "PRED_SPIKES_45" ] = 96,
    ( *__receptor_type )[ "PRED_SPIKES_46" ] = 97,
    ( *__receptor_type )[ "PRED_SPIKES_47" ] = 98,
    ( *__receptor_type )[ "PRED_SPIKES_48" ] = 99,
    ( *__receptor_type )[ "PRED_SPIKES_49" ] = 100,
    ( *__d )[ "receptor_types" ] = __receptor_type;

  (*__d)[nest::names::recordables] = recordablesMap_.get_list();
}

inline void state_neuron_nestml::set_status(const DictionaryDatum &__d)
{
  // parameters
  double tmp_kp = get_kp();
  nest::updateValueParam<double>(__d, nest::state_neuron_nestml_names::_kp, tmp_kp, this);
  // Resize vectors
  if (tmp_kp != get_kp())
  {
  }
  bool tmp_pos = get_pos();
  nest::updateValueParam<bool>(__d, nest::state_neuron_nestml_names::_pos, tmp_pos, this);
  // Resize vectors
  if (tmp_pos != get_pos())
  {
  }
  double tmp_base_rate = get_base_rate();
  nest::updateValueParam<double>(__d, nest::state_neuron_nestml_names::_base_rate, tmp_base_rate, this);
  // Resize vectors
  if (tmp_base_rate != get_base_rate())
  {
  }
  double tmp_buffer_size = get_buffer_size();
  nest::updateValueParam<double>(__d, nest::state_neuron_nestml_names::_buffer_size, tmp_buffer_size, this);
  // Resize vectors
  if (tmp_buffer_size != get_buffer_size())
  {
  }
  long tmp_simulation_steps = get_simulation_steps();
  nest::updateValueParam<long>(__d, nest::state_neuron_nestml_names::_simulation_steps, tmp_simulation_steps, this);
  // Resize vectors
  if (tmp_simulation_steps != get_simulation_steps())
  {
  }
  long tmp_N_fbk = get_N_fbk();
  nest::updateValueParam<long>(__d, nest::state_neuron_nestml_names::_N_fbk, tmp_N_fbk, this);
  // Resize vectors
  if (tmp_N_fbk != get_N_fbk())
  {
    std::vector< double >  _tmp_current_fbk_input = get_current_fbk_input();
    _tmp_current_fbk_input.resize(tmp_N_fbk, 0.);
    set_current_fbk_input(_tmp_current_fbk_input);
    std::vector< double >  _tmp_fbk_counts = get_fbk_counts();
    _tmp_fbk_counts.resize(tmp_N_fbk, 0.);
    set_fbk_counts(_tmp_fbk_counts);
  }
  long tmp_N_pred = get_N_pred();
  nest::updateValueParam<long>(__d, nest::state_neuron_nestml_names::_N_pred, tmp_N_pred, this);
  // Resize vectors
  if (tmp_N_pred != get_N_pred())
  {
    std::vector< double >  _tmp_current_pred_input = get_current_pred_input();
    _tmp_current_pred_input.resize(tmp_N_pred, 0.);
    set_current_pred_input(_tmp_current_pred_input);
    std::vector< double >  _tmp_pred_counts = get_pred_counts();
    _tmp_pred_counts.resize(tmp_N_pred, 0.);
    set_pred_counts(_tmp_pred_counts);
  }
  long tmp_fbk_bf_size = get_fbk_bf_size();
  nest::updateValueParam<long>(__d, nest::state_neuron_nestml_names::_fbk_bf_size, tmp_fbk_bf_size, this);
  // Resize vectors
  if (tmp_fbk_bf_size != get_fbk_bf_size())
  {
    std::vector< double >  _tmp_fbk_buffer = get_fbk_buffer();
    _tmp_fbk_buffer.resize(tmp_fbk_bf_size, 0.);
    set_fbk_buffer(_tmp_fbk_buffer);
  }
  long tmp_pred_bf_size = get_pred_bf_size();
  nest::updateValueParam<long>(__d, nest::state_neuron_nestml_names::_pred_bf_size, tmp_pred_bf_size, this);
  // Resize vectors
  if (tmp_pred_bf_size != get_pred_bf_size())
  {
    std::vector< double >  _tmp_pred_buffer = get_pred_buffer();
    _tmp_pred_buffer.resize(tmp_pred_bf_size, 0.);
    set_pred_buffer(_tmp_pred_buffer);
  }
  double tmp_time_wait = get_time_wait();
  nest::updateValueParam<double>(__d, nest::state_neuron_nestml_names::_time_wait, tmp_time_wait, this);
  // Resize vectors
  if (tmp_time_wait != get_time_wait())
  {
  }
  double tmp_time_trial = get_time_trial();
  nest::updateValueParam<double>(__d, nest::state_neuron_nestml_names::_time_trial, tmp_time_trial, this);
  // Resize vectors
  if (tmp_time_trial != get_time_trial())
  {
  }

  // initial values for state variables in ODE or kernel
  double tmp_in_rate = get_in_rate();
  nest::updateValueParam<double>(__d, nest::state_neuron_nestml_names::_in_rate, tmp_in_rate, this);
  double tmp_out_rate = get_out_rate();
  nest::updateValueParam<double>(__d, nest::state_neuron_nestml_names::_out_rate, tmp_out_rate, this);
  long tmp_spike_count_out = get_spike_count_out();
  nest::updateValueParam<long>(__d, nest::state_neuron_nestml_names::_spike_count_out, tmp_spike_count_out, this);
  std::vector< double >  tmp_current_fbk_input = get_current_fbk_input();
  updateValue<std::vector< double > >(__d, nest::state_neuron_nestml_names::_current_fbk_input, tmp_current_fbk_input);
   
  // Check if the new vector size matches its original size
  if ( tmp_current_fbk_input.size() != tmp_N_fbk )
  {
    std::stringstream msg;
    msg << "The vector \"current_fbk_input\" does not match its size: " << tmp_N_fbk;
    throw nest::BadProperty(msg.str());
  }
  std::vector< double >  tmp_current_pred_input = get_current_pred_input();
  updateValue<std::vector< double > >(__d, nest::state_neuron_nestml_names::_current_pred_input, tmp_current_pred_input);
   
  // Check if the new vector size matches its original size
  if ( tmp_current_pred_input.size() != tmp_N_pred )
  {
    std::stringstream msg;
    msg << "The vector \"current_pred_input\" does not match its size: " << tmp_N_pred;
    throw nest::BadProperty(msg.str());
  }
  std::vector< double >  tmp_fbk_buffer = get_fbk_buffer();
  updateValue<std::vector< double > >(__d, nest::state_neuron_nestml_names::_fbk_buffer, tmp_fbk_buffer);
   
  // Check if the new vector size matches its original size
  if ( tmp_fbk_buffer.size() != tmp_fbk_bf_size )
  {
    std::stringstream msg;
    msg << "The vector \"fbk_buffer\" does not match its size: " << tmp_fbk_bf_size;
    throw nest::BadProperty(msg.str());
  }
  std::vector< double >  tmp_pred_buffer = get_pred_buffer();
  updateValue<std::vector< double > >(__d, nest::state_neuron_nestml_names::_pred_buffer, tmp_pred_buffer);
   
  // Check if the new vector size matches its original size
  if ( tmp_pred_buffer.size() != tmp_pred_bf_size )
  {
    std::stringstream msg;
    msg << "The vector \"pred_buffer\" does not match its size: " << tmp_pred_bf_size;
    throw nest::BadProperty(msg.str());
  }
  std::vector< double >  tmp_fbk_counts = get_fbk_counts();
  updateValue<std::vector< double > >(__d, nest::state_neuron_nestml_names::_fbk_counts, tmp_fbk_counts);
   
  // Check if the new vector size matches its original size
  if ( tmp_fbk_counts.size() != tmp_N_fbk )
  {
    std::stringstream msg;
    msg << "The vector \"fbk_counts\" does not match its size: " << tmp_N_fbk;
    throw nest::BadProperty(msg.str());
  }
  std::vector< double >  tmp_pred_counts = get_pred_counts();
  updateValue<std::vector< double > >(__d, nest::state_neuron_nestml_names::_pred_counts, tmp_pred_counts);
   
  // Check if the new vector size matches its original size
  if ( tmp_pred_counts.size() != tmp_N_pred )
  {
    std::stringstream msg;
    msg << "The vector \"pred_counts\" does not match its size: " << tmp_N_pred;
    throw nest::BadProperty(msg.str());
  }
  long tmp_tick = get_tick();
  nest::updateValueParam<long>(__d, nest::state_neuron_nestml_names::_tick, tmp_tick, this);
  long tmp_position_count = get_position_count();
  nest::updateValueParam<long>(__d, nest::state_neuron_nestml_names::_position_count, tmp_position_count, this);
  double tmp_mean_fbk = get_mean_fbk();
  nest::updateValueParam<double>(__d, nest::state_neuron_nestml_names::_mean_fbk, tmp_mean_fbk, this);
  double tmp_mean_pred = get_mean_pred();
  nest::updateValueParam<double>(__d, nest::state_neuron_nestml_names::_mean_pred, tmp_mean_pred, this);
  double tmp_var_fbk = get_var_fbk();
  nest::updateValueParam<double>(__d, nest::state_neuron_nestml_names::_var_fbk, tmp_var_fbk, this);
  double tmp_var_pred = get_var_pred();
  nest::updateValueParam<double>(__d, nest::state_neuron_nestml_names::_var_pred, tmp_var_pred, this);
  double tmp_CV_fbk = get_CV_fbk();
  nest::updateValueParam<double>(__d, nest::state_neuron_nestml_names::_CV_fbk, tmp_CV_fbk, this);
  double tmp_CV_pred = get_CV_pred();
  nest::updateValueParam<double>(__d, nest::state_neuron_nestml_names::_CV_pred, tmp_CV_pred, this);
  double tmp_total_CV = get_total_CV();
  nest::updateValueParam<double>(__d, nest::state_neuron_nestml_names::_total_CV, tmp_total_CV, this);
  double tmp_lambda_poisson = get_lambda_poisson();
  nest::updateValueParam<double>(__d, nest::state_neuron_nestml_names::_lambda_poisson, tmp_lambda_poisson, this);

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
  set_N_fbk(tmp_N_fbk);
  set_N_pred(tmp_N_pred);
  set_fbk_bf_size(tmp_fbk_bf_size);
  set_pred_bf_size(tmp_pred_bf_size);
  set_time_wait(tmp_time_wait);
  set_time_trial(tmp_time_trial);
  set_in_rate(tmp_in_rate);
  set_out_rate(tmp_out_rate);
  set_spike_count_out(tmp_spike_count_out);
  set_current_fbk_input(tmp_current_fbk_input);
  set_current_pred_input(tmp_current_pred_input);
  set_fbk_buffer(tmp_fbk_buffer);
  set_pred_buffer(tmp_pred_buffer);
  set_fbk_counts(tmp_fbk_counts);
  set_pred_counts(tmp_pred_counts);
  set_tick(tmp_tick);
  set_position_count(tmp_position_count);
  set_mean_fbk(tmp_mean_fbk);
  set_mean_pred(tmp_mean_pred);
  set_var_fbk(tmp_var_fbk);
  set_var_pred(tmp_var_pred);
  set_CV_fbk(tmp_CV_fbk);
  set_CV_pred(tmp_CV_pred);
  set_total_CV(tmp_total_CV);
  set_lambda_poisson(tmp_lambda_poisson);





  // recompute internal variables in case they are dependent on parameters or state that might have been updated in this call to set_status()
  recompute_internal_variables();
};



#endif /* #ifndef STATE_NEURON_NESTML */
