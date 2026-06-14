
/**
 *  tracking_neuron_planner_nestml.h
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
 *  Generated from NESTML at time: 2024-10-11 08:49:25.216976
**/
#ifndef TRACKING_NEURON_PLANNER_NESTML
#define TRACKING_NEURON_PLANNER_NESTML

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
namespace tracking_neuron_planner_nestml_names
{
    const Name _out_rate( "out_rate" );
    const Name _lambda( "lambda" );
    const Name _spike_count_out( "spike_count_out" );
    const Name _current_step( "current_step" );
    const Name _curr_traj( "curr_traj" );
    const Name _kp( "kp" );
    const Name _pos( "pos" );
    const Name _base_rate( "base_rate" );
    const Name _simulation_steps( "simulation_steps" );
    const Name _traj( "traj" );
}
}




#include "nest_time.h"
  typedef size_t nest_port_t;
  typedef size_t nest_rport_t;

/* BeginDocumentation
  Name: tracking_neuron_planner_nestml

  Description:

    

  Parameters:
  The following parameters can be set in the status dictionary.
kp [real]  Gain parameter
pos [boolean]  Sensitivity of neuron to positive or negative values
base_rate [real]  Base firing rate
simulation_steps [integer]  Simulation steps -> computed before calling the neuron model as the length of the time vector
traj [real]  Trajectory vector (defined by function GetDesiredTrajectory)


  Dynamic state variables:


  Sends: nest::SpikeEvent

  Receives: Spike,  DataLoggingRequest
*/

// Register the neuron model
void register_tracking_neuron_planner_nestml( const std::string& name );

class tracking_neuron_planner_nestml : public nest::StructuralPlasticityNode
{
public:
  /**
   * The constructor is only used to create the model prototype in the model manager.
  **/
  tracking_neuron_planner_nestml();

  /**
   * The copy constructor is used to create model copies and instances of the model.
   * @node The copy constructor needs to initialize the parameters and the state.
   *       Initialization of buffers and interal variables is deferred to
   *       @c init_buffers_() and @c pre_run_hook() (or calibrate() in NEST 3.3 and older).
  **/
  tracking_neuron_planner_nestml(const tracking_neuron_planner_nestml &);

  /**
   * Destructor.
  **/
  ~tracking_neuron_planner_nestml() override;

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

  inline double get_out_rate() const
  {
    return S_.out_rate;
  }

  inline void set_out_rate(const double __v)
  {
    S_.out_rate = __v;
  }

  inline double get_lambda() const
  {
    return S_.lambda;
  }

  inline void set_lambda(const double __v)
  {
    S_.lambda = __v;
  }

  inline long get_spike_count_out() const
  {
    return S_.spike_count_out;
  }

  inline void set_spike_count_out(const long __v)
  {
    S_.spike_count_out = __v;
  }

  inline long get_current_step() const
  {
    return S_.current_step;
  }

  inline void set_current_step(const long __v)
  {
    S_.current_step = __v;
  }

  inline double get_curr_traj() const
  {
    return S_.curr_traj;
  }

  inline void set_curr_traj(const double __v)
  {
    S_.curr_traj = __v;
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

  inline long get_simulation_steps() const
  {
    return P_.simulation_steps;
  }

  inline void set_simulation_steps(const long __v)
  {
    P_.simulation_steps = __v;
  }

  inline std::vector< double >  get_traj() const
  {
    return P_.traj;
  }

  inline void set_traj(const std::vector< double >  __v)
  {
    P_.traj = __v;
  }


  // -------------------------------------------------------------------------
  //   Getters/setters for internals
  // -------------------------------------------------------------------------

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
    INH_SPIKES = 0,
    EXC_SPIKES = 1,
    SPIKES = 2,
    MAX_SPIKE_RECEPTOR = 3
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
  friend class nest::RecordablesMap<tracking_neuron_planner_nestml>;
  friend class nest::UniversalDataLogger<tracking_neuron_planner_nestml>;

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
    //!  Gain parameter
    double kp;
    //!  Sensitivity of neuron to positive or negative values
    bool pos;
    //!  Base firing rate
    double base_rate;
    //!  Simulation steps -> computed before calling the neuron model as the length of the time vector
    long simulation_steps;
    //!  Trajectory vector (defined by function GetDesiredTrajectory)
    std::vector< double >  traj;

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
    double out_rate;
    double lambda;
    long spike_count_out;
    long current_step;
    double curr_traj;

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
    Buffers_(tracking_neuron_planner_nestml &);
    Buffers_(const Buffers_ &, tracking_neuron_planner_nestml &);

    /**
     * Logger for all analog data
    **/
    nest::UniversalDataLogger<tracking_neuron_planner_nestml> logger_;

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
  static nest::RecordablesMap<tracking_neuron_planner_nestml> recordablesMap_;
  nest::normal_distribution normal_dev_; //!< random deviate generator
  nest::poisson_distribution poisson_dev_; //!< random deviate generator

}; /* neuron tracking_neuron_planner_nestml */

inline nest_port_t tracking_neuron_planner_nestml::send_test_event(nest::Node& target, nest_rport_t receptor_type, nest::synindex, bool)
{
  // You should usually not change the code in this function.
  // It confirms that the target of connection @c c accepts @c nest::SpikeEvent on
  // the given @c receptor_type.
  nest::SpikeEvent e;
  e.set_sender(*this);
  return target.handles_test_event(e, receptor_type);
}

inline nest_port_t tracking_neuron_planner_nestml::handles_test_event(nest::SpikeEvent&, nest_port_t receptor_type)
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

inline nest_port_t tracking_neuron_planner_nestml::handles_test_event(nest::DataLoggingRequest& dlr, nest_port_t receptor_type)
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

inline void tracking_neuron_planner_nestml::get_status(DictionaryDatum &__d) const
{
  // parameters
  def< double >(__d, nest::tracking_neuron_planner_nestml_names::_kp, get_kp());
  def< bool >(__d, nest::tracking_neuron_planner_nestml_names::_pos, get_pos());
  def< double >(__d, nest::tracking_neuron_planner_nestml_names::_base_rate, get_base_rate());
  def< long >(__d, nest::tracking_neuron_planner_nestml_names::_simulation_steps, get_simulation_steps());
  def< std::vector< double >  >(__d, nest::tracking_neuron_planner_nestml_names::_traj, get_traj());

  // initial values for state variables in ODE or kernel
  def< double >(__d, nest::tracking_neuron_planner_nestml_names::_out_rate, get_out_rate());
  def< double >(__d, nest::tracking_neuron_planner_nestml_names::_lambda, get_lambda());
  def< long >(__d, nest::tracking_neuron_planner_nestml_names::_spike_count_out, get_spike_count_out());
  def< long >(__d, nest::tracking_neuron_planner_nestml_names::_current_step, get_current_step());
  def< double >(__d, nest::tracking_neuron_planner_nestml_names::_curr_traj, get_curr_traj());

  StructuralPlasticityNode::get_status( __d );

  (*__d)[nest::names::recordables] = recordablesMap_.get_list();
}

inline void tracking_neuron_planner_nestml::set_status(const DictionaryDatum &__d)
{
  // parameters
  double tmp_kp = get_kp();
  nest::updateValueParam<double>(__d, nest::tracking_neuron_planner_nestml_names::_kp, tmp_kp, this);
  // Resize vectors
  if (tmp_kp != get_kp())
  {
  }
  bool tmp_pos = get_pos();
  nest::updateValueParam<bool>(__d, nest::tracking_neuron_planner_nestml_names::_pos, tmp_pos, this);
  // Resize vectors
  if (tmp_pos != get_pos())
  {
  }
  double tmp_base_rate = get_base_rate();
  nest::updateValueParam<double>(__d, nest::tracking_neuron_planner_nestml_names::_base_rate, tmp_base_rate, this);
  // Resize vectors
  if (tmp_base_rate != get_base_rate())
  {
  }
  long tmp_simulation_steps = get_simulation_steps();
  nest::updateValueParam<long>(__d, nest::tracking_neuron_planner_nestml_names::_simulation_steps, tmp_simulation_steps, this);
  // Resize vectors
  if (tmp_simulation_steps != get_simulation_steps())
  {
    std::vector< double >  _tmp_traj = get_traj();
    _tmp_traj.resize(tmp_simulation_steps, 0.);
    set_traj(_tmp_traj);
  }
  std::vector< double >  tmp_traj = get_traj();
  updateValue<std::vector< double > >(__d, nest::tracking_neuron_planner_nestml_names::_traj, tmp_traj);
  // Resize vectors
  if (tmp_traj != get_traj())
  {
  }
   
  // Check if the new vector size matches its original size
  if ( tmp_traj.size() != tmp_simulation_steps )
  {
    std::stringstream msg;
    msg << "The vector \"traj\" does not match its size: " << tmp_simulation_steps;
    throw nest::BadProperty(msg.str());
  }

  // initial values for state variables in ODE or kernel
  double tmp_out_rate = get_out_rate();
  nest::updateValueParam<double>(__d, nest::tracking_neuron_planner_nestml_names::_out_rate, tmp_out_rate, this);
  double tmp_lambda = get_lambda();
  nest::updateValueParam<double>(__d, nest::tracking_neuron_planner_nestml_names::_lambda, tmp_lambda, this);
  long tmp_spike_count_out = get_spike_count_out();
  nest::updateValueParam<long>(__d, nest::tracking_neuron_planner_nestml_names::_spike_count_out, tmp_spike_count_out, this);
  long tmp_current_step = get_current_step();
  nest::updateValueParam<long>(__d, nest::tracking_neuron_planner_nestml_names::_current_step, tmp_current_step, this);
  double tmp_curr_traj = get_curr_traj();
  nest::updateValueParam<double>(__d, nest::tracking_neuron_planner_nestml_names::_curr_traj, tmp_curr_traj, this);

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  StructuralPlasticityNode::set_status(__d);

  // if we get here, temporaries contain consistent set of properties
  set_kp(tmp_kp);
  set_pos(tmp_pos);
  set_base_rate(tmp_base_rate);
  set_simulation_steps(tmp_simulation_steps);
  set_traj(tmp_traj);
  set_out_rate(tmp_out_rate);
  set_lambda(tmp_lambda);
  set_spike_count_out(tmp_spike_count_out);
  set_current_step(tmp_current_step);
  set_curr_traj(tmp_curr_traj);





  // recompute internal variables in case they are dependent on parameters or state that might have been updated in this call to set_status()
  recompute_internal_variables();
};



#endif /* #ifndef TRACKING_NEURON_PLANNER_NESTML */
