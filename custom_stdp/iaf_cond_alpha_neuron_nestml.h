
/**
 *  iaf_cond_alpha_neuron_nestml.h
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
#ifndef IAF_COND_ALPHA_NEURON_NESTML
#define IAF_COND_ALPHA_NEURON_NESTML

#ifndef HAVE_LIBLTDL
#error "NEST was compiled without support for dynamic loading. Please install libltdl and recompile NEST."
#endif

// C++ includes:
#include <cmath>

#include "config.h"

#ifndef HAVE_GSL
#error "The GSL library is required for the Runge-Kutta solver."
#endif

// External includes:
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// Includes from nestkernel:
#include "extended_post_history_archiving_node.h"
#include "connection.h"
#include "dict_util.h"
#include "event.h"
#include "nest_types.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"

// Includes from sli:
#include "dictdatum.h"

// uncomment the next line to enable printing of detailed debug information
// #define DEBUG

namespace nest
{
namespace iaf_cond_alpha_neuron_nestml_names
{
    const Name _V_m( "V_m" );
    const Name _refr_t( "refr_t" );
    const Name _g_inh__X__inh_spikes__DOT__weight( "g_inh__X__inh_spikes__DOT__weight" );
    const Name _g_inh__X__inh_spikes__DOT__weight__d( "g_inh__X__inh_spikes__DOT__weight__d" );
    const Name _g_exc__X__exc_spikes__DOT__weight( "g_exc__X__exc_spikes__DOT__weight" );
    const Name _g_exc__X__exc_spikes__DOT__weight__d( "g_exc__X__exc_spikes__DOT__weight__d" );
    const Name _C_m( "C_m" );
    const Name _g_L( "g_L" );
    const Name _E_L( "E_L" );
    const Name _refr_T( "refr_T" );
    const Name _V_th( "V_th" );
    const Name _V_reset( "V_reset" );
    const Name _E_exc( "E_exc" );
    const Name _use_offset_trick( "use_offset_trick" );
    const Name _E_inh( "E_inh" );
    const Name _tau_syn_exc( "tau_syn_exc" );
    const Name _tau_syn_inh( "tau_syn_inh" );
    const Name _I_e( "I_e" );

    const Name gsl_abs_error_tol("gsl_abs_error_tol");
    const Name gsl_rel_error_tol("gsl_rel_error_tol");
}
}



/**
 * Function computing right-hand side of ODE for GSL solver.
 * @note Must be declared here so we can befriend it in class.
 * @note Must have C-linkage for passing to GSL. Internally, it is
 *       a first-class C++ function, but cannot be a member function
 *       because of the C-linkage.
 * @note No point in declaring it inline, since it is called
 *       through a function pointer.
 * @param void* Pointer to model neuron instance.
 *
 * Integrate the variables: refr_t
**/
extern "C" inline int iaf_cond_alpha_neuron_nestml_dynamics_refr_t( double, const double ode_state[], double f[], void* pnode );
/**
 * Function computing right-hand side of ODE for GSL solver.
 * @note Must be declared here so we can befriend it in class.
 * @note Must have C-linkage for passing to GSL. Internally, it is
 *       a first-class C++ function, but cannot be a member function
 *       because of the C-linkage.
 * @note No point in declaring it inline, since it is called
 *       through a function pointer.
 * @param void* Pointer to model neuron instance.
 *
 * Integrate the variables: V_m
**/
extern "C" inline int iaf_cond_alpha_neuron_nestml_dynamics_V_m( double, const double ode_state[], double f[], void* pnode );
/**
 * Function computing right-hand side of ODE for GSL solver.
 * @note Must be declared here so we can befriend it in class.
 * @note Must have C-linkage for passing to GSL. Internally, it is
 *       a first-class C++ function, but cannot be a member function
 *       because of the C-linkage.
 * @note No point in declaring it inline, since it is called
 *       through a function pointer.
 * @param void* Pointer to model neuron instance.
 *
 * Integrate the variables:
**/
extern "C" inline int iaf_cond_alpha_neuron_nestml_dynamics( double, const double ode_state[], double f[], void* pnode );

#include "nest_time.h"

typedef size_t nest_port_t;
typedef size_t nest_rport_t;

// Register the neuron model
void register_iaf_cond_alpha_neuron_nestml( const std::string& name );

class iaf_cond_alpha_neuron_nestml : public nest::ExtendedPostHistoryArchivingNode
{
public:
  /**
   * The constructor is only used to create the model prototype in the model manager.
  **/
  iaf_cond_alpha_neuron_nestml();

  /**
   * The copy constructor is used to create model copies and instances of the model.
   * @node The copy constructor needs to initialize the parameters and the state.
   *       Initialization of buffers and interal variables is deferred to
   *       @c init_buffers_() and @c pre_run_hook() (or calibrate() in NEST 3.3 and older).
  **/
  iaf_cond_alpha_neuron_nestml(const iaf_cond_alpha_neuron_nestml &);

  /**
   * Destructor.
  **/
  ~iaf_cond_alpha_neuron_nestml() override;

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
  void handle(nest::CurrentEvent &) override;      //! accept input current

  void handle(nest::DataLoggingRequest &) override;//! allow recording with multimeter
  nest_port_t handles_test_event(nest::SpikeEvent&, nest_port_t) override;
  nest_port_t handles_test_event(nest::CurrentEvent&, nest_port_t) override;
  nest_port_t handles_test_event(nest::DataLoggingRequest&, nest_port_t) override;

  // -------------------------------------------------------------------------
  //   Functions for getting/setting parameters and state values.
  // -------------------------------------------------------------------------

  void get_status(DictionaryDatum &) const override;
  void set_status(const DictionaryDatum &) override;


  // -------------------------------------------------------------------------
  //   Getters/setters for state block
  // -------------------------------------------------------------------------

  inline double get_V_m() const
  {
    return S_.ode_state[State_::V_m];
  }

  inline void set_V_m(const double __v)
  {
    S_.ode_state[State_::V_m] = __v;
  }

  inline double get_refr_t() const
  {
    return S_.ode_state[State_::refr_t];
  }

  inline void set_refr_t(const double __v)
  {
    S_.ode_state[State_::refr_t] = __v;
  }

  inline double get_g_inh__X__inh_spikes__DOT__weight() const
  {
    return S_.ode_state[State_::g_inh__X__inh_spikes__DOT__weight];
  }

  inline void set_g_inh__X__inh_spikes__DOT__weight(const double __v)
  {
    S_.ode_state[State_::g_inh__X__inh_spikes__DOT__weight] = __v;
  }

  inline double get_g_inh__X__inh_spikes__DOT__weight__d() const
  {
    return S_.ode_state[State_::g_inh__X__inh_spikes__DOT__weight__d];
  }

  inline void set_g_inh__X__inh_spikes__DOT__weight__d(const double __v)
  {
    S_.ode_state[State_::g_inh__X__inh_spikes__DOT__weight__d] = __v;
  }

  inline double get_g_exc__X__exc_spikes__DOT__weight() const
  {
    return S_.ode_state[State_::g_exc__X__exc_spikes__DOT__weight];
  }

  inline void set_g_exc__X__exc_spikes__DOT__weight(const double __v)
  {
    S_.ode_state[State_::g_exc__X__exc_spikes__DOT__weight] = __v;
  }

  inline double get_g_exc__X__exc_spikes__DOT__weight__d() const
  {
    return S_.ode_state[State_::g_exc__X__exc_spikes__DOT__weight__d];
  }

  inline void set_g_exc__X__exc_spikes__DOT__weight__d(const double __v)
  {
    S_.ode_state[State_::g_exc__X__exc_spikes__DOT__weight__d] = __v;
  }


  // -------------------------------------------------------------------------
  //   Getters/setters for parameters
  // -------------------------------------------------------------------------

  inline double get_C_m() const
  {
    return P_.C_m;
  }

  inline void set_C_m(const double __v)
  {
    P_.C_m = __v;
  }

  inline double get_g_L() const
  {
    return P_.g_L;
  }

  inline void set_g_L(const double __v)
  {
    P_.g_L = __v;
  }

  inline double get_E_L() const
  {
    return P_.E_L;
  }

  inline void set_E_L(const double __v)
  {
    P_.E_L = __v;
  }

  inline double get_refr_T() const
  {
    return P_.refr_T;
  }

  inline void set_refr_T(const double __v)
  {
    P_.refr_T = __v;
  }

  inline double get_V_th() const
  {
    return P_.V_th;
  }

  inline void set_V_th(const double __v)
  {
    P_.V_th = __v;
  }

  inline double get_V_reset() const
  {
    return P_.V_reset;
  }

  inline void set_V_reset(const double __v)
  {
    P_.V_reset = __v;
  }

  inline double get_E_exc() const
  {
    return P_.E_exc;
  }

  inline void set_E_exc(const double __v)
  {
    P_.E_exc = __v;
  }

  inline bool get_use_offset_trick() const
  {
    return P_.use_offset_trick;
  }

  inline void set_use_offset_trick(const bool __v)
  {
    P_.use_offset_trick = __v;
  }

  inline double get_E_inh() const
  {
    return P_.E_inh;
  }

  inline void set_E_inh(const double __v)
  {
    P_.E_inh = __v;
  }

  inline double get_tau_syn_exc() const
  {
    return P_.tau_syn_exc;
  }

  inline void set_tau_syn_exc(const double __v)
  {
    P_.tau_syn_exc = __v;
  }

  inline double get_tau_syn_inh() const
  {
    return P_.tau_syn_inh;
  }

  inline void set_tau_syn_inh(const double __v)
  {
    P_.tau_syn_inh = __v;
  }

  inline double get_I_e() const
  {
    return P_.I_e;
  }

  inline void set_I_e(const double __v)
  {
    P_.I_e = __v;
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
  inline double get___P__refr_t__refr_t() const
  {
    return V_.__P__refr_t__refr_t;
  }

  inline void set___P__refr_t__refr_t(const double __v)
  {
    V_.__P__refr_t__refr_t = __v;
  }
  inline double get___P__g_inh__X__inh_spikes__DOT__weight__g_inh__X__inh_spikes__DOT__weight() const
  {
    return V_.__P__g_inh__X__inh_spikes__DOT__weight__g_inh__X__inh_spikes__DOT__weight;
  }

  inline void set___P__g_inh__X__inh_spikes__DOT__weight__g_inh__X__inh_spikes__DOT__weight(const double __v)
  {
    V_.__P__g_inh__X__inh_spikes__DOT__weight__g_inh__X__inh_spikes__DOT__weight = __v;
  }
  inline double get___P__g_inh__X__inh_spikes__DOT__weight__g_inh__X__inh_spikes__DOT__weight__d() const
  {
    return V_.__P__g_inh__X__inh_spikes__DOT__weight__g_inh__X__inh_spikes__DOT__weight__d;
  }

  inline void set___P__g_inh__X__inh_spikes__DOT__weight__g_inh__X__inh_spikes__DOT__weight__d(const double __v)
  {
    V_.__P__g_inh__X__inh_spikes__DOT__weight__g_inh__X__inh_spikes__DOT__weight__d = __v;
  }
  inline double get___P__g_inh__X__inh_spikes__DOT__weight__d__g_inh__X__inh_spikes__DOT__weight() const
  {
    return V_.__P__g_inh__X__inh_spikes__DOT__weight__d__g_inh__X__inh_spikes__DOT__weight;
  }

  inline void set___P__g_inh__X__inh_spikes__DOT__weight__d__g_inh__X__inh_spikes__DOT__weight(const double __v)
  {
    V_.__P__g_inh__X__inh_spikes__DOT__weight__d__g_inh__X__inh_spikes__DOT__weight = __v;
  }
  inline double get___P__g_inh__X__inh_spikes__DOT__weight__d__g_inh__X__inh_spikes__DOT__weight__d() const
  {
    return V_.__P__g_inh__X__inh_spikes__DOT__weight__d__g_inh__X__inh_spikes__DOT__weight__d;
  }

  inline void set___P__g_inh__X__inh_spikes__DOT__weight__d__g_inh__X__inh_spikes__DOT__weight__d(const double __v)
  {
    V_.__P__g_inh__X__inh_spikes__DOT__weight__d__g_inh__X__inh_spikes__DOT__weight__d = __v;
  }
  inline double get___P__g_exc__X__exc_spikes__DOT__weight__g_exc__X__exc_spikes__DOT__weight() const
  {
    return V_.__P__g_exc__X__exc_spikes__DOT__weight__g_exc__X__exc_spikes__DOT__weight;
  }

  inline void set___P__g_exc__X__exc_spikes__DOT__weight__g_exc__X__exc_spikes__DOT__weight(const double __v)
  {
    V_.__P__g_exc__X__exc_spikes__DOT__weight__g_exc__X__exc_spikes__DOT__weight = __v;
  }
  inline double get___P__g_exc__X__exc_spikes__DOT__weight__g_exc__X__exc_spikes__DOT__weight__d() const
  {
    return V_.__P__g_exc__X__exc_spikes__DOT__weight__g_exc__X__exc_spikes__DOT__weight__d;
  }

  inline void set___P__g_exc__X__exc_spikes__DOT__weight__g_exc__X__exc_spikes__DOT__weight__d(const double __v)
  {
    V_.__P__g_exc__X__exc_spikes__DOT__weight__g_exc__X__exc_spikes__DOT__weight__d = __v;
  }
  inline double get___P__g_exc__X__exc_spikes__DOT__weight__d__g_exc__X__exc_spikes__DOT__weight() const
  {
    return V_.__P__g_exc__X__exc_spikes__DOT__weight__d__g_exc__X__exc_spikes__DOT__weight;
  }

  inline void set___P__g_exc__X__exc_spikes__DOT__weight__d__g_exc__X__exc_spikes__DOT__weight(const double __v)
  {
    V_.__P__g_exc__X__exc_spikes__DOT__weight__d__g_exc__X__exc_spikes__DOT__weight = __v;
  }
  inline double get___P__g_exc__X__exc_spikes__DOT__weight__d__g_exc__X__exc_spikes__DOT__weight__d() const
  {
    return V_.__P__g_exc__X__exc_spikes__DOT__weight__d__g_exc__X__exc_spikes__DOT__weight__d;
  }

  inline void set___P__g_exc__X__exc_spikes__DOT__weight__d__g_exc__X__exc_spikes__DOT__weight__d(const double __v)
  {
    V_.__P__g_exc__X__exc_spikes__DOT__weight__d__g_exc__X__exc_spikes__DOT__weight__d = __v;
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
  const nest_port_t MAX_SPIKE_RECEPTOR = 2;



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
  friend class nest::RecordablesMap<iaf_cond_alpha_neuron_nestml>;
  friend class nest::UniversalDataLogger<iaf_cond_alpha_neuron_nestml>;

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
    //!  Membrane capacitance
    double C_m;
    //!  Leak conductance
    double g_L;
    //!  Leak reversal potential (aka resting potential)
    double E_L;
    //!  Duration of refractory period
    double refr_T;
    //!  Spike threshold potential
    double V_th;
    //!  Reset potential
    double V_reset;
    //!  Excitatory reversal potential
    double E_exc;
    bool use_offset_trick;
    //!  Inhibitory reversal potential
    double E_inh;
    //!  Synaptic time constant of excitatory synapse
    double tau_syn_exc;
    //!  Synaptic time constant of inhibitory synapse
    double tau_syn_inh;
    //!  constant external input current
    double I_e;

    double __gsl_abs_error_tol;
    double __gsl_rel_error_tol;

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

    // non-ODE state variables
    //! Symbolic indices to the elements of the state vector y
    enum StateVecElems
    {
      V_m,
      refr_t,
      g_inh__X__inh_spikes__DOT__weight,
      g_inh__X__inh_spikes__DOT__weight__d,
      g_exc__X__exc_spikes__DOT__weight,
      g_exc__X__exc_spikes__DOT__weight__d,
      // moved state variables from synapse (numeric)
      // moved state variables from synapse (analytic)
      // final entry to easily get the vector size
      STATE_VEC_SIZE
    };

    //! state vector, must be C-array for GSL solver
    double ode_state[STATE_VEC_SIZE];

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
    double __P__refr_t__refr_t;
    double __P__g_inh__X__inh_spikes__DOT__weight__g_inh__X__inh_spikes__DOT__weight;
    double __P__g_inh__X__inh_spikes__DOT__weight__g_inh__X__inh_spikes__DOT__weight__d;
    double __P__g_inh__X__inh_spikes__DOT__weight__d__g_inh__X__inh_spikes__DOT__weight;
    double __P__g_inh__X__inh_spikes__DOT__weight__d__g_inh__X__inh_spikes__DOT__weight__d;
    double __P__g_exc__X__exc_spikes__DOT__weight__g_exc__X__exc_spikes__DOT__weight;
    double __P__g_exc__X__exc_spikes__DOT__weight__g_exc__X__exc_spikes__DOT__weight__d;
    double __P__g_exc__X__exc_spikes__DOT__weight__d__g_exc__X__exc_spikes__DOT__weight;
    double __P__g_exc__X__exc_spikes__DOT__weight__d__g_exc__X__exc_spikes__DOT__weight__d;
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
    Buffers_(iaf_cond_alpha_neuron_nestml &);
    Buffers_(const Buffers_ &, iaf_cond_alpha_neuron_nestml &);

    /**
     * Logger for all analog data
    **/
    nest::UniversalDataLogger<iaf_cond_alpha_neuron_nestml> logger_;

    // -----------------------------------------------------------------------
    //   Spike buffers and sums of incoming spikes/currents per timestep
    // -----------------------------------------------------------------------

    // input port: exc_spikes
    nest::RingBuffer spike_input_exc_spikes__DOT__weight_;
    double spike_input_exc_spikes__DOT__weight_grid_sum_;
    nest::RingBuffer spike_input_exc_spikes_spike_input_received_;
    double spike_input_exc_spikes_spike_input_received_grid_sum_;

    // input port: inh_spikes
    nest::RingBuffer spike_input_inh_spikes__DOT__weight_;
    double spike_input_inh_spikes__DOT__weight_grid_sum_;
    nest::RingBuffer spike_input_inh_spikes_spike_input_received_;
    double spike_input_inh_spikes_spike_input_received_grid_sum_;

    // -----------------------------------------------------------------------
    //   Continuous-input buffers
    // -----------------------------------------------------------------------



    nest::RingBuffer
     I_stim;   //!< Buffer for input (type: pA)

    inline nest::RingBuffer& get_I_stim() {
        return I_stim;
    }

    double I_stim_grid_sum_;

    // -----------------------------------------------------------------------
    //   GSL ODE solver data structures
    // -----------------------------------------------------------------------

    gsl_odeiv_step* __s;    //!< stepping function
    gsl_odeiv_control* __c; //!< adaptive stepsize control function
    gsl_odeiv_evolve* __e;  //!< evolution function
    gsl_odeiv_system __sys; //!< struct describing system

    // __integration_step should be reset with the neuron on ResetNetwork,
    // but remain unchanged during calibration. Since it is initialized with
    // step_, and the resolution cannot change after nodes have been created,
    // it is safe to place both here.
    double __step;             //!< step size in ms
    double __integration_step; //!< current integration time step, updated by GSL
  };

  // -------------------------------------------------------------------------
  //   Getters/setters for inline expressions
  // -------------------------------------------------------------------------

  inline double get_I_syn_exc() const
  {
    return S_.ode_state[State_::g_exc__X__exc_spikes__DOT__weight] * (S_.ode_state[State_::V_m] - P_.E_exc);
  }

  inline double get_I_syn_inh() const
  {
    return S_.ode_state[State_::g_inh__X__inh_spikes__DOT__weight] * (S_.ode_state[State_::V_m] - P_.E_inh);
  }

  inline double get_I_leak() const
  {
    return P_.g_L * (S_.ode_state[State_::V_m] - P_.E_L);
  }



  // -------------------------------------------------------------------------
  //   Getters/setters for spike input buffers
  // -------------------------------------------------------------------------

  // input port: exc_spikes
  inline nest::RingBuffer& get_spike_input_exc_spikes__DOT__weight()
  {
    return B_.spike_input_exc_spikes__DOT__weight_;
  }

  inline double get_spike_input_exc_spikes__DOT__weight_grid_sum_()
  {
    return B_.spike_input_exc_spikes__DOT__weight_grid_sum_;
  }

  inline nest::RingBuffer& get_spike_input_exc_spikes_spike_input_received_()
  {
    return B_.spike_input_exc_spikes_spike_input_received_;
  }

  inline double get_spike_input_exc_spikes_spike_input_received_grid_sum_()
  {
    return B_.spike_input_exc_spikes_spike_input_received_grid_sum_;
  }

  // input port: inh_spikes
  inline nest::RingBuffer& get_spike_input_inh_spikes__DOT__weight()
  {
    return B_.spike_input_inh_spikes__DOT__weight_;
  }

  inline double get_spike_input_inh_spikes__DOT__weight_grid_sum_()
  {
    return B_.spike_input_inh_spikes__DOT__weight_grid_sum_;
  }

  inline nest::RingBuffer& get_spike_input_inh_spikes_spike_input_received_()
  {
    return B_.spike_input_inh_spikes_spike_input_received_;
  }

  inline double get_spike_input_inh_spikes_spike_input_received_grid_sum_()
  {
    return B_.spike_input_inh_spikes_spike_input_received_grid_sum_;
  }


  // -------------------------------------------------------------------------
  //   Getters/setters for continuous-time input buffers
  // -------------------------------------------------------------------------



  inline nest::RingBuffer& get_I_stim() {
      return B_.get_I_stim();
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
  static nest::RecordablesMap<iaf_cond_alpha_neuron_nestml> recordablesMap_;
  friend int iaf_cond_alpha_neuron_nestml_dynamics_refr_t( double, const double ode_state[], double f[], void* pnode );
  friend int iaf_cond_alpha_neuron_nestml_dynamics_V_m( double, const double ode_state[], double f[], void* pnode );
  friend int iaf_cond_alpha_neuron_nestml_dynamics( double, const double ode_state[], double f[], void* pnode );

}; /* neuron iaf_cond_alpha_neuron_nestml */

inline nest_port_t iaf_cond_alpha_neuron_nestml::send_test_event(nest::Node& target, nest_rport_t receptor_type, nest::synindex, bool)
{
  // You should usually not change the code in this function.
  // It confirms that the target of connection @c c accepts @c nest::SpikeEvent on
  // the given receptor_type.
  nest::SpikeEvent e;
  e.set_sender(*this);
  return target.handles_test_event(e, receptor_type);
}

inline nest_port_t iaf_cond_alpha_neuron_nestml::handles_test_event(nest::SpikeEvent&, nest_port_t receptor_type)
{
    if ( receptor_type < 0 or receptor_type >= MAX_SPIKE_RECEPTOR )
    {
      throw nest::UnknownReceptorType( receptor_type, get_name() );
    }
    return receptor_type;
}

inline nest_port_t iaf_cond_alpha_neuron_nestml::handles_test_event(nest::CurrentEvent&, nest_port_t receptor_type)
{
  // You should usually not change the code in this function.
  // It confirms to the connection management system that we are able
  // to handle a CurrentEvent on port 0. You need to extend the function
  // if you want to differentiate between input ports.
  if (receptor_type != 0)
  {
    throw nest::UnknownReceptorType(receptor_type, get_name());
  }
  return 0;
}

inline nest_port_t iaf_cond_alpha_neuron_nestml::handles_test_event(nest::DataLoggingRequest& dlr, nest_port_t receptor_type)
{
  // You should usually not change the code in this function.
  // It confirms to the connection management system that we are able
  // to handle a DataLoggingRequest on port 0.
  // The function also tells the built-in UniversalDataLogger that this node
  // is recorded from and that it thus needs to collect data during simulation.
  if (receptor_type != 0)
  {
    throw nest::UnknownReceptorType(receptor_type, get_name());
  }

  return B_.logger_.connect_logging_device(dlr, recordablesMap_);
}

inline void iaf_cond_alpha_neuron_nestml::get_status(DictionaryDatum &__d) const
{
  // parameters
  def< double >(__d, nest::iaf_cond_alpha_neuron_nestml_names::_C_m, get_C_m());
  def< double >(__d, nest::iaf_cond_alpha_neuron_nestml_names::_g_L, get_g_L());
  def< double >(__d, nest::iaf_cond_alpha_neuron_nestml_names::_E_L, get_E_L());
  def< double >(__d, nest::iaf_cond_alpha_neuron_nestml_names::_refr_T, get_refr_T());
  def< double >(__d, nest::iaf_cond_alpha_neuron_nestml_names::_V_th, get_V_th());
  def< double >(__d, nest::iaf_cond_alpha_neuron_nestml_names::_V_reset, get_V_reset());
  def< double >(__d, nest::iaf_cond_alpha_neuron_nestml_names::_E_exc, get_E_exc());
  def< double >(__d, nest::iaf_cond_alpha_neuron_nestml_names::_E_inh, get_E_inh());
  def< bool >(__d, nest::iaf_cond_alpha_neuron_nestml_names::_use_offset_trick, get_use_offset_trick());
  def< double >(__d, nest::iaf_cond_alpha_neuron_nestml_names::_tau_syn_exc, get_tau_syn_exc());
  def< double >(__d, nest::iaf_cond_alpha_neuron_nestml_names::_tau_syn_inh, get_tau_syn_inh());
  def< double >(__d, nest::iaf_cond_alpha_neuron_nestml_names::_I_e, get_I_e());

  // initial values for state variables in ODE or kernel
  def< double >(__d, nest::iaf_cond_alpha_neuron_nestml_names::_V_m, get_V_m());
  def< double >(__d, nest::iaf_cond_alpha_neuron_nestml_names::_refr_t, get_refr_t());
  def< double >(__d, nest::iaf_cond_alpha_neuron_nestml_names::_g_inh__X__inh_spikes__DOT__weight, get_g_inh__X__inh_spikes__DOT__weight());
  def< double >(__d, nest::iaf_cond_alpha_neuron_nestml_names::_g_inh__X__inh_spikes__DOT__weight__d, get_g_inh__X__inh_spikes__DOT__weight__d());
  def< double >(__d, nest::iaf_cond_alpha_neuron_nestml_names::_g_exc__X__exc_spikes__DOT__weight, get_g_exc__X__exc_spikes__DOT__weight());
  def< double >(__d, nest::iaf_cond_alpha_neuron_nestml_names::_g_exc__X__exc_spikes__DOT__weight__d, get_g_exc__X__exc_spikes__DOT__weight__d());

  ExtendedPostHistoryArchivingNode::get_status( __d );
  DictionaryDatum __receptor_type = new Dictionary();
( *__receptor_type )[ "EXC_SPIKES" ] = 1;
( *__receptor_type )[ "INH_SPIKES" ] = 2;

    ( *__d )[ "receptor_types" ] = __receptor_type;

  (*__d)[nest::names::recordables] = recordablesMap_.get_list();
  def< double >(__d, nest::iaf_cond_alpha_neuron_nestml_names::gsl_abs_error_tol, P_.__gsl_abs_error_tol);
  if ( P_.__gsl_abs_error_tol <= 0. ){
    throw nest::BadProperty( "The gsl_abs_error_tol must be strictly positive." );
  }
  def< double >(__d, nest::iaf_cond_alpha_neuron_nestml_names::gsl_rel_error_tol, P_.__gsl_rel_error_tol);
  if ( P_.__gsl_rel_error_tol < 0. ){
    throw nest::BadProperty( "The gsl_rel_error_tol must be zero or positive." );
  }
}

inline void iaf_cond_alpha_neuron_nestml::set_status(const DictionaryDatum &__d)
{
  // parameters
  double tmp_C_m = get_C_m();
  nest::updateValueParam<double>(__d, nest::iaf_cond_alpha_neuron_nestml_names::_C_m, tmp_C_m, this);
  double tmp_g_L = get_g_L();
  nest::updateValueParam<double>(__d, nest::iaf_cond_alpha_neuron_nestml_names::_g_L, tmp_g_L, this);
  double tmp_E_L = get_E_L();
  nest::updateValueParam<double>(__d, nest::iaf_cond_alpha_neuron_nestml_names::_E_L, tmp_E_L, this);
  double tmp_refr_T = get_refr_T();
  nest::updateValueParam<double>(__d, nest::iaf_cond_alpha_neuron_nestml_names::_refr_T, tmp_refr_T, this);
  double tmp_V_th = get_V_th();
  nest::updateValueParam<double>(__d, nest::iaf_cond_alpha_neuron_nestml_names::_V_th, tmp_V_th, this);
  double tmp_V_reset = get_V_reset();
  nest::updateValueParam<double>(__d, nest::iaf_cond_alpha_neuron_nestml_names::_V_reset, tmp_V_reset, this);
  double tmp_E_exc = get_E_exc();
  nest::updateValueParam<double>(__d, nest::iaf_cond_alpha_neuron_nestml_names::_E_exc, tmp_E_exc, this);
  bool tmp_use_offset_trick = get_use_offset_trick();
  nest::updateValueParam<bool>(__d, nest::iaf_cond_alpha_neuron_nestml_names::_use_offset_trick, tmp_use_offset_trick, this);
  double tmp_E_inh = get_E_inh();
  nest::updateValueParam<double>(__d, nest::iaf_cond_alpha_neuron_nestml_names::_E_inh, tmp_E_inh, this);
  double tmp_tau_syn_exc = get_tau_syn_exc();
  nest::updateValueParam<double>(__d, nest::iaf_cond_alpha_neuron_nestml_names::_tau_syn_exc, tmp_tau_syn_exc, this);
  double tmp_tau_syn_inh = get_tau_syn_inh();
  nest::updateValueParam<double>(__d, nest::iaf_cond_alpha_neuron_nestml_names::_tau_syn_inh, tmp_tau_syn_inh, this);
  double tmp_I_e = get_I_e();
  nest::updateValueParam<double>(__d, nest::iaf_cond_alpha_neuron_nestml_names::_I_e, tmp_I_e, this);

  // initial values for state variables in ODE or kernel
  double tmp_V_m = get_V_m();
  nest::updateValueParam<double>(__d, nest::iaf_cond_alpha_neuron_nestml_names::_V_m, tmp_V_m, this);
  double tmp_refr_t = get_refr_t();
  nest::updateValueParam<double>(__d, nest::iaf_cond_alpha_neuron_nestml_names::_refr_t, tmp_refr_t, this);
  double tmp_g_inh__X__inh_spikes__DOT__weight = get_g_inh__X__inh_spikes__DOT__weight();
  nest::updateValueParam<double>(__d, nest::iaf_cond_alpha_neuron_nestml_names::_g_inh__X__inh_spikes__DOT__weight, tmp_g_inh__X__inh_spikes__DOT__weight, this);
  double tmp_g_inh__X__inh_spikes__DOT__weight__d = get_g_inh__X__inh_spikes__DOT__weight__d();
  nest::updateValueParam<double>(__d, nest::iaf_cond_alpha_neuron_nestml_names::_g_inh__X__inh_spikes__DOT__weight__d, tmp_g_inh__X__inh_spikes__DOT__weight__d, this);
  double tmp_g_exc__X__exc_spikes__DOT__weight = get_g_exc__X__exc_spikes__DOT__weight();
  nest::updateValueParam<double>(__d, nest::iaf_cond_alpha_neuron_nestml_names::_g_exc__X__exc_spikes__DOT__weight, tmp_g_exc__X__exc_spikes__DOT__weight, this);
  double tmp_g_exc__X__exc_spikes__DOT__weight__d = get_g_exc__X__exc_spikes__DOT__weight__d();
  nest::updateValueParam<double>(__d, nest::iaf_cond_alpha_neuron_nestml_names::_g_exc__X__exc_spikes__DOT__weight__d, tmp_g_exc__X__exc_spikes__DOT__weight__d, this);

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  ExtendedPostHistoryArchivingNode::set_status(__d);

  // if we get here, temporaries contain consistent set of properties
  set_C_m(tmp_C_m);
  set_g_L(tmp_g_L);
  set_E_L(tmp_E_L);
  set_refr_T(tmp_refr_T);
  set_V_th(tmp_V_th);
  set_V_reset(tmp_V_reset);
  set_E_exc(tmp_E_exc);
  set_use_offset_trick(tmp_use_offset_trick);
  set_E_inh(tmp_E_inh);
  set_tau_syn_exc(tmp_tau_syn_exc);
  set_tau_syn_inh(tmp_tau_syn_inh);
  set_I_e(tmp_I_e);
  set_V_m(tmp_V_m);
  set_refr_t(tmp_refr_t);
  set_g_inh__X__inh_spikes__DOT__weight(tmp_g_inh__X__inh_spikes__DOT__weight);
  set_g_inh__X__inh_spikes__DOT__weight__d(tmp_g_inh__X__inh_spikes__DOT__weight__d);
  set_g_exc__X__exc_spikes__DOT__weight(tmp_g_exc__X__exc_spikes__DOT__weight);
  set_g_exc__X__exc_spikes__DOT__weight__d(tmp_g_exc__X__exc_spikes__DOT__weight__d);




  updateValue< double >(__d, nest::iaf_cond_alpha_neuron_nestml_names::gsl_abs_error_tol, P_.__gsl_abs_error_tol);
  if ( P_.__gsl_abs_error_tol <= 0. )
  {
    throw nest::BadProperty( "The gsl_abs_error_tol must be strictly positive." );
  }
  updateValue< double >(__d, nest::iaf_cond_alpha_neuron_nestml_names::gsl_rel_error_tol, P_.__gsl_rel_error_tol);
  if ( P_.__gsl_rel_error_tol < 0. )
  {
    throw nest::BadProperty( "The gsl_rel_error_tol must be zero or positive." );
  }

  // recompute internal variables in case they are dependent on parameters or state that might have been updated in this call to set_status()
  recompute_internal_variables();
};



#endif /* #ifndef IAF_COND_ALPHA_NEURON_NESTML */
