
/**
 *  eglif_cond_alpha_multisyn.h
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
#ifndef EGLIF_COND_ALPHA_MULTISYN
#define EGLIF_COND_ALPHA_MULTISYN

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
#include "archiving_node.h"
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
namespace eglif_cond_alpha_multisyn_names
{
    const Name _V_m( "V_m" );
    const Name _I_dep( "I_dep" );
    const Name _I_adap( "I_adap" );
    const Name _r( "r" );
    const Name _syn_kernel2__X__syn2_spike( "syn_kernel2__X__syn2_spike" );
    const Name _syn_kernel2__X__syn2_spike__d( "syn_kernel2__X__syn2_spike__d" );
    const Name _syn_kernel4__X__syn4_spike( "syn_kernel4__X__syn4_spike" );
    const Name _syn_kernel4__X__syn4_spike__d( "syn_kernel4__X__syn4_spike__d" );
    const Name _syn_kernel1__X__syn1_spike( "syn_kernel1__X__syn1_spike" );
    const Name _syn_kernel1__X__syn1_spike__d( "syn_kernel1__X__syn1_spike__d" );
    const Name _syn_kernel3__X__syn3_spike( "syn_kernel3__X__syn3_spike" );
    const Name _syn_kernel3__X__syn3_spike__d( "syn_kernel3__X__syn3_spike__d" );
    const Name _I_syn( "I_syn" );
    const Name _C_m( "C_m" );
    const Name _tau_m( "tau_m" );
    const Name _E_L( "E_L" );
    const Name _t_ref( "t_ref" );
    const Name _I_e( "I_e" );
    const Name _V_min( "V_min" );
    const Name _V_th( "V_th" );
    const Name _lambda_0( "lambda_0" );
    const Name _tau_V( "tau_V" );
    const Name _V_reset( "V_reset" );
    const Name _k_1( "k_1" );
    const Name _k_2( "k_2" );
    const Name _k_adap( "k_adap" );
    const Name _A1( "A1" );
    const Name _A2( "A2" );
    const Name _E_rev1( "E_rev1" );
    const Name _tau_syn1( "tau_syn1" );
    const Name _E_rev2( "E_rev2" );
    const Name _tau_syn2( "tau_syn2" );
    const Name _E_rev3( "E_rev3" );
    const Name _tau_syn3( "tau_syn3" );
    const Name _E_rev4( "E_rev4" );
    const Name _tau_syn4( "tau_syn4" );
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
**/
extern "C" inline int eglif_cond_alpha_multisyn_dynamics( double, const double y[], double f[], void* pnode );


#include "nest_time.h"
  typedef size_t nest_port_t;
  typedef size_t nest_rport_t;

/* BeginDocumentation
  Name: eglif_cond_alpha_multisyn

  Description:

    """
  eglif_cond_alpha_multisyn - Conductance based extended generalized leaky integrate and fire neuron model#######################################################################################################

  Description
  +++++++++++

  eglif_cond_alpha_multisyn is the generalized leaky integrate and fire neuron according to
  Geminiani et al. (2018) [#geminiani_2018]_, with post-synaptic conductances in the form of a
  bi-exponential ("alpha") function.

  The membrane potential is given by the following differential equation:

  .. math::
     \begin{cases}
         C_m \dfrac{dV_m}{dt} = \dfrac{C_m}{\tau_m}(V_m-E_L)
                                  - I_{adap} + I_{dep} + I_e + I_{stim} + I_{syn}\\
         \dfrac{dI_{adap}}{dt}=k_{adap}(V_m-E_L)-k_2 \cdot I_{adap}\\
         \dfrac{dI_{dep}}{dt} =-k_1 \cdot I_{dep}\\
     \end{cases}

  where the synaptic current :math:`I_{syn}` integrates the input from 4 postsynaptic receptors:

  .. math::
         I_{syn} = \sum_{i=1} ^{4} g_i (V - E_{rev,i})

  Here, the synapse `i` is excitatory or inhibitory depending on the value of :math:`E_{rev,i}`.

  Neuron produces spikes stochastically according to a point process with the firing intensity:

  .. math::
    \lambda = \lambda_0 \cdot e^{\dfrac{V_m -V_{th}}{\tau_V}}

  In case of spike emission, the spike-triggered adaptation currents :math:`I_{adap}` and
  :math:`I_{dep}` are respectively increased and set by their respective constant
  (which can be positive or negative):

  .. math::
     V_m &= V_{reset}\\
     I_{dep} &= A1\\
     I_{adap} &= I_{adap} + A2

  References
  ++++++++++

  .. start-references

  .. [#geminiani_2018] Geminiani, A., Casellato, C., Locatelli, F., Prestori, F., Pedrocchi, A.,
     & D'Angelo, E. (2018). Complex dynamics in simplified neuronal models:
     reproducing Golgi cell electroresponsiveness. Frontiers in neuroinformatics, 12, 88.
     https://doi.org/10.3389/fninf.2018.00088


  See also
  ++++++++

  aeif_cond_alpha, aeif_cond_exp
  """


  Parameters:
  The following parameters can be set in the status dictionary.
C_m [pF]  Membrane Capacitance
tau_m [ms]  Membrane time constant
E_L [mV]  Leak reversal Potential (aka resting potential)
t_ref [ms]  Refractory period
I_e [pA]  Constant endogenous current
V_min [mV]  Minimal membrane potential value
V_th [mV]  Spike Threshold
lambda_0 [real]  Escape rate parameter
tau_V [mV]  Escape rate parameter
V_reset [mV]  Reset potential
k_1 [real]  Decay rate
k_2 [real]  Adaptation constant
k_adap [real]  Adaptation constant
A1 [pA]  Current update constant
A2 [pA]  Current update constant
E_rev1 [mV]  Postsynaptic receptor reversal Potential
tau_syn1 [ms]  Postsynaptic receptor time constant
E_rev2 [mV]  Postsynaptic receptor reversal Potential
tau_syn2 [ms]  Postsynaptic receptor time constant
E_rev3 [mV]  Postsynaptic receptor reversal Potential
tau_syn3 [ms]  Postsynaptic receptor time constant
E_rev4 [mV]  Postsynaptic receptor reversal Potential
tau_syn4 [ms]  Postsynaptic receptor time constant


  Dynamic state variables:
V_m [mV]  Membrane potential in mV
I_dep [pA]  Depolarizing spike-triggered current
I_adap [pA]  Adaptation current
r [integer]  Refractory state


  Sends: nest::SpikeEvent

  Receives: Spike, Current, DataLoggingRequest
*/

// Register the neuron model
void register_eglif_cond_alpha_multisyn( const std::string& name );

class eglif_cond_alpha_multisyn : public nest::ArchivingNode
{
public:
  /**
   * The constructor is only used to create the model prototype in the model manager.
  **/
  eglif_cond_alpha_multisyn();

  /**
   * The copy constructor is used to create model copies and instances of the model.
   * @node The copy constructor needs to initialize the parameters and the state.
   *       Initialization of buffers and interal variables is deferred to
   *       @c init_buffers_() and @c pre_run_hook() (or calibrate() in NEST 3.3 and older).
  **/
  eglif_cond_alpha_multisyn(const eglif_cond_alpha_multisyn &);

  /**
   * Destructor.
  **/
  ~eglif_cond_alpha_multisyn() override;

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

  inline double get_I_dep() const
  {
    return S_.ode_state[State_::I_dep];
  }

  inline void set_I_dep(const double __v)
  {
    S_.ode_state[State_::I_dep] = __v;
  }

  inline double get_I_adap() const
  {
    return S_.ode_state[State_::I_adap];
  }

  inline void set_I_adap(const double __v)
  {
    S_.ode_state[State_::I_adap] = __v;
  }

  inline long get_r() const
  {
    return S_.r;
  }

  inline void set_r(const long __v)
  {
    S_.r = __v;
  }

  inline double get_syn_kernel2__X__syn2_spike() const
  {
    return S_.ode_state[State_::syn_kernel2__X__syn2_spike];
  }

  inline void set_syn_kernel2__X__syn2_spike(const double __v)
  {
    S_.ode_state[State_::syn_kernel2__X__syn2_spike] = __v;
  }

  inline double get_syn_kernel2__X__syn2_spike__d() const
  {
    return S_.ode_state[State_::syn_kernel2__X__syn2_spike__d];
  }

  inline void set_syn_kernel2__X__syn2_spike__d(const double __v)
  {
    S_.ode_state[State_::syn_kernel2__X__syn2_spike__d] = __v;
  }

  inline double get_syn_kernel4__X__syn4_spike() const
  {
    return S_.ode_state[State_::syn_kernel4__X__syn4_spike];
  }

  inline void set_syn_kernel4__X__syn4_spike(const double __v)
  {
    S_.ode_state[State_::syn_kernel4__X__syn4_spike] = __v;
  }

  inline double get_syn_kernel4__X__syn4_spike__d() const
  {
    return S_.ode_state[State_::syn_kernel4__X__syn4_spike__d];
  }

  inline void set_syn_kernel4__X__syn4_spike__d(const double __v)
  {
    S_.ode_state[State_::syn_kernel4__X__syn4_spike__d] = __v;
  }

  inline double get_syn_kernel1__X__syn1_spike() const
  {
    return S_.ode_state[State_::syn_kernel1__X__syn1_spike];
  }

  inline void set_syn_kernel1__X__syn1_spike(const double __v)
  {
    S_.ode_state[State_::syn_kernel1__X__syn1_spike] = __v;
  }

  inline double get_syn_kernel1__X__syn1_spike__d() const
  {
    return S_.ode_state[State_::syn_kernel1__X__syn1_spike__d];
  }

  inline void set_syn_kernel1__X__syn1_spike__d(const double __v)
  {
    S_.ode_state[State_::syn_kernel1__X__syn1_spike__d] = __v;
  }

  inline double get_syn_kernel3__X__syn3_spike() const
  {
    return S_.ode_state[State_::syn_kernel3__X__syn3_spike];
  }

  inline void set_syn_kernel3__X__syn3_spike(const double __v)
  {
    S_.ode_state[State_::syn_kernel3__X__syn3_spike] = __v;
  }

  inline double get_syn_kernel3__X__syn3_spike__d() const
  {
    return S_.ode_state[State_::syn_kernel3__X__syn3_spike__d];
  }

  inline void set_syn_kernel3__X__syn3_spike__d(const double __v)
  {
    S_.ode_state[State_::syn_kernel3__X__syn3_spike__d] = __v;
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

  inline double get_tau_m() const
  {
    return P_.tau_m;
  }

  inline void set_tau_m(const double __v)
  {
    P_.tau_m = __v;
  }

  inline double get_E_L() const
  {
    return P_.E_L;
  }

  inline void set_E_L(const double __v)
  {
    P_.E_L = __v;
  }

  inline double get_t_ref() const
  {
    return P_.t_ref;
  }

  inline void set_t_ref(const double __v)
  {
    P_.t_ref = __v;
  }

  inline double get_I_e() const
  {
    return P_.I_e;
  }

  inline void set_I_e(const double __v)
  {
    P_.I_e = __v;
  }

  inline double get_V_min() const
  {
    return P_.V_min;
  }

  inline void set_V_min(const double __v)
  {
    P_.V_min = __v;
  }

  inline double get_V_th() const
  {
    return P_.V_th;
  }

  inline void set_V_th(const double __v)
  {
    P_.V_th = __v;
  }

  inline double get_lambda_0() const
  {
    return P_.lambda_0;
  }

  inline void set_lambda_0(const double __v)
  {
    P_.lambda_0 = __v;
  }

  inline double get_tau_V() const
  {
    return P_.tau_V;
  }

  inline void set_tau_V(const double __v)
  {
    P_.tau_V = __v;
  }

  inline double get_V_reset() const
  {
    return P_.V_reset;
  }

  inline void set_V_reset(const double __v)
  {
    P_.V_reset = __v;
  }

  inline double get_k_1() const
  {
    return P_.k_1;
  }

  inline void set_k_1(const double __v)
  {
    P_.k_1 = __v;
  }

  inline double get_k_2() const
  {
    return P_.k_2;
  }

  inline void set_k_2(const double __v)
  {
    P_.k_2 = __v;
  }

  inline double get_k_adap() const
  {
    return P_.k_adap;
  }

  inline void set_k_adap(const double __v)
  {
    P_.k_adap = __v;
  }

  inline double get_A1() const
  {
    return P_.A1;
  }

  inline void set_A1(const double __v)
  {
    P_.A1 = __v;
  }

  inline double get_A2() const
  {
    return P_.A2;
  }

  inline void set_A2(const double __v)
  {
    P_.A2 = __v;
  }

  inline double get_E_rev1() const
  {
    return P_.E_rev1;
  }

  inline void set_E_rev1(const double __v)
  {
    P_.E_rev1 = __v;
  }

  inline double get_tau_syn1() const
  {
    return P_.tau_syn1;
  }

  inline void set_tau_syn1(const double __v)
  {
    P_.tau_syn1 = __v;
  }

  inline double get_E_rev2() const
  {
    return P_.E_rev2;
  }

  inline void set_E_rev2(const double __v)
  {
    P_.E_rev2 = __v;
  }

  inline double get_tau_syn2() const
  {
    return P_.tau_syn2;
  }

  inline void set_tau_syn2(const double __v)
  {
    P_.tau_syn2 = __v;
  }

  inline double get_E_rev3() const
  {
    return P_.E_rev3;
  }

  inline void set_E_rev3(const double __v)
  {
    P_.E_rev3 = __v;
  }

  inline double get_tau_syn3() const
  {
    return P_.tau_syn3;
  }

  inline void set_tau_syn3(const double __v)
  {
    P_.tau_syn3 = __v;
  }

  inline double get_E_rev4() const
  {
    return P_.E_rev4;
  }

  inline void set_E_rev4(const double __v)
  {
    P_.E_rev4 = __v;
  }

  inline double get_tau_syn4() const
  {
    return P_.tau_syn4;
  }

  inline void set_tau_syn4(const double __v)
  {
    P_.tau_syn4 = __v;
  }


  // -------------------------------------------------------------------------
  //   Getters/setters for internals
  // -------------------------------------------------------------------------

  inline long get_RefractoryCounts() const
  {
    return V_.RefractoryCounts;
  }

  inline void set_RefractoryCounts(const long __v)
  {
    V_.RefractoryCounts = __v;
  }
  inline double get___h() const
  {
    return V_.__h;
  }

  inline void set___h(const double __v)
  {
    V_.__h = __v;
  }
  inline double get___P__I_dep__I_dep() const
  {
    return V_.__P__I_dep__I_dep;
  }

  inline void set___P__I_dep__I_dep(const double __v)
  {
    V_.__P__I_dep__I_dep = __v;
  }
  inline double get___P__syn_kernel2__X__syn2_spike__syn_kernel2__X__syn2_spike() const
  {
    return V_.__P__syn_kernel2__X__syn2_spike__syn_kernel2__X__syn2_spike;
  }

  inline void set___P__syn_kernel2__X__syn2_spike__syn_kernel2__X__syn2_spike(const double __v)
  {
    V_.__P__syn_kernel2__X__syn2_spike__syn_kernel2__X__syn2_spike = __v;
  }
  inline double get___P__syn_kernel2__X__syn2_spike__syn_kernel2__X__syn2_spike__d() const
  {
    return V_.__P__syn_kernel2__X__syn2_spike__syn_kernel2__X__syn2_spike__d;
  }

  inline void set___P__syn_kernel2__X__syn2_spike__syn_kernel2__X__syn2_spike__d(const double __v)
  {
    V_.__P__syn_kernel2__X__syn2_spike__syn_kernel2__X__syn2_spike__d = __v;
  }
  inline double get___P__syn_kernel2__X__syn2_spike__d__syn_kernel2__X__syn2_spike() const
  {
    return V_.__P__syn_kernel2__X__syn2_spike__d__syn_kernel2__X__syn2_spike;
  }

  inline void set___P__syn_kernel2__X__syn2_spike__d__syn_kernel2__X__syn2_spike(const double __v)
  {
    V_.__P__syn_kernel2__X__syn2_spike__d__syn_kernel2__X__syn2_spike = __v;
  }
  inline double get___P__syn_kernel2__X__syn2_spike__d__syn_kernel2__X__syn2_spike__d() const
  {
    return V_.__P__syn_kernel2__X__syn2_spike__d__syn_kernel2__X__syn2_spike__d;
  }

  inline void set___P__syn_kernel2__X__syn2_spike__d__syn_kernel2__X__syn2_spike__d(const double __v)
  {
    V_.__P__syn_kernel2__X__syn2_spike__d__syn_kernel2__X__syn2_spike__d = __v;
  }
  inline double get___P__syn_kernel4__X__syn4_spike__syn_kernel4__X__syn4_spike() const
  {
    return V_.__P__syn_kernel4__X__syn4_spike__syn_kernel4__X__syn4_spike;
  }

  inline void set___P__syn_kernel4__X__syn4_spike__syn_kernel4__X__syn4_spike(const double __v)
  {
    V_.__P__syn_kernel4__X__syn4_spike__syn_kernel4__X__syn4_spike = __v;
  }
  inline double get___P__syn_kernel4__X__syn4_spike__syn_kernel4__X__syn4_spike__d() const
  {
    return V_.__P__syn_kernel4__X__syn4_spike__syn_kernel4__X__syn4_spike__d;
  }

  inline void set___P__syn_kernel4__X__syn4_spike__syn_kernel4__X__syn4_spike__d(const double __v)
  {
    V_.__P__syn_kernel4__X__syn4_spike__syn_kernel4__X__syn4_spike__d = __v;
  }
  inline double get___P__syn_kernel4__X__syn4_spike__d__syn_kernel4__X__syn4_spike() const
  {
    return V_.__P__syn_kernel4__X__syn4_spike__d__syn_kernel4__X__syn4_spike;
  }

  inline void set___P__syn_kernel4__X__syn4_spike__d__syn_kernel4__X__syn4_spike(const double __v)
  {
    V_.__P__syn_kernel4__X__syn4_spike__d__syn_kernel4__X__syn4_spike = __v;
  }
  inline double get___P__syn_kernel4__X__syn4_spike__d__syn_kernel4__X__syn4_spike__d() const
  {
    return V_.__P__syn_kernel4__X__syn4_spike__d__syn_kernel4__X__syn4_spike__d;
  }

  inline void set___P__syn_kernel4__X__syn4_spike__d__syn_kernel4__X__syn4_spike__d(const double __v)
  {
    V_.__P__syn_kernel4__X__syn4_spike__d__syn_kernel4__X__syn4_spike__d = __v;
  }
  inline double get___P__syn_kernel1__X__syn1_spike__syn_kernel1__X__syn1_spike() const
  {
    return V_.__P__syn_kernel1__X__syn1_spike__syn_kernel1__X__syn1_spike;
  }

  inline void set___P__syn_kernel1__X__syn1_spike__syn_kernel1__X__syn1_spike(const double __v)
  {
    V_.__P__syn_kernel1__X__syn1_spike__syn_kernel1__X__syn1_spike = __v;
  }
  inline double get___P__syn_kernel1__X__syn1_spike__syn_kernel1__X__syn1_spike__d() const
  {
    return V_.__P__syn_kernel1__X__syn1_spike__syn_kernel1__X__syn1_spike__d;
  }

  inline void set___P__syn_kernel1__X__syn1_spike__syn_kernel1__X__syn1_spike__d(const double __v)
  {
    V_.__P__syn_kernel1__X__syn1_spike__syn_kernel1__X__syn1_spike__d = __v;
  }
  inline double get___P__syn_kernel1__X__syn1_spike__d__syn_kernel1__X__syn1_spike() const
  {
    return V_.__P__syn_kernel1__X__syn1_spike__d__syn_kernel1__X__syn1_spike;
  }

  inline void set___P__syn_kernel1__X__syn1_spike__d__syn_kernel1__X__syn1_spike(const double __v)
  {
    V_.__P__syn_kernel1__X__syn1_spike__d__syn_kernel1__X__syn1_spike = __v;
  }
  inline double get___P__syn_kernel1__X__syn1_spike__d__syn_kernel1__X__syn1_spike__d() const
  {
    return V_.__P__syn_kernel1__X__syn1_spike__d__syn_kernel1__X__syn1_spike__d;
  }

  inline void set___P__syn_kernel1__X__syn1_spike__d__syn_kernel1__X__syn1_spike__d(const double __v)
  {
    V_.__P__syn_kernel1__X__syn1_spike__d__syn_kernel1__X__syn1_spike__d = __v;
  }
  inline double get___P__syn_kernel3__X__syn3_spike__syn_kernel3__X__syn3_spike() const
  {
    return V_.__P__syn_kernel3__X__syn3_spike__syn_kernel3__X__syn3_spike;
  }

  inline void set___P__syn_kernel3__X__syn3_spike__syn_kernel3__X__syn3_spike(const double __v)
  {
    V_.__P__syn_kernel3__X__syn3_spike__syn_kernel3__X__syn3_spike = __v;
  }
  inline double get___P__syn_kernel3__X__syn3_spike__syn_kernel3__X__syn3_spike__d() const
  {
    return V_.__P__syn_kernel3__X__syn3_spike__syn_kernel3__X__syn3_spike__d;
  }

  inline void set___P__syn_kernel3__X__syn3_spike__syn_kernel3__X__syn3_spike__d(const double __v)
  {
    V_.__P__syn_kernel3__X__syn3_spike__syn_kernel3__X__syn3_spike__d = __v;
  }
  inline double get___P__syn_kernel3__X__syn3_spike__d__syn_kernel3__X__syn3_spike() const
  {
    return V_.__P__syn_kernel3__X__syn3_spike__d__syn_kernel3__X__syn3_spike;
  }

  inline void set___P__syn_kernel3__X__syn3_spike__d__syn_kernel3__X__syn3_spike(const double __v)
  {
    V_.__P__syn_kernel3__X__syn3_spike__d__syn_kernel3__X__syn3_spike = __v;
  }
  inline double get___P__syn_kernel3__X__syn3_spike__d__syn_kernel3__X__syn3_spike__d() const
  {
    return V_.__P__syn_kernel3__X__syn3_spike__d__syn_kernel3__X__syn3_spike__d;
  }

  inline void set___P__syn_kernel3__X__syn3_spike__d__syn_kernel3__X__syn3_spike__d(const double __v)
  {
    V_.__P__syn_kernel3__X__syn3_spike__d__syn_kernel3__X__syn3_spike__d = __v;
  }


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
    SYN1_SPIKE = 1,
    SYN2_SPIKE = 2,
    SYN3_SPIKE = 3,
    SYN4_SPIKE = 4,
    MAX_SPIKE_RECEPTOR = 5
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
  friend class nest::RecordablesMap<eglif_cond_alpha_multisyn>;
  friend class nest::UniversalDataLogger<eglif_cond_alpha_multisyn>;

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
    //!  Membrane Capacitance
    double C_m;
    //!  Membrane time constant
    double tau_m;
    //!  Leak reversal Potential (aka resting potential)
    double E_L;
    //!  Refractory period
    double t_ref;
    //!  Constant endogenous current
    double I_e;
    //!  Minimal membrane potential value
    double V_min;
    //!  Spike Threshold
    double V_th;
    //!  Escape rate parameter
    double lambda_0;
    //!  Escape rate parameter
    double tau_V;
    //!  Reset potential
    double V_reset;
    //!  Decay rate
    double k_1;
    //!  Adaptation constant
    double k_2;
    //!  Adaptation constant
    double k_adap;
    //!  Current update constant
    double A1;
    //!  Current update constant
    double A2;
    //!  Postsynaptic receptor reversal Potential
    double E_rev1;
    //!  Postsynaptic receptor time constant
    double tau_syn1;
    //!  Postsynaptic receptor reversal Potential
    double E_rev2;
    //!  Postsynaptic receptor time constant
    double tau_syn2;
    //!  Postsynaptic receptor reversal Potential
    double E_rev3;
    //!  Postsynaptic receptor time constant
    double tau_syn3;
    //!  Postsynaptic receptor reversal Potential
    double E_rev4;
    //!  Postsynaptic receptor time constant
    double tau_syn4;

    double __gsl_error_tol;

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
//!  Refractory state
long r;
    //! Symbolic indices to the elements of the state vector y
    enum StateVecElems
    {
      V_m,
      I_dep,
      I_adap,
      syn_kernel2__X__syn2_spike,
      syn_kernel2__X__syn2_spike__d,
      syn_kernel4__X__syn4_spike,
      syn_kernel4__X__syn4_spike__d,
      syn_kernel1__X__syn1_spike,
      syn_kernel1__X__syn1_spike__d,
      syn_kernel3__X__syn3_spike,
      syn_kernel3__X__syn3_spike__d,
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
    //!  refractory time in steps
    long RefractoryCounts;
    double __h;
    double __P__I_dep__I_dep;
    double __P__syn_kernel2__X__syn2_spike__syn_kernel2__X__syn2_spike;
    double __P__syn_kernel2__X__syn2_spike__syn_kernel2__X__syn2_spike__d;
    double __P__syn_kernel2__X__syn2_spike__d__syn_kernel2__X__syn2_spike;
    double __P__syn_kernel2__X__syn2_spike__d__syn_kernel2__X__syn2_spike__d;
    double __P__syn_kernel4__X__syn4_spike__syn_kernel4__X__syn4_spike;
    double __P__syn_kernel4__X__syn4_spike__syn_kernel4__X__syn4_spike__d;
    double __P__syn_kernel4__X__syn4_spike__d__syn_kernel4__X__syn4_spike;
    double __P__syn_kernel4__X__syn4_spike__d__syn_kernel4__X__syn4_spike__d;
    double __P__syn_kernel1__X__syn1_spike__syn_kernel1__X__syn1_spike;
    double __P__syn_kernel1__X__syn1_spike__syn_kernel1__X__syn1_spike__d;
    double __P__syn_kernel1__X__syn1_spike__d__syn_kernel1__X__syn1_spike;
    double __P__syn_kernel1__X__syn1_spike__d__syn_kernel1__X__syn1_spike__d;
    double __P__syn_kernel3__X__syn3_spike__syn_kernel3__X__syn3_spike;
    double __P__syn_kernel3__X__syn3_spike__syn_kernel3__X__syn3_spike__d;
    double __P__syn_kernel3__X__syn3_spike__d__syn_kernel3__X__syn3_spike;
    double __P__syn_kernel3__X__syn3_spike__d__syn_kernel3__X__syn3_spike__d;
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
    Buffers_(eglif_cond_alpha_multisyn &);
    Buffers_(const Buffers_ &, eglif_cond_alpha_multisyn &);

    /**
     * Logger for all analog data
    **/
    nest::UniversalDataLogger<eglif_cond_alpha_multisyn> logger_;

    // -----------------------------------------------------------------------
    //   Buffers and sums of incoming spikes/currents per timestep
    // -----------------------------------------------------------------------
    // Buffer containing the incoming spikes
    

inline std::vector< nest::RingBuffer >& get_spike_inputs_()
{
    return spike_inputs_;
}
std::vector< nest::RingBuffer > spike_inputs_;

    // Buffer containing the sum of all the incoming spikes
    

inline std::vector< double >& get_spike_inputs_grid_sum_()
{
    return spike_inputs_grid_sum_;
}
std::vector< double > spike_inputs_grid_sum_;

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
  inline double get_I_syn1() const
  {
    return S_.ode_state[State_::syn_kernel1__X__syn1_spike] * 1.0 * (P_.E_rev1 - S_.ode_state[State_::V_m]);
  }

  inline double get_I_syn2() const
  {
    return S_.ode_state[State_::syn_kernel2__X__syn2_spike] * 1.0 * (P_.E_rev2 - S_.ode_state[State_::V_m]);
  }

  inline double get_I_syn3() const
  {
    return S_.ode_state[State_::syn_kernel3__X__syn3_spike] * 1.0 * (P_.E_rev3 - S_.ode_state[State_::V_m]);
  }

  inline double get_I_syn4() const
  {
    return S_.ode_state[State_::syn_kernel4__X__syn4_spike] * 1.0 * (P_.E_rev4 - S_.ode_state[State_::V_m]);
  }

  inline double get_I_syn() const
  {
    return (S_.ode_state[State_::syn_kernel1__X__syn1_spike] * 1.0 * (P_.E_rev1 - S_.ode_state[State_::V_m])) + (S_.ode_state[State_::syn_kernel2__X__syn2_spike] * 1.0 * (P_.E_rev2 - S_.ode_state[State_::V_m])) + (S_.ode_state[State_::syn_kernel3__X__syn3_spike] * 1.0 * (P_.E_rev3 - S_.ode_state[State_::V_m])) + (S_.ode_state[State_::syn_kernel4__X__syn4_spike] * 1.0 * (P_.E_rev4 - S_.ode_state[State_::V_m]));
  }



  // -------------------------------------------------------------------------
  //   Getters/setters for input buffers
  // -------------------------------------------------------------------------

  // Buffer containing the incoming spikes
  

inline std::vector< nest::RingBuffer >& get_spike_inputs_()
{
    return B_.get_spike_inputs_();
}

  

inline std::vector< double >& get_spike_inputs_grid_sum_()
{
    return B_.get_spike_inputs_grid_sum_();
}
  
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
  static nest::RecordablesMap<eglif_cond_alpha_multisyn> recordablesMap_;
  friend int eglif_cond_alpha_multisyn_dynamics( double, const double y[], double f[], void* pnode );


}; /* neuron eglif_cond_alpha_multisyn */

inline nest_port_t eglif_cond_alpha_multisyn::send_test_event(nest::Node& target, nest_rport_t receptor_type, nest::synindex, bool)
{
  // You should usually not change the code in this function.
  // It confirms that the target of connection @c c accepts @c nest::SpikeEvent on
  // the given @c receptor_type.
  nest::SpikeEvent e;
  e.set_sender(*this);
  return target.handles_test_event(e, receptor_type);
}

inline nest_port_t eglif_cond_alpha_multisyn::handles_test_event(nest::SpikeEvent&, nest_port_t receptor_type)
{
    assert( B_.spike_inputs_.size() == NUM_SPIKE_RECEPTORS );
    if ( receptor_type < MIN_SPIKE_RECEPTOR or receptor_type >= MAX_SPIKE_RECEPTOR )
    {
      throw nest::UnknownReceptorType( receptor_type, get_name() );
    }
    return receptor_type - MIN_SPIKE_RECEPTOR;
}

inline nest_port_t eglif_cond_alpha_multisyn::handles_test_event(nest::CurrentEvent&, nest_port_t receptor_type)
{
  // You should usually not change the code in this function.
  // It confirms to the connection management system that we are able
  // to handle @c CurrentEvent on port 0. You need to extend the function
  // if you want to differentiate between input ports.
  if (receptor_type != 0)
  {
    throw nest::UnknownReceptorType(receptor_type, get_name());
  }
  return 0;
}

inline nest_port_t eglif_cond_alpha_multisyn::handles_test_event(nest::DataLoggingRequest& dlr, nest_port_t receptor_type)
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

inline void eglif_cond_alpha_multisyn::get_status(DictionaryDatum &__d) const
{
  // parameters
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_C_m, get_C_m());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_tau_m, get_tau_m());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_E_L, get_E_L());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_t_ref, get_t_ref());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_I_e, get_I_e());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_V_min, get_V_min());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_V_th, get_V_th());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_lambda_0, get_lambda_0());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_tau_V, get_tau_V());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_V_reset, get_V_reset());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_k_1, get_k_1());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_k_2, get_k_2());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_k_adap, get_k_adap());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_A1, get_A1());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_A2, get_A2());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_E_rev1, get_E_rev1());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_tau_syn1, get_tau_syn1());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_E_rev2, get_E_rev2());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_tau_syn2, get_tau_syn2());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_E_rev3, get_E_rev3());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_tau_syn3, get_tau_syn3());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_E_rev4, get_E_rev4());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_tau_syn4, get_tau_syn4());

  // initial values for state variables in ODE or kernel
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_V_m, get_V_m());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_I_dep, get_I_dep());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_I_adap, get_I_adap());
  def<long>(__d, nest::eglif_cond_alpha_multisyn_names::_r, get_r());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_syn_kernel2__X__syn2_spike, get_syn_kernel2__X__syn2_spike());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_syn_kernel2__X__syn2_spike__d, get_syn_kernel2__X__syn2_spike__d());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_syn_kernel4__X__syn4_spike, get_syn_kernel4__X__syn4_spike());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_syn_kernel4__X__syn4_spike__d, get_syn_kernel4__X__syn4_spike__d());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_syn_kernel1__X__syn1_spike, get_syn_kernel1__X__syn1_spike());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_syn_kernel1__X__syn1_spike__d, get_syn_kernel1__X__syn1_spike__d());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_syn_kernel3__X__syn3_spike, get_syn_kernel3__X__syn3_spike());
  def<double>(__d, nest::eglif_cond_alpha_multisyn_names::_syn_kernel3__X__syn3_spike__d, get_syn_kernel3__X__syn3_spike__d());

  ArchivingNode::get_status( __d );
  DictionaryDatum __receptor_type = new Dictionary();
    ( *__receptor_type )[ "SYN1_SPIKE" ] = 1;
    ( *__receptor_type )[ "SYN2_SPIKE" ] = 2;
    ( *__receptor_type )[ "SYN3_SPIKE" ] = 3;
    ( *__receptor_type )[ "SYN4_SPIKE" ] = 4;
    ( *__d )[ "receptor_types" ] = __receptor_type;

  (*__d)[nest::names::recordables] = recordablesMap_.get_list();
  def< double >(__d, nest::names::gsl_error_tol, P_.__gsl_error_tol);
  if ( P_.__gsl_error_tol <= 0. ){
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }
}

inline void eglif_cond_alpha_multisyn::set_status(const DictionaryDatum &__d)
{
  // parameters
  double tmp_C_m = get_C_m();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_C_m, tmp_C_m, this);
  double tmp_tau_m = get_tau_m();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_tau_m, tmp_tau_m, this);
  double tmp_E_L = get_E_L();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_E_L, tmp_E_L, this);
  double tmp_t_ref = get_t_ref();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_t_ref, tmp_t_ref, this);
  double tmp_I_e = get_I_e();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_I_e, tmp_I_e, this);
  double tmp_V_min = get_V_min();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_V_min, tmp_V_min, this);
  double tmp_V_th = get_V_th();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_V_th, tmp_V_th, this);
  double tmp_lambda_0 = get_lambda_0();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_lambda_0, tmp_lambda_0, this);
  double tmp_tau_V = get_tau_V();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_tau_V, tmp_tau_V, this);
  double tmp_V_reset = get_V_reset();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_V_reset, tmp_V_reset, this);
  double tmp_k_1 = get_k_1();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_k_1, tmp_k_1, this);
  double tmp_k_2 = get_k_2();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_k_2, tmp_k_2, this);
  double tmp_k_adap = get_k_adap();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_k_adap, tmp_k_adap, this);
  double tmp_A1 = get_A1();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_A1, tmp_A1, this);
  double tmp_A2 = get_A2();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_A2, tmp_A2, this);
  double tmp_E_rev1 = get_E_rev1();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_E_rev1, tmp_E_rev1, this);
  double tmp_tau_syn1 = get_tau_syn1();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_tau_syn1, tmp_tau_syn1, this);
  double tmp_E_rev2 = get_E_rev2();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_E_rev2, tmp_E_rev2, this);
  double tmp_tau_syn2 = get_tau_syn2();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_tau_syn2, tmp_tau_syn2, this);
  double tmp_E_rev3 = get_E_rev3();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_E_rev3, tmp_E_rev3, this);
  double tmp_tau_syn3 = get_tau_syn3();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_tau_syn3, tmp_tau_syn3, this);
  double tmp_E_rev4 = get_E_rev4();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_E_rev4, tmp_E_rev4, this);
  double tmp_tau_syn4 = get_tau_syn4();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_tau_syn4, tmp_tau_syn4, this);

  // initial values for state variables in ODE or kernel
  double tmp_V_m = get_V_m();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_V_m, tmp_V_m, this);
  double tmp_I_dep = get_I_dep();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_I_dep, tmp_I_dep, this);
  double tmp_I_adap = get_I_adap();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_I_adap, tmp_I_adap, this);
  long tmp_r = get_r();
  nest::updateValueParam<long>(__d, nest::eglif_cond_alpha_multisyn_names::_r, tmp_r, this);
  double tmp_syn_kernel2__X__syn2_spike = get_syn_kernel2__X__syn2_spike();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_syn_kernel2__X__syn2_spike, tmp_syn_kernel2__X__syn2_spike, this);
  double tmp_syn_kernel2__X__syn2_spike__d = get_syn_kernel2__X__syn2_spike__d();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_syn_kernel2__X__syn2_spike__d, tmp_syn_kernel2__X__syn2_spike__d, this);
  double tmp_syn_kernel4__X__syn4_spike = get_syn_kernel4__X__syn4_spike();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_syn_kernel4__X__syn4_spike, tmp_syn_kernel4__X__syn4_spike, this);
  double tmp_syn_kernel4__X__syn4_spike__d = get_syn_kernel4__X__syn4_spike__d();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_syn_kernel4__X__syn4_spike__d, tmp_syn_kernel4__X__syn4_spike__d, this);
  double tmp_syn_kernel1__X__syn1_spike = get_syn_kernel1__X__syn1_spike();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_syn_kernel1__X__syn1_spike, tmp_syn_kernel1__X__syn1_spike, this);
  double tmp_syn_kernel1__X__syn1_spike__d = get_syn_kernel1__X__syn1_spike__d();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_syn_kernel1__X__syn1_spike__d, tmp_syn_kernel1__X__syn1_spike__d, this);
  double tmp_syn_kernel3__X__syn3_spike = get_syn_kernel3__X__syn3_spike();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_syn_kernel3__X__syn3_spike, tmp_syn_kernel3__X__syn3_spike, this);
  double tmp_syn_kernel3__X__syn3_spike__d = get_syn_kernel3__X__syn3_spike__d();
  nest::updateValueParam<double>(__d, nest::eglif_cond_alpha_multisyn_names::_syn_kernel3__X__syn3_spike__d, tmp_syn_kernel3__X__syn3_spike__d, this);

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  ArchivingNode::set_status(__d);

  // if we get here, temporaries contain consistent set of properties
  set_C_m(tmp_C_m);
  set_tau_m(tmp_tau_m);
  set_E_L(tmp_E_L);
  set_t_ref(tmp_t_ref);
  set_I_e(tmp_I_e);
  set_V_min(tmp_V_min);
  set_V_th(tmp_V_th);
  set_lambda_0(tmp_lambda_0);
  set_tau_V(tmp_tau_V);
  set_V_reset(tmp_V_reset);
  set_k_1(tmp_k_1);
  set_k_2(tmp_k_2);
  set_k_adap(tmp_k_adap);
  set_A1(tmp_A1);
  set_A2(tmp_A2);
  set_E_rev1(tmp_E_rev1);
  set_tau_syn1(tmp_tau_syn1);
  set_E_rev2(tmp_E_rev2);
  set_tau_syn2(tmp_tau_syn2);
  set_E_rev3(tmp_E_rev3);
  set_tau_syn3(tmp_tau_syn3);
  set_E_rev4(tmp_E_rev4);
  set_tau_syn4(tmp_tau_syn4);
  set_V_m(tmp_V_m);
  set_I_dep(tmp_I_dep);
  set_I_adap(tmp_I_adap);
  set_r(tmp_r);
  set_syn_kernel2__X__syn2_spike(tmp_syn_kernel2__X__syn2_spike);
  set_syn_kernel2__X__syn2_spike__d(tmp_syn_kernel2__X__syn2_spike__d);
  set_syn_kernel4__X__syn4_spike(tmp_syn_kernel4__X__syn4_spike);
  set_syn_kernel4__X__syn4_spike__d(tmp_syn_kernel4__X__syn4_spike__d);
  set_syn_kernel1__X__syn1_spike(tmp_syn_kernel1__X__syn1_spike);
  set_syn_kernel1__X__syn1_spike__d(tmp_syn_kernel1__X__syn1_spike__d);
  set_syn_kernel3__X__syn3_spike(tmp_syn_kernel3__X__syn3_spike);
  set_syn_kernel3__X__syn3_spike__d(tmp_syn_kernel3__X__syn3_spike__d);




  updateValue< double >(__d, nest::names::gsl_error_tol, P_.__gsl_error_tol);
  if ( P_.__gsl_error_tol <= 0. )
  {
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }

  // recompute internal variables in case they are dependent on parameters or state that might have been updated in this call to set_status()
  recompute_internal_variables();
};



#endif /* #ifndef EGLIF_COND_ALPHA_MULTISYN */
