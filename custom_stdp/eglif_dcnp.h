
/**
 *  eglif_dcnp.h
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
 *  Generated from NESTML at time: 2024-10-21 06:27:55.250947
**/
#ifndef EGLIF_DCNP
#define EGLIF_DCNP

#ifndef HAVE_LIBLTDL
//#error "NEST was compiled without support for dynamic loading. Please install libltdl and recompile NEST."
#endif

// C++ includes:
#include <cmath>

#include "config.h"

// Includes for random number generator
#include <random>

#ifndef HAVE_GSL
//error "The GSL library is required for the Runge-Kutta solver."
#endif

// External includes:
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// Includes from nestkernel:
#include "archiving_node.h"
#include "extended_post_history_archiving_node.h"
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
   void register_eglif_dcnp( const std::string& name );
}
namespace eglif_dcnp_names
{
    const Name _V_m( "V_m" );
    const Name _I_adap( "I_adap" );
    const Name _I_dep( "I_dep" );
    const Name _r( "r" );
    const Name _lambda( "lambda" );
    const Name _offset( "offset" );
    const Name _tick( "tick" );
    const Name _cf_buffer( "cf_buffer" );
    const Name _gr_buffer( "gr_buffer" );
    const Name _last_io( "last_io" );
    const Name _g4__X__rec4( "g4__X__rec4" );
    const Name _g4__X__rec4__d( "g4__X__rec4__d" );
    const Name _g2__X__rec2( "g2__X__rec2" );
    const Name _g2__X__rec2__d( "g2__X__rec2__d" );
    const Name _g1__X__rec1( "g1__X__rec1" );
    const Name _g1__X__rec1__d( "g1__X__rec1__d" );
    const Name _g3__X__rec3( "g3__X__rec3" );
    const Name _g3__X__rec3__d( "g3__X__rec3__d" );
    const Name _I_syn( "I_syn" );
    const Name _I_tot( "I_tot" );
    const Name _C_m( "C_m" );
    const Name _tau_m( "tau_m" );
    const Name _E_L( "E_L" );
    const Name _t_ref( "t_ref" );
    const Name _V_reset( "V_reset" );
    const Name _V_th( "V_th" );
    const Name _Vmin( "Vmin" );
    const Name _I_e( "I_e" );
    const Name _Vinit( "Vinit" );
    const Name _lambda_0( "lambda_0" );
    const Name _tau_V( "tau_V" );
    const Name _kadap( "kadap" );
    const Name _k2( "k2" );
    const Name _k1( "k1" );
    const Name _A1( "A1" );
    const Name _A2( "A2" );
    const Name _E_rev1( "E_rev1" );
    const Name _E_rev2( "E_rev2" );
    const Name _E_rev3( "E_rev3" );
    const Name _E_rev4( "E_rev4" );
    const Name _tau_syn1( "tau_syn1" );
    const Name _tau_syn2( "tau_syn2" );
    const Name _tau_syn3( "tau_syn3" );
    const Name _tau_syn4( "tau_syn4" );
    const Name _complex_flag( "complex_flag");

}
extern "C" inline int eglif_dcnp_dynamics( double, const double ode_state[], double f[], void* pnode );

#include "nest_time.h"
  typedef size_t nest_port_t;
  typedef size_t nest_rport_t;


class eglif_dcnp : public nest::ExtendedPostHistoryArchivingNode
{
public:
  /**
   * The constructor is only used to create the model prototype in the model manager.
  **/
  eglif_dcnp();

  /**
   * The copy constructor is used to create model copies and instances of the model.
   * @node The copy constructor needs to initialize the parameters and the state.
   *       Initialization of buffers and interal variables is deferred to
   *       @c init_buffers_() and @c pre_run_hook() (or calibrate() in NEST 3.3 and older).
  **/
  eglif_dcnp(const eglif_dcnp &);

  /**
   * Destructor.
  **/
  ~eglif_dcnp() override;

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

   bool
    is_off_grid() const override
    {
    return true;
    }
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

  inline double get_I_adap() const
  {
    return S_.ode_state[State_::I_adap];
  }

  inline void set_I_adap(const double __v)
  {
    S_.ode_state[State_::I_adap] = __v;
  }

  inline double get_I_dep() const
  {
    return S_.ode_state[State_::I_dep];
  }

  inline void set_I_dep(const double __v)
  {
    S_.ode_state[State_::I_dep] = __v;
  }

  inline long get_r() const
  {
    return S_.r;
  }

  inline void set_r(const long __v)
  {
    S_.r = __v;
  }

  inline double get_lambda() const
  {
    return S_.lambda;
  }

  inline void set_lambda(const double __v)
  {
    S_.lambda = __v;
  }

  inline double get_cf_buffer() const
  {
    return S_.cf_buffer;
  }
  
  inline void set_cf_buffer(const double __v)
  {
    S_.cf_buffer = __v;
  }

  inline double get_gr_buffer() const
  {
    return S_.gr_buffer;
  }

  inline void set_gr_buffer(const double __v)
  {
    S_.gr_buffer = __v;
  }

  inline double get_last_io() const
  {
    return S_.last_io;
  }

  inline void set_last_io(const double __v)
  {
    S_.last_io = __v;
  }

  inline double get_tick() const
  {
    return S_.tick;
  }

  inline void set_tick(const double __v)
  {
    S_.tick = __v;
  }

  inline double get_complex_flag() const
  {
    return S_.complex_flag;
  }

  inline void set_complex_flag(const double __v)
  {
    S_.complex_flag = __v;
  }

  inline double get_g4__X__rec4() const
  {
    return S_.ode_state[State_::g4__X__rec4];
  }

  inline void set_g4__X__rec4(const double __v)
  {
    S_.ode_state[State_::g4__X__rec4] = __v;
  }

  inline double get_g4__X__rec4__d() const
  {
    return S_.ode_state[State_::g4__X__rec4__d];
  }

  inline void set_g4__X__rec4__d(const double __v)
  {
    S_.ode_state[State_::g4__X__rec4__d] = __v;
  }

  inline double get_g2__X__rec2() const
  {
    return S_.ode_state[State_::g2__X__rec2];
  }

  inline void set_g2__X__rec2(const double __v)
  {
    S_.ode_state[State_::g2__X__rec2] = __v;
  }

  inline double get_g2__X__rec2__d() const
  {
    return S_.ode_state[State_::g2__X__rec2__d];
  }

  inline void set_g2__X__rec2__d(const double __v)
  {
    S_.ode_state[State_::g2__X__rec2__d] = __v;
  }

  inline double get_g1__X__rec1() const
  {
    return S_.ode_state[State_::g1__X__rec1];
  }

  inline void set_g1__X__rec1(const double __v)
  {
    S_.ode_state[State_::g1__X__rec1] = __v;
  }

  inline double get_g1__X__rec1__d() const
  {
    return S_.ode_state[State_::g1__X__rec1__d];
  }

  inline void set_g1__X__rec1__d(const double __v)
  {
    S_.ode_state[State_::g1__X__rec1__d] = __v;
  }

  inline double get_g3__X__rec3() const
  {
    return S_.ode_state[State_::g3__X__rec3];
  }

  inline void set_g3__X__rec3(const double __v)
  {
    S_.ode_state[State_::g3__X__rec3] = __v;
  }

  inline double get_g3__X__rec3__d() const
  {
    return S_.ode_state[State_::g3__X__rec3__d];
  }

  inline void set_g3__X__rec3__d(const double __v)
  {
    S_.ode_state[State_::g3__X__rec3__d] = __v;
  }


  // -------------------------------------------------------------------------
  //   Getters/setters for parameters
  // -------------------------------------------------------------------------
  inline double get_offset() const
  {
    return P_.offset;
  }

  inline void set_offset(const double __v)
  {
    P_.offset = __v;
  }
  
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

  inline double get_V_reset() const
  {
    return P_.V_reset;
  }

  inline void set_V_reset(const double __v)
  {
    P_.V_reset = __v;
  }

  inline double get_V_th() const
  {
    return P_.V_th;
  }

  inline void set_V_th(const double __v)
  {
    P_.V_th = __v;
  }

  inline double get_Vmin() const
  {
    return P_.Vmin;
  }

  inline void set_Vmin(const double __v)
  {
    P_.Vmin = __v;
  }

  inline double get_I_e() const
  {
    return P_.I_e;
  }

  inline void set_I_e(const double __v)
  {
    P_.I_e = __v;
  }

  inline double get_Vinit() const
  {
    return P_.Vinit;
  }

  inline void set_Vinit(const double __v)
  {
    P_.Vinit = __v;
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

  inline double get_kadap() const
  {
    return P_.kadap;
  }

  inline void set_kadap(const double __v)
  {
    P_.kadap = __v;
  }

  inline double get_k2() const
  {
    return P_.k2;
  }

  inline void set_k2(const double __v)
  {
    P_.k2 = __v;
  }

  inline double get_k1() const
  {
    return P_.k1;
  }

  inline void set_k1(const double __v)
  {
    P_.k1 = __v;
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

  inline double get_E_rev2() const
  {
    return P_.E_rev2;
  }

  inline void set_E_rev2(const double __v)
  {
    P_.E_rev2 = __v;
  }

  inline double get_E_rev3() const
  {
    return P_.E_rev3;
  }

  inline void set_E_rev3(const double __v)
  {
    P_.E_rev3 = __v;
  }

  inline double get_E_rev4() const
  {
    return P_.E_rev4;
  }

  inline void set_E_rev4(const double __v)
  {
    P_.E_rev4 = __v;
  }

  inline double get_tau_syn1() const
  {
    return P_.tau_syn1;
  }

  inline void set_tau_syn1(const double __v)
  {
    P_.tau_syn1 = __v;
  }

  inline double get_tau_syn2() const
  {
    return P_.tau_syn2;
  }

  inline void set_tau_syn2(const double __v)
  {
    P_.tau_syn2 = __v;
  }

  inline double get_tau_syn3() const
  {
    return P_.tau_syn3;
  }

  inline void set_tau_syn3(const double __v)
  {
    P_.tau_syn3 = __v;
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
  inline double get___P__g4__X__rec4__g4__X__rec4() const
  {
    return V_.__P__g4__X__rec4__g4__X__rec4;
  }

  inline void set___P__g4__X__rec4__g4__X__rec4(const double __v)
  {
    V_.__P__g4__X__rec4__g4__X__rec4 = __v;
  }
  inline double get___P__g4__X__rec4__g4__X__rec4__d() const
  {
    return V_.__P__g4__X__rec4__g4__X__rec4__d;
  }

  inline void set___P__g4__X__rec4__g4__X__rec4__d(const double __v)
  {
    V_.__P__g4__X__rec4__g4__X__rec4__d = __v;
  }
  inline double get___P__g4__X__rec4__d__g4__X__rec4() const
  {
    return V_.__P__g4__X__rec4__d__g4__X__rec4;
  }

  inline void set___P__g4__X__rec4__d__g4__X__rec4(const double __v)
  {
    V_.__P__g4__X__rec4__d__g4__X__rec4 = __v;
  }
  inline double get___P__g4__X__rec4__d__g4__X__rec4__d() const
  {
    return V_.__P__g4__X__rec4__d__g4__X__rec4__d;
  }

  inline void set___P__g4__X__rec4__d__g4__X__rec4__d(const double __v)
  {
    V_.__P__g4__X__rec4__d__g4__X__rec4__d = __v;
  }
  inline double get___P__g2__X__rec2__g2__X__rec2() const
  {
    return V_.__P__g2__X__rec2__g2__X__rec2;
  }

  inline void set___P__g2__X__rec2__g2__X__rec2(const double __v)
  {
    V_.__P__g2__X__rec2__g2__X__rec2 = __v;
  }
  inline double get___P__g2__X__rec2__g2__X__rec2__d() const
  {
    return V_.__P__g2__X__rec2__g2__X__rec2__d;
  }

  inline void set___P__g2__X__rec2__g2__X__rec2__d(const double __v)
  {
    V_.__P__g2__X__rec2__g2__X__rec2__d = __v;
  }
  inline double get___P__g2__X__rec2__d__g2__X__rec2() const
  {
    return V_.__P__g2__X__rec2__d__g2__X__rec2;
  }

  inline void set___P__g2__X__rec2__d__g2__X__rec2(const double __v)
  {
    V_.__P__g2__X__rec2__d__g2__X__rec2 = __v;
  }
  inline double get___P__g2__X__rec2__d__g2__X__rec2__d() const
  {
    return V_.__P__g2__X__rec2__d__g2__X__rec2__d;
  }

  inline void set___P__g2__X__rec2__d__g2__X__rec2__d(const double __v)
  {
    V_.__P__g2__X__rec2__d__g2__X__rec2__d = __v;
  }
  inline double get___P__g1__X__rec1__g1__X__rec1() const
  {
    return V_.__P__g1__X__rec1__g1__X__rec1;
  }

  inline void set___P__g1__X__rec1__g1__X__rec1(const double __v)
  {
    V_.__P__g1__X__rec1__g1__X__rec1 = __v;
  }
  inline double get___P__g1__X__rec1__g1__X__rec1__d() const
  {
    return V_.__P__g1__X__rec1__g1__X__rec1__d;
  }

  inline void set___P__g1__X__rec1__g1__X__rec1__d(const double __v)
  {
    V_.__P__g1__X__rec1__g1__X__rec1__d = __v;
  }
  inline double get___P__g1__X__rec1__d__g1__X__rec1() const
  {
    return V_.__P__g1__X__rec1__d__g1__X__rec1;
  }

  inline void set___P__g1__X__rec1__d__g1__X__rec1(const double __v)
  {
    V_.__P__g1__X__rec1__d__g1__X__rec1 = __v;
  }
  inline double get___P__g1__X__rec1__d__g1__X__rec1__d() const
  {
    return V_.__P__g1__X__rec1__d__g1__X__rec1__d;
  }

  inline void set___P__g1__X__rec1__d__g1__X__rec1__d(const double __v)
  {
    V_.__P__g1__X__rec1__d__g1__X__rec1__d = __v;
  }
  inline double get___P__g3__X__rec3__g3__X__rec3() const
  {
    return V_.__P__g3__X__rec3__g3__X__rec3;
  }

  inline void set___P__g3__X__rec3__g3__X__rec3(const double __v)
  {
    V_.__P__g3__X__rec3__g3__X__rec3 = __v;
  }
  inline double get___P__g3__X__rec3__g3__X__rec3__d() const
  {
    return V_.__P__g3__X__rec3__g3__X__rec3__d;
  }

  inline void set___P__g3__X__rec3__g3__X__rec3__d(const double __v)
  {
    V_.__P__g3__X__rec3__g3__X__rec3__d = __v;
  }
  inline double get___P__g3__X__rec3__d__g3__X__rec3() const
  {
    return V_.__P__g3__X__rec3__d__g3__X__rec3;
  }

  inline void set___P__g3__X__rec3__d__g3__X__rec3(const double __v)
  {
    V_.__P__g3__X__rec3__d__g3__X__rec3 = __v;
  }
  inline double get___P__g3__X__rec3__d__g3__X__rec3__d() const
  {
    return V_.__P__g3__X__rec3__d__g3__X__rec3__d;
  }

  inline void set___P__g3__X__rec3__d__g3__X__rec3__d(const double __v)
  {
    V_.__P__g3__X__rec3__d__g3__X__rec3__d = __v;
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
    REC1 = 1,
    REC2 = 2,
    REC3 = 3,
    REC4 = 4,
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
  friend class nest::RecordablesMap<eglif_dcnp>;
  friend class nest::UniversalDataLogger<eglif_dcnp>;

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
    //!  Neurophysiological quantities
    //!  Membrane potential
    double C_m;
    //!  Membrane time constant
    double tau_m;
    //!  Membrane resting potential
    double E_L;
    //!  Refractory period duration
    double t_ref;
    //!  Reset potential
    double V_reset;
    //!  Spike generation threshold potential
    double V_th;
    //!  Lower bound on membrane voltage potential
    double Vmin;
    //!  Endogenous current
    double I_e;
    double Vinit;
    //!  Spike generation parameters
    double lambda_0;
    double tau_V;
    //!  Functional parameters to be optimized
    //!  Adaptation threshold
    double kadap;
    //!  Adaptation threshold
    double k2;
    //!  Idep decay rate
    double k1;
    //!  Update parameter for Idep
    double A1;
    //!  Update parameter for Iadap
    double A2;
    //!  Synaptic parameters -> 4 receptors
    double E_rev1;
    double E_rev2;
    double E_rev3;
    double E_rev4;
    double tau_syn1;
    double tau_syn2;
    double tau_syn3;
    double tau_syn4;

    double offset;

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
//!  Counter of ticks during refractory period
long r;
//!  Stochasticity function for spike generation
double lambda;
//!  Flag for simple vs complex spikes
double cf_buffer;
double gr_buffer;
double last_io;
double tick;
double complex_flag;
    //! Symbolic indices to the elements of the state vector y
    enum StateVecElems
    {
      V_m,
      I_dep,
      I_adap,
      g4__X__rec4,
      g4__X__rec4__d,
      g2__X__rec2,
      g2__X__rec2__d,
      g1__X__rec1,
      g1__X__rec1__d,
      g3__X__rec3,
      g3__X__rec3__d,
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
    //!  Duration of the refractory period in simulation steps
    long RefractoryCounts;
    double __h;
    double __P__I_dep__I_dep;
    double __P__g4__X__rec4__g4__X__rec4;
    double __P__g4__X__rec4__g4__X__rec4__d;
    double __P__g4__X__rec4__d__g4__X__rec4;
    double __P__g4__X__rec4__d__g4__X__rec4__d;
    double __P__g2__X__rec2__g2__X__rec2;
    double __P__g2__X__rec2__g2__X__rec2__d;
    double __P__g2__X__rec2__d__g2__X__rec2;
    double __P__g2__X__rec2__d__g2__X__rec2__d;
    double __P__g1__X__rec1__g1__X__rec1;
    double __P__g1__X__rec1__g1__X__rec1__d;
    double __P__g1__X__rec1__d__g1__X__rec1;
    double __P__g1__X__rec1__d__g1__X__rec1__d;
    double __P__g3__X__rec3__g3__X__rec3;
    double __P__g3__X__rec3__g3__X__rec3__d;
    double __P__g3__X__rec3__d__g3__X__rec3;
    double __P__g3__X__rec3__d__g3__X__rec3__d;
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
    Buffers_(eglif_dcnp &);
    Buffers_(const Buffers_ &, eglif_dcnp &);

    /**
     * Logger for all analog data
    **/
    nest::UniversalDataLogger<eglif_dcnp> logger_;

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
  inline double get_I_syn() const
  {
    return S_.ode_state[State_::g1__X__rec1] * 1.0 * (P_.E_rev1 - S_.ode_state[State_::V_m]) + S_.ode_state[State_::g2__X__rec2] * 1.0 * (P_.E_rev2 - S_.ode_state[State_::V_m]) + S_.ode_state[State_::g3__X__rec3] * 1.0 * (P_.E_rev3 - S_.ode_state[State_::V_m]) + S_.ode_state[State_::g4__X__rec4] * 1.0 * (P_.E_rev4 - S_.ode_state[State_::V_m]);
  }

  inline double get_I_tot() const
  {
    return S_.ode_state[State_::I_dep] - S_.ode_state[State_::I_adap] + P_.I_e + B_.I_stim_grid_sum_ + (S_.ode_state[State_::g1__X__rec1] * 1.0 * (P_.E_rev1 - S_.ode_state[State_::V_m]) + S_.ode_state[State_::g2__X__rec2] * 1.0 * (P_.E_rev2 - S_.ode_state[State_::V_m]) + S_.ode_state[State_::g3__X__rec3] * 1.0 * (P_.E_rev3 - S_.ode_state[State_::V_m]) + S_.ode_state[State_::g4__X__rec4] * 1.0 * (P_.E_rev4 - S_.ode_state[State_::V_m]));
  }



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
  static nest::RecordablesMap<eglif_dcnp> recordablesMap_;
  friend int eglif_dcnp_dynamics( double, const double ode_state[], double f[], void* pnode );
  nest::normal_distribution normal_dev_; //!< random deviate generator
  nest::poisson_distribution poisson_dev_; //!< random deviate generator

}; /* neuron eglif_dcnp */

inline nest_port_t eglif_dcnp::send_test_event(nest::Node& target, nest_rport_t receptor_type, nest::synindex, bool)
{
  // You should usually not change the code in this function.
  // It confirms that the target of connection @c c accepts @c nest::SpikeEvent on
  // the given @c receptor_type.
  nest::SpikeEvent e;
  e.set_sender(*this);
  return target.handles_test_event(e, receptor_type);
}

inline nest_port_t eglif_dcnp::handles_test_event(nest::SpikeEvent&, nest_port_t receptor_type)
{
    assert( B_.spike_inputs_.size() == NUM_SPIKE_RECEPTORS );
    if ( receptor_type < MIN_SPIKE_RECEPTOR or receptor_type >= MAX_SPIKE_RECEPTOR )
    {
      throw nest::UnknownReceptorType( receptor_type, get_name() );
    }
    return receptor_type - MIN_SPIKE_RECEPTOR;
}

inline nest_port_t eglif_dcnp::handles_test_event(nest::CurrentEvent&, nest_port_t receptor_type)
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

inline nest_port_t eglif_dcnp::handles_test_event(nest::DataLoggingRequest& dlr, nest_port_t receptor_type)
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

inline void eglif_dcnp::get_status(DictionaryDatum &__d) const
{
  // parameters
  def< double >(__d, eglif_dcnp_names::_C_m, get_C_m());
  def< double >(__d, eglif_dcnp_names::_tau_m, get_tau_m());
  def< double >(__d, eglif_dcnp_names::_E_L, get_E_L());
  def< double >(__d, eglif_dcnp_names::_t_ref, get_t_ref());
  def< double >(__d, eglif_dcnp_names::_V_reset, get_V_reset());
  def< double >(__d, eglif_dcnp_names::_V_th, get_V_th());
  def< double >(__d, eglif_dcnp_names::_Vmin, get_Vmin());
  def< double >(__d, eglif_dcnp_names::_I_e, get_I_e());
  def< double >(__d, eglif_dcnp_names::_Vinit, get_Vinit());
  def< double >(__d, eglif_dcnp_names::_lambda_0, get_lambda_0());
  def< double >(__d, eglif_dcnp_names::_tau_V, get_tau_V());
  def< double >(__d, eglif_dcnp_names::_kadap, get_kadap());
  def< double >(__d, eglif_dcnp_names::_k2, get_k2());
  def< double >(__d, eglif_dcnp_names::_k1, get_k1());
  def< double >(__d, eglif_dcnp_names::_A1, get_A1());
  def< double >(__d, eglif_dcnp_names::_A2, get_A2());
  def< double >(__d, eglif_dcnp_names::_E_rev1, get_E_rev1());
  def< double >(__d, eglif_dcnp_names::_E_rev2, get_E_rev2());
  def< double >(__d, eglif_dcnp_names::_E_rev3, get_E_rev3());
  def< double >(__d, eglif_dcnp_names::_E_rev4, get_E_rev4());
  def< double >(__d, eglif_dcnp_names::_tau_syn1, get_tau_syn1());
  def< double >(__d, eglif_dcnp_names::_tau_syn2, get_tau_syn2());
  def< double >(__d, eglif_dcnp_names::_tau_syn3, get_tau_syn3());
  def< double >(__d, eglif_dcnp_names::_tau_syn4, get_tau_syn4());
  def< double >(__d, eglif_dcnp_names::_offset, get_offset());

  // initial values for state variables in ODE or kernel
  def< double >(__d, eglif_dcnp_names::_V_m, get_V_m());
  def< double >(__d, eglif_dcnp_names::_I_adap, get_I_adap());
  def< double >(__d, eglif_dcnp_names::_I_dep, get_I_dep());
  def< long >(__d, eglif_dcnp_names::_r, get_r());
  def< double >(__d, eglif_dcnp_names::_lambda, get_lambda());
  def< double >(__d, eglif_dcnp_names::_cf_buffer, get_cf_buffer());
  def< double >(__d, eglif_dcnp_names::_gr_buffer, get_gr_buffer());
  def< double >(__d, eglif_dcnp_names::_last_io, get_last_io());
  def< double >(__d, eglif_dcnp_names::_tick, get_tick());
  def< double >(__d, eglif_dcnp_names::_complex_flag, get_complex_flag());
  def< double >(__d, eglif_dcnp_names::_g4__X__rec4, get_g4__X__rec4());
  def< double >(__d, eglif_dcnp_names::_g4__X__rec4__d, get_g4__X__rec4__d());
  def< double >(__d, eglif_dcnp_names::_g2__X__rec2, get_g2__X__rec2());
  def< double >(__d, eglif_dcnp_names::_g2__X__rec2__d, get_g2__X__rec2__d());
  def< double >(__d, eglif_dcnp_names::_g1__X__rec1, get_g1__X__rec1());
  def< double >(__d, eglif_dcnp_names::_g1__X__rec1__d, get_g1__X__rec1__d());
  def< double >(__d, eglif_dcnp_names::_g3__X__rec3, get_g3__X__rec3());
  def< double >(__d, eglif_dcnp_names::_g3__X__rec3__d, get_g3__X__rec3__d());

  ExtendedPostHistoryArchivingNode::get_status( __d );
  DictionaryDatum __receptor_type = new Dictionary();
    ( *__receptor_type )[ "REC1" ] = 1;
    ( *__receptor_type )[ "REC2" ] = 2;
    ( *__receptor_type )[ "REC3" ] = 3;
    ( *__receptor_type )[ "REC4" ] = 4;
    ( *__d )[ "receptor_types" ] = __receptor_type;

  (*__d)[nest::names::recordables] = recordablesMap_.get_list();
  def< double >(__d, nest::names::gsl_error_tol, P_.__gsl_error_tol);
  if ( P_.__gsl_error_tol <= 0. ){
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }
}

inline void eglif_dcnp::set_status(const DictionaryDatum &__d)
{
  // parameters
  double tmp_C_m = get_C_m();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_C_m, tmp_C_m, this);
  double tmp_tau_m = get_tau_m();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_tau_m, tmp_tau_m, this);
  double tmp_E_L = get_E_L();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_E_L, tmp_E_L, this);
  double tmp_t_ref = get_t_ref();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_t_ref, tmp_t_ref, this);
  double tmp_V_reset = get_V_reset();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_V_reset, tmp_V_reset, this);
  double tmp_V_th = get_V_th();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_V_th, tmp_V_th, this);
  double tmp_Vmin = get_Vmin();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_Vmin, tmp_Vmin, this);
  double tmp_I_e = get_I_e();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_I_e, tmp_I_e, this);
  double tmp_Vinit = get_Vinit();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_Vinit, tmp_Vinit, this);
  double tmp_lambda_0 = get_lambda_0();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_lambda_0, tmp_lambda_0, this);
  double tmp_tau_V = get_tau_V();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_tau_V, tmp_tau_V, this);
  double tmp_kadap = get_kadap();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_kadap, tmp_kadap, this);
  double tmp_k2 = get_k2();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_k2, tmp_k2, this);
  double tmp_k1 = get_k1();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_k1, tmp_k1, this);
  double tmp_A1 = get_A1();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_A1, tmp_A1, this);
  double tmp_A2 = get_A2();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_A2, tmp_A2, this);
  double tmp_E_rev1 = get_E_rev1();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_E_rev1, tmp_E_rev1, this);
  double tmp_E_rev2 = get_E_rev2();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_E_rev2, tmp_E_rev2, this);
  double tmp_E_rev3 = get_E_rev3();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_E_rev3, tmp_E_rev3, this);
  double tmp_E_rev4 = get_E_rev4();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_E_rev4, tmp_E_rev4, this);
  double tmp_tau_syn1 = get_tau_syn1();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_tau_syn1, tmp_tau_syn1, this);
  double tmp_tau_syn2 = get_tau_syn2();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_tau_syn2, tmp_tau_syn2, this);
  double tmp_tau_syn3 = get_tau_syn3();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_tau_syn3, tmp_tau_syn3, this);
  double tmp_tau_syn4 = get_tau_syn4();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_tau_syn4, tmp_tau_syn4, this);
  double tmp_offset = get_offset();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_offset, tmp_offset, this);
  // initial values for state variables in ODE or kernel
  double tmp_V_m = get_V_m();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_V_m, tmp_V_m, this);
  double tmp_I_adap = get_I_adap();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_I_adap, tmp_I_adap, this);
  double tmp_I_dep = get_I_dep();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_I_dep, tmp_I_dep, this);
  long tmp_r = get_r();
  nest::updateValueParam<long>(__d, eglif_dcnp_names::_r, tmp_r, this);
  double tmp_lambda = get_lambda();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_lambda, tmp_lambda, this);
  double tmp_cf_buffer = get_cf_buffer();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_cf_buffer, tmp_cf_buffer, this);
  double tmp_gr_buffer = get_gr_buffer();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_gr_buffer, tmp_gr_buffer, this);
  double tmp_last_io = get_last_io();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_last_io, tmp_last_io, this);
  double tmp_tick = get_tick();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_tick, tmp_tick, this);
  double tmp_complex_flag = get_complex_flag();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_complex_flag, tmp_complex_flag, this);
  double tmp_g4__X__rec4 = get_g4__X__rec4();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_g4__X__rec4, tmp_g4__X__rec4, this);
  double tmp_g4__X__rec4__d = get_g4__X__rec4__d();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_g4__X__rec4__d, tmp_g4__X__rec4__d, this);
  double tmp_g2__X__rec2 = get_g2__X__rec2();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_g2__X__rec2, tmp_g2__X__rec2, this);
  double tmp_g2__X__rec2__d = get_g2__X__rec2__d();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_g2__X__rec2__d, tmp_g2__X__rec2__d, this);
  double tmp_g1__X__rec1 = get_g1__X__rec1();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_g1__X__rec1, tmp_g1__X__rec1, this);
  double tmp_g1__X__rec1__d = get_g1__X__rec1__d();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_g1__X__rec1__d, tmp_g1__X__rec1__d, this);
  double tmp_g3__X__rec3 = get_g3__X__rec3();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_g3__X__rec3, tmp_g3__X__rec3, this);
  double tmp_g3__X__rec3__d = get_g3__X__rec3__d();
  nest::updateValueParam<double>(__d, eglif_dcnp_names::_g3__X__rec3__d, tmp_g3__X__rec3__d, this);

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  ExtendedPostHistoryArchivingNode::set_status(__d);

  // if we get here, temporaries contain consistent set of properties
  set_C_m(tmp_C_m);
  set_tau_m(tmp_tau_m);
  set_E_L(tmp_E_L);
  set_t_ref(tmp_t_ref);
  set_V_reset(tmp_V_reset);
  set_V_th(tmp_V_th);
  set_Vmin(tmp_Vmin);
  set_I_e(tmp_I_e);
  set_Vinit(tmp_Vinit);
  set_lambda_0(tmp_lambda_0);
  set_tau_V(tmp_tau_V);
  set_kadap(tmp_kadap);
  set_k2(tmp_k2);
  set_k1(tmp_k1);
  set_A1(tmp_A1);
  set_A2(tmp_A2);
  set_E_rev1(tmp_E_rev1);
  set_E_rev2(tmp_E_rev2);
  set_E_rev3(tmp_E_rev3);
  set_E_rev4(tmp_E_rev4);
  set_tau_syn1(tmp_tau_syn1);
  set_tau_syn2(tmp_tau_syn2);
  set_tau_syn3(tmp_tau_syn3);
  set_tau_syn4(tmp_tau_syn4);
  set_V_m(tmp_V_m);
  set_I_adap(tmp_I_adap);
  set_I_dep(tmp_I_dep);
  set_r(tmp_r);
  set_lambda(tmp_lambda);
  set_offset(tmp_offset);
  set_cf_buffer(tmp_cf_buffer);
  set_gr_buffer(tmp_gr_buffer);
  set_last_io(tmp_last_io);
  set_tick(tmp_tick);
  set_complex_flag(tmp_complex_flag);
  set_g4__X__rec4(tmp_g4__X__rec4);
  set_g4__X__rec4__d(tmp_g4__X__rec4__d);
  set_g2__X__rec2(tmp_g2__X__rec2);
  set_g2__X__rec2__d(tmp_g2__X__rec2__d);
  set_g1__X__rec1(tmp_g1__X__rec1);
  set_g1__X__rec1__d(tmp_g1__X__rec1__d);
  set_g3__X__rec3(tmp_g3__X__rec3);
  set_g3__X__rec3__d(tmp_g3__X__rec3__d);




  updateValue< double >(__d, nest::names::gsl_error_tol, P_.__gsl_error_tol);
  if ( P_.__gsl_error_tol <= 0. )
  {
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }

  // recompute internal variables in case they are dependent on parameters or state that might have been updated in this call to set_status()
  recompute_internal_variables();
};



#endif /* #ifndef EGLIF_DCNP */
