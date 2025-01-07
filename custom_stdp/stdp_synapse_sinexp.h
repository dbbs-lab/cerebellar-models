/*
 *  stdp_synapse.h
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
 */

#ifndef STDP_SYNAPSE_SINEXP_H
#define STDP_SYNAPSE_SINEXP_H

// C++ includes:
#include <cmath>

// Includes from nestkernel:
#include "common_synapse_properties.h"
#include "connection.h"
#include "connector_model.h"
#include "event.h"

// Includes from sli:
#include "dictdatum.h"
#include "dictutils.h"

#include "histentry.h"
#include "archiving_node.h"
#include "extended_hist_entry.h"
#include "extended_post_history_archiving_node.h"

namespace nest
{

/* BeginUserDocs: synapse, spike-timing-dependent plasticity

Short description
+++++++++++++++++

Synapse type for spike-timing dependent plasticity

Description
+++++++++++

``stdp_synapse`` is a connector to create synapses with spike time
dependent plasticity (as defined in [1]_). Here the weight dependence
exponent can be set separately for potentiation and depression.

.. warning::

   This synaptic plasticity rule does not take
   :ref:`precise spike timing <sim_precise_spike_times>` into
   account. When calculating the weight update, the precise spike time part
   of the timestamp is ignored.

See also [2]_, [3]_, [4]_.

Parameters
++++++++++

========= =======  ======================================================
 tau_plus  ms      Time constant of STDP window, potentiation
                   (tau_minus defined in postsynaptic neuron)
 lambda    real    Step size
 alpha     real    Asymmetry parameter (scales depressing increments as
                   alpha*lambda)
 mu_plus   real    Weight dependence exponent, potentiation
 mu_minus  real    Weight dependence exponent, depression
 Wmax      real    Maximum allowed weight
========= =======  ======================================================

Transmits
+++++++++

SpikeEvent

References
++++++++++

.. [1] Guetig et al. (2003). Learning input correlations through nonlinear
       temporally asymmetric hebbian plasticity. Journal of Neuroscience,
       23:3697-3714 DOI: https://doi.org/10.1523/JNEUROSCI.23-09-03697.2003
.. [2] Rubin J, Lee D, Sompolinsky H (2001). Equilibrium
       properties of temporally asymmetric Hebbian plasticity. Physical Review
       Letters, 86:364-367. DOI: https://doi.org/10.1103/PhysRevLett.86.364
.. [3] Song S, Miller KD, Abbott LF (2000). Competitive Hebbian learning
       through spike-timing-dependent synaptic plasticity. Nature Neuroscience
       3(9):919-926.
       DOI: https://doi.org/10.1038/78829
.. [4] van Rossum MCW, Bi G-Q, Turrigiano GG (2000). Stable Hebbian learning
       from spike timing-dependent plasticity. Journal of Neuroscience,
       20(23):8812-8821.
       DOI: https://doi.org/10.1523/JNEUROSCI.20-23-08812.2000

See also
++++++++

tsodyks_synapse, static_synapse

Examples using this model
+++++++++++++++++++++++++

.. listexamples:: stdp_synapse

EndUserDocs */

// connections are templates of target identifier type (used for pointer /
// target index addressing) derived from generic connection template

void register_stdp_synapse_sinexp( const std::string& name );

template < typename targetidentifierT >
class stdp_synapse_sinexp : public Connection< targetidentifierT >
{

public:
  typedef CommonSynapseProperties CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;

  static constexpr ConnectionModelProperties properties = ConnectionModelProperties::HAS_DELAY
    | ConnectionModelProperties::IS_PRIMARY | ConnectionModelProperties::SUPPORTS_HPC
    | ConnectionModelProperties::SUPPORTS_LBL;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  stdp_synapse_sinexp();


  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  stdp_synapse_sinexp( const stdp_synapse_sinexp& ) = default;
  stdp_synapse_sinexp& operator=( const stdp_synapse_sinexp& ) = default;

  // Explicitly declare all methods inherited from the dependent base
  // ConnectionBase. This avoids explicit name prefixes in all places these
  // functions are used. Since ConnectionBase depends on the template parameter,
  // they are not automatically found in the base class.
  using ConnectionBase::get_delay;
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;

  /**
   * Get all properties of this connection and put them into a dictionary.
   */
  void get_status( DictionaryDatum& d ) const;

  /**
   * Set properties of this connection from the values given in dictionary.
   */
  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param cp common properties of all synapses (empty).
   */
  bool send( Event& e, size_t t, const CommonSynapseProperties& cp );


  class ConnTestDummyNode : public ConnTestDummyNodeBase
  {
  public:
    // Ensure proper overriding of overloaded virtual functions.
    // Return values from functions are ignored.
    using ConnTestDummyNodeBase::handles_test_event;
    size_t
    handles_test_event( SpikeEvent&, size_t ) override
    {
      return invalid_port;
    }
  };

  void
  check_connection( Node& s, Node& t, size_t receptor_type, const CommonPropertiesType& )
  {
    ConnTestDummyNode dummy_target;

    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );

    t.register_stdp_connection( t_lastspike_ - get_delay(), get_delay() );
  }

  void
  set_weight( double w )
  {
    weight_ = w;
  }
  std::vector<double> simple_post_spikes_;

private:
  double
  facilitate_()
  {
    double k = Aplus_;
    return k;
  }

  double
  depress_(double tempo)
  {
    double k = (1/C_) * std::exp((-tempo+t0_)/tau_) * std::pow(sin(2*3.14159*(tempo-t0_)/tau_), 20);
    return k * Aminus_;
  }

  // data members of each connection
  double weight_;
  double Wmin_;
  double Wmax_;
  double tau_;
  double C_;
  double t0_;
  double Aplus_;
  double Aminus_;
  

  double t_lastspike_;
};

template < typename targetidentifierT >
constexpr ConnectionModelProperties stdp_synapse_sinexp< targetidentifierT >::properties;

/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param t The thread on which this connection is stored.
 * \param cp Common properties object, containing the stdp parameters.
 */
template < typename targetidentifierT >
inline bool
stdp_synapse_sinexp< targetidentifierT >::send( Event& e, size_t t, const CommonSynapseProperties& )
{

  // synapse STDP depressing/facilitation dynamics
  const double t_spike = e.get_stamp().get_ms();

  // use accessor functions (inherited from Connection< >) to obtain delay and
  // target
  ExtendedPostHistoryArchivingNode* target =dynamic_cast<ExtendedPostHistoryArchivingNode*>(get_target( t ));
  double dendritic_delay = get_delay();

  // get spike history in relevant range (t1, t2] from postsynaptic neuron
  std::deque< ExtendedHistEntry >::iterator start;
  std::deque< ExtendedHistEntry >::iterator finish;

  // For a new synapse, t_lastspike_ contains the point in time of the last
  // spike. So we initially read the
  // history(t_last_spike - dendritic_delay, ..., T_spike-dendritic_delay]
  // which increases the access counter for these entries.
  // At registration, all entries' access counters of
  // history[0, ..., t_last_spike - dendritic_delay] have been
  // incremented by ArchivingNode::register_stdp_connection(). See bug #218 for
  // details.
  target->get_history( t_lastspike_ - dendritic_delay, t_spike - dendritic_delay, &start, &finish );
  // facilitation due to postsynaptic spikes since last pre-synaptic spike
  double minus_dt;
  //double post = 0;
  double LTD = 0;
  //double last_post_spike = 0;
  while ( start != finish )
  {
    minus_dt = t_lastspike_ - ( start->t_ + dendritic_delay );
    const double offset = start->offset_;
    // get_history() should make sure that
    // start->t_ > t_lastspike - dendritic_delay, i.e. minus_dt < 0
    std::cout << "Custom STDP synapse: processing post spike with offset = " << offset << std::endl;
    std::cout << "VALUE OF MINUS_DT = " << minus_dt << std::endl;
    //simple_post_spikes_.push_back(t_lastspike_);
    //assert( minus_dt < -1.0 * kernel().connection_manager.get_stdp_eps() );
  //  if (offset != 0)
  //  {
    // le righe di codice commentate sono obsolete !!!
    //  for (unsigned int GR=0; GR < simple_post_spikes_.size(); GR++) {
    //    double sd = simple_post_spikes_[GR] - minus_dt;
        //std::cout << "sddddddddddddd: " << sd;
        //double sd_minus = simple_post_spikes_[GR-1] - minus_dt;
    //    if (sd < 0 && sd > -200 && post != 0){
    if (offset != 0 && minus_dt <= 0 && minus_dt >= -200){
      LTD += depress_(minus_dt);
    }
    //    }
    //    else if (sd < 0 && sd > -200 && post == 0){
    //      LTD += 2*depress_(sd);
          //LTD += depress_(sd_minus);
    //      ++post;
    //    }
    //    last_post_spike = sd;
    //  }
    weight_ = weight_ + LTD;
    if (weight_ <= Wmin_)
    {
      weight_ = Wmin_;
    }

    //  if (last_post_spike <= -200){
    //    post = 0;
    //  }
  //  }
    ++start;
  }


  weight_ = weight_ + facilitate_();
//  std::cout << "t_lastspike STA QUIIIIIII = " << t_lastspike_ << std::endl;
  if (weight_ >= Wmax_){
    weight_ = Wmax_;
  }
  //  while(simple_post_spikes_[0] < t_spike - 200){
  //    simple_post_spikes_.erase(simple_post_spikes_.begin());
  //  }
  //  std::cout << "vector: " << simple_post_spikes_[0] << std::endl;

  e.set_receiver( *target );
  e.set_weight( weight_ );
    // use accessor functions (inherited from Connection< >) to obtain delay in
    // steps and rport
  e.set_delay_steps( get_delay_steps() );
  e.set_rport( get_rport() );
  e();
  t_lastspike_ = t_spike;
//  }

  return true;
}


template < typename targetidentifierT >
stdp_synapse_sinexp< targetidentifierT >::stdp_synapse_sinexp()
  : ConnectionBase()
  , weight_( 0.14 )
  , Wmin_( 0.0 )
  , Wmax_( 0.6 )
  , tau_(1000)
  , C_( 1.2848 )
  , t0_( 100 )
  , Aplus_( 0.00005 )
  , Aminus_( -0.005 )

  , t_lastspike_( 0.0 )
{
}

template < typename targetidentifierT >
void
stdp_synapse_sinexp< targetidentifierT >::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );
  def< double >( d, names::Wmin, Wmin_ );
  def< double >( d, names::Wmax, Wmax_ );
  def< double >( d, names::tau, tau_ );
  def< double >( d, names::Aplus, Aplus_ );
  def< double >( d, names::Aminus, Aminus_ );
  def< long >( d, names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
stdp_synapse_sinexp< targetidentifierT >::set_status( const DictionaryDatum& d, ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
  updateValue< double >( d, names::Wmin, Wmin_ );
  updateValue< double >( d, names::Wmax, Wmax_ );
  updateValue< double >( d, names::tau, tau_ );
  updateValue< double >( d, names::Aplus, Aplus_ );
  updateValue< double >( d, names::Aminus, Aminus_ );

  // check if weight_ and Wmax_ has the same sign
  if ( not( ( ( weight_ >= 0 ) - ( weight_ < 0 ) ) == ( ( Wmax_ >= 0 ) - ( Wmax_ < 0 ) ) ) )
  {
    throw BadProperty( "Weight and Wmax must have same sign." );
  }

  if ( Aplus_ < 0 )
  {
    throw BadProperty( "Aplus must be non-negative." );
  }
  
  if ( Aminus_ > 0 )
  {
    throw BadProperty( "Aminus must be non-positive." );
  }
}

} // of namespace nest

#endif // of #ifndef STDP_SYNAPSE_SINEXP_H
