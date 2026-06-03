.. include-nestml:: ../../../../cerebellar_models/nest_models/eglif_multirec.nestml
    :end-before: References

.. warning::
   The model described here is not matching the other LIF based models because of the
   sign in the membrane potential equations: the leak current should drive the membrane potential
   towards the resting state and not the opposite.

Neuron parameters
+++++++++++++++++

Spike Stochasticity
-------------------
The generation of action potentials can be configured to behave either deterministically
or stochastically through the boolean parameter ``stochastic_spiking`` (default: ``false``).

* **Deterministic Mode** (``stochastic_spiking = false``): A spike is generated
  subitaneously whenever the membrane potential crosses the fixed threshold value:

  .. math::
     V_m \ge V_{th}

* **Stochastic Mode** (``stochastic_spiking = true``): The neuron produces spikes
  stochastically according to a point process. The instantaneous firing intensity
  (escape rate) :math:`\lambda` depends exponentially on the distance between the
  current membrane potential and the threshold:

  .. math::
     \lambda = \lambda_0 \cdot e^{\dfrac{V_m - V_{th}}{\tau_V}}

  Where :math:`\lambda_0` is the base escape rate and :math:`\tau_V` is the voltage
  scale parameter. The probability of emitting a spike within a simulation time step
  :math:`\Delta t` is given by:

  .. math::
     P(\text{spike}) = 1 - e^{-\lambda \cdot \Delta t}

`In-vitro` state
----------------
The automatic tuning procedure was applied to the parameters of the E-GLIF model
(:math:`A_1`, :math:`A_2`, :math:`k_1`, :math:`k_2`, :math:`I_e`, :math:`k_{adap}`), yielding the cell-specific
values. The rest of the biophysical parameters were defined in accordance with the electrophysiological studies for each
neuronal population [#geminiani_2019]_.

For each cell type, the goal was to preserve the main electrophysiological features of the
corresponding multicompartmental reference model within the reduced point-neuron description.
Reference responses were generated through depolarizing and hyperpolarizing current-injection protocols
in NEURON. The same protocols were then simulated with the E-GLIF model, and electrophysiological
descriptors were extracted from both model responses.

These descriptors were used to define cell-specific objective (fitness) functions, yielding a
multi-objective minimization problem that captured the :math:`f-I` relationship together with
relevant dynamical features such as:

* Rheobase
* Tonic firing rate
* Coefficient of Variation (CV)
* :math:`f-I` response
* :math:`f-I` curve
* Post-inhibitory pauses (depending on the cell type)
* Rebound responses

See Table 1 for details in [#degrazia_2026]_.

The parameter search was implemented in a custom Python workflow based on the DEAP evolutionary
framework, using the NSGA-II algorithm as a multi-objective optimizer. This approach was chosen to
handle the simultaneous calibration of multiple electrophysiological features over a nonlinear
parameter space, a setting in which evolutionary methods have been widely used for neuronal model tuning.

Since NSGA-II returns a Pareto set of non-dominated solutions, the final parameter set for each cell
type was selected using an achievement scalarizing function, identifying a balanced trade-off among
the different objectives [#degrazia_2026]_.

The resulting parameters for all optimized cells are summarized in Table 2 and fitness errors are reported in
Table 3 of the reference paper [#degrazia_2026]_.

Synaptic parameters
+++++++++++++++++++

Static Synapse
--------------

`In-vitro` state
^^^^^^^^^^^^^^^^

For the static synapse configuration, synaptic transmission integrates receptor-specific kinetics
for AMPA, NMDA, and GABA receptors. To capture multi-timescale synaptic profiles, the total
conductance :math:`g_{X}(t)` for each receptor is modeled as the sum of a fast and a slow
exponential contribution:

.. math::
   g_{X}(t) = g_{X,\text{fast}}(t) + g_{X,\text{slow}}(t)

where :math:`X \in \{\text{AMPA}, \text{NMDA}, \text{GABA}\}`.

The underlying kernel parameters regulating these kinetic profiles, i.e., specifically the initial
conductance :math:`g_{\text{init}}`, the rise scaling factor :math:`A_r`, the decay scaling factors
(:math:`A_{d1}`, :math:`A_{d2}`), and the respective time constants (:math:`\tau_r`, :math:`\tau_{d1}`,
:math:`\tau_{d2}`) were calibrated using a multi-objective genetic algorithm (NSGA-II).

For the mathematical details underlying AMPA, NMDA, GABA synaptic kernels, see Equations 4–7 in [#degrazia_2026]_.

The optimization was targeted against reference synaptic conductance traces extracted from
multicompartmental NEURON simulations under voltage-clamp protocols. The resulting time
constants and scaling factors for each connection and postsynaptic cell type are illustrated
in Figure 5 of the reference paper [#degrazia_2026]_.

Furthermore, to compensate for the loss of dendritic attenuation inherent to the reduction from
spatially extended structures to point-neuron representations, a network-level weight calibration
was performed. The final effective static synaptic weights optimized for each pathway are summarized in Table 4 of [#degrazia_2026]_.

Tsodyks Markram Synapse
-----------------------

Short-term plasticity (STP) is incorporated into the network using the Tsodyks-Markram
formalism, with parameters assigned according to each specific connection type.
Circuits based on the De Grazia et al. model [#degrazia_2026]_ leverage the
`tsodyks_synapse <https://nest-simulator.readthedocs.io/en/latest/models/tsodyks_synapse.html>`_
version of the model.

`In-vitro` state
^^^^^^^^^^^^^^^^

The connection-specific parameters regulating these dynamic synapses
are defined in accordance with S1 Table of the reference paper [#degrazia_2026]_.

To preserve consistency between static and plastic network simulations, an automated
synaptic weight rescaling procedure was applied. A rescaling factor relative to the corresponding
static weight was computed for each connection and receptor type, ensuring that the area under
the synaptic conductance curve matches that of the static model under low-frequency stimulation (5 Hz),
where STP effects are assumed to be negligible [#degrazia_2026]_.

The final calibrated effective STP synaptic weights resulting from this connection-specific
rescaling are reported in Table 4 of the reference paper [#degrazia_2026]_.


References
++++++++++

.. include-nestml:: ../../../../cerebellar_models/nest_models/eglif_multirec.nestml
    :start-after: start-references
    :end-before: See also

.. [#geminiani_2019] Geminiani, A., Casellato, C., D’Angelo, E., & Pedrocchi, A. (2019).
   Complex electroresponsive dynamics in olivocerebellar neurons represented with extended-generalized
   leaky integrate and fire models. Frontiers in Computational Neuroscience, 13, 35.
   https://doi.org/10.3389/fncom.2019.00035