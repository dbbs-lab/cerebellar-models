.. include-nestml:: ../../../../cerebellar_models/nest_models/eglif_cond_alpha_multisyn.nestml
    :end-before: References

.. warning::
   The model described here is not matching the other LIF based models because of the
   sign in the membrane potential equations: the leak current should drive the membrane potential
   towards the resting state and not the opposite.

   We have not being able to reproduce all the results of the Geminiani et al. (2018 and 2019)
   [#geminiani_2018]_ [#geminiani_2019]_, so we advise you to be careful when using it.

Neuron parameters
+++++++++++++++++

`In-vitro` state
----------------

The parameters for the EGLIF models were extracted from Table 2 and Table 3 in Geminiani et al.
(2019) [#geminiani_2019]_. We did the following modifications to match better the paper's results:

- The :math:`k_2` parameter should be greater than :math:`\dfrac{1}{\tau_m}` to prevent unstable
  oscillations of the membrane potential (see proof of Geminiani et al. [#geminiani_2018]_)
  but the authors seemed to have rounded down the values which resulted in an unstable behavior
  for the GrC. We therefore rounded up this value, since the GrC should not produce spikes without stimulus.
- For the PC, we modified also the :math:`I_e` value so that the tonic firing rate of PC is ~45 Hz
  [#telgkamp_2002]_ but maintained the F/I curve slope from the paper.
- The parameters of DCNp and DCNi populations were slightly changed so that
  their tonic firing rates remains around ~10Hz [#moscato_2019]_.

.. warning::
   - It is not clear how the spiking parameters (i.e :math:`\lambda_0` and :math:`\tau_V` and initial :math:`V_m`)
     are obtained in the Geminiani et al. (2019) paper [#geminiani_2019]_ .
     These parameters were manually set to reproduce the F/I curves from the Figure 4 and Figure 3 from
     respectively Geminiani et al. (2018 and 2019) papers [#geminiani_2018]_ [#geminiani_2019]_.
   - On a side note, in the optimization section of the Geminiani et al. (2018) paper
     [#geminiani_2018]_, the authors wrote that the :math:`k_2` parameter should not be optimized but
     set to :math:`\dfrac{1}{\tau_m}` to have stable oscillations but this is not the case for most of
     the :math:`k_2` parameters listed in  Geminiani et al. (2019) paper [#geminiani_2019]_ .

The postsynaptic currents are integrated to the soma with alpha exponential functions. Each function
is defined with a reversal potential parameter :math:`E_{rev}` and a time constant :math:`\tau_{syn}`.
These parameters depend on the connection types. In NEST, they are defined in the neuron equations.

The postsynaptic receptor parameters are listed in Table 2 of Geminiani et al. (2019b)
[#geminiani_2019b]_ .

For the Unipolar Brush cells, the LIF parameters of the EGLIF model were extracted from
Locatelli et al. (2013) [#locatelli_2013]_, Subramaniyam et al. (2014) [#subramaniyam_2014]_
and Russo et al. (2007). The rest of the EGLIF parameters were optimized to match results of
Locatelli et al. (2013) [#locatelli_2013]_ using the Geminiani et al. (2018) [#geminiani_2018]_ method.

Awake state
-----------

The parameters for the awake state are the same as the `in-vitro` state because for most of the cells
we do not have the data to fit our model to. However, we changed for the following:

- The endogenous current :math:`I_e` of PC was set to 700 pA and :math:`\lambda_0` :math:`\tau_V` were
  changed to increase the F/I curve slope. We targeted here ~80 Hz of tonic firing rate to match the
  range of Table 1 from Geminiani et al. 2024 [#geminiani_2024]_.
- The spiking parameters (i.e :math:`\lambda_0` and :math:`\tau_V`) were tuned for DCNp to match the
  Geminiani et al. 2024 [#geminiani_2024]_.

Synaptic parameters
+++++++++++++++++++

Static Synapse
--------------

`In-vitro` state
^^^^^^^^^^^^^^^^

The synaptic parameters used for the `canonical circuit` corresponds to the one listed in Table B of
supplementary document 1 in Geminiani et al. (2024) [#geminiani_2024]_. The receptor id corresponds
to the postsynaptic receptor used for the connection.
However, we made the following adjustments to obtain results close to the paper and literature:

- In our experiments, we decreased the weights for the pf-SC, pf-BC so that the activity of MLI lies around
  ~15 Hz for both BC and SC [#kim_2021]_. Then aa-PC, pf-PC were decreased to maintain the PC in a stable
  low activity ~50Hz [#telgkamp_2002]_.
- Synaptic parameter values of DCNp, DCNi and IO were manually adjusted through trial and error to ensure a
  reasonable excitation/inhibition activity ratio.
- Finally, the SC-PC was scaled to take into account the increase of synapses from the connectivity rule.

.. warning::
   It is currently unclear from the paper, how the synaptic parameters were optimized, or which features were targeted.

Awake state
^^^^^^^^^^^

The parameters for the awake state are the same as the `in-vitro` state, except for the connections
GrC(pf)-PC, GrC(aa)-PC, GrC(pf)-SC, GrC(pf)-BC, BC-PC for which we adapted the weights to match the Table 1
from Geminiani et al. 2024 [#geminiani_2024]_.


Tsodyks Markram Synapse
-----------------------

Circuits based on the Geminiani et al. model [#geminiani_2018]_ leverage the
`tsodyks2_synapse <https://nest-simulator.readthedocs.io/en/latest/models/tsodyks2_synapse.html>`_
version of the model.

For each synapse of the `canonical circuit`, the initial value of ``u`` was set to ``U`` and ``x`` to ``1.0``.

`In-vitro` state
^^^^^^^^^^^^^^^^

The synaptic parameters used for the canonical circuit correspond to those
listed in the table below, obtained from Masoli et al. (2022) [#masoli_2022]_ .
The receptor ID corresponds to the postsynaptic receptor used for the connection.
The weights have been rescaled under the assumption
that the first peak of the postsynaptic conductance (:math:`g_{syn_0}`) for the
Tsodyks–Markram synapse must have the same amplitude as the ones obtained with a static_synapse model.

:math:`weight_{tsodyks} = \dfrac{{weight_{static}}^2}{g_{syn_0}}`

.. Note::
   The connections mf-glom and GoC-GoC are both considered static since, for these two connections,
   we do not have Tsodyks-Markram parameters.
   Moreover, for the pf-SC connection, the weight was adjusted manually to keep the firing rate
   within the desired range.


Awake state
^^^^^^^^^^^

The parameters for the awake state are the same as the in-vitro state, except for the following connections:
GrC(pf)-PC, GrC(aa)-PC, GrC(pf)-SC, GrC(pf)-BC, BC-PC for which we adapted the weights to match the Table 1
from Geminiani et al. 2024 [#geminiani_2024]_.

.. Note::
   For the simulations using Tsodyks-Markram synapse, the mean firing rates and mean interspike intervals (ISI)
   obtained for each neuron population from both in-vitro and awake states are expected to be the same,
   as the ones obtained with static synapses.
   For pf-SC connection weight was adjusted manually to keep the firing rate in the desired range.

References
++++++++++

.. include-nestml:: ../../../../cerebellar_models/nest_models/eglif_cond_alpha_multisyn.nestml
    :start-after: start-references
    :end-before: See also

.. [#geminiani_2019] Geminiani, A., Casellato, C., D’Angelo, E., & Pedrocchi, A. (2019).
   Complex electroresponsive dynamics in olivocerebellar neurons represented with extended-generalized
   leaky integrate and fire models. Frontiers in Computational Neuroscience, 13, 35.
   https://doi.org/10.3389/fncom.2019.00035
.. [#geminiani_2019b] Geminiani, A., Pedrocchi, A., D’Angelo, E., & Casellato, C. (2019). Response
   dynamics in an olivocerebellar spiking neural network with non-linear neuron properties.
   Frontiers in computational neuroscience, 13, 68.
   https://doi.org/10.3389/fncom.2019.00068
.. [#geminiani_2024] Geminiani, A., Casellato, C., Boele, H. J., Pedrocchi, A., De Zeeuw, C. I., &
   D’Angelo, E. (2024). Mesoscale simulations predict the role of synergistic cerebellar plasticity
   during classical eyeblink conditioning. PLOS Computational Biology, 20(4), e1011277.
   https://doi.org/10.1371/journal.pcbi.1011277
.. [#moscato_2019] Moscato, L., Montagna, I., De Propris, L., Tritto, S., Mapelli, L., & D’Angelo, E. (2019).
   Long-lasting response changes in deep cerebellar nuclei in vivo correlate with low-frequency oscillations.
   Frontiers in Cellular Neuroscience, 13, 433625.
   https://doi.org/10.3389/fncel.2019.00084
.. [#locatelli_2013] Locatelli, F., Bottà, L., Prestori, F., Masetto, S., & D’Angelo, E. (2013).
   "Late-onset bursts evoked by mossy fibre bundle stimulation in unipolar brush cells: Evidence for the involvement
   of H- and TRP-currents."
   The Journal of Physiology, 591(4), 899–918.
   https://doi.org/10.1113/jphysiol.2012.242180
.. [#subramaniyam_2014] Subramaniyam, S., Solinas, S., Perin, P., Locatelli, F., Masetto, S., & D’Angelo, E. (2014).
   "Computational modeling predicts the ionic mechanism of late-onset responses in unipolar brush cells."
   Frontiers in Cellular Neuroscience, 8.
   https://doi.org/10.3389/fncel.2014.00237
.. [#telgkamp_2002] Telgkamp, P., & Raman, I. M. (2002).
   Depression of inhibitory synaptic transmission between Purkinje cells and neurons of the cerebellar nuclei.
   Journal of Neuroscience, 22(19), 8447-8457.
   https://doi.org/10.1523/JNEUROSCI.22-19-08447.2002.
.. [#kim_2021] Kim, J., & Augustine, G. J. (2021).
   Molecular layer interneurons: key elements of cerebellar network computation and behavior.
   Neuroscience, 462, 22-35.
   https://doi.org/10.1016/j.neuroscience.2020.10.008
.. [#masoli_2022] Masoli, S., Rizza, M. F., Tognolina, M., Prestori, F., & D’Angelo, E. (2022).
   Computational models of neurotransmission at cerebellar synapses unveil the impact on network computation.
   Frontiers in Computational Neuroscience, 16, 1006989.
   https://doi.org/10.3389/fncom.2022.1006989
