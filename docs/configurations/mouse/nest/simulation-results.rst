Nest Simulation Results
~~~~~~~~~~~~~~~~~~~~~~~

We present here the results of the `canonical circuit` simulations with the different protocol presented in :doc:`nest`
(in-vitro and awake states).
Neurons are here represented as `eglif_cond_alpha_multisyn` and are connected with `static synapses`.

Basal activity
##############

The results are the same with or without the Deep cerebellar nuclei regions as their activity does not effect
the rest of the circuit.
No basal activity changes are observed in the cerebellar network with Inferior Olive cells, because these cells
presents no auto-rhythm [#de_gruijl_2012]_ [#lefler_2013]_.
We also include the results of the Unipolar brush cell simulations (only in-vitro), as their effect on the other
populations is limited.

For this simulation, the mean firing rates and mean ISI obtained for each neuron population are as follows (expressed
as mean :math:`\pm` standard deviation):

`In-vitro` state
++++++++++++++++

.. csv-table:: Results of the canonical circuit in basal activity
   :header-rows: 1
   :delim: ;

    Cell name;Mean Firing rate (Hz); Mean ISI (ms)
    Mossy cell; :math:`4.0 \pm 1.4`; :math:`250 \pm 140`
    Granule cell; :math:`3.4 \pm 3.4`; :math:`330 \pm 270`
    Golgi cell;:math:`11 \pm 5.4`; :math:`120 \pm 94`
    Purkinje cell;:math:`47 \pm 1.1`; :math:`22 \pm 0.48`
    Basket cell;:math:`16 \pm 7.7`; :math:`77 \pm 45`
    Stellate cell;:math:`13 \pm 11`; :math:`180 \pm 240`
    DCNp; :math:`11 \pm 7.0`; :math:`97 \pm 110`
    DCNi; :math:`11 \pm 1.2`; :math:`90 \pm 11`
    Unipolar brush cell; :math:`0.21 \pm 0`; :math:`1100 \pm 91`

Awake state
+++++++++++

.. csv-table:: Results of the canonical circuit in basal activity
   :header-rows: 1
   :delim: ;

    Cell name;Mean Firing rate (Hz); Mean ISI (ms)
    Mossy cell; :math:`4.0 \pm 1.4`; :math:`250 \pm 140`
    Granule cell; :math:`3.4 \pm 3.4`; :math:`330 \pm 270`
    Golgi cell;:math:`11 \pm 5.4`; :math:`120 \pm 94`
    Purkinje cell;:math:`91 \pm 2.2`; :math:`11 \pm 0.27`
    Basket cell;:math:`21 \pm 9.0`; :math:`56 \pm 25`
    Stellate cell;:math:`18 \pm 14`; :math:`120 \pm 140`
    DCNp; :math:`45 \pm 7.0`; :math:`23 \pm 1.8`
    DCNi; :math:`13 \pm 1.2`; :math:`81 \pm 4.1`


Mossy fiber stimulus
####################

Similar to the ``basal activity``, we include the spiking results of the UBC (only in-vitro) and DCN populations
as their activity has little to no effect on the rest of the circuit. IO cells are here too silent.

For this simulation, **during the stimulus**, the mean firing rates and mean ISI obtained for each
neuron population are as follows (expressed in mean :math:`\pm` standard deviation):

`In-vitro` state
++++++++++++++++

.. csv-table:: Results of the canonical circuit during stimulus of the mossy
   :header-rows: 1
   :delim: ;

    Cell name;Mean Firing rate (Hz); Mean ISI (ms)
    Mossy cell; :math:`45 \pm 74`; :math:`6.5 \pm 3.`
    Granule cell; :math:`21 \pm 46`; :math:`9.1 \pm 7.3`
    Golgi cell;:math:`52 \pm 45`; :math:`17 \pm 12`
    Purkinje cell;:math:`59 \pm 11`; :math:`17 \pm 3.7`
    Basket cell;:math:`59 \pm 88`; :math:`12 \pm 7.3`
    Stellate cell;:math:`52 \pm 47`; :math:`14 \pm 7.5`
    DCNp; :math:`24 \pm 11`; :math:`30.0 \pm 4.2`
    DCNi; :math:`9.4 \pm 10`; not enough spikes per neuron
    Unipolar brush cell; :math:`16 \pm 7.9`; not enough spikes per neuron

Awake state
+++++++++++

.. csv-table:: Results of the canonical circuit during stimulus of the mossy
   :header-rows: 1
   :delim: ;

    Cell name;Mean Firing rate (Hz); Mean ISI (ms)
    Mossy cell; :math:`45 \pm 74`; :math:`6.5 \pm 3.`
    Granule cell; :math:`21 \pm 46`; :math:`9.1 \pm 7.3`
    Golgi cell;:math:`52 \pm 45`; :math:`17 \pm 12`
    Purkinje cell;:math:`140 \pm 33`; :math:`7.1 \pm 1.9`
    Basket cell;:math:`77 \pm 54`; :math:`8.7 \pm 5.9`
    Stellate cell;:math:`71 \pm 55`; :math:`9.6 \pm 6.2`
    DCNp; :math:`57 \pm 16`; :math:`20 \pm 6.2`
    DCNi; :math:`10 \pm 10`; not enough spikes per neuron


Eyeblink Classical Conditioning
###############################

`In-vitro` state
++++++++++++++++

.. csv-table:: Results of the canonical circuit in `in-vitro` state with DCN and IO during stimulus of the mf and the IO
   :header-rows: 1
   :delim: ;

    Cell name;Mean Firing rate (Hz) [CS]; Mean ISI (ms) [CS]; Mean Firing rate (Hz) [CS+US]; Mean ISI (ms) [CS+US]
    Mossy cell; :math:`44 \pm 13`; :math:`22 \pm 8.0`; :math:`38 \pm 60`; :math:`3.1 \pm 2.3`
    Granule cell; :math:`21 \pm 30`; :math:`36 \pm 36`; :math:`20 \pm 48`; :math:`6.1 \pm 1.1`
    Golgi cell;:math:`42 \pm 17`; :math:`28 \pm 13`; :math:`26 \pm 50`; :math:`8.7 \pm 0.2`
    Purkinje cell;:math:`55 \pm 5.9`; :math:`18.0 \pm 2.0`; :math:`140 \pm 100`; :math:`3.1 \pm 2.8`
    Basket cell;:math:`66 \pm 24`; :math:`19 \pm 9.5`; :math:`61 \pm 49`; not enough spikes per neuron
    Stellate cell;:math:`46 \pm 43`; :math:`31 \pm 31`; :math:`48 \pm 57`; :math:`7.2 \pm 1.1`
    DCNp; :math:`30 \pm 4.7`; :math:`33 \pm 5.2`; :math:`33 \pm 47`; not enough spikes per neuron
    DCNi; :math:`11 \pm 2.3`; :math:`93 \pm 11`;  :math:`12 \pm 32`; not enough spikes per neuron
    IO; 0; no spikes; :math:`270 \pm 120`; :math:`1.9 \pm 0.5`

Awake state
+++++++++++

.. csv-table:: Results of the canonical circuit in awake state with DCN and IO during stimulus of the mf and the IO
   :header-rows: 1
   :delim: ;

    Cell name;Mean Firing rate (Hz) [CS]; Mean ISI (ms) [CS]; Mean Firing rate (Hz) [CS+US]; Mean ISI (ms) [CS+US]
    Mossy cell; :math:`43 \pm 13`; :math:`23 \pm 8.1`; :math:`46 \pm 68`; :math:`2.8 \pm 1.8`
    Granule cell; :math:`22 \pm 35`; :math:`34 \pm 37`; :math:`22 \pm 50`; :math:`6.4 \pm 1.2`
    Golgi cell;:math:`42 \pm 21`; :math:`29 \pm 15`; :math:`22 \pm 45`; :math:`7.0 \pm 0.0`
    Purkinje cell;:math:`140 \pm 14`; :math:`7.1 \pm 0.71`; :math:`300 \pm 150`; :math:`2.7 \pm 1.9`
    Basket cell;:math:`83 \pm 49`; :math:`19 \pm 21`; :math:`55 \pm 54`; :math:`6.3 \pm 1.1`
    Stellate cell;:math:`74 \pm 68`; :math:`22 \pm 24`; :math:`57 \pm 81`; :math:`5.6 \pm 1.3`
    DCNp; :math:`59 \pm 7.6`; :math:`17 \pm 2.1`; :math:`65 \pm 56`; :math:`6.0 \pm 2.7`
    DCNi; :math:`11 \pm 2.0`; :math:`90 \pm 13`; :math:`9.1 \pm 29`; not enough spikes per neuron
    IO; 0; no spikes; :math:`350 \pm 150`; :math:`1.9 \pm 0.78`

References
##########

.. [#de_gruijl_2012] De Gruijl, J. R., Bazzigaluppi, P., de Jeu, M. T., & De Zeeuw, C. I. (2012).
   "Climbing fiber burst size and olivary sub-threshold oscillations in a network setting."
   PLoS computational biology, 8(12), e1002814.
   https://doi.org/10.1371/journal.pcbi.1002814
.. [#lefler_2013] Lefler, Y., Torben-Nielsen, B., & Yarom, Y. (2013).
   "Oscillatory activity, phase differences, and phase resetting in the inferior olivary nucleus."
   Frontiers in systems neuroscience, 7, 22.
   https://doi.org/10.3389/fnsys.2013.00022
