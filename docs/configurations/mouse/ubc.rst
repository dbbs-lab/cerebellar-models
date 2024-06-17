Unipolar Brush Cells
~~~~~~~~~~~~~~~~~~~~

Unipolar Brush cells (UBC) are excitatory interneurons located in the granular layer that acts as relay cells for
glomerulus [#mugnaini_2011]_. The default configuration for UBCs is implemented in
`ubc.yaml <https://github.com/dbbs-lab/cerebellum/blob/master/configurations/mouse/ubc/ubc.yaml>`_.

Cellular density
^^^^^^^^^^^^^^^^

The density of UBC depends on the cerebellar cortex regions. Sekerková et al. (2014) [#sekerkova_2014]_ have reported
the following densities for mouse brain:

.. csv-table::
   :header-rows: 1
   :delim: ;

    Brain region name; Density (:math:`mm^{-3}`)
    Lingula (I); 6000
    Lobule II; 6000
    Central Lobule; 6000
    Culmen; 6000
    Declive (VI); 6000
    Folium-tuber vermis (VII); 6000
    Pyramus (VIII); 6000
    Uvula (IX); 25600
    Nodulus (X); 98400
    Simple lobule; 6000
    Ansiform lobule; 6000
    Paramedian lobule; 25600
    Copula pyramidis; 25600
    Flocculus; 78000

Remember that these values need to be converted into cell / :math:`\mu m^{-3}` to be used in BSB configurations.

Connectivity
^^^^^^^^^^^^

UBC are here represented with one dendrite (max. ``50`` :math:`\mu m` long), that connects to a single glomerulus and
produces ``4`` axons (max. ``500`` :math:`\mu m` long), each having one UBC glomerulus terminals.

.. note::
   UBC cells can have more than one dendrite but this subpopulation is not represented in the current circuit.
   It is also unclear if only a single UBC can connect to a glomerulus so this constraint is at the moment not
   implemented in the provided configurations.

Mf and UBC based glomeruli act similarly in the network: GrC, GoC and UBC can connect to them simultaneously.
UBC are also inhibited by Golgi cells that sends axons to their presynaptic glomerulus. One third of the UBC population
is reported to receive input from mf glomerulus while the rest receives input from UBC glomerulus. See review of
Mugnaini et al. (2011) [#mugnaini_2011]_ for more information. The connection id of the first column corresponds to
the numbers reported in :numref:`fig-network` and overrides the table :ref:`table-connectivity` of the canonical circuit.

.. csv-table::
   :header-rows: 1
   :delim: ;

   #; Source Name; Source Branch; Target Name; Target Branch; Strategy; References
   2; glom & ubc_glon ; /; GrC; dendrites; :ref:`glom_grc`; Houston et al. (2017) [#houston_2017]_
   3; glom & ubc_glon; /; GoC; basal dendrites; :ref:`glom_goc`; Kanichay and Silver (2008) [#kanichay_2008]_
   16; UBC ; /; ubc_glom ; /; :ref:`ubc_glom`; Mugnaini et al. (2011) [#mugnaini_2011]_
   17; glom & ubc_glon ; /; UBC ; /; :ref:`glom_ubc`; Mugnaini et al. (2011) [#mugnaini_2011]_
   18; Golgi ; axon; UBC ; /; :ref:`goc_glom`; Mugnaini et al. (2011) [#mugnaini_2011]_

NEST simulation
^^^^^^^^^^^^^^^

Neuron parameters
+++++++++++++++++

The UBC neurons were represented as a EGLIF point neuron models (see :doc:`NEST section <nest>`).

The following LIF parameters for the UBC cells were extracted from Locatelli et al. (2013) [#locatelli_2013]_,
Subramaniyam et al. (2014) [#subramaniyam_2014]_ and Russo et al. (2007):

.. csv-table:: UBC LIF neuron parameters
   :header-rows: 1
   :delim: ;

   :math:`C_m\ (pF)`;:math:`\tau_m\ (ms)`;:math:`E_L\ (mV)`;:math:`t_{ref}\ (ms)`;:math:`V_{reset}\ (mV)`;:math:`V_{th}\ (mV)`
   16.5; 13.2;-66.7;1.67;-76,7;-55.8

EGLIF parameters were optimized to match results of Locatelli et al. (2013) [#locatelli_2013]_ using the
Geminiani et al. (2018) [#geminiani_2018]_ method:

.. csv-table:: UBC EGLIF neuron parameters
   :header-rows: 1
   :delim: ;

    :math:`k_{adap}\ (nS \cdot ms^{-1})`;:math:`k_1\ (ms^{-1})`;:math:`k_2\ (ms^{-1})`;:math:`A_1\ (pA)`;:math:`A_2\ (pA)`;:math:`I_e\ (pA)`
    1.17; 0.14; 0.76; 0.0001; 0.0001; 1.0

It is not clear how the spiking parameters (i.e :math:`\lambda_0` and :math:`\tau_V`) are obtained in the Geminiani et
al. (2018) paper [#geminiani_2018]_. The values were extracted from a BSB configuration provided by the authors.

The postsynaptic receptors are defined according to the following table:

.. _ubc-table-receptor:
.. csv-table:: UBC Postsynaptic receptor parameters
   :header-rows: 1
   :delim: ;

   Receptor id; :math:`E_{rev,i}\ (mV)`; :math:`\tau_{syn,i}\ (ms)`; Type
   1; 0; 0.2; exc.
   2; -80; 2.0; inh.
   3; 0; 1.2; exc.

Synapse parameters
++++++++++++++++++

UBC connections are represented as ``static synapses`` (see :doc:`NEST section <nest>`).
The receptor id corresponds to the postsynaptic receptor used for the connection (see table :ref:`ubc-table-receptor`).
However, it is currently unclear how these parameters were optimized, or which features were targeted:

.. csv-table:: Presynaptic parameters
   :header-rows: 1
   :delim: ;

    Source-Target;:math:`weight \ (nS)`;:math:`delay \ (ms)`; Receptor id
    UBC-ubc_glom;1;1;1
    glom-GrC;0.23;1;1
    glom-GoC;0.24;1;1
    GoC-GrC;0.24;2;2

Simulation paradigms
++++++++++++++++++++

The `ubc_nest.yaml <https://github.com/dbbs-lab/cerebellum/blob/master/configurations/mouse/ubc/ubc_nest.yaml>`_ are
including all the simulation paradigms described in the :doc:`NEST section <nest>`) but include the UBC cells in the
circuit.

Basal activity
##############

For this simulation, the mean firing rates and mean ISI obtained for each neuron population are as
follows (expressed in mean :math:`\pm` standard deviation):

.. csv-table:: Results of the canonical circuit with UBC in basal activity
   :header-rows: 1
   :delim: ;

    Cell name;Mean Firing rate (Hz); Mean ISI (ms)
    Mossy cell; :math:`4.1 \pm 0.93`; :math:`240 \pm 61`
    Granule cell; :math:`2.9 \pm 2.1`; :math:`480 \pm 440`
    Unipolar brush cell; :math:`0.2 \pm 0`; not enough spikes per neuron
    Golgi cell;:math:`7.4 \pm 4.2`; :math:`220 \pm 270`
    Purkinje cell;:math:`48 \pm 1.6`; :math:`21 \pm 0.73`
    Basket cell;:math:`19 \pm 8.4`; :math:`70 \pm 47`
    Stellate cell;:math:`24 \pm 13`; :math:`80 \pm 130`

Mossy fiber stimulus
####################

For this simulation, **during the stimulus**, the mean firing rates and mean ISI obtained for each
neuron population are as follows (expressed in mean :math:`\pm` standard deviation):

.. csv-table:: Results of the canonical circuit with UBC during stimulus of the mossy
   :header-rows: 1
   :delim: ;

    Cell name;Mean Firing rate (Hz); Mean ISI (ms)
    Mossy cell; :math:`46 \pm 73`; :math:`6.3 \pm 2.8`
    Granule cell; :math:`19 \pm 36`; :math:`9.2 \pm 4.8`
    Unipolar brush cell; :math:`16 \pm 8.1`; not enough spikes per neuron
    Golgi cell;:math:`32 \pm 16`; :math:`7.2 \pm 3.5`
    Purkinje cell;:math:`59 \pm 15`; :math:`15 \pm 2.9`
    Basket cell;:math:`73 \pm 53`; :math:`8.4 \pm 4.1`
    Stellate cell;:math:`98 \pm 68`; :math:`7.9 \pm 5.3`

References
^^^^^^^^^^

.. [#mugnaini_2011] Mugnaini, E., Sekerková, G., & Martina, M. (2011). The unipolar brush cell: a remarkable neuron finally
   receiving deserved attention. Brain research reviews, 66(1-2), 220-245.
   https://doi.org/10.1016/j.brainresrev.2010.10.001
.. [#sekerkova_2014] Sekerková, G., Watanabe, M., Martina, M., & Mugnaini, E. (2014). Differential distribution of
   phospholipase C beta isoforms and diaglycerol kinase-beta in rodents cerebella corroborates the division of unipolar
   brush cells into two major subtypes. Brain Structure and Function, 219, 719-749.
   https://doi.org/10.1007/s00429-013-0531-9
.. [#houston_2017] Houston, C. M., Diamanti, E., Diamantaki, M., Kutsarova, E., Cook, A., Sultan, F.,
   & Brickley, S. G. (2017). Exploring the significance of morphological diversity for cerebellar
   granule cell excitability. Scientific Reports, 7(1), 1-16. https://doi.org/10.1038/srep46147
.. [#kanichay_2008] Kanichay, R. T., & Silver, R. A. (2008). Synaptic and cellular properties of the
   feedforward inhibitory circuit within the input layer of the cerebellar cortex. Journal of
   Neuroscience, 28(36), 8955-8967. https://doi.org/10.1523/JNEUROSCI.5469-07.2008
.. [#locatelli_2013] Locatelli, F., Bottà, L., Prestori, F., Masetto, S., & D’Angelo, E. (2013). Late-onset bursts
   evoked by mossy fibre bundle stimulation in unipolar brush cells: Evidence for the involvement of H- and
   TRP-currents. The Journal of Physiology, 591(4), 899–918. https://doi.org/10.1113/jphysiol.2012.242180
.. [#subramaniyam_2014] Subramaniyam, S., Solinas, S., Perin, P., Locatelli, F., Masetto, S., & D’Angelo, E. (2014).
   Computational modeling predicts the ionic mechanism of late-onset responses in unipolar brush cells.
   Frontiers in Cellular Neuroscience, 8. https://doi.org/10.3389/fncel.2014.00237
.. [#geminiani_2018] Geminiani, A., Casellato, C., Locatelli, F., Prestori, F., Pedrocchi, A., & D'Angelo, E. (2018).
   Complex dynamics in simplified neuronal models: reproducing Golgi cell electroresponsiveness. Frontiers in
   neuroinformatics, 12, 88. https://doi.org/10.3389/fninf.2018.00088
