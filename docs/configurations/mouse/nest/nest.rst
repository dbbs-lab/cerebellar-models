NEST simulations
----------------

NEST Installation
^^^^^^^^^^^^^^^^^
To reproduce the experiments presented below, you should install the NEST simulator (see
instructions :doc:`here </getting-started/installation>`).

Nest modules are automatically compiled with BSB as ``components`` and deployed as ``cerebmodule``
through the provided configurations:

.. code-block:: yaml

    components:
      - cerebellar_models.nest_models.build_models
    simulations:
        simulation_name:
            simulator: nest
            modules:
                - cerebmodule

Alternatively, you can manually compile them, running the
`build_models.py <https://github.com/dbbs-lab/cerebellar-models/blob/master/cerebellar_models/nest_models/build_models.py>`_
python script.

.. toctree::
    :maxdepth: 2
    :caption: NEST simulations

    parameters
    eglif_cond_alpha_multisyn
    simulation-results
