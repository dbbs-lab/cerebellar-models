Installation
------------
We advise you to implement the ``cerebellar-models`` package and the related simulators within a Conda
or Python virtual environment.

The ``cerebellar-models`` package requires Python 3.10+.

Create a virtual environment in the current folder and activate it:

.. code-block:: bash

   python3 -m venv cereb_env
   # activate the environment
   # this command can be added to your .bashrc, adapting the path
   source cereb_venv/bin/activate
   # if you need to leave the environment
   deactivate

You can install the ``cerebellar-models`` package trough `pip`:

.. code-block:: bash

    pip install cerebellar-models

To include simulation support for NEST and/or NEURON:

.. code-block:: bash

    pip install cerebellar-models[nest]  # For NEST
    pip install cerebellar-models[neuron]  # For NEURON

Alternatively, you can install it from sources:

.. code-block:: bash

    git clone git@github.com:dbbs-lab/cerebellar-models
    cd cerebellar-models
    pip install -e .
    # do not forget to add the optional simulator supports required.
    pip install -e .[nest]  # for NEST

If you wish to contribute to the cerebellar-models repository, please also install
the :doc:`developers' packages <for-developers>`.

NEST simulator
~~~~~~~~~~~~~~

``cerebellar-models`` compiles its own NESTML-based cell models (``cerebellar-models
build-nestml``) against your NEST installation. This requires a NEST install with a real,
existing install prefix and its C++ headers available, which not every installation method
below provides.

From source
^^^^^^^^^^^

Download and install the NEST simulator in the current folder within your virtual environment:

.. code-block:: bash

   pip install cython
   git clone https://github.com/nest/nest-simulator
   cd nest-simulator
   mkdir build/ && cd build
   # The minimal cmake instruction is `cmake ..`, but if you wish to run simulation with MPI
   cmake .. -Dwith-mpi=ON
   # To speed up the process, you can set the number of process to compile NEST, here 4
   make install -j4

NEST 3.9 and later additionally require the Boost headers to build against
(``libboost-dev`` on Debian/Ubuntu, ``boost`` via Homebrew, ...).

With conda or mamba
^^^^^^^^^^^^^^^^^^^^

NEST is also available on conda-forge, with a working install prefix that
``cerebellar-models build-nestml`` can compile against:

.. code-block:: bash

   mamba install -c conda-forge nest-simulator cmake compilers libboost-headers

``libboost-headers`` is required alongside ``nest-simulator``: NEST's conda-forge package only
depends on Boost at *build* time, so it is not pulled in automatically even though NEST's
installed headers need it at compile time.

.. warning::
   Do not use ``pip install nest-simulator``. The PyPI wheel bundles a self-contained NEST
   binary without a usable install prefix, and its ``nest-config`` script reports hardcoded
   paths from the wheel's own CI build that do not exist on your machine. As a result,
   ``cerebellar-models build-nestml`` (and NESTML in general) cannot build custom extension
   modules against it — see `nest/nest-simulator#3737
   <https://github.com/nest/nest-simulator/issues/3737>`_. Use the source or conda/mamba routes
   above instead.

Compiling the NESTML models
~~~~~~~~~~~~~~~~~~~~~~~~~~~

``cerebellar-models`` defines its cell and synapse models in NESTML file format
and compiles them into a NEST extension module (``cerebmodule``) on top of the
NEST installation set up above.

This happens automatically the first time it is needed: ``bsb compile``/``bsb simulate`` compile
and deploy ``cerebmodule`` for you when your configuration lists the
``cerebellar_models.nest_models.build_models`` component, e.g.:

.. code-block:: yaml

    components:
      - cerebellar_models.nest_models.build_models
    simulations:
        simulation_name:
            simulator: nest
            modules:
                - cerebmodule

Alternatively, you can trigger (or force) the compilation manually with the
:ref:`cerebellar-models build-nestml <cli-nestml>` CLI command.

What next:
~~~~~~~~~~

If you want more information on the Biological context, head over to this
:doc:`section <biological-context>`.

Otherwise, you can continue to the :doc:`Contents section <content>`.
