For developers
---------------
The ``cerebellar-models`` package requires Python 3.10+.

Developers best use pip’s *editable* install. This creates a live link between the installed package
and the local git repository:

.. code-block:: bash

    git clone git@github.com:dbbs-lab/cerebellar-models
    cd cerebellar-models
    pip install -e .[dev]
    pre-commit install

Documentation
~~~~~~~~~~~~~

You can produce the documentation for this code using the following commands, assuming you are in
the ``cerebellar-models`` folder:

.. code-block:: bash

   pip install -e .[docs]
   cd docs
   sphinx-build -n -b html . _build/html

Testing
~~~~~~~
You can also run the unit tests and functional tests using the following commands, assuming you are
in the ``cerebellar-models`` folder:

.. code-block:: bash

   pip install -e .[tests]
   python -m unittest discover -v -s ./tests

Configuration file layout
~~~~~~~~~~~~~~~~~~~~~~~~~

The ``configurations/`` directory holds all YAML fragments that the CLI assembles into a final BSB
configuration. A normal user never touches these files — the ``cerebellar-models configure``
command does the merging — but developers need to know where each kind of parameter lives when
checking or changing values.

The tree is organised in three nested levels: **species → state → simulator**.

.. code-block:: text

   configurations/
   └── <species>/                              # e.g. mouse/
       ├── <species>_cerebellar_cortex.yaml    # base circuit (network geometry, placements,
       │                                       # connectivity rules, no simulator parameters)
       ├── cell_types/
       │   ├── dcn.yaml                        # optional: adds DCN cell types
       │   ├── io.yaml                         # optional: adds IO cell types
       │   └── ubc.yaml                        # optional: adds UBC cell types
       └── <state>/                            # e.g. awake/ or in-vitro/
           └── <simulator>/                    # e.g. nest/
               ├── <simulation_scenario>.yaml  # simulation settings and recording devices
               └── cell_models/
                   └── <model_name>.yaml       # neuron model constants + synaptic weights

Layer 1 — base circuit
''''''''''''''''''''''

``<species>_cerebellar_cortex.yaml`` defines the minimal cerebellar cortex: network size, layer
geometry, cell-type densities, placement strategies, and connectivity rules. It contains no
simulator-specific keys. All other fragments are merged on top of it by the CLI.

**Where to edit:** network geometry, cell density / radius, placement or connectivity strategy
names (i.e. anything that applies regardless of simulation state or tool).

Layer 2 — optional cell types
''''''''''''''''''''''''''''''

Each YAML in ``cell_types/`` patches the base circuit with an additional population (DCN, IO, UBC).
The CLI deep-merges the user-selected patches into the base config before proceeding.

**Where to edit:** the definition of an optional population — its density, its connectivity to the
rest of the circuit, and any geometry adjustments (the ``network.z`` extension needed to fit the
new layer).

Layer 3 — state folder
'''''''''''''''''''''''

A *state* (e.g. ``awake``, ``in-vitro``) groups parameter regimes that differ between experimental
conditions. Adding a new sibling folder here is enough to make it selectable in the CLI — no code
change is required.

Inside each state folder, a subfolder per simulator (currently only ``nest/``) contains:

**Simulation scenario files** (``basal_awake.yaml``, ``mf_stimulus.yaml``, …)
  Each file corresponds to one simulation the user can select. It specifies:

  - ``simulations.<name>``: duration, resolution, seed, and the NEST ``devices`` (stimulators and
    recorders) for that scenario.
  - ``packages`` and ``components``: BSB packages and build scripts required at runtime.

  The ``awake`` scenario files use BSB's ``$import`` directive to pull their content from the
  matching ``in-vitro`` files, then override only the parameters that differ between states.
  This avoids duplicating large blocks of YAML. When editing ``awake`` parameters, first check
  whether the key you need is defined locally or inherited from ``in-vitro``.

  **Where to edit:** simulation duration/seed, recording-device targetting, stimulation rates,
  or the list of packages required for a scenario.

**Cell model files** (``cell_models/<model_name>.yaml``)
  Each file corresponds to one neuron model that the user can select during the interactive
  survey. It is structured in three sections:

  - ``cell_models``: per-cell-type NEST ``constants`` (membrane, threshold, reset voltages,
    synaptic time constants, …) and the ``model`` key that names the NEST model to instantiate.
  - ``<synapse_type>_connection_models``: synaptic weight and delay parameters for every
    projection, grouped by synapse model (e.g. ``static_synapse_connection_models``,
    ``tsodyks_synapse_connection_models``). The CLI prompts the user to choose one group; the
    chosen group is renamed to ``connection_models`` in the output and the others are discarded.

  Currently two neuron models are available for NEST:

  - ``eglif_cond_alpha_multisyn`` — conductance-based E-GLIF with alpha-shaped synaptic currents
    (legacy model, available in both ``awake`` and ``in-vitro``).
  - ``eglif_multirec`` — E-GLIF with multi-receptor support (available in ``in-vitro`` only;
    the ``awake`` variant is not yet provided).

  **Where to edit:** a cell type's intrinsic parameters (``cell_models.<cell_type>.constants``),
  synaptic weights or delays (``<synapse_type>_connection_models.<projection_name>``).

Adding a new neuron or synapse model
'''''''''''''''''''''''''''''''''''''

To expose a new NEST model to the CLI, create a YAML file under
``configurations/<species>/<state>/nest/cell_models/``. The file must follow the structure
described above: a ``cell_models`` section with per-cell-type constants and ``model``, plus one or
more ``<synapse_type>_connection_models`` sections. The filename (without ``.yaml``) becomes the
choice label in the survey — no code change is needed.

To add a new synapse model variant within an existing cell-model file, add a new
``<synapse_type>_connection_models`` top-level key alongside the existing ones. The CLI will list
it as an additional synapse-model choice for that neuron model.
