Customising configurations with BSB imports
-------------------------------------------

The ``cerebellar-models`` configurations are designed to be composable. Rather than copying a full
configuration file and editing it, BSB's ``$import`` directive lets you pull specific sections from
an existing file and override only what you need. This keeps your custom configurations short and
ensures that upstream fixes propagate automatically.

This page explains how ``$import`` works and walks through the patterns used in this repository.

How ``$import`` works
~~~~~~~~~~~~~~~~~~~~~

``$import`` is a special YAML/JSON key recognised by BSB at load time. It merges selected sections
from another configuration file into the document that contains it. The basic shape is:

.. code-block:: yaml

    $import:
      ref: <path-to-file>#<json-pointer>
      values:
        - <key1>
        - <key2>

``ref``
  Path to the source file, optionally followed by a ``#`` and a
  `JSON Pointer <https://datatracker.ietf.org/doc/html/rfc6901>`_ that navigates to the node
  within that file to import from. Use ``#/`` to import from the root of the file, or
  ``#/simulations/basal_activity`` to start from a nested node.

``values``
  List of keys to copy from the referenced node into the current document location. Only the listed
  keys are imported; everything else in the source file is ignored.

**Merge behaviour**: ``$import`` is applied as a recursive merge. If a key already exists in the
current document at the same location, the local definition takes precedence — the imported value
does not overwrite it. This makes it safe to define overrides *before* the ``$import`` block.

``$import`` can appear at the top level of a file or nested inside any mapping, including inside a
simulation definition.

Pattern 1 — importing top-level sections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simplest use case is pulling a set of top-level sections (``packages``, ``components``,
``simulations``, …) wholesale from another file. The ``awake`` simulation files do this to inherit
from the ``in-vitro`` equivalents:

.. code-block:: yaml

    # configurations/mouse/awake/nest/basal_awake.yaml
    $import:
      ref: ../../in-vitro/nest/basal_vitro.yaml#/
      values:
        - packages
        - simulations
        - components

The ``#/`` pointer targets the root of ``basal_vitro.yaml``, so ``packages``, ``simulations`` and
``components`` are imported as-is. The awake file itself adds nothing extra — it is just a thin
alias that can be extended later.

Pattern 2 — importing a subsection into a nested context
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``$import`` can be placed inside a mapping to pull content from a *nested* location of the source
file into the *current* mapping. For example, ``mf_stimulus.yaml`` defines a new simulation that
reuses the cell models and recording devices of the basal activity:

.. code-block:: yaml

    # configurations/mouse/in-vitro/nest/mf_stimulus.yaml
    simulations:
      mf_stimulus:
        simulator: nest
        resolution: 0.1
        duration: 2000
        seed: 1234
        devices:
          mf_stimulus:           # only the new stimulus device is defined here
            device: poisson_generator
            rate: 150
            start: 1200
            stop: 1260
            targetting:
              strategy: sphere
              radius: 90
              origin: [150, 65, 100]
              cell_models:
                - mossy_fibers
            weight: 1
            delay: 0.1
        $import:
          ref: ./basal_vitro.yaml#/simulations/basal_activity
          values:
            - cell_models   # import all cell model definitions from basal_activity
            - devices       # merge the basal recording devices into the new simulation
    $import:
      ref: ./basal_vitro.yaml#/
      values:
        - packages
        - components

Here two ``$import`` blocks work together:

- The nested one (inside ``mf_stimulus``) targets
  ``basal_vitro.yaml#/simulations/basal_activity`` and merges ``cell_models`` and ``devices``
  directly into the ``mf_stimulus`` simulation mapping.
- The top-level one imports ``packages`` and ``components`` from the same file.

The local ``mf_stimulus`` device key defined before the nested import is not overwritten by the
imported ``devices`` block, because local values take precedence.

Pattern 3 — overriding imported values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To change specific values from an imported file, define them locally *before* the ``$import``
block. The local definition wins over the imported one. The awake version of the eyeblink
conditioning stimulus demonstrates this: it inherits the full in-vitro stimulus configuration but
raises the climbing-fiber weight:

.. code-block:: yaml

    # configurations/mouse/awake/nest/mf_cf_stimulus.yaml
    simulations:
      mf_cf_stimulus:
        devices:
          cf_stimulus:
            weight: 100.    # override: higher weight for awake state

    $import:
      ref: ../../in-vitro/nest/mf_cf_stimulus.yaml#/
      values:
        - packages
        - simulations
        - components

The in-vitro file sets ``cf_stimulus.weight`` to ``55``; the local definition of ``100.`` is
already present when the import is applied, so it is preserved.

.. important::
    Note that list in the configuration, will be fully overwritten by the imported values,
    unlike dictionaries which are merged together. If you wish to overwrite a list in the
    configuration to add a element, then you will need to define it entirely in the final
    configuration.

Creating a custom simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As a worked example, suppose you want a mossy fiber stimulus that lasts ``3000 ms`` and targets a
different sphere. Create a new file next to the provided ones:

.. code-block:: yaml

    # my_custom_stimulus.yaml
    simulations:
      my_stimulus:
        simulator: nest
        resolution: 0.1
        duration: 3000          # custom duration
        seed: 42
        devices:
          my_mf_stimulus:
            device: poisson_generator
            rate: 200           # custom rate
            start: 500
            stop: 600
            targetting:
              strategy: sphere
              radius: 50
              origin: [100, 50, 80]
              cell_models:
                - mossy_fibers
            weight: 1
            delay: 0.1
        $import:
          ref: ./basal_vitro.yaml#/simulations/basal_activity
          values:
            - cell_models   # reuse neuron model definitions unchanged
            - devices       # add the basal recording devices

    $import:
      ref: ./basal_vitro.yaml#/
      values:
        - packages
        - components

Further reading
~~~~~~~~~~~~~~~

BSB's configuration system is described in detail in the
:ref:`BSB configuration documentation <bsb:cfg_file_ref>`.
