Configuring reconstruction and simulation randomness
------------------------------------------------------

Placement, connectivity, post-processing hooks, and
simulations all draw random numbers. BSB configures all of that randomness through a single ``rng``
block in the configuration, rather than through separate, strategy-specific ``seed`` attributes.
This page covers how that block is used across ``cerebellar-models``' configurations; for the
full explanation of the mechanism itself see `BSB's randomness guide
<https://bsb.readthedocs.io/en/latest/core/randomness.html>`_.

Why this matters
~~~~~~~~~~~~~~~~

Two things are needed from a stochastic reconstruction, and they pull in opposite directions:

- Running the same configuration repeatedly should give **technical replicates** — networks that
  differ only in their randomness, so that averaging over several runs is meaningful.
- Any one of those runs should be **reproducible** afterwards, so a specific result can be
  checked or re-analysed later.

The ``rng`` block gives both by making an unset seed mean *draw one and write it down*: leave it
unset and every run is an independent replicate; the seed that was actually used is recorded back
into the configuration stored with the run's output, so feeding that configuration back
reproduces the run exactly.

.. warning::
    Never treat several runs that share a fixed seed as a sample. Repeating one seed repeats the
    exact same draws, so it is one run counted several times, not several independent
    observations — see `BSB's randomness guide
    <https://bsb.readthedocs.io/en/latest/core/randomness.html>`_ for a worked example. Leave
    ``seed`` unset whenever you need a set of runs to average over.

Where is it set in the configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``rng`` sits at the root of a configuration, alongside ``network``, ``storage`` and ``packages`` —
not nested inside ``simulations`` or any single strategy. Only
``configurations/mouse/mouse_cerebellar_cortex.yaml`` sets it:

.. code-block:: yaml

    rng:
      seed: 1234

None of the NEST scenario files (``basal_vitro.yaml``, ``mf_stimulus.yaml``, …) declare their own
``rng``, and they do not need to. The ``cerebellar-models configure`` CLI always merges a scenario
file's ``simulations`` (and any other key it defines) *onto* the base circuit configuration, and
``rng`` is an ordinary config node — a merge leaves a key it does not itself define untouched — so
every generated circuit inherits ``mouse_cerebellar_cortex.yaml``'s ``rng`` unless a scenario file
overrides it on purpose.

This is unlike ``packages``: that one *has* to be redeclared in full in every scenario file that
needs an extra package such as ``bsb-nest``, because a list cannot be merged item-by-item — only
replaced wholesale. ``rng``, like ``components``, carries no such restriction.

A scenario file can still pin its own, independent seed by declaring ``rng`` locally — see
:ref:`nest-kernel-seeding` for the kernel case, or "Holding one part of the model fixed" above for
the general mechanism.

Getting a fresh replicate every run
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To let a configuration vary from run to run — for example while collecting a sample of networks
to average structural or spiking measures over — remove the ``seed`` line (or set it to
``null``):

.. code-block:: yaml

    rng:
      seed: null

A fresh seed is then drawn when the network is compiled, and it is written back into the
configuration stored alongside that run's output, so any single run from the sample can still be
reproduced afterwards by copying that recorded seed back into ``rng.seed``.

Holding one part of the model fixed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``rng`` can also declare named ``generators``. An entry with its own ``seed`` stays fixed across
runs; a component (a placement or connectivity strategy) names the entry it draws from through
its own ``rng`` attribute:

.. code-block:: yaml

    rng:
      seed: null
      generators:
        structure:
          seed: 42
    placement:
      granule_layer_placement:
        strategy: bsb.placement.RandomPlacement
        rng: structure
        # ...

Here every run reconstructs the exact same placement of granule cells, while connectivity, any
other placement strategy, and the simulation still vary between runs. This is the difference
between a technical replicate of the whole network and one that holds a chosen part of it fixed.

.. _nest-kernel-seeding:

NEST kernel seeding
~~~~~~~~~~~~~~~~~~~~

NEST simulations used to take their own ``seed`` attribute directly on ``simulations.<name>``.
That attribute no longer exists — see :doc:`the NEST parameters page
</configurations/mouse/nest/parameters>`. By default, the kernel's master seed is derived from the
network's ``rng.seed``, keyed on the simulation's name, so:

- pinning ``rng.seed`` on the base circuit config also pins every NEST simulation built from it,
  and
- two different simulations of the same network still get independent NEST streams from one
  another.

To pin the kernel to a specific seed independently of the rest of ``rng``, declare a
``rng.settings`` entry of strategy ``nest`` and name it on the simulation:

.. code-block:: yaml

    rng:
      seed: 1234
      settings:
        kernel:
          strategy: nest
          seed: 999
    simulations:
      basal_activity:
        rng: kernel
        # ...

Further reading
~~~~~~~~~~~~~~~~

`BSB's randomness guide <https://bsb.readthedocs.io/en/latest/core/randomness.html>`_ covers the
full mechanism: how streams are derived from what is being drawn for rather than from the MPI rank
(so a run reproduces regardless of how many ranks it runs on), and how to write a generator of
your own.
