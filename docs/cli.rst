######################
Command-Line Interface
######################

Open up your favorite terminal and enter the ``cerebellar-models --help`` command
to verify you correctly installed the software.

The list of ``cerebellar-models`` commands are listed below.

.. note::
  Parameters included between angle brackets are example values, parameters between square
  brackets are optional, leave off the brackets in the actual command.


Retrieve a canonical circuit configuration
==========================================

.. warning::
    For the time being, this command has to be run inside the `cerebellar-models` folder.

.. code-block:: bash

  cerebellar-models configure [--output_folder <./path>] [--species <mouse>] [--extension <yaml>] [--microzones]

This command will construct a BSB configuration file based on the canonical cerebellar circuit
developed by the DBBS.

You will be presented with a list of forms to fill directly in the terminal.
Use the up and down arrows of your keyboard to select an option that you want to modify, then
press tab to select it.

Options can either be a list or an entry to manually fill (for text or number). For lists,
you would have to select one or many items. An item can be chosen with the up and down arrow,
selected by pressing the right arrow and validated with the enter key.


* ``--species``: Species to build the configuration from.
* ``--output_folder``: Path to the output folder where the configuration will be stored.
* ``--extension``: File extension for the configuration. Can be either ``json`` or ``yaml``.
* ``--microzones``: Flag to split your circuit into 2 separated microzones (not set by default).

.. note::
    When selecting a ``NEST`` simulation via the CLI, the output configuration will also
    contain an additional simulation based on the associated stimulation protocol
    (see :ref:`NEST paradigms <nest-paradigms>`).


Re-compile NESTML neuron models
================================

.. code-block:: bash

  cerebellar-models build-nestml [--model_dir <path>] [--build_dir <path>] [--module_name <name>]

This command forces a full re-compilation of the NESTML neuron models used by the NEST
simulator. It is equivalent to calling ``_build_nest_models(redo=True)`` programmatically:
the existing build cache is cleared and PyNESTML re-generates the C++ source code and
reinstalls the NEST module ``cerebmodule``.

Running this command is useful when:

* you have updated or added a ``.nestml`` model file and want to force the rebuild without
  clearing the cache manually;
* the cached module has become stale after a NEST or PyNESTML upgrade.

.. note::
    This command requires a working NEST and PyNESTML installation
    (install with ``pip install cerebellar-models[nest]``).

* ``--model_dir``: Directory that contains the ``.nestml`` source files.
  Defaults to the ``nest_models`` folder inside the installed package.
* ``--build_dir``: Directory where the compiled NEST module will be written.
  Defaults to the user-level cache directory managed by the package
  (``appdirs.user_cache_dir("cerebellar_models")``).
