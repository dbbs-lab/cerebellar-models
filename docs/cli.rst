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
