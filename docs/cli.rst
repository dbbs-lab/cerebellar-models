######################
Command-Line Interface
######################

Open up your favorite terminal and enter the ``cerebellum --help`` command
to verify you correctly installed the software.

The list of cerebellum commands are listed below.

.. note::
  Parameters included between angle brackets are example values, parameters between square
  brackets are optional, leave off the brackets in the actual command.


Retrieve a canonical circuit configuration
==========================================

.. warning::
    For the time being, this command has to be run inside the cerebellum folder.

.. code-block:: bash

  cerebellum configure [--output_folder <./path>] [--species <mouse>] [--extension <yaml>]

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
* ``--path``: File Extension for the configuration. Can be either ``json`` or ``yaml``.
