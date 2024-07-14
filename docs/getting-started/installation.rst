Installation
------------
We advise you to implement the ``cerebellum`` package and the related simulators within a Conda
or Python virtual environment.

Create a virtual environment in the current folder and activate it:

.. code-block:: bash

   python3 -m venv cereb_env
   # activate the environment
   # this command can be added to your .bashrc, adapting the path
   source cereb_venv/bin/activate
   # if you need to leave the environment
   deactivate

The ``cerebellum`` package requires Python 3.9+.

.. code-block:: bash

    git clone git@github.com:dbbs-lab/cerebellum
    cd cerebellum
    pip install -e .

If you wish to contribute to the cerebellum repository, please also install
the :doc:`developers' packages <for-developers>`.

NEST simulator
~~~~~~~~~~~~~~

Download and install the NEST simulator in the current folder within your virtual environment:

.. code-block:: bash

   pip install cython cmake
   git clone https://github.com/nest/nest-simulator
   cd nest-simulator
   mkdir build/ && cd build
   # The minimal cmake instruction is `cmake ..`, but if you wish to run simulation with MPI
   cmake .. -Dwith-mpi=ON
   # To speed up the process, you can set the number of process to compile NEST, here 4
   make install -j4



What next:
~~~~~~~~~~

If you want more information on the Biological context, head over to this
:doc:`section <biological-context>`.

Otherwise, you can continue to the :doc:`Contents section <content>`.
