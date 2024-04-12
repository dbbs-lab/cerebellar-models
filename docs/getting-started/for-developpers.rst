For developpers
---------------
The ``cerebellum`` package requires Python 3.9+.

Developers best use pip’s *editable* install. This creates a live link between the installed package
and the local git repository:

.. code-block:: bash

    git clone git@github.com:dbbs-lab/cerebellum
    cd cerebellum
    pip install -e .[dev]
    pre-commit install

Documentation
~~~~~~~~~~~~~

You can produce the documentation for this code using the following commands, assuming you are in
the ``cerebellum`` folder:

.. code-block:: bash

   pip install -e .[docs]
   cd docs
   sphinx-build -n -b html . _build/html

Testing
~~~~~~~
You can also run the unit tests and functional tests using the following commands, assuming you are
in the ``cerebellum`` folder:

.. code-block:: bash

   pip install -e .[tests]
   pytest test/
