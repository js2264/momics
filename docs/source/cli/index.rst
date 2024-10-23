.. _cli:

CLI Reference
=============

.. toctree::
  :maxdepth: 1

Contents
--------

.. program:: momics
.. code-block:: shell

    momics [OPTIONS] COMMAND [ARGS]...

.. list-table::
    :widths: 25 100
    :align: left
    :header-rows: 1

    * - Data I/O
      -
    * - `momics ingest <#momics-ingest>`_
      - Ingest datasets
    * - `momics cp <#momics-cp>`_
      - Copy datasets
    * - `momics remove <#momics-remove>`_
      - Remove datasets


.. list-table::
    :widths: 25 100
    :align: left
    :header-rows: 1

    * - Data query engine
      -
    * - `momics query <#momics-query>`_
      - Query a repository


.. list-table::
    :widths: 25 100
    :align: left
    :header-rows: 1

    * - Cloud configuration
      -
    * - `momics config <#momics-config>`_
      - Create a repository


.. list-table::
    :widths: 25 100
    :align: left
    :header-rows: 1

    * - Repository management
      -
    * - `momics create <#momics-create>`_
      - Create a repository
    * - `momics delete <#momics-delete>`_
      - Delete a repository
    * - `momics ls <#momics-ls>`_
      - List tables
    * - `momics tree <#momics-tree>`_
      - List repository content
    * - `momics consolidate <#momics-consolidate>`_
      - Consolidate a repository to optimize storage and query performance
    * - `momics manifest <#momics-manifest>`_
      - Generate a manifest of the repository configuration and timestamps


.. list-table::
    :widths: 25 100
    :align: left
    :header-rows: 1

    * - Momics utils
      -
    * - `momics binnify <#momics-binnify>`_
      - Binnify a genome
    * - `momics version <#momics-version>`_
      - Get the version of the momics package

.. rubric:: Options

.. option:: -v, --version

    Show the version and exit.

.. option:: -h, --help

    Show help menu and exit.


----

.. click:: momics.cli:cli
  :prog: momics
  :nested: full
