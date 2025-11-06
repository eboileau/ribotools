.. _api_sample:

get-sample-table
================

Create a sample table

.. code-block:: bash

   usage: get-sample-table [--help] [--ribo] [--rna] [--dea] config


Positional Arguments
--------------------
``config``
    Yaml config file (same as used for ``run-htseq-workflow``)

Flags
-----
-h, --help  show this help message and exit
-r, --ribo  Prepare for Ribo-seq samples only
--rna       Prepare for RNA-seq samples only
-d, --dea   Use config key ``dea_data`` (default: ``tea_data``)
