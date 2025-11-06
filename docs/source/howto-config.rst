.. _howto_config:

How to add a config file
========================

Create separate YAML configuration files, or a single YAML configuration file, for Ribo- and RNA-seq. For Ribo-seq, consult `How to add a config file <https://rp-bp.readthedocs.io/en/latest/howto-config.html>`_.

.. note::

   Required keys such as ``fasta`` and ``ribosomal_fasta`` are required only if you use **Rp-Bp** to create annotations, see :ref:`howto_annotation`.


For RNA-seq, additional keys are required:

* ``rnaseq_samples`` *(required, input)* A dictionary *key: value*, where *key* is used to construct filenames, and *value* is the full path to the FASTQ.gz file for a given sample. The *key* should not contain spaces or special characters.
* ``matching_samples`` *(optional, input)* A dictionary *key: value*, where *key* is the same as ``rnaseq_samples`` *key*, and *value* is the same as ``riboseq_samples`` *key* (matched samples). Required when using ``--trim-rna-to-max-fragment-size``.
* ``rna_adapter_file`` *(optional, input)* Path to adapter sequences to be removed (FASTA).
* ``rnaseq_data`` *(required, output)* The base output location for all created files.
* ``rnaseq_sample_name_map`` *(optional, output)* A dictionary *key: value*, where *key* is the same as ``rnaseq_samples`` *key*, and *value* is a fancy name for *key* to use in downstream analyses.

For TE or DE analyses, the following keys are required:

* ``contrasts`` *(required, input)* A dictionary *key: value*, where *key* is a name for the contrast to be tested, and *value* contains 2 items, the first item is the condition to be tested against the second (reference).
* ``tea_data`` *(required, output)* The base output location for all created files (TE).
* ``dea_data`` *(required, output)* The base output location for all created files (DE).
* ``sample_table`` *(optional, input)* The path to a sample table *e.g.* if running the analysis from existing data.
* ``count_table`` *(optional, input)* The path to a count table *e.g.* if running the analysis from existing data.

.. attention::

   Use ``sample_table`` and/or ``count_table`` only if these files are already available, *e.g.* to estimate TE or DE using
   existing sample and/or count tables. To automatically create these tables from the output of the pipeline, do not use these keys.

To download an example configuration file, check the test dataset included with the :ref:`all_tutorials`, or consult the examples :ref:`prep_tables_te` or :ref:`prep_tables_de`. To change the default parameters, see :ref:`defaults`.
