.. _user_guide:

User guide
==========

The **Ribotools** box
---------------------

**Ribotools** can be used in combination with `Rp-Bp <http://rp-bp.readthedocs.io/en/latest/>`_, or directly from raw reads, to generate count tables for Ribo-seq, RNA-seq, or both. If you have run **Rp-Bp** on a set of Ribo-seq samples, you can use the output to generate Ribo-seq count tables. In fact, it can also be used with existing alignments (Ribo- and/or RNA-seq), provided the files follow the naming convention as described in the **Rp-Bp** package, see *e.g.* `How to use existing alignment files <https://rp-bp.readthedocs.io/en/latest/existing-alignments.html>`_.

.. important::

    If using existing output or alignment files, do not use the ``--overwrite`` option!

Once Ribo- and RNA-seq count tables are available, **Ribotools** can be used to estimate translation efficiency.

.. _top:
.. use with `back to top <#top>`_

How to prepare the configuration file
-------------------------------------

A single YAML configuration file can be used for both Ribo- and RNA-seq. For Ribo-seq, consult the **Rp-Bp** documentation, in particular `How to prepare the configuration file <https://rp-bp.readthedocs.io/en/latest/user-guide.html#how-to-prepare-the-configuration-file>`_. For RNA-seq, additional keys are required:


* ``rnaseq_samples`` *(required, input)* A dictionary *key: value*, where *key* is used to construct filenames, and *value* is the full path to the FASTQ.gz file for a given sample. The *key* should not contain spaces or special characters.
* ``matching_samples`` *(required, input)* A dictionary *key: value*, where *key* is the same as ``rnaseq_samples`` *key*, and *value* is the same as ``riboseq_samples`` *key* (matched samples)

* ``rna_adapter_file`` *(optional, input)* Path to adapter sequences to be removed (FASTA).

* ``rnaseq_data`` *(required, output)* The base output location for all created files.

* ``rnaseq_sample_name_map`` *(optional, output)* A dictionary *key: value*, where *key* is the same as ``rnaseq_samples`` *key*, and *value* is a fancy name for *key* to use in downstream analyses.


.. tip::

    Use keywords such as *Ribo* and *RNA* in ``riboseq_sample_name_map`` and ``rnaseq_sample_name_map`` values, respectively, to facilitate data integration for TE analysis.


.. _alignment_workflow:

How to estimate abundance
-------------------------

To get started, you need reference annotations (GTF), the ribosomal `Bowtie 2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ index, and the `STAR <https://github.com/alexdobin/STAR>`_ index. Protocol-specific or general adapter sequences to be removed are also required for `Flexbar <https://github.com/seqan/flexbar/wiki/Manual>`_.


.. tip::

    **Rp-Bp** is installed as a dependency, and can be used to prepare the annotations, see `How to prepare genome indices and annotations <https://rp-bp.readthedocs.io/en/latest/user-guide.html#how-to-prepare-genome-indices-and-annotations>`_.


For default Flexbar and STAR parameters, consult the **Rp-Bp** documentation, in particular `Default parameters and options <https://rp-bp.readthedocs.io/en/latest/user-guide.html#default-parameters-and-options>`_. These are overridden via command line arguments.


.. important::

    All Ribo-seq samples (including biological replicates) in the configuration file must be from the same organism and use the same ``genome_base_path``, ``star_index``, ``ribosomal_index``, *etc.* Samples from different organisms or using different annotations must be "split" into different configuration files, and run separately.


.. note::

    Ribo-seq abundance is estimated from ribosome-protected RNA fragments (ribosome footprints) using periodic reads only, unless the ``--skip-periodicity-estimation`` option is used. **Ribotools** does not currently handle paired-end RNA-seq. In the config file, specify only one pair (typically the first) for ``rnaseq_samples``. Use the option ``--trim-rna-to-max-fragment-size`` to trim RNA-seq reads to the maximum periodic fragment length of a matched Ribo-seq sample, to minimize mapping bias. Matching samples are specified in the config file via the ``matching_samples`` key.


.. note::

    Default STAR options are for Ribo-seq and trimmed RNA-seq reads. Unless you use ``--trim-rna-to-max-fragment-size``, or to change the default parameters, use ``--star-options`` to override defaults, see `Default parameters and options <https://rp-bp.readthedocs.io/en/latest/user-guide.html#default-parameters-and-options>`_ for more details.


.. _ribotools_usage:

General usage
^^^^^^^^^^^^^

**Ribotools** can be called for Ribo-seq only with the ``ribo`` option, RNA-seq only with the ``rna`` option, or for both Ribo- and RNA-seq at once, using the ``ribo --run-all`` option.

.. code-block:: bash

    # Ribo-seq only (estimate periodicity and generate count tables)
    run-htseq-workflow ribo <config> [options]

    # RNA-seq only (generate count tables)
    # with matching Ribo-seq data use [--trim-rna-to-max-fragment-size] [--ribo-config RIBO_CONFIG]
    run-htseq-workflow rna <config> [options]

    # Ribo- and RNA-seq at once, i.e. one after the other
    run-htseq-workflow ribo <config> --run-all [--trim-rna-to-max-fragment-size] [--ribo-config RIBO_CONFIG] [options]


- examples

HTSeq options:
  --htseq-options [HTSEQ_OPTIONS ...]
                        Optional argument: a space-delimited list of options to pass to htseq-count. Each option must be quoted separately as in "--htseqOption value",
                        using soft quotes, where '--htseqOption' is the long parameter name from htseq-count and 'value' is the value given to this parameter. If
                        specified, htseq-count options will override default settings. (default: None)


  --stranded {yes,reverse,no}
                        Optional argument: library strandedness for RNA-seq. This option is passed to htseq-count and overrides the same option passed via [--htseq-
                        options] and used for Ribo-seq. Unless given, the default value will be used. (default: no)
  --trim-rna-to-max-fragment-size
                        Flag: trim RNA post adapter removal using max fragment size from the matching Ribo-seq sample. Note* At least the "periodic-offsets" file must be
                        available. The config file must also include "matching_samples" and the path to the Ribo-seq config must be given [--ribo-config]) (default:
                        False)
  --ribo-config RIBO_CONFIG
                        Optional argument: the Ribo-seq config file if using [--trim-rna-to-max-fragment-size]. (default: None)
  --rna-config RNA_CONFIG
                        Optional argument: the RNA-seq config file if using [--run-all]. (default: None)




**Ribotools** can be run with the `SLURM <http://slurm.schedmd.com>`_ scheduler. For all options, use ``run-htseq-workflow -h``.


.. tip::

    You can use **Rp-Bp** to perform read filtering quality control, use the ``-k/--keep-intermediate-files`` option. Intermediate files *e.g.* Flexbar, or Bowtie2 output can be deleted afterwards. See `Visualization and quality control <https://rp-bp.readthedocs.io/en/latest/apps.html>`_.



- output files same as Rp-bp + count-tables, as in docs


----


Default parameters and options
------------------------------

The parameters and options decribed below are all optional. All parameters and options have default values that do not normally need to be modified.


.. important::

    **Rp-Bp** parameters can be changed via the configuration file, and options for external programs (Flexbar, STAR) are handled via command line arguments.
    You do not need to include **Rp-Bp** parameters in the configuration file, unless you wish to change their values.


check Rp-Bp docs, add htseq


Flexbar and STAR options
^^^^^^^^^^^^^^^^^^^^^^^^

Default options for external programs (Flexbar, STAR) are overridden via command line using ``--flexbar-options`` or ``--star-options``. Currently, no options can be passed to Bowtie2.

Flexbar
"""""""

* ``max-uncalled`` Default: 1.
* ``pre-trim-left`` Default: 0.
* ``qtrim-format`` Default: sanger.
* ``qtrim`` Default: TAIL.
* ``qtrim-threshold`` Default: 10.
* ``zip-output`` Default: GZ.


check Rp-Bp docs

Rp-Bp parameters
^^^^^^^^^^^^^^^^
