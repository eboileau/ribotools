.. _user_guide:

User guide
==========

The **Ribotools** box
---------------------

**Ribotools** can be used in combination with `Rp-Bp <http://rp-bp.readthedocs.io/en/latest/>`_, or directly from raw reads, to generate count tables for Ribo-seq, RNA-seq, or both. If you have run **Rp-Bp** on a set of Ribo-seq samples, you can use the output to generate Ribo-seq count tables. In fact, it can also be used with existing alignments (Ribo- and/or RNA-seq), provided the files follow the naming convention as described in the **Rp-Bp** package, see *e.g.* `How to use existing alignment files <https://rp-bp.readthedocs.io/en/latest/existing-alignments.html>`_.

.. important::

    If using existing output or alignment files, do not use the ``--overwrite`` option!

Once Ribo- and RNA-seq count tables are available, **Ribotools** can be used to estimate translation efficiency. This is described in `How to estimate TE <estimate-te.html>`_. **Ribotools** can also be used to estimate differential expression based on Ribo-seq periodic fragment lengths, *i.e.* determine which features, genes or Ribo-seq ORFs, are differentially regulated. This is described in `How to estimate DE <estimate-de.html>`_.

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


For TE or DE analyses, the following keys are required:

* ``contrasts`` *(required, input)* A dictionary *key: value*, where *key* is a name for the contrast to be tested, and *value* contains 2 items, the first item is the condition to be tested against the second (reference).

* ``tea_data`` *(required, output)* The base output location for all created files, wheter you are interested in TE or in DE.

* ``sample_table`` *(optional, input)* The path to a sample table *e.g.* if only running the analysis from existing data.
* ``count_table`` *(optional, input)* The path to a count table *e.g.* if only running the analysis from existing data.


.. A *template* configuration file is available to download with the Tutorials.

.. tip::

    Use keywords such as *Ribo* and *RNA* in ``riboseq_sample_name_map`` and ``rnaseq_sample_name_map`` values, respectively, to facilitate data integration for the analysis.


.. _prepare_genome:


How to prepare genome indices and annotations
---------------------------------------------

To get started, you need reference annotations (GTF), the ribosomal `Bowtie 2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ index, and the `STAR <https://github.com/alexdobin/STAR>`_ index. Protocol-specific or general adapter sequences to be removed are also required for `Flexbar <https://github.com/seqan/flexbar/wiki/Manual>`_. We assume that annotations and indices are already available.

.. tip::

    **Rp-Bp** is installed as a dependency, and can be used to prepare the annotations, see `How to prepare genome indices and annotations <https://rp-bp.readthedocs.io/en/latest/user-guide.html#how-to-prepare-genome-indices-and-annotations>`_. Use the same configuration file as described above, with the necessary keys to run ``prepare-rpbp-genome``.


.. _alignment_workflow:

How to estimate abundance
-------------------------

For default Flexbar and STAR parameters, consult the **Rp-Bp** documentation, in particular `Default parameters and options <https://rp-bp.readthedocs.io/en/latest/user-guide.html#default-parameters-and-options>`_. Note that providing ``--post-trim-length`` as a ``--flexbar-options`` will overwrite ``--trim-rna-to-max-fragment-size``. Default STAR options are for Ribo-seq and trimmed RNA-seq reads. Unless you use ``--trim-rna-to-max-fragment-size``, or to change the mapping parameters, use ``--star-options`` to override defaults, and run Ribo-seq and RNA-seq separately.


.. important::

    All Ribo-seq samples (including biological replicates) in the configuration file must be from the same organism and use the same ``genome_base_path``, ``star_index``, ``ribosomal_index``, *etc.* Samples from different organisms or using different annotations must be "split" into different configuration files, and run separately.


.. note::

    Ribo-seq abundance is estimated from ribosome-protected RNA fragments (RPFs) using periodic reads only, unless the ``--skip-periodicity-estimation`` option is used. **Ribotools** does not currently handle paired-end RNA-seq. In the config file, specify only one pair (typically the first) for ``rnaseq_samples``. But you should specify whether the data is from a strand-specific assay for counting reads, see below for further details. Use the option ``--trim-rna-to-max-fragment-size`` to trim RNA-seq reads to the maximum periodic fragment length of a matched Ribo-seq sample, to minimize mapping bias. Matching samples are specified in the config file via the ``matching_samples`` key.


.. tip::

    If some Ribo-seq samples fail to pass the periodicity estimation quality control, you can still recover matching RNA-seq samples, provided there is enough replicates to run the TE analysis. Simply replace the corresponding Ribo-seq sample name for the matching RNA-seq sample (*value* from ``matching_samples``) with an integer *e.g.* ``!!int 31``. The value will be used for ``--trim-rna-to-max-fragment-size``, and can correspond to the mean or median periodic fragment size across all available libraries. Alternatively, you can use `fixed lengths and offsets <https://rp-bp.readthedocs.io/en/latest/user-guide.html#fixed-lengths-and-offsets>`_, although this is not recommended.



.. _ribotools_usage:

General usage
^^^^^^^^^^^^^

**Ribotools** can be called for Ribo-seq only with the ``ribo`` option, RNA-seq only with the ``rna`` option, or for both Ribo- and RNA-seq at once, using the ``ribo --run-all`` option. **Ribotools** can also be run with the `SLURM <http://slurm.schedmd.com>`_ scheduler. For all options, use ``run-htseq-workflow -h``.

.. code-block:: bash

    # Ribo-seq only (estimate periodicity and generate count tables)
    # e.g. unpaired protocol (default)
    run-htseq-workflow ribo <config> [options]

    # RNA-seq only (generate count tables)
    # e.g. pair 1 (reverse protocol) for RNA
    # with matching Ribo-seq data use [--trim-rna-to-max-fragment-size] [--ribo-config RIBO_CONFIG]
    run-htseq-workflow rna <config> --htseq-options "--stranded reverse" [options]

    # Ribo- and RNA-seq at once, i.e. one after the other
    # e.g. stranded protocol for Ribo passed via [--htseq-options] (default is no),
    # and pair 1 (reverse protocol) for RNA (passed via [--stranded], since we use [--run-all])
    run-htseq-workflow ribo <config> --run-all [--trim-rna-to-max-fragment-size] [--rna-config RNA_CONFIG] --htseq-options "--stranded yes" --stranded reverse [options]


If Ribo-seq ORFs are available from **Rp-Bp**, TE or DE can be estimated for Ribo-seq ORFs, instead of genes (CDS by default, or exon). In this case, you need to prepare a GTF file with Ribo-seq ORFs before, and ``run-htseq-workflow`` with additional options, see `How to estimate TE using Ribo-seq ORFs <ribo-seq-orfs.html>`_ for details.


.. important::

    If using **Ribotools** with a *de novo* assembly generated with **Rp-Bp**, specifying ``--htseq-options --type=exon``, or type other than ``CDS`` can have unexpected results! This is because the GTF file created under ``genome_base_path`` is a concatenation of ``gtf`` and ``de_novo_gtf``, and possibly contains repeated features (see `How to prepare genome indices and annotations <https://rp-bp.readthedocs.io/en/latest/user-guide.html#how-to-prepare-genome-indices-and-annotations>`_). For mapping this is not a problem. For abundance estimation, however, this can be problematic. Unless this GTF file is manually curated, only CDS features should be used (default).


.. tip::

    Use ``--star-options "--quantMode GeneCounts"`` to get count tables. You can check counts for unstranded data (column 2), counts for the 1st read strand (column 3, htseq-count -s yes), and counts for the 2nd read strand (column 4, htseq-count -s reverse). The stranded column (3 or 4) with the lowest *N_noFeature* count should correspond to the correct strand option.


.. tip::

    You can use **Rp-Bp** to perform read filtering quality control, use the ``-k/--keep-intermediate-files`` option. Intermediate files *e.g.* Flexbar, or Bowtie2 output can be deleted afterwards. See `Visualization and quality control <https://rp-bp.readthedocs.io/en/latest/apps.html>`_.


Output files
^^^^^^^^^^^^

Except for *orf_profiles*, all output files follow the conventions described in **Rp-Bp** `output files <https://rp-bp.readthedocs.io/en/latest/user-guide.html#id10>`_. For RNA-seq, we follow the same conventions and nomenclature. Count tables are written to *<riboseq_data>/count-tables* and *<rnaseq_data>/count-tables*.


Default parameters and options
------------------------------

The parameters and options decribed below are all optional. All parameters and options have default values that do not normally need to be modified.


.. note::

    **Rp-Bp** parameters can be changed via the configuration file, and options for external programs (Flexbar, STAR) are handled via command line arguments.
    You do not need to include **Rp-Bp** parameters in the configuration file, unless you wish to change their values.


Check `Default parameters and options <https://rp-bp.readthedocs.io/en/latest/user-guide.html#default-parameters-and-options>`_.


HTSeq
^^^^^

Default options are overridden via command line using ``--htseq-options``.


* ``format`` Default: bam.
* ``stranded`` Default: no.
* ``type`` Default: CDS.
* ``idattr`` Default: gene_id.
* ``additional-attr`` Default: gene_name.
* ``mode`` Default: intersection-nonempty.
* ``secondary-alignments`` Default: ignore.
* ``supplementary-alignments`` Default: ignore.
