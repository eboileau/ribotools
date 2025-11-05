.. _running_htseq_workflow:

How to estimate abundance
=========================

See :ref:`ribotools_usage` for a short description of required input and expected output. Output from Flexbar, Bowtie2, and STAR are written in FASTQ or `BAM <https://samtools.github.io/hts-specs/>`_ formats.


.. attention::

    All Ribo-seq samples (including biological replicates) in the configuration file must be from the same organism and use the same ``genome_base_path``, ``star_index``, ``ribosomal_index``, *etc.* Samples from different organisms or using different annotations must be "split" into different configuration files, and run separately.


.. note::

    Ribo-seq abundance is estimated from ribosome-protected RNA fragments (RPFs) using periodic reads only, unless the ``--skip-periodicity-estimation`` option is used. **Ribotools** does not currently handle paired-end RNA-seq. In the config file, specify only one pair (typically the first) for ``rnaseq_samples``. But you should specify whether the data is from a strand-specific assay for counting reads, see below for further details. Use the option ``--trim-rna-to-max-fragment-size`` to trim RNA-seq reads to the maximum periodic fragment length of a matched Ribo-seq sample, to minimize mapping bias. Matching samples are specified in the config file via the ``matching_samples`` key.


.. tip::

    If some Ribo-seq samples fail to pass the periodicity estimation quality control, you can still recover matching RNA-seq samples, provided there is enough replicates to run the TE analysis. Simply replace the corresponding Ribo-seq sample name for the matching RNA-seq sample (*value* from ``matching_samples``) with an integer *e.g.* ``!!int 31``. The value will be used for ``--trim-rna-to-max-fragment-size``, and can correspond to the mean or median periodic fragment size across all available libraries. Alternatively, you can use `fixed lengths and offsets <https://rp-bp.readthedocs.io/en/latest/defaults.html#fixed-lengths-and-offsets>`_, although this is not recommended.

.. _ribotools_usage:

General usage
-------------

**Ribotools** can be called for Ribo-seq only with the ``ribo`` option, RNA-seq only with the ``rna`` option, or for both Ribo- and RNA-seq at once, using the ``ribo --run-all`` option.

.. code-block:: bash

    # Ribo-seq only (estimate periodicity and generate count tables)
    # e.g. stranded protocol (default is "no", i.e. not from a strand-specific protocol)
    run-htseq-workflow ribo --htseq-options "--stranded yes" [options] config

    # RNA-seq only (generate count tables)
    # e.g. pair 1 (reverse protocol) for RNA
    # with matching Ribo-seq data use [--trim-rna-to-max-fragment-size] [--ribo-config RIBO_CONFIG]
    run-htseq-workflow rna --htseq-options "--stranded reverse" [options] config

    # Ribo- and RNA-seq at once, i.e. one after the other
    # e.g. trim RNA post adapter removal using max fragment size from matching Ribo-seq samples
    run-htseq-workflow ribo --run-all --trim-rna-to-max-fragment-size --ribo-config CONFIG --rna-config CONFIG --htseq-options "--stranded yes" --rna-stranded reverse [options] config

For all options, consult the API for :ref:`api_workflow`. Even if you use the same configuration file for Ribo- and RNA-seq, you may have to pass it multiple times depending on the selected options, *e.g.* via ``--ribo-config`` and/or ``--rna-config``, in addition to the required positional argument ``config``. See also :ref:`howto_config`.

.. caution::

    If using **Ribotools** with a *de novo* assembly generated with **Rp-Bp**, specifying ``--htseq-options --type=exon``, or type other than ``CDS`` can have unexpected results! This is because the GTF file created under ``genome_base_path`` is a concatenation of ``gtf`` and ``de_novo_gtf``, and possibly contains repeated features (see `More about de novo ORF discovery <https://rp-bp.readthedocs.io/en/latest/rpbp-genome.html#more-about-de-novo-orf-discovery>`_). For mapping this is not a problem. For abundance estimation, however, this can be problematic. Unless this GTF file is manually curated, only CDS features should be used (default).


.. hint::

   To estimate abundance for Ribo-seq ORFs instead of genes (CDS by default, or exon), you need to prepare a GTF file with Ribo-seq ORFs before, and then ``run-htseq-workflow`` with additional options, see :ref:`using_riboseq_orfs`.

.. tip::

    To perform Ribo-seq read filtering quality control (QC), use the ``-k/--keep-intermediate-files`` option, and the
    **Rp-Bp** profile construction dashboard, see `Visualization and QC <https://rp-bp.readthedocs.io/en/latest/howto-qc.html>`_.

.. tip::

    Use ``--star-options "--quantMode GeneCounts"`` to get count tables. You can check counts for unstranded data (column 2), counts for the 1st read strand (column 3, htseq-count -s yes), and counts for the 2nd read strand (column 4, htseq-count -s reverse). The stranded column (3 or 4) with the lowest *N_noFeature* count should correspond to the correct strand option.

Required input
^^^^^^^^^^^^^^

All the input files are given in the configuration file.

Expected output
^^^^^^^^^^^^^^^

Except for *orf_profiles*, all output files follow the conventions described in **Rp-Bp** `ORF profile construction <https://rp-bp.readthedocs.io/en/latest/howto-run.html#orf-profile-construction>`_. For RNA-seq, the same conventions and nomenclature are used. Count tables are written to *<riboseq_data>/count-tables* and *<rnaseq_data>/count-tables*.
