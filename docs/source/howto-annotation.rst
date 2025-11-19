.. _howto_annotation:

How to prepare annotations
==========================

For abundance estimation only. Run this once for any given reference genome, as long as the package version and its dependencies remain the same.

.. _genome_usage:

General usage
-------------

If you are not interested in periodicity estimation, you only need reference annotations (GTF), the ribosomal `Bowtie 2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ index, and the `STAR <https://github.com/alexdobin/STAR>`_ index. Protocol-specific or general adapter sequences to be removed are also required for `Flexbar <https://github.com/seqan/flexbar/wiki/Manual>`_.

To prepare the Bowtie2 index

.. code-block:: bash

   bowtie2-build-s [options] ribosomal_fasta ribosomal_index

where ``ribosomal_fasta`` and ``ribosomal_index`` must match the values from the configuration file.

To prepare the STAR index

.. code-block:: bash

   STAR --runMode genomeGenerate [options] --genomeDir star_index --genomeFastaFiles fasta

where ``fasta`` and ``star_index`` must match the values from the configuration file. See :ref:`howto_config`.

For periodicity estimation (default), a BED12+ file with annotated transcripts is also required. You can create this file with

.. code-block:: bash

   gtf-to-bed12 [options] gtf bed

where ``gtf`` is the reference GTF annotation (must match the value from the configuration file), and ``bed`` is the BED12+ file to generate. To use this file, it must be placed under the directory specified by ``genome_base_path`` and its name must match the pattern *<genome_name>.annotated.bed.gz*.

.. tip::

    Use **Rp-Bp** to prepare all annotations and indices at once. Use the same configuration file with the necessary keys to
    run ``prepare-rpbp-genome``. See `How to prepare annotations <https://rp-bp.readthedocs.io/en/latest/howto-annotation.html>`_ for
    more information. This is a must if you want to QC your Ribo-seq data, see :ref:`ribotools_qc`.
