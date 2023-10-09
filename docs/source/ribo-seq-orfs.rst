How to estimate TE using Ribo-seq ORFs
======================================

By default, TE is estimated based on RPFs and RNA abundance for ``type=CDS``, *i.e.* all other feature types (3rd column in GTF file) are ignored. If **Rp-Bp** has been used for Ribo-seq ORFs discovery, this set of Ribo-seq ORFs can be used to generate a new GTF file, and TE can be estimated for Ribo-seq ORFs, instead of genes (CDS or exon). This also works with a *de novo* assembly, see *e.g.* `More about prepare-rpbp-genome <https://rp-bp.readthedocs.io/en/latest/rpbp-genome.html>`_.

You first need to locate the output of ``summarize-rpbp-predictions``, in particular the BED12+ file *<project_name>[.note][-unique][.filtered].predicted-orfs.bed.gz*, see `Summarizing the Rp-Bp predictions <https://rp-bp.readthedocs.io/en/latest/apps.html#summarizing-the-rp-bp-predictions>`_. This file contains the combined predicted translation events, or Ribo-seq ORFs, from all samples and replicates, with additional columns. To prepare a GTF file with the Ribo-seq ORFs

.. code-block:: bash

    get-gtf-from-predictions bed out [options]

where ``bed`` is the path to *<project_name>[.note][-unique][.filtered].predicted-orfs.bed.gz*, and ``out`` if the path to the GTF file to be created, *e.g.* *<project_name>[.note][-unique][.filtered].predicted-orfs.gtf*.

This file contain all the Ribos-eq ORFs, with CDS and exon features, and additional attributes, such as ORF type.

You are now ready to call

.. code-block:: bash

    run-htseq-workflow ribo <config> --run-all [--trim-rna-to-max-fragment-size] [--rna-config RNA_CONFIG] --htseq-options "--idattr orf_id" "--additional-attr orf_type" "--additional-attr gene_name" --gtf HTSEQ_GTF [options]

The important options are

* ``--htseq-options "--idattr orf_id"``, this will indicate ``htseq-count`` to use the ORF id as feature for counting. You can add additional options, such as ``"--additional-attr orf_type"``, *etc.*
* ``--gtf HTSEQ_GTF``, where ``HTSEQ_GTF`` is the newly created GTF file *<project_name>[.note][-unique][.filtered].predicted-orfs.gtf*, this will indicate the program to suse the Ribo-seq ORFs for abundance estimation, and not the existing GTF file.


Finally, you need to check the ``htseq-count`` output tables, to determine which column contains which attributes. Currently, ``htseq-count`` outputs ``"--additional-attr`` in the order as they are given, but this behaviour is not guaranteed. In particular if using multiple additional attributes, we need to make sure the right column is used, and they can be specified with ``<-symbolCol COLUMN>`` and ``<-orfCol COLUMN_NUMBER>``, where ``symbolCol`` is the ``gene_name`` column, and ``orfCol`` is the ORF type column, if used.


.. code-block:: bash

    run-tea -config [CONFIG] <-method LRT/deltaTE> <-lfcThreshold L2FC> <-alpha ALPHA> <-symbolCol COLUMN> <-orfCol COLUMN_NUMBER> <-delim TAB/CSV> <-batch> <-filter>


.. note::

    By default, **DESeq2** does not specify the field separator for ``read.table``, *i.e.* is uses the default white space. If some GTF attributes are missing, *e.g.* ``gene_name`` with a *de novo* assembly for novel Ribo-seq ORFs, ``htseq-count`` will omit the attribute, and the count file will not be read correctly. In this case, you should specify the field separator using ``-delim TAB``.
