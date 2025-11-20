Ribotools
=========

Introduction
------------

**Ribotools** is a toolbox for the analysis of matched ribosome profiling (Ribo-seq) and RNA sequencing (RNA-seq) data. **Ribotools** can be used in combination with `Rp-Bp <http://rp-bp.readthedocs.io/en/latest/>`_, or as a standalone software to generate count tables from raw Ribo-seq reads, raw RNA-seq reads, or both, and to calculate translation efficiency (TE) and differential expression (DE).

**Ribotools** uses `Rp-Bp <http://rp-bp.readthedocs.io/en/latest/>`_ for periodicity estimation, and follows the same alignment workflow, directory structure, and naming convention. It uses `HTSeq <https://htseq.readthedocs.io/en/master/>`_ for abundance estimation.

**Ribotools** can find TE and/or DE regulated Ribo-seq ORFs across conditions, using Ribo-seq periodic fragment lengths only. Both approaches use `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_ to accommodate for sample-to-sample variance and/or confounding factors.


.. toctree::
   :titlesonly:

   getting-started
   installation
   user-guide
   tutorial
   faq
   api
