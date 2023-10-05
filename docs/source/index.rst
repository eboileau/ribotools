Ribotools
=========

Introduction
------------

**Ribotools** is a toolbox for the analysis of matched ribosome profiling (Ribo-seq) and RNA sequencing (RNA-seq) data. It uses `Rp-Bp <http://rp-bp.readthedocs.io/en/latest/>`_ for periodicity estimation, and follows the same alignment workflow, directory structure, and naming convention.

**Ribotools** can be used in combination with `Rp-Bp <http://rp-bp.readthedocs.io/en/latest/>`_, or directly from raw reads, to generate count tables for Ribo-seq, RNA-seq, or both. It uses `HTSeq <https://htseq.readthedocs.io/en/master/>`_ for abundance estimation.

**Ribotools** include scripts to calculate translation efficiency (TE). Translation efficiency can be defined as the number of ribosomes per gene, normalized to transcript abundance. A gene can be regulated transcriptionally and/or translationally, resulting in several different regulatory profiles. It uses `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_ to accommodate for sample-to-sample variance, complex experimental designs, and/or confounding factors.


.. toctree::
   :titlesonly:

   getting-started
   user-guide
   estimate-te
..    installation
