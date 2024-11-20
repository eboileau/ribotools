Ribotools
=========

Introduction
------------

**Ribotools** is a toolbox for the analysis of matched ribosome profiling (Ribo-seq) and RNA sequencing (RNA-seq) data. **Ribotools** can be used in combination with `Rp-Bp <http://rp-bp.readthedocs.io/en/latest/>`_, or as a standalone software to generate count tables for raw Ribo-seq reads, raw RNA-seq reads, or both, and to calculate translation efficiency (TE) and differential expression (DE).

**Ribotools** uses `Rp-Bp <http://rp-bp.readthedocs.io/en/latest/>`_ for periodicity estimation, and follows the same alignment workflow, directory structure, and naming convention. It uses `HTSeq <https://htseq.readthedocs.io/en/master/>`_ for abundance estimation.

Translation efficiency can be defined as the number of ribosomes per feature, typically a gene, normalized to abundance. A feature can be regulated transcriptionally and/or translationally, resulting in several different regulatory profiles. Differential expression is based on Ribo-seq periodic fragment lengths and can be used to find regulated Ribo-seq ORFs across conditions. Both approaches use `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_ to accommodate for sample-to-sample variance, complex experimental designs, and/or confounding factors.


.. toctree::
   :titlesonly:

   getting-started
   installation
   user-guide
   estimate-te
   estimate-de
   ribo-seq-orfs
