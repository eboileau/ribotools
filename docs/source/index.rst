Ribotools
=========

Introduction
------------

**Ribotools** is a toolbox for the analysis of matched ribosome profiling (Ribo-seq) and RNA sequencing (RNA-seq) data. **Ribotools** can be used in standalone mode or in combination with `Rp-Bp <http://rp-bp.readthedocs.io/en/latest/>`_ for Ribo-seq quality control, open reading frames (ORFs) discovery, translation efficiency (TE), and differential expression (DE) analyses.

**Ribotools** uses `Rp-Bp <http://rp-bp.readthedocs.io/en/latest/>`_ for periodicity estimation and follows the same alignment workflow, directory structure, and naming convention. It uses `HTSeq <https://htseq.readthedocs.io/en/master/>`_ for abundance estimation and `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_ for statistical testing.


.. toctree::
   :titlesonly:

   getting-started
   installation
   user-guide
   tutorial
   faq
   api
