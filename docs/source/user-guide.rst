.. _user_guide:

User guide
==========

You want to generate count tables from matched Ribo- and RNA-seq data? First, you need to prepare genome indices and annotations for your organism. This has to be done once for any given reference genome and annotation. Consult :ref:`howto_annotation`.

You can estimate translation efficiency, you can find regulated Ribo-seq ORFs across conditions, without matched RNA-seq data, and much more. Consult :ref:`running_htseq_workflow`, :ref:`running_te` or :ref:`running_de`.

.. hint::

   You can use the output of **Rp-Bp** as input to **Ribotools**. In fact, **Ribotools** can be used with existing alignments (Ribo- and/or
   RNA-seq), provided the files follow the naming convention as described in the **Rp-Bp** package, see
   `How to use existing alignment files <https://rp-bp.readthedocs.io/en/latest/existing-alignments.html>`_. If using existing output or
   alignment files, do not use the ``--overwrite`` option!

.. toctree::
   :maxdepth: 1

   howto-config
   howto-annotation
   howto-run
   howto-te
   howto-de
   howto-riboseq-orfs
   defaults
