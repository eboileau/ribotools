Frequently asked questions
==========================

* :ref:`q1`
* :ref:`q2`
* :ref:`q3`
* :ref:`q4`
* :ref:`q5`
* :ref:`q6`

.. _q1:

I don't want/I can't install it on my computer. Can I still use **Ribotools** without installing the package?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Yes, you can use a Docker or a Singularity container. Simply pull, and you're done! See :ref:`installation_full` for instructions.
Example calls are also given in the user guide and tutorials.

.. _q2:

I have alignments, can I use **Ribotools**?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The short answer is yes. The pipeline is designed to handle all steps from raw FASTQ files up to the final count tables, but you can start the pipeline from any step. Check the `Rp-Bp tutorial on how to use existing alignment files <https://rp-bp.readthedocs.io/en/latest/existing-alignments.html>`_.

.. _q3:

I have count tables, can I use **Ribotools** for TE/DE only?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Yes, check :ref:`prep_tables_te_general` (TE) or :ref:`prep_tables_de_general` (DE).

.. _q4:

Can I quality control (QC) my data?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For Ribo-seq, **Ribotools** estimates periodicity by default, so you know which samples and read lengths are usable for downstream analyses.
If you use **Rp-Bp** to create index files and run the pipeline with the ``--create-orf-profiles`` option, you can use the `Rp-Bp profile construction dashboard <https://rp-bp.readthedocs.io/en/latest/howto-qc.html#summarizing-the-profile-construction>`_ to facilitate visualization and QC of your Ribo-seq data. See :ref:`ribotools_qc` for more details. For RNA-seq, there are already plenty of QC tools.

.. _q5:

I called ``run-htseq-workflow`` with the ``--run-all`` flag, why are there no RNA count tables?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If one of your Ribo-seq sample failed to complete, *e.g.* due to low quality (no periodic reads were found), all RNA samples
will fail to run. In such cases, you can try to troubleshoot the problem by looking at the logs. You can decide to drop failed
samples, or use fixed lengths and offsets, and re-run the workflow. The pipeline will skip existing files and continue were it failed,
unless you use the ``--overwrite`` flag. This behavior may be different whether you are running samples sequentially or in parallel, *e.g.*
with the ``--use-slurm`` option.

.. _q6:

For TE/DE estimation, I get errors about undefined columns or NA values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you get errors such as ``undefined columns selected``, ``Error in DESeqDataSet: NA values are not allowed``, or ``Gene IDs (first column) differ between files``, this is generally due to not specifying ``--orfCol``, ``--symbolCol``, and/or the correct ``--delim`` for the count table. Check :ref:`using_orfs_tede` for some hints.
