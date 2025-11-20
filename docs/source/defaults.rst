.. _defaults:

Default parameters
==================

The parameters and options decribed below are all optional. All parameters and options have default values that do not normally need to be modified.

.. important::

    **Ribotools** parameters are changed via the configuration file, but options for programs such as Flexbar and STAR are handled via command line arguments. You do not need to include **Ribotools** parameters in the configuration file, unless you wish to change their values.


Flexbar and STAR options
------------------------

For default parameters, consult the **Rp-Bp** documentation, in particular `Flexbar and STAR options <https://rp-bp.readthedocs.io/en/latest/defaults.html#flexbar-and-star-options>`_.

.. caution::

   Providing ``--post-trim-length`` as a ``--flexbar-options`` will overwrite ``--trim-rna-to-max-fragment-size``.

.. note::

   Default STAR options are used for both Ribo-seq and trimmed RNA-seq reads. Unless you use ``--trim-rna-to-max-fragment-size``, or to
   change the mapping parameters, use ``--star-options`` to override defaults, and run Ribo-seq and RNA-seq separately.

Ribotools parameters
--------------------

For default parameters, consult the **Rp-Bp** documentation, in particular `Rp-Bp parameters <https://rp-bp.readthedocs.io/en/latest/defaults.html#rp-bp-parameters>`_ (general, shared MCMC, and metagene and periodicity estimation parameters).

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
