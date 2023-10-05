Getting started
===============

What is **Ribotools**?
------------------

**Ribotools** is a toolbox for the analysis of matched ribosome profiling (Ribo-seq) and RNA sequencing (RNA-seq) data. It can be used to generate count tables for Ribo-seq, RNA-seq, or both, and perform translation efficiency (TE) analysis.

Ribo-seq abundance is estimated from ribosome-protected RNA fragments (ribosome footprints) using periodic reads only. Optionally, RNA-seq reads can be trimmed to the maximum periodic fragment length of a matched Ribo-seq sample, to minimize mapping bias.

To get started, you need

* Ribo-seq data (FASTQ) - or output from `Rp-Bp <http://rp-bp.readthedocs.io/en/latest/>`_
* matched RNA-seq data (FASTQ) - if performing TE analysis

* annotation for your organism (GTF)
* indices for `Bowtie 2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_  and `STAR <https://github.com/alexdobin/STAR>`_
* protocol-specific or general adapter sequences to be removed (FASTA)

.. _getting_started:

Installation
------------

Install with

.. code-block:: bash

    # currently only local install
    # see https://github.com/mamba-org/mamba/issues/633
    mamba create -n ribotools
    mamba activate ribotools
    git clone https://github.com/eboileau/ribotools.git && cd ribotools
    mamba env update -n ribotools --file environment.yml
    pip --verbose install . 2>&1 | tee install.log


.. For detailed installation instructions, refer to `Installation <installation.html>`_.


**Ribotools** quickstart
------------------------

**Ribotools** alignment workflow wraps calls to `Flexbar <https://github.com/seqan/flexbar/wiki/Manual>`_, `Bowtie 2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_, and `STAR <https://github.com/alexdobin/STAR>`_. Indices must be available.

To generate count tables, simply call

.. code-block:: bash

    run-htseq-workflow <{rna,ribo}> <config> [options]

To perform TE analysis, call

.. code-block:: bash

    run-tea -config CONFIG


For more information and guidelines on how to prepare the configuration file, refer to the `User guide <user-guide.html>`_.


How to report issues
--------------------

Bugs and issues should be reported in the `bug tracker <https://github.com/eboileau/ribotools/issues>`_. Follow the instructions and guidelines given in the template.


How to contribute
-----------------

Contributions are welcome! New code should follow `Black <https://black.readthedocs.io/en/stable/>`_ and `flake8 <https://flake8.pycqa.org/en/latest/>`_. Install development dependencies inside a virtual environment. A typical development workflow would include *(i)* forking the repository, *(ii)* creating a new branch for your PR, *(iii)* adding features or bug fixes, *(iv)* making sure all tests are passing, *(v)* building the documentation if necessary, and *(vi)* opening a PR back to the main repository. If you're fixing a bug, add a test. Run it first to confirm it fails, then fix the bug, and run it again to confirm it's fixed. If adding a new feature, add a test, or first open an issue to discuss the idea.

Running the tests
^^^^^^^^^^^^^^^^^

Dependencies can be installed with ``pip install -e .[tests]``.

Building the docs
^^^^^^^^^^^^^^^^^

Dependencies for building the documentation can be installed with ``pip install -e .[docs]``.

Semantic versioning
^^^^^^^^^^^^^^^^^^^

We try to follow `semantic versioning <https://semver.org/>`_.


License
-------

The MIT License (MIT). Copyright (c) 2019 Etienne Boileau.
