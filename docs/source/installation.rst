.. _installation_full:

Installation
============

Containers
----------

To use a container (Docker or Singularity) with **Ribotools** pre-installed, simply pull, and you're done!

.. code-block:: bash

   # docker or...
   docker pull quay.io/biocontainers/ribotools:<tag>
   # ...singularity
   singularity pull ribotools.sif docker://quay.io/biocontainers/ribotools:<tag>

There is no *latest* tag, you need to specify the version tag. See `ribotools/tags <https://quay.io/repository/biocontainers/ribotools?tab=tags>`_ for valid values for <tag>.

..
  Check the `Tutorials <tutorial.html>`_ on how to use the containers.


.. _conda_install:

Conda installation
------------------

If required, set up the conda channels as described `here <https://bioconda.github.io/#usage>`_, and install with

.. code-block:: bash

    # preferably install in some conda environment...
    conda install ribotools

or create an environment, called *ribotools*, containing the **Ribotools** package

.. code-block:: bash

    conda create -n ribotools ribotools

.. tip::

    `Mamba <https://mamba.readthedocs.io/en/latest/user_guide/mamba.html#mamba>`_ can be used as a drop-in replacement, you can swap almost all commands between conda and mamba.

.. _pypi_install:

Contributing to **Ribotools**
-----------------------------

To install the local VCS project in development mode

.. code-block:: bash

    # create a conda environment...
    mamba env create -n ribotools
    # ... activate...
    mamba activate ribotools
    # ... and install dependencies...
    mamba install --only-deps ribotools
    # ... clone the git repository and install the package
    git clone https://github.com/eboileau/ribotools.git && cd ribotools
    pip install --no-deps --editable . 2>&1 | tee install.log

Alternatively, clone the git repository and install from the `yaml spec file <https://mamba.readthedocs.io/en/latest/user_guide/mamba.html#conda-yaml-spec-files>`_

.. code-block:: bash

   mamba create -f environment.yml
   # ... activate environment...
   mamba activate ribotools
   pip install --no-deps --editable . 2>&1 | tee install.log

Finally, install test dependencies.

PyPI installation
^^^^^^^^^^^^^^^^^

To install the package from `PyPI <https://pypi.org/project/ribotools>`_

.. code-block:: bash

    # create a virtual environment...
    python3 -m venv ribotools
    # ... activate...
    source ribotools/bin/activate
    # ... and install the package
    pip install ribotools

**Required dependencies:** `Flexbar <https://github.com/seqan/flexbar>`_, `Bowtie 2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_, `STAR <https://github.com/alexdobin/STAR>`_, `Samtools <http://www.htslib.org>`_, and `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_. You also need a working `R <https://www.r-project.org/>`_ installation with additional packages.

.. warning::

    A PyPI installation only installs the python package. You need to install required dependencies separately. Executables or binaries must be in your ``$PATH``.

.. _uninstall:

Uninstallation
--------------

Remove the conda environment

.. code-block:: bash

    mamba env remove --name ribotools

or remove the package installed in another environment

.. code-block:: bash

    # remove the ribotools package from myenv environment...
    mamba remove -n myenv ribotools

To remove **Ribotools** if installed with pip

.. code-block:: bash

    pip uninstall ribotools

If the package is installed in a dedicated python virtual environment, remove this environment.
