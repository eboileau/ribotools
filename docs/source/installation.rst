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
    mamba create -n ribotools
    # ...activate it...
    mamba activate ribotools
    # ... and install dependencies
    git https://github.com/eboileau/ribotools.git && cd ribotools
    mamba env update -n ribotools --file environment.yml
    # install
    pip --verbose install --editable .[tests,docs] 2>&1 | tee install.log


PyPI installation
^^^^^^^^^^^^^^^^^

We do not recommend to install **Ribotools** directly from `PyPI <https://pypi.org/project/ribotools>`_.
However, if you already have the required dependencies installed on your system, to install

.. code-block:: bash

    # create a virtual environment...
    python3 -m venv ribotools_pypi
    # ... activate it ...
    source ribotools_pypi/bin/activate
    # ... and install Ribotools (ribotools_pypi is now activated)
    pip install ribotools


**Required dependencies:** `Flexbar <https://github.com/seqan/flexbar>`_, `Bowtie 2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_, `STAR <https://github.com/alexdobin/STAR>`_, `Samtools <http://www.htslib.org>`_.

.. warning::

    Conda installation or containers include all dependencies. With a PyPI installation, you need to install required dependencies. Executables or binaries must be in your ``$PATH``.


.. _uninstall:

Uninstallation
--------------

Remove the conda environment

.. code-block:: bash

    conda env remove --name ribotools

or remove the package installed in another environment

.. code-block:: bash

    # remove the ribotools package from myenv environment...
    conda remove -n myenv ribotools


To remove **Ribotools** if installed with pip

.. code-block:: bash

    pip uninstall ribotools

If the package is installed in a dedicated python virtual environment, this environment can also be removed.
