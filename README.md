# Description

This repository contains various scripts, wrappers or snippets of Python3 code that are useful
for pre-/post-processing data associated with the `Rp-Bp` and/or `B-tea` pipelines, but that have not
yet been integrated into these projects. This could be due to the fact that they are associated 
with recent (under development) data processing and/or visualisation or have a limited 
usage/application, or else that there are similar functionalities available, and they have not yet 
been superseded.

## Getting Started

There is currently no documentation available.

### Prerequisites

Pinned version of selected dependencies are listed in the `requirements.txt` file for reproducible installation.
In particular, this requires a running installation of [pbio](https://github.com/dieterich-lab/pybio-utils) and [rpbp](https://github.com/dieterich-lab/rp-bp/tree/1.1.12).

Note: There may be some issues with dependencies if this package is not installed concurrently with `Rp-Bp`!

### Installation

To install the local VCS project in development mode, use the `--editable` or `-e` option, otherwise
this flag can be ignored. 

To install `pproc` and dependencies, first create a virtual environment:
 
```
python3 -m venv /path/to/virtual/environment
```

For information about Python virtual environments, see the [venv](https://docs.python.org/3/library/venv.html) documentation.
To activate the new virtual environment and install `pproc`:

```
# Activate the new virtual environment.
source /path/to/virtual/environment/bin/activate

# If necessary, upgrade pip and wheel or additional packages (such as setuptools if installing in editable mode).
pip install --upgrade pip setuptools wheel

# Clone the git repository
git clone https://github.com/eboileau/pyproc-utils
cd pyproc-utils

# The period is required, it is the local project path (pybio-utils)
pip --verbose install -r requirements.txt [-e] . 2>&1 | tee install.log

```

#### Anaconda installation

The package can also be installed within an [anaconda](https://www.continuum.io/) environment. 

```
# Create the anaconda environment.
conda create -n my_new_environment python=3.6 anaconda

# Activate the new environment.
source activate my_new_environment

# Clone the git repository
git clone https://github.com/eboileau/pyproc-utils
cd pyproc-utils

pip --verbose install -r requirements.txt [-e] . 2>&1 | tee install.log
```

## Uninstallation

To remove the `pproc` package:

```
pip uninstall pproc
```

If the package is installed in a dedicated virtual environment, this environment can also be cleared or removed.

## Running the tests

## Contributing

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

Some scripts are adapted from the `Rp-Bp` and/or `B-tea` pipeline, authored by Brandon Malone.

