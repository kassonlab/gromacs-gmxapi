This is a fork of the [main Gromacs project](http://www.gromacs.org/) in which interface, API, and extensibility issues are being investigated.
The forked project lives on GitHub at [https://github.com/kassonlab/gromacs-gmxapi](https://github.com/kassonlab/gromacs-gmxapi/)

[![Build Status](https://travis-ci.org/kassonlab/gromacs-gmxapi.svg?branch=master)](https://travis-ci.org/kassonlab/gromacs-gmxapi)

In addition to a regular GROMACS installation, this fork provides `libgmxapi` for
high-level C++ access to GROMACS MD simulation.
It exists primarily to support the [`gmxapi`](https://github.com/gmxapi/gmxapi) companion project that provides a Python module and bindings.

This README.md file supplants the main README file to avoid merge conflicts while providing convenient documentation to the repository browser.

# Installation

Install as you would a regular copy of GROMACS. The following example downloads the source into a directory named `gromacs`,
creates a parallel (out-of-source) `build` directory, configures, builds, and installs. Use e.g. `make -j10 install` to build in parallel with 10 processes.

    $ git clone https://github.com/kassonlab/gromacs-gmxapi gromacs
    $ mkdir build
    $ cd build
    $ cmake ../gromacs -DCMAKE_INSTALL_PREFIX=/path/to/where/i/want/gromacs -DGMX_THREAD_MPI=ON -DGMX_GPU=OFF
    $ make install

You may then either source the GMXRC file (as usual for GROMACS use) or export the environment variable
`gmxapi_DIR=/path/to/where/i/want/gromacs` to help `gmxapi` clients such as the Python 
package or your own CMake project to find
what it needs to build against the gmxapi library.

# Documentation

To build additional documentation, use the additional build targets `gmxapi_cppdocs` and `gmxapi_cppdocs_dev`.
You will need to have `doxygen` installed.
If you would like to write code that uses `libgmxapi`, use `make gmxapi_cppdocs`.
For more detail, or if you would like to extend or contribute to the API, `make gmxapi_cppdocs_dev`.

Then refer either to `docs/html/doxygen/api-user/index.html` or
`docs/html/doxygen/api-dev/index.html` in the build directory.

Also, please use the issue tracking system or feel free to suggest other modes of communication.
