This is a fork of the [main Gromacs project](http://www.gromacs.org/) in which interface, API, and extensibility issues are being investigated.
The forked project lives on Bitbucket at http://bitbucket.org/kassonlab/gromacs

In addition to a regular GROMACS installation, this fork provides `libgmxapi` for
high-level C++ access to GROMACS MD simulation.
It exists primarily to support the `gmxpy` companion project that provides a [Python module](http://bitbucket.org/kassonlab/gmxpy/)

This README.md file supplants the main README file to avoid merge conflicts while providing convenient documentation to the BitBucket repository browser.

# Installation

Install as you would a regular copy of GROMACS.

    $ git clone git clone https://bitbucket.org/kassonlab/gromacs.git
    $ mkdir build
    $ cd build
    $ cmake ../gromacs -DCMAKE_INSTALL_PREFIX=/path/to/where/i/want/gromacs
    $ make install

You may then either source the gmxrc file as usual or export the environment variable
`gmxapi_DIR=/path/to/where/i/want/gromacs` to help `gmxpy` or your own CMake project to find
what it needs to build against the gmxapi library.

# Documentation

To build additional documentation, use the additional build targets `gmxapi_cppdocs` and `gmxapi_cppdocs_dev`.
You will need to have `doxygen` installed.
If you would like to write code that uses `libgmxapi`, `make gmxapi_cppdocs`.
For more detail, or if you would like to extend or contribute to the API, `make gmxapi_cppdocs_dev`.

Then refer either to `docs/html/doxygen/api-user/index.html` or
`docs/html/doxygen/api-dev/index.html` in the build directory.

Also, please use the issue tracking system on Bitbucket or feel free to suggest other modes of communication.
