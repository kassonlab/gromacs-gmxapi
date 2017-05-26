For web-based documentation, please visit http://gromacs-api.readthedocs.io/

This is a fork of the main Gromacs project in which interface, API, and extensibility issues are being investigated.
This README.md file supplants the main README file to avoid merge conflicts while providing convenient documentation to the BitBucket repository browser.
The usual Gromacs README is still here without the .md file extension.

The first few commits in this repository after forking from master are just documentation built in Sphinx.
To build and view, do a normal Gromacs build with the ``GMX_API=ON`` cmake advanced option. E.g.

    git clone https://bitbucket.org/kassonlab/gromacs.git
    mkdir build; cd build
    cmake ../gromacs/ -DGMX_API=ON
    make

Then open `src/api/docs/html/index.html`

Please use the Bitbucket issue tracking system and other features to make comments, feature requests,
pull requests, etcetera, or email M. Eric Irrgang.

Please note this repository is in migration and will rapidly evolve over the week or two following the 2017 summer Gromacs workshop at MPI.

updated 26 May, 2017