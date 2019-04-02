#!/bin/bash
set -e

# Set up environment for "testing" user here.
export PATH=$HOME/scripts:$PATH
source /usr/local/gromacs/bin/GMXRC
source $HOME/testing/bin/activate

exec "$@"