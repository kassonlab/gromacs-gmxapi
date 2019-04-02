#!/bin/bash
set -e

# Set up environment for "testing" user here.
cd $HOME
export PATH=$HOME/scripts:$PATH
source $HOME/testing/bin/activate

exec "$@"
