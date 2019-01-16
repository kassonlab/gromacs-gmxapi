#!/bin/bash
set -ev

export GMX_DOUBLE=OFF

if [ ${CI_MPI} -eq 0 ] ; then
    export GMX_MPI=OFF
    export GMX_THREAD_MPI=ON
else
    export GMX_MPI=ON
    export GMX_THREAD_MPI=OFF
fi

rm -rf build
mkdir build
pushd build
    cmake -DGMX_BUILD_HELP=OFF \
         -DGMX_ENABLE_CCACHE=ON \
         -DCMAKE_CXX_COMPILER=$CXX \
         -DCMAKE_C_COMPILER=$CC \
         -DGMX_DOUBLE=$GMX_DOUBLE \
         -DGMX_MPI=$GMX_MPI \
         -DGMX_THREAD_MPI=$GMX_THREAD_MPI \
         -DGMXAPI=ON \
         -DCMAKE_INSTALL_PREFIX=$HOME/install/gromacs \
         ..
    make -j2 check
    make -j2 install
popd
ccache -s
