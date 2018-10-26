#!/bin/bash
set -ev

SOURCEDIR=$HOME/mpich-3.2.1
BUILDDIR=/tmp/mpich-build
INSTALLDIR=$MPICH_DIR

rm -rf $SOURCEDIR
rm -rf $BUILDDIR
rm -rf $INSTALLDIR

pushd $HOME
    wget http://www.mpich.org/static/downloads/3.2.1/mpich-3.2.1.tar.gz
    tar zxf mpich-3.2.1.tar.gz
    mkdir $BUILDDIR
    pushd $BUILDDIR
        mkdir $INSTALLDIR
        $SOURCEDIR/configure \
            --prefix=$INSTALLDIR \
            --enable-shared \
            --disable-fortran \
            CC=`which gcc-${GCC}` \
            CXX=`which g++-${GCC}`
        make -j2 install
    popd
popd
