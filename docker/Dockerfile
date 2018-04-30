# Build GROMACS with modifications for gmxapi
# See https://github.com/kassonlab/gromacs-gmxapi for more information.
# This image serves as a base for integration with the gmxapi Python tools and sample code.

FROM jupyter/scipy-notebook

# Allow the build type to be specified with `docker build --build-arg TYPE=something`
ARG TYPE=Release
# Allow build for an arbitrary branch or tag, but default to the tip of `master`
ARG BRANCH=master

# jupyter/scipy-notebook sets the USER to jovyan. Temporarily switch it back to root to install more stuff.
USER root

RUN \
    apt-get update && \
    apt-get -y upgrade && \
    apt-get install -y software-properties-common && \
    apt-get install -y byobu curl git htop man unzip vim wget

RUN DEBIAN_FRONTEND=noninteractive apt-get install -y libopenblas-dev

RUN \
  DEBIAN_FRONTEND=noninteractive apt-get install -y \
  gdb doxygen mscgen graphviz openmpi-bin libopenmpi-dev libfftw3-dev vim

# Remove cache to save space in the image, but only after all apt-get activity is done.
RUN rm -rf /var/lib/apt/lists/*

USER jovyan
RUN conda install pytest tox cmake

RUN cd /home/jovyan && \
    wget https://github.com/kassonlab/gromacs-gmxapi/archive/$BRANCH.tar.gz && \
    tar xvf $BRANCH.tar.gz && \
    mv gromacs-gmxapi-$BRANCH gromacs-source && \
    rm $BRANCH.tar.gz

# build shared object library with OpenMPI support.
# Note: AVX2 instructions not available in older docker engines.
RUN (mkdir -p /home/jovyan/gromacs-build && \
    cd /home/jovyan/gromacs-build && \
    rm -f CMakeCache.txt && \
    cmake /home/jovyan/gromacs-source \
        -DCMAKE_INSTALL_PREFIX=/home/jovyan/install/gromacs \
        -DGMX_THREAD_MPI=ON \
        -DGMX_BUILD_HELP=OFF \
        -DGMX_SIMD=AVX_256 \
        -DGMX_USE_RDTSCP=OFF \
        -DGMX_HWLOC=OFF \
        -DBUILD_TESTING=ON \
        -DCMAKE_BUILD_TYPE=$TYPE && \
    make -j4)

RUN (cd /home/jovyan/gromacs-build && make -j4 tests)

RUN (cd /home/jovyan/gromacs-build && make -j4 install)

RUN (cd /home/jovyan/gromacs-build && make check)

RUN rm -rf /home/jovyan/gromacs-source && \
    rm -rf /home/jovyan/gromacs-build

# We should also try tests with an MPI-connected set of docker containers

# To be able to step through with gdb, run with something like the following, replacing
# 'imagename' with the name of the docker image built with this recipe.
# docker run --rm -ti --security-opt seccomp=unconfined imagename bash
