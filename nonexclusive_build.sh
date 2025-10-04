#!/bin/bash

export PATH=$HOME/gcc-13.2/bin:$PATH
export LD_LIBRARY_PATH=$HOME/gcc-13.2/lib64:$LD_LIBRARY_PATH

# if passed use the first argument as number of threads, otherwise use 6
if [[ $# -eq 2  ]]; then
    BUILD_DIR=$1
    NUM_THREADS=$2
else
    BUILD_DIR=build
    NUM_THREADS=6
fi

if [ -d ${BUILD_DIR} ]; then
    echo "Build directory ${BUILD_DIR} already exists: moving to ${BUILD_DIR}_last"      
    if [ -d ${BUILD_DIR}_last ]; then
        echo "Removing old ${BUILD_DIR}_last"
        rm -rf ${BUILD_DIR}_last/*
    fi
    mv ${BUILD_DIR} ${BUILD_DIR}_last
fi
mkdir ${BUILD_DIR}

cd ${BUILD_DIR} || { echo "Build directory ${BUILD_DIR} could not be created."; exit 1; }
# nonexclusive cmake .. --preset=dev -DKAHYPAR_CI_BUILD=ON -DKAHYPAR_DOWNLOAD_TBB=On # -DKAHYPAR_DOWNLOAD_GMP=On

nonexclusive cmake .. -DCMAKE_C_COMPILER=${HOME}/gcc-13.2/bin/gcc \
        -DCMAKE_CXX_COMPILER=${HOME}/gcc-13.2/bin/g++ \
        -DCMAKE_CXX_FLAGS="-static-libstdc++ -static-libgcc -L${HOME}/gcc-13.2/lib64 -Wl,-rpath,${HOME}/gcc-13.2/lib64" \
        -DKAHYPAR_DOWNLOAD_TBB=On \
        -DCMAKE_BUILD_TYPE=RELEASE


NUM_THREADS=6
if [[ $# -eq 1  ]]; then NUM_THREADS=$1; fi
 
nonexclusive make -j"$NUM_THREADS" MtKaHyPar

