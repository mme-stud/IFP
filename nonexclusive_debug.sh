#!/bin/bash

export PATH=$HOME/gcc-13.2/bin:$PATH
export LD_LIBRARY_PATH=$HOME/gcc-13.2/lib64:$LD_LIBRARY_PATH

# enter the name of build directory. It should begin with build_debug_

read -p "Enter the name of build directory (default: build_debug): build_debug" build_debug
BUILD_DIR=${build_debug:-$BUILD_DIR}
BUILD_DIR=build_debug${BUILD_DIR}

rm -rf ${BUILD_DIR}_last/*
mv ${BUILD_DIR} ${BUILD_DIR}_last

mkdir ${BUILD_DIR}
cd ${BUILD_DIR}
nonexclusive cmake .. --preset=dev -DKAHYPAR_CI_BUILD=ON -DKAHYPAR_DOWNLOAD_TBB=On -DCMAKE_BUILD_TYPE=Debug # -DKAHYPAR_DOWNLOAD_GMP=On

NUM_THREADS=6
if [[ $# -eq 1  ]]; then NUM_THREADS=$1; fi
 
nonexclusive make -j"$NUM_THREADS" MtKaHyPar

cd ../_debug

#nonexclusive ./run_prompt.sh
#nonexclusive ./summary_results.sh
