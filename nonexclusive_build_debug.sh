#!/bin/bash

export PATH=$HOME/gcc-13.2/bin:$PATH
export LD_LIBRARY_PATH=$HOME/gcc-13.2/lib64:$LD_LIBRARY_PATH

rm -rf build_debug_last/*
mv build_debug build_debug_last
mkdir build_debug
cd build_debug
nonexclusive cmake .. --preset=dev -DKAHYPAR_CI_BUILD=ON -DKAHYPAR_DOWNLOAD_TBB=On -DCMAKE_BUILD_TYPE=Debug # -DKAHYPAR_DOWNLOAD_GMP=On


NUM_THREADS=6
if [[ $# -eq 1  ]]; then NUM_THREADS=$1; fi
 
nonexclusive make -j"$NUM_THREADS" MtKaHyPar

