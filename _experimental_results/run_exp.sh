#! /bin/bash

export PATH=$HOME/gcc-13.2/bin:$PATH
export LD_LIBRARY_PATH=$HOME/gcc-13.2/lib64:$LD_LIBRARY_PATH


if [ "$#" -ne 7 ]; then
    echo "Illegal number of parameters. Usage: run_exp.sh <build_dir> <preset> <version name> <mode> <scenario> <num_threads> <num_graphs>"
    exit 1
fi
BUILD_DIR=$1
echo "Using build directory: $BUILD_DIR"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HOME}/mt-kahypar/IFP/${BUILD_DIR}/lib

preset=$2
version=$3
mode=$4
scenario=$5
num_threads=$6
instances=$7

if [ -d survey_benchmark/${mode}/${scenario}/${version} ]; then
    echo "removing directory survey_benchmark/${mode}/${scenario}/${version}"
    rm -rf survey_benchmark/${mode}/${scenario}/${version}
fi

mkdir survey_benchmark/${mode}/${scenario}/${version}

./partition.out -p $preset -v $version -m $mode -s $scenario -t $num_threads -instances $instances


# Args:
# $2 - version / config name
# $3 - mode
# $4 - scenario
# $5 - num_threads
# $6 - graphs (no spaces, comma-separated)