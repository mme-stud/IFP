#! /bin/bash

export PATH=$HOME/gcc-13.2/bin:$PATH
export LD_LIBRARY_PATH=$HOME/gcc-13.2/lib64:$LD_LIBRARY_PATH


if [ "$#" -ne 7 ]; then
    echo "Illegal number of parameters. Usage: <build_dir> <preset> <version> <data_root_dir> <mode> <scenario> <graphs>"
    echo "Example: run_AON_HMLL_experiments.sh build AON_HMLL_clique survey_benchmark DCHSBM scenA1 "1,2,3""
    exit 1
fi

JULIA=~/HyperModularity.jl/julia-1.11.4/bin/julia
BUILD_DIR=$1
echo "Using build directory: $BUILD_DIR"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/mt-kahypar/IFP/${BUILD_DIR}/lib

preset=$2
version=$3
data_root_dir=$4
mode=$5
scenario=$6
graphs=$7


if [ -d ${data_root_dir}/${mode}/${scenario}/${version} ]; then
    echo "removing directory ${data_root_dir}/${mode}/${scenario}/${version}"
    rm -rf ${data_root_dir}/${mode}/${scenario}/${version}
fi

mkdir ${data_root_dir}/${mode}/${scenario}/${version}

rm -rf ${data_root_dir}/${mode}/${scenario}/${version}.results

${JULIA} AON_HMLL_exp.jl -v "$version" -d "$data_root_dir/" -m "$mode/" -s "$scenario/" -i "$graphs"

echo "Analyzing AON_HMLL results for version=$version, data_root_dir=$data_root_dir, mode=$mode, scenario=$scenario, graphs=$graphs"

python3 analyze_results_via_metrics.py --version $version --data_root_dir $data_root_dir --mode $mode --scenario $scenario --graphs $graphs
./conductance.out -p $preset -v $version -m $mode -s $scenario -instances $graphs
