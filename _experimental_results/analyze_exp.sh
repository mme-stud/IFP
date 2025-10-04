#! /bin/bash

export PATH=$HOME/gcc-13.2/bin:$PATH
export LD_LIBRARY_PATH=$HOME/gcc-13.2/lib64:$LD_LIBRARY_PATH


if [ "$#" -ne 8 ]; then
    echo "Illegal number of parameters. Usage: analyze_exp.sh <build_dir> <preset> <version> <data_root_dir> <mode> <scenario> <num_threads> <graphs>"
    echo "Example: analyze_exp.sh build test cluster survey_benchmark DCHSBM scenA1 3 "1,2,3""
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
num_threads=$7
graphs=$8


rm  ${data_root_dir}/${mode}/${scenario}/${version}.results

echo "Analyzing results for version=$version, data_root_dir=$data_root_dir, mode=$mode, scenario=$scenario, num_threads=$num_threads, graphs=$graphs"

python3 analyze_results_via_metrics.py --version $version --data_root_dir $data_root_dir --mode $mode --scenario $scenario --graphs $graphs
${JULIA} analyze_results_via_modularity.jl -v "$version" -d "$data_root_dir/" -m "$mode/" -s "$scenario/" -t $num_threads -i "$graphs"
./conductance.out -p $preset -v $version -m $mode -s $scenario -instances $graphs
