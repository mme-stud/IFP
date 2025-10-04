#! /bin/bash

export PATH=$HOME/gcc-13.2/bin:$PATH
export LD_LIBRARY_PATH=$HOME/gcc-13.2/lib64:$LD_LIBRARY_PATH

function print_usage {
    echo "Usage: $0 [build_directory]"
    echo "If no build_directory is provided, 'build_lib' will be used by default."
}

function compile_library {
    local build_dir="${HOME}/mt-kahypar/IFP/${1}"
    cd $build_dir || { echo "Build directory $build_dir does not exist."; exit 1; }
    make -j8 install-mtkahypar
    cd - || { echo "Failed to return to previous directory."; exit 1; }
}

function compile_build {
    local root_dir="${HOME}/mt-kahypar/IFP"
    cd $root_dir || { echo "Root directory $root_dir does not exist."; exit 1; }
    local build_dir=$1
    ./build.sh $build_dir 8
    cd - || { echo "Failed to return to previous directory."; exit 1; }
}

if [ "$#" -eq 1 ]; then
    BUILD_DIR=$1
else
    BUILD_DIR=build_lib
fi
echo "Using build directory: $BUILD_DIR"


read -p "Do you want to compile the library? (y/n) " COMPILE_LIB
if [ "$COMPILE_LIB" = "y" ]; then
    read -p "Do you want to compile the build? (y/n) " COMPILE_WHOLE
    if [ "$COMPILE_WHOLE" = "y" ]; then
        echo "Compiling the build..."
        compile_build $BUILD_DIR
    fi
    echo "Compiling the library..."
    compile_library $BUILD_DIR
else
    echo "Skipping library compilation."
fi

LIB_DIR=${HOME}/mt-kahypar/IFP/${BUILD_DIR}/lib
echo "Using library directory: ${LIB_DIR}"

g++ -DNDEBUG -o partition.out -O3 -std=c++20 run_experiments_write_clusters_times_conductances.cc ${LIB_DIR}/libmtkahypar.so -I ${HOME}/mt-kahypar/IFP/include -lpthread
g++ -DNDEBUG -o conductance.out -O3 -std=c++20 analyze_results_via_conductance.cc ${LIB_DIR}/libmtkahypar.so  -I ${HOME}/mt-kahypar/IFP/include -lpthread