#!/bin/bash

LOG_OUTPUT=n
LOG_FILE=log_experiments.txt
if [ "$#" -eq 1 ]; then
    LOG_OUTPUT=y
    LOG_FILE=$1
    echo "Logging output to $LOG_FILE"
fi

############################### Functions #################################
################# Running functions #################

function run_scenario {
    local scenario=$1
    local mode=$2

    if [ "$ONLY_AON_HMLL" = "n" ]; then
        echo "Running scenario: $scenario"
        ./run_exp.sh $BUILD_DIR $PRESET $VERSION $mode $scenario $NUM_THREADS $GRAPHS
        echo "Experiments for $scenario completed."

        ./analyze_exp.sh $BUILD_DIR $PRESET $VERSION survey_benchmark $mode $scenario $NUM_THREADS $GRAPHS
        echo "Analysis for $scenario completed."
    fi

    if [ "$RUN_AON_HMLL" = "y" ]; then
        echo "Running AON_HMLL and analyzing its results"
        ./run_AON_HMLL_experiments.sh $BUILD_DIR $PRESET $AON_VERSION survey_benchmark $mode $scenario $GRAPHS
        echo "Analysis for AON_HMLL in $scenario completed."
    fi
}

# Loop through scenarios in the directory of the specified mode. Scenarios are the names of subdirectories.

function run_mode {
    local mode=$1

    local scenario_dir="survey_benchmark/${mode}"
    echo "Running $GRAPHS from each scenario in mode: $mode"
    for scenario in $(ls $scenario_dir); do
        if [[ $scenario != scen* ]]; then
            continue
        fi
        run_scenario $scenario $mode
        if [ "$LOG_OUTPUT" = "y" ]; then
            echo "Finished scenario: $scenario" | tee -a $LOG_FILE
        fi
    done
}

function run_all_modes {

    local modes=("DCHSBM" "hABCD" "HyperSBM")
    for mode in "${modes[@]}"; do
        if [ "$LOG_OUTPUT" = "y" ]; then
            echo "Running mode $mode:" | tee -a $LOG_FILE
        fi
        run_mode $mode
        if [ "$LOG_OUTPUT" = "y" ]; then
            echo "----------------------------------------" | tee -a $LOG_FILE
        fi
    done

    echo "All modes completed."

}

########### User Input Functions ###########

function set_graphs {
    local use_preset
    local preset_name
    read -p "Use graphs from a preset? (y/n, default: n): " use_preset
    use_preset=${use_preset:-n}
    if [ "$use_preset" != "n" ]; then
        if [ "$use_preset" = "y" ]; then
            read -p "Enter preset name (e.g., all, test, experiment): " preset_name
        else
            preset_name=${use_preset}
            echo "Using preset name: $preset_name"
        fi
        if [ -f "presets/${preset_name}_graphs.txt" ]; then
            GRAPHS=$(<presets/${preset_name}_graphs.txt)
            echo "Using graphs from preset '${preset_name}': $GRAPHS"
        else
            echo "Preset file 'presets/${preset_name}_graphs.txt' not found."
            use_preset="n"
        fi
    fi
    if [ "$use_preset" = "n" ]; then
        read -p "Enter ids of graphs (comma-separated string): " GRAPHS
        GRAPHS=${GRAPHS}
    fi
}

function set_mode {
    read -p "Enter mode (1=DCHSBM, 2=hABCD, 3=HyperSBM): " MODE
    if [[ "$MODE" == "1" ]]; then
        MODE="DCHSBM"
    elif [[ "$MODE" == "2" ]]; then
        MODE="hABCD"
    elif [[ "$MODE" == "3" ]]; then
        MODE="HyperSBM"
    else
        echo "Invalid mode selection."
        exit 1
    fi
}

function set_scenario {
    local all_scenarios
    read -p "Do you want to run all scenarios? (y/n): " all_scenarios
    if [ "$all_scenarios" = "y" ]; then
        SCENARIO=all
    elif [ "$all_scenarios" = "n" ]; then
        read -p "Enter scenario (e.g., scenA1): " SCENARIO
    else
        SCENARIO=${all_scenarios}
        echo "Using scenario: $SCENARIO"
    fi
    echo $SCENARIO
}

################################## Input ##################################


read -p "Enter version name (default: test): " VERSION
VERSION=${VERSION:-test}

read -p "Enter preset name (default: cluster): " PRESET
PRESET=${PRESET:-cluster}

read -p "Enter build directory (default: build): build" BUILD_DIR
BUILD_DIR=build${BUILD_DIR}

GRAPHS=1,2,3
set_graphs

read -p "Enter number of threads (default: 3): " NUM_THREADS
NUM_THREADS=${NUM_THREADS:-3}

# Select mode and scenario
read -p "Do you want to run all modes and scenarios? (y/n, default: y): " MODE
MODE=${MODE:-y}
if [ "$MODE" = "y" ]; then
    MODE=all
    SCENARIO=all
elif [ "$MODE" = "n" ]; then
    set_mode
    set_scenario
else
    echo "Using input as mode: $MODE"
    set_scenario
fi

read -p "Do you want to also run AON_HMLL for comparison? (y/n, default: n): " RUN_AON_HMLL
RUN_AON_HMLL=${RUN_AON_HMLL:-n}
AON_VERSION="AON_HMLL_clique"
ONLY_AON_HMLL="n"
if [ "$RUN_AON_HMLL" = "y" ]; then
    read -p "Do you want to ONLY run AON_HMLL and skip main experiments? (y/n, default: n): " ONLY_AON_HMLL
    ONLY_AON_HMLL=${ONLY_AON_HMLL:-n}
fi

read -p "Do you want to plot combined results? (y/n, default: y): " PLOT_COMBINED
PLOT_COMBINED=${PLOT_COMBINED:-y}
if [ "$PLOT_COMBINED" = "y" ]; then
    read -p "Enter versions to compare with (comma-separated): $(VERSION)" VERSIONS_TO_PLOT
    VERSIONS_TO_PLOT=$VERSION$VERSIONS_TO_PLOT
fi

################################## Build ##################################
if [ "$ONLY_AON_HMLL" != "y" ]; then
    ./build_exp.sh $BUILD_DIR
fi

############################# Run Experiments #############################

# Execute experiments based on selections
if [ "$MODE" = "all" ]; then
    run_all_modes
    echo "All modes and scenarios completed."
elif [ "$SCENARIO" = "all" ]; then
    run_mode $MODE
    echo "All scenarios completed."
else 
    run_scenario $SCENARIO $MODE
fi

echo "All selected experiments completed."
############################# Plot Combined Results #############################
if [ "$PLOT_COMBINED" = "y" ]; then
    echo "Plotting combined results for versions: $VERSIONS_TO_PLOT"
    ./plot_all.sh $VERSIONS_TO_PLOT survey_benchmark
    echo "Plotting completed."
fi
