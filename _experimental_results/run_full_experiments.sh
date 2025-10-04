#! /bin/bash

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters. Usage: run_full_experiments.sh <prompt_name>"
    exit 1
fi

# today's date
date=$(date +%Y-%m-%d)

prompt_name=$1
echo "Starting full experiments with prompt: $prompt_name" | tee log_full_experiments.txt
# Build mt-kahypar library in build directory, build binary for experiments
./build_exp.sh build

for ip_short in AON_BAY AON SING; do
# for ip_short in AON AON_BAY SING; do
    name=""
    for vc in 3 1 0; do
        for ct in 1 2 5; do
            for rlp in ON OFF; do
                preset_name="IP=${ip_short}_VC=${vc}_CT=${ct}_RLP=${rlp}"
                if [ -n "$name" ]; then
                    name="${name},${preset_name}"
                else
                    name="${preset_name}"
                fi
                echo "Running experiments for preset: $preset_name"
                prompt=$(./prompt_preset_experiment.sh $preset_name $prompt_name)
                echo "Prompt to be used:"
                echo -e "$prompt"
                echo -e "$prompt" | ./experiments.sh log_full_experiments.txt
                echo "Finished experiments for preset: $preset_name" | tee -a log_full_experiments.txt
                echo "----------------------------------------" | tee -a log_full_experiments.txt
            done
        done
    done
    ./plot_all $name survey_benchmark ${ip_short}_${date}_${prompt_name}
done
echo "All experiments for prompt: $prompt_name completed." | tee -a log_full_experiments.txt
echo "----------------------------------------" | tee -a log_full_experiments.txt