#! /bin/bash

if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters. Usage: prompt_preset_experiment.sh  <preset_name> <experiment_prompt>"
    echo "Example: prompt_preset_experiment.sh dummy prompt_all"
    exit 1
fi

preset_name=$1
experiment_prompt=$2

if [ ! -f prompts/${experiment_prompt}.sh ]; then
    echo "Preset config file prompts/${experiment_prompt}.sh doesn't exist"
    exit 1
fi

if [ ! -f config/${preset_name}_preset.ini ]; then
    echo "Preset config file config/${preset_name}.ini doesn't exist"
    exit 1
fi

# version name
#echo $preset_name
# preset type
#echo $preset_name
# rest
#./prompts/${experiment_prompt}.sh

rest=$(./prompts/${experiment_prompt}.sh)

echo -e "$preset_name\n$preset_name\n$rest"