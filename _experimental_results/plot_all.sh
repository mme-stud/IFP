#! /bin/bash

NAME=""

if [ "$#" -eq 3 ]; then
    NAME=$3
elif [ "$#" -lt 2 ]; then
    echo "Illegal number of parameters. Usage: plot_all.sh <comma-separated versions> <data_root_dir>"
    exit 1
fi

VERSIONS=$1
DATA_ROOT_DIR=$2
python3 plot_all_data.py --versions $VERSIONS --data_root_dir $DATA_ROOT_DIR --output $NAME