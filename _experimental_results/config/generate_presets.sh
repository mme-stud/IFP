#! /bin/bash

# Generate altered presets for experiments
# Base preset is dummy_preset.ini
# New presets are generated from it by making following changes:

# Type one of the change: change IP:
# 1. Use aon_hypermodularity: [IP=AON]
# Was: i-enabled-ip-algos=0 # aon_hypermodularity
# Will be: i-enabled-ip-algos=1 # aon_hypermodularity
# 2. Use aon_hypermodularity_bayesian: [IP=AON_BAY]
# Was: i-enabled-ip-algos=0 # aon_hypermodularity_bayesian
# Will be: i-enabled-ip-algos=1 # aon_hypermodularity_bay
# 3. Use singleton: [IP=SING]
# Was: i-enabled-ip-algos=0 # singleton
# Will be: i-enabled-ip-algos=1 # singleton

# Type two of the change: change number of V-cycles:
# 1. Use 0 V-cycles: [VC=0]
# Was: num-vcycles=0
# Will be: num-vcycles=0
# 2. Use 1 V-cycle: [VC=1]
# Was: num-vcycles=0
# Will be: num-vcycles=1
# 3. Use 3 V-cycles: [VC=3]
# Was: num-vcycles=0
# Will be: num-vcycles=3

# Type three of the change: change c-t
# 1. Use c-t=5: [CT=5]
# Was: c-t=20
# Will be: c-t=5
# 2. Use c-t=20: [CT=20]
# Was: c-t=20
# Will be: c-t=20
# 3. Use c-t=80: [CT=80]
# Was: c-t=20
# Will be: c-t=80

# Type four of the change: use i-r-lp
# 1. Use r-lp:[RLP=ON]
# Was: 
# i-r-lp-type=do_nothing
# i-r-lp-maximum-iterations=0
# Will be: 
# i-r-lp-type=label_propagation
# i-r-lp-maximum-iterations=5
# 2. Use no r-lp: [RLP=OFF]
# Was: 
# i-r-lp-type=do_nothing
# i-r-lp-maximum-iterations=0
# Will be: 
# i-r-lp-type=do_nothing
# i-r-lp-maximum-iterations=0

function change_ip {
    local preset_file=$1
    local new_ip=$2
    # new_ip is one of aon_hypermodularity, aon_hypermodularity_bayesian, singleton
    # don't ignore # in the line
    sed -i "s/i-enabled-ip-algos=0 # $new_ip$/i-enabled-ip-algos=1 # $new_ip/" $preset_file
}

function change_vc {
    local preset_file=$1
    local new_vc=$2
    # new_vc is one of 0, 1, 3
    sed -i "s/num-vcycles=[0-9]*/num-vcycles=$new_vc/" $preset_file
}

function change_ct {
    local preset_file=$1
    local new_ct=$2
    # new_ct is one of 5, 20, 80
    sed -i "s/c-t=[0-9]*/c-t=$new_ct/" $preset_file
}

function change_rlp {
    local preset_file=$1
    local new_rlp=$2
    # new_rlp is one of ON, OFF
    if [ "$new_rlp" == "ON" ]; then
        sed -i "s/i-r-lp-type=do_nothing/i-r-lp-type=label_propagation/" $preset_file
        sed -i "s/i-r-lp-maximum-iterations=0/i-r-lp-maximum-iterations=5/" $preset_file
    else
        sed -i "s/i-r-lp-type=label_propagation/i-r-lp-type=do_nothing/" $preset_file
        sed -i "s/i-r-lp-maximum-iterations=[0-9]*/i-r-lp-maximum-iterations=0/" $preset_file
    fi
}

for ip in aon_hypermodularity aon_hypermodularity_bayesian singleton; do
    for vc in 0 1 3; do
        for ct in 1 2 5; do
            for rlp in ON OFF; do
                ip_short=""
                if [ $ip == aon_hypermodularity ]; then
                    ip_short=AON
                elif [ $ip == aon_hypermodularity_bayesian ]; then
                    ip_short=AON_BAY
                else
                    ip_short=SING
                fi
                new_preset="IP=${ip_short}_VC=${vc}_CT=${ct}_RLP=${rlp}_preset.ini"
                # remove existing file
                if [ -f config/$new_preset ]; then
                    rm config/$new_preset
                fi
                cp dummy_preset.ini $new_preset
                change_ip $new_preset $ip
                change_vc $new_preset $vc
                change_ct $new_preset $ct
                change_rlp $new_preset $rlp
                echo "Generated preset: $new_preset"
            done
        done
    done
done