#!/bin/bash
set -e
#parameter
directory="/home/henry/mininet-optical/examples"
execfile="unilinear2.py"
cmd="python3 "

for boost_target_gain in 16 17 18; do
    for numAmp in 1 2 3 4 5 6 7 8 9 10; do
        sudo $cmd $directory$execfile $boost_target_gain $numAmp<<EOF
        config
        reset
        exit
EOF
    done
done