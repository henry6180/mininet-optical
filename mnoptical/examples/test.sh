#!/bin/bash
set -e

# cd .. && make clean && make dist && make force && cd examples

execfile="unilinear2.py"
sudo="sudo python3 "

max_boost_gain=30
max_numAmp=30
for boost_target_gain in $(seq 1 $max_boost_gain); do
    for numAmp in $(seq 1 $max_numAmp); do
        $sudo $execfile $boost_target_gain $numAmp<<EOF
        config
        reset
        exit
EOF
    done
done
exit 0