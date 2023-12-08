#!/bin/bash
set -e

# cd .. && make clean && make dist && make force && cd examples
# directory="/home/henry/mininet-optical/examples/"
execfile="unilinear2.py"
# path=$directory$execfile
sudo="sudo python3 "

for boost_target_gain in 19 20 21 22 23 24 25 26 27 28 29; do
    for numAmp in 1 2 3 4 5 6 7 8; do
        $sudo $execfile $boost_target_gain $numAmp<<EOF
        config
        reset
        exit
EOF
    done
done