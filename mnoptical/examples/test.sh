#!/bin/bash
set -e

# cd .. && make clean && make dist && make force && cd examples
execfile="unilinear2.py"
sudo="sudo python3 "

length=100
roadm_insertion_loss=17
if [ $# -ge 1 ]; then
    length=$(($1))
fi
if [ $# -ge 2 ]; then
    roadm_insertion_loss=$(($2))
fi

minB=$(($roadm_insertion_loss-5))
maxB=$(($roadm_insertion_loss+5))
max_numAmp=20

for boost_target_gain in $(seq $minB $maxB); do
    for numAmp in $(seq 1 $max_numAmp); do
        $sudo $execfile $length $roadm_insertion_loss $numAmp $boost_target_gain<<EOF
        config
        reset
        exit
EOF
    done
done

exit 0