#!/bin/bash
set -e

execfile="unilinear2.py"
simfile="simulate.py"
sudo="sudo python3 "

length=100
roadm_insertion_loss=17
numAmp=2
boost_target_gain=17
if [ $# -ge 1 ]; then
    length=$(($1))
fi
if [ $# -ge 2 ]; then
    roadm_insertion_loss=$(($2))
fi
if [ $# -ge 3 ]; then
    numAmp=$(($3))
fi
if [ $# -ge 4 ]; then
    boost_target_gain=$(($4))
fi

$sudo $execfile $length $roadm_insertion_loss $numAmp $boost_target_gain<<EOF
    config
    reset
    exit
EOF
$sudo $simfile $length $roadm_insertion_loss
exit 0