#!/bin/bash
set -e
# Usage: ./test2.sh length roadm_insertion_loss numAmp boost_target_gain ber

execfile="unilinear2.py"
simfile="simulate.py"
sudo="sudo python3 "

length=100
roadm_insertion_loss=17
numAmp=2
boost_target_gain=17
ber="1"
if [ $# -ge 1 ]; then
    length=$1
fi
if [ $# -ge 2 ]; then
    roadm_insertion_loss=$2
fi
if [ $# -ge 3 ]; then
    numAmp=$3
fi
if [ $# -ge 4 ]; then
    boost_target_gain=$4
fi
if [ $# -ge 5 ]; then
    ber=$5
fi

# input of execfile is length roadm_insertion_loss numAmp boost_target_gain
# input of simfile is length roadm_insertion_loss ber

# output format of execfile in result1.txt is boost_target_gain numAmp t1-gosnr t2-gosnr t1-ber t2-ber
# output of simfile in result2.txt are those cases whose t1-ber and t2-ber are both less than ber and
# the format of result2.txt is the same with result1.txt

# In result2.txt t1-ber(or t2-ber) == 1 if its gosnr is negative(such that out of the domain of get_ber())


$sudo $execfile $length $roadm_insertion_loss $numAmp $boost_target_gain<<EOF
    config
    reset
    exit
EOF
$sudo $simfile $length $roadm_insertion_loss $ber
exit 0