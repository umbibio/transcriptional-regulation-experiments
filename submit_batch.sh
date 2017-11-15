#!/bin/bash

MBS=$1
BDIST="geometric"

if [ -z $MBS ] || [ ! $MBS -gt 0 ]; then
    echo "ERROR - Need to specify first argument: Mean Burst Size"
    exit 1
fi

for i in {1..5}
do
    bsub -n 8 -R rusage[mem=128] -R span[hosts=1]\
        -W 2:00 -q short -J $(printf "%02d" $MBS)${BDIST:0:2}$i \
        python sim.py \
        -kb 10 -um 2 -kp 15 -up 0.05 -mbs $MBS \
        -n 20000 -bsd $BDIST -d 150 \
        -nfp -q -nj 8 -i $i \
        -rf results.dat
done