#!/bin/bash

delayUntil=${1:-"None"}

if [ ! $delayUntil == "None" ];
then
    startTime=$(date +%s)
    endTime=$(date -d "$delayUntil" +%s)
    if [ $? -ne 0 ];
    then
        echo "Time provided is not valid"
        exit 1
    fi
    timeToWait=$(($endTime- $startTime))
    sleep $timeToWait
fi

for i in {1..50}
do
    for mbs in {11..15}
    do
        echo -e "\nRound #: $i, mRNA Burst Size: $mbs\n"
        ./sim.py \
            --results-file results.dat \
            --skip-final-plot \
            --redo-simulation \
            --timeseries-framestep 1000 \
            --n-cell 10000 \
            --mean-burst-size $mbs
    done
done
