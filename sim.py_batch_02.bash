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

for mbs in {1..5}
do
    for i in {1..5}
    do
        echo -e "\nmRNA Burst Size: $mbs, Round #: $i\n"
        ./sim.py \
            --results-file results.dat \
            --skip-final-plot \
            --redo-simulation \
            --experiment-duration 5 \
            --n-cell 100000 \
            --burst-size-distribution conditional_geometric \
            --mean-burst-size $mbs
    done
done
