#!/bin/bash

for i in {1..50}
do
    for mbs in {11..15}
    do
        echo "\nRound #: $i, mRNA Burst Size: $mbs\n"
        ./sim.py \
            --results-file results.dat \
            --skip-final-plot \
            --redo-simulation \
            --timeseries-framestep 1000 \
            --n-cell 10000 \
            --mean-burst-size $mbs
    done
done
