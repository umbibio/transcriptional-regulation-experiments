#!/bin/bash

for i in {1..500}
do
    for mbs in {1..100}
    do
        ./sim.py \
            --results-file results.dat \
            --skip-final-plot \
            --redo-simulation \
            --timeseries-framestep 1000 \
            --n-cell 10000 \
            --mean-burst-size $mbs \
            > /dev/null
    done
done
