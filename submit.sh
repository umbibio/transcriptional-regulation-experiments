#!/bin/bash

#BSUB -n 4
#BSUB -R rusage[mem=512] # this is per job slot
#BSUB -W 2:00
#BSUB -q short # which queue we want to run in
#BSUB -R "span[hosts=1]" # All job slots on the same node (needed for threaded applications)

module load python/2.7.9_packages/matplotlib/1.4.3

python sim.py -kb 10 -um 2 -kp 15 -up 0.05 -mbs 4 -n 20000 -bsd geometric -d 150 -nfp -q -nj 4 -i 1
