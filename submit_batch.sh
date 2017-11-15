#!/bin/bash

# experiment parameters
MBS=$1
BDIST="geometric"
KB=10
UM=2
KP=15
UP=0.05
EXPD=150

# jobs parameters
HL=180      # 3 hours
NCPU=16

# load proper modules
module load python/2.7.9_packages/matplotlib/1.4.3


if [ -z $MBS ] || [ ! $MBS -gt 0 ]; then
    echo "ERROR - Need to specify first argument: Mean Burst Size"
    exit 1
fi

for i in {1..5}
do
    bsub -n $NCPU -R rusage[mem=64] -R span[hosts=1]\
        -W $HL -q short -J $(printf "%02d" $MBS)${BDIST:0:2}i$i \
        python sim.py \
        -kb $KB -um $UM -kp $KP -up $UP -mbs $MBS \
        -n 20000 -bsd $BDIST -d $EXPD \
        -nfp -q -nj $NCPU -i $i \
        -rf results_kb${KB}-um${UM}-kp${KP}-up${UP}.dat \
        -l $((($HL-3)*60))
done
