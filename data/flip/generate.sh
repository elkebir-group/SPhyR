#!/bin/bash
if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <perturb_executable>" >&2
    exit 1
fi

k=1
n=50
m=50
loss=0.2
for s in {1..5}
do
    for seed in {1,2,3}
    do
        for alpha in {0.001,0.01,0.05}
        do
            for beta in {0.1,0.2,0.3}
            do
                filename=m${m}_n${n}_${s}_k${k}_seed${seed}_loss${loss}
                $1 -a $alpha -b $beta ../k_dollo/$filename.B ${filename}_a${alpha}_b${beta}.B
            done
        done
    done
done

