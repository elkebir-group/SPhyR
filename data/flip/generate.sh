#!/bin/bash
if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <perturb_executable>" >&2
    exit 1
fi

m=50
alpha=0.001
for n in {25,50,100}
do
    for s in {1..5}
    do
        for seed in {1,2,3}
        do
            for k in {1,2}
            do
                for loss in {0.1,0.2,0.4}
                do
                    for beta in {0.1,0.2,0.3}
                    do
                        filename=m${m}_n${n}_${s}_k${k}_seed${seed}_loss${loss}
                        $1 -a $alpha -b $beta ../k_dollo/$filename.B ${filename}_a${alpha}_b${beta}.B
                    done
                done
            done
        done
    done
done
