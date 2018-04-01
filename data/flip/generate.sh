#!/bin/bash
if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <perturb_executable>" >&2
    exit 1
fi

m=50
n=50
alpha=0.001
beta=0.2
k=1
for loss in {0.1,0.2,0.4}
do
    filename=m${m}_n${n}_s${s}_k${k}_loss${loss}
    $1 -a $alpha -b $beta ../k_dollo/$filename.B ${filename}_a${alpha}_b${beta}.B
done
