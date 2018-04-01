#!/bin/bash
if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <SPhyR>" >&2
    exit 1
fi

exec=$1
alpha=0.0152
beta=0.0789
s=10
t=15

$exec ../../data/CRC/CRC1.SPhyR.input -k 0 -a $alpha -b $beta -s $s -t $t > SPhyR_k0.A
$exec ../../data/CRC/CRC1.SPhyR.input -k 1 -a $alpha -b $beta -s $s -t $t > SPhyR_k1.A
