#!/bin/bash
for f in flip_SCITE/*ml0.gv
do
    g=$(basename $f)
    m=`echo $g | sed -e "s/^m\([0-9]*\)_n\([0-9]*\)_s\([0-9]*\)_k\([0-9]*\)_loss\(0\.[0-9]*\)_a\(0\.[0-9]*\)_b\(0\.[0-9]*\).*/\1/g"`
    n=`echo $g | sed -e "s/^m\([0-9]*\)_n\([0-9]*\)_s\([0-9]*\)_k\([0-9]*\)_loss\(0\.[0-9]*\)_a\(0\.[0-9]*\)_b\(0\.[0-9]*\).*/\2/g"`
    python processSCITE.py $f $m $n > $(dirname $f)/$(basename $f _ml0.gv).A
done
