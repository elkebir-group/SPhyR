#!/bin/bash

for f in *inA*.log
do
    tail -n 1 $f | cut -d' ' -f 4 > $(basename $f .log)_mlTree.newick
done
