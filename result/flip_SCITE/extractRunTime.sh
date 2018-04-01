#!/bin/bash
for f in *inA*.log
do
    (echo "$(tail -n 1 $f | cut -d' ' -f 3) / 1000" | bc) > $(basename $f .log).time
done
