#!/bin/bash
for f in *.log
do
    grep "elapsed time" $f | tail -n 1 | cut -d' ' -f 6 > $(basename $f .log).time
done
