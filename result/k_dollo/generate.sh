#!/bin/bash
if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <kDP_executable>" >&2
    exit 1
fi

data="../../data/k_dollo/"
for n in {25,50,75,100,250}
do
    for m in {50,100}
    do
	for s in {1..20}
	do
	    for loss in {0.1,0.2,0.4}
	    do
		for k in {1,2}
		do
		    filename=m${m}_n${n}_s${s}_k${k}_seed${seed}_loss${loss}
		    echo Running $filename...
		    #$1 -k $k -t 1 -T 1200 ${data}/${filename}.B ${filename}.A -v 2> ${filename}.log
		done
	    done
	done
    done
done

