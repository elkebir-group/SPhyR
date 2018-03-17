#!/bin/bash
if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <simulate_executable>" >&2
    exit 1
fi

for n in {25,50,100,250}
do
    for m in {50,100}
    do
	for s in {1..5}
	do
	    for loss in {0.1,0.2,0.4}
	    do
		for k in {1,2}
		do
		    for seed in {1,2,3}
		    do
			filename=m${m}_n${n}_${s}_k${k}_seed${seed}_loss${loss}
			$1 -loss $loss -k $k -s $seed -A $filename.A -B $filename.B -dot $filename.dot ../perfect_phylogeny/m${m}_n${n}_s${s}.txt
		    done
		done
	    done
	done
    done
done
