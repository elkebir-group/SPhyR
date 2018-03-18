#!/bin/bash

if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <ms_executable>" >&2
    echo "To generate the same files, make sure that ms is compiled using rand1.c" >&2
    exit 1
fi

seeds="7369 217"

for n in {25,50,75,100,250,500}
do
    for m in {50,100,150,200}
    do
	for s in {1..20}
	do
	    filename=m${m}_n${n}_s${s}.txt
	    echo "$m # cells" > $filename
	    echo "$n # mutations" >> $filename
	    $1 $m 1 -s $n -seeds $seeds $s | tail -n $m >> $filename
	    sed -i -e '3,$s/\([0-9]\)/\1 /g' $filename
	    sed -i -e '3,$s/ $//g' $filename
	done
    done
done

rm *.txt-e
