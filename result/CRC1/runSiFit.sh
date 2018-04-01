#!/bin/bash
if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <SiFit.jar>" >&2
    exit 1
fi

exec=$1
fp=0.0152
java -jar $exec -r 100 -m 178 -n 16 -fp $fp -fn 0.0789 -iter 10000 -df 0 -ipMat ../../data/CRC/CRC1.SiFit.input > CRC1.log
fn=$(tail -n 5 CRC1.log | head -n 1 | cut -d' ' -f 6)
deletion=$(tail -n 4 CRC1.log | head -n 1 | cut -d' ' -f 5)
LOH=$(tail -n 3 CRC1.log | head -n 1 | cut -d' ' -f 5)

echo $fn
echo $deletion
echo $LOH

java -cp $exec SiFit.algorithm.InferAncestralStates -fp $fp -fn $fn -w $LOH -d $deletion -df 0 -ipMat ../../data/CRC/CRC1.SiFit.input -tree CRC1.SiFit.input_mlTree.newick -geneNames CRC1_SNV.SiFit.labels -cellNames CRC1_cell.SiFit.labels -expectedMatrix CRC1.leaves > CRC1.tree
python ../parseSiFit.py CRC1.tree CRC1.leaves 16 > CRC1.SPhyR.tree

