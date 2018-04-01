#!/bin/bash
if [ ! $# -eq 2 ]
then
    echo "Usage: $0 <SiFit.jar> <DIR>" >&2
    exit 1
fi

#java -cp ~/Projects/sifit/SiFit.jar SiFit.algorithm.InferAncestralStates -fp 0.001 -fn 0.095746 -w 0.966823 -d 0.706015 -df 0 -ipMat ../data/flip/m50_n50_1_k1_seed1_loss0.2_a0.001_b0.1.SiFit -tree flip_SiFit/m50_n50_1_k1_seed1_loss0.2_a0.001_b0.1.SiFit_mlTree.newick -geneNames SiFit.characters.txt -cellNames SiFit.taxa.txt -expectedMatrix test.txt

#flip_SiFit/m50_n50_1_k1_seed1_loss0.2_a0.001_b0.1.SiFit.log

for f in flip_${2}/*.log
do
    g=`basename $f`
    m=`echo $g | sed -e "s/^m\([0-9]*\)_n\([0-9]*\)_s\([0-9]*\)_k\([0-9]*\)_loss\(0\.[0-9]*\)_a\(0\.[0-9]*\)_b\(0\.[0-9]*\).*/\1/g"`
    n=`echo $g | sed -e "s/^m\([0-9]*\)_n\([0-9]*\)_s\([0-9]*\)_k\([0-9]*\)_loss\(0\.[0-9]*\)_a\(0\.[0-9]*\)_b\(0\.[0-9]*\).*/\2/g"`
    s=`echo $g | sed -e "s/^m\([0-9]*\)_n\([0-9]*\)_s\([0-9]*\)_k\([0-9]*\)_loss\(0\.[0-9]*\)_a\(0\.[0-9]*\)_b\(0\.[0-9]*\).*/\3/g"`
    k=`echo $g | sed -e "s/^m\([0-9]*\)_n\([0-9]*\)_s\([0-9]*\)_k\([0-9]*\)_loss\(0\.[0-9]*\)_a\(0\.[0-9]*\)_b\(0\.[0-9]*\).*/\4/g"`
    loss=`echo $g | sed -e "s/^m\([0-9]*\)_n\([0-9]*\)_s\([0-9]*\)_k\([0-9]*\)_loss\(0\.[0-9]*\)_a\(0\.[0-9]*\)_b\(0\.[0-9]*\).*/\5/g"`
    a=`echo $g | sed -e "s/^m\([0-9]*\)_n\([0-9]*\)_s\([0-9]*\)_k\([0-9]*\)_loss\(0\.[0-9]*\)_a\(0\.[0-9]*\)_b\(0\.[0-9]*\).*/\6/g"`
    b=`echo $g | sed -e "s/^m\([0-9]*\)_n\([0-9]*\)_s\([0-9]*\)_k\([0-9]*\)_loss\(0\.[0-9]*\)_a\(0\.[0-9]*\)_b\(0\.[0-9]*\).*/\7/g"`
    inA=`echo $g | sed -e "s/^m\([0-9]*\)_n\([0-9]*\)_s\([0-9]*\)_k\([0-9]*\)_loss\(0\.[0-9]*\)_a\(0\.[0-9]*\)_b\(0\.[0-9]*\).SiFit_inA\(0\.[0-9]*\)_inB\(0\.[0-9]*\).*/\8/g"`
    inB=`echo $g | sed -e "s/^m\([0-9]*\)_n\([0-9]*\)_s\([0-9]*\)_k\([0-9]*\)_loss\(0\.[0-9]*\)_a\(0\.[0-9]*\)_b\(0\.[0-9]*\).SiFit_inA\(0\.[0-9]*\)_inB\(0\.[0-9]*\).*/\9/g"`
    fp=$inA
    fn=$(tail -n 5 $f | head -n 1 | cut -d' ' -f 6)
    deletion=$(tail -n 4 $f | head -n 1 | cut -d' ' -f 5)
    LOH=$(tail -n 3 $f | head -n 1 | cut -d' ' -f 5)

    input=m${m}_n${n}_s${s}_k${k}_loss${loss}_a0.001_b0.2.SiFit

    java -cp $1 SiFit.algorithm.InferAncestralStates -fp $fp -fn $fn -w $LOH -d $deletion -df 0 -ipMat ../data/flip/${input} -tree $(dirname $f)/$(basename $g .log)_mlTree.newick -geneNames SiFit.characters.txt -cellNames SiFit.taxa.txt -expectedMatrix $(dirname $f)/$(basename $g .log).leaves > $(dirname $f)/$(basename $g .log).tree
    python parseSiFit.py $(dirname $f)/$(basename $g .log).tree $(dirname $f)/$(basename $g .log).leaves $n > $(dirname $f)/$(basename $g .log).SPhyR.tree
done
