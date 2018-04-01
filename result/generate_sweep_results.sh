#!/bin/bash
if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <ANALYZE_EXEC>" >&2
    exit 1
fi

echo "method,m,n,s,k,loss,alpha,beta,time,RF,norm_RF,anc_loss_recall,inc_loss_recall,cls_loss_recall,anc_recall,inc_recall,cls_recall,taxa_RI,taxa_recall,taxa_precision,char_RI,char_recall,char_precision,L,back_mut_inf,par_evo_inf,back_mut_true,par_evo_true,loss_recall,loss_precision,loss_F1,flip01_correct,flip01_incorrect,flip10_correct,flip10_incorrect,TN,FN,FP,TP,ones,zeros,inTN,inFN,inFP,inTP"
#for method in {SPhyR_lT10_lC10,SPhyR_lT10_lC15,SPhyR_lT10_lC20,SPhyR_lT15_lC10,SPhyR_lT15_lC15,SPhyR_lT15_lC20,SCITE,SiFit}
for method in {SCITE,SiFit}
do
    ./processSweep.sh $method $1
done
for k in {0,1,2}
do
    ./processSweep.sh SPhyR_inK${k} $1
done
