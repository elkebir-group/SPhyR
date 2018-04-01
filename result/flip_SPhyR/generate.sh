#!/bin/bash
if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <kDPFC_executable>" >&2
    exit 1
fi

exec="~/src/SPhyR/build/kDPFC"
data="~/src/SPhyR/data/flip/"
result="~/src/SPhyR/result/flip/"

k=1
n=50
m=50
alpha=0.001
beta=0.2
lT=10
lC=35

for input_alpha in {0.0001,0.001,0.01}
do
    for input_beta in {0.1,0.2,0.3}
    do
        for loss in {0.1,0.2,0.4}
        do
            for input_k in {0,1,2}
            do
                pbs_filename=swp_a${input_alpha}_b${input_beta}_loss${loss}_k${input_k}.pbs
                sh_filename=swp_a${input_alpha}_b${input_beta}_loss${loss}_k${input_k}.sh
                echo "#PBS -l nodes=1:ppn=20" > $pbs_filename
                echo "#PBS -l walltime=10:00:00" >> $pbs_filename
                echo "#PBS -N a${input_alpha}_b${input_beta}_loss${loss}_k${input_k}" >> $pbs_filename
                echo "#PBS -e a${input_alpha}_b${input_beta}_loss${loss}_k${input_k}.err" >> $pbs_filename
                echo "#PBS -o a${input_alpha}_b${input_beta}_loss${loss}_k${input_k}.out" >> $pbs_filename
                echo "cd ~/scratch/sphyr" >> $pbs_filename
                echo "aprun -n 1 -d 20 /bin/bash $result/$sh_filename" >> $pbs_filename
                echo "#!/bin/bash" > $sh_filename
                for s in {1..20}
                do
                    filename=m${m}_n${n}_s${s}_k${k}_loss${loss}_a${alpha}_b${beta}
                    out_filename=${filename}_inA${input_alpha}_inB${input_beta}_inK${input_k}
                    echo "($exec -k $input_k -T 25000 -lC $lC -lT $lT -N 100 -a $input_alpha -b $input_beta $data/${filename}.B ${out_filename}.A 2> ${out_filename}.log && mv ${out_filename}.* $result/) &" >> $sh_filename
                done
                echo "wait" >> $sh_filename
                chmod a+x $sh_filename
            done
        done
    done
done

