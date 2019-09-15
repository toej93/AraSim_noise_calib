#/bin/bash
#PBS -l nodes=1:ppn=28
#PBS -j oe
#PBS -A PAS0654
#PBS -o logs_pbs/
#PBS -m n
#PBS -l mem=128000MB
#PBS -l walltime=60:20:00

source /users/PCON0003/cond0068/.bash_profile_pitzer
#module load modules/au2016

cd $RUN_DIR

j=101
while [ $j -lt 201 ]
do
    END=$[$j+28]
    for i in $(seq $j $END) #3785
    do
        echo $i
        ./AraSim $INPUTFILE $i $TMPDIR &
    done
    wait
    j=$[$j+28]
    pbsdcp $TMPDIR/'*' $OUTPUT_DIR
done


