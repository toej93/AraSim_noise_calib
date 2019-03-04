#/bin/bash
#PBS -l nodes=1:ppn=40
#PBS -j oe
#PBS -A PCON0003
#PBS -m e
#PBS -l mem=128000MB
#PBS -l walltime=60:20:00

source /users/PCON0003/cond0068/.bash_profile_pitzer
#module load modules/au2016

cd $RUN_DIR

j=600
while [ $j -lt 1000 ]
do
    END=$[$j+40]
    for i in $(seq $j $END) #3785
    do
        echo $i
        ./AraSim $INPUTFILE $i $TMPDIR &
    done
    wait
    j=$[$j+40]
    pbsdcp $TMPDIR/'*' $OUTPUT_DIR
done


