#/bin/bash
#PBS -N desA_e20
#PBS -l nodes=1:ppn=28
#PBS -j oe
#PBS -A PAS0654
#PBS -m e
#PBS -l mem=64000MB
#PBS -l walltime=10:20:00

source /users/PCON0003/cond0068/.bash_profile_owens
module load modules/au2016

cd $RUN_DIR

for i in {1..25}
do
	./AraSim $INPUTFILE $i $OUTPUT_DIR & 
done
wait

for i in {25..50}
do
	./AraSim $INPUTFILE $i $OUTPUT_DIR & 
done
wait

for i in {51..75}
do
	./AraSim $INPUTFILE $i $OUTPUT_DIR & 
done
wait

for i in {76..100}
do
	./AraSim $INPUTFILE $i $OUTPUT_DIR & 
done
wait

#pbsdcp $TMPDIR/'*' $OUTPUT_DIR
