#/bin/bash
#PBS -N RayleighA2
#PBS -l nodes=1:ppn=3
#PBS -j oe
#PBS -A PCON0003
#PBS -m e
#PBS -l mem=64000MB
#PBS -l walltime=02:00:00

source /users/PCON0003/cond0068/.bash_profile_pitzer

cd /users/PAS0654/osu8354/ARA/AraSim/trunk/

for i in {0..2}
do
    ./AraSim setup.txt &
done
wait
