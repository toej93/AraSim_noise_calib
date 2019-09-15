#!/bin/bash

declare -a arr=("4.5" "6.0" "6.2")
#declare -a arr=("5.70" "5.72" "5.74" "5.76" "5.78" "5.80" "5.82" "5.84" "5.86" "5.88" "5.90")

#declare -a arr=("5.8")

for i in "${arr[@]}"
do
    AraSimDir='/users/PAS0654/osu8354/ARA/AraSim/trunk/'
    SetUpFile='/users/PAS0654/osu8354/ARA/AraSim/trunk/SETUP/setup_'$i'_th.txt'
    #SetUpFile='/users/PAS0654/osu8354/ARA/AraSim/trunk/setup_'$i'_th_stdnoise.txt'
    #OutputDir='/users/PAS0654/osu8354/ARA/AraSim/trunk/outputs/'
    OutputDir='/fs/scratch/PAS0654/jorge/sim_results/'
    export AraSimDir
    export SetUpFile
    export OutputDir
    
    echo ""
    echo "--------------------------------------------"
    echo "----- Preparing to batch submit AraSim "
    echo "----- "
    echo "----- AraSimDir: " $AraSimDir
    echo "----- SetUpFile: " $SetUpFile
    echo "----- outputDir: " $OutputDir
    echo "----- "
    echo "--------------------------------------------"
    echo ""
    
    #	read -p "Press Enter to run! " RUNNOW
    
    qsub -v INPUTFILE=$SetUpFile,RUN_DIR=$AraSimDir,OUTPUT_DIR=$OutputDir -N AraSim_noise_"$i" run_desA_e17.sh
    
    
    #sed -i '$ a RAY_TRACE_ICE_MODEL_PARAMS=2' des* #append to last line
    #sed -i 's/old-text/new-text/g' input.txt
    
done
