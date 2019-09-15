#!/bin/bash
declare -a var_name

#for i in $(seq 5.7 0.02 5.9)

declare -a thres=("4.5" "5.0" "5.5" "6.0" "6.5")
for i in "${thres[@]}"
do
    for j in $(seq 0.4 0.1 1.0)
    do
	#echo $i
	#var_name[$i]=$i
	cp setup_4.5_th.txt setup_"$i"_th_sigma_"$j".txt
	sed -i 's/POWERTHRESHOLD=-4.5/POWERTHRESHOLD=-'$i'/g' setup_"$i"_th_sigma_"$j".txt
	sed -i 's/SIGMA_THRES=0/SIGMA_THRES='$j'/g' setup_"$i"_th_sigma_"$j".txt
    done
done
