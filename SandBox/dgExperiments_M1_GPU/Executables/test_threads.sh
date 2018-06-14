#!/bin/bash

arr=(1 2 4 8)

source ~/thornado/Workflow/SetEnvironment.sh mcarpe21

for thread in "${arr[@]}"
do

	for sec_thread in "${arr[@]}"
	do

		export OMP_NUM_THREADS=$thread

		export MKL_NUM_THREADS=$sec_thread
	
		./StreamingSineWave_mcarpe21 > data.txt

		python sort_data.py $thread $sec_thread

	done

done


