#!/bin/bash -l
#this script is used to submit the task analysis script for each of teh permutations used in the task stability analysis

SUB=${1}
permnum=${2} #max number of permutations


MINVAR=$(seq 1 1 ${permnum}); 

for PERM in ${MINVAR}; 
do 
sbatch sbatch_first_level_task_analysis_script_stability_par_inside.sh ${SUB} ${PERM}
done
