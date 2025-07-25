#!/bin/bash -l

#SBATCH -J nilearn_task
#SBATCH --ntasks=6
#SBATCH --tmp=40gb
#SBATCH --mem=40gb
#SBATCH -t 01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=NONE
#SBATCH -p msismall
#SBATCH -o output_logs/pycharm_%A_%a.out
#SBATCH -e output_logs/pycharm_%A_%a.err

SUB=${1}
ACQ=${2}

source /mypath/envs/miniconda3/load_miniconda3.sh; 
conda activate cabinet

python3 /mypath/task_analysis/surface-based_first-level_analysis_ION_combined_glover.py ${SUB} ${ACQ}
