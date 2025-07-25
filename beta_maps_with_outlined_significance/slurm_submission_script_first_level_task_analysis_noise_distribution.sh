#!/bin/bash -l

#SBATCH -J nilearn_task
#SBATCH --ntasks=6
#SBATCH --tmp=80gb
#SBATCH --mem=80gb
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH -p msismall
#SBATCH -o output_logs/pycharm_%A_%a.out
#SBATCH -e output_logs/pycharm_%A_%a.err

SUB=${1}
ACQ=${2}
permnum=1000

work_dir=/tmp/sub-${SUB}/acq-${ACQ}/noise

if [ ! -d ${work_dir} ]; then
	mkdir -p ${work_dir}

fi

source /mypath/envs/miniconda3/load_miniconda3.sh; 
conda activate cabinet

MINVAR=$(seq 1 1 ${permnum}); 

for PERM in ${MINVAR}; 
do
python3 /mypath/task_analysis/surface-based_first-level_analysis_ION_noise_distribution.py ${SUB} ${ACQ} ${PERM}
done
#save outputs to matrix
module load matlab
matlab -nodisplay -nosplash -r "addpath('/mypath/task_analysis'); noise_distribution('${work_dir}','${SUB}', '${ACQ}'); exit;";
