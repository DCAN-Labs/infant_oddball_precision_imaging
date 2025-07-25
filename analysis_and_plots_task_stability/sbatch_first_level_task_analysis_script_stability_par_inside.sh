#!/bin/bash -l

#SBATCH -J nilearn_task
#SBATCH --ntasks=6
#SBATCH --tmp=20gb
#SBATCH --mem=20gb
#SBATCH -t 00:07:00
#SBATCH --mail-type=NONE
#SBATCH -p msismall


SUB=${1}
SES=MENORDIC
TASK=oddball
ACQ='3T2mm'
#PERM=1 # permutation number
PERM=${2}


#define list of included runs
#create concatinated timeseries that only include these runs
#create mask to mask design matrix

in_dir=/mypath/XCP-D_derivatives_task/ION${SUB}_${SES}_combined_task/sub-${SUB}/ses-${SES}/func
dtseries_in=${in_dir}/sub-${SUB}_ses-${SES}_task-${TASK}_acq-${ACQ}_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_SMOOTHED_2.25.dtseries.nii
outlier_file=${in_dir}/sub-${SUB}_ses-${SES}_task-${TASK}_acq-${ACQ}_outliers.tsv



module load matlab

source /mypath/envs/miniconda3/load_miniconda3.sh; 
conda activate cabinet


work_dir=/tmp/sub-${SUB}/acq-${ACQ}/perm${PERM}

if [ ! -d ${work_dir} ]; then
	mkdir -p ${work_dir}

fi

matlab -nodisplay -nosplash -r "addpath('/mypath/task_analysis'); select_runs_save_timeseries_split_half('${SUB}', '${SES}', '${TASK}', '${ACQ}', ${PERM}, '${work_dir}', '${dtseries_in}', '${outlier_file}'); exit;";

#python task analysis script to create beta maps for half 1
python3 /home/yaco0006/shared/projects/Baby7T/code/task_analysis/task_stability/surface-based_first-level_analysis_ION_stability.py ${SUB} ${ACQ} ${PERM};
#python task analysis script to create beta maps for half 2
python3 /home/yaco0006/shared/projects/Baby7T/code/task_analysis/task_stability/surface-based_first-level_analysis_ION_stability2.py ${SUB} ${ACQ} ${PERM};



