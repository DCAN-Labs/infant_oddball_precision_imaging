
import nibabel as nb
import nilearn
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from nilearn.plotting import plot_design_matrix

import os
#from docopt import docopt
import argparse

def task_code_ION(sub, acq, perm):
    # define TR
    if "7T" in acq:
        t_r = 1.768
    else: 
        t_r = 1.761
    
    ses='MENORDIC'
    task='oddball'
    slice_time_ref = 0  # ??
    # read in cifti data and transpose it to expected shape
    defined_path = f"/mypath/XCP-D_derivatives_task/ION{sub}_{ses}_combined_task/sub-{sub}/ses-{ses}/func/"
    fname = os.path.join(defined_path,
                         f"sub-{sub}_ses-{ses}_task-{task}_acq-{acq}_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_SMOOTHED_2.25.dtseries.nii")
    cifti = nb.load(fname)
    cifti_data = cifti.get_fdata(dtype='f4')
    texture = cifti_data.transpose()

    # read in event file
    event_path = '/mypath/task_analysis/'
    if "7T" in acq:
        ename = os.path.join(event_path, f"ION{sub}_7T_event_file.csv")
    else: 
        ename = os.path.join(event_path, f"ION{sub}_event_file.csv")
 
    events = pd.read_table(ename, delimiter=',')

    # read in motion trace
    mname = os.path.join(defined_path, f'sub-{sub}_ses-{ses}_task-{task}_acq-{acq}_motion.tsv')
    motion = pd.read_table(mname, delimiter='\t')
    motion_reg = motion[['framewise_displacement', 'rot_x', 'rot_y', 'rot_z', 'trans_x', 'trans_y', 'trans_z']]


    # define number of frames and their timing
    n_scans = texture.shape[1]
    frame_times = t_r * (np.arange(n_scans) + 0)  # relates to slice time ref? # initially 0.5

    # Create the design matrix.
    # We specify an HRF model containing the Glover model and its time derivative The drift model is implicitly a cosine basis with a period cutoff at 128s.

    from nilearn.glm.first_level import make_first_level_design_matrix

    initial_design_matrix = make_first_level_design_matrix(frame_times,
                                                           events=events,
                                                           hrf_model='glover + derivative',
                                                           add_regs=motion_reg,
                                                           drift_model=None
                                                           )
    
    # read in mask to remove data that is not included in the permutation and outliers
    defined_path2 = f"/tmp/sub-{sub}/acq-{acq}/perm{perm}/"
    oname = os.path.join(defined_path2, f'sub-{sub}_ses-{ses}_acq-{acq}_mask_runs2.txt')
    outliers = pd.read_table(oname, delimiter='\t',header=None)
    outlier_logical = outliers.astype(bool)

    design_matrix = initial_design_matrix[outlier_logical.values]
   
    # load cleaned cifti data
    fname2 = os.path.join(defined_path2,
                          f'sub-{sub}_ses-{ses}_task-{task}_acq-{acq}_bold_shuffled_timeseries_scrubbed2.dtseries.nii')
    cifti = nb.load(fname2)
    cifti_data = cifti.get_fdata(dtype='f4')
    texture = cifti_data.transpose()

    # Setup and fit GLM.
    # Note that the output consists in 2 variables: labels and fit.
    # labels tags voxels according to noise autocorrelation. estimates contains the parameter estimates.
    # We keep them for later contrast computation.

    from nilearn.glm.first_level import run_glm

    labels, estimates = run_glm(texture.T, design_matrix.values)

    # Estimate contrasts
    # Specify the contrasts.

    # For practical purpose, we first generate an identity matrix whose size is the number of columns of the design matrix.

    contrast_matrix = np.eye(design_matrix.shape[1])

    # At first, we create basic contrasts.

    basic_contrasts = {column: contrast_matrix[i]
                       for i, column in enumerate(design_matrix.columns)}

    # Next, we add some intermediate contrasts and one contrast adding all conditions with some auditory parts.

    contrasts = {
        'oddball': (
            basic_contrasts['Noise']+
            basic_contrasts['Noise_derivative'])
    }
    # Let's estimate the contrasts by iterating over them.

    from nilearn.glm.contrasts import compute_contrast

    for index, (contrast_id, contrast_val) in enumerate(contrasts.items()):
        print(f"  Contrast {index + 1:1} out of {len(contrasts)}: "
              f"{contrast_id}")
        # compute contrast-related statistics
        contrast = compute_contrast(labels, estimates, contrast_val, 't')
        # we present the Z-transform of the t map
        z_score = contrast.z_score()
        # read in dscalar example to write contrast z_score values on
        cifti_dscalar = nb.load('/mypath/test.dscalar.nii')
        # reshape data to fit expected format
        inp_dat = np.reshape(z_score, (1, -1))
        # create object with header
        cifti_data_dscalar3 = nb.cifti2.cifti2.Cifti2Image(dataobj=inp_dat, header=cifti_dscalar.header)
        # save cifti
        nb.cifti2.save(cifti_data_dscalar3, f"/mypath/oddball_task/stability/half2/ION{sub}/sub-{sub}_acq-{acq}_contrast_{contrast_id}_zscored_perm{perm}_h2.dscalar.nii")
 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="task analysis script ION stability")
    parser.add_argument("sub", help="subject ID")
    parser.add_argument("acq", help="acquisition type (3T or 7T and resolution)")
    parser.add_argument("perm", help="permutation number for stability analysis. Determines the runs included in a given calculation")

    args = parser.parse_args()
    task_code_ION(args.sub, args.acq, args.perm)

#if __name__ == '__main__':
 #   args = docopt(__doc__)
  #  test_code_nelson(args['<sub>'])
