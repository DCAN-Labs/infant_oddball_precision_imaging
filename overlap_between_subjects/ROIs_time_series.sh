#!/bin/bash


#   wb_command -surface-geodesic-rois
#      <surface> - the surface to draw on
#      <limit> - geodesic distance limit from vertex, in mm
#      <vertex-list-file> - a text file containing the vertices to draw ROIs
#         around
#      <metric-out> - output - the output metric
SUB=${1}
vertex=${2}
side=${3}
ACQ=3T2mm
SES=MENORDIC
TASK=oddball

module load workbench/1.5.0

#save vertex as txt file as function expects this input
temp_dir=/tmp/vertex_list
mkdir -p ${temp_dir}
echo ${vertex} > ${temp_dir}/list_vertex.txt

derivs_dir=/mypath/XCP-D_derivatives_task/ION${SUB}_${SES}_combined_task/sub-${SUB}/ses-${SES};
results_dir=/mypath/oddball_task

#select appropriate surface
if [ "${side}" == "right" ]; then
selected_surface=${derivs_dir}/anat/sub-${SUB}_ses-${SES}_hemi-R_space-fsLR_den-32k_desc-hcp_midthickness.surf.gii
elif [ "${side}" == "left" ]; then
selected_surface=${derivs_dir}/anat/sub-${SUB}_ses-${SES}_hemi-L_space-fsLR_den-32k_desc-hcp_midthickness.surf.gii

fi

#select size
limit=10
list_vertex=${temp_dir}/list_vertex.txt
metric_out=${results_dir}/sub-${SUB}_vertex-${vertex}_ROI_${side}.func.gii

#create ROI based on surface mesh
wb_command -surface-geodesic-rois $selected_surface $limit $list_vertex $metric_out

out_ts=${results_dir}/sub-${SUB}_vertex-${vertex}_${side}_average_time_series.txt
dtseries=${derivs_dir}/func/sub-${SUB}_ses-${SES}_task-${TASK}_acq-${ACQ}_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_SMOOTHED_2.25.dtseries.nii

# average tiem series over this ROI
wb_command -cifti-roi-average $dtseries $out_ts -${side}-roi $metric_out

#   wb_command -cifti-create-label
#      <cifti-out> - output - the output cifti file
#      [-left-label] - label file for left surface
#         <label> - the label file
#         [-roi-left] - roi of vertices to use from left surface
#            <roi-metric> - the ROI as a metric file



#$wbc -cifti-create-label $cifti_out -left-label 1 -roi-left $metric_out
