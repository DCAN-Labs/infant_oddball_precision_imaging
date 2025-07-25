function select_runs_save_timeseries_split_half(SUB, SES, TASK, ACQ,  permnum, outfolder, infile, outlierfile)

addpath(genpath('/mypath/utilities/cifti-matlab'));

% example inputs: 
%SUB='0011'
%SES='MENORDIC'
%TASK='oddball';
%ACQ='3T2mm'
%permnum=3;
%infolder=['/mypath/XCP-D_derivatives/ION0011_MENORDIC_combined_task/sub-' SUB '/ses-' SES '/func'];
%infile=[infolder '/sub-' SUB '_ses-' SES '_task-' TASK '_acq-' ACQ '_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_SMOOTHED_2.25.dtseries.nii']
%outlierfile=[infolder '/sub-' SUB '_ses-' SES '_task-' TASK '_acq-' ACQ '_outliers.tsv']

%% step 1, load random run order
load(['permuted_runs/sub-' SUB '_run_combinations_perm.mat']) %variable is called 'B'
run_order=B(permnum,:);
sizehalf=floor(size(run_order,2)/2); %take half and round to floor
%runs half 1
run_order1=sort(run_order(1,1:sizehalf));
%runs half 2
run_order2=sort(run_order(1,sizehalf+1:sizehalf*2));
%% step 2, load concatenated dtseries

%load concatenated time series
concatenated_timeseries=cifti_read([infile]);
% list of frames per run - important in case of incomplete runs to make sure run boarders are detrmined accurately
load(['sub-' SUB '_framelist.mat']);

%create array with starting and stopping point of frame list
framelist_start=1;
for i=2:size(framelist,1)
    framelist_start(i,1)=framelist(i-1,1)+framelist_start(i-1,1);
end

for i=1:size(framelist,1)
    framelist_stop(i,1)=framelist_start(i,1)+framelist(i,1)-1;
end

%split data by run according to the frames per run
for i=1:size(framelist,1)
    data_by_run{i}=concatenated_timeseries.cdata(:,framelist_start(i):framelist_stop(i));
end

%% step 3, load motion file and create masks
%load motion file
outliers=readtable([outlierfile], "FileType","text",'Delimiter', '\t');

retained_frames=abs(table2array(outliers)-1);
%% step 3.1, sort motion file according to new run order

%split data by run according to the frames per run
for i=1:size(framelist,1)
    fd_by_run{i}=retained_frames(framelist_start(i):framelist_stop(i));
end
% for run order 1
for i=1:size(framelist,1)
    fd_by_run_empty{i}=zeros(size(fd_by_run{i}));
end
%fill in numbers for relevant runs
for k=1:size(run_order1,2)
    fd_by_run_empty{run_order1(1,k)}=fd_by_run{run_order1(1,k)};
end
%concatenate data again
k=1;
rearranged_fd=[fd_by_run_empty{1}];
while k<size(fd_by_run_empty,2)
rearranged_fd=[rearranged_fd; fd_by_run_empty{k+1}];
k=k+1;
end
% and for run order 2
for i=1:size(framelist,1)
    fd_by_run_empty2{i}=zeros(size(fd_by_run{i}));
end
%fill in numbers for relevant runs
for k=1:size(run_order2,2)
    fd_by_run_empty2{run_order2(1,k)}=fd_by_run{run_order2(1,k)};
end
%concatenate data again
k=1;
rearranged_fd2=[fd_by_run_empty2{1}];
while k<size(fd_by_run_empty2,2)
rearranged_fd2=[rearranged_fd2; fd_by_run_empty2{k+1}];
k=k+1;
end
%% step 3.3, create mask and write them to tmp space

writematrix(rearranged_fd, [outfolder '/sub-' SUB '_ses-' SES '_acq-' ACQ '_mask_runs.txt']);
writematrix(rearranged_fd2, [outfolder '/sub-' SUB '_ses-' SES '_acq-' ACQ '_mask_runs2.txt']);

%% make sure motion censoring is also applied to time series

for i=1:size(framelist,1)
    data_by_run{i}=data_by_run{i}(:,logical(fd_by_run{i}));
end
%for run order 1
k=1;
rearranged_data=[data_by_run{1,run_order1(k)}];
while k<size(run_order1,2)
rearranged_data=[rearranged_data, data_by_run{1,run_order1(k+1)}];
k=k+1;
end

new_timeseries=concatenated_timeseries;
new_timeseries.cdata=rearranged_data;
new_timeseries.diminfo{1,2}.length=size(rearranged_data,2);


cifti_write(new_timeseries, [outfolder '/sub-' SUB '_ses-' SES '_task-' TASK '_acq-' ACQ '_bold_shuffled_timeseries_scrubbed.dtseries.nii']);
%for run order 2
k=1;
rearranged_data2=[data_by_run{1,run_order2(k)}];
while k<size(run_order2,2)
rearranged_data2=[rearranged_data2, data_by_run{1,run_order2(k+1)}];
k=k+1;
end

new_timeseries2=concatenated_timeseries;
new_timeseries2.cdata=rearranged_data2;
new_timeseries2.diminfo{1,2}.length=size(rearranged_data2,2);

cifti_write(new_timeseries2, [outfolder '/sub-' SUB '_ses-' SES '_task-' TASK '_acq-' ACQ '_bold_shuffled_timeseries_scrubbed2.dtseries.nii']);


end

