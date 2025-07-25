clear all

addpath(genpath('/mypath/cifti-matlab'));
%
folder=['/mypath/oddball_task'];
output='/mypath/oddball_task/';
%filelist=dir([folder '/*glover_pval.dscalar*']); %same code used for overlap between significant values and overlap between 75% strongest beta values
filelist=dir([folder '/*25_per*.dscalar*']);
%%
for n=1:size(filelist)
    b=cifti_read([folder '/' filelist(n,1).name]);
    data=b.cdata;
    datasign=find(data>0); %input data is already thresholded. So everything below the threshold is 0 (or negative for negative significant values)
    data_matrix=nan(size(data));
    data_matrix(datasign)=data(datasign);
    data_all(:,n)=data_matrix;
end
% overlap will be counted from 0-10
for p=1:size(data_all,1)
    row=data_all(p,:);
    rowidx=find(~isnan(row));
    rowcount=size(rowidx,2);
    overlap(p,1)=rowcount;
end

b.cdata=overlap;
cifti_write(b, [output 'sub-overlap_acq-3T2mm_25_percentile_pos.dscalar.nii']);

%% add on analysis: how does the percentage of overlap change when using all, early or late oddballs
a=cifti_read([output 'sub-overlap_acq-3T2mm_pval_0.01_pos_noise_dist.dscalar.nii']);
e=cifti_read([output 'habituation/sub-overlap_acq-3T2mm_pval_0.01_pos_noise_dist_habi_early.dscalar.nii']);
l=cifti_read([output 'habituation/sub-overlap_acq-3T2mm_pval_0.01_pos_noise_dist_habi_late.dscalar.nii']);

ao=find(a.cdata>=5);
pcta=size(ao,1)/size(a.cdata,1)*100;

eo=find(e.cdata>=5);
pcte=size(eo,1)/size(e.cdata,1)*100;

lo=find(l.cdata>=5);
pctl=size(lo,1)/size(l.cdata,1)*100;
