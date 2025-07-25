clear all

addpath(genpath('/mypath/utilities/cifti-matlab'));
%
sublist=['0001'; '0002';'0003';'0004';'0005';'0006';'0007';'0008';'0010';'0011'];
for x=1:size(sublist,1)
    SUB=sublist(x,:)
    output='/mypath/oddball_task';
    %%
    %betafile=['sub-' SUB '_acq-3T2mm_contrast_oddball_glover_zscored.dscalar.nii']
    betafile=['habituation/sub-' SUB '_acq-3T2mm_contrast_oddball_zscored_habi_late.dscalar.nii']
    
    load([output '/noise_distribution/sub-' SUB '_acq-3T2mmnoise_distribution_99.5pctl_habi_late.mat']);
    load([output '/noise_distribution/sub-' SUB '_acq-3T2mmnoise_distribution_0.5pctl_habi_late.mat']);
    
    pctlhigh=pctlvhigh; %for p<0.01
    pctllow=pctlvlow;
    
    
    b=cifti_read([output '/' betafile]);
    data=b.cdata;
    
    %seperate positive and negative values to make sure they are seperately
    %considered and not mixed up
    dataposidx=find(data>0);
    datapos=zeros(size(data));
    datapos(dataposidx)=data(dataposidx);
    
    datanegidx=find(data<0);
    dataneg=zeros(size(data));
    dataneg(datanegidx)=data(datanegidx);
    %quick check to make sure that all percentile values are indeed
    %positive/negative
    
    dataposidxtest=find(pctlhigh>0);
    datanegidxtest=find(pctllow<0);
    
    % find values that are greater than threshold
    datahigh=find((datapos-pctlhigh)>0);
    datalow=find((dataneg-pctllow)<0);
    
    %fill in values into data matrix
    data_matrix=zeros(size(data));
    data_matrix(datahigh)=data(datahigh);
    data_matrix(datalow)=data(datalow);
    
    % save data as dscalar
    
    b.cdata=data_matrix;
    cifti_write(b, [output '/sub-' SUB '_acq-3T2mm_thresholded_p0.01_from_noise_distr.dscalar.nii']);
    

end
