clear all

addpath(genpath('/mypath/utilities/cifti-matlab'));
%
sublist=['0001'; '0002';'0003';'0004';'0005';'0006';'0007';'0008';'0010';'0011'];
for x=1:size(sublist,1)
    SUB=sublist(x,:)
    output='/mypath/oddball_task';
    %%
    betafile=['sub-' SUB '_acq-3T2mm_contrast_oddball_glover_zscored.dscalar.nii']
       
    b=cifti_read([output '/' betafile]);
    data=b.cdata;
    
    %seperate positive and negative values to make sure they are seperately
    %considered and not mixed up
    thresh_high=prctile((data(data>0)),75); %25 percent highest betas
    thresh_low=prctile(data(data<0),25); % 25 percent lowest betas
    datahigh=find(data>=thresh_high);
    datalow=find(data<=thresh_low);

    datapos=zeros(size(data));
    datapos(datahigh)=data(datahigh);
    
    dataneg=zeros(size(data));
    dataneg(datalow)=data(datalow);
    
    %fill in values into data matrix
    data_matrix=datapos+dataneg;

    
    % save data as dscalar
    
    b.cdata=data_matrix;
    cifti_write(b, [output '/sub-' SUB '_acq-3T2mm_thresholded_at25_percent.dscalar.nii']);
    
    threshlistH(x,1)=thresh_high;
    threshlistL(x,1)=thresh_low;


end

thresholdstab=table(str2num(sublist), threshlistH, threshlistL);
thresholdstab=renamevars(thresholdstab, 'Var1', 'subject');

writetable(thresholdstab, [output '/threshold_values_25percentile_all_sub.csv'])
