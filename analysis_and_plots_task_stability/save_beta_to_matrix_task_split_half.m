clear all
addpath(genpath('/mypath/utilities/cifti-matlab'));
addpath(genpath('/mypath/utilities/gifti/'));

SUB='0011' % example subject
permnum=462;
permlist=[1:permnum];

%half1 - rad in all beta dscalars and combine them in amatrix to facilitate further analysis
stabilitypath=['/mypath/oddball_task/stability/half1/ION' SUB '/'];
betamat=zeros(91282,permnum); %number of grayordinates
for i=1:size(permlist,2)
    filename=['sub-' SUB '_acq-3T2mm_contrast_oddball_zscored_perm' num2str(permlist(i)) '_h1.dscalar.nii'];
    betaval=cifti_read([stabilitypath filename]);
    betamat(:,i)=betaval.cdata;
end

save([stabilitypath 'sub-' SUB '_matrix_all_perm_h1.mat'], 'betamat')

%half2
stabilitypath2=['/mypath/oddball_task/stability/half2/ION' SUB '/'];
betamat2=zeros(91282,permnum);
for i=1:size(permlist,2)
    filename=['sub-' SUB '_acq-3T2mm_contrast_oddball_zscored_perm' num2str(permlist(i)) '_h2.dscalar.nii'];
    betaval=cifti_read([stabilitypath2 filename]);
    betamat2(:,i)=betaval.cdata;
end

save([stabilitypath2 'sub-' SUB '_matrix_all_perm_h2.mat'], 'betamat2')
