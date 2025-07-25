function noise_distribution(work_dir, SUB, ACQ)
%example inputs:
%work_dir='/tmp/sub-0011/acq-3T2mm/noise'
%SUB='0011'
%ACQ='3T2mm'
addpath(genpath('/mypath/cifti-matlab'));
%%
output='/mypath/oddball_task/noise_distribution';
%%
filelist=dir([work_dir '/*zscore*dscalar*']);
%
for n=1:size(filelist)
    b=cifti_read([work_dir '/' filelist(n,1).name]);
    data=b.cdata;
    data_all(:,n)=data;
end

pctlhigh=prctile(data_all, 97.5,2);
pctllow=prctile(data_all, 2.5,2);

pctlvhigh=prctile(data_all, 99.5,2);
pctlvlow=prctile(data_all, 0.5,2);

save([output '/sub-' SUB '_acq-' ACQ 'noise_distribution.mat'], 'data_all')
save([output '/sub-' SUB '_acq-' ACQ 'noise_distribution_97.5pctl.mat'], 'pctlhigh')
save([output '/sub-' SUB '_acq-' ACQ 'noise_distribution_2.5pctl.mat'], 'pctllow')
save([output '/sub-' SUB '_acq-' ACQ 'noise_distribution_99.5pctl.mat'], 'pctlvhigh')
save([output '/sub-' SUB '_acq-' ACQ 'noise_distribution_0.5pctl.mat'], 'pctlvlow')

end

