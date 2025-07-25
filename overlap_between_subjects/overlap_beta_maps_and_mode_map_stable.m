clear all
%ad relavant paths
addpath(genpath('/mypath/utilities/cifti-matlab'));
wb_command='/mypath/utilities/workbench/1.4.2/workbench/bin_rh_linux64/wb_command';


output='/mypath/oddball_task/stability';
filelist=dir([output '/*25p*.dscalar.nii']); %list contains file for each subject with percentage of overlap of overlapping split halves
%%
% read in all betas and determine threshold

for n=1:size(filelist)
    b=cifti_read([output '/' filelist(n,1).name]);
    data=b.cdata;
    thresh_high=50; %50% of cases show overlap
    datahigh=find(data>=thresh_high);
    data_matrix=zeros(size(data));
    data_matrix(datahigh)=data(datahigh);
    data_all(:,n)=data_matrix;
end

% read in mode map for each participant
sublist=['0001'; '0002';'0003';'0004';'0005';'0006';'0007';'0008';'0010';'0011'];
for j=1:size(sublist,1)
    SUB=sublist(j,:);
    map_path=['/mypath/TM_results/infant_template/no_task/surfonly/'];
    map_file=['sub-' SUB '_ses-MENORDIC_task-oddball_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_template_matched_Zscored_recolored.dscalar.nii'];
    try
        mode_map=cifti_read([map_path map_file]);
    catch
        map_file=['sub-' SUB '_ses-MENORDIC_task-oddball_acq-3T2mm_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_template_matched_Zscored_recolored.dscalar.nii'];
        mode_map=cifti_read([map_path map_file]);
    end
    mode_map_all(:,j)=mode_map.cdata;
      %count network vertices
    for k=1:18
        idxk=find(mode_map.cdata==k);
        netwsizesub(k,1)=size(idxk,1);
    end
    netwsize_all(:,j)=netwsizesub;
end


% get network name list (network_names)
% read in eta file for the purpose of having network labels in matlab
eta_files=dir([map_path '*_template_matched_Zscored.mat']);
load([map_path eta_files(1,1).name])
namelist=network_names;
namelist([4,6])=[] %remove empty names

%cut data all to size of surface vertices as infant template only has hese
data_all=data_all(1:size(mode_map_all,1), :);

% find number of network vertices with overlap and count
% number of vertives in each network

for j=1:size(data_all,2)
    idxpos=find(data_all(:,j)>0);
    posnetlist=mode_map_all(idxpos, j);
    for k=1:18
        idxp=find(posnetlist==k);
        netwproppos(k,j)=size(idxp,1);
    end
    netwpctpos(:,j)=netwproppos(:,j)*100./netwsize_all(:,j);
end

netwpctpos([4,6, 17, 18],:)=[]; %remove empties

% create color map for histograms
colormap = [246, 34, 47
    48, 40, 172
    252, 255, 64
    55, 199, 54
    55, 177, 177
    39, 41, 38
    110, 38, 169 
    97, 233, 229
    251, 149, 65
    182, 87, 251
    44, 78, 114
    105, 247, 103
    53, 39, 252
    252, 253, 212]; 
colormap=colormap./255;
%'DMN'	'Vis'	'FP'	'DAN'	'VAN'	'Sal'	'CO'	'SMd'	'SMl'	'Aud'	'Tpole'	'MTL'	'PMN'	'PON'


pbsublist=['PB0007'; 'PB0009';'PB0010';'PB0011';'PB0015';'PB0016';'PB0017';'PB0018';'PB0019';'PB0020'];


%%
addpath(genpath('/mypath/spider_plot/'))
newDefaultColors = ([43 66 49
    34 136 51
    147 157 92
    220 155 65
    202 91 72
    225 151 144
    170 51 119
    56 37 133
    86 180 233
    187 187 187])./255;

%spiderplot with all subjects in one plot
P=netwpctpos';
figure
spider_plot(P,...
    'AxesLimits', [0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
    50 50 50 50 50 50 50 50 50 50 50 50 50 50],...
    'AxesLabels', namelist, ...
    'Color', newDefaultColors,...
    'AxesInterval', 5,...
    'AxesDisplay', 'one',...
    'AxesPrecision', 0,...
    'AxesLabelsRotate', 'on',...
    'AxesLabelsOffset', 0.1,...
    'AxesRadial', 'off');
legend(pbsublist, 'Location', 'eastoutside')

%spiderplot with all subjects individually
figure
for i=1:10
subplot(2,5,i)
P=netwpctpos(:,i)';

spider_plot(P,...
    'AxesLimits', [0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
    50 50 50 50 50 50 50 50 50 50 50 50 50 50],...
    'AxesLabels', namelist, ...
    'Color', newDefaultColors(i,:),...
    'AxesInterval', 5,...
    'AxesDisplay', 'one',...
    'AxesPrecision', 0,...
    'AxesLabelsRotate', 'on',...
    'AxesLabelsOffset', 0.1,...
    'AxesRadial', 'off');
title(sprintf(pbsublist(i,:)),...
    'FontSize', 12);
ax = gca;
ax.TitleHorizontalAlignment = 'left';

end

writematrix(netwpctpos, [output '/values_networks_and_significant_stable_betas.csv'])
