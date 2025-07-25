clear all
close all

addpath(genpath('/mypath/utilities/cifti-matlab'));
addpath(genpath('/mypath/utilities/gifti/'));

runnumlist=[24; 23; 13; 11; 20; 12; 13; 15; 11; 13]; %number of runs for each subject
pbsublist=['PB0007'; 'PB0009';'PB0010';'PB0011';'PB0015';'PB0016';'PB0017';'PB0018';'PB0019';'PB0020'];
sublist=['0001'; '0002';'0003';'0004';'0005';'0006';'0007';'0008';'0010';'0011'];
for n=1:size(sublist,1)
    SUB=sublist(n,:);
    runnum=runnumlist(n,1);
    %special case for PB0019 as one run is missing
    if SUB=='0010'
        runvec=['01'; '02'; '03'; '04'; '05'; '06'; '08'; '09'; '10'; '11'; '12'];
    else
        runvec=['01'; '02'; '03'; '04'; '05'; '06';  '07'; '08'; '09'; '10'; '11'; '12'; '13'; '14'; '15'; '16'; '17'; '18'; '19'; '20'; '21'; '22'; '23'; '24'];
    end

%% load vector for dice overlap and identify those with low dice
    dicevec=load(['/mypath/oddball_task/stability/sub-' SUB '_dice_vector_top_25pct.csv']);
    dicevec_all(:,n)=dicevec;
    [max_dice, maxdiceidx(n,1)]=max(dicevec);
    [min_dice, mindiceidx(n,1)]=min(dicevec);
    % load rund combinations
    codepath=['/mypath/task_analysis/task_stability/permuted_runs/'];
    load([codepath 'sub-' SUB '_run_combinations_perm.mat']); %variable is called B
    B=B(1:size(dicevec,1),:);
    allB{n}=B;


%% load motion files for each run
    derivspath=['/mypath/XCP-D_derivatives_task/ION' SUB '_MENORDIC_combined_task/sub-' SUB '/ses-MENORDIC/func/'];
    
    % read in motion files and output mean FD for each run
    for k=1:runnum
        try
            motionfile=readtable([derivspath 'sub-' SUB '_ses-MENORDIC_task-oddball_acq-3T2mm_run-' runvec(k,:) '_motion.tsv'], 'FileType','text', 'Delimiter', '\t');
        catch
            motionfile=readtable([derivspath 'sub-' SUB '_ses-MENORDIC_task-oddball_run-' runvec(k,:) '_motion.tsv'], 'FileType','text', 'Delimiter', '\t');
        end
        motionfile=table2array(motionfile);
        FDavg(k,1)=mean(motionfile(:,end),1);
    end


%create FD matrix for all possible run combos
    FDrunmatrix_all=zeros(size(B));
    for j=1:size(B,1)
        FDrunmatrix_all(j,:)=FDavg(B(j,:));
    end
    
    % calculate mean FD per half
    runnum_h=floor(runnum/2);
    FDh1_all=mean(FDrunmatrix_all(:,1:runnum_h),2);
    FDh2_all=mean(FDrunmatrix_all(:,runnum_h+1:runnum_h*2),2);
    FDdiff(:,n)=abs(FDh2_all-FDh1_all);
    FDh1_sub(:,n)=FDh1_all;
    FDh2_sub(:,n)=FDh2_all;
    FD_all_runs{n}=FDrunmatrix_all(1,:);
    FDavg_all(1,n)=mean(FDrunmatrix_all(1,:));

end

%define plot colors
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

% plot dice values in relation to difference in FD between halves
f3=figure
set(f3,'Color','w')
for j=1:10
    subplot(2,5,j)
    plot(FDdiff(:,j), dicevec_all(:,j), '.', 'Color', newDefaultColors(j,:))
    xlabel('FD difference')
    ylabel('dice coefficient')
    ylim([0 0.6])
    fontsize(12,"points")
    title(pbsublist(j,:))
    hold on
    p = polyfit(FDdiff(:,j), dicevec_all(:,j), 1);
    yfit = polyval(p, FDdiff(:,j));
    plot(FDdiff(:,j), yfit, 'k:', 'LineWidth', 2);
end
box off

% test fit
for j=1:10
    [r,p]=corrcoef(FDdiff(:,j), dicevec_all(:,j));
    rvec(j,1)=r(1,2);
    pvec(j,1)=p(1,2);
end

%% pick particulalry high/low FD runs

for j=1:10
    min_FDh1(:,j)=find(FDh1_sub(:,j)==min(FDh1_sub(:,j)));
    max_FDh1(:,j)=find(FDh1_sub(:,j)==max(FDh1_sub(:,j)));
end


%% calculate dice overlap between all of the runs and the lalve with highest and lowest motion
sublist=['0001'; '0002';'0003';'0004';'0005';'0006';'0007';'0008';'0010';'0011'];
for l=1:10
    SUB=sublist(l,:)
    subnum=l;
    folder=['/mypath/oddball_task/'];
    output='/mypath/oddball_task/stability/';
    
    %load .mat file that contains betas
    load([output 'half1/ION' SUB '/sub-' SUB '_matrix_all_perm_h1.mat']);
    
    % read in betas for all data
    b=cifti_read([folder 'sub-' SUB '_acq-3T2mm_contrast_oddball_glover_zscored.dscalar.nii']);
    
    thresh_all=prctile((b.cdata(b.cdata>0)),75); %25 percent highest betas
    datahigh_all=find(b.cdata>=thresh_all); %find data above 25%
    data_matrix_all=nan(size(b.cdata));% create Nan matrix
    data_matrix_all(datahigh_all)=b.cdata(datahigh_all);
    binary_datamat_all=zeros(size(data_matrix_all));
    binary_datamat_all(~isnan(data_matrix_all),1)=1;
    
	%same for betas from high motion half
    high_motion_avg=betamat(:,max_FDh1(:,subnum));
    
    thresh_high=prctile((high_motion_avg(high_motion_avg>0)),75); %25 percent highest betas
    datahigh_hFD=find(high_motion_avg>=thresh_high); %find data above 25%
    data_matrix_hFD=nan(size(high_motion_avg));% create Nan matrix
    data_matrix_hFD(datahigh_hFD)=high_motion_avg(datahigh_hFD);
    binary_datamat_hFD=zeros(size(data_matrix_hFD));
    binary_datamat_hFD(~isnan(data_matrix_hFD),1)=1;
    
	%same for betas from low motion half
    low_motion_avg=betamat(:,min_FDh1(:,subnum));

    thresh_low=prctile((low_motion_avg(low_motion_avg>0)),75); %25 percent highest betas
    datahigh_lFD=find(low_motion_avg>=thresh_low); %find data above 25%
    data_matrix_lFD=nan(size(low_motion_avg));% create Nan matrix
    data_matrix_lFD(datahigh_lFD)=low_motion_avg(datahigh_lFD);
    binary_datamat_lFD=zeros(size(data_matrix_lFD));
    binary_datamat_lFD(~isnan(data_matrix_lFD),1)=1;

    overlap_dice_low(l,1)=dice(binary_datamat_all, binary_datamat_lFD);
    overlap_dice_high(l,1)=dice(binary_datamat_all, binary_datamat_hFD);
    overlap_dice_mixed(l,1)=dice(binary_datamat_lFD, binary_datamat_hFD);

end

% plot dice for high and low motion halves
f3=figure
set(f3,'Color','w')
plot([1:10],overlap_dice_high, 'r.', 'MarkerSize', 40)
hold on
plot([1:10],overlap_dice_low, 'b.', 'MarkerSize', 40)
plot([1:10],overlap_dice_mixed, 'k.', 'MarkerSize', 40)
xticklabels(pbsublist)
xlim([0.5 10.5])
fontsize(12,"points")
ylabel('dice coefficient')
legend('dice full data and high FD runs', 'dice full data and low FD runs', 'dice high and low FD runs', 'Location', 'southeast')
box off

% test id difference is significant
[h,p,CI,STATS]=ttest(overlap_dice_high, overlap_dice_low)


%% boxplot for FD by run
FD_boxplot=cell2mat(FD_all_runs);
group_boxplot=[zeros(1,size(FD_all_runs{1},2)), ones(1,size(FD_all_runs{2},2)), 2*ones(1,size(FD_all_runs{3},2)), 3*ones(1,size(FD_all_runs{4},2)), 4*ones(1,size(FD_all_runs{5},2)), ...
    5*ones(1,size(FD_all_runs{6},2)), 6*ones(1,size(FD_all_runs{7},2)), 7*ones(1,size(FD_all_runs{8},2)), 8*ones(1,size(FD_all_runs{9},2)), 9*ones(1,size(FD_all_runs{10},2))];

f2=figure
set(0, 'CurrentFigure', f2)
boxplot(FD_boxplot, group_boxplot, 'Color', newDefaultColors)   
colormap(newDefaultColors)
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),get(h(j),'Color'),'FaceAlpha',.3);
end
ylabel('mean FD by run')
xticks([1:10])
xticklabels(pbsublist)
set(f2,'Color','w')
xlim padded
fontsize(16,"points")
box off
hold on
xCenter = 1:numel(FD_all_runs); 
spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(FD_all_runs)
    plot(rand(size(FD_all_runs{i}))*spread -(spread/2) + xCenter(i), FD_all_runs{i}, 'o','linewidth', 1, 'Color', newDefaultColors(i,:))
end

%% outliers by run instead of motion by run
for n=1:size(sublist,1)
    SUB=sublist(n,:);
    runnum=runnumlist(n,1);
    if SUB=='0010'
        runvec=['01'; '02'; '03'; '04'; '05'; '06'; '08'; '09'; '10'; '11'; '12'];
    else
        runvec=['01'; '02'; '03'; '04'; '05'; '06';  '07'; '08'; '09'; '10'; '11'; '12'; '13'; '14'; '15'; '16'; '17'; '18'; '19'; '20'; '21'; '22'; '23'; '24'];
    end
    derivspath=['/mypath/XCP-D_derivatives_task/ION' SUB '_MENORDIC_combined_task/sub-' SUB '/ses-MENORDIC/func/'];
    
    for k=1:runnum
        try
            outliersfile=readtable([derivspath 'sub-' SUB '_ses-MENORDIC_task-oddball_acq-3T2mm_run-' runvec(k,:) '_outliers.tsv'], 'FileType','text', 'Delimiter', '\t');
        catch
            outliersfile=readtable([derivspath 'sub-' SUB '_ses-MENORDIC_task-oddball_run-' runvec(k,:) '_outliers.tsv'], 'FileType','text', 'Delimiter', '\t');
        end
        outliersfile=table2array(outliersfile);
        outliers_total(1,k)=sum(outliersfile);
    end
outliers_all{n}=outliers_total;
 clear outliers_total
end

%% same plot for outliers by run
ouliers_boxplot=cell2mat(outliers_all);

f4=figure
set(0, 'CurrentFigure', f4)
boxplot(ouliers_boxplot, group_boxplot, 'Color', newDefaultColors)   
colormap(newDefaultColors)
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),get(h(j),'Color'),'FaceAlpha',.3);
end
ylabel('scrubbed frames by run')
xticks([1:10])
xticklabels(pbsublist)
set(f4,'Color','w')
xlim padded
fontsize(16,"points")
box off
hold on
xCenter = 1:numel(outliers_all); 
spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(outliers_all)
    plot(rand(size(outliers_all{i}))*spread -(spread/2) + xCenter(i), outliers_all{i}, 'o','linewidth', 1, 'Color', newDefaultColors(i,:))
end

%% figure for number of runs
f2=figure
set(0, 'CurrentFigure', f2)
scatter([1:10], runnumlist, 100, newDefaultColors, 'filled')
colormap(newDefaultColors)
ylabel('number of runs')
xticklabels(pbsublist)
ylim([0, 25])
set(f2,'Color','w')
xlim padded
box off


%% figure for overall average FD
f2=figure
set(0, 'CurrentFigure', f2)
scatter([1:10], FDavg_all, 100, newDefaultColors, 'filled')
colormap(newDefaultColors)
ylabel('average FD')
xticklabels(pbsublist)
%ylim([0, 25])
set(f2,'Color','w')
xlim padded
box off
