clear all

addpath(genpath('/mypath/utilities/cifti-matlab'));

folder=['/mypath/analysis/oddball_task/'];
output='/mypath/analysis/oddball_task/stability/';
%filelist=dir([folder '*zscored*dscalar*']);
sublist=['0001'; '0002';'0003';'0004';'0005';'0006';'0007';'0008';'0010';'0011']; % subject identifiers in BIDS data are not eh same as PB numbers
for l=1:10
    SUB=sublist(l,:)
    %load .mat file with accumulated data from permutations
    load([output 'half1/ION' SUB '/sub-' SUB '_matrix_all_perm_h1.mat']);
    load([output 'half2/ION' SUB '/sub-' SUB '_matrix_all_perm_h2.mat']);
    % read in one random dscalar for cifti header info
    b=cifti_read([folder 'sub-' SUB '_acq-3T2mm_contrast_oddball_glover_zscored.dscalar.nii']);
    
    %% find maps that overlap between the two halves and then look at the overlap between those overlaps
    for n=1:size(betamat,2)
    
        data1=betamat(:,n);
        data2=betamat2(:,n);
        thresh_high1=prctile((data1(data1>0)),75); %25 percent highest betas
        datahigh1=find(data1>=thresh_high1);
        data_matrix1=nan(size(data1));
        data_matrix1(datahigh1)=data1(datahigh1);
    
        thresh_high2=prctile((data2(data2>0)),75); %25 percent highest betas
        datahigh2=find(data2>=thresh_high2);
        data_matrix2=nan(size(data2));
        data_matrix2(datahigh2)=data2(datahigh2);

        binary_datamat1=zeros(size(data_matrix1)); % binarize matrix for dice overlap
        binary_datamat1(~isnan(data_matrix1),1)=1;
        binary_datamat2=zeros(size(data_matrix2));
        binary_datamat2(~isnan(data_matrix2),1)=1;
        overlap_dice(n,1)=dice(binary_datamat1, binary_datamat2);
    
    % use overlap across split-halves to also look at overlap of this overlap across all permutations
        for p=1:size(data_matrix1,1)
            row=[data_matrix1(p,1),data_matrix2(p,1)];
            rowidx=find(~isnan(row));
            rowpct=size(rowidx,2)/size(row,2);
            if rowpct<1
                overlap=0;
            elseif rowpct==1
                overlap=1;
            end
            halves_overlap(p,1)=overlap;
        end
        overlap_all(:,n)=halves_overlap;
    end
    %expless the overlap in percent - percent of permutations that contain this overlap pattern
    for p=1:size(overlap_all,1)
        row=overlap_all(p,:);
        rowidx=find(row>0);
        rowpct=size(rowidx,2)*100/size(row,2);
        pct_stabilitypos(p,1)=rowpct;
    end
    %write out matrix that is used for violon plot
    writematrix(overlap_dice, [output 'sub-' SUB '_dice_vector_top_25pct.csv'])

%write out dscalar for display of percent overlap
b.cdata=pct_stabilitypos;
cifti_write(b, [output 'sub-' SUB '_acq-3T2mm_25percentile_stability_map_split_half.dscalar.nii']);
end

%% dice coefficient to halves of other participants

sublist=['0001'; '0002';'0003';'0004';'0005';'0006';'0007';'0008';'0010';'0011'];
for j=1:size(sublist)
    SUB=sublist(j,:)
    stabilitypath=['/mypath/analysis/oddball_task/stability/half1/ION' SUB '/'];
    load([stabilitypath 'sub-' SUB '_matrix_all_perm_h1.mat'])
    betamat_all(:,:,j)=betamat;
end

% dice overlap between hall first halves of all subjects
for l=1:size(betamat_all,2)
    for j=1:size(sublist)
    data1=squeeze(betamat_all(:,l,j));
    thresh_high1=prctile((data1(data1>0)),75); %25 percent highest betas
    datahigh1=find(data1>=thresh_high1);
    data_matrix1=nan(size(data1));
    data_matrix1(datahigh1)=data1(datahigh1);
    binary_datamat1=zeros(size(data_matrix1));
    binary_datamat1(~isnan(data_matrix1),1)=1;
    binary_data_mat_all_subj(:,j)=binary_datamat1;
    end
    for j=1:size(sublist)
        for m=1:size(sublist)
             dice_overlap_vec(l,j, m)=dice(binary_data_mat_all_subj(:,j), binary_data_mat_all_subj(:,m));
        end
    end
end




%% save files
%write out matrix that is used for violon plot

sublist=['0001'; '0002';'0003';'0004';'0005';'0006';'0007';'0008';'0010';'0011'];
for j=1:size(sublist)
    SUB=sublist(j,:)
    diceall=squeeze(dice_overlap_vec(:,j,:));
    diceall(:,j)=[];
    stabilitypath=['/mypath/analysis/oddball_task/stability/half1/ION' SUB '/'];
    writematrix(diceall, [stabilitypath 'sub-' SUB '_dice_overlap_vector_all_subj.csv'])
end

