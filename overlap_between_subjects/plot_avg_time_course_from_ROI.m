clear all 
close all
clc

cd('/mypath/task_analysis')
outpath='/mypath/oddball_task/';
sublist={'0001','0002', '0003', '0004', '0005', '0006', '0007', '0008', '0010', '0011'};
pbsublist={'PB0007','PB0009', 'PB0010', 'PB0011', 'PB0015', 'PB0016', 'PB0017', 'PB0018', 'PB0019', 'PB0020'};

TR=1.761;
num_frames=14; % 25 sec
xaxistics=([0:num_frames]*TR);
vertex='4368'; % each average time series is labled by the vertex that is the center of the region of interest (see README)
side='right';

%%
for k=1:size(sublist,2)

    SUB=sublist{k}
    time_series=load([outpath 'sub-' SUB '_vertex-' vertex '_' side '_average_time_series.txt']);
    events=readtable([ '/mypath/task_analysis/ION' SUB '_event_file_habi.csv']);
    event_timecourse=zeros(size(time_series));
    event_frames=round(events.onset/TR)+1;
    event_timecourse(event_frames)=1;
    % rad in outlier fiel to censor these frames in the time course
    derivspath=['/mypath/XCP-D_derivatives_task/ION' SUB '_MENORDIC_combined_task/sub-' SUB '/ses-MENORDIC/func/'];
    outliers=readtable([derivspath 'sub-' SUB '_ses-MENORDIC_task-oddball_acq-3T2mm_outliers.tsv'], "FileType","text",'Delimiter', '\t');
    outliers=logical(table2array(outliers));
    time_series(outliers)=NaN;

	%identify event snipets from tiem course
    for i=1:sum(event_timecourse)
        beginning=event_frames(i);
        ending=beginning+num_frames;
        event_snippets(i,:)=time_series(beginning:ending,1);
        event_onsets(i,:)=event_timecourse(beginning:ending,1);
    end

	%average across all event snipets
    evoked_response(k,:)=nanmean(event_snippets);

end

%% Figure
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


f3=figure
set(0, 'CurrentFigure', f3)
plot(xaxistics, evoked_response, 'LineWidth', 3)
set(gca, 'ColorOrder', newDefaultColors)
hold on
yline(0, '--')
xline(xaxistics(find((sum(event_onsets)))))
%title('all subjects')
set(f3,'Color','w')
legend(pbsublist, 'Location','eastoutside', 'FontSize',14)
ylabel('BOLD signal change')
xlabel('time')
ylim([-18 26])
ax = gca;
ax.FontSize = 14; 
box off

% all subjects - look at when first peak occurs
narrowwindow=evoked_response(:,1:7);
[m, maxtime]=max(narrowwindow, [], 2);



