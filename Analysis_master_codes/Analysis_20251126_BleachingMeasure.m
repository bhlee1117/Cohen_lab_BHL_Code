%% Load file path
clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:AA31');

ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,9),'UniformOutput',false);

oblique_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,10),'UniformOutput',false);
PeriSoma_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,11),'UniformOutput',false);
basal_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false);
apical_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,13),'UniformOutput',false);
distal_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,14),'UniformOutput',false);

fpath=raw(:,1)';
StructureData=raw(:,8);
BadROI=cellfun(@(x) (str2num(num2str(x))),raw(:,17),'UniformOutput',false);
EndFrame=cell2mat(raw(:,15));
ifmotionReject=cell2mat(raw(:,16));
ifdirtRemov=cell2mat(raw(:,18));
Pixelsize=cell2mat(raw(:,6));
NSeesawComponent=cell2mat(raw(:,25));
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
PlaceFieldList=cellfun(@(x) (str2num(num2str(x))),raw(:,21),'UniformOutput',false);
PlaceFieldBin=cellfun(@(x) (str2num(num2str(x))),raw(:,22),'UniformOutput',false);
set(0,'DefaultFigureWindowStyle','docked')
%foi=[1 4 5 6 8 10 11 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27];
foi=[1 4 5 6 8 10 11 15 16 17 18 19 20 21 22 23 24 25 26 27];
%%
tc=[];
figure(2); clf;
for f=foi
    load(fullfile(fpath{f}, 'PC_Result.mat'), 'Result')
    if ~isfield(Result,'lwpass_fit')
        Ftrace=movmedian(-mean(Result.traces,1,'omitnan'),1000,2,'omitnan');
    else
        Ftrace=-Result.lwpass_fit;
    end

    t_fit=find(~isnan(Ftrace));
    [y_fit tc(f)]=expfitDM_2(t_fit',Ftrace(t_fit)',[1:600000]',100000);
    nexttile([1 1])
    plot(Ftrace); hold all
    plot(y_fit,'r');
    title(f); drawnow;
    Bleachingat10min(f)=y_fit(end)/y_fit(1);
end
%%
BleachingM=mean(Bleachingat10min(foi),'omitnan');
BleachingS=std(Bleachingat10min(foi),'omitnan');
disp(['Bleaching at 10 min, mean: ' num2str(BleachingM,3) ', s.d.: ' num2str(BleachingS,3)]);

