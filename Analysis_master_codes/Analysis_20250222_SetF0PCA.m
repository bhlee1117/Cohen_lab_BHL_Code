%
% clear;
% cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism';
% [~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
%     'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'C5:Q175');
%
% save_dir='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Figures/invivoPrism/FigureOptopatch';
% fpath=raw(:,1);
% Mouse=cell2mat(raw(:,2));
% NeuronInd=cell2mat(raw(:,5));
% CamType=raw(:,3);
% StructureData=raw(:,10);
% StimROI=raw(:,6);
% StimWfn=raw(:,7);
% isGoodCell=cell2mat(raw(:,11));
% PixelSize=cell2mat(cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false));
% refROI=cellfun(@(x) (str2num(num2str(x))),raw(:,14),'UniformOutput',false);
% maintrunkROI=cellfun(@(x) (str2num(num2str(x))),raw(:,15),'UniformOutput',false);
% place_bin=150; time_segment=15000; overlap=200;
% alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
% bound=6;
% title_str={'Basal','Apical','Peri-Soma'};
% [~, unqInd] = unique([Mouse NeuronInd] ,'row');
% set(0,'DefaultFigureWindowStyle','docked')

%%

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
save_figto='/Volumes/BHL18TB_D2/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
PlaceFieldList=cellfun(@(x) (str2num(num2str(x))),raw(:,21),'UniformOutput',false);
PlaceFieldBin=cellfun(@(x) (str2num(num2str(x))),raw(:,22),'UniformOutput',false);
set(0,'DefaultFigureWindowStyle','docked')
%foi=[1 4 5 6 8 10 11 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27];
foi=[1 4 5 6 8 10 11 15 16 17 18 19 20 21 22 23 24 25 26 27];
%foi=23;

%%
[a, unqInd] = unique([Mouse NeuronInd] ,'row');

for i=unqInd([57:end])'
    cat_trace=[];
    cat_spike=[];
    cat_ind=[];
    i
    SameCellInd=find(Mouse==Mouse(i) & NeuronInd==NeuronInd(i));
    for j=SameCellInd'
        j
        load(fullfile(fpath{j},'OP_Result.mat'),'Result')
        cat_ind=[cat_ind repmat(j,1,size(Result.traces,2))];
        if isfield(Result,'spike')
            cat_spike=[cat_spike Result.spike(1,:)];
        else
            try
                tr=Result.traces_bvMask-movmedian(Result.traces_bvMask,100,2);
                Result.normTraces=Result.traces_bvMask./get_threshold(Result.traces_bvMask,1);
            catch
                tr=Result.traces-movmedian(Result.traces,100,2);
                Result.normTraces=Result.traces./get_threshold(Result.traces,1);
            end
            cat_spike=[cat_spike find_spike_bh(tr(1,:)./get_threshold(tr(1,:),1),5,3)];
        end

        if isfield(Result,'normTraces')
            tr=Result.normTraces;
        else if isfield(Result,'traces_bvMask')
                if ~isempty(Result.traces_bvMask)
                    tr=Result.traces_bvMask./get_threshold(Result.traces_bvMask,1);
                else
                    tr=Result.traces./get_threshold(Result.traces,1);
                end
        else
            tr=Result.traces./get_threshold(Result.traces,1);
        end
        end
        tr=tr-movprc(tr,1000,30,2);
        cat_trace=[cat_trace tr];
        save([fpath{j} '/OP_Result.mat'],'Result','-v7.3')
    end
    t_spike=find(cat_spike(1,:)>0);
    [~, spamp_bin] = zscore_binning({[t_spike' cat_trace(1,t_spike)']}, [1:250:size(cat_spike,2)]);
    validInd=~isnan(spamp_bin{1}(:,2));
    y_fit=expfitDM_2(spamp_bin{1}(validInd,1),spamp_bin{1}(validInd,2),[1:size(cat_spike,2)]',1000);
    y_fit=y_fit./mean(y_fit(1:10));
    %cat_trace=cat_trace./y_fit';

    cat_sub=get_subthreshold(cat_trace,cat_spike(1,:),7,25);
    [F0_PCA N_used]=get_F0PCA(cat_sub);
    if N_used==1
        [F0_PCA]=get_F0PCA(cat_sub,3);
    end

    for j=SameCellInd'
        load(fullfile(fpath{j},'OP_Result.mat'),'Result')
        %Result.normTraces=Result.normTraces./y_fit(cat_ind==j)';
        Result.F0_PCA=F0_PCA;
        save([fpath{j} '/OP_Result.mat'],'Result','-v7.3')
    end
end
normTr=cat_trace./F0_PCA;
imagesc(normTr(Result.dist_order,:))

%% synchronize

isResult=zeros(1,length(fpath));
isTrace=zeros(1,length(fpath));
for f=1:length(fpath)
    f
    resultfilelist=dir(fullfile(fpath{f},['OP_Result*']));
    if ~isempty(resultfilelist)
        clear Result
        load(fullfile(fpath{f},resultfilelist(1).name),'Result')
        % save(fullfile(fpath{f},'OP_Result.mat'),'Result','-v7.3')
        if contains(fpath{f}, 'BHL18TB_D2')
            backupServer(fpath{f},'BHL18TB_D2','cohen_lab/Lab/Labmembers/Byung Hun Lee/Data','OP_Result.mat')
        else
            backupServer(fpath{f},'cohen_lab/Lab/Labmembers/Byung Hun Lee/Data','BHL18TB_D2','OP_Result.mat')
        end
        isResult(f)=1;
        if isfield(Result,'traces')
            isTrace(f)=1;
        end
    end
end

%% Calculate F0PCA (PC data)
for f=24
    load(fullfile(fpath{f},'PC_Result.mat'),'Result')
    normTr=Result.normTraces;
    SubTr=get_subthreshold(normTr,Result.spike(1,:),7,25);
    [F0_PCA N_used]=get_F0PCA(SubTr,[1:6]);
    if N_used==1
        [F0_PCA]=get_F0PCA(cat_sub,3);
    end

    for j=SameCellInd'
        load(fullfile(fpath{j},'OP_Result.mat'),'Result')
        %Result.normTraces=Result.normTraces./y_fit(cat_ind==j)';
        Result.F0_PCA=F0_PCA;
        save([fpath{j} '/OP_Result.mat'],'Result','-v7.3')
    end
end
normTr=cat_trace./F0_PCA;
imagesc(normTr(Result.dist_order,:))

%% synchronize