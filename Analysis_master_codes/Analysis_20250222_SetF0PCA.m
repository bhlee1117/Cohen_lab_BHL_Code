
clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'C5:Q175');

save_dir='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Figures/invivoPrism/FigureOptopatch';
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);
StimROI=raw(:,6);
StimWfn=raw(:,7);
isGoodCell=cell2mat(raw(:,11));
PixelSize=cell2mat(cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false));
refROI=cellfun(@(x) (str2num(num2str(x))),raw(:,14),'UniformOutput',false);
maintrunkROI=cellfun(@(x) (str2num(num2str(x))),raw(:,15),'UniformOutput',false);
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
[~, unqInd] = unique([Mouse NeuronInd] ,'row');
set(0,'DefaultFigureWindowStyle','docked')

%%
[a, unqInd] = unique([Mouse NeuronInd] ,'row');

for i=unqInd(46:end)'
    cat_trace=[];
    cat_spike=[];
    i
    SameCellInd=find(Mouse==Mouse(i) & NeuronInd==NeuronInd(i));
    for j=SameCellInd'
        load(fullfile(fpath{j},'OP_Result.mat'),'Result')
        if isfield(Result,'spike')
        cat_spike=[cat_spike Result.spike];
        else
            tr=Result.traces_bvMask-movmedian(Result.traces_bvMask,100,2);
        cat_spike=[cat_spike find_spike_bh(tr./get_threshold(tr,1),5,3)];    
        end
        if isfield(Result,'normTraces')
        cat_trace=[cat_trace Result.normTraces];
        else if isfield(Result,'traces_bvMask')
                if ~isempty(Result.traces_bvMask)
        cat_trace=[cat_trace Result.traces_bvMask./get_threshold(Result.traces_bvMask,1)];
                else
        cat_trace=[cat_trace Result.traces./get_threshold(Result.traces,1)];                
                end
        else
        cat_trace=[cat_trace Result.traces./get_threshold(Result.traces,1)];    
        end
        end
    end
    
    cat_sub=get_subthreshold(cat_trace,cat_spike(1,:),7,15);
    F0_PCA=get_F0PCA(cat_sub);

    for j=SameCellInd'
        load(fullfile(fpath{j},'OP_Result.mat'),'Result')
        Result.F0_PCA=F0_PCA;
        save([fpath{j} '/OP_Result.mat'],'Result','-v7.3')
    end

end


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