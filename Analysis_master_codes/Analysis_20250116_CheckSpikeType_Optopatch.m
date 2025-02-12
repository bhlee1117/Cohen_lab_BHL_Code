
clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'B5:P175');

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

for f=75
    load(fullfile(fpath{f},'OP_Result.mat'),'Result');

    if sum(ismember(refROI{f},1))==0
        refROI{f}=[refROI{f} 1];        
    end

    noi=[1:size(Result.ftprnt,3)];
    noi_dist=ismember(Result.dist_order,noi);

    normTr=Result.normTraces./Result.F0_PCA;
    normTr=normTr(Result.dist_order(noi_dist),:);
    SpMat=Result.SpClass(1:2,:);
    SpMat_modified=interactiveFrameCheck(normTr,SpMat,[200]);

    Result.SpClass(1:2,:)=SpMat_modified;
    CS_trace_add=tovec(find(SpMat_modified(2,:))+[-2:50]');
    Result.CStrace(CS_trace_add)=1;
    unmarkspike=find((SpMat(1,:)-SpMat_modified(1,:))==1 & Result.CStrace==0);
    markspike=find((SpMat_modified(1,:)-SpMat(1,:))==1);

    Result.spike(refROI{f},unmarkspike)=0;
    Result.spike(refROI{f},markspike)=1;
    Result.spike(refROI{f},find(max(Result.SpClass(1:2,:),[],1)>0))=1;

    % remove CStrace without CS
    bwCS=bwlabel(Result.CStrace);
    ValidCS=unique(bwCS.*Result.SpClass(2,:));
    ValidCS=ValidCS(2:end);
    Result.CStrace=ismember(bwCS,ValidCS);

    bwCS=bwlabel(Result.CStrace);
    newCS=zeros(1,size(bwCS,2));
    for b=1:max(bwCS)
        frmb=find(bwCS==b);
        CS1st=find(Result.SpClass(2,:)==1 & bwCS==b,1);
        newCS(CS1st:frmb(end))=1;
    end
    Result.CStrace=newCS;

    % Remove consecutive spikes
    tr_ref=Result.normTraces(refROI{f}(1),:);
    for sp_type=1:2
        sp=max(Result.SpClass(sp_type,:),[],1);
        [~, shift]=max(reshape(tr_ref(find(sp>0)+[-1:0]'),2,[]),[],1);
        shift=shift-2;
        sp_time = find(sp>0)+shift;
        sp_new=zeros(1,length(tr_ref));
        sp_new(sp_time)=1;
        Result.SpClass(sp_type,:)=sp_new;
    end
    notdSP=find(Result.SpClass(3,:)>0 & max(Result.SpClass([1:2],:),[],1)>0);
    Result.SpClass(3,notdSP)=0;

    newSpclassSp=find((Result.SpClass(1,:)-Result.spike(refROI{f}(1),:))==1);
    Result.spike(refROI{f},newSpclassSp)=1;
    SpMat=Result.spike(1,:);
    SpMat_modified=interactiveFrameCheck(normTr,SpMat,[100]);

    tr_ref=Result.normTraces(refROI{f}(1),:);
    sp=SpMat_modified;
    [~, shift]=max(reshape(tr_ref(find(sp>0)+[-1:0]'),2,[]),[],1);
    shift=shift-2;
    sp_time = find(sp>0)+shift;
    sp_new=zeros(1,length(tr_ref));
    sp_new(sp_time)=1;
    Result.spike(refROI{f},:)=repmat(sp_new,length(refROI{f}),1);

    sp=max(Result.spike(refROI{f},:),[],1);
    [~, shift]=max(reshape(tr_ref(find(sp>0)+[-1:0]'),2,[]),[],1);
    shift=shift-2;
    sp_time = find(sp>0)+shift;
    sp_new=zeros(1,length(tr_ref));
    sp_new(sp_time)=1;
    Result.spike(refROI{f},:)=repmat(sp_new,length(refROI{f}),1);

    clf; tiledlayout(3,1)
    ax1=nexttile([2 1]);
    imagesc(normTr,[-0.003 0.015]);
    colormap(turbo);
    ax2=nexttile([1 1]);
    plot(Result.SpClass'); hold all
    plot(Result.spike(refROI{f}(1),:)*1.5,'r'); hold all
    plot(Result.CStrace);
    linkaxes([ax1 ax2 ax3],'x')
    legend({'SS','CS','dS','CStr'})

    save(fullfile(fpath{f},'OP_Result.mat'),'Result','-v7.3');
    load(fullfile(fpath{f},'OP_Result.mat'),'Result');
    backupServer(fpath{f},'BHL18TB_D2','ByungHun2TB/Prism_PlaceCell_backup','OP_Result.mat')
end

