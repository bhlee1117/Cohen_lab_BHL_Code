
clear
clc;
cd '/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'C5:Q196');

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

%% Get F0_PCA
[a, unqInd] = unique([Mouse NeuronInd] ,'row');

for i=unqInd([8])'
    i
    cat_trace=[];
    cat_spike=[];
    SameCellInd=find(Mouse==Mouse(i) & NeuronInd==NeuronInd(i));
    for j=SameCellInd'
        disp([num2str(j) 'th path is loaded'])
        load(fullfile(fpath{j},'OP_Result.mat'),'Result')
        if isfield(Result,'spike')
        cat_spike=[cat_spike Result.spike];
        else
            tr=Result.traces_bvMask-movmedian(Result.traces_bvMask,100,2);
        cat_spike=[cat_spike find_spike_bh(tr./get_threshold(tr,1),5,3)];    
        end
        if isfield(Result,'normTraces')
        cat_trace=[cat_trace Result.normTraces];
        else 
        Result.normTraces=Result.traces_bvMask./get_threshold(Result.traces_bvMask,1);
        save([fpath{j} '/OP_Result.mat'],'Result','-v7.3')
        cat_trace=[cat_trace Result.normTraces];
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

%% Spike detection & checking

CS_thres=[7 1.5];
Sp_thres=[4 2];
exclude_frq2=[55.5 56]; %motion
freq_lowhigh2=exclude_frq2/(1000/2);
[b2, a2] = butter(4, freq_lowhigh2, 'stop');

for f=[84]
    f
    load(fullfile(fpath{f},'OP_Result.mat'),'Result');

    if sum(ismember(refROI{f},1))==0
        refROI{f}=[refROI{f} 1];        
    end

    Result.SpClass=[];
    Result.CStrace=[];

    noi=[1:size(Result.ftprnt,3)];
    noi_dist=ismember(Result.dist_order,noi);
    [nROI nTime]=size(Result.traces_bvMask);
    Result.spike=zeros(nROI,nTime);
 
    % Spike detection
    Result.normTraces=Result.traces_bvMask-prctile(Result.traces_bvMask,25,2);
    Result.normTraces=Result.normTraces./get_threshold(Result.normTraces,1);
    for n=1:size(Result.ftprnt,3)
        Result.normTraces(n,:) = filtfilt(b2, a2, Result.normTraces(n,:));
    end
    normTr=Result.normTraces./Result.F0_PCA;
    normTr=normTr(Result.dist_order(noi_dist),:);

    tr_ref=Result.normTraces(refROI{f},:);
    tr_ref_mean=mean(tr_ref,1,'omitnan');
    tr_ref_mean=tr_ref_mean./get_threshold(tr_ref_mean,1);
    sp_soma=find_spike_bh(tr_ref_mean-movprc(tr_ref_mean,200,30,2),Sp_thres(1),Sp_thres(2));

    sp_ref_tmp=sum(sp_soma,1)>0;
    sp_ref_1st_ind=find((sp_ref_tmp(2:end)-sp_ref_tmp(1:end-1))==1)+1;
    sp_ref_sumcon=zeros(1,size(sp_soma,2));
    sp_ref_sumcon(1,sp_ref_1st_ind)=sum(sum(reshape(sp_soma(:,sp_ref_1st_ind'+[-1:0]),size(sp_soma,1),[],2),3),1);
    sp_ref=sp_ref_sumcon;

    [~, shift]=max(reshape(tr_ref_mean(find(sp_ref>=1)+[-1:0]'),2,[]),[],1);
    shift=shift-2;
    sp_time_Soma = find(sp_ref>=1)+shift;
    sp_soma=zeros(1,nTime);
    sp_soma(sp_time_Soma)=1;
    sp_soma=[0 (sp_soma(2:end)-sp_soma(1:end-1))==1]; %remove consecutive spikes
    Result.spike(refROI{f},:)=repmat(sp_soma,length(refROI{f}),1);
    
    SpMat=Result.spike(1,:);
    SpMat_modified=interactiveFrameCheck(normTr,SpMat,[200]); %check spike detection
    sp=SpMat_modified;
    [~, shift]=max(reshape(tr_ref_mean(find(sp>0)+[-1:0]'),2,[]),[],1);
    shift=shift-2;
    sp_time = find(sp>0)+shift;
    sp_new=zeros(1,length(tr_ref_mean));
    sp_new(sp_time)=1;
    Result.spike(refROI{f},:)=repmat(sp_new,length(refROI{f}),1);
    sp_soma=Result.spike(1,:);

    % Find complex spike
    tr_sub=mean(tr_ref_mean,1)-movprc(mean(tr_ref_mean,1),200,20,2);
    tr_sub=get_subthreshold(tr_sub,sp_soma,7,17);

    [trans tr_trace]=detect_transient2(tr_sub,CS_thres,sp_soma,20);
    if isempty(trans.amp)
    CS_trace=zeros(1,nTime);
    CS_spike_time=[];
    else
    transcand=cell2mat(cellfun(@(x) length(x)>1,trans.ISI,'UniformOutput',false));
    meanISI_frnt=cellfun(@(x) mean(x(1:1)),trans.ISI(transcand));
    meanISI_first3=zeros(1,length(trans.length));
    meanISI_first3(transcand)=meanISI_frnt;

    %CS_ind=find(trans.spike_number>2 & trans.mean_ISI<15);
    CS_ind=find(trans.spike_number>2 & meanISI_first3<20);
    CS_trace=ismember(tr_trace,CS_ind);
    CS_spike=sp_soma.*bwlabel(CS_trace);
    [~, CS_spike_time]=unique(CS_spike);
    end

    sp_total=max([sp_soma; sp(2:end,:)],[],1);
    bAP_ind=zeros(1,nTime);
    bAP_ind(unique(find(sp_soma)'+[0:3]))=1;

    SpikeClassMat=zeros(3,nTime);
    SpikeClassMat(1,:)=sp_soma.*(1-CS_trace); %bAPs
    if ~isempty(CS_spike_time)
    SpikeClassMat(2,CS_spike_time(2:end))=1; %Complex spikes
    else
        SpikeClassMat(2,:)=0;
    end
    SpikeClassMat(3,:)=sp_total.*(1-bAP_ind); %dSpikes

    Result.SpClass=SpikeClassMat;
    SpMat=Result.SpClass(1:2,:);
    SpMat_modified=interactiveFrameCheck(normTr,SpMat,[200]);
    Result.SpClass(1:2,:)=SpMat_modified;

    Result.CStrace=CS_trace;
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
    for sp_type=1:2
        sp=max(Result.SpClass(sp_type,:),[],1);
        [~, shift]=max(reshape(tr_ref_mean(find(sp>0)+(-1:0)'),2,[]),[],1);
        shift=shift-2;
        sp_time = find(sp>0)+shift;
        sp_new=zeros(1,length(tr_ref_mean));
        sp_new(sp_time)=1;
        Result.SpClass(sp_type,:)=sp_new;
    end
    notdSP=find(Result.SpClass(3,:)>0 & max(Result.SpClass([1:2],:),[],1)>0);
    Result.SpClass(3,notdSP)=0;

    newSpclassSp=find((Result.SpClass(1,:)-Result.spike(refROI{f}(1),:))==1);
    Result.spike(refROI{f},newSpclassSp)=1;

    figure(f);
    clf; tiledlayout(3,1)
    ax1=nexttile([2 1]);
    imagesc(normTr,[-0.005 0.025]);
    colormap(turbo);
    ax2=nexttile([1 1]);
    plot(Result.SpClass'); hold all
    plot(Result.spike(refROI{f}(1),:)*1.5,'r'); hold all
    plot(Result.CStrace);
    linkaxes([ax1 ax2],'x')
    legend({'SS','CS','dS','CStr'})

    save(fullfile(fpath{f},'OP_Result.mat'),'Result','-v7.3');
    % load(fullfile(fpath{f},'OP_Result.mat'),'Result');
    % backupServer(fpath{f},'BHL18TB_D2','cohen_lab/Lab/Labmembers/Byung Hun Lee/Data','OP_Result.mat')
end
%%

f=96;
 load(fullfile(fpath{f},'OP_Result.mat'),'Result');
  figure(f);
    clf; tiledlayout(4,1)
    ax1=nexttile([2 1]);
        noi=[1:size(Result.ftprnt,3)];
    noi_dist=ismember(Result.dist_order,noi);
    normTr=Result.normTraces./Result.F0_PCA;
    normTr=normTr(Result.dist_order(noi_dist),:);
    imagesc(normTr,[-0.005 0.025]);
    colormap(turbo);
    ax2=nexttile([1 1]);
    plot(Result.SpClass'); hold all
    plot(Result.spike(refROI{f}(1),:)*1.5,'r'); hold all
    plot(Result.CStrace);

    ax3=nexttile([1 1]);
    plot(Result.Blue);
    linkaxes([ax1 ax2 ax3],'x')
    legend({'SS','CS','dS','CStr'})
