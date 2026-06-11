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
%foi=23;
%% Load Result files
win_pre = 30; win_post = 20;
clear PCresult

for f = [foi]
    f
    load(fullfile(fpath{f}, 'PC_Result.mat'), 'Result')

    if ~isempty(find(ismember(find(max(Result.SpClass([1 2],:),[],1)>0),find(Result.spike(1,:)>0))==0))
        diffind=find(ismember(find(max(Result.SpClass([1 2],:),[],1)>0),find(Result.spike(1,:)>0))==0);
        ss=find(max(Result.SpClass([1 2],:),[],1)>0);
        disp([num2str(ss(diffind)) ' th frame is not matching']);
        error('Spclass and spike mat is not matching');
    end

    Result.spike = Result.spike > 0;
    Result.SpClass = Result.SpClass > 0;
    PCresult{f}.Dist_order = Result.BrancDist_order;
    PCresult{f}.BranchLabel = Result.BranchLabel;
    PCresult{f}.interDendDist = Result.interDendDist * Pixelsize(f);
    PCresult{f}.pixelsize = Pixelsize(f);
    PCresult{f}.swc=Result.SWC;
    Dsign = ones(1, size(PCresult{f}.interDendDist, 1));
    Dsign(PCresult{f}.Dist_order(1:find(PCresult{f}.Dist_order == 1) - 1)) = -1;
    perisomaROI = setdiff(find(PCresult{f}.interDendDist(1, :) < 40), BadROI{f});
    noi = setdiff(1:size(Result.ftprnt, 3), BadROI{f});
    PCresult{f}.noi_dist = ismember(PCresult{f}.Dist_order, noi);

    PCresult{f}.sortDist=PCresult{f}.Dist_order(PCresult{f}.noi_dist);
    dendaxis = PCresult{f}.interDendDist(1,:) .* Dsign;
    PCresult{f}.dendaxis = dendaxis(PCresult{f}.sortDist);

    roisD(f,:)={basal_ROI{f},PeriSoma_ROI{f},apical_ROI{f},oblique_ROI{f},distal_ROI{f}}; %set the ROIs
    for dClass=1:size(roisD,2) % 1:basal, 2:soma, 3:trunk, 4:obliqu, 5:distal
        g=1;
        if ~isnan(roisD{f,dClass})
            for d=roisD{f,dClass} %branch label
                dind=setdiff(find(Result.BranchLabel==d),BadROI{f});
                PCresult{f}.roisD_order{dClass,g}=ismember(PCresult{f}.sortDist,dind)*d;
                PCresult{f}.roi_dClass(ismember(PCresult{f}.sortDist,dind))=dClass;
                PCresult{f}.branch_dClass(ismember(PCresult{f}.sortDist,dind))=dClass;
                g=g+1;
            end
        end
    end

    NormalizedTrace = Result.normTraces ./ Result.F0_PCA;
    [bAP_STA MatbAP]= get_STA(NormalizedTrace, Result.spike(1,:), win_pre, win_post);
    bAP_STA = bAP_STA - prctile(bAP_STA, 30, 2);
    PCresult{f}.MatbAP=permute(MatbAP(PCresult{f}.sortDist,:,:),[1 3 2]);
    [PCresult{f}.SpikeHeight] = max(mean(bAP_STA(perisomaROI,:), 1), [], 2);
    [~, PCresult{f}.SpikeDelay] = max(bAP_STA, [], 2);
    PCresult{f}.SpikeDelay = PCresult{f}.SpikeDelay - (win_pre + 1);

    NormalizedTrace_dirt = NormalizedTrace;
    NormalizedTrace_dirt(:, Result.motionReject > 0) = NaN;
    NormTrCheck = cellfun(@(x) x ./ Result.F0_PCA, Result.norm_trace_check, 'UniformOutput', false);
    NormTrCheck{1}(:, Result.motionReject > 0) = NaN;
    NormTrCheck{2}(:, Result.motionReject > 0) = NaN;
    if ifdirtRemov(f)
        NormalizedTrace_dirt(Result.dirtTrace > 0) = NaN;
        NormTrCheck{1}(Result.dirtTrace > 0) = NaN;
        NormTrCheck{2}(Result.dirtTrace > 0) = NaN;
    end
    NormalizedTrace_dirt = NormalizedTrace_dirt(PCresult{f}.sortDist, :);
    NormTrCheck{1} = NormTrCheck{1}(PCresult{f}.sortDist, :);
    NormTrCheck{2} = NormTrCheck{2}(PCresult{f}.sortDist, :);
    [nROI nTime]=size(NormalizedTrace_dirt);

    PCresult{f}.NormalizedTrace_dirt = NormalizedTrace_dirt;
    PCresult{f}.NormalizedTrace_ch = NormTrCheck;
    PCresult{f}.Subthreshold = get_subthreshold(NormalizedTrace_dirt, max(Result.spike(1,:), [], 1) > 0, 7, 17);
    PCresult{f}.Subthreshold(isnan(NormalizedTrace_dirt)) = NaN;

    sptime = find(Result.spike(1,:));
    brst = bwlabel((sptime(2:end) - sptime(1:end-1)) < 20);
    SpClass = Result.SpClass;
    BS_trace = zeros(1, nTime);
    CS_trace= Result.CStrace;
    bwCStrace= bwlabel(CS_trace);
    CS_fistSp=find(SpClass(2,:)>0);
    if any(CS_trace(CS_fistSp)==0)
        CS_trace(CS_fistSp(find(CS_trace(CS_fistSp)==0))'+[0:50])=1;
    end

    SubthInterp=interpolateNaN(mean(PCresult{f}.Subthreshold(roisD{f,2},:),1,'omitnan'))./PCresult{f}.SpikeHeight; %soma subthreshold
    bwCScand=bwlabel(SubthInterp>0.45); % Subthreshold level higher than 0.3 spike height is included to complex spike interval
    CS_trace_new=CS_trace;
    for b=1:max(bwCStrace)
        candInd=unique(bwCScand(bwCStrace==b));
        candInd(candInd==0)=[];
        CS_trace_new(ismember(bwCScand,candInd))=1;
    end
    % brst_fistSp=[];
    % for b=1:max(brst)
    % brst_fistSp(b)=sptime(find(brst==b,1));
    % end
    CS_trace_new(CS_fistSp-1)=0;    %do not connect two CSs
    %CS_trace_new(brst_fistSp-1)=0;    %do not connect CSs and BSs

    bwCStrace= bwlabel(CS_trace_new);
    for b=1:max(bwCStrace)
        if any(SpClass(2,bwCStrace==b)>0)
        else
            CS_trace_new(bwCStrace==b)=0; %remove the CS traces without spikes
        end
    end
    CS_trace=CS_trace_new;
    bwCStrace= bwlabel(CS_trace);

    for b = 1:max(brst)
        bind=find(brst == b);
        bwn = sptime([bind]);
        if any(CS_trace(bwn)) % part of CS?
            CSind=unique(bwCStrace(bwn));
            CSind(CSind==0)=[];
            if length(CSind)>1
                error(['Spike train in more than one CS, check frame #' num2str(bwn)]);
            end
            if SpClass(2,bwn(1))==0 & CS_trace(bwn(1))==0 % the first spike is not assigned to CS?
                originalCSinterval=find(bwCStrace==CSind);
                SpClass(:,bwn(1):originalCSinterval(end))=0;
                SpClass(2,bwn(1))=1;
                CS_trace(bwn(1):originalCSinterval(end))=1;
            end
        else
            SpClass(1, bwn(1):sptime(bind(end)+1)) = 0;
            SpClass(4, bwn(1)) = 1; % assign first spike of BS
            BS_trace(1, [bwn(1):sptime(bind(end)+1)+7]) = b;
        end
    end

    BS_trace=BS_trace(1:nTime); CS_trace=CS_trace(1:nTime);
    SpClass = SpClass([1 2 4 3], :); %SS, CS, BS ,dSP
    Classvec = get_Class2index(SpClass);
    SpikeClassvec = Classvec .* Result.spike(1,:);

    PCresult{f}.SpikeMat = double(Result.spike>0);
    PCresult{f}.SpikeMat(:, Result.motionReject > 0) = NaN;
    if ifdirtRemov(f)
        PCresult{f}.SpikeMat(Result.dirtTrace > 0) = NaN;
    end
    PCresult{f}.SpikeClassMat = SpClass; % first spike of class
    PCresult{f}.SpikeClassVecMat = SpikeClassvec; % all the spikes in class

    BrstOrder = get_BurstOrder(Result.spike(1,:), 20) - SpClass(1,:);
    BrstOrder(SpClass(3,:) > 0) = 1;
    ComplexSpikeOrder = get_spikeOrder(CS_trace, Result.spike(1,:));
    PCresult{f}.SpikeOrder=max([BrstOrder; ComplexSpikeOrder; SpClass(1,:)],[],1);

    ActivityTr=max([Result.spike(1,:); BS_trace; CS_trace],[],1)>0;
    PCresult{f}.Subthreshold_larger = get_subthreshold(NormalizedTrace_dirt, ActivityTr, 40, 500);
    PCresult{f}.Subthreshold_larger(isnan(NormalizedTrace_dirt)) = NaN;

    for ch=1:2
        PCresult{f}.Subthreshold_larger_ch{ch}= get_subthreshold(PCresult{f}.NormalizedTrace_ch{ch}, ActivityTr, 40, 500);
        PCresult{f}.Subthreshold_larger_ch{ch}(isnan(NormalizedTrace_dirt)) = NaN;
    end

    % for stype = 1:3
    %     s_list = find(SpClass(stype,:) > 0);
    %     [~, PCresult{f}.SpClassTrSp{stype}] = get_STA(Result.spike(1,:), SpClass(stype,:), -nTau(1), nTau(end));
    %     [~, PCresult{f}.SpClassTrBlue{stype}] = get_STA(Result.Blue, SpClass(stype,:), -nTau(1), nTau(end));
    %     [~, PCresult{f}.SpClassTrCStr{stype}] = get_STA(Result.CStrace, SpClass(stype,:), -nTau(1), nTau(end));
    %     [~, PCresult{f}.SpClassTrOrder{stype}] = get_STA(BrstOrder, SpClass(stype,:), -nTau(1), nTau(end));
    %
    %     switch stype
    %         case 1
    %             [~, PCresult{f}.SpClassTrV{stype,1} sptime]=get_STA(NormalizedTrace_dirt,SpClass(stype,:),-nTau(1),nTau(end));
    %             [~, PCresult{f}.SpClassTrSub{stype,1}] = get_STA(PCresult{f}.Subthreshold, SpClass(stype,:), -nTau(1), nTau(end));
    %             [~, PCresult{f}.SpClassTrSubLar{stype,1}] = get_STA(PCresult{f}.Subthreshold_larger, SpClass(stype,:), -nTau(1), nTau(end));
    %             PCresult{f}.SpClassTrList{stype,1}=sptime;
    %         case 2
    %             for ns=1:5
    %                 [~, PCresult{f}.SpClassTrV{stype,ns} sptime]=get_STA(NormalizedTrace_dirt,(ComplexSpikeOrder)==ns,-nTau(1),nTau(end));
    %                 [~, PCresult{f}.SpClassTrSub{stype,ns}] = get_STA(PCresult{f}.Subthreshold,(ComplexSpikeOrder)==ns,-nTau(1), nTau(end));
    %                 [~, PCresult{f}.SpClassTrSubLar{stype,ns}] = get_STA(PCresult{f}.Subthreshold_larger, (ComplexSpikeOrder)==ns, -nTau(1), nTau(end));
    %                 PCresult{f}.SpClassTrList{stype,ns}=sptime;
    %             end
    %         case 3
    %             for ns=1:5
    %                 [~, PCresult{f}.SpClassTrV{stype,ns} sptime]=get_STA(NormalizedTrace_dirt,(BrstOrder.*(1-Result.CStrace))==ns,-nTau(1),nTau(end));
    %                 [~, PCresult{f}.SpClassTrSub{stype,ns}] = get_STA(PCresult{f}.Subthreshold,(BrstOrder.*(1-Result.CStrace))==ns,-nTau(1), nTau(end));
    %                 [~, PCresult{f}.SpClassTrSubLar{stype,ns}] = get_STA(PCresult{f}.Subthreshold_larger,(BrstOrder.*(1-Result.CStrace))==ns, -nTau(1), nTau(end));
    %                 PCresult{f}.SpClassTrList{stype,ns}=sptime;
    %             end
    %     end
    % end
    %
    % PCresult{f}.SpClassTrV=cellfun(@(x) permute(x,[1 3 2]),PCresult{f}.SpClassTrV,'UniformOutput',false);
    % PCresult{f}.SpClassTrSp=cellfun(@(x) permute(x,[1 3 2]),PCresult{f}.SpClassTrSp,'UniformOutput',false);
    % PCresult{f}.SpClassTrBlue=cellfun(@(x) permute(x,[1 3 2]),PCresult{f}.SpClassTrBlue,'UniformOutput',false);
    % PCresult{f}.SpClassTrSub=cellfun(@(x) permute(x,[1 3 2]),PCresult{f}.SpClassTrSub,'UniformOutput',false);
    % PCresult{f}.SpClassTrSubLar=cellfun(@(x) permute(x,[1 3 2]),PCresult{f}.SpClassTrSubLar,'UniformOutput',false);
    % PCresult{f}.SpClassTrOrder=cellfun(@(x) permute(x,[1 3 2]),PCresult{f}.SpClassTrOrder,'UniformOutput',false);
    % PCresult{f}.SpClassTrCStr=cellfun(@(x) permute(x,[1 3 2]),PCresult{f}.SpClassTrCStr,'UniformOutput',false);

    PCresult{f}.BlueStim = Result.Blue;
    PCresult{f}.CStrace = CS_trace;
    PCresult{f}.BStrace = BS_trace;
    PCresult{f}.avgImg = Result.ref_im;
    PCresult{f}.Ftprnts = Result.ftprnt;
    PCresult{f}.VR = Result.VR;

    % Detect peaks
    StopFreq=[12 200]; %filter high frequency
    [nROI nTime]=size(PCresult{f}.Subthreshold);
    Subthreshold_int=interpolateNaN(PCresult{f}.Subthreshold);
    FilterTr=[];
    PCresult{f}.peakvec=zeros(nROI,nTime); PCresult{f}.troughvec=zeros(nROI,nTime);
    Npeak=[]; Ntrough=[];
    % figure(f+400); clf;
    for b=unique(Result.BranchLabel)
        branchLabel=Result.BranchLabel(PCresult{f}.sortDist);
        branchInd=find(branchLabel==b);
        %[PhaseTrace BasalSubFilt BasalthetaPower] = get_phase(Subthreshold_int(n,:), 1000, FilterFreq);
        perispikefrm=unique(find(PCresult{f}.SpikeMat(1,:))'+[-5:50]);
        perispikefrm(perispikefrm<1 | perispikefrm>nTime)=[];
        FilterTr(b,:)=get_bandstop(mean(Subthreshold_int(branchInd,:),1,'omitnan')./PCresult{f}.SpikeHeight,1000,StopFreq);
        FilterTr(b,:)=FilterTr(b,:)-movmedian(FilterTr(b,:),300,2);
        [pks, locs] = findpeaks(FilterTr(b,:),'MinPeakHeight', 0.5,'MinPeakDistance',100, ...
            'MinPeakProminence', 0.2);  % Prominence can also be tuned
        [trough, loc_trgh] = findpeaks(-FilterTr(b,:),'MinPeakHeight', 0.5,'MinPeakDistance',100, ...
            'MinPeakProminence', 0.2);  % Prominence can also be tuned
        % Turn to 0 during blue Stim and peri-spike frame
        PCresult{f}.peakvec(branchInd,:)=repmat(ind2vec(nTime,locs,1),length(branchInd),1);
        PCresult{f}.peakvec(branchInd,PCresult{f}.BlueStim>0 | ind2vec(nTime,perispikefrm,1))=0;
        PCresult{f}.troughvec(branchInd,:)=repmat(ind2vec(nTime,loc_trgh,1),length(branchInd),1);
        PCresult{f}.troughvec(branchInd,PCresult{f}.BlueStim>0 | ind2vec(nTime,perispikefrm,1))=0;
    end

    if ~isempty(PlaceFieldList{f})
        binTrack = ceil(PCresult{f}.VR(5,:) / (115 / 150));
        PFvec = zeros(1, nTime);
        for p = 1:length(PlaceFieldBin{f})/2
            if PlaceFieldBin{f}(2*p-1) > PlaceFieldBin{f}(2*p)
                Pvec = ~(binTrack < PlaceFieldBin{f}(2*p-1) & binTrack > PlaceFieldBin{f}(2*p));
            else
                Pvec = binTrack > PlaceFieldBin{f}(2*p-1) & binTrack < PlaceFieldBin{f}(2*p);
            end
            Lapvec = PCresult{f}.VR(8,:) > PlaceFieldList{f}(2*p-1) & PCresult{f}.VR(8,:) < PlaceFieldList{f}(2*p);
            PFvec = PFvec | (Lapvec & Pvec);
        end
        PCresult{f}.PFvec = PFvec;
    end

    PCresult{f}.runVec = double(imdilate(movmean(PCresult{f}.VR(end,:),200)>0.002, ones(1,2001)));
end

%% Measuring bAP amplitude

x_bin_edges=[-300 -50 50 140 460];
x_bin_edges2=[-300 -220 -150 -50 50 100 200 300 400];

nTau = [-200:250]; % SS, CS, Brst
PerisomaRange=[-70 70]; presub_readtime=[-10:-3];
distalpeakdetectionwindow=[10 40];
AlignedbAPall=[]; AlignedEvntall=[]; bAPPropsMat=[]; EventPropsMat=[]; nTau_EV=[-200 500]; AlignedEvntcell_bin=[]; AlignedbAPcell_bin=[];
AUCbAPall=[];
for f=[foi]%[8 10 11 15 16 17 18 19 20 21 22 23 24 25 26 27];
    f
    bAPPropsMat{f,1}=[]; bAPPropsMat{f,2}=[];
    EventPropsMat{f,1}=[]; EventPropsMat{f,2}=[];
    AUCbAPall{f,1}=[]; AUCbAPall{f,2}=[];
    AlignedEvntall{f,1}=[]; AlignedEvntall{f,2}=[]; AlignedEvntall{f,3}=[];
    AlignedbAPall{f,1}=[]; AlignedbAPall{f,2}=[]; LabelMat=[]; AUCall=[]; Suball=[];
    dax=PCresult{f}.dendaxis;
    [nROI nTime]=size(PCresult{f}.Subthreshold);

    BaselineTrace=PCresult{f}.Subthreshold_larger;
    EventTrace=zeros(3,nTime);
    EventTrace(1,find(PCresult{f}.SpikeClassMat(1,:))'+[-4:3])=1;
    EventTrace(3,:)=PCresult{f}.BStrace>0;
    EventTrace(2,:)=PCresult{f}.CStrace;
    EventTrace(4,:)=PCresult{f}.SpikeMat(1,:);
    %BaselineTrace=get_subthreshold(PCresult{f}.NormalizedTrace_dirt,PCresult{f}.SpikeMat(1,:)>0,50,300);
    [~, AlignedbAP SpikeList]=get_STA(PCresult{f}.NormalizedTrace_dirt,PCresult{f}.SpikeMat(1,:)>0,-nTau(1),nTau(end));
    [~, AlignedbAPch]=get_STA(PCresult{f}.NormalizedTrace_ch{1},PCresult{f}.SpikeMat(1,:)>0,-nTau(1),nTau(end));
    [~, AlignedbAPch2]=get_STA(PCresult{f}.NormalizedTrace_ch{2},PCresult{f}.SpikeMat(1,:)>0,-nTau(1),nTau(end));
    [~, AlignedBaseline]=get_STA(BaselineTrace,PCresult{f}.SpikeMat(1,:)>0,-nTau(1),nTau(end));
    [~, AlignedBaselinech]=get_STA(PCresult{f}.Subthreshold_larger_ch{1},PCresult{f}.SpikeMat(1,:)>0,-nTau(1),nTau(end));
    [~, AlignedBaselinech2]=get_STA(PCresult{f}.Subthreshold_larger_ch{2},PCresult{f}.SpikeMat(1,:)>0,-nTau(1),nTau(end));
    [~, AlignedEV EVList]=get_STA(PCresult{f}.NormalizedTrace_dirt,max(PCresult{f}.SpikeClassMat(1:3,:),[],1)>0,-nTau_EV(1),nTau_EV(end));
    [~, AlignedSub]=get_STA(PCresult{f}.Subthreshold,max(PCresult{f}.SpikeClassMat(1:3,:),[],1)>0,-nTau_EV(1),nTau_EV(end));
    [~, AlignedEVBaseline]=get_STA(BaselineTrace,max(PCresult{f}.SpikeClassMat(1:3,:),[],1)>0,-nTau_EV(1),nTau_EV(end));
    [~, AlignedEVCSTrace]=get_STA(EventTrace>0,max(PCresult{f}.SpikeClassMat(1:3,:),[],1)>0,-nTau_EV(1),nTau_EV(end));
    AlignedEVCSTrace=permute(AlignedEVCSTrace,[2 3 1]);
    AlignedbAPall{f,1}=permute(AlignedbAP,[1 3 2]); %raw
    AlignedbAPall{f,2}=permute(AlignedbAP-AlignedBaseline,[1 3 2]); %subracted
    AlignedbAPall{f,3}=permute(AlignedbAPch-AlignedBaselinech,[1 3 2]); %subracted
    AlignedbAPall{f,4}=permute(AlignedbAPch2-AlignedBaselinech2,[1 3 2]); %subracted

    AlignedEvntall{f,1}=cat(3,AlignedEvntall{f,1},permute(AlignedEV,[1 3 2])); %raw
    AlignedEvntall{f,2}=cat(3,AlignedEvntall{f,2},permute(AlignedEV-AlignedEVBaseline,[1 3 2])); %subracted
    AlignedEvntall{f,3}=cat(3,AlignedEvntall{f,3},permute(AlignedSub,[1 3 2])); %subthreshold

    SpikeOrder=PCresult{f}.SpikeOrder(SpikeList);
    IsCS=PCresult{f}.CStrace(SpikeList)>0;
    IsBlue=PCresult{f}.BlueStim(SpikeList)>0;
    IsNA=sum(isnan(PCresult{f}.NormalizedTrace_dirt(:,SpikeList)),1)>0;
    IsPF=PCresult{f}.PFvec(SpikeList);
    IsRun=PCresult{f}.runVec(SpikeList);
    IsIsolated=(movsum([diff([SpikeList])>20 1],2)>=2);
    delayMat=NaN(size(AlignedbAPall{f,1},1),1,size(AlignedbAPall{f,1},3));

    for spo=unique(SpikeOrder)
        STAbAP=mean(AlignedbAPall{f,1}(:,:,SpikeOrder==spo),3,'omitnan');
        [~, tmax]=max(STAbAP(:,-nTau(1)+[-2:4]),[],2);
        tmax=tmax-nTau(1)-3;
        delayMat(:,:,SpikeOrder==spo)=repmat(tmax,1,1,sum(SpikeOrder==spo));
    end
    [AUCbAP, AUCrawbAP kinkAmp kink_raw]=get_AUC(AlignedbAPall{f,2},delayMat,2,3);
    [AUCbAP_short, AUCrawbAP_short]=get_AUC(AlignedbAPall{f,2},delayMat,1,1);

    [~, AUCrawbAP_ch1]=get_AUC(AlignedbAPall{f,3},delayMat,2,3); %ch1
    [~, AUCrawbAP_short_ch1]=get_AUC(AlignedbAPall{f,3},delayMat,1,1);
    [~, AUCrawbAP_ch2]=get_AUC(AlignedbAPall{f,4},delayMat,2,3); %ch2
    [~, AUCrawbAP_short_ch2]=get_AUC(AlignedbAPall{f,4},delayMat,1,1);

    AUCbAP_cell = num2cell(AUCrawbAP, 1); % Change this line to use AUC vs AUC raw for Tx rate
    AUCbAP_cell = cellfun(@(x) [dax' x(:,1)],AUCbAP_cell,'UniformOutput',false);
    AmpbAP_cell = num2cell(AUCrawbAP_short, 1);
    AmpbAP_cell = cellfun(@(x) [dax' x(:,1)],AmpbAP_cell,'UniformOutput',false);
    AlignedEvntcell = cellfun(@(x) [dax' x(:,1)],num2cell(AlignedEvntall{f,2},1),'UniformOutput',false);
    AlignedbAPcell = cellfun(@(x) [dax' x(:,1)],num2cell(AlignedbAPall{f,2},1),'UniformOutput',false);

    [~, AUCbAPbin_cell, dcenter] = zscore_binning(AUCbAP_cell, x_bin_edges);
    [~, AmpbAPbin_cell, dcenter] = zscore_binning(AmpbAP_cell, x_bin_edges);
    [~, AUCbAPbin_cell2, dcenter2] = zscore_binning(AUCbAP_cell, x_bin_edges2);
    emptycell=cellfun(@isempty,AUCbAPbin_cell);
    [~, EVbin, dcenter] = zscore_binning(reshape(AlignedEvntcell,1,[]), x_bin_edges);
    EVbin=reshape(EVbin,size(AlignedEvntall{f,2},2),size(AlignedEvntall{f,2},3));
    EVbin=cellfun(@(x) x(:,2),EVbin,'UniformOutput',false);
    AlignedEvntcell_bin=reshape(cell2mat(EVbin),length(dcenter),diff(nTau_EV)+1,[]);
    [~, bAPbin, dcenter] = zscore_binning(reshape(AlignedbAPcell,1,[]), x_bin_edges);
    bAPbin=reshape(bAPbin,size(AlignedbAPall{f,2},2),size(AlignedbAPall{f,2},3));
    bAPbin=cellfun(@(x) x(:,2),bAPbin,'UniformOutput',false);
    AlignedbAPcell_bin=reshape(cell2mat(bAPbin),length(dcenter),-nTau(1)+nTau(end)+1,[]);

    centroid_bAPAUC=sum(dax'.*AUCrawbAP,1,'omitnan')./sum(AUCrawbAP,1,'omitnan');
    trnsmit_bAPAUC=NaN(length(SpikeList),1);
    trnsmit_bAPAUC(find(~emptycell))=cellfun(@(x,y) x(4,2)/y(2,2),AUCbAPbin_cell(~emptycell),AmpbAPbin_cell(~emptycell)); %Divided by short AUC
    trnsmit_bAPAUC2=NaN(length(SpikeList),1);
    trnsmit_bAPAUC2(find(~emptycell))=cellfun(@(x) x(4,2)/x(2,2),AUCbAPbin_cell(~emptycell)); %Divided by same length AUC
    bAPAUC_Apical=cellfun(@(x) x(4,2),AUCbAPbin_cell(~emptycell));
    bAPAUC_Soma=cellfun(@(x) x(2,2),AmpbAPbin_cell(~emptycell));
    ISI=[NaN diff(SpikeList)]';
    [SpType, ~]=find(PCresult{f}.SpikeClassVecMat(:,SpikeList));

    catlab=[SpikeList' SpikeOrder' IsBlue' IsCS' IsNA' IsPF' IsRun' IsIsolated' centroid_bAPAUC' trnsmit_bAPAUC bAPAUC_Apical bAPAUC_Soma ISI trnsmit_bAPAUC2 SpType];
    bAPPropsMat{f,1}=[catlab];
    bAPPropsMat{f,2}={AUCbAPbin_cell,AUCbAPbin_cell2,AmpbAPbin_cell,AUCbAP,AUCrawbAP,AUCrawbAP_short,AlignedbAPcell_bin, AUCrawbAP_ch1, AUCrawbAP_ch2, AUCrawbAP_short_ch1, AUCrawbAP_short_ch2};

    fieldName={'SpikeFrame', 'SpikeOrder','IsBlue','IsCS','IsNA','IsPF','IsRun','IsIsolated','AUC_centroid','TransmissionRate','AUC_apical','AUC_soma','ISI','TransmissionRate_equal','Spike_Type'};
    bAPPropsMat{f,1}=array2table(bAPPropsMat{f,1},'VariableNames',fieldName);

    fieldName2={'AUC_binned','AUC_binned_fine','Amp_binned','AUCbAP','AUCbAPraw','AUCbAPrawShort','AlignedbAP_binned'...
               ,'AUCbAPraw_ch1','AUCbAPraw_ch2','AUCbAPrawShort_ch1','AUCbAPrawShort_ch2'};
    bAPPropsMat{f,2}=array2table(bAPPropsMat{f,2},'VariableNames',fieldName2);
    %run Analysis_20250608_STAMovGen.m

    % Event props
    IntervalMat=[]; AUCrawEV=[]; ISIEV=[]; kink_rawEV=[]; Evtype=[]; NSpike=[]; TRevMax=[]; TRevFirst=[]; TransRatebAPinEvent=[];
    PostSub_Basal=[]; PostSub_Distal=[]; ApicalAUCFirst=[];
    PostAUC_Soma=[]; PostAUC_Distal=[]; PostPeak_Soma=[]; PostPeak_Distal=[];
    post_idx = (-nTau_EV(1)+1+distalpeakdetectionwindow(1)) : (-nTau_EV(1)+1+distalpeakdetectionwindow(2));
    somaROI   = PCresult{f}.roi_dClass == 2;
    distalROI = ismember(PCresult{f}.roi_dClass,[5]);
    for ev=1:length(EVList) %Event
        Evtype(ev)=find(PCresult{f}.SpikeClassMat([1:3],EVList(ev)));
        bwOn=bwlabel(AlignedEVCSTrace(ev,:,Evtype(ev))>0);
        EvOnInterval=find(bwOn==bwOn(-nTau_EV(1)+2));
        IntervalMat(ev)=length(EvOnInterval);
        NSpike(ev)=sum(AlignedEVCSTrace(ev,-nTau_EV(1)+[1:IntervalMat(ev)],4),2);

        EventFrame=EVList(ev)+[0:IntervalMat(ev)];
        bAPinEvent=ismember(bAPPropsMat{f,1}.SpikeFrame,EventFrame);
        TransRatebAPinEvent{ev}=bAPPropsMat{f,1}.TransmissionRate(bAPinEvent);
        ApicalAUCbAPinEvent{ev}=bAPPropsMat{f,1}.AUC_apical(bAPinEvent);
        if isempty(TransRatebAPinEvent{ev})
            TRevFirst(ev)=NaN; TRevMax(ev)=NaN;
            ApicalAUCFirst(ev)=NaN;
        else
            TRevFirst(ev)=TransRatebAPinEvent{ev}(1);
            TRevMax(ev)=max(TransRatebAPinEvent{ev});
            ApicalAUCFirst(ev)=ApicalAUCbAPinEvent{ev}(1);
        end

        [AUCEV, AUCrawEV(:,ev) kinkAmpEV kink_rawEV(:,ev)]=get_AUC(AlignedEvntall{f,2}(:,:,ev),repmat(round(-nTau_EV(1)+IntervalMat(ev)/2),nROI,1),IntervalMat(ev)+2,IntervalMat(ev));
        Interv=EVList(ev)-find(PCresult{f}.SpikeMat(1,:)>0);
        if isempty(find(Interv>0)) %First spike
            ISIEV(ev)=NaN;
        else
            ISIEV(ev)=min(Interv(find(Interv>0)));
        end
        if any((EvOnInterval(end)+[1:8])>size(AlignedEvntall{f,3},2))
        PostSub_Basal(ev,1)=NaN;
        PostSub_Distal(ev,1)=NaN;
        else
        PostSub_Basal(ev,1)=squeeze(mean(AlignedEvntall{f,3}(ismember(PCresult{f}.roi_dClass,1),EvOnInterval(end)+[1:8],ev),[1 2],'omitnan'));
        PostSub_Distal(ev,1)=squeeze(mean(AlignedEvntall{f,3}(ismember(PCresult{f}.roi_dClass,5),EvOnInterval(end)+[1:8],ev),[1 2],'omitnan'));
        end

        % Post-spike raw ΔF/F: AUC and peak amplitude (20-80 ms window)
        v_soma   = mean(AlignedEvntall{f,3}(somaROI,   post_idx, ev), 1, 'omitnan')';  % [61 x 1]
        v_distal = mean(AlignedEvntall{f,3}(distalROI, post_idx, ev), 1, 'omitnan')';

        if any(somaROI)
            vSort = sort(v_soma, 'descend');
            PostPeak_Soma(ev,1) = mean(vSort(1:min(5,end)), 'omitnan');
            PostAUC_Soma(ev,1)  = sum(v_soma, 'omitnan');
        else
            PostPeak_Soma(ev,1) = NaN;  PostAUC_Soma(ev,1)  = NaN;
        end
        if any(distalROI)
            vSort = sort(v_distal, 'descend');
            PostPeak_Distal(ev,1) = mean(vSort(1:min(5,end)), 'omitnan');
            PostAUC_Distal(ev,1)  = sum(v_distal, 'omitnan');
        else
            PostPeak_Distal(ev,1) = NaN;  PostAUC_Distal(ev,1) = NaN;
        end
    end

    PreSub=squeeze(mean(AlignedEvntall{f,3}(:,-nTau_EV(1)+presub_readtime,:),2,'omitnan'));
    PreSub_Basal=squeeze(mean(AlignedEvntall{f,3}(ismember(PCresult{f}.roi_dClass,1),-nTau_EV(1)+presub_readtime,:),[1 2],'omitnan'));
    PreSub_Soma=squeeze(mean(AlignedEvntall{f,3}(ismember(PCresult{f}.roi_dClass,2),-nTau_EV(1)+presub_readtime,:),[1 2],'omitnan'));
    PreSub_Distal=squeeze(mean(AlignedEvntall{f,3}(ismember(PCresult{f}.roi_dClass,5),-nTau_EV(1)+presub_readtime,:),[1 2],'omitnan'));
    PreSub_cell = num2cell(PreSub, 1);
    PreSub_cell = cellfun(@(x) [dax' x(:,1)],PreSub_cell,'UniformOutput',false);
    [~, PreSubbin_cell] = zscore_binning(PreSub_cell, x_bin_edges);

    AUCEV_cell = num2cell(AUCrawEV, 1);
    AUCEV_cell = cellfun(@(x) [dax' x(:,1)],AUCEV_cell,'UniformOutput',false);
    [~, AUCEVbin_cell, dcenter] = zscore_binning(AUCEV_cell, x_bin_edges);
    emptycellEV=cellfun(@isempty,AUCEVbin_cell);

    AmpEV_cell = num2cell(kink_rawEV, 1);
    AmpEV_cell = cellfun(@(x) [dax' x(:,1)],AmpEV_cell,'UniformOutput',false);
    [~, AmpEVbin_cell, dcenter] = zscore_binning(AmpEV_cell, x_bin_edges);

    centroid_EVAUC=sum(dax'.*AUCrawEV,1,'omitnan')./sum(AUCrawEV,1,'omitnan');
    trnsmit_EVAUC=NaN(length(EVList),1);
    trnsmit_EVAUC(find(~emptycellEV))=cellfun(@(x) x(4,2)/x(2,2),AUCEVbin_cell(~emptycellEV));
    EVAUC_Apical=cellfun(@(x) x(4,2),AUCEVbin_cell(~emptycellEV));
    EVAUC_Soma=cellfun(@(x) x(2,2),AUCEVbin_cell(~emptycellEV));
    EVAmp_Apical=cellfun(@(x) x(4,2),AmpEVbin_cell(~emptycellEV));
    PreSub_Apical=cellfun(@(x) x(4,2),PreSubbin_cell(~emptycellEV));

    Evlab=[EVList' PCresult{f}.BlueStim(EVList)'>0 sum(isnan(PCresult{f}.NormalizedTrace_dirt(:,EVList)),1)'>0 PCresult{f}.PFvec(EVList)'...
        PCresult{f}.runVec(EVList)' IntervalMat' centroid_EVAUC' trnsmit_EVAUC EVAUC_Apical EVAUC_Soma ISIEV' EVAmp_Apical Evtype' NSpike' PreSub_Soma PreSub_Apical TRevMax' TRevFirst' PreSub_Distal PreSub_Basal...
        PostSub_Basal PostSub_Distal ApicalAUCFirst' Evtype' ...
        PostAUC_Soma PostAUC_Distal PostPeak_Soma PostPeak_Distal];
    EventPropsMat{f,1}=[Evlab];
    EventPropsMat{f,2}={TransRatebAPinEvent,AUCrawEV,AUCEVbin_cell,AlignedEvntcell_bin};

    fieldName={'SpikeFrame','IsBlue','IsNA','IsPF','IsRun','Length','AUC_centroid','TransmissionRate','AUC_apical','AUC_soma',...
               'ISI','Amp_apical','Spike_type','Nspike','PreSub_Soma','PreSub_Apical','TRmax_bAP','TRfrst_bAP','PreSub_Distal',...
               'PreSub_Basal','PostSub_Basal','PostSub_Distal','AUCapical_frst','SpClass',...
               'PostAUC_Soma','PostAUC_Distal','PostPeak_Soma','PostPeak_Distal'};
    EventPropsMat{f,1}=array2table(EventPropsMat{f,1},'VariableNames',fieldName);
    EventPropsMat{f,2}=array2table(EventPropsMat{f,2},'VariableNames',{'bAPtransmissionRate','AUC','AUC_binned','AlignedEv_binned'});
end

% f=23;
% figure(23); clf; tiledlayout(2,2);
% [~, spsort]=sort(bAPPropsMat{f,1}.TransmissionRate,'ascend');
% nexttile(1,[1 1])
% imagesc(squeeze(max(AlignedbAPall{f,1}(:,-nTau(1)+[-1:3],spsort),[],2)));
% xlabel('Spike ID'); ylabel('Basal to apical');
% nexttile(3,[1 1])
% plot(bAPPropsMat{f,1}.SpikeOrder(spsort),'.')
% xlabel('Spike ID'); ylabel('Spike order');
for f=foi
    showInd=find(~bAPPropsMat{f,1}.IsBlue & ~bAPPropsMat{f,1}.IsNA);
    showDat=[log10(bAPPropsMat{f,1}.ISI(showInd)) bAPPropsMat{f,1}.TransmissionRate(showInd)];
    showDat((showDat(:,2)>50 | showDat(:,2)<-5),:)=[];
    showDatAA=[log10(bAPPropsMat{f,1}.ISI(showInd)) bAPPropsMat{f,1}.AUC_apical(showInd)];
    showDatAA((showDatAA(:,2)>50 | showDatAA(:,2)<-5),:)=[];

    Normconst=mean(showDat(showDat(:,1)>log10(300),2),1,'omitnan');
    NormconstAA=mean(showDatAA(showDatAA(:,1)>log10(300),2),1,'omitnan');

    NormTXratio=bAPPropsMat{f,1}.TransmissionRate/Normconst;
    NormTXratio_event_frst=EventPropsMat{f,1}.TRfrst_bAP/Normconst;
    NormTXratio_event_max=EventPropsMat{f,1}.TRmax_bAP/Normconst;

    NormAA=bAPPropsMat{f,1}.AUC_apical/NormconstAA;
    NormAA_event_frst=EventPropsMat{f,1}.AUCapical_frst/NormconstAA;

    bAPPropsMat{f,1}.TransmissionRate_norm=NormTXratio;
    bAPPropsMat{f,1}.AUC_apical_norm=NormAA;
    EventPropsMat{f,1}.TRfrst_bAP_norm=NormTXratio_event_frst;
    EventPropsMat{f,1}.TRmax_bAP_norm=NormTXratio_event_max;
    EventPropsMat{f,1}.AUCapical_frst_norm=NormAA_event_frst;
end

%% Compare somatic and apical distal voltage trace: SS vs CS vs BS

foi2show=20;
PeakAmp_spclass=[];
cmap_sc = [0.20 0.45 0.80; 0.85 0.20 0.20; 0.20 0.70 0.30];  % SS, CS, BS
figure(22); clf;
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for sc=1:3
    ev2show=find(EventPropsMat{foi2show}.SpClass==sc);
    if sc==2
    %ev2show([5 6])=[];
    end
    PeakAmp_spclass{sc}=[EventPropsMat{foi2show}.PostPeak_Soma(ev2show) EventPropsMat{foi2show}.PostPeak_Distal(ev2show) ev2show ones(length(ev2show),1)*sc];
    PeakAUC_spclass{sc}=[EventPropsMat{foi2show}.PostAUC_Soma(ev2show) EventPropsMat{foi2show}.PostAUC_Distal(ev2show) ev2show ones(length(ev2show),1)*sc];

    nexttile(1,[1 1])
    scatter(PeakAmp_spclass{sc}(:,1),PeakAmp_spclass{sc}(:,2),20,cmap_sc(sc,:),'filled'); hold all;
    xlabel('Peak amplitude at soma (Z-score)'); ylabel('Peak amplitude at distal dendrite (Z-score)');

    nexttile(2,[1 1])
    scatter(PeakAUC_spclass{sc}(:,1),PeakAUC_spclass{sc}(:,2),20,cmap_sc(sc,:),'filled'); hold all;
    xlabel('AUC at soma (Z-score)'); ylabel('AUC at dendrite (Z-score)');
end
legend({'SS','CS','BS'},'Location','northwest')
set_fontsize(15);

SomROI = PCresult{foi2show}.roi_dClass ==2;
DendROI = PCresult{foi2show}.roi_dClass ==5;

figure(200); clf; l=[];
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for sc = 1:3
        
        SomTraces2show=mean(AlignedEvntall{foi2show,2}(SomROI,:,PeakAmp_spclass{sc}(:,3)),1,'omitnan');
        DendTraces2show=mean(AlignedEvntall{foi2show,2}(DendROI,:,PeakAmp_spclass{sc}(:,3)),1,'omitnan');
        N_spike=length(PeakAmp_spclass{sc}(:,3));

        MSomTraces2show=squeeze(mean(SomTraces2show,3,'omitnan')); SSomTraces2show=squeeze(std(SomTraces2show,0,3,'omitnan'))./sqrt(N_spike);
        MDendTraces2show=squeeze(mean(DendTraces2show,3,'omitnan')); SDendTraces2show=squeeze(std(DendTraces2show,0,3,'omitnan'))./sqrt(N_spike);

        t2show=[nTau_EV(1):nTau_EV(2)];
        nexttile(1,[1 1]); hold on;
        errorbar_shade(t2show,MSomTraces2show,SSomTraces2show,cmap_sc(sc,:));       
        xlabel('Time from spike (ms)');
        ylabel('Voltage (Z-score)');
        xline(0, '--k', 'LineWidth', 1);                          % spike time
        xlim([-50 150]);
        title(['Somatic dendrite voltage (mean ± SEM)'])

        nexttile(2,[1 1]); hold on;
        l(sc)=errorbar_shade(t2show,MDendTraces2show,SDendTraces2show,cmap_sc(sc,:));       
        xlabel('Time from spike (ms)');
        ylabel('Voltage (Z-score)');
        xline(0, '--k', 'LineWidth', 1);                          % spike time
        xlim([-50 150]);
        title(['Distal dendrite voltage (mean ± SEM)'])
    end
legend(l,{'SS','CS','BS'})
set_fontsize(15)%%
figure(24); clf;
PeakAmpBin=[-0.5 0.5 1 5]*5;
PeakAUCBin=[-30 30 100];
clear CSProbhist
CSProbhist{1}=NaN(max(foi),length(PeakAmpBin)-1);
CSProbhist{2}=NaN(max(foi),length(PeakAUCBin)-1);
probcmap=copper(length(PeakAmpBin)-1); probcmapAUC=copper(length(PeakAUCBin)-1);
PeakAmpAUC_Cat=cell(3,1);
for f=foi(1:end-1)
    for sc=1:3
        ev2show=find(EventPropsMat{f}.SpClass==sc);
        PeakAmp_spclass{sc}=[EventPropsMat{f}.PostPeak_Soma(ev2show) EventPropsMat{f}.PostPeak_Distal(ev2show) ev2show ones(length(ev2show),1)*sc];
        PeakAUC_spclass{sc}=[EventPropsMat{f}.PostAUC_Soma(ev2show) EventPropsMat{f}.PostAUC_Distal(ev2show) ev2show ones(length(ev2show),1)*sc];
        PeakAmpAUC_Cat{sc}=[PeakAmpAUC_Cat{sc}; [EventPropsMat{f}.PostPeak_Distal(ev2show) EventPropsMat{f}.PostAUC_Distal(ev2show)]./PCresult{f}.SpikeHeight];
    end

    PeakAmpCat=cell2mat(PeakAmp_spclass');
    PeakAmpAll=PeakAmpCat(:,2);%./PCresult{f}.SpikeHeight;
    SpikeClassAll=PeakAmpCat(:,4);
    [~,~,~,bin_ind]=binning_data({repmat(PeakAmpAll,1,2)},PeakAmpBin);
    CSProbhist{1}(f,:)=cellfun(@(x) sum(x & SpikeClassAll==2)./sum(x),bin_ind);
    CSProbhist{1}(f,cellfun(@sum,bin_ind)'<11)=NaN;

    PeakAUCCat=cell2mat(PeakAUC_spclass');
    PeakAUCpAll=PeakAUCCat(:,2)./PCresult{f}.SpikeHeight;
    [~,~,~,bin_ind]=binning_data({repmat(PeakAUCpAll,1,2)},PeakAUCBin);
    CSProbhist{2}(f,:)=cellfun(@(x) sum(x & SpikeClassAll==2)./sum(x),bin_ind);
    CSProbhist{2}(f,cellfun(@sum,bin_ind)'<4)=NaN;
end
nexttile(1,[1 1])
Boxplot_wPoints2(CSProbhist{1},probcmap);
ylim([0 1.2]);
set(gca,'XTick',[1:length(PeakAmpBin)-1],'XTickLabel',edge2bin(PeakAmpBin));
nexttile(2,[1 1])
Boxplot_wPoints2(CSProbhist{2},probcmapAUC);
set(gca,'XTick',[1:length(PeakAUCBin)-1],'XTickLabel',edge2bin(PeakAUCBin));
ylim([0 1.2]);

Ampbin=[-0.5:0.1:2.5]; AUCbin=[-30:4:90];
for sc=1:3
    nexttile(3,[1 1])
    histogram(PeakAmpAUC_Cat{sc}(:,1),Ampbin,'FaceColor',cmap_sc(sc,:),'FaceAlpha',0.7,'Normalization','probability'); hold all
    box off;
    nexttile(4,[1 1])
    histogram(PeakAmpAUC_Cat{sc}(:,2),AUCbin,'FaceColor',cmap_sc(sc,:),'FaceAlpha',0.7,'Normalization','probability'); hold all
    box off;
end
nexttile(3,[1 1])
xlabel('Distal peak amplitude (somatic spike height)'); ylabel('Probability');
legend({'SS','CS','BS'});
nexttile(4,[1 1])
xlabel('Distal AUC (A.U.)'); ylabel('Probability');




% %% 
% %
% % Extracts traces from 20 to 80 ms after the triggering spike for:
% %   - Compartment 1: Peri-Soma       (roi_dClass == 2)
% %   - Compartment 2: Apical Distal   (roi_dClass == 5)
% % Event types from EventPropsMat{f,1}.SpClass: 1 = SS, 2 = CS, 3 = BS
% %
% % Parameters computed per event × compartment:
% %   PeakAmp  : mean of the top-5 time points in the window (ΔF/F)
% %   AUC      : trapz of the positive part of the raw trace in the window
% %   LenRaw   : number of frames (ms) the RAW trace is above 0 in the window
% %   LenSub   : number of frames (ms) the SUBTHRESHOLD trace is above 0 in the window
% %
% % Assumes 1 kHz sampling (1 frame = 1 ms), nTau_EV = [-200 500].
% 
% % ---------- Window definition ----------
% post_win= [20 80];              % ms after trigger
% show_win= [-50 200];              % display window (ms), wider for context
% trig_idx = -nTau_EV(1) + 1;      % 1-indexed trigger frame = 201
% win_idx      = (trig_idx + post_win(1))  : (trig_idx + post_win(2));   % 221:281 – parameter window
% show_wind_inx = (trig_idx + show_win(1)) : (trig_idx + show_win(2));   % 151:401 – display window
% 
% SpClassLabels = {'SS', 'CS', 'BS'};
% compLabels    = {'Soma', 'Apical Distal'};
% compDClass    = [2, 5];           % roi_dClass values for each compartment
% 
% foi2show=20;
% rawAll = AlignedEvntall{foi2show,2};   % [nROI x nTime x nEv]  baseline corrected ΔF/F
% subAll = AlignedEvntall{foi2show,3};   % [nROI x nTime x nEv]  subthreshold
% 
% spClass_ev = EventPropsMat{f,1}.SpClass;   % [nEv x 1], values 1/2/3
% PeakAmp_ftn=@(v) mean(topkrows(v,5),'omitnan');
% AUC_ftn=  @(v) trapz(max(v, 0));
% 
% for sc=1:3
% 
% end
% 
% 
% 
% 
% 
% % ---------- Initialize accumulators ----------
% for sc = 1:3
%     for comp = 1:2
%         SpProp(sc,comp).PeakAmp  = [];
%         SpProp(sc,comp).AUC      = [];
%         SpProp(sc,comp).LenRaw   = [];
%         SpProp(sc,comp).LenSub   = [];
%         SpProp(sc,comp).traceRaw = [];  % [nWin x nEvents] for grand-mean trace
%         SpProp(sc,comp).traceSub = [];
%     end
% end
% 
% % ---------- Aggregate across files ----------
% for f = foi
%     if isempty(EventPropsMat{f,1}), continue; end
%     spClass_ev = EventPropsMat{f,1}.SpClass;   % [nEv x 1], values 1/2/3
% 
%     rawAll = AlignedEvntall{f,2};   % [nROI x nTime x nEv]  baseline corrected ΔF/F
%     subAll = AlignedEvntall{f,3};   % [nROI x nTime x nEv]  subthreshold
% 
%     for sc = 1:3
%         evIdx = find(spClass_ev == sc);
%         if isempty(evIdx), continue; end
% 
%         for comp = 1:2
%             roiMask = PCresult{f}.roi_dClass == compDClass(comp);
%             if ~any(roiMask), continue; end
% 
%             % ---- Display-window traces (-50 to 200 ms) for plotting ----
%             % mean over ROI dim → [1 x nShow x nEvType] → squeeze → [nShow x nEvType]
%             trRaw = squeeze(mean(rawAll(roiMask, show_wind_inx, evIdx), 1, 'omitnan'));
%             trSub = squeeze(mean(subAll(roiMask, show_wind_inx, evIdx), 1, 'omitnan'));
%             if length(evIdx) == 1
%                 trRaw = trRaw(:);
%                 trSub = trSub(:);
%             end
%             if size(trRaw, 1) ~= length(show_wind_inx)
%                 trRaw = trRaw';
%                 trSub = trSub';
%             end
%             % trRaw / trSub : [nShow x nEvType]
% 
%             % ---- Parameter-window traces (20–80 ms) for scalar metrics ----
%             trRaw_p = squeeze(mean(rawAll(roiMask, win_idx, evIdx), 1, 'omitnan'));
%             trSub_p = squeeze(mean(subAll(roiMask, win_idx, evIdx), 1, 'omitnan'));
%             if length(evIdx) == 1
%                 trRaw_p = trRaw_p(:);
%                 trSub_p = trSub_p(:);
%             end
%             if size(trRaw_p, 1) ~= length(win_idx)
%                 trRaw_p = trRaw_p';
%                 trSub_p = trSub_p';
%             end
%             % trRaw_p / trSub_p : [nParamWin x nEvType]
% 
%             nEvType = length(evIdx);
%             for e = 1:nEvType
%                 v = trRaw_p(:, e);
%                 s = trSub_p(:, e);
% 
%                 % Peak amplitude: mean of top-5 time points in 20-80 ms window
%                 vSort   = sort(v, 'descend');
%                 peakAmp = mean(vSort(1:min(5, length(vSort))), 'omitnan');
% 
%                 % AUC: trapz of positive part of raw trace (20-80 ms)
%                 auc = trapz(max(v, 0));
% 
%                 % Length above 0 within 20-80 ms window, in ms (1 frame = 1 ms)
%                 lenRaw = sum(v > 0);
%                 lenSub = sum(s > 0);
% 
%                 SpProp(sc,comp).PeakAmp(end+1,1)  = peakAmp;
%                 SpProp(sc,comp).AUC(end+1,1)      = auc;
%                 SpProp(sc,comp).LenRaw(end+1,1)   = lenRaw;
%                 SpProp(sc,comp).LenSub(end+1,1)   = lenSub;
%             end
% 
%             % Store display-window traces for mean ± SEM plots
%             SpProp(sc,comp).traceRaw = [SpProp(sc,comp).traceRaw, trRaw];
%             SpProp(sc,comp).traceSub = [SpProp(sc,comp).traceSub, trSub];
%         end
%     end
% end
% 
% % ---------- Plot 1: Mean ± SEM voltage traces (-50 to 200 ms) ----------
% colors_sc = [0.20 0.45 0.80;   % SS  (blue)
%              0.85 0.20 0.20;   % CS  (red)
%              0.20 0.70 0.30];  % BS  (green)
% tplot = (show_win(1) : show_win(2))';   % -50:200 ms axis
% 
% figure(200); clf;
% tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
% for comp = 1:2
%     for sc = 1:3
%         nexttile; hold on;
%         tr = SpProp(sc,comp).traceRaw;
%         if ~isempty(tr)
%             mTr   = mean(tr, 2, 'omitnan');
%             semTr = std(tr, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(tr), 2));
%             fill([tplot; flipud(tplot)], [mTr+semTr; flipud(mTr-semTr)], ...
%                  colors_sc(sc,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%             plot(tplot, mTr, 'Color', colors_sc(sc,:), 'LineWidth', 2);
%         end
%         nEv = size(SpProp(sc,comp).traceRaw, 2);
%         title([compLabels{comp} ' – ' SpClassLabels{sc} '  (n=' num2str(nEv) ')']);
%         xlabel('Time from spike (ms)');
%         ylabel('\DeltaF/F');
%         xline(0, '--k', 'LineWidth', 1);                          % spike time
%         xline(post_win(1), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);  % param window start
%         xline(post_win(2), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);  % param window end
%         xlim(show_win);
%         box off;
%     end
% end
% sgtitle('Mean ± SEM Voltage Trace (–50 to 200 ms post-spike) by Spike Class & Compartment');
% 
% % ---------- Plot 2: Parameter comparison (scatter + mean) ----------
% paramLabels = {'Peak Amplitude (\DeltaF/F)', 'AUC (a.u.)', ...
%                'Length Raw (ms above 0)', 'Length Subthreshold (ms above 0)'};
% paramFields = {'PeakAmp', 'AUC', 'LenRaw', 'LenSub'};
% 
% figure(201); clf;
% tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% for p = 1:4
%     for comp = 1:2
%         nexttile; hold on;
%         allData  = [];
%         groupVec = [];
%         for sc = 1:3
%             dat = SpProp(sc,comp).(paramFields{p});
%             dat = dat(~isnan(dat));
%             if isempty(dat), continue; end
%             allData  = [allData;  dat];
%             groupVec = [groupVec; sc * ones(length(dat), 1)];
%             xjitter  = sc + (rand(length(dat), 1) - 0.5) * 0.35;
%             scatter(xjitter, dat, 12, colors_sc(sc,:), 'filled', 'MarkerFaceAlpha', 0.35);
%             plot([sc-0.3, sc+0.3], [mean(dat) mean(dat)], '-', ...
%                  'Color', colors_sc(sc,:), 'LineWidth', 2.5);
%         end
%         % Kruskal-Wallis across spike classes
%         if length(unique(groupVec)) > 1 && length(allData) > 2
%             [pKW, ~, ~] = kruskalwallis(allData, groupVec, 'off');
%             title([compLabels{comp} '  (KW p = ' sprintf('%.3f', pKW) ')']);
%         else
%             title(compLabels{comp});
%         end
%         set(gca, 'XTick', 1:3, 'XTickLabel', SpClassLabels);
%         ylabel(paramLabels{p});
%         box off;
%     end
% end
% sgtitle('Spike-class comparison: raw trace parameters, 20–80 ms post-spike');
% 
% % ---------- Plot 3: SS vs CS vs BS – overlaid per compartment ----------
% % Each panel is one compartment; all three spike classes overlaid (mean ± SEM)
% figure(202); clf;
% tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% for comp = 1:2
%     nexttile; hold on;
%     hLeg = gobjects(0);
%     legStr = {};
%     for sc = 1:3
%         tr = SpProp(sc,comp).traceRaw;
%         if isempty(tr), continue; end
%         mTr   = mean(tr, 2, 'omitnan');
%         semTr = std(tr, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(tr), 2));
%         fill([tplot; flipud(tplot)], [mTr+semTr; flipud(mTr-semTr)], ...
%              colors_sc(sc,:), 'FaceAlpha', 0.18, 'EdgeColor', 'none');
%         h = plot(tplot, mTr, '-', 'Color', colors_sc(sc,:), 'LineWidth', 2.5);
%         hLeg(end+1) = h;
%         legStr{end+1} = [SpClassLabels{sc} '  (n=' num2str(size(tr,2)) ')'];
%     end
%     xline(0, '--k', 'LineWidth', 1);                                       % spike time
%     xline(post_win(1), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);      % param window
%     xline(post_win(2), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
%     xlim(show_win);
%     xlabel('Time from spike (ms)');
%     ylabel('\DeltaF/F');
%     title(compLabels{comp});
%     legend(hLeg, legStr, 'Location', 'best');
%     box off;
% end
% sgtitle('SS vs CS vs BS mean trace (±SEM) — dashed lines mark parameter window (20–80 ms)');
% 
% % ---------- Summary table printed to Command Window ----------
% fprintf('\n===== Parameter Summary: 20–80 ms post-spike =====\n');
% for comp = 1:2
%     fprintf('\n--- %s ---\n', compLabels{comp});
%     fprintf('%-6s  %10s  %10s  %10s  %10s\n', ...
%             'Class', 'PeakAmp', 'AUC', 'LenRaw(ms)', 'LenSub(ms)');
%     for sc = 1:3
%         n = length(SpProp(sc,comp).PeakAmp(~isnan(SpProp(sc,comp).PeakAmp)));
%         fprintf('%-6s  %10.4f  %10.4f  %10.1f  %10.1f  (n=%d)\n', ...
%             SpClassLabels{sc}, ...
%             mean(SpProp(sc,comp).PeakAmp, 'omitnan'), ...
%             mean(SpProp(sc,comp).AUC,     'omitnan'), ...
%             mean(SpProp(sc,comp).LenRaw,  'omitnan'), ...
%             mean(SpProp(sc,comp).LenSub,  'omitnan'), n);
%     end
% end
% fprintf('\n');
% 
% % Each panel is one compartment; all three spike classes overlaid (mean ± SEM)
% figure(204); clf;
% tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
% for sc = 1:3
%     PeakSom = SpProp(sc,1).PeakAmp;
%     PeakDend = SpProp(sc,2).PeakAmp;
%     scatter(PeakSom,PeakDend,8,colors_sc(sc,:),'filled'); hold all
% end
% xlabel('Peak post-spike amplitude at soma (Z-score)');
% ylabel('Peak post-spike at distal dendrite (Z-score)');
% legend({'SS','CS','BS'});
% box off;