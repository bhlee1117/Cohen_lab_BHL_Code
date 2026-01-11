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
%foi=[1 4 5 6 8 10 11 15 16 17 18 19 20 21 22 23 24 25 26 27];
foi=22;
%% Load Result files
win_pre = 30; win_post = 20;
nTau = [-200:50]; % SS, CS, Brst
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

PerisomaRange=[-70 70]; presub_readtime=[-10:-3];
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
    end

    PreSub=squeeze(mean(AlignedEvntall{f,3}(:,-nTau_EV(1)+presub_readtime,:),2,'omitnan'));
    PreSub_Basal=squeeze(mean(AlignedEvntall{f,3}(ismember(PCresult{f}.roi_dClass,1),-nTau_EV(1)+presub_readtime,:),[1 2],'omitnan'));
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
    PreSub_Soma=cellfun(@(x) x(2,2),PreSubbin_cell(~emptycellEV));

    Evlab=[EVList' PCresult{f}.BlueStim(EVList)'>0 sum(isnan(PCresult{f}.NormalizedTrace_dirt(:,EVList)),1)'>0 PCresult{f}.PFvec(EVList)'...
        PCresult{f}.runVec(EVList)' IntervalMat' centroid_EVAUC' trnsmit_EVAUC EVAUC_Apical EVAUC_Soma ISIEV' EVAmp_Apical Evtype' NSpike' PreSub_Soma PreSub_Apical TRevMax' TRevFirst' PreSub_Distal PreSub_Basal...
        PostSub_Basal PostSub_Distal ApicalAUCFirst'];
    EventPropsMat{f,1}=[Evlab];
    EventPropsMat{f,2}={TransRatebAPinEvent,AUCrawEV,AUCEVbin_cell,AlignedEvntcell_bin};

    fieldName={'SpikeFrame','IsBlue','IsNA','IsPF','IsRun','Length','AUC_centroid','TransmissionRate','AUC_apical','AUC_soma',...
        'ISI','Amp_apical','Spike_type','Nspike','PreSub_Apical','PreSub_Soma','TRmax_bAP','TRfrst_bAP','PreSub_Distal',...
        'PreSub_Basal','PostSub_Basal','PostSub_Distal','AUCapical_frst'};
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

%% dSpikes
nTauDelay=[-3:4]; % Short time window to estimate delay.
nTauAlign=[20 15];
dclass_str={'Basal','Soma','Trunk','Oblique','Distal'};
for f=[foi]
    [nROI nTime]=size(PCresult{f}.Subthreshold);
    nROI_basal=find(ismember(PCresult{f}.branch_dClass,[1]));
    nROI_apical=find(ismember(PCresult{f}.branch_dClass,[5]));
    bAP_delay=PCresult{f}.SpikeDelay(PCresult{f}.sortDist);
    validTimeLength=sum(sum(isnan(PCresult{f}.Subthreshold),1)==0);
    branchInd=PCresult{f}.BranchLabel(PCresult{f}.sortDist);
    SortedFtprnt=PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist);
    SpikeMatbranch=NaN(max(branchInd),nTime); SpikeMatdClass=NaN(5,nTime);
    dSpikeTA_Branch=zeros(nROI,length(nTau),max(branchInd)); dSpikeTA_dClass=zeros(nROI,length(nTau),5);
    dSpikeTAsp_Branch=zeros(max(branchInd),length(nTau)); dSpikeTAsp_dClass=zeros(5,length(nTau));
    FtprntBranch=[]; FtprntdClass=[];

    SomaSpikeFrame=find(max(PCresult{f}.SpikeMat(PCresult{f}.BranchLabel==1,:),[],1)>0);
    SpikeMatsorted=PCresult{f}.SpikeMat;
    SpikeMatsorted=SpikeMatsorted(PCresult{f}.sortDist,:);
    bAPFrame=unique(SomaSpikeFrame'+[0:3]);

    SpikeVec=(max(SpikeMatsorted,[],1)>0);
    SpikeVec=SpikeVec.*(PCresult{f}.BlueStim==0);
    [~, SpikeTrV ,SpikeFrame]=get_STA(PCresult{f}.NormalizedTrace_dirt-PCresult{f}.Subthreshold_larger,SpikeVec,10,10); %align every spike
    SpikeTrV=permute(SpikeTrV,[1 3 2]);

    tmax_mat=repmat(11,nROI,1,length(SpikeFrame)); % n X t X R;
    tmax_mat(:,:,find(ismember(SpikeFrame,SomaSpikeFrame)))=repmat(bAP_delay+11,1,1,sum(ismember(SpikeFrame,SomaSpikeFrame)));

    [AUCbAP AUCrawbAP KinkbAP Kinkraw]=get_AUC(SpikeTrV,tmax_mat,6,6);
    dSpikePropsMat{f}=[];

    for sp=SpikeFrame %Spikes
        brOI=unique(branchInd(find(SpikeMatsorted(:,sp)>0))); % branch that have spike
        spInd=find(SpikeFrame==sp);
        for b=brOI %branches
            nROI_brOI=find(ismember(branchInd,b'));
            nROI_Soma=find(ismember(branchInd,1));
            DendriteClass=unique(PCresult{f}.branch_dClass(nROI_brOI));
            distfromsoma=mean(PCresult{f}.dendaxis(nROI_brOI));
            if length(DendriteClass)>1; error('There is more than one dendrite class in your ROI'); end

            nearestbackSoma=SomaSpikeFrame(find(SomaSpikeFrame>=sp,1));
            nearestFrontSoma=find(SomaSpikeFrame<=sp);
            if isempty(find(SomaSpikeFrame>sp,1)) | isempty(nearestFrontSoma) %Spike at the first and end of trace
                nearestbackSoma=NaN; nearestFrontSoma=NaN;
                props2cat=[sp b DendriteClass distfromsoma ismember(sp,bAPFrame) PCresult{f}.CStrace(sp) NaN(1,4) ...
                    mean(AUCrawbAP(nROI_brOI,spInd),1,'omitnan') mean(AUCrawbAP(nROI_Soma,spInd),1,'omitnan') mean(AUCrawbAP(nROI_apical,spInd),1,'omitnan') mean(AUCrawbAP(nROI_basal,spInd),1,'omitnan') PCresult{f}.SpikeOrder(sp) NaN];
            else
                nearestFrontSoma=SomaSpikeFrame(nearestFrontSoma(end));
                if ismember(sp,bAPFrame)
                    nearestbackSoma=SomaSpikeFrame(find(SomaSpikeFrame>sp,1));
                end

                VoltageTrace_periSpike=[mean(PCresult{f}.NormalizedTrace_dirt(nROI_Soma,nearestFrontSoma + nTauDelay),1,'omitnan'); ...
                    mean(PCresult{f}.NormalizedTrace_dirt(nROI_Soma,nearestbackSoma + nTauDelay),1,'omitnan'); ...
                    mean(PCresult{f}.NormalizedTrace_dirt(nROI_brOI,sp + nTauDelay),1,'omitnan')]; %Soma_front, Soma_back, dendrite
                SpikeTimeInterp = get_delay(VoltageTrace_periSpike, 20, -nTauDelay(1)+[-1 2]) + [nearestFrontSoma nearestbackSoma sp]' + nTauDelay(1)-1;
                dt=[sp-nearestFrontSoma, nearestbackSoma-sp];
                [~, shift]=min(abs(dt));
                if shift==1;
                    spord=PCresult{f}.SpikeOrder(nearestFrontSoma);
                    nearestSomSP=nearestFrontSoma;
                else
                    spord=PCresult{f}.SpikeOrder(nearestbackSoma);
                    nearestSomSP=nearestbackSoma;
                end
                dt_interp=[SpikeTimeInterp(3)-SpikeTimeInterp(1) SpikeTimeInterp(2)-SpikeTimeInterp(3)];
                props2cat=[sp b DendriteClass distfromsoma ismember(sp,bAPFrame) PCresult{f}.CStrace(sp) dt dt_interp ...
                    mean(AUCrawbAP(nROI_brOI,spInd),1,'omitnan') mean(AUCrawbAP(nROI_Soma,spInd),1,'omitnan') mean(AUCrawbAP(nROI_apical,spInd),1,'omitnan') mean(AUCrawbAP(nROI_basal,spInd),1,'omitnan') ...
                    spord nearestSomSP];
                % Spike frame, brach Index, dClass, is bAP, distance from soma, time from preceding soma
                % spike, time from following spike, interpolated time pre,
                % post, AUC at the branch, AUC at the soma, AUC at Basal, Spike
                % Order, nearest somAP
            end
            dSpikePropsMat{f}=[dSpikePropsMat{f}; props2cat];
        end
    end
    fieldName={'Spike_frame','Branch_Index','dClass','Distance_from_soma','IsbAP','IsCS','dt_pre_somAP','dt_post_somAP'...
        'dt_pre_somAP_int','dt_post_somAP_int','AUC_branch','AUC_soma','AUC_apical','AUC_basal','SpikeOrder','Nearest_somAP'};
    dSpikePropsMat{f}=array2table(dSpikePropsMat{f},'VariableNames',fieldName);

    % Define somAP and dAP
    conductionspeed=230; %um/ms
    dt_error=1; %ms
    bAPthreshold=@(x) conductionspeed*(x+dt_error);
    dAPCandidates=find(~dSpikePropsMat{f}.IsbAP & ~dSpikePropsMat{f}.IsCS & dSpikePropsMat{f}.dClass~=2 & dSpikePropsMat{f}.dt_post_somAP<2.5);
    bAPsCandidates=find(dSpikePropsMat{f}.IsbAP & ~dSpikePropsMat{f}.IsCS & dSpikePropsMat{f}.dClass~=2 & dSpikePropsMat{f}.dt_pre_somAP<4);
    catCandidates=[bAPsCandidates;dAPCandidates];
    dt2show=[dSpikePropsMat{f}.dt_pre_somAP_int(bAPsCandidates); -dSpikePropsMat{f}.dt_post_somAP_int(dAPCandidates)];
    dst2show=[dSpikePropsMat{f}.Distance_from_soma([bAPsCandidates;dAPCandidates])];
    dAPCandidates_threshold=catCandidates(find((dst2show-bAPthreshold(dt2show))>0));

    SomAPInd=dSpikePropsMat{f}.Branch_Index==1;
    CSframeInd=dSpikePropsMat{f}.Branch_Index==1 & dSpikePropsMat{f}.dClass~=0 & dSpikePropsMat{f}.IsCS==1 & dSpikePropsMat{f}.SpikeOrder==1;
    SSframeInd=dSpikePropsMat{f}.Branch_Index==1 & dSpikePropsMat{f}.dClass~=0 & dSpikePropsMat{f}.IsCS==0;

    FiringRate{f}.DendIntSpFrame=unique(dSpikePropsMat{f}.Nearest_somAP(dAPCandidates_threshold));
    FiringRate{f}.SomaSpikeFrame=setdiff(dSpikePropsMat{f}.Spike_frame(SomAPInd),FiringRate{f}.DendIntSpFrame);
    FiringRate{f}.ComplexSpikeFrame=dSpikePropsMat{f}.Spike_frame(CSframeInd);
    FiringRate{f}.SimpleSpikeFrame=dSpikePropsMat{f}.Spike_frame(SSframeInd);
    FiringRate{f}.ValidPeriod=validTimeLength;

    DendOnlySpikes=find(dSpikePropsMat{f}.Branch_Index~=1 & dSpikePropsMat{f}.dClass~=0 & ~dSpikePropsMat{f}.IsbAP & ~dSpikePropsMat{f}.IsCS & dSpikePropsMat{f}.dt_pre_somAP>=4 & dSpikePropsMat{f}.dt_post_somAP>=2.5);
    dSpikeFrame=unique(dSpikePropsMat{f}.Spike_frame(DendOnlySpikes));
    dSpikeFrame(find(diff(dSpikeFrame)==1)+1)=[];
    FiringRate{f}.DendOnlySpFrame=setdiff(dSpikeFrame,FiringRate{f}.DendIntSpFrame);

    DendIntSpClassMat=NaN(5,length(FiringRate{f}.DendIntSpFrame));
    for d=1:length(FiringRate{f}.DendIntSpFrame)
        dint_ind=dAPCandidates_threshold(find(dSpikePropsMat{f}.Nearest_somAP(dAPCandidates_threshold)==FiringRate{f}.DendIntSpFrame(d)));
        DendIntSpClassMat(:,d)=ind2vec(5,dSpikePropsMat{f}.dClass(dint_ind)',1);
    end

    dSpikeClassMat=NaN(5,length(FiringRate{f}.DendOnlySpFrame));
    for d=1:length(FiringRate{f}.DendOnlySpFrame)
        dSpikeClassMat(:,d)=ind2vec(5,dSpikePropsMat{f}.dClass(dSpikePropsMat{f}.Spike_frame==FiringRate{f}.DendOnlySpFrame(d) & dSpikePropsMat{f}.Branch_Index~=1)',1);
    end
    dSpikeClassMat(1,max(dSpikeClassMat([3 4 5],:),[],1)>0 & max(dSpikeClassMat(1,:),[],1)>0)=0;

    dSpikeBranch=NaN(max(branchInd),2);
    for b=1:max(branchInd)
        nROI_brOI=find(ismember(branchInd,b'));
        DendriteClass=unique(PCresult{f}.branch_dClass(nROI_brOI));
        if DendriteClass~=0
            dSpikeBranch(b,:)=[sum(dSpikePropsMat{f}.Branch_Index==b)/validTimeLength*1000 DendriteClass];
        end
    end

    FiringRate{f}.SomAP=[length(FiringRate{f}.SomaSpikeFrame)]/validTimeLength*1000;
    FiringRate{f}.CS=[length(FiringRate{f}.ComplexSpikeFrame)]/validTimeLength*1000;
    FiringRate{f}.SS=[length(FiringRate{f}.SimpleSpikeFrame)]/validTimeLength*1000;
    FiringRate{f}.dAP= length(FiringRate{f}.DendIntSpFrame)/validTimeLength*1000;
    FiringRate{f}.dAPdClass=[sum(DendIntSpClassMat(1,:)>0),...
        sum(max(DendIntSpClassMat([3 4 5],:),[],1)>0)]/size(DendIntSpClassMat,2);
    FiringRate{f}.dSpike=[length(FiringRate{f}.DendOnlySpFrame)]/validTimeLength*1000;
    FiringRate{f}.dSpdClass=[sum(dSpikeClassMat(1,:)>0),...
        sum(max(dSpikeClassMat([3 4 5],:),[],1)>0)]/size(dSpikeClassMat,2); %basal, Apical
    FiringRate{f}.dSpikeBasal=FiringRate{f}.DendOnlySpFrame(dSpikeClassMat(1,:)>0);
    FiringRate{f}.dSpikeApical=FiringRate{f}.DendOnlySpFrame(max(dSpikeClassMat([3 4 5],:),[],1)>0);

    [~, AligndSpiketr{f}]=get_STA(PCresult{f}.NormalizedTrace_dirt,dSpikePropsMat{f}.Spike_frame,nTauAlign(1),nTauAlign(2));
    AligndSpiketr{f}=permute(AligndSpiketr{f},[1 3 2]);
end

%% Show dAP neuron image (Figure S_dAP)
f=22;
conductionspeed=230; %um/ms
dt_error=1; %ms
bAPthreshold=@(x) conductionspeed*(x+dt_error);
dAPCandidates=find(~dSpikePropsMat{f}.IsbAP & ~dSpikePropsMat{f}.IsCS & dSpikePropsMat{f}.dClass~=2 & dSpikePropsMat{f}.dt_post_somAP<2.5);
bAPsCandidates=find(dSpikePropsMat{f}.IsbAP & ~dSpikePropsMat{f}.IsCS & dSpikePropsMat{f}.dClass~=2 & dSpikePropsMat{f}.dt_pre_somAP<4);
catCandidates=[bAPsCandidates;dAPCandidates];
dt2show=[dSpikePropsMat{f}.dt_pre_somAP_int(bAPsCandidates); -dSpikePropsMat{f}.dt_post_somAP_int(dAPCandidates)];
dst2show=[dSpikePropsMat{f}.Distance_from_soma([bAPsCandidates;dAPCandidates])];
dAPCandidates_threshold=catCandidates(find((dst2show-bAPthreshold(dt2show))>0));

dAPSTA=get_STA(PCresult{f}.NormalizedTrace_dirt,FiringRate{f}.DendIntSpFrame,10,20);
somAPSTA=get_STA(PCresult{f}.NormalizedTrace_dirt,setdiff(FiringRate{f}.SimpleSpikeFrame,FiringRate{f}.DendIntSpFrame),10,20);
nROI=size(PCresult{f}.NormalizedTrace_dirt,1);
dAPSTA=dAPSTA-median(dAPSTA(:,1:5),2);
somAPSTA=somAPSTA-median(somAPSTA(:,1:5),2);

figure(33); clf; tiledlayout(3,2);
ax1=nexttile([2 1]); cax=[0 10];
imshow2(PCresult{f}.avgImg,[]);
drawScaleBar(100/PCresult{f}.pixelsize,'horizontal','color',[1 1 1]);
nexttile([2 1]);
scatter(dt2show,dst2show,20,[0 0 0],'filled'); hold all
f2=plot([-1:2],bAPthreshold([-1:2]),'r--','linewidth',2);
xlabel('\Deltat (ms)'); ylabel('Distance from soma (\mum)'); axis tight;
ax3=nexttile([1 1]);
imagesc([-10:20],[1:nROI],dAPSTA,cax);
set_kymoYtick(PCresult{f}.dendaxis)
xlabel('Time (ms)'); ylabel('Distance (\mum)'); axis tight;
title('Average of spikes inititated from dendrite (SomAPs)')
ax2=nexttile([1 1]);
imagesc([-10:20],[1:nROI],somAPSTA,cax);
set_kymoYtick(PCresult{f}.dendaxis)
title('Average of spikes inititated from soma (SomAPs)')
xlabel('Time (ms)'); ylabel('Distance (\mum)'); axis tight;
colormap(ax2,turbo);
colormap(ax3,turbo);
cb=colorbar;
cb.Label.String='Z score';
set_fontsize(13);

%% Save STA movies
time_segment=15000; bound=5; overlap=200; nTauSTA=[150 150]; Dirtthreshold=0.3;
for f=[foi]
    f
    cd(fpath{f})
    SpikeLabelMat=[]; MovSTAinfo=[];
    [nROI nTime]=size(PCresult{f}.Subthreshold);
    load(fullfile(fpath{f},'PC_Result.mat'),'Result');
    load([fpath{f} '/output_data.mat'])
    load([fpath{f} '/MovdAPSTAinfo.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    frm_end=EndFrame(f);

    f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
    f_seg_real=[f_seg(1:end-1)' f_seg(2:end)'];
    f_seg_real(1:end-1,2)=f_seg_real(1:end-1,2)+overlap;
    f_seg_real(2:end,1)=f_seg_real(2:end,1)+overlap+1;

    take_window=repmat([1 time_segment],length(f_seg)-1,1);
    take_window(2:end,1)=take_window(2:end,1)+overlap; take_window(1:end-1,2)=take_window(1:end-1,2)+overlap;
    take_window(end)=mod(f_seg(end),time_segment);
    take_window(take_window==0)=time_segment;

    perispike_time=unique(find(PCresult{f}.SpikeMat(1,:))'+[-5:30]); perispike_time(perispike_time<1 | perispike_time>nTime)=[];
    periblue_time=unique(find(PCresult{f}.BlueStim>0)'+[-5:30]); periblue_time(periblue_time<1 | periblue_time>nTime)=[];
    t_fit= (ind2vec(size(Result.traces_bvMask,2),periblue_time,1)==0) & (ind2vec(size(Result.traces_bvMask,2),perispike_time,1)==0);
    [bleaching_fit] = expfitDM_2(find(t_fit)',-mean(Result.traces_bvMask(:,t_fit))',[1:size(Result.traces_bvMask,2)]',[100000 10000]);

    isset_template=1;
    SpikeHeight=Result.SpikeHeight_fit;
    Mov_SpikeTA=zeros(sz(2)*sz(1),sum(nTauSTA)+1,1);
    Movfilt_SpikeTA=zeros(sz(2)*sz(1),sum(nTauSTA)+1,1);
    N_added=zeros(sz(2)*sz(1),sum(nTauSTA)+1,1);
    SpClassVec=ind2vec(nTime,FiringRate{f}.DendIntSpFrame,1);

    for j=1:length(f_seg)-1

        [nInd fInd]=find(SpClassVec(:,[f_seg_real(j):f_seg_real(j+1)]));
        [fInd_sp]=find(Result.spike(1,[f_seg_real(j):f_seg_real(j+1)]));

        if ~isempty(nInd)
            mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
            load([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat']);

            mov_mc=mov_mc(:,:,[take_window(j,1):take_window(j,2)]);
            %mc= movmean(mcTrace.xymean-movmedian(mcTrace.xymean,500,1),3,1);
            mc=mcTrace.xymean([take_window(j,1):take_window(j,2)],:);

            bkg = zeros(1, size(mov_mc,3));
            bkg(1,:) = bleaching_fit(f_seg_real(j,1):f_seg_real(j,2));  % bleaching

            mov_res= mov_mc-median(mov_mc,3);
            mov_res = SeeResiduals(mov_res,mc);
            mov_res = SeeResiduals(mov_res,mc.^2);
            mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
            mov_res = SeeResiduals(mov_res,bkg,1);
            %mov_res=tovec(mov_res);
            mov_res= mov_res./reshape(SpikeHeight(f_seg_real(j,1):f_seg_real(j,2)),1,1,[]);

            if ~isfield(MovSTAinfo,'template_basis') & isset_template % set templates
                mov_res_sub=toimg(get_subthreshold(tovec(mov_res),ind2vec(size(mov_res,3),fInd_sp,1),11,25),sz(2),sz(1));
                mov_res_sub=maskEdge(mov_res_sub,5,0);
                mov_STA_tmp=get_STA(tovec(mov_res),fInd_sp,3,3);
                mov_STA_tmp=toimg(mov_STA_tmp,sz(2),sz(1));
                %[~,~,Vsub]=get_eigvector(tovec(imgaussfilt(mov_res_sub,1.5))',10);
                %[~,~,Vsub]=get_eigvector(tovec(mov_res_sub)',20);
                [~,~,Vsub]=get_eigvector(tovec(imgaussfilt(mov_res_sub(:,:,1050:7200),1))',20);
                imshow2_patch(toimg(Vsub,sz(2),sz(1)))
                sub2keep=input('subthreshold component to keep \n');
                [~,~,Vspike]=get_eigvector(tovec(mov_STA_tmp)',7);
                imshow2_patch(toimg(Vspike,sz(2),sz(1)))
                sp2keep=input('spike component to keep \n');
                template_basis=cat(3,toimg(Vsub(:,sub2keep),sz(2),sz(1)),toimg(Vspike(:,sp2keep),sz(2),sz(1)));
                isset_template=0;
            else
                template_basis=MovSTAinfo.template_basis;
                isset_template=0;
            end
            mov_res_filt=pcafilt_template(imresize(mov_res,0.6),imresize(template_basis,0.6));
            mov_res_filt=imresize(mov_res_filt,[sz(2) sz(1)]);

            mov_res=tovec(mov_res);
            mov_res_filt=tovec(mov_res_filt);

            for k=1:length(fInd) %frames
                f2add=fInd(k);
                n2add=nInd(k);
                if f2add-nTauSTA(1)>1 & f2add+nTauSTA(2)<size(mov_mc,3)
                    mov2add=mov_res(:,f2add+[-nTauSTA(1):nTauSTA(2)]);
                    movfilt2add=mov_res_filt(:,f2add+[-nTauSTA(1):nTauSTA(2)]);
                    if isfield(Result,'dirtTrace')
                        validFrm=sum(Result.dirtTrace(:,f2add+[-nTauSTA(1):nTauSTA(2)])>0,1)/size(Result.dirtTrace,1)<Dirtthreshold;
                    else
                        validFrm=ones(1,size(mov2add,2));
                    end
                    mov2add(:,setdiff([1:sum(nTauSTA)+1],find(validFrm)))=NaN;
                    movfilt2add(:,setdiff([1:sum(nTauSTA)+1],find(validFrm)))=NaN;
                    Mov_SpikeTA(:,:,n2add)=sum(cat(3,Mov_SpikeTA(:,:,n2add),mov2add),3,'omitnan');
                    Movfilt_SpikeTA(:,:,n2add)=sum(cat(3,Movfilt_SpikeTA(:,:,n2add),movfilt2add),3,'omitnan');
                    N_added(:,:,n2add)=N_added(:,:,n2add)+(~isnan(mov2add));
                    SpikeLabelMat=[SpikeLabelMat; [f2add n2add]];
                end
            end
            disp([num2str(j) ' th segment is stacked']);
        else
            disp([num2str(j) ' has no valid index']);
        end
    end
    Mov_SpikeTA=-Mov_SpikeTA./N_added;
    Range_MovSTA=[min(Mov_SpikeTA(:)) max(Mov_SpikeTA(:))];

    Mov_SpikeTA=vm(mat2gray(double(Mov_SpikeTA))*10000);
    MovTAname=['/dAPSpikeTriggeredMatMovie_Type.bin'];
    Mov_SpikeTA.transpose.savebin([fpath{f} MovTAname])

    Movfilt_SpikeTA=-Movfilt_SpikeTA./N_added;
    Range_MovFiltSTA=[min(Movfilt_SpikeTA(:)) max(Movfilt_SpikeTA(:))];
    Movfilt_SpikeTA=vm(mat2gray(double(Movfilt_SpikeTA))*10000);
    MovTAfiltname=['/dAPSpikeTriggeredMatMovieFilt_Type.bin'];
    Movfilt_SpikeTA.transpose.savebin([fpath{f} MovTAfiltname])

    fieldName={'StackedFrame','StackedStype'};
    SpikeLabelMat=array2table(SpikeLabelMat,'VariableNames',fieldName);
    MovSTAinfo.StackedSpike=SpikeLabelMat;
    MovSTAinfo.Range_MovFiltSTA=Range_MovFiltSTA;
    MovSTAinfo.Range_MovSTA=Range_MovSTA;
    MovSTAinfo.N_added=N_added;
    MovSTAinfo.template_basis=template_basis;
    save([fpath{f} '/MovdAPSTAinfo.mat'],'MovSTAinfo','-v7.3')
end
%% Convert to dF/F and color
for f=[foi]
    cd(fpath{f}); F0img=[];
    load(fullfile(fpath{f},'PC_Result.mat'),'Result');
    load([fpath{f} '/output_data.mat'])
    load([fpath{f} '/MovdAPSTAinfo.mat'])
    %load(fullfile(fpath{f},'coloredSTA.mat'));

    [nROI nTime]=size(Result.normTraces);

    if isempty(F0img)
        sz=double(Device_Data{1, 3}.ROI([2 4]));
        frm_end=EndFrame(f);
        f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
        f_seg_real=[f_seg(1:end-1)' f_seg(2:end)'];
        f_seg_real(1:end-1,2)=f_seg_real(1:end-1,2)+overlap;
        f_seg_real(2:end,1)=f_seg_real(2:end,1)+overlap+1;

        take_window=repmat([1 time_segment],length(f_seg)-1,1);
        take_window(2:end,1)=take_window(2:end,1)+overlap; take_window(1:end-1,2)=take_window(1:end-1,2)+overlap;
        take_window(end)=mod(f_seg(end),time_segment);
        take_window(take_window==0)=time_segment;
        perispike_time=unique(find(PCresult{f}.SpikeMat(1,:))'+[-5:30]); perispike_time(perispike_time<1 | perispike_time>nTime)=[];
        periblue_time=unique(find(PCresult{f}.BlueStim>0)'+[-5:30]); periblue_time(periblue_time<1 | periblue_time>nTime)=[];
        t_fit= (ind2vec(size(Result.traces_bvMask,2),periblue_time,1)==0) & (ind2vec(size(Result.traces_bvMask,2),perispike_time,1)==0);
        [bleaching_fit] = expfitDM_2(find(t_fit)',-mean(Result.traces_bvMask(:,t_fit))',[1:size(Result.traces_bvMask,2)]',[100000 10000]);
        SpikeHeight=Result.SpikeHeight_fit;

        MovSTA=double(readBinMov('dAPSpikeTriggeredMatMovie_Type.bin',sz(2)*sz(1),301));
        MovSTA=reshape(MovSTA,sz(2),sz(1),301,[]);
        MovSTAfilt=double(readBinMov('dAPSpikeTriggeredMatMovieFilt_Type.bin',sz(2)*sz(1),301));
        MovSTAfilt=reshape(MovSTAfilt,sz(2),sz(1),301,[]);

        oldRange = [min(MovSTA(:)) max(MovSTA(:))];
        NewRange = [min(MovSTAinfo.Range_MovSTA(:)) max(MovSTAinfo.Range_MovSTA(:))];
        MovSTA = (MovSTA - oldRange(1)) * diff(NewRange) / diff(oldRange) + NewRange(1);

        oldRange = [min(MovSTAfilt(:)) max(MovSTAfilt(:))];
        NewRange = [min(MovSTAinfo.Range_MovFiltSTA(:)) max(MovSTAinfo.Range_MovFiltSTA(:))];
        MovSTAfilt = (MovSTAfilt - oldRange(1)) * diff(NewRange) / diff(oldRange) + NewRange(1);

        j=3;
        mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
        load([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat']);
        mov_mc=mov_mc(:,:,[take_window(j,1):take_window(j,2)]);
        %mc= movmean(mcTrace.xymean-movmedian(mcTrace.xymean,500,1),3,1);
        mc=mcTrace.xymean([take_window(j,1):take_window(j,2)],:);

        bkg = zeros(1, size(mov_mc,3));
        bkg(1,:) = bleaching_fit(f_seg_real(j,1):f_seg_real(j,2));  % bleaching

        mov_res= mov_mc-median(mov_mc,3);
        mov_res = SeeResiduals(mov_res,mc);
        mov_res = SeeResiduals(mov_res,mc.^2);
        mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
        mov_res = SeeResiduals(mov_res,bkg,1);
        %mov_res=tovec(mov_res);
        mov_res= mov_res./reshape(SpikeHeight(f_seg_real(j,1):f_seg_real(j,2)),1,1,[]);
        mov_res_sub=movmedian(mov_res,20,3);
        mov_res_sub=maskEdge(mov_res_sub,5,0);
        F0img=get_F0img_PCA(imresize(mov_res_sub,0.7),[3000:8000]);
        %F0img2=get_F0img(imresize(mov_res_sub,0.7),[3000:8000]);
        F0img=imresize(F0img,[sz(2) sz(1)]);
        %F0img2=imresize(F0img2,[sz(2) sz(1)]);
    end

    cellMask=point2img(Result.SWC(:,[1 2]),5,size(Result.ref_im));

    avgImg=Result.ref_im-100;%./imgaussfilt(Result.ref_im,50);
    strImgMasked=avgImg;
    strImgMasked(cellMask==1 | max(Result.ftprnt>0.5,[],3)>0)=NaN;
    strImgbkg = medfilt2nan(strImgMasked, ones(1,2)*30);
    [~, idx]=bwdist(~isnan(strImgbkg));
    strImgbkg(isnan(strImgbkg)) = strImgbkg(idx(isnan(strImgbkg)));  % assign nearest values
    strImgbkg = imgaussfiltnan(strImgbkg, 5);

    figure(f); clf; ax1=[];
    ax1=[ax1 nexttile([1 1])];  imshow2(Result.ref_im,[])
    ax1=[ax1 nexttile([1 1])];  imshow2(strImgbkg,[])
    ax1=[ax1 nexttile([1 1])];  imshow2(avgImg./strImgbkg,[])
    avgImg=avgImg./strImgbkg;
    DR=[prctile(avgImg(cellMask==0),70),prctile(avgImg(cellMask==1),95)];
    strImg_bin2=grs2rgb(avgImg,colormap('gray'),DR(1),DR(2));
    strImg_bin2=strImg_bin2(:,:,1); strImg_bin2(strImg_bin2<0.2)=0;
    ax1=[ax1 nexttile([1 1])];
    imshow2(strImg_bin2>0,[])
    linkaxes(ax1)

    if isfield(Result,'Structure')
        strImg=Result.Structure.*point2img(Result.SWC(:,[1 2]),7,size(Result.ref_im));
        strImg_bin=strImg>0;
    else
        strImg=avgImg;
        strImg_bin=strImg_bin2;
    end

    MovSTAfilt=MovSTAfilt-prctile(MovSTAfilt(:,:,1:nTauSTA(1)/2,:),60,3);
    MovSTA=MovSTA-prctile(MovSTA(:,:,1:nTauSTA(1)/2,:),60,3);

    save(fullfile(fpath{f},'dAPSTA.mat'),'MovSTAfilt','MovSTA','F0img','F0img2','-v7.3');
end
%%
f=22;
load(fullfile(fpath{f},'PC_Result.mat'),'Result');
load(fullfile(fpath{f},'dAPSTA.mat'));

STAmoviedF=MovSTA(:,:,nTauSTA(1)+[-20:20]);
STAmoviedF=STAmoviedF-median(MovSTA(:,:,1:50),3);
mask=max(strImg_bin,[],3)>0.01;
cellMask=point2img(Result.SWC(:,[1 2]),5,size(Result.ref_im));

[dtimg] = get_dtimg(STAmoviedF);
STAmoviedF=maskEdge(STAmoviedF,2,NaN);
AmpImg= max(imgaussfilt_NaN(STAmoviedF,1.5),[],3)./F0img.*mask;
StrImg=grs2rgb(strImg,colormap('gray'),1,5);
StrImg=rgb2gray(StrImg);
dtimg2=dtimg;
dtimg2(isnan(dtimg2))=median(dtimg2(:),'omitnan');
dtimg2(cellMask==0)=NaN;
dtimg2=medfilt2nan(dtimg2,[5 5]);
%%
bound=4;
% === Build SNAPT Gaussian-flash movie ===
subframeT = 0.01; % ms per subframe
initialT  = -3;   % ms
finalT    =  1;   % ms
sigmaT    = 0.1; % ms flash width

[ysize, xsize, ~] = size(STAmoviedF);
times   = initialT:subframeT:finalT;
nFrames = numel(times);

stdimgNorm = mat2gray(StrImg);
AmpNorm=grs2rgb(AmpImg,colormap('gray'),0,20);
AmpNorm=imgaussfilt(rgb2gray(AmpNorm),2);

GaussPeaksmov = zeros(ysize,xsize,nFrames);
for q = 1:nFrames
    GaussPeaksmov(:,:,q) = exp(-(dtimg2-times(q)).^2/(2*sigmaT^2));
end
%superlocmov = GaussPeaksmov .* repmat(stdimgNorm,[1 1 nFrames]) .* AmpNorm;
superlocmov = GaussPeaksmov.* AmpNorm;
superlocmov(isnan(superlocmov)) = 0;

% === Convert to RGB overlay on structural image ===
superlocColormov = zeros(ysize,xsize,3,nFrames);
StrImgGray = grs2rgb(double(stdimgNorm),colormap(gray));
for j = 1:nFrames
    ColorLayer = grs2rgb(double(superlocmov(:,:,j)), colormap("hot"),0, 0.8).*strImg_bin2;
    superlocColormov(:,:,:,j) = StrImgGray + ColorLayer*2;
end

superlocColormov=superlocColormov(bound:end-bound,bound:end-50,:,:);

figure(20); clf;
fpath2read=fpath{f};
v = VideoWriter([fpath2read '/dAPSNAPT_movie'],'MPEG-4');
%v = VideoWriter([fpath2read '/SNAPT_movie'],'Uncompressed AVI');
v.FrameRate = 25;  %can adjust this, 5 - 10 works well for me
v.Quality= 100;
open(v);

for j = 1:length(times)
    clf;
    %set(gca,'units','pixels','position',[200 0 1000 800])
    imshow2(superlocColormov(:,:,:,j),[0 1])
    pbaspect([size(double(superlocColormov(:,:,:,j)),2) size(double(superlocColormov(:,:,:,j)),1) 1]),colormap(gray)
    hold all
    axis off
    text(505,10,[num2str(times(j)+1.8) ' ms'], 'FontSize', 15, 'color', [0.99 0.99 0.99])% the value 1. is to adjust timing by eyes
    drawScaleBar(100/1.17,'horizontal','color',[1 1 1],'Linewidth',3,'Position',[120 160]);
    text(20,145,['100 \mum'], 'FontSize', 15, 'color', [0.99 0.99 0.99])% the value 1. is to adjust timing by eyes
    pause(0.1)
    set(gcf,'color','w')    % Sets background to white
    frame = getframe(gcf);
    writeVideo(v,frame);
    pause(0.1);

end;
close(v);

%%
figure(22); clf;
for j = 1:36:length(times)
    nexttile([1 1])
    %set(gca,'units','pixels','position',[200 0 1000 800])
    imshow2(superlocColormov(:,:,:,j),[0 1])
    pbaspect([size(double(superlocColormov(:,:,:,j)),2) size(double(superlocColormov(:,:,:,j)),1) 1]),colormap(gray)
    hold all
    axis off
    text(570,20,[num2str(times(j)+1.8) ' ms'], 'FontSize', 13, 'color', [0.99 0.99 0.99],'HorizontalAlignment','right')% the value 1. is to adjust timing by eyes
    set(gcf,'color','w')    % Sets background to white
end;

drawScaleBar(100/0.936,'horizontal','color',[1 1 1],'Linewidth',3,'Position',[120 160]);
    text(20,135,['100 \mum'], 'FontSize', 13, 'color', [0.99 0.99 0.99])% the value 1. is to adjust timing by eyes