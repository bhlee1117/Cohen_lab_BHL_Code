%% Load file path
clear
clc;
fprintf('Section: Load file paths and metadata\n');
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
save_figto='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism/';
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
%% Load Result files11
% All voltage traces are normalized to dF/F with the robust dF/F fit scale
% factor (dFFrobust.ScaleFactor), matching Analysis_20260522_Figure2all.m.
% Everything downstream reads PCresult{f}, so a single normalization choice
% here propagates to every Figure 3/4 panel.
fprintf('Section: Load Result files and build PCresult\n');
win_pre = 30; win_post = 20;
nTau = [-200:50]; % SS, CS, Brst
clear PCresult

for f = [foi]
    f
    load(fullfile(fpath{f}, 'PC_Result.mat'), 'Result')
    dFFrobust=importdata(fullfile(fpath{f},'RobustdFFfit.mat'));

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

    PCresult{f}.Blue=Result.Blue;
    PCresult{f}.traces_bvMask=Result.traces_bvMask;
    PCresult{f}.SpikeHeight_fit=Result.SpikeHeight_fit;
    PCresult{f}.ref_im=Result.ref_im;
    PCresult{f}.bvMask=Result.bvMask;
    PCresult{f}.normTraces=Result.normTraces;
    PCresult{f}.motionReject=Result.motionReject;
    PCresult{f}.F0_PCA=Result.F0_PCA;

    PCresult{f}.sortDist=PCresult{f}.Dist_order(PCresult{f}.noi_dist);
    dendaxis = PCresult{f}.interDendDist(1,:) .* Dsign;
    PCresult{f}.dendaxis_unsorted = dendaxis;
    PCresult{f}.R2= vertcat(dFFrobust.dFFresultsInteractive.Rsq);
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

    %NormalizedTrace = Result.normTraces ./ Result.F0_PCA;
    NormalizedTrace = Result.normTraces ./ dFFrobust.ScaleFactor;
    [bAP_STA MatbAP]= get_STA(NormalizedTrace, Result.spike(1,:), win_pre, win_post);
    bAP_STA = bAP_STA - prctile(bAP_STA, 30, 2);
    PCresult{f}.MatbAP=permute(MatbAP(PCresult{f}.sortDist,:,:),[1 3 2]);
    [PCresult{f}.SpikeHeight] = max(mean(bAP_STA(perisomaROI,:), 1), [], 2);
    [~, PCresult{f}.SpikeDelay] = max(bAP_STA, [], 2);
    PCresult{f}.SpikeDelay = PCresult{f}.SpikeDelay - (win_pre + 1);

    NormalizedTrace_dirt = NormalizedTrace;
    NormalizedTrace_dirt(:, Result.motionReject > 0) = NaN;
    NormTrCheck = cellfun(@(x) x ./ dFFrobust.ScaleFactor, Result.norm_trace_check, 'UniformOutput', false);
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
    PCresult{f}.normTraces=Result.normTraces;
    PCresult{f}.scalingfactor=dFFrobust.ScaleFactor;

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

    % dendrite-class label per sorted ROI (1 basal, 2 soma, 3 trunk, 4 oblique, 5 distal)
    roisD_order_ind = cellfun(@find, PCresult{f}.roisD_order, 'UniformOutput', false);
    PCresult{f}.labelvec = NaN(1, nROI);
    for dClass = 1:5
        PCresult{f}.labelvec(cell2mat(roisD_order_ind(dClass,:)')) = dClass;
    end

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
        binTrack = ceil(PCresult{f}.VR(5,:) / (115 / 150)); %bin # = 150, VR length = 115
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
%% bAP gallery
fprintf('Section: bAP gallery\n');
for f=foi
[STAtr]= get_STA(PCresult{f}.normTraces, PCresult{f}.SpikeMat(1,:),50,30);
STAtr=STAtr-median(STAtr(:,1:40),2);
STAtr=STAtr./PCresult{f}.scalingfactor;
figure(f); clf; tiledlayout(1,2,'Padding','compact')
nexttile([1 1]);
[~, dsort]=sort(PCresult{f}.dendaxis_unsorted);
imagesc(STAtr(dsort,:)); colormap(turbo(256));
set(gca,'ytick',[1:size(STAtr,1)],'YTickLabel',num2str(dsort'))
nexttile([1 1])
scatter(PCresult{f}.dendaxis_unsorted,max(STAtr,[],2),30,'filled'); hold all;
text(PCresult{f}.dendaxis_unsorted,max(STAtr,[],2),num2str([1:length(PCresult{f}.dendaxis_unsorted)]'));
end

%% Measuring bAP amplitude
fprintf('Section: Measuring bAP amplitude\n');

% Cache: this section is expensive, so its results are saved to disk and reloaded
% on subsequent runs. Delete bAPampDataFile to force recomputation. PCresult is
% NOT cached here (it holds the full traces), so the "Load Result files" section
% must have run first on either path.
bAPampDataFile = fullfile(save_figto, 'Figure3_bAPamplitude_data.mat');
if isfile(bAPampDataFile)
    fprintf('  Loading cached results (skipping recompute): %s\n', bAPampDataFile);
    load(bAPampDataFile, 'bAPPropsMat', 'EventPropsMat', 'AlignedbAPall', ...
        'AlignedEvntall', 'AUCbAPall', 'SpikeTable', 'EventTable', 'NeuronTable', ...
        'nTau_EV', 'x_bin_edges', 'x_bin_edges2', ...
        'presub_readtime', 'PerisomaRange');
else

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
        PostSub_Basal PostSub_Distal ApicalAUCFirst'];
    EventPropsMat{f,1}=[Evlab];
    EventPropsMat{f,2}={TransRatebAPinEvent,AUCrawEV,AUCEVbin_cell,AlignedEvntcell_bin};

    fieldName={'SpikeFrame','IsBlue','IsNA','IsPF','IsRun','Length','AUC_centroid','TransmissionRate','AUC_apical','AUC_soma',...
               'ISI','Amp_apical','Spike_type','Nspike','PreSub_Soma','PreSub_Apical','TRmax_bAP','TRfrst_bAP','PreSub_Distal',...
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

% Consolidate every per-spike / per-event descriptor into sorted, tidy tables.
% SpikeTable (one row per somatic spike), EventTable (one row per SS/CS/BS event)
% and NeuronTable (one row per neuron) supersede the scattered bAPPropsMat /
% EventPropsMat variables. The full ROI x time aligned traces stay in
% AlignedbAPall / AlignedEvntall; slice them with (Neuron, LocalIdx).
fprintf('  Building unified SpikeTable / EventTable ...\n');
[SpikeTable, EventTable, NeuronTable] = build_spike_event_tables(PCresult, bAPPropsMat, EventPropsMat, foi);
fprintf('    %d spikes, %d events, %d neurons\n', height(SpikeTable), height(EventTable), height(NeuronTable));

fprintf('  Saving bAP results to %s\n', bAPampDataFile);
save(bAPampDataFile, 'bAPPropsMat', 'EventPropsMat', 'AlignedbAPall', ...
    'AlignedEvntall', 'AUCbAPall', 'SpikeTable', 'EventTable', 'NeuronTable', ...
    'nTau_EV', 'x_bin_edges', 'x_bin_edges2', ...
    'presub_readtime', 'PerisomaRange', 'foi', '-v7.3');
end

%% Coorelation mat of spike AUC
fprintf('Section: Correlation matrix of spike AUC\n');
figure(53); clf;
for f=foi
    St = SpikeTable(SpikeTable.Neuron==f,:);
    somaROI = ismember(PCresult{f}.roi_dClass,2);
    AUCraw    = cell2mat(St.AUC_perROI');          % nROI x nSpike
    AUCshort  = cell2mat(St.AUCshort_perROI');
    AUCraw1   = cell2mat(St.AUC_perROI_ch1');
    AUCraw2   = cell2mat(St.AUC_perROI_ch2');
    AUCshort1 = cell2mat(St.AUCshort_perROI_ch1');
    AUCshort2 = cell2mat(St.AUCshort_perROI_ch2');
    bAPmat  = AUCraw ./mean(AUCshort(somaROI,:),1,'omitnan');
    bAPmat1 = AUCraw1./mean(AUCshort1(somaROI,:),1,'omitnan');
    bAPmat2 = AUCraw2./mean(AUCshort2(somaROI,:),1,'omitnan');
    bAPCorr{f}=get_corrMat(bAPmat1,bAPmat2);
    nexttile([1 1])
    imagesc(bAPCorr{f},[-0.2 0.8])
end

%% Correlation between pre-spike sub and tx ratio
fprintf('Section: Correlation between pre-spike subthreshold and TX ratio\n');
%figure(46); clf;
f=23;
Et = EventTable(EventTable.Neuron==f,:);
evName = ["SS" "CS" "BS"];
clf; %tiledlayout(4,2);
sgtitle(['Subthreshold from ' num2str(presub_readtime(1)) ' ms to ' num2str(presub_readtime(end)) ' ms'])

cmapEv=hsv(3);
nexttile([1 1]);
scatter(Et.GlobalAmp, Et.SeesawAmp,20,vec2cmap(Et.TRfrst_bAP,'jet',[0 3]),'filled');
xlabel('Basal + distal'); ylabel('Basal - distal');
title('Pre-spike subthreshold (color: TX ratio)');
xlim([-3 5]/50); ylim([-4 4]/50);
colormap(jet);
cb=colorbar;
cb.Label.String='TX ratio';
cb.Ticks=[0 1];
cb.TickLabels={'0','3'};

nexttile([1 1]);
scatter(Et.PreSub_Basal, Et.PreSub_Distal,20,vec2cmap(Et.TRfrst_bAP,'jet',[0 3]),'filled');
xlabel('Basal'); ylabel('Distal');
title('Basal vs distal pre-sub (color: TX ratio)');
xlim([-3 5]/50); ylim([-4 4]/50);
colormap(jet);
cb=colorbar;
cb.Label.String='TX ratio';
cb.Ticks=[0 1];
cb.TickLabels={'0','3'};

nexttile([1 1]);
for sptype=1:3
    isc=Et.EventClass==evName(sptype);
    scatter(Et.PostSub_Basal(isc),Et.PostSub_Distal(isc),20,cmapEv(sptype,:),'filled'); hold all;
end
xlim([-2 5]/50); ylim([-3 6]/50);
xlabel('Post-spike subthreshold, basal'); ylabel('Post-spike subthreshold, distal');
title('Post-spike subthreshold by event type');
legend({'SS','CS','BS'});

nexttile([1 1]);
SeesawAmp=Et.SeesawAmp;
GlobalAmp=Et.GlobalAmp;
for sptype=1:3
    isc=Et.EventClass==evName(sptype);
    scatter(Et.AUCapical_frst(isc),SeesawAmp(isc),20,cmapEv(sptype,:),'filled'); hold all;
end
validInd = ~isnan(Et.AUCapical_frst) & ~isnan(SeesawAmp);
[r p]=corr(Et.AUCapical_frst(validInd),SeesawAmp(validInd));
xlabel('First bAP apical AUC'); ylabel('Basal - distal');
xlim([-4 17]/50); ylim([-5 5]/50);
title(sprintf('Coor: %.2f, p: %.2d',r,p))
legend({'SS','CS','BS'});

nexttile([1 1]);
for sptype=1:3
    isc=Et.EventClass==evName(sptype);
    scatter(Et.TRfrst_bAP(isc),SeesawAmp(isc),20,cmapEv(sptype,:),'filled'); hold all;
end
validInd = ~isnan(Et.TRfrst_bAP) & ~isnan(SeesawAmp);
[r p]=corr(Et.TRfrst_bAP(validInd),SeesawAmp(validInd));
xlabel('TX ratio'); ylabel('Basal - distal');
title(sprintf('Coor: %.2f, p: %.2d',r,p));
xlim([-0.5 5]); ylim([-5 5]/50);
legend({'SS','CS','BS'});

nexttile([1 1]);
for sptype=1:3
    isc=Et.EventClass==evName(sptype);
    scatter(Et.AUCapical_frst(isc),GlobalAmp(isc),20,cmapEv(sptype,:),'filled'); hold all;
end
validInd = ~isnan(Et.AUCapical_frst) & ~isnan(GlobalAmp);
[r p]=corr(Et.AUCapical_frst(validInd),GlobalAmp(validInd));
xlabel('First bAP apical AUC'); ylabel('Basal + distal');
xlim([-4 17]/50); ylim([-2 5]/50);
title(sprintf('Coor: %.2f, p: %.2d',r,p))
legend({'SS','CS','BS'});

nexttile([1 1]);
for sptype=1:3
    isc=Et.EventClass==evName(sptype);
    scatter(Et.TRfrst_bAP(isc),GlobalAmp(isc),20,cmapEv(sptype,:),'filled'); hold all;
end
validInd = ~isnan(Et.TRfrst_bAP) & ~isnan(GlobalAmp);
[r p]=corr(Et.TRfrst_bAP(validInd),GlobalAmp(validInd));
xlabel('TX ratio'); ylabel('Basal + distal');
title(sprintf('Coor: %.2f, p: %.2d',r,p))
xlim([-0.5 5]); ylim([-2 5]/50);
legend({'SS','CS','BS'});

nexttile([1 1]);
for sptype=1:3
    isc=Et.EventClass==evName(sptype);
    scatter(Et.TRfrst_bAP(isc),Et.PreSub_Basal(isc),20,cmapEv(sptype,:),'filled'); hold all;
end
validInd = ~isnan(Et.TRfrst_bAP) & ~isnan(Et.PreSub_Basal);
[r p]=corr(Et.AUCapical_frst(validInd),Et.PreSub_Basal(validInd));
xlabel('TX ratio'); ylabel('Basal');
xlim([-0.5 5]); ylim([-2 5]/50);
title(sprintf('Coor: %.2f, p: %.2d',r,p))
legend({'SS','CS','BS'});

nexttile([1 1]);
for sptype=1:3
    isc=Et.EventClass==evName(sptype);
    scatter(Et.TRfrst_bAP(isc),Et.PreSub_Distal(isc),20,cmapEv(sptype,:),'filled'); hold all;
end
validInd = ~isnan(Et.TRfrst_bAP) & ~isnan(Et.PreSub_Distal);
[r p]=corr(Et.TRfrst_bAP(validInd),Et.PreSub_Distal(validInd));
xlabel('TX ratio'); ylabel('Distal');
title(sprintf('Coor: %.2f, p: %.2d',r,p))
xlim([-0.5 5]); ylim([-2 5]/50);
legend({'SS','CS','BS'});

%% Sort by TX ratio
fprintf('Section: Sort events by TX ratio\n');
f=23;
Et = EventTable(EventTable.Neuron==f,:);
evName = ["SS" "CS" "BS"];
SeesawAmp=Et.SeesawAmp;
GlobalAmp=Et.GlobalAmp;
stype_str={'SS','CS','BS'};
prctbin=[0:0.2:1]*100; cax=[-1 6]/50; ax2=[]; show_interval=[-50:50]; minISI=50;
for sptype=1:3
    figure(47+sptype); clf;
    tiledlayout(5,4);
    for b=1:length(prctbin)-1
        nROI=size(AlignedEvntall{f,1},1);
        pctLabel=sprintf('%g-%g%%',prctbin(b),prctbin(b+1));
        sptypeInd=[Et.EventClass==evName(sptype)];
        STAkymo=mean(AlignedEvntall{f,1}(:,-nTau_EV(1)+show_interval,sptypeInd),3,'omitnan');
        prcInd=[SeesawAmp>=prctile(SeesawAmp,prctbin(b)) & SeesawAmp<prctile(SeesawAmp,prctbin(b+1)) ...
              & Et.EventClass==evName(sptype) & Et.ISI>minISI];
        nexttile((b-1)*4+1,[1 1]);
        imagesc(show_interval,[1:nROI],mean(AlignedEvntall{f,2}(:,-nTau_EV(1)+show_interval,prcInd),3,'omitnan'),cax)
        set_kymoYtick(PCresult{f}.dendaxis);
        xlabel('Time from event (ms)'); ylabel('Distance from soma (\mum)'); title(['Seesaw ' pctLabel]);
        ax2=[ax2 nexttile((b-1)*4+2,[1 1])];
        imagesc(show_interval,[1:nROI],mean(AlignedEvntall{f,2}(:,-nTau_EV(1)+show_interval,prcInd)-STAkymo,3,'omitnan'),[-2 2]/50)
        set_kymoYtick(PCresult{f}.dendaxis);
        xlabel('Time from event (ms)'); ylabel('Distance from soma (\mum)'); title(['Seesaw ' pctLabel ', \Delta mean']);

        prcInd=[GlobalAmp>=prctile(GlobalAmp,prctbin(b)) & GlobalAmp<prctile(GlobalAmp,prctbin(b+1)) ...
              & Et.EventClass==evName(sptype) & Et.ISI>minISI];
        nexttile((b-1)*4+3,[1 1]);
        imagesc(show_interval,[1:nROI],mean(AlignedEvntall{f,2}(:,-nTau_EV(1)+show_interval,prcInd),3,'omitnan'),cax)
        set_kymoYtick(PCresult{f}.dendaxis);
        xlabel('Time from event (ms)'); ylabel('Distance from soma (\mum)'); title(['Global ' pctLabel]);
        ax2=[ax2 nexttile((b-1)*4+4,[1 1])];
        imagesc(show_interval,[1:nROI],mean(AlignedEvntall{f,2}(:,-nTau_EV(1)+show_interval,prcInd)-STAkymo,3,'omitnan'),[-2 2]/50)
        set_kymoYtick(PCresult{f}.dendaxis);
        xlabel('Time from event (ms)'); ylabel('Distance from soma (\mum)'); title(['Global ' pctLabel ', \Delta mean']);
    end
    colormap(turbo);
    sgtitle(sprintf('%s events: subtracted kymographs by pre-spike subthreshold percentile', stype_str{sptype}));
end
arrayfun(@(x) colormap(x,gen_colormap([0 0.4 1;1 1 1;1 0 0])),ax2);

figure(51); clf;
tiledlayout(5,6); stype_str={'SS','CS','BS'}; ax2=[];
for sptype=1:3
    for b=1:length(prctbin)-1
        pctLabel=sprintf('%g-%g%%',prctbin(b),prctbin(b+1));
        prcInd=[Et.PreSub_Distal>=prctile(Et.PreSub_Distal,prctbin(b)) & Et.PreSub_Distal<prctile(Et.PreSub_Distal,prctbin(b+1))...
               & Et.EventClass==evName(sptype) & Et.ISI>minISI];
        sptypeInd=[Et.EventClass==evName(sptype)];
        nexttile((b-1)*6+(sptype-1)*2+1,[1 1]);
        imagesc(show_interval,[1:nROI],mean(AlignedEvntall{f,1}(:,-nTau_EV(1)+show_interval,prcInd),3,'omitnan'),cax)
        set_kymoYtick(PCresult{f}.dendaxis);
        xlabel('Time from event (ms)'); ylabel('Distance from soma (\mum)');
        title([stype_str{sptype} ', distal ' pctLabel]);
        STAkymo=mean(AlignedEvntall{f,1}(:,-nTau_EV(1)+show_interval,sptypeInd),3,'omitnan');
        ax2=[ax2 nexttile((b-1)*6+(sptype-1)*2+2,[1 1])];
        imagesc(show_interval,[1:nROI],mean(AlignedEvntall{f,1}(:,-nTau_EV(1)+show_interval,prcInd)-STAkymo,3,'omitnan'),[-2 2]/50)
        set_kymoYtick(PCresult{f}.dendaxis);
        xlabel('Time from event (ms)'); ylabel('Distance from soma (\mum)');
        title([stype_str{sptype} ', distal ' pctLabel ', \Delta mean']);
    end
    colormap(turbo);
end
arrayfun(@(x) colormap(x,gen_colormap([0 0.4 1;1 1 1;1 0 0])),ax2);
sgtitle('Raw kymographs by distal pre-spike subthreshold percentile');

%% Plot correlation between distal pre-sub and bAP apical AUC (figure 4);
fprintf('Section: Figure 4 - distal pre-sub vs bAP apical AUC\n');
f=23; cax=[-0.02 0.1]; minISI=50; rmvROI2show=4;
Et = EventTable(EventTable.Neuron==f,:);
vec2threshold=Et.AUCapical_frst_norm;
threshold2showLow=[prctile(vec2threshold,1), prctile(vec2threshold,15)];
threshold2showHigh=[prctile(vec2threshold,75), prctile(vec2threshold,99)];
LowInd=[vec2threshold<threshold2showLow(2) & vec2threshold>=threshold2showLow(1) & Et.ISI>minISI & Et.EventClass=='SS'];
HighInd=[vec2threshold<threshold2showHigh(2) & vec2threshold>=threshold2showHigh(1) & Et.ISI>minISI & Et.EventClass=='SS'];
dax=PCresult{f}.dendaxis; nROI=length(PCresult{f}.sortDist);
[daxsorted daxsortedInd]=sort(dax,'ascend');
LowTXSTA=mean(AlignedEvntall{f,1}(daxsortedInd,:,LowInd),3,'omitnan');
HighTXSTA=mean(AlignedEvntall{f,1}(daxsortedInd,:,HighInd),3,'omitnan');
% clf; imagesc(LowTXSTA(:,200+[-50:50]),cax)
%LowTXSTA=LowTXSTA-median(LowTXSTA(:,1:20),2);
%HighTXSTA=HighTXSTA-median(HighTXSTA(:,1:20),2);

figure(54); clf; tiledlayout(7,2,'Padding','tight','TileSpacing','tight'); 
ax3=[]; roi2show=setdiff([1:nROI],rmvROI2show);
ax3=[ax3 nexttile([4 1])];
nTau_baseline=[-100:-30];
imagesc([nTau_EV(1):nTau_EV(2)],[1:nROI],LowTXSTA(roi2show,:)-median(LowTXSTA(roi2show,-nTau_EV(1)+nTau_baseline),2),cax)
set_kymoYtick(daxsorted); ylabel('Distance from soma (µm)');  xlabel('Time (ms)');
title(sprintf('Average of \nlow transmission bAPs'));
ax3=[ax3 nexttile([4 1])];
imagesc([nTau_EV(1):nTau_EV(2)],[1:nROI],HighTXSTA(roi2show,:)-median(HighTXSTA(roi2show,-nTau_EV(1)+nTau_baseline),2),cax);
set(gca,'ytick',[]); xlabel('Time (ms)');
%set_kymoYtick(daxsorted);
title(sprintf('Average of \nhigh transmission bAPs'));
colormap(turbo); cb=colorbar; cb.Label.String='∆F/F';

ax3=[ax3 nexttile([3 1])];
cmap_AP=gen_colormap(Plasma,10); 
LowTXSTAtr=mean(cat(3,Et.AlignedEv_binned{LowInd}),3,'omitnan');
HighTXSTAtr=mean(cat(3,Et.AlignedEv_binned{HighInd}),3,'omitnan');
LowTXSTAtr=LowTXSTAtr-median(LowTXSTAtr(:,-nTau_EV(1)+nTau_baseline),2);
HighTXSTAtr=HighTXSTAtr-median(HighTXSTAtr(:,-nTau_EV(1)+nTau_baseline),2);

l=plot([nTau_EV(1):nTau_EV(2)],[LowTXSTAtr(2,:); LowTXSTAtr([4],:)]','LineWidth',1.5);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap_AP([1 10],:),2)); axis off tight; ylim([-0.04 0.14]);

ax3=[ax3 nexttile([3 1])];
l=plot([nTau_EV(1):nTau_EV(2)],[HighTXSTAtr(2,:); HighTXSTAtr([4],:)]','LineWidth',1.5);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap_AP([1 10],:),2)); axis off tight; ylim([-0.04 0.14]);
drawScaleBar(0.1,'vertical','position',[19.5 0.04],'color','k');
linkaxes(ax3,'x')
xlim([-50 21])
set_font('Arial'); set_fontsize(22); 
set_figsize(265,175);
% Correlation of first-bAP apical AUC vs pre-spike distal/basal subthreshold,
% per session (corrTXnPreSub) and resolved in time before the spike
% (corrTXnPreSub_timelapse). See get_presub_apicalAUC_corr at end of file.

[corrTXnPreSub, corrTXnPreSub_timelapse, t2avg] = get_presub_apicalAUC_corr(EventTable, AlignedEvntall, PCresult, foi, nTau_EV);

figure(55); clf;
tiledlayout(1,8,'padding','tight');
nexttile([1 3]);
f=23;
Et = EventTable(EventTable.Neuron==f,:);
validInd = ~isnan(Et.AUCapical_frst_norm) & ~isnan(Et.PreSub_Distal);
scatter(Et.PreSub_Distal,Et.AUCapical_frst_norm,20,[0 0 0],'filled','MarkerFaceAlpha',1); hold all;
%scatter_density(Et.PreSub_Distal,Et.TRfrst_bAP); hold all;
validInd = ~isnan(Et.AUCapical_frst_norm) & ~isnan(Et.PreSub_Distal);
[R, P] = corr(Et.AUCapical_frst_norm(validInd), Et.PreSub_Apical(validInd), 'Type', 'Pearson');
% least-squares fit + 95% confidence bounds, light red
lightRed=[1 0.6 0.6];
xq=linspace(min(Et.PreSub_Distal(validInd)),max(Et.PreSub_Distal(validInd)),100)';
mdlfit=fitlm(Et.PreSub_Distal(validInd),Et.AUCapical_frst_norm(validInd));
[yq,yqci]=predict(mdlfit,xq);
plot(xq,yqci,'--','Color',lightRed,'LineWidth',1);
plot(xq,yq,'-','Color',lightRed,'LineWidth',1.5);
ylabel('bAP apical AUC'); xlabel(['Pre-spike subthreshold' newline 'in distal dendrite']);
title(sprintf('R= %.2f, P: %.2d',R,P))
xlim([-0.1 0.1]); ylim([-0.1 2]);

nexttile([1 2]);
p=Boxplot_wPoints2(corrTXnPreSub(foi,[2 1]),cmap_AP([1 10],:)); ylim([-0.5 1]);
set(gca,'xtick',[1:2],'XTickLabel',{'Basal','Distal'}); box off;
ylabel(['Correlation between' newline 'bAP apical AUC and pre-spike subthreshold']);

nexttile([1 3]);
tax=mean([t2avg(2:end); t2avg(1:end-1)],1);
corrTXnPreSub_t_mean=squeeze(mean(corrTXnPreSub_timelapse(foi,:,:),1,'omitnan'));
corrTXnPreSub_t_std=squeeze(std(corrTXnPreSub_timelapse(foi,:,:),0,1,'omitnan'));
corrTXnPreSub_t_N=squeeze(sum(~isnan(corrTXnPreSub_timelapse(foi,:,:)),1));
corrTXnPreSub_t_sem=corrTXnPreSub_t_std./sqrt(corrTXnPreSub_t_N);
l=errorbar(tax,corrTXnPreSub_t_mean([2 1],:),corrTXnPreSub_t_sem([2 1],:),'linewidth',1.5)';
set(gca,'xdir','normal','YAxisLocation','right')
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap_AP([1 10],:),2)); box off;
ylabel(['Correlation between' newline 'bAP apical AUC and pre-spike subthreshold']);
xlabel('Time before spike (ms)')
set_font('Arial'); set_fontsize(16);
set_figsize(290,110);

%% Plot correlation between distal pre-sub and TX ratio (figure 4a);
fprintf('Section: Figure 4a - distal pre-sub vs TX ratio\n');
f=23; cax=[-0.02 0.1]; minISI=50; cax_tr=[-0.04 0.15];
Et = EventTable(EventTable.Neuron==f,:);
vec2threshold=Et.TRfrst_bAP_norm;
% threshold2showLow=[prctile(vec2threshold,1), prctile(vec2threshold,25)];
% threshold2showHigh=[prctile(vec2threshold,75), prctile(vec2threshold,99)];
% LowInd=[vec2threshold<threshold2showLow(2) & vec2threshold>=threshold2showLow(1) & Et.ISI>minISI & Et.EventClass=='SS'];
% HighInd=[vec2threshold<threshold2showHigh(2) & vec2threshold>=threshold2showHigh(1) & Et.ISI>minISI & Et.EventClass=='SS'];
ValidInd=find([Et.ISI>minISI & Et.EventClass=='SS']);
threshold2showLow=[prctile(vec2threshold(ValidInd),1), prctile(vec2threshold(ValidInd),20)];
threshold2showHigh=[prctile(vec2threshold(ValidInd),80), prctile(vec2threshold(ValidInd),99)];
LowInd=ValidInd([vec2threshold(ValidInd)<threshold2showLow(2) & vec2threshold(ValidInd)>=threshold2showLow(1)]); 
HighInd=ValidInd([vec2threshold(ValidInd)<threshold2showHigh(2) & vec2threshold(ValidInd)>=threshold2showHigh(1)]); 

omitROI=[2 22 33 38]; 
showROI=setdiff([1:size(AlignedEvntall{f,2},1)],omitROI);
daxshowROI=PCresult{f}.dendaxis(showROI);
[sorteddax, dsort]=sort(daxshowROI,'ascend');

nTau_baseline=[-50:-30];

LowTXSTA=mean(AlignedEvntall{f,1}(showROI,:,LowInd),3,'omitnan');
HighTXSTA=mean(AlignedEvntall{f,1}(showROI,:,HighInd),3,'omitnan');
LowTXSTA=LowTXSTA-median(LowTXSTA(:,-nTau_EV(1)+nTau_baseline),2);
HighTXSTA=HighTXSTA-median(HighTXSTA(:,-nTau_EV(1)+nTau_baseline),2);

figure(54); clf; tiledlayout(7,2,'TileSpacing','compact'); ax3=[]; 
ax3=[ax3 nexttile([4 1])];
imagesc([nTau_EV(1):nTau_EV(2)],[1:length(showROI)],LowTXSTA(dsort,:),cax)
set_kymoYtick(sorteddax); ylabel('Distance from soma (µm)');
title(sprintf('Average of \nlow transmission bAPs'));

ax3=[ax3 nexttile([4 1])];
imagesc([nTau_EV(1):nTau_EV(2)],[1:length(showROI)],HighTXSTA(dsort,:),cax);
set(gca,'ytick',[]); xlabel('Time (ms)');
title(sprintf('Average of \nhigh transmission bAPs'));
colormap(turbo); cb=colorbar; cb.Label.String='∆F/F';

ax3=[ax3 nexttile([3 1])];
cmap_AP=gen_colormap(Plasma,10);
LowTXSTAtr=mean(cat(3,Et.AlignedEv_binned{LowInd}),3,'omitnan');
HighTXSTAtr=mean(cat(3,Et.AlignedEv_binned{HighInd}),3,'omitnan');
LowTXSTAtr=LowTXSTAtr-median(LowTXSTAtr(:,-nTau_EV(1)+nTau_baseline),2);
HighTXSTAtr=HighTXSTAtr-median(HighTXSTAtr(:,-nTau_EV(1)+nTau_baseline),2);

l=plot([nTau_EV(1):nTau_EV(2)],[LowTXSTAtr(2,:); LowTXSTAtr([4],:)]','LineWidth',1.5);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap_AP([1 10],:),2)); axis off tight; ylim(cax_tr);

ax3=[ax3 nexttile([3 1])];
l=plot([nTau_EV(1):nTau_EV(2)],[HighTXSTAtr(2,:); HighTXSTAtr([4],:)]','LineWidth',1.5);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap_AP([1 10],:),2)); axis off tight; ylim(cax_tr);
linkaxes(ax3,'x')
%drawScaleBar(50,'horizontal','Position',[19 -0.5]); drawScaleBar(5,'vertical','Position',[19 -0.5]); 
xlim([-50 20]); 
set_font('Arial'); set_fontsize(22); 
set_figsize(265,175);


%% Correlation of first-bAP apical AUC vs pre-spike distal/basal subthreshold,
% per session (corrTXnPreSub) and resolved in time before the spike
% (corrTXnPreSub_timelapse). See get_presub_apicalAUC_corr at end of file.
[corrTXnPreSub, corrTXnPreSub_timelapse, t2avg] = get_presub_apicalAUC_corr(EventTable, AlignedEvntall, PCresult, foi, nTau_EV);

figure(55); clf;
tiledlayout(1,8,'padding','tight','TileSpacing','tight');
nexttile([1 3]);
f=23;
Et = EventTable(EventTable.Neuron==f,:);
validInd = ~isnan(Et.AUCapical_frst_norm) & ~isnan(Et.PreSub_Distal);
scatter(Et.PreSub_Distal,Et.AUCapical_frst_norm,20,[0 0 0],'filled','MarkerFaceAlpha',1); hold all;
%scatter_density(Et.PreSub_Distal,Et.TRfrst_bAP); hold all;
validInd = ~isnan(Et.AUCapical_frst_norm) & ~isnan(Et.PreSub_Distal);
[R, P] = corr(Et.AUCapical_frst_norm(validInd), Et.PreSub_Distal(validInd), 'Type', 'Pearson');
% least-squares fit + 95% confidence bounds, light red
lightRed=[1 0.6 0.6];
xq=linspace(min(Et.PreSub_Distal(validInd)),max(Et.PreSub_Distal(validInd)),100)';
mdlfit=fitlm(Et.PreSub_Distal(validInd),Et.AUCapical_frst_norm(validInd));
[yq,yqci]=predict(mdlfit,xq);
plot(xq,yqci,'--','Color',lightRed,'LineWidth',1);
plot(xq,yq,'-','Color',lightRed,'LineWidth',1.5);
ylabel('TX ratio'); xlabel(['Pre-spike subthreshold' newline 'in distal dendrite (∆F/F)']);
title(sprintf('\\rho = %.2f, \\it{P} = %.2d',R,P))
xlim([-0.09 0.09]); ylim([0 1.9]);

nexttile([1 2]);
p=Boxplot_wPoints2(corrTXnPreSub(foi,[2 1]),cmap_AP([1 10],:)); 
drawPValueLines(p, 0, 'TextYOffset', 0.05); box off;
ylim([-0.5 1]);
set(gca,'xtick',[1:2],'XTickLabel',{'Basal','Distal'}); box off;
ylabel(['Correlation between' newline 'TX ratio and pre-spike subthreshold']);

nexttile([1 3]);
tax=mean([t2avg(2:end); t2avg(1:end-1)],1);
corrTXnPreSub_t_mean=squeeze(mean(corrTXnPreSub_timelapse(foi,:,:),1,'omitnan'));
corrTXnPreSub_t_std=squeeze(std(corrTXnPreSub_timelapse(foi,:,:),0,1,'omitnan'));
corrTXnPreSub_t_N=squeeze(sum(~isnan(corrTXnPreSub_timelapse(foi,:,:)),1));
corrTXnPreSub_t_sem=corrTXnPreSub_t_std./sqrt(corrTXnPreSub_t_N);
l=errorbar(tax(2:end),corrTXnPreSub_t_mean([2 1],2:end),corrTXnPreSub_t_sem([2 1],2:end),'linewidth',1.5)';
p_t=NaN(length(tax),1);
for tt=2:length(tax)
pp=get_pValue([squeeze(corrTXnPreSub_timelapse(foi,1,tt)),squeeze(corrTXnPreSub_timelapse(foi,2,tt))],1);
p_t(tt)=pp(1,2);
end
set(gca,'xdir','normal','YAxisLocation','right')
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap_AP([1 10],:),2)); box off;

% significance star slightly above the distal (orange) errorbar where p_t < 0.05
% distal mean/sem = row 1 of corrTXnPreSub_t_*  (basal = row 2)
hold on; yl=ylim; yoff=0.05*diff(yl); maxStarY=yl(2);
for tt=2:length(tax)
    if p_t(tt)<0.05
        if p_t(tt)<0.001, star='***'; elseif p_t(tt)<0.01, star='**'; else, star='*'; end
        yStar=corrTXnPreSub_t_mean(1,tt)+corrTXnPreSub_t_sem(1,tt)+yoff;
        text(tax(tt),yStar,star,'Color',[0 0 0],'FontWeight','bold','FontSize',12, ...
             'HorizontalAlignment','center','VerticalAlignment','bottom');
        maxStarY=max(maxStarY,yStar);
    end
end
ylim([yl(1) max(yl(2),maxStarY+yoff)]);
xlim([-150 0]);

ylabel(['Correlation between' newline 'TX ratio and pre-spike subthreshold']);
xlabel('Time before spike (ms)');
set_font('Arial'); set_fontsize(19);
set_figsize(320,140);

%% Get variance in bAP success/fail
fprintf('Section: Variance in bAP success/fail\n');
% figure(53); clf; tiledlayout('padding','tight');
Varfraction=[]; 
for f=foi([1:end])
% branch-normalized bAP AUC matrix (soma rows use short AUC), non-blue spikes only
[bAPmat, NonBlueInd] = build_branchNorm_bAPmat(SpikeTable, PCresult, f);
St = SpikeTable(SpikeTable.Neuron==f,:);

sortedbranchInd=PCresult{f}.BranchLabel(PCresult{f}.sortDist);
bAPmat_branch=NaN(max(sortedbranchInd),size(bAPmat,2));
ftprnts_branch=NaN(size(PCresult{f}.Ftprnts,1),size(PCresult{f}.Ftprnts,2),max(sortedbranchInd));
ftprnts_branch_sorted=PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist);
for b=1:max(sortedbranchInd)
    bindtmp=find(sortedbranchInd==b);
if ~isempty(bindtmp)
bAPmat_branch(b,:)=mean(bAPmat(bindtmp,:),1,'omitnan');    
ftprnts_branch(:,:,b)=max(ftprnts_branch_sorted(:,:,bindtmp),[],3);
end
end
omitInds=[find(sum(isnan(bAPmat_branch),2)>size(bAPmat_branch,2)/2)];
bAPmat_branch(omitInds,:)=[];
ftprnts_branch(:,:,omitInds)=[];

TXthres=[prctile(St.TransmissionRate,75) prctile(St.TransmissionRate,25)];
ind2Successtemplate=St.TransmissionRate(NonBlueInd)>TXthres(1); ind2Failtemplate=St.TransmissionRate(NonBlueInd)<TXthres(2);
%ind2Successtemplate=eigt_p(1,:)>TXthres(1); ind2Failtemplate=eigt_p(1,:)<TXthres(2);
Successtemplate=mean(bAPmat_branch(:,ind2Successtemplate),2,'omitnan')-mean(bAPmat_branch(:,ind2Failtemplate),2,'omitnan');
Successtemp_trace=Successtemplate'*bAPmat_branch;

[V D eigt validI] =get_eigvector(bAPmat_branch);

swc2show=PCresult{f}.swc;
swc2show(:,4)=2; swc2show(1,4)=15;
% figure(5); clf;
% nexttile([1 1]);
% showScaleScatter(Successtemplate,swc2show,ftprnts_branch,turbo(256),[prctile(Successtemplate(:),5) prctile(Successtemplate(:),95)])
% title(sprintf('TX template'))
% for pc=1:11
%     nexttile([1 1]);
%     showScaleScatter(V(:,pc),swc2show,ftprnts_branch,turbo(256),[prctile(V(:),5) prctile(V(:),95)])
%     title(sprintf('PC: %1.0f, Fraction: %.2f %',pc,D(pc)/sum(D)))
% end

dax=PCresult{f}.dendaxis;

validI=find(validI);
TXratio_sub=St.TransmissionRate(NonBlueInd);
TXratio_sub=TXratio_sub(validI);
%validIJ=find(abs(zscore(TXratio_sub./eigt(:,1)))<4);
validIJ=[1:length(validI)];

[Varfraction(f,1)] = get_variance(bAPmat_branch(:,validI(validIJ)),eigt(validIJ,1)');
Varfraction(f,5)=D(1)/sum(D);
[Varfraction(f,2)] = get_variance(bAPmat_branch(:,validI(validIJ)),St.AUC_apical(validI(validIJ))');
[Varfraction(f,3)] = get_variance(bAPmat_branch(:,validI(validIJ)),St.TransmissionRate(validI(validIJ))');
[Varfraction(f,4)] = get_variance(bAPmat_branch(:,validI(validIJ)),Successtemp_trace(validI(validIJ)));

% AUCprnt=bAPmat(:,validI(validIJ))'\bAPPropsMat{f,1}.AUC_apical(validI(validIJ));
% TXprnt=bAPmat(:,validI(validIJ))'\bAPPropsMat{f,1}.TransmissionRate(validI(validIJ));
% 
% [Varfraction(f,1)] = get_variance(bAPmat(:,validI(validIJ))',V(:,1)');
% [Varfraction(f,2)] = get_variance(bAPmat(:,validI(validIJ))',AUCprnt');
% [Varfraction(f,3)] = get_variance(bAPmat(:,validI(validIJ))',TXprnt');
% [Varfraction(f,4)] = get_variance(bAPmat(:,validI(validIJ))',Successtemplate');
end

figure(52); clf;
%Boxplot_wPoints2(Varfraction(foi,[1 4 2 3]),[1 0 0])
% set(gca,'xtick',[1:4],'XTickLabel',{'PC1','High TX template','AUC apical','TX ratio'})
Boxplot_wPoints2(Varfraction(foi,4),[0 0 0])
set(gca,'xtick',[],'FontSize',13); ylim([0 1]); xlim([0.5 3.5]);
ylabel('Variance explained'); box off;
fprintf('Variance explained (branch averaged), %2.2f ± %2.2f, mean ± s.d.\n', mean(Varfraction(foi,4)),std(Varfraction(foi,4)/sqrt(length(foi))))

corrApical=NaN(1,max(foi));
omitN=[24 27]; dlimit=250;
for f=setdiff(foi,[])

% branch-normalized bAP AUC matrix (soma rows use short AUC), non-blue spikes only
[bAPmat, NonBlueInd] = build_branchNorm_bAPmat(SpikeTable, PCresult, f);
dcclass=PCresult{f}.roi_dClass;
dax=PCresult{f}.dendaxis;

bAPCorrMat=get_corrMat(bAPmat,bAPmat);
uppertri_ind=triu(ones(size(bAPCorrMat)),1);
label_ind=ismember(dcclass,[5])'*ismember(dcclass,[5]);
Ndend(f)=max(sum(uppertri_ind & label_ind,2));
Ndendpair(f)=sum(uppertri_ind & label_ind,[1 2]);
%label_ind=(dax>dlimit)'*(dax>dlimit);
cax=tovec(bAPCorrMat(uppertri_ind & label_ind));
corrApical(f)=mean(cax);
end
fprintf('Apical bAP correlation: %2.2f ± %2.2f, mean ± s.e.m.\n', mean(corrApical(foi),'omitnan'),std(corrApical(foi),'omitnan')/sqrt(length(foi)))

%% Show TX ratio STAs (Figure 4F)
fprintf('Section: Show TX ratio STAs\n');
figure(16); clf;
St = SpikeTable(SpikeTable.Neuron==f,:);
prcax=[0:20:100]; ax2=[];
tiledlayout(3,length(prcax)-1,'Padding','tight','TileSpacing','tight')
SWCpoints=PCresult{f}.swc;
SWCpoints(1,4)=40;
SWCpoints(:,4)=SWCpoints(:,4)+6;
cmap_ExTr=gen_colormap(Plasma,10);
isSS = St.SpikeClass=='SS';
for p=1:length(prcax)-1
    rangeVal=[prctile(St.TransmissionRate(isSS),prcax(p)) prctile(St.TransmissionRate(isSS),prcax(p+1))];
    showsp=find(St.TransmissionRate>rangeVal(1) & St.TransmissionRate<rangeVal(2) & isSS);
    nexttile(p,[2 1])
    ShowAmp=mean(AlignedbAPall{f,2}(:,:,showsp),3,'omitnan');
    %ShowAmp=ShowAmp-median(ShowAmp(:,1:200),2,'omitnan');
    %ShowAmp=max(ShowAmp(:,200:205),[],2);
    ShowAUC=cell2mat(St.AUC_perROI(showsp)');%./mean(AUCbAPall{f,1}(PCresult{f}.branch_dClass==2,showsp),1,'omitnan');
    ShowAUC=mean(ShowAUC,2,'omitnan');
    showScaleScatter(ShowAUC, SWCpoints, PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist),'turbo',[0 0.3]);
    axis equal tight off
    % nexttile([1 1])

    %ShowbAPTrace=ShowbAPTrace-median(ShowbAPTrace(1:30,:),1,'omitnan');
    % plot(ShowbAPTrace(:,:))
    if p==1
        bar100um=100/PCresult{f}.pixelsize;
        h = drawScaleBar(bar100um, 'vertical','color','k','position',[150 60]);
    end
    ax2=[ax2 nexttile(2*(length(prcax)-1)+p,[1 1])];
    ShowTr=mean(cat(3,St.AlignedbAP_binned{showsp}),3,'omitnan'); ShowTr=ShowTr(:,150:250);
    ShowTr=ShowTr-median(ShowTr(:,1:40),2);
    l=plot(ShowTr([2 4],20:80)','linewidth',1.5);
    arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap_ExTr([1 10],:),2));
    % if p==1   % label the two traces by the mean ROI distance from soma in each bin
    %     dax_f=PCresult{f}.dendaxis;   % same bin rule as zscore_binning: edge(k) <= d < edge(k+1)
    %     binMeanDist=arrayfun(@(k) mean(dax_f(dax_f>=x_bin_edges(k) & dax_f<x_bin_edges(k+1)),'omitnan'), ...
    %                          1:numel(x_bin_edges)-1);
    %     legend(arrayfun(@(d) [num2str(round(d)) ' \mum'], binMeanDist([2 4]), ...
    %            'UniformOutput',false),'Box','off','Location','northeastoutside');
    % end
    if p==length(prcax)-1
        drawScaleBar(10,'horizontal','color','k','position',[60, 0.05]);
        drawScaleBar(0.1,'vertical','color','k','position',[60, 0.05]);
    end
    axis off;
end
linkaxes(ax2); ylim([-0.02 0.15]);
nexttile(length(prcax)-1,[2 1])
cb=colorbar;
cb.Label.String = 'AUC (A.U.)';
set_font('Arial'); set_fontsize(16);
set_figsize(275,150);

%% Event props, Linear discrimination analysis
fprintf('Section: Event properties - linear discrimination analysis\n');
f=23;
Et = EventTable(EventTable.Neuron==f,:);
A = Et.AUCapical_frst_norm;
% Scan every numeric predictor column against A (skip identifiers, flags and the
% target itself). isnumeric drops the categorical / cell columns automatically.
excludeVar = {'EventUID','Neuron','LocalIdx','Frame','Time_s','Length','Nspike', ...
              'IsBlue','IsNA','IsPF','IsRun','IsValid','FirstSpikeUID','AUCapical_frst_norm'};
isNumCol = varfun(@(v) isnumeric(v) && size(v,2)==1, Et, 'OutputFormat','uniform');
VariableName = Et.Properties.VariableNames(isNumCol & ~ismember(Et.Properties.VariableNames, excludeVar));
features = Et{:, VariableName};
[n, t] = size(features);  % features is N x T matrix
rsq = zeros(1, t);
for i = 1:t
    mdl = fitlm(features(:,i), A);
    rsq(i) = mdl.Rsquared.Ordinary;
end
[~, bestIdx] = max(rsq);
fprintf('Best variable: %s (R^2 = %.3f)\n', VariableName{bestIdx}, rsq(bestIdx));

figure(19); clf;
bar(rsq)
set(gca,'xtick',[1:size(features,2)],'xticklabel',VariableName,'TickLabelInterpreter', 'none')
ylabel('R-squared')


%% Plot bAP amplitude versus spike order (figure 4, representative case)
fprintf('Section: Figure 4 - bAP amplitude vs spike order\n');
f=23;
St = SpikeTable(SpikeTable.Neuron==f,:);
Et = EventTable(EventTable.Neuron==f,:);
TrRateEdge=[-0.005:0.002:0.036]; CShist=[]; CSpoints=[]; BSpoints=[];
showInd=find(~St.IsBlue & ~St.IsNA & St.IsIsolated);
SSpoints={St.AUC_apical_norm(showInd), St.AUC_soma(showInd), St.TransmissionRate_norm(showInd)};
Ind2show=[];
for Spikeorder2show=[1:6]
    Ind2show{Spikeorder2show}=find(~St.IsBlue & ~St.IsNA & St.SpikeOrder==Spikeorder2show & ~St.IsIsolated & St.TransmissionRate<5.8);
    BSpoints{Spikeorder2show,1}=St.AUC_apical_norm(Ind2show{Spikeorder2show});
    BSpoints{Spikeorder2show,2}=St.AUC_soma(Ind2show{Spikeorder2show});
    BSpoints{Spikeorder2show,3}=St.TransmissionRate_norm(Ind2show{Spikeorder2show});
end
cmap_spikeorder=gray(10);
figure(10); clf; tiledlayout(1,2);
ax1=[];
nexttile([1 1]);
CSind=Et.EventClass=='CS'; SSind=Et.EventClass=='SS';
Boxplot_wPoints2({St.TransmissionRate_norm(showInd),Et.TRfrst_bAP_norm(CSind)},hsv(2))
%Boxplot_wPoints2({Et.TRfrst_bAP_norm(SSind),Et.TRfrst_bAP_norm(CSind)},hsv(2))
ylim([-1 5]/3)
set(gca,'XTick', [1 2],'XTickLabel',{'Isolated Spike','1st spike of CS'},'Box','off');
ylabel('Transmission ratio')
nexttile([1 1]);
show_data=[SSpoints{1,3} BSpoints(:,3)'];
show_data=cellfun(@(x) subsasgn(x, substruct('()', {find(x>2.6)}), NaN), show_data, 'UniformOutput', false);
p=Violin_wPoints(show_data,cmap_spikeorder(1:7,:));
orderstring=counting_string(1:5);
set(gca,'xtick',[1:6],'XTickLabel',[{'Isolated'} orderstring]);
ylabel('Normalized transmission ratio'); xlim([0.5 6.5]); ylim([0 3.2]);

% ---- significance vs the Isolated and 1st-spike baselines (stacked stars) ----
% show_data groups: 1 Isolated, 2 1st, 3 2nd, 4 3rd, 5 4th, 6 5th. Two compact
% rows of stars sit above each 2nd-5th violin (top = vs Isolated, bottom = vs
% 1st), labeled once on the left; Isolated-vs-1st is a short bracket. See
% pair_p / sig_star at end of file.
hold on;
testX = 3:6;  yA = 3.05;  yB = 2.82;
text(2.6, yA, 'vs Isolated','HorizontalAlignment','right','FontSize',9,'Color',[0 0 0]);
text(2.6, yB, 'vs 1st',     'HorizontalAlignment','right','FontSize',9,'Color',[0.45 0.45 0.45]);
for x = testX
    text(x, yA, sig_star(pair_p(show_data{1},show_data{x})), ...
        'HorizontalAlignment','center','FontWeight','bold','FontSize',11,'Color',[0 0 0]);
    text(x, yB, sig_star(pair_p(show_data{2},show_data{x})), ...
        'HorizontalAlignment','center','FontWeight','bold','FontSize',11,'Color',[0.45 0.45 0.45]);
end
yC = 2.35;
plot([1 1 2 2], yC+[0 0.06 0.06 0], 'k','LineWidth',1);
text(1.5, yC+0.08, sig_star(pair_p(show_data{1},show_data{2})), ...
    'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10);
set_font('Arial'); set_fontsize(16);
set_figsize(150,110);
%% Distribution of normalized TX ratio across neurons
fprintf('Section: Distribution of normalized TX ratio across neurons\n');
figure(12); clf;
TXratiobins=[0:0.1:2.5]
distTXratio_norm=NaN(max(foi),length(TXratiobins)-1);
for f=foi
St = SpikeTable(SpikeTable.Neuron==f,:);
showInd=find(~St.IsBlue & ~St.IsNA & St.TransmissionRate<5.5);
distTXratio_norm(f,:)=histcounts(St.TransmissionRate_norm(showInd),TXratiobins,'Normalization','probability');
xlabel('Normalized Transmission ratio'); ylabel('Probability'); xlim([0.5 7.5]/3);
end
plot(edge2bin(TXratiobins),distTXratio_norm','color',[0.6 0.6 0.6]); hold all;
errorbar_shade(edge2bin(TXratiobins),mean(distTXratio_norm,1,'omitnan'),std(distTXratio_norm,0,1,'omitnan'),[0 0 0]);
xlabel(sprintf('Normalized \ntransmission ratio')); ylabel('Probability');
xlim([0 2.5]); box off;
ylim([0 0.25]);
set_font('Arial'); set_fontsize(16);
set_figsize(80,117);

%% Representative cases (figure 3F)
fprintf('Section: Representative cases (figure 3f)\n');
f=23; %load(fullfile(fpath{f},'PC_Result.mat')) % load the result file
St = SpikeTable(SpikeTable.Neuron==f,:);
figure(17); clf; t=tiledlayout(2,4,'TileSpacing','compact');
ShowInd=[602 329 530 365 525 479 785 705]-3; show_t=[190:210]; omitROI=[1 2 3 5 21 37 13 12 11]; cax=[-0.05 0.15];

showROI=setdiff([1:size(AlignedEvntall{f,2},1)],omitROI);
daxshowROI=PCresult{f}.dendaxis(showROI);
[sorteddax, dsort]=sort(daxshowROI,'ascend');
g=1;
for sp=ShowInd
    showTr=St.AlignedbAP_binned{sp}([2 4],:);
    showTr=showTr-movmedian(showTr,50,2);
    % l=plot(show_t+nTau(1),showTr(:,show_t)','LineWidth',1.5); box off;
    % arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap_ExTr([1 10],:),2));
    %title(['Tx rate: ' num2str(St.TransmissionRate(sp))]);
    nexttile(g,[1 1])
    %filteredkymo=AlignedbAPall{f,2}(showROI,:,sp);
    ROIvec=Result.maintrunk;
    ROIvec_distorder=find_index_bh(PCresult{f}.sortDist,ROIvec);
    dax=PCresult{f}.dendaxis(ROIvec_distorder); tax=nTau;
    dendriteaxis_bin=linspace(min(dax),max(dax),17);

    filteredkymo=pcafilterTrace(AlignedbAPall{f,2}(showROI,:,sp),[1:20]);
    %filteredkymo=AlignedbAPall{f,2}(ROIvec_distorder,:,sp);
    %filteredkymo=pcafilterTrace(filteredkymo,[1:6]);

    [~, sortind]=sort(dax,'ascend');
    [Xq, Yq] = meshgrid([tax], dendriteaxis_bin);
    filteredkymo_main = interp2([tax], dax(sortind), filteredkymo(sortind,:), Xq, Yq, 'linear');

    imagesc(show_t+nTau(1),[1:length(showROI)],filteredkymo(dsort,show_t),cax);
    %imagesc(show_t+nTau(1),dendriteaxis_bin,filteredkymo_main(:,show_t),cax);
    if ismember(g,[1:4]);
        set(gca,'xtick',[]);
    end
    if ismember(g,[1 5]);
        set_kymoYtick(sorteddax);
        set(gca,'FontSize',12);
    else
    colormap(turbo); set(gca,'ytick',[],'FontSize',12);
    end
    %set(gca,'ytick',[1:length(dsort)],'YTickLabel',dsort)
    g=g+1;
end
% single colorbar on the right, labeled dF/F
cb=colorbar; cb.Layout.Tile='east'; cb.Label.String='∆F/F';
cb.Label.FontSize=16; 
% distance-from-soma y ticks + row titles on the left column
% shared axis labels: ylabel only on the left, xlabel under the bottom row
ylabel(t,'Distance from soma (µm)','FontSize',16,'FontName','Arial');
xlabel(t,'Time (ms)','fontsize',16);
% tick font size
set_font('Arial'); %set_fontsize(16);
set_figsize(192,107);

%% Representative cases (figure 4H)
fprintf('Section: Representative cases (figure 4H)\n');
figure(18); clf; t=tiledlayout(2,3,'padding','tight');
ShowInd=[223 335 55 330 16 144]-1;
show_t=[180:250]; omitROI=[2 3 13 16 21 37 33]; cax=[-0.04 0.15]; g=1; 
showROI=setdiff([1:size(AlignedEvntall{f,2},1)],omitROI);

daxshowROI=PCresult{f}.dendaxis(showROI);
[sorteddax, dsort]=sort(daxshowROI,'ascend');
for sp=ShowInd
    nexttile(g,[1 1])
    % showTr=bAPPropsMat{f,2}.AlignedbAP_binned{1}([2 4],:,sp);
    % showTr=showTr-movmedian(showTr,30,2);
    % l=plot([-10:55],showTr(:,190:255)','LineWidth',1.5); box off;
    % arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap_ExTr([1 10],:),2));
    % showTr=bAPPropsMat{f,2}.AlignedbAP_binned{1}([2 4],:,sp);
    filteredkymo=pcafilterTrace(AlignedEvntall{f,2}(showROI,:,sp),[1:20]);
    imagesc(show_t+nTau_EV(1),[1:length(showROI)],filteredkymo(dsort,show_t),cax);
    set(gca,'ytick',[]);
    %set(gca,'ytick',[1:length(showROI)],'yticklabel',dsort);
    colormap(turbo); 
    g=g+1;
end
% left column: distance-from-soma y ticks
nexttile(1); set_kymoYtick(sorteddax);
nexttile(4); set_kymoYtick(sorteddax);
% remove x ticks on the top row
for tl=1:3, nexttile(tl); set(gca,'xtick',[]); end
% one title per row (centered over the middle tile)
nexttile(2); title('Doublets');
nexttile(5); title('Complex spikes'); xlabel('Time (ms)');
% single shared y-axis label for both rows
ylabel(t,'Distance from soma (µm)','fontsize',20,'FontName','Arial');
% single colorbar on the far-right, labeled dF/F
axr=nexttile(6); cb=colorbar(axr); cb.Layout.Tile='east'; cb.Label.String='\DeltaF/F';
cb.Ticks=[-0.04 0 0.04 0.08 0.12];
set_font('Arial'); set_fontsize(20);
set_figsize(300, 170);
%% Multiple neurons (Figure 4G)
fprintf('Section: Figure 4G - multiple neurons\n');

figure(11); clf; tiledlayout('horizontal');
showDat=[]; showDatall=[]; bAPdecayConst=[];
M=[]; N=[]; S=[]; t_edge=[0.5:0.25:1.6 2 2.6 3.2 3.8]; R2=[];

% Per-neuron interactive exponential fits of the raw (unbinned) ISI vs bAP
% apical AUC data, replacing expfitDM_2. One fit per neuron, cached so the GUIs
% only run once; delete isiFitFile to refit. Model has no offset (decays to 0),
% so fit AUC-1 (add 1 back for plotting); X is ISI in ms.
isiFitFile = fullfile(save_figto,'Figure4G_ISIfit.mat');
haveFit = isfile(isiFitFile);
if haveFit, Sload=load(isiFitFile,'isiFit'); isiFit=Sload.isiFit; else, isiFit=cell(1,max(foi)); end

for f=foi
    St = SpikeTable(SpikeTable.Neuron==f,:);
    showInd=find(~St.IsBlue & ~St.IsNA);
    showDat=[log10(St.ISI(showInd)) St.AUC_apical_norm(showInd)];
    showDat((showDat(:,2)>50/3 | showDat(:,2)<-5/3),:)=[];
    %[M(f,:) S(f,:) t_center Ntmp]=binning_data_median({showDat},t_edge);
    bR=binning_data({showDat},t_edge);
    M(f,:)=bR.mean;
    S(f,:)=bR.std;
    t_center=bR.centers;
    showDatall=[showDatall; showDat];

    if ~haveFit    % interactively fit this neuron's raw showDat
        Xf=10.^showDat(:,1); Yf=showDat(:,2)-1; v=~isnan(Xf)&~isnan(Yf);
        if nnz(v)>=5
            optFit=struct('t_guess',20,'tau_lb',0,'tau_ub',Inf,'predefined_idx',[(Xf(v)<250) & (Xf(v)>2)]);
            fprintf('  Interactive ISI fit, neuron %d ...\n', f);
            isiFit{f}=interactive_expfit(Xf(v),Yf(v),optFit);
        end
        if isiFit{f}.aborted
            break;
        end
    end
    if ~isempty(isiFit{f}) && ~isempty(isiFit{f}.tau) && ~all(isnan(isiFit{f}.tau))
        bAPdecayConst(f)=max(isiFit{f}.tau);
        R2(f)=isiFit{f}.Rsq;
    else
        bAPdecayConst(f)=NaN; R2(f)=NaN;
    end
    %N(f,:)=sqrt(cellfun(@sum,Ntmp)');
end
if ~haveFit
    save(isiFitFile,'isiFit');
    fprintf('  Saved %d per-neuron ISI fits: %s\n', numel(foi), isiFitFile);
end

% 
% nexttile([1 1]);
% dsc=densityScatterChart(showDatall(:,1),showDatall(:,2));
% [Show_M Show_S t_center Show_N]=binning_data({showDatall},t_edge);
% dsc.ColorbarLabel = 'Probability';
% [~, ax2, scat]=unmanage(dsc); scat.SizeData=10; hold all;
% set(ax2, 'XTick', [0 1 2 3 4],'XTickLabel',10.^[0 1 2 3 4],'Box','off');
% errorbar(ax2, t_center, Show_M, Show_S ./ sqrt(cellfun(@sum,Show_N)'), 'Color', [1 0 0], 'LineWidth', 1.5);
% xlabel('Time since previous spike (ms)'); ylabel('Normalized transmission rate');
% ylim([-2 10]);

nexttile([1 1]);
plot(t_center,M(foi,:)','color',[.6 .6 .6]); hold all
Sall=sem(M(foi,:),1);
errorbar(t_center,mean(M(foi,:),1,'omitnan'),Sall,'color','r','LineWidth',2.5)
% dotted black line = mean of the per-neuron interactive fits, reconstructed on a
% common ISI grid from amp/tau (+offset); fit is on AUC-1, so add 1 back
xcommon=logspace(1,3,200)';
Yhat=nan(numel(xcommon),max(foi));
for f=foi
    r=isiFit{f};
    if ~isempty(r) && ~isempty(r.tau) && ~all(isnan(r.tau))
        off=0; if isfield(r,'offset') && ~isempty(r.offset) && ~isnan(r.offset), off=r.offset; end
        Yhat(:,f)=off + sum(r.amp(:)'.*exp(-xcommon./r.tau(:)'),2);
    end
end
if any(any(~isnan(Yhat(:,foi))))
    plot(log10(xcommon), mean(Yhat(:,foi),2,'omitnan')+1, 'k--', 'LineWidth', 2);
end
set(gca,'XTick', [0 1 2 3 4],'XTickLabel',10.^[0 1 2 3 4],'Box','off');
xlabel('Time since previous spike (ms)'); ylabel('Normalized bAP apical AUC');
ylim([0.7 1.8]); xlim([0.5 3.6]);
set(gca,'xdir','reverse','YAxisLocation','right')

% pairwise significance vs the pooled long-ISI baseline (bins 7:8), stars over bins 1:6
meanM = mean(M(foi,:),1,'omitnan');
pooledRef = M(foi,7:8); pooledRef = pooledRef(:);   % baseline group (2 x nNeuron values)
p_bin = nan(1,6);
for b=1:6
    pp = get_pValue({pooledRef, M(foi,b)}, 0);      % unpaired (different group sizes)
    p_bin(b) = pp(1,2);
end
yl=ylim; yoff=0.03*diff(yl);
for b=1:4
    if p_bin(b)<0.05
        if p_bin(b)<0.001, star='***'; elseif p_bin(b)<0.01, star='**'; else, star='*'; end
        text(t_center(b), meanM(b)+Sall(b)+yoff, star, 'Color','k','FontWeight','bold', ...
             'FontSize',12,'HorizontalAlignment','center','VerticalAlignment','bottom');
    end
end

fprintf('Decay cont; mean: %.2f, std: %.2f\n',mean(bAPdecayConst(foi)),std(bAPdecayConst(foi)))
set_font('Arial'); set_fontsize(16);
set_figsize(150,110);

%% ISI vs TX ratio (Figure 4)
fprintf('Section: Figure 4 - ISI vs TX ratio\n');

figure(11); clf; tiledlayout(2,2);
maxspikeN=max(LabelMat(:,1)); cmap=hsv(maxspikeN);
nexttile([1 1]);
for Spikeorder2show=[1:maxspikeN]
    showInd=find(LabelMat(:,5)==0 & LabelMat(:,2)==0 & LabelMat(:,1)==Spikeorder2show & emptycell'==0);
    ShowDat=cell2mat(cellfun(@(x) x([2 4],2),AUCbAPbin_cell(showInd),'UniformOutput',false));
    scatter(ShowDat(1,:),ShowDat(2,:),20,cmap(Spikeorder2show,:),'filled'); hold all
end
xlabel('Soma AUC'); ylabel('Apical AUC');
legend(counting_string(1:maxspikeN),'Location','northwest')

cmapst=[0 0 0; 1 0 0; 0 0.5 1];
nexttile([1 1]);
for stype=[3 1 2]
    showInd=find(LabelMat(:,5)==0 & LabelMat(:,2)==0 & SpType==stype & LabelMat(:,1)==1 & emptycell'==0);
    ShowDat=cell2mat(cellfun(@(x) x([2 4],2),AUCbAPbin_cell(showInd),'UniformOutput',false));
    scatter(ShowDat(1,:),ShowDat(2,:),20,cmapst(stype,:),'filled'); hold all
end
xlabel('Soma AUC'); ylabel('Apical AUC');
legend({'BS 1st','SS','CS 1st'},'Location','northwest')

cmapst=[0 0 0; 1 0 0; 0 0.5 1];
nexttile([1 1]);
for stype=[3 2]
    showInd=find(LabelMat(:,5)==0 & LabelMat(:,2)==0 & SpType==stype & LabelMat(:,1)~=1 & emptycell'==0);
    ShowDat=cell2mat(cellfun(@(x) x([2 4],2),AUCbAPbin_cell(showInd),'UniformOutput',false));
    scatter(ShowDat(1,:),ShowDat(2,:),20,cmapst(stype,:),'filled'); hold all
end
xlabel('Soma AUC'); ylabel('Apical AUC');
legend({'BS','CS'},'Location','northwest')

nexttile([1 1]);
showInd=find(LabelMat(:,5)==0 & LabelMat(:,2)==0 & emptycell'==0);
ShowDat=cell2mat(cellfun(@(x) x([4],2),AUCbAPbin_cell(showInd),'UniformOutput',false));
dsc=densityScatterChart(log10(ISIvec(LabelMat(showInd,6))),ShowDat(1,:));
dsc.ColorbarLabel = 'Probability';
[~, ax2, scat]=unmanage(dsc); scat.SizeData=40;
set(ax2, 'XTick', [0 1 2 3 4],'XTickLabel',10.^[0 1 2 3 4],'Box','off');
xlabel('Time since previous spike (ms)'); ylabel('Apical AUC');

%Eigenspace of bAP amplitude
validInd=find(LabelMat(:,5)==0 & LabelMat(:,2)==0 & sum(isnan(AUCall),1)'==0);
eigAmp=NaN(nROI,size(AUCall,2));
[V D E]=get_eigvector(AUCall(:,validInd));
eigAmp(:,validInd)=E';

figure(13); clf; tiledlayout(5,3); ax1=[];
for i=[1:5]
    showsp_low=find(eigAmp(i,:)<prctile(eigAmp(i,validInd),10));
    showsp_hi=find(eigAmp(i,:)>prctile(eigAmp(i,validInd),90));

    nexttile([1 1])
    showScaleImage(PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist)>0,V(:,i),colormap(turbo));
    title(['PC #' num2str(i) ', fraction: ' num2str(D(i)/sum(D),2)])

    ax1=[ax1 nexttile([1 1])];
    imagesc(nTau,[1:nROI],mean(AlignedbAPall{1}(:,:,showsp_low,1),3,'omitnan'),cax)
    set(gca,'YTick',[1 find(dax==0) length(dax)],'YTickLabel',num2str([min(dax) 0 max((dax))]',3))
    title(['PC #' num2str(i) ', low 10% spike-triggered average'])

    ax1=[ax1 nexttile([1 1])];
    imagesc(nTau,[1:nROI],mean(AlignedbAPall{1}(:,:,showsp_hi,1),3,'omitnan'),cax)
    set(gca,'YTick',[1 find(dax==0) length(dax)],'YTickLabel',num2str([min(dax) 0 max((dax))]',3))
    title(['PC #' num2str(i) ', top 10% spike-triggered average'])
end
linkaxes(ax1)
xlim([-50 30])

figure(14); clf; tiledlayout(2,2);
nexttile([1 1])
for Spikeorder2show=1:maxspikeN
    showInd=find(LabelMat(:,5)==0 & LabelMat(:,2)==0 & LabelMat(:,1)==Spikeorder2show & emptycell'==0);
    scatter(eigAmp(1,showInd),eigAmp(2,showInd),10,cmap(LabelMat(showInd,1),:),'filled'); hold all;
end
xlabel('PC 1'); ylabel('PC 2');
legend(counting_string(1:maxspikeN),'Location','northeast')

nexttile([1 1]); cmapst=[0 0 0; 1 0 0; 0 0.5 1];
for stype=1:3
    showInd=find(LabelMat(:,5)==0 & LabelMat(:,2)==0 & SpType==stype & LabelMat(:,1)==1 & emptycell'==0);
    scatter(eigAmp(1,showInd),eigAmp(2,showInd),10,cmapst(stype,:),'filled'); hold all;
end
xlabel('PC 1'); ylabel('PC 2');
legend({'SS','CS 1st','BS 1st'},'Location','northwest')

nexttile([1 1]); cmapst=[0 0 0; 1 0 0; 0 0.5 1];
for stype=[3 2]
    showInd=find(LabelMat(:,5)==0 & LabelMat(:,2)==0 & SpType==stype & LabelMat(:,1)==2 & emptycell'==0);
    scatter(eigAmp(1,showInd),eigAmp(2,showInd),10,cmapst(stype,:),'filled'); hold all;
end
xlabel('PC 1'); ylabel('PC 2');
legend({'BS 2nd','CS 2nd'},'Location','northwest')

nexttile([1 1]); cmapst=[0 0 0; 1 0 0; 0 0.5 1];
for stype=[3 2]
    showInd=find(LabelMat(:,5)==0 & LabelMat(:,2)==0 & SpType==stype & LabelMat(:,1)~=1 & emptycell'==0);
    scatter(eigAmp(1,showInd),eigAmp(2,showInd),10,cmapst(stype,:),'filled'); hold all;
end
xlabel('PC 1'); ylabel('PC 2');
legend({'BS','CS'},'Location','northwest')

%Subthreshold vs bAP amplitude

validInd=find(LabelMat(:,5)==0 & LabelMat(:,2)==0 & LabelMat(:,1)==1 & sum(isnan(AUCall),1)'==0);
Subz_cell=cell(1,size(Suball,2));
Sub_cell = num2cell(Suball(:, validInd), 1);
Sub_cell = cellfun(@(x) [dax' x(:,1)],Sub_cell,'UniformOutput',false);
[~, ztmp] = zscore_binning(Sub_cell, x_bin_edges);
Subz_cell(validInd)=ztmp;

figure(15); clf; tiledlayout(4,4);
showInd=find(LabelMat(:,5)==0 & LabelMat(:,2)==0 & LabelMat(:,1)==1 & emptycell'==0);
for dxi=1:4
    for dxj=dxi:4
        nexttile(4*(dxi-1)+dxj,[1 1]);
        ShowDatX=(cellfun(@(x) x(dxi,2),Subz_cell(showInd)));
        ShowDatY=(cellfun(@(x) x(dxj,2),AUCbAPbin_cell(showInd)));
        [R P]=corr(ShowDatX',ShowDatY');
        dsc=densityScatterChart(ShowDatX,ShowDatY);
        title_string = sprintf('Sub: from %d %s m to %d %s m, \n bAP: from %d %s m to %d %s m, \n R= %.3f, P= %.3f', ...
            x_bin_edges(dxi), char(181), x_bin_edges(dxi+1), char(181), ...
            x_bin_edges(dxj), char(181), x_bin_edges(dxj+1), char(181), R, P);
        title(title_string);
        xlabel('Pre-spike subthreshold (Z score)')
        ylabel('AUC')
    end
end

figure(16); clf; tiledlayout(2,1); cmap=hsv(5);
nexttile([1 1]);
caxAUC=[prctile(AUCall(:),5) prctile(AUCall(:),97)];
showInd=find(LabelMat(:,5)==0 & LabelMat(:,2)==0 & emptycell'==0);
[~, spEigsort]=sort(eigAmp(1,showInd));
imagesc([1:length(showInd)],[1:nROI],AUCall(:,showInd(spEigsort)),caxAUC)
xlabel('Spike ID'); ylabel('Distance from soma (\mum)');
set(gca,'YTick',[1 find(dax==0) length(dax)],'YTickLabel',num2str([min(dax) 0 max((dax))]',3))

nexttile([1 1]);
for Spikeorder2show=1:5
    showInd=find(LabelMat(:,5)==0 & LabelMat(:,2)==0 & LabelMat(:,1)==Spikeorder2show & emptycell'==0);
    ShowDat=cell2mat(cellfun(@(x) x(:,2),AUCbAPbin_cell2(showInd),'UniformOutput',false));
    ShowDat(:,sum(ShowDat<-1 | ShowDat>3,1)>0)=[];
    plot(dcenter2,ShowDat,'color',cmap(Spikeorder2show,:)); hold all
end
xlabel('Distance from soma (\mum)'); ylabel('Normalized AUC'); axis tight;
legend(counting_string(1:maxspikeN),'Location','northwest')

%% dSpikes
fprintf('Section: dSpikes\n');
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
    % conduction-speed threshold separating back-propagating APs from
    % dendritically-initiated spikes (see classify_dAP_candidates at end of file)
    [dAPCandidates_threshold, dt2show, dst2show, bAPthreshold] = classify_dAP_candidates(dSpikePropsMat{f});

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
fprintf('Section: Figure S_dAP - dAP neuron image\n');
f=22;
[dAPCandidates_threshold, dt2show, dst2show, bAPthreshold] = classify_dAP_candidates(dSpikePropsMat{f});

dAPSTA=get_STA(PCresult{f}.NormalizedTrace_dirt,FiringRate{f}.DendIntSpFrame,10,20);
somAPSTA=get_STA(PCresult{f}.NormalizedTrace_dirt,setdiff(FiringRate{f}.SimpleSpikeFrame,FiringRate{f}.DendIntSpFrame),10,20);
nROI=size(PCresult{f}.NormalizedTrace_dirt,1);
dAPSTA=dAPSTA-median(dAPSTA(:,1:5),2);
somAPSTA=somAPSTA-median(somAPSTA(:,1:5),2);

figure(33); clf; tiledlayout(3,2);
ax1=nexttile([2 1]); cax=[-0.02 0.2];
imshow2(PCresult{f}.avgImg,[]);
drawScaleBar(100/PCresult{f}.pixelsize,'horizontal','color',[1 1 1]);
nexttile([2 1]);
scatter(dt2show,dst2show,20,[0 0 0],'filled'); hold all
f2=plot([-1:2],bAPthreshold([-1:2]),'r--','linewidth',2);
xlabel('\Deltat (ms)'); ylabel('Distance from soma (µm)'); axis tight;
ax3=nexttile([1 1]);
imagesc([-10:20],[1:nROI],dAPSTA,cax);
set_kymoYtick(PCresult{f}.dendaxis)
xlabel('Time (ms)'); ylabel('Distance (µm)'); axis tight;
title(sprintf('Average of spikes \ninititated from dendrite (dAPs)'))
ax2=nexttile([1 1]);
imagesc([-10:20],[1:nROI],somAPSTA,cax);
set_kymoYtick(PCresult{f}.dendaxis)
title(sprintf('Average of spikes \ninititated from soma (SomAPs)'))
xlabel('Time (ms)'); ylabel('Distance (µm)'); axis tight;
colormap(ax2,turbo);
colormap(ax3,turbo);
cb=colorbar;
cb.Label.String='∆F/F';
set_fontsize(16); set_font('Arial');
set_figsize(310,200)
%% Show neuron image and example (figure 3)
fprintf('Section: Show neuron image and example (figure 3)\n');
f=18;
%load(fullfile(fpath{f}, 'PC_Result.mat'), 'Result')
omitROI=[3]; scaleoffset=0.35; cmap_spclass=[0 0 0; 1 0.57 0.11; 1 0.57 0.11; 0.9 0.46 0.1];
[nROI nTime]=size(PCresult{f}.NormalizedTrace_dirt);
ShowClassMat=[ind2vec(nTime, FiringRate{f}.SomaSpikeFrame,1); ind2vec(nTime, FiringRate{f}.DendIntSpFrame,1); ind2vec(nTime, FiringRate{f}.DendOnlySpFrame,1)];
ShowROI=setdiff([1:nROI],omitROI);  subShow_ROI=[2 4 6 8 13 15 16]; %subShow_ROI=[4 6 12 13 14 17 18]; 
sz=size(Result.ref_im);
cmap_ROIs=gen_colormap(Plasma,length(ShowROI));
ftprnts=PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist);
str=PCresult{f}.swc;

figure(11); clf;
nexttile([2 4]);
l=plot(PCresult{f}.NormalizedTrace_dirt(ShowROI,:)'-scaleoffset*[1:length(ShowROI)]);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap_ROIs,2)); axis off; hold all;
[spc spt]=find(ShowClassMat); g=1;
for spclass=1:4
    if ~isempty(find(spc==spclass))
        plot([spt(spc==spclass) spt(spc==spclass)]',(g-1)*scaleoffset+[0 scaleoffset*2/3]','color',cmap_spclass(spclass,:),'linewidth',1)
        g=g+1;
    end
end
plot(rescale(PCresult{f}.VR(5,:))*scaleoffset+(g-1)*scaleoffset,'color',[0 0 0]);
drawScaleBar(60000,'horizontal','color','k','position',[1.01*nTime -nROI*scaleoffset],'linewidth',2);
drawScaleBar(0.2,'vertical','color','k','position',[1.01*nTime -nROI*scaleoffset],'linewidth',2)
drawScaleBar(scaleoffset,'vertical','color','k','position',[1.01*nTime scaleoffset*2],'linewidth',2)
text(ones(length(ShowROI),1)*(-5000),-scaleoffset*[1:length(ShowROI)],...
    num2str(PCresult{f}.dendaxis(ShowROI)','%2.0f µm'),'fontsize',12,'horizontalalignment','right','fontname','Arial');
text(ones(1,1)*(-5000),scaleoffset*[0]+scaleoffset*0.25,...
    {'SomAP'},'fontsize',14,'horizontalalignment','right');
text(ones(1,1)*(-5000),scaleoffset*[1]+scaleoffset*0.25,...
    {'dSpike'},'fontsize',14,'horizontalalignment','right','color',cmap_spclass(2,:));
text(ones(1,1)*(-5000),scaleoffset*[2]+scaleoffset*0.25,...
    {'VR'},'fontsize',14,'horizontalalignment','right');
set_font('Arial');
set_figsize(170,220);

figure(12); clf; tiledlayout(2,1,'padding','tight','TileSpacing','tight'); frame2show=[105426 401970]; length2show=1300; ax1=[];
scaleoffset=0.22;
for fr=frame2show
    ax1=[ax1 nexttile([1 1])];
    x_show=[(fr):fr+length2show];
    tr2show=PCresult{f}.NormalizedTrace_dirt(ShowROI,x_show);
    tr2show=pcafilterTrace(tr2show,[1:5]);
    l=plot(tr2show(subShow_ROI,:)'-scaleoffset*[1:length(subShow_ROI)],'linewidth',1);
    arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap_ROIs(subShow_ROI,:),2)); axis off; hold all;
    [spc spt]=find(ShowClassMat(:,x_show)); g=1;
    for spclass=1:3
        if ~isempty(find(spc==spclass))
            scatter([spt(spc==spclass)],repmat((g-1)*scaleoffset/3,1,sum(spc==spclass)),40,cmap_spclass(spclass,:),'linewidth',2,'marker','|')
            g=g+1;
        end
    end
    text(ones(length(subShow_ROI),1)*(-50),-scaleoffset*[1:length(subShow_ROI)]+scaleoffset*0.5,...
        num2str(PCresult{f}.dendaxis(ShowROI(subShow_ROI))','%2.0f µm'),'fontname','Arial','fontsize',11)
    text(ones(1,1)*(35),scaleoffset*[0],...
        {'SomAP'},'fontsize',11,'horizontalalignment','right');
    text(ones(1,1)*(35),scaleoffset*[0.4],...
        {'dSpike'},'fontsize',11,'horizontalalignment','right','color',cmap_spclass(2,:));
end
linkaxes(ax1,'xy');
ylim([-1.7 0.1]);
drawScaleBar(100,'horizontal','color','k','position',[1310 0]);
drawScaleBar(0.1,'vertical','color','k','position',[1310 0])
xlim([0 1320])
set_font('Arial'); set_fontsize(12);
set_figsize(220,190);
%%
figure(19); clf; %show footprint
theta = 167.3;  % degrees
Rmat = [cosd(theta) -sind(theta) 0;  sind(theta)  cosd(theta) 0; 0 0 1];
str_rot=[[transformPointsForward(affine2d(Rmat),str(:,[1 2]))] str(:,3)];
scatter(str_rot(:,2),str_rot(:,1),str_rot(:,3),[0.6 0.6 0.6],'filled'); axis equal tight off; hold all;

for r=1:length(ShowROI)
    DMDmask=ftprnts(:,:,ShowROI(r))>0;
    str_ROI=str(find(str(:,1)>0.5 & str(:,1)<(sz(2)-0.5) & str(:,2)>0.5 & str(:,2)<(sz(1)-0.5)),:);
    str_ROI=str_ROI(DMDmask(sub2ind(sz,round(str_ROI(:,2)),round(str_ROI(:,1)))),:);
    str_ROI_rot=[[transformPointsForward(affine2d(Rmat),str_ROI(:,[1 2]))] str_ROI(:,3)];
    scatter(str_ROI_rot(:,2),str_ROI_rot(:,1),str_ROI_rot(:,3),cmap_ROIs(r,:),'filled');
end
%set(gca,'ydir','reverse','xdir','reverse')
bar100um=100/PCresult{f}.pixelsize;
plot(prctile(str_rot(:,2),25)-[0 bar100um],prctile(str_rot(:,1),10)*[1 1],'color',[0 0 0],'LineWidth',3)

%% AUC scatter plot
figure(23); clf;
dclass2show=[3 4 5];
nexttile([1 1]);
SomAP2show = find(ismember(dSpikePropsMat{f}.Spike_frame,FiringRate{f}.SomaSpikeFrame) & dSpikePropsMat{f}.Branch_Index==1 & dSpikePropsMat{f}.SpikeOrder==1);
dAP2show = find(ismember(dSpikePropsMat{f}.Spike_frame,FiringRate{f}.DendIntSpFrame) & ismember(dSpikePropsMat{f}.dClass,dclass2show));
dSpike2show = find(ismember(dSpikePropsMat{f}.Spike_frame,FiringRate{f}.DendOnlySpFrame) & ismember(dSpikePropsMat{f}.dClass,dclass2show));
dSpike2show_B = find(ismember(dSpikePropsMat{f}.Spike_frame,FiringRate{f}.dSpikeBasal));
dSpike2show_A = find(ismember(dSpikePropsMat{f}.Spike_frame,FiringRate{f}.dSpikeApical));

inBurst=find(PCresult{f}.CStrace(dSpikePropsMat{f}.Spike_frame) | PCresult{f}.BStrace(dSpikePropsMat{f}.Spike_frame));
inBurst=inBurst(ismember(inBurst,[SomAP2show; dAP2show; dSpike2show]));
dSpike2show2=dSpike2show;
dSpike2show2(ismember(dSpike2show,inBurst))=[];
dSpike2show_B(ismember(dSpike2show_B,inBurst))=[];
dSpike2show_A(ismember(dSpike2show_A,inBurst))=[];

scatter(dSpikePropsMat{f}.AUC_soma(SomAP2show),dSpikePropsMat{f}.AUC_apical(SomAP2show),25,'filled'); hold all
scatter(dSpikePropsMat{f}.AUC_soma(dAP2show),dSpikePropsMat{f}.AUC_apical(dAP2show),35,'filled');
scatter(dSpikePropsMat{f}.AUC_soma(dSpike2show),dSpikePropsMat{f}.AUC_apical(dSpike2show),25,'filled');
%scatter(dSpikePropsMat{f}(inBurst,10),dSpikePropsMat{f}(inBurst,11),35,'co');
xlabel('AUC soma'); ylabel('AUC apical');
title(['Neuron # ' num2str(f)])
legend({'SomAP','dAP','dSpike'});
%% dSpike Kymograph
fprintf('Section: dSpike kymograph\n');
%ROIvec=Result.maintrunk;
%ROIvec_distorder=find_index_bh(PCresult{f}.sortDist,ROIvec);
ROIvec_distorder=ShowROI;
dax=PCresult{f}.dendaxis(ROIvec_distorder);
SomAPkymo=mean(AligndSpiketr{f}(:,:,SomAP2show),3,'omitnan');
SomAPkymo_main=SomAPkymo(ROIvec_distorder,:);
dAPkymo=mean(AligndSpiketr{f}(:,:,dAP2show),3,'omitnan');
dSpikekymo=mean(AligndSpiketr{f}(:,:,dSpike2show2),3,'omitnan');
dSpikekymo_main=dSpikekymo(ROIvec_distorder,:);
dSpikeBasal_kymo=mean(AligndSpiketr{f}(:,:,dSpike2show_B),3,'omitnan');
dSpikeApical_kymo=mean(AligndSpiketr{f}(:,:,dSpike2show_A),3,'omitnan');
showKymoCell={SomAPkymo,dAPkymo,dSpikeApical_kymo};
showKymoCell_main={SomAPkymo_main,dAPkymo,dSpikekymo_main};

% coordinate rotation replacing view(-155,90) so the vertical scalebar is not rotated
Rv = viewmtx(-155,90); Rv = Rv(1:2,1:2);
cax=[-0.02 0.13];
SWCpoints=PCresult{f}.swc;
[dax_sorted, sortind]=sort(dax,'ascend');

% SomAP (k=1) and dSpike (k=3) are now split into separate figures
kList=[1 3]; figNums=[24 25]; titleStr={'SomAP','dSpike'};
for gi=1:numel(kList)
    k=kList(gi);
    showKymoCell{k}=showKymoCell{k}-mean(showKymoCell{k}(:,[1:10]),2);

    if gi==1
        figure(figNums(gi)); clf; tiledlayout(2,6,'padding','compact','TileSpacing','tight'); 
        sgtitle(titleStr{gi},'FontSize',23);
    else
        figure(figNums(gi)); clf; tiledlayout(2,9,'padding','compact','TileSpacing','tight'); sgtitle(titleStr{gi},'FontSize',23);
        axm2=nexttile(4,[2 3]);
        dSpikeRate=sum(PCresult{f}.SpikeMat(PCresult{f}.sortDist,:).*ind2vec(nTime,unique(FiringRate{f}.DendOnlySpFrame+[-2:2]),1),2,'omitnan')./FiringRate{f}.ValidPeriod;
        showScaleScatter_rot(dSpikeRate*1000, SWCpoints, PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist),'copper',[0 max(dSpikeRate*1000)],Rv);
        cb=colorbar;
        colormap(axm2,'copper');
        cb.Label.String = 'dSpike rate (Hz)';
        cb.Location="southoutside";
    end

    % morphology with rotated coordinates (vertical scalebar, no view())
    axm=nexttile(1,[2 3]);
    ShowAmp=max(showKymoCell{k}(:,nTauAlign(1)+[-2:3]),[],2);
    showScaleScatter_rot(ShowAmp, SWCpoints, PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist),'turbo',cax,Rv);
    drawScaleBar(100/Pixelsize(f),'vertical','color','k','position',[-100 -330]); axis tight;
    cbm=colorbar(axm,'southoutside'); cbm.Label.String='∆F/F';

    % kymograph
    ShowKymo_main=showKymoCell_main{k}(sortind,:);
    ShowKymo_main=ShowKymo_main-mean(ShowKymo_main(:,1:5),2);
    axk=nexttile([1 3]);
    imagesc([-nTauAlign(1):nTauAlign(2)],[1:size(showKymoCell_main{k},1)],ShowKymo_main,cax);
    [yt yl]=set_kymoYtick(dax_sorted);
    ylabel('Distance (µm)');
    colormap(axk,turbo);
    cb1=colorbar(axk); cb1.Label.String='∆F/F';

    % per-ROI traces
    axt=nexttile([1 3]);
    l=plot([-nTauAlign(1):nTauAlign(2)],ShowKymo_main,'LineWidth',1.5);
    arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap_ROIs,2));
    xlabel('Time (ms)'); ylabel('∆F/F'); box off; ylim([-0.03 0.14]);
    colormap(axt,cmap_ROIs);
    cb2=colorbar(axt); cb2.Ticks=yt/max(yt); cb2.TickLabels=yl;
    cb2.Label.String='Distance from soma (µm)';
    cb2.Label.FontName='Arial';

    linkaxes([axk axt],'x'); 
    if gi==1
    xlim(axt,[-10 5]);    
    else
    xlim(axt,[-12 10]);
    end
    set_font('Arial'); set_fontsize(18);
    if gi==1
    set_figsize(200,170);
    else
    colormap(axm2,'copper');    
    set_figsize(270,170);
    end
end

%% Plot SomAP vs dAP vs dSpike Rate
fprintf('Section: SomAP vs dAP vs dSpike rate\n');
figure(22); clf; tiledlayout(4,2); cmap_FR=nebula(5); cmap_region=[0 0.4 1; 0.4 1 0; 1 0 0];
SomAP_FR=cellfun(@(x) x.SomAP,FiringRate(foi));
dAP_FR=cellfun(@(x) x.dAP,FiringRate(foi));
CS_FR=cellfun(@(x) x.CS,FiringRate(foi));
SS_FR=cellfun(@(x) x.SS,FiringRate(foi));
dSpike_FR=cellfun(@(x) x.dSpike,FiringRate(foi));
Ratio_FR=[SomAP_FR; dAP_FR; dSpike_FR]./sum([SomAP_FR; dAP_FR; dSpike_FR],1);
dAPClass_FR=cell2mat(cellfun(@(x) x.dAPdClass*x.dAP,FiringRate(foi),'UniformOutput',false)');
Ratio_dClass=[SomAP_FR' dAPClass_FR]./(SomAP_FR+dAP_FR)';
Ratio_dSpdClass=cell2mat(cellfun(@(x) x.dSpdClass,FiringRate(foi),'UniformOutput',false)');

nexttile([2 1])
p=Boxplot_wPoints2([SomAP_FR; dAP_FR; dSpike_FR]',cmap_FR([3 2 4],:));
set(gca,'xtick',[1:4],'xticklabel',{'SomAP','dAP','dSpike'});
drawPValueLines(p,0,'TextYOffset',0.5,'StepHeight',1)
ylabel('Firing rate (Hz)'); box off;

nexttile([2 1])
p=Boxplot_wPoints2([SS_FR; CS_FR]',cmap_FR([1 5],:));
set(gca,'xtick',[1:4],'xticklabel',{'SS','CS'});
drawPValueLines(p,0,'TextYOffset',0.5,'StepHeight',1)
ylabel('Firing rate (Hz)'); box off;

nexttile([1 1])
p=piechart(mean(Ratio_FR,2),{'SomAP','dAP','dSpike'});
p.ExplodedWedges=[2 3];
p.ColorOrder=cmap_FR([3 2 4],:);
p.EdgeColor=[1 1 1];
nexttile([1 1])
d = donutchart(mean([Ratio_FR([1 2],:); Ratio_FR(3,:).*Ratio_dSpdClass'],2),{'','','Basal','Apical'});
d.ColorOrder(1:4,:)=[1 1 1; 1 1 1; cmap_FR(4,:)/4; cmap_FR(4,:)/4*3];
d.EdgeColor=[1 1 1];
d.InnerRadius = 0.8;

% Plot ratio of initiation region
figure(24); clf;
nexttile([1 1])
p=Boxplot_wPoints2(Ratio_dClass,cmap_region([2 1 3],:));
set(gca,'xtick',[1 2 3],'xticklabel',{'Soma','Basal','Apical'});
drawPValueLines(p,0.05,'TextYOffset',0.05,'StepHeight',0.1)
xlim([0.5 3.5]); ylabel('Fraction of AP'); xlabel('Initiation region'); box off;

nexttile([1 1])
p=Boxplot_wPoints2(Ratio_dSpdClass,cmap_region([1 3],:));
set(gca,'xtick',[1 2],'xticklabel',{'Basal','Apical'});
drawPValueLines(p,0,'TextYOffset',0.05,'StepHeight',0.1)
xlim([0.5 2.5]); ylabel('Fraction of dSpike'); xlabel('Initiation region'); box off;

%% Spikes on the VR position
fprintf('Section: Spikes on the VR position\n');
figure(27); clf;
tiledlayout('flow','Padding','compact');
for f=foi
    StayTime=histcounts(PCresult{f}.VR,[0:5:120]);
    CShist=histcounts(PCresult{f}.VR(5,FiringRate{f}.ComplexSpikeFrame),[0:5:120]);
    SShist=histcounts(PCresult{f}.VR(5,FiringRate{f}.SimpleSpikeFrame),[0:5:120]);
    dSPhist=histcounts(PCresult{f}.VR(5,FiringRate{f}.DendOnlySpFrame),[0:5:120]);
    apicalROI=find(ismember(PCresult{f}.branch_dClass,[3 4 5]));
    PeakFrame=find(max(PCresult{f}.peakvec(apicalROI,:),[],1)>0);
    TroughFrame=find(max(PCresult{f}.troughvec(apicalROI,:),[],1)>0);
    dSPhist=histcounts(PCresult{f}.VR(5,FiringRate{f}.DendOnlySpFrame),[0:5:120]);
    Peakhist=histcounts(PCresult{f}.VR(5,PeakFrame),[0:5:120]);
    Troughhist=histcounts(PCresult{f}.VR(5,TroughFrame),[0:5:120]);
    nexttile([1 1]);
    plot([2.5:5:117.5],rescale2([CShist; SShist; dSPhist; Peakhist; Troughhist]./StayTime,2)');
    nexttile([1 1]);
    dSpike2show_A = find(ismember(dSpikePropsMat{f}.Spike_frame,FiringRate{f}.dSpikeApical));
    scatter(PCresult{f}.VR(5,dSpikePropsMat{f}.Spike_frame(dSpike2show_A))',PCresult{f}.VR(8,dSpikePropsMat{f}.Spike_frame(dSpike2show_A))',5,'filled');  hold all;
    scatter(PCresult{f}.VR(5,PeakFrame)',PCresult{f}.VR(8,PeakFrame)',5,'filled');
    scatter(PCresult{f}.VR(5,TroughFrame)',PCresult{f}.VR(8,TroughFrame)',5,'filled');
    title(f)
end
legend({'CS','SS','dSP','Peak','Trough','dSP','Peak','Trough'})

%% dSpike-triggered average of other event types and behavioral state
fprintf('Section: dSpike-triggered average (event types & state)\n');
% Align to every dSpike (DendOnlySpFrame) and average the occurrence of other event
% types (SS, BS, CS, dAP, apical peaks) and the behavioral-state occupancy
% (running/resting, in/out place field) as a function of lag from the dSpike.
% This shows directly whether a dSpike tends to precede/follow each spike type,
% and in which state dSpikes occur.
Twin=[300 300]; lags=-Twin(1):Twin(2);   % ms before/after (1 kHz -> 1 ms/frame)
smthW=5;                                  % ms smoothing for display
Pwin=50;                                 % ms window for before/after event probability
evNames={'SS','BS','CS','dAP','Peak'};
stateNames={'Running','Resting','In PF','Out PF'};
sumEv  = zeros(numel(evNames),numel(lags));
sumRun = zeros(2,numel(lags));  nDsp_all=0;   % running/resting (all neurons)
sumPF  = zeros(2,numel(lags));  nDsp_PF =0;   % in/out place field (PF neurons only)
% P(>=1 event of each type within Pwin ms) before / after each dSpike
cntBefore=zeros(1,numel(evNames)); cntAfter=zeros(1,numel(evNames));
colAfter  = Twin(1)+1 + (1:Pwin);      % lag columns  +1 .. +Pwin
colBefore = Twin(1)+1 - (Pwin:-1:1);   % lag columns  -Pwin .. -1
% shuffle null: circularly shift the dSpike times nShuffle times
nShuffle=1000; rng(0);                  % rng seed for reproducible shuffles
shufBefore=zeros(numel(evNames),nShuffle); shufAfter=zeros(numel(evNames),nShuffle);
nNeuron_dSp=0;
foi2analyze=foi;
for i=1:length(foi2analyze)
    f=foi2analyze(i);
    nTime=size(PCresult{f}.NormalizedTrace_dirt,2);
    somaAP=find(max(PCresult{f}.SpikeMat(PCresult{f}.BranchLabel==1,:),[],1)>0);
    apicalROI=find(ismember(PCresult{f}.branch_dClass,[3 4 5]));

    % event trains (1 x nTime binary)
    evTrains=[ ind2vec(nTime,FiringRate{f}.SimpleSpikeFrame,1); ...             % SS
               ind2vec(nTime,somaAP(PCresult{f}.BStrace(somaAP)>0),1); ...      % BS (somatic AP in burst)
               ind2vec(nTime,FiringRate{f}.ComplexSpikeFrame,1); ...           % CS
               ind2vec(nTime,FiringRate{f}.DendIntSpFrame,1); ...              % dAP
               ind2vec(nTime,find(max(PCresult{f}.peakvec(apicalROI,:),[],1)>0),1) ]; % apical peaks

    dsp=FiringRate{f}.DendOnlySpFrame(:)';
    dsp=dsp(dsp>Twin(1) & dsp<=nTime-Twin(2));  % keep only windows fully in range
    if isempty(dsp); continue; end
    W=dsp'+lags;                                 % nDsp x nLag frame indices
    for e=1:numel(evNames)
        tr=evTrains(e,:);                        % row vector; single-subscript keeps size(W)
        M=tr(W);                                 % nDsp x nLag
        sumEv(e,:)=sumEv(e,:)+sum(M,1);
        cntAfter(e) =cntAfter(e) +sum(any(M(:,colAfter )>0,2)); % dSpikes with event after
        cntBefore(e)=cntBefore(e)+sum(any(M(:,colBefore)>0,2)); % dSpikes with event before
    end

    % --- circular-shift shuffle null for the before/after probability ---
    baseBA = dsp' + [(-Pwin:-1) (1:Pwin)];        % nDsp x 2*Pwin window frame indices
    for sh=1:nShuffle
        s=randi(nTime);                            % random circular shift of dSpike times
        idxS = mod(baseBA + s - 1, nTime) + 1;     % shifted window frames (wrap around ends)
        for e=1:numel(evNames)
            tr=evTrains(e,:);
            Msh=tr(idxS);                          % nDsp x 2*Pwin
            shufBefore(e,sh)=shufBefore(e,sh)+sum(any(Msh(:,1:Pwin)>0,2));
            shufAfter(e,sh) =shufAfter(e,sh) +sum(any(Msh(:,Pwin+1:end)>0,2));
        end
    end

    runV=PCresult{f}.runVec>0;
    sumRun(1,:)=sumRun(1,:)+sum(runV(W),1);
    sumRun(2,:)=sumRun(2,:)+sum(~runV(W),1);
    nDsp_all=nDsp_all+length(dsp);
    nNeuron_dSp=nNeuron_dSp+1;
    if isfield(PCresult{f},'PFvec')
        pfV=PCresult{f}.PFvec>0;
        sumPF(1,:)=sumPF(1,:)+sum(pfV(W),1);
        sumPF(2,:)=sumPF(2,:)+sum(~pfV(W),1);
        nDsp_PF=nDsp_PF+length(dsp);
    end
end
% event types -> instantaneous rate (Hz); states -> occupancy fraction
evRate  = movmean(sumEv/nDsp_all*1000, smthW, 2);
runFrac = movmean(sumRun/nDsp_all,      smthW, 2);
pfFrac  = movmean(sumPF/nDsp_PF,        smthW, 2);
% probability of each spike type within Pwin ms before / after a dSpike
probBefore=cntBefore/nDsp_all; probAfter=cntAfter/nDsp_all;
% shuffle null distributions (nEv x nShuffle) and one-sided p (obs >= chance)
probShufBefore=shufBefore/nDsp_all; probShufAfter=shufAfter/nDsp_all;
pBefore=mean(probShufBefore>=probBefore(:),2)'; pAfter=mean(probShufAfter>=probAfter(:),2)';
fprintf('dSpike-triggered average: %d dSpikes from %d neurons (%d in PF neurons)\n',nDsp_all,nNeuron_dSp,nDsp_PF);
fprintf('P(event within %d ms of dSpike), obs vs %d-shuffle null (mean), one-sided p:\n',Pwin,nShuffle);
fprintf('   type   obs[before|after]   shuffle[before|after]   p[before|after]\n');
for e=1:numel(evNames)
    fprintf('   %-5s  %.3f | %.3f       %.3f | %.3f          %.3g | %.3g\n', evNames{e}, ...
        probBefore(e),probAfter(e), mean(probShufBefore(e,:)),mean(probShufAfter(e,:)), pBefore(e),pAfter(e));
end

cmap_ev=[0 0 0; 0.2 0.5 1; 0.9 0.2 0.2; 1 0.6 0.1; 0.4 0.7 0.3];   % SS BS CS dAP Peak
cmap_st=[0 0.6 0.2; 0.6 0.6 0.6; 1 0.4 0.1; 0.3 0.3 0.3];          % Run Rest inPF outPF
figure(280); clf; tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
nexttile([1 1]);
set(gca,'ColorOrder',cmap_ev,'NextPlot','replacechildren');
plot(lags,evRate','LineWidth',1.5); hold all;
xline(0,'k--'); xlabel('Time from dSpike (ms)'); ylabel('Event rate (Hz)'); box off;
legend(evNames,'Location','northwest'); title('dSpike-triggered event rate'); axis tight;
% probability of each spike type within Pwin ms, before vs after the dSpike
nexttile([1 1]);
b=bar([probBefore; probAfter]'*100); hold all;
b(1).FaceColor=[0.55 0.55 0.55]; b(2).FaceColor=[0.85 0.33 0.10];
% overlay shuffle null: mean marker + 95% CI whiskers on each bar
mB=mean(probShufBefore,2)'*100; mA=mean(probShufAfter,2)'*100;
loB=prctile(probShufBefore,2.5,2)'*100; hiB=prctile(probShufBefore,97.5,2)'*100;
loA=prctile(probShufAfter, 2.5,2)'*100; hiA=prctile(probShufAfter, 97.5,2)'*100;
errorbar(b(1).XEndPoints,mB,mB-loB,hiB-mB,'LineStyle','none','Color','k','Marker','_','CapSize',6,'LineWidth',1);
eS=errorbar(b(2).XEndPoints,mA,mA-loA,hiA-mA,'LineStyle','none','Color','k','Marker','_','CapSize',6,'LineWidth',1);
set(gca,'xtick',1:numel(evNames),'xticklabel',evNames);
ylabel(sprintf('P(spike within %d ms) (%%)',Pwin)); box off;
legend([b(1) b(2) eS],{'before','after','shuffle 95% CI'},'Location','northwest');
title('Spike probability around dSpike');
nexttile([1 1]);
set(gca,'ColorOrder',cmap_st,'NextPlot','replacechildren');
plot(lags,[runFrac; pfFrac]','LineWidth',1.5); hold all;
xline(0,'k--'); xlabel('Time from dSpike (ms)'); ylabel('State occupancy (fraction)'); box off;
legend(stateNames,'Location','best'); title('dSpike-triggered state occupancy'); axis tight;
set_font('Arial'); set_fontsize(14);

%% Somatic AUC histogram: dSpike vs simple spike (cf. AUC scatter above)
fprintf('Section: Somatic AUC histogram (dSpike vs simple spike)\n');
dclass2show=[3 4 5];
AUCsoma_SS_all=[]; AUCsoma_dSp_all=[];
nNeuron_SS=0; nNeuron_dSp=0;
for f=foi
    SSrows  = find(ismember(dSpikePropsMat{f}.Spike_frame,FiringRate{f}.SimpleSpikeFrame) & dSpikePropsMat{f}.Branch_Index==1 & dSpikePropsMat{f}.SpikeOrder==1);
    dSprows = find(ismember(dSpikePropsMat{f}.Spike_frame,FiringRate{f}.DendOnlySpFrame) & ismember(dSpikePropsMat{f}.dClass,dclass2show));
    AUCsoma_SS_all =[AUCsoma_SS_all;  dSpikePropsMat{f}.AUC_soma(SSrows)];
    AUCsoma_dSp_all=[AUCsoma_dSp_all; dSpikePropsMat{f}.AUC_soma(dSprows)];
    nNeuron_SS =nNeuron_SS +~isempty(SSrows);
    nNeuron_dSp=nNeuron_dSp+~isempty(dSprows);
end
fprintf('AUC histogram: %d simple spikes (%d neurons), %d dSpikes (%d neurons), of %d neurons total\n', ...
    sum(~isnan(AUCsoma_SS_all)),nNeuron_SS,sum(~isnan(AUCsoma_dSp_all)),nNeuron_dSp,numel(foi));

fExample=18; % same example neuron as the AUC scatter above
SSrows_e  = find(ismember(dSpikePropsMat{fExample}.Spike_frame,FiringRate{fExample}.SimpleSpikeFrame) & dSpikePropsMat{fExample}.Branch_Index==1 & dSpikePropsMat{fExample}.SpikeOrder==1);
dSprows_e = find(ismember(dSpikePropsMat{fExample}.Spike_frame,FiringRate{fExample}.DendOnlySpFrame) & ismember(dSpikePropsMat{fExample}.dClass,dclass2show));
AUCsoma_SS_e  = dSpikePropsMat{fExample}.AUC_soma(SSrows_e);
AUCsoma_dSp_e = dSpikePropsMat{fExample}.AUC_soma(dSprows_e);

cmap_hist=[0 0.45 0.74; 0.85 0.33 0.10]; % SS / dSpike
figure(281); clf; tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
% pooled across neurons
edges_p=linspace(prctile([AUCsoma_SS_all; AUCsoma_dSp_all]/13,1),prctile([AUCsoma_SS_all; AUCsoma_dSp_all]/13,99),40);
histogram(AUCsoma_SS_all/13 ,edges_p,'Normalization','probability','FaceColor',cmap_hist(1,:),'EdgeColor','none'); hold all;
histogram(AUCsoma_dSp_all/13,edges_p,'Normalization','probability','FaceColor',cmap_hist(2,:),'EdgeColor','none');
xlabel('Mean somatic voltage (∆F/F)'); ylabel('Probability'); box off;
p_auc=ranksum(AUCsoma_SS_all(~isnan(AUCsoma_SS_all)),AUCsoma_dSp_all(~isnan(AUCsoma_dSp_all)));
legend({'Simple spike','dSpike'}); title(sprintf('All neurons (ranksum p = %.2e)',p_auc));
set_font('Arial'); set_fontsize(14);
set_figsize(150,100);

%% Spike types in place field & during running vs resting
fprintf('Section: spike types place field and locomotion state\n');
% For each spike type (SS, BS, CS, dAP, dSpike) split by locomotion state (runVec)
% and place field (PFvec):
%   - occupancy-normalized RATE (Hz) in each condition
%   - PROBABILITY that a spike of that type is detected in vs out of the place field
%     ( P = fraction of that type's spikes falling in/out PF; the two sum to 1 ).
% Neurons without a defined place field contribute NaN (skipped by the paired stats).
% column layout: [type1_A .. typeN_A  type1_B .. typeN_B]  (A/B = the two states)
spkNames={'SS','BS','CS','dAP','dSpike'};
nSpk=numel(spkNames);
runRate=NaN(length(foi),2*nSpk);   % A=running  B=resting            (Hz)
PFRate =NaN(length(foi),2*nSpk);   % A=in PF    B=out of PF          (Hz)
runProb=NaN(length(foi),2*nSpk);   % A=running  B=resting            (fraction of that type)
PFprob =NaN(length(foi),2*nSpk);   % A=in PF    B=out of PF          (fraction of that type)
occRun =NaN(length(foi),1);        % running occupancy fraction (chance level for P(run))
occPF  =NaN(length(foi),1);        % PF occupancy fraction (chance level for P(in PF))
for i=1:length(foi)
    f=foi(i);
    somaAP=find(max(PCresult{f}.SpikeMat(PCresult{f}.BranchLabel==1,:),[],1)>0);
    frames={ FiringRate{f}.SimpleSpikeFrame(:)', ...                  % SS
             somaAP(PCresult{f}.BStrace(somaAP)>0), ...               % BS (somatic AP in burst)
             FiringRate{f}.ComplexSpikeFrame(:)', ...                 % CS
             FiringRate{f}.DendIntSpFrame(:)', ...                    % dAP
             FiringRate{f}.DendOnlySpFrame(:)' };                     % dSpike

    runV=PCresult{f}.runVec>0;
    Trun=sum(runV)/1000; Trest=sum(~runV)/1000;                       % seconds (1 kHz)
    occRun(i)=Trun/(Trun+Trest);
    for s=1:nSpk
        fr=frames{s};
        runRate(i,s)      = sum(runV(fr))/Trun;
        runRate(i,s+nSpk) = sum(~runV(fr))/Trest;
        nTot=numel(fr);
        if nTot>0
            runProb(i,s)      = sum(runV(fr))/nTot;   % P(run  | spike type)
            runProb(i,s+nSpk) = sum(~runV(fr))/nTot;  % P(rest | spike type)
        end
    end

    if isfield(PCresult{f},'PFvec')
        pfV=PCresult{f}.PFvec>0;
        Tpf=sum(pfV)/1000; Tout=sum(~pfV)/1000;
        occPF(i)=Tpf/(Tpf+Tout);
        for s=1:nSpk
            fr=frames{s};
            PFRate(i,s)      = sum(pfV(fr))/Tpf;
            PFRate(i,s+nSpk) = sum(~pfV(fr))/Tout;
            nTot=numel(fr);
            if nTot>0
                PFprob(i,s)      = sum(pfV(fr))/nTot;    % P(in PF  | spike type)
                PFprob(i,s+nSpk) = sum(~pfV(fr))/nTot;   % P(out PF | spike type)
            end
        end
    end
end

% Significance vs chance: per-neuron observed fraction vs per-neuron occupancy,
% paired across neurons (Wilcoxon signed-rank). Stars are drawn only when p<0.05.
p_run=nan(1,nSpk); p_pf=nan(1,nSpk);
for s=1:nSpk
    m=~isnan(runProb(:,s)) & ~isnan(occRun);
    if sum(m)>=3; p_run(s)=signrank(runProb(m,s),occRun(m)); end
    m=~isnan(PFprob(:,s)) & ~isnan(occPF);
    if sum(m)>=3; p_pf(s)=signrank(PFprob(m,s),occPF(m)); end
end

% Figure 282: Rows 1-2 = occupancy-normalized firing rate (Hz); Row 3 = fraction of
% each spike type by state (stacked to 1), dashed = chance = time-occupancy fraction,
% stars = paired signed-rank of the observed fraction vs chance (only if p<0.05).
figure(282); clf; tiledlayout(3,nSpk,'Padding','compact','TileSpacing','compact');

% Row 1: firing rate (Hz), running vs resting  (paired t-test, labelled above)
for s=1:nSpk
    nexttile([1 1]);
    p=Boxplot_wPoints2(runRate(:,[s s+nSpk]),[0 0.6 0.2; 0.6 0.6 0.6]); hold all;
    set(gca,'xtick',[1 2],'xticklabel',{'Run','Rest'});
    ylabel('Rate (Hz)'); title(spkNames{s}); box off;
    yl=ylim; r=yl(2)-yl(1); yb=yl(2)+r*0.02;
    plot([1 1 2 2],[yb yb+r*0.03 yb+r*0.03 yb],'k','LineWidth',1);
    text(1.5,yb+r*0.06,sig_star(p(1,2)),'HorizontalAlignment','center','FontSize',12,'FontWeight','bold');
    ylim([yl(1) yb+r*0.16]);
end
% Row 2: firing rate (Hz), in PF vs out PF  (paired t-test, labelled above)
for s=1:nSpk
    nexttile([1 1]);
    p=Boxplot_wPoints2(PFRate(:,[s s+nSpk]),[1 0.4 0.1; 0.6 0.6 0.6]); hold all;
    set(gca,'xtick',[1 2],'xticklabel',{'In PF','Out PF'});
    ylabel('Rate (Hz)'); title(spkNames{s}); box off;
    yl=ylim; r=yl(2)-yl(1); yb=yl(2)+r*0.02;
    plot([1 1 2 2],[yb yb+r*0.03 yb+r*0.03 yb],'k','LineWidth',1);
    text(1.5,yb+r*0.06,sig_star(p(1,2)),'HorizontalAlignment','center','FontSize',12,'FontWeight','bold');
    ylim([yl(1) yb+r*0.16]);
end

% Row 3, left span: running vs resting fraction of each spike type
nexttile([1 2]); hold all;
mRun=[mean(runProb(:,1:nSpk),1,'omitnan')' mean(runProb(:,nSpk+1:end),1,'omitnan')']; % nSpk x [run rest]
hb=bar(mRun,'stacked','FaceColor','flat');
hb(1).CData=repmat([0 0.6 0.2],nSpk,1); hb(2).CData=repmat([0.75 0.75 0.75],nSpk,1);
% SEM (across neurons) of the bottom fraction, drawn at the segment boundary
semRun=std(runProb(:,1:nSpk),0,1,'omitnan')./sqrt(sum(~isnan(runProb(:,1:nSpk)),1));
errorbar(1:nSpk,mRun(:,1)',semRun,'k','LineStyle','none','CapSize',6,'LineWidth',1);
for s=1:nSpk; hc=plot([s-0.45 s+0.45],mean(occRun,'omitnan')*[1 1],'k--','LineWidth',1.2); end
for s=1:nSpk; if p_run(s)<0.05; text(s,1.06,sig_star(p_run(s)),'HorizontalAlignment','center','FontSize',12,'FontWeight','bold'); end; end
set(gca,'xtick',1:nSpk,'xticklabel',spkNames); ylim([0 1.18]); box off;
ylabel('Fraction of spikes'); title('Running vs resting (fraction)');
legend([hb(1) hb(2) hc],{'Running','Resting','chance'},'Location','eastoutside');

% Row 3, right span: in vs out of place field fraction of each spike type
nexttile([1 2]); hold all;
mPF=[mean(PFprob(:,1:nSpk),1,'omitnan')' mean(PFprob(:,nSpk+1:end),1,'omitnan')'];   % nSpk x [inPF outPF]
hb2=bar(mPF,'stacked','FaceColor','flat');
hb2(1).CData=repmat([1 0.5 0.1],nSpk,1); hb2(2).CData=repmat([0.75 0.75 0.75],nSpk,1);
% SEM (across neurons) of the bottom fraction, drawn at the segment boundary
semPF=std(PFprob(:,1:nSpk),0,1,'omitnan')./sqrt(sum(~isnan(PFprob(:,1:nSpk)),1));
errorbar(1:nSpk,mPF(:,1)',semPF,'k','LineStyle','none','CapSize',6,'LineWidth',1);
for s=1:nSpk; hc2=plot([s-0.45 s+0.45],mean(occPF,'omitnan')*[1 1],'k--','LineWidth',1.2); end
for s=1:nSpk; if p_pf(s)<0.05; text(s,1.06,sig_star(p_pf(s)),'HorizontalAlignment','center','FontSize',12,'FontWeight','bold'); end; end
set(gca,'xtick',1:nSpk,'xticklabel',spkNames); ylim([0 1.18]); box off;
ylabel('Fraction of spikes'); title('Place field in vs out (fraction)');
legend([hb2(1) hb2(2) hc2],{'In PF','Out PF','chance'},'Location','eastoutside');

sgtitle('Rows 1-2: firing rate (Hz).   Row 3: fraction of spikes (stacked to 1), dashed = chance, * = vs chance');
set_font('Arial'); set_fontsize(11);

%% dSpike across laps and track position (PlaceTrigger_average)
fprintf('Section: dSpike laps x position map\n');
% Bin each dSpike train by lap x track position with PlaceTrigger_average to see
% where along the track (and in which laps) dSpikes occur.
place_bin=150; vel_thresh=0.002; lap_dist=115;
LapdSp=cell(1,max(foi)); LapdSpFR=cell(1,max(foi)); LapdSpN=cell(1,max(foi));

figure(283); clf; tiledlayout('flow','Padding','compact','TileSpacing','compact');
for f=foi
    nTime=size(PCresult{f}.NormalizedTrace_dirt,2);
    dSpikeTrain=ind2vec(nTime,FiringRate{f}.DendOnlySpFrame,1);
    [LapdSpFR{f}, LapdSp{f}, LapdSpN{f}]=PlaceTrigger_average(dSpikeTrain,place_bin,PCresult{f}.VR,vel_thresh,lap_dist);
    nexttile; hold all;
    nLaps=size(LapdSp{f},1);
    imagesc([1:place_bin]/place_bin*lap_dist,[1:nLaps],LapdSp{f});
    colormap(turbo); caxis([0 1]);
    xlabel('Track position (cm)'); ylabel('Lap'); title(sprintf('Neuron %d',f)); axis tight;
    % shade place field as a rectangle: [start lap, end lap] x [start position, end position]
    if ~isempty(PlaceFieldBin{f})
        for p=1:length(PlaceFieldBin{f})/2
            pb=PlaceFieldBin{f}(2*(p-1)+[1 2]);            % [start end] position bin
            pl=PlaceFieldList{f}(2*(p-1)+[1 2]);           % [start end] lap
            pl=[max(pl(1),1) min(pl(2),nLaps)];            % clamp laps to the map
            if pb(1)<=pb(2), xseg={[pb(1) pb(2)]};
            else,            xseg={[pb(1) place_bin],[1 pb(2)]}; % field wraps past teleport
            end
            for xs=1:numel(xseg)
                xc=[xseg{xs}(1) xseg{xs}(2) xseg{xs}(2) xseg{xs}(1)]/place_bin*lap_dist;
                yc=[pl(1) pl(1) pl(2) pl(2)];
                patch(xc,yc,[1 1 1],'FaceAlpha',0.18,'EdgeColor','w','LineWidth',1.2);
            end
        end
    end
end
sgtitle('dSpike count per lap \times position');
set_font('Arial'); set_fontsize(11);

% Example neuron: lap x position map + position tuning (rate averaged over laps)
fExample=18;
figure(284); clf; tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
posAxis=[1:place_bin]/place_bin*lap_dist;
nexttile([1 1]);
imagesc(posAxis,[1:size(LapdSp{fExample},1)],LapdSp{fExample}); colormap(turbo); caxis([0 1]);
ylabel('Lap'); title(sprintf('Neuron %d  dSpike (count)',fExample)); axis tight;
if ~isempty(PlaceFieldBin{fExample})
    for p=1:length(PlaceFieldBin{fExample})/2
        xl=sort(PlaceFieldBin{fExample}(2*(p-1)+[1 2]))/place_bin*lap_dist;
        xline(xl(1),'w-','LineWidth',1.5); xline(xl(2),'w-','LineWidth',1.5);
    end
end
nexttile([1 1]);
plot(posAxis,mean(LapdSpFR{fExample},1,'omitnan')*1000,'k','LineWidth',1.5);
xlabel('Track position (cm)'); ylabel('dSpike rate (Hz)'); box off; axis tight;
set_font('Arial'); set_fontsize(13);

%% Is the dSpike position-tuned? (spatial tuning + shuffle test)
fprintf('Section: dSpike position tuning\n');
% Occupancy-normalized spatial tuning curve of dSpikes (Hz per track-position bin),
% and a Skaggs spatial-information significance test against a circular-shift null:
% shift the dSpike train by a random offset (1000x), recompute the tuning + info,
% and compare the observed info to the null. p = (1+#{null>=obs})/(nShuffle+1).
nShuffle_pt=1000; rng(1);
smooth_bin=5;                           % circular moving-mean window (bins) for tuning curves
dSpikeTuning=NaN(max(foi),place_bin);   % Hz per position bin, per neuron (smoothed)
SI_obs=NaN(1,max(foi)); SI_p=NaN(1,max(foi));
for f=foi
    nTime=size(PCresult{f}.NormalizedTrace_dirt,2);
    VR=PCresult{f}.VR;
    not_running = VR(end,:) < vel_thresh;
    bin_dist=ceil(VR(5,:)/(lap_dist/place_bin)); bin_dist=min(max(bin_dist,1),place_bin);
    dTrain=ind2vec(nTime,FiringRate{f}.DendOnlySpFrame,1);
    occSum=sum(LapdSpN{f},1,'omitnan'); occSum(occSum==0)=NaN;   % valid running frames per bin
    spk=dTrain; spk(not_running)=0;                              % running-only dSpikes
    cntObs=accumarray(bin_dist(spk>0)',1,[place_bin 1])';
    dSpikeTuning(f,:)=ringmovMean(cntObs./occSum*1000,smooth_bin);   % Hz, circularly smoothed
    SI_obs(f)=skaggs_info(dSpikeTuning(f,:),occSum);
    SI_sh=zeros(1,nShuffle_pt);
    for sh=1:nShuffle_pt
        shifted=circshift(dTrain,randi(nTime)); shifted(not_running)=0;
        c=accumarray(bin_dist(shifted>0)',1,[place_bin 1])';
        SI_sh(sh)=skaggs_info(ringmovMean(c./occSum*1000,smooth_bin),occSum); % same smoothing as obs
    end
    SI_p(f)=(1+sum(SI_sh>=SI_obs(f)))/(nShuffle_pt+1);
end
sigMask=SI_p(foi)<0.05;
fprintf('Position-tuned dSpike neurons: %d/%d (Skaggs info, shuffle p<0.05)\n',sum(sigMask),numel(foi));

% normalized tuning curves across neurons, sorted by peak position
Tn=dSpikeTuning(foi,:); Tn=Tn./max(Tn,[],2,'omitnan');     % peak-normalized per neuron
[~,pk]=max(Tn,[],2,'omitnan'); [~,ord]=sort(pk);
figure(285); clf; tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
nexttile([1 1]);
imagesc(posAxis,1:numel(foi),Tn(ord,:),[0 1]); colormap(turbo);
xlabel('Track position (cm)'); ylabel('Neuron (sorted by peak)');
title('Normalized dSpike tuning'); cb=colorbar; cb.Label.String='Norm. rate';
% mark significant neurons on the y-axis
hold all; sig_ord=find(sigMask(ord));
plot(zeros(numel(sig_ord),1)+posAxis(1),sig_ord,'w*','MarkerSize',6);

% example neuron tuning with shuffle 95% band
nexttile([1 1]); hold all;
f=fExample; nTime=size(PCresult{f}.NormalizedTrace_dirt,2);
VR=PCresult{f}.VR; not_running=VR(end,:)<vel_thresh;
bin_dist=ceil(VR(5,:)/(lap_dist/place_bin)); bin_dist=min(max(bin_dist,1),place_bin);
dTrain=ind2vec(nTime,FiringRate{f}.DendOnlySpFrame,1);
occSum=sum(LapdSpN{f},1,'omitnan'); occSum(occSum==0)=NaN;
shufT=zeros(nShuffle_pt,place_bin);
for sh=1:nShuffle_pt
    shifted=circshift(dTrain,randi(nTime)); shifted(not_running)=0;
    shufT(sh,:)=ringmovMean(accumarray(bin_dist(shifted>0)',1,[place_bin 1])'./occSum*1000,smooth_bin);
end
band=prctile(shufT,[2.5 97.5],1);
fill([posAxis fliplr(posAxis)],[band(1,:) fliplr(band(2,:))],[0.8 0.8 0.8],'EdgeColor','none','FaceAlpha',0.6);
plot(posAxis,dSpikeTuning(f,:),'k','LineWidth',1.5);
if ~isempty(PlaceFieldBin{f})
    for p=1:length(PlaceFieldBin{f})/2
        xl=sort(PlaceFieldBin{f}(2*(p-1)+[1 2]))/place_bin*lap_dist;
        xline(xl(1),'r--'); xline(xl(2),'r--');
    end
end
xlabel('Track position (cm)'); ylabel('dSpike rate (Hz)'); box off; axis tight;
title(sprintf('Neuron %d (p=%.3f)',fExample,SI_p(fExample)));
legend({'shuffle 95%','observed'},'Location','best');

% observed spatial info; significant (p<0.05) neurons marked with a star above the bar
nexttile([1 1]); hold all;
bar(1:numel(foi),SI_obs(foi),'FaceColor',[0.6 0.6 0.6]);
sigIdx=find(sigMask);
text(sigIdx,SI_obs(foi(sigMask)),repmat('*',numel(sigIdx),1),...
    'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14,'FontWeight','bold');
set(gca,'xtick',1:numel(foi),'xticklabel',foi); xtickangle(90);
xlabel('Neuron'); ylabel('Spatial info (bits/spike)'); box off;
title(sprintf('%d/%d position-tuned (p<0.05)',sum(sigMask),numel(foi)));
set_font('Arial'); set_fontsize(12);

%% bAP amplitude in place fields
fprintf('Section: bAP amplitude in place fields\n');
nTau_short=[1 3];
nTau_show=[15 5];

placefield_bAPAMP=[];
placefield_SubAmp=[];
placeTfield_bAPAMP=[];
placefield_spClassVec=[];
PF_center=[];
Placefield_LapFR=[];

t_bin=[-200 200:400:4000];

for f=foi
    % --- source everything from robust-dF/F PCresult{f} (adapter for legacy names) ---
    NT_dirt   = PCresult{f}.NormalizedTrace_dirt;   % robust dF/F, distance-sorted ROIs
    Sub_f     = PCresult{f}.Subthreshold;           % robust dF/F subthreshold
    VR_f      = PCresult{f}.VR;                      % row 5 = track position, row 8 = lap
    dax_f     = PCresult{f}.dendaxis;               % distance-from-soma axis (sorted)
    blue_f    = PCresult{f}.BlueStim;
    somaSpike = PCresult{f}.SpikeMat(1,:);           % somatic spike train
    anySpike  = max(PCresult{f}.SpikeMat,[],1);      % spike on any ROI
    [nROI, nTime] = size(NT_dirt);

    roisD_order_ind=cellfun(@find,PCresult{f}.roisD_order,'UniformOutput',false);
    if isempty(cell2mat(roisD_order_ind(1,:)'))
        basalind=cell2mat(roisD_order_ind(2,:)'); %if there is no basal, use soma
    else
        basalind=cell2mat(roisD_order_ind(1,:)');
    end
    apicalind=cell2mat(roisD_order_ind(5,:)'); %apical = distal dend
    somaind=cell2mat(roisD_order_ind(2,:)');

    [bAP_STA, bAP_MAT, sp_time]=get_STA(NT_dirt,somaSpike.*double(blue_f==0),nTau_short(1),nTau_short(2));
    [bAP_Amp, shift]=max(bAP_MAT,[],3);
    subMat=tovec(permute(bAP_Amp,[1 3 2]));

    % Lap x position-bin arrays (reconstructed from PCresult; replace legacy LapSub/LapSpclassVec).
    % LapSub: [lap, 150 position bins, ROI] mean robust-dF/F subthreshold.
    % LapSpc: [lap, 150 position bins, class] spike counts (1 SS, 2 CS, 3 BS).
    posBin = ceil(VR_f(5,:)/(115/150));            % track position -> bin 1..150
    lapIdx = VR_f(8,:);
    nLap   = max(lapIdx);
    okPB   = posBin>=1 & posBin<=150 & lapIdx>=1;
    LapSub = NaN(nLap,150,nROI);
    for roi=1:nROI
        LapSub(:,:,roi)=accumarray([lapIdx(okPB)' posBin(okPB)'], Sub_f(roi,okPB)', ...
            [nLap 150], @(x) mean(x,'omitnan'), NaN);
    end
    LapSpc = zeros(nLap,150,3);
    spcVec = PCresult{f}.SpikeClassVecMat;          % 4 x nTime, row = class (1 SS,2 CS,3 BS,4 dSP)
    for c=1:3
        sel = okPB & spcVec(c,:)>0;
        if any(sel)
            LapSpc(:,:,c)=accumarray([lapIdx(sel)' posBin(sel)'], 1, [nLap 150]);
        end
    end

    if ~isempty(PlaceFieldList{f}) % in place field
        binTrack=(ceil(VR_f(5,:)/((115)/150)));
        PFvec=zeros(1,nTime);
        for p=1:length(PlaceFieldBin{f})/2
            if PlaceFieldBin{f}(2*(p-1)+1)>PlaceFieldBin{f}(2*(p-1)+2) %place field includes teleport
                Pvec=~(binTrack<(PlaceFieldBin{f}(2*(p-1)+1)) & binTrack>(PlaceFieldBin{f}(2*(p-1)+2)));
            else
                Pvec=(binTrack>(PlaceFieldBin{f}(2*(p-1)+1)) & binTrack<(PlaceFieldBin{f}(2*(p-1)+2)));
            end
            Lapvec=(VR_f(8,:)>PlaceFieldList{f}(2*(p-1)+1) & VR_f(8,:)<PlaceFieldList{f}(2*(p-1)+2));
            PFvec=PFvec| (Lapvec & Pvec);

            bin_dist=ceil(VR_f(5,:)/(115/50));
            PF_rescale=[floor(PlaceFieldBin{f}(2*(p-1)+[1])/150*50) ceil(PlaceFieldBin{f}(2*(p-1)+[2])/150*50)];

            bin_STA=[];
            for b=1:max(bin_dist)
                bin_STA(:,:,b)=get_STA(NT_dirt,anySpike.*double(blue_f==0).*double(bin_dist==b).*double(PFvec),nTau_short(1),nTau_short(2));
            end

            tbin_STA=[];
            first_PFspike_time=find(get_spikeOrder(PFvec,anySpike)==1);
            for t=1:length(t_bin)-1
                t_vec=ind2vec(nTime,tovec(first_PFspike_time'+[t_bin(t):t_bin(t+1)-1]),1);
                tbin_STA(:,:,t)=get_STA(NT_dirt,anySpike.*double(blue_f==0).*double(t_vec),nTau_short(1),nTau_short(2));
            end

            bin_STA_Amp=squeeze(max(bin_STA,[],2));
            tbin_STA_Amp=squeeze(max(tbin_STA,[],2));
            tmp=repmat(bin_STA_Amp,1,2); tmp2=repmat(LapSub,1,2); tmp3=repmat(LapSpc,1,2);

            Placefield_LapFR{f,p}=sum(LapSpc([PlaceFieldList{f}(2*(p-1)+1):PlaceFieldList{f}(2*(p-1)+2)],:,:),3,'omitnan');
            if PlaceFieldBin{f}(2*(p-1)+1)>PlaceFieldBin{f}(2*(p-1)+2) %place field includes teleport
                bin_rescale=([PF_rescale(1):PF_rescale(2)+50]-PF_rescale(1))/50*2;
                placefield_bAPAMP{f,p}=[[NaN bin_rescale]; [dax_f' tmp(:,PF_rescale(1):PF_rescale(2)+50)]];
                placefield_SubAmp{f,p}=permute(mean(tmp2([PlaceFieldList{f}(2*(p-1)+1):PlaceFieldList{f}(2*(p-1)+2)],[PlaceFieldBin{f}(2*(p-1)+1):PlaceFieldBin{f}(2*(p-1)+2)+150],:),1,'omitnan'),[3 2 1]);
                placefield_SubAmp{f,p}=[[NaN [PlaceFieldBin{f}(2*(p-1)+1):PlaceFieldBin{f}(2*(p-1)+2)+150]-PlaceFieldBin{f}(2*(p-1)+1)]; [dax_f' placefield_SubAmp{f,p}]];
                placefield_spClassVec{f,p}=tmp3([PlaceFieldList{f}(2*(p-1)+1):PlaceFieldList{f}(2*(p-1)+2)],[PlaceFieldBin{f}(2*(p-1)+1):PlaceFieldBin{f}(2*(p-1)+2)+150],:);
                PF_center(f,p)=mod(mean([PlaceFieldBin{f}(2*(p-1)+1) PlaceFieldBin{f}(2*(p-1)+2)+150]),150)/150*2;
            else
                bin_rescale=([PF_rescale(1):PF_rescale(2)]-PF_rescale(1))/50*2;
                placefield_bAPAMP{f,p}=[[NaN bin_rescale]; [dax_f' tmp(:,PF_rescale(1):PF_rescale(2))]];
                placefield_SubAmp{f,p}=permute(mean(tmp2([PlaceFieldList{f}(2*(p-1)+1):PlaceFieldList{f}(2*(p-1)+2)],[PlaceFieldBin{f}(2*(p-1)+1):PlaceFieldBin{f}(2*(p-1)+2)],:),1,'omitnan'),[3 2 1]);
                placefield_SubAmp{f,p}=[[NaN [PlaceFieldBin{f}(2*(p-1)+1):PlaceFieldBin{f}(2*(p-1)+2)]-PlaceFieldBin{f}(2*(p-1)+1)]; [dax_f' placefield_SubAmp{f,p}]];
                placefield_spClassVec{f,p}=tmp3([PlaceFieldList{f}(2*(p-1)+1):PlaceFieldList{f}(2*(p-1)+2)],[PlaceFieldBin{f}(2*(p-1)+1):PlaceFieldBin{f}(2*(p-1)+2)],:);
                PF_center(f,p)=mod(mean([PlaceFieldBin{f}(2*(p-1)+1) PlaceFieldBin{f}(2*(p-1)+2)]),150)/150*2;
            end
            meanAmp= max(bAP_STA,[],2);
            Normcst= [mean(meanAmp(somaind,:),1,'omitnan'); mean(meanAmp(apicalind,:),1,'omitnan')];
            placeTfield_bAPAMP{f,p}=[mean(tbin_STA_Amp(somaind,:),1,'omitnan'); mean(tbin_STA_Amp(apicalind,:),1,'omitnan')]./Normcst;
        end
    end
end
PF_center(PF_center==0)=NaN;

v=cellfun(@(x) ~isempty(x),placefield_bAPAMP);
PFbAPdat=placefield_bAPAMP(v);

PFbAPdat_normPF=[]; dbin=80;
for p=1:length(PFbAPdat)
    PFbAPdat_normPF{p}=PFbAPdat{p};
    PFbAPdat_normPF{p}(1,2:end)=rescale(PFbAPdat_normPF{p}(1,2:end));
end

figure(31); clf;
imshow_patch(PFbAPdat_normPF,[0.5 2]);
colormap(turbo)
figure(36); clf;
PF_show=[1 2 3 4 7 8 9 10 11 12 13 14];
show3Dbinning(PFbAPdat_normPF(PF_show),[0:0.07:1],[-dbin*3/2:dbin:400],'image')
shading flat; set(gca,'YDir','reverse')
colormap(turbo)
caxis([0.6 1.3])
xlabel('Normalized place field')
ylabel('Distance from soma (\mum)')
cb=colorbar;
cb.Label.String='bAP amplitude';

figure(32); clf;
imshow_patch(placefield_SubAmp(v),[-0.2 0.7]);
colormap(turbo)
PFSubdat=placefield_SubAmp(v);

figure(37); clf; PFSubdat_normPF=[];
for p=1:length(PFSubdat)
    PFSubdat_normPF{p}=PFSubdat{p};
    PFSubdat_normPF{p}(1,2:end)=rescale(PFSubdat_normPF{p}(1,2:end));
end
show3Dbinning(PFSubdat_normPF(PF_show),[0 0.03:0.05:1],[-dbin*3/2:dbin:400],'image')
shading flat; set(gca,'YDir','reverse')
colormap(turbo)
caxis([-0.05 0.25])
xlabel('Normalized place field')
ylabel('Distance from soma (\mum)')
cb=colorbar;
cb.Label.String='Subthreshold';

figure(35); clf;
%imshow_patch(cellfun(@(x) im_merge(ringmovMean(x(:,:,[1 2 3]),3),[1 1 1;1 0 0;1 1 1]), placefield_spClassVec(v),'UniformOutput',false));
tmppf=placefield_spClassVec(v);
for p=1:size(tmppf,1)
    nexttile([1 1])
    l=plot(squeeze(mean(tmppf{p},1,'omitnan')));
    arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell([0 0 0;1 0 0;0 0.5 1],2));
    xlabel('VR position')
    ylabel('Firing rate');
end
legend({'SS','CS','BS'})


figure(33); clf;
% [binnedZ, centerX, centerY]=show3Dbinning(placeTfield_bAPAMP(v),t_bin,[-100:70:450],'image');
% colormap(turbo); caxis([0.5 2]);
% shading flat
% set(gca,'YDir','reverse')
somaTfield=cell2mat(cellfun(@(x) x(1,:),placeTfield_bAPAMP(v),'UniformOutput',false));
apicalTfield=cell2mat(cellfun(@(x) x(2,:),placeTfield_bAPAMP(v),'UniformOutput',false));

dat=apicalTfield;%./max(apicalTfield,[],2);
dat_som=somaTfield;%./somaTfield(:,1);
T=mean([t_bin(2:end); t_bin(1:end-1)])/1000;
M=mean(dat,1,'omitnan'); S=std(dat,0,1,'omitnan'); N=sum(~isnan(dat),1);
Ms=mean(dat_som,1,'omitnan'); Ss=std(dat_som,0,1,'omitnan'); Ns=sum(~isnan(dat_som),1);
plot(T,dat','color',[1 0.7 0.7]);
hold all
plot(T,dat_som','color',[0.7 0.7 1]);
errorbar_shade(T,M,S/sqrt(N),[1 0 0])
errorbar_shade(T,Ms,Ss/sqrt(Ns),[0 0 1])

%% Somatic spike vs dendritic spike, how many spikes are preceding? -> No spikes
fprintf('Section: Somatic vs dendritic spike, preceding spikes\n');
nTau_short=[2 3];
nTau_show=[15 5];
FractionMat=[];
figure(38); clf;
for f=foi
    % --- source everything from robust-dF/F PCresult{f} ---
    NT_dirt   = PCresult{f}.NormalizedTrace_dirt;   % robust dF/F, distance-sorted
    anyClassSpike = max(PCresult{f}.SpikeClassMat,[],1);  % first spike of any class
    nTime = size(PCresult{f}.Subthreshold,2);

    roisD_order_ind=cellfun(@find,PCresult{f}.roisD_order,'UniformOutput',false);
    if isempty(cell2mat(roisD_order_ind(1,:)'))
        basalind=cell2mat(roisD_order_ind(2,:)'); %if there is no basal, use soma
    else
        basalind=cell2mat(roisD_order_ind(1,:)');
    end
    apicalind=cell2mat(roisD_order_ind(5,:)'); %apical = distal dend
    somaind=cell2mat(roisD_order_ind(2,:)');

    [bAP_STA, bAP_MAT sp_time]=get_STA(NT_dirt,anyClassSpike,nTau_short(1),nTau_short(2));
    [~, Sp_STA]=get_STA(PCresult{f}.SpikeMat,anyClassSpike,nTau_short(1),nTau_short(2));
    Sp_STA=max(Sp_STA(PCresult{f}.sortDist,:,1:nTau_short(1)),[],3);   % reorder to distance-sorted ROIs
    [bAP_Amp, shift]=max(bAP_MAT,[],3);
    bAP_STA_Amp=max(bAP_STA,[],2);

    apicalIntSp=find(max(Sp_STA(apicalind,:),[],1)>0);
    basalIntSp=find(max(Sp_STA(basalind,:),[],1)>0);

    apicalIntSp_vec=ind2vec(nTime,sp_time(apicalIntSp),1);
    basalIntSp_vec=ind2vec(nTime,sp_time(basalIntSp),1);

    [bAP_STA_show ]=get_STA(NT_dirt,anyClassSpike,nTau_show(1),nTau_show(2));
    [~, AP_INT_MAT]=get_STA(NT_dirt,apicalIntSp_vec,nTau_show(1),nTau_show(2));
    [~, BS_INT_MAT]=get_STA(NT_dirt,basalIntSp_vec,nTau_show(1),nTau_show(2));

    FractionMat(f,:)=[sum(basalIntSp_vec), sum(apicalIntSp_vec), sum(basalIntSp_vec & apicalIntSp_vec)]/sum(anyClassSpike);

    nexttile([1 1]);
    imagesc([squeeze(mean(BS_INT_MAT,2,'omitnan')) squeeze(mean(AP_INT_MAT,2,'omitnan')) bAP_STA_show])
    title(num2str(f))
end
colormap(turbo)

%% Orthogonal-axis subthreshold correlation
fprintf('Section: Orthogonal-axis subthreshold correlation\n');
% NOTE: in the source this cell referenced variables that were never defined
% (filteredSubthres_orthAxis / _2, max_xcorrT, FilterFreq) and looped a single
% neuron (f=26). Reconstructed here from robust-dF/F PCresult: the per-channel
% subthreshold is projected onto the dendritic orthogonal axes and cross-
% correlated across sessions. The cross-channel (ch1 x ch2) autocorrelation on
% the diagonal removes independent shot noise, as in Figure 2. Validate the
% projection/constants against the intended design before using in a figure.
run_threshold=0.002;
OrthAxis =[1 0 1 0 1; 1 0 -1 0 -1; 0 1 0 0 0]; %basal, soma, trunk, oblique, distal
OrthAxis2=[1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 0 1]; %basal, soma, trunk, distal
max_xcorrT=500;      % cross-correlation half-window (frames)
FilterFreq=[4 10];   % theta band for get_phase (Hz)

nAx1=size(OrthAxis,1); nAx2=size(OrthAxis2,1);
OrthaxisXcorr_silent=NaN(nAx2,nAx2,2*max_xcorrT+1,max(foi));
OrthaxisXcorr       =NaN(nAx2,nAx2,2*max_xcorrT+1,max(foi));
STmat_SubthOrth=cell(max(foi),4); PhaseTr=cell(1,max(foi)); ThetaPow=cell(1,max(foi)); psd_norm=cell(1,max(foi));

for f=foi
    f
    NT_dirt  = PCresult{f}.NormalizedTrace_dirt;            % robust dF/F, distance-sorted
    Sub_c    = PCresult{f}.Subthreshold_larger;             % combined subthreshold
    Sub_ch1  = PCresult{f}.Subthreshold_larger_ch{1};       % per-channel subthreshold
    Sub_ch2  = PCresult{f}.Subthreshold_larger_ch{2};
    blue_f   = PCresult{f}.BlueStim;
    spkClass = PCresult{f}.SpikeClassMat;
    somaSpike= PCresult{f}.SpikeMat(1,:);
    CS_f     = PCresult{f}.CStrace;
    [nROI, nTime]=size(NT_dirt);

    % orthogonal-axis weight vectors (basal/distal seesaw, global, etc.)
    roisD_order_ind=cellfun(@find,PCresult{f}.roisD_order,'UniformOutput',false);
    OrthAxis_vec=zeros(nAx1,nROI); OrthAxis_vec2=zeros(nAx2,nROI);
    for d=1:5
        dclassinds=cell2mat(roisD_order_ind(d,:)');
        for v=1:nAx1
            if abs(OrthAxis(v,d))>0 & isempty(dclassinds)
                OrthAxis_vec(v,:)=NaN;
            else
                suborthvec=ind2vec(nROI,dclassinds,OrthAxis(v,d)/length(dclassinds));
                OrthAxis_vec(v,abs(suborthvec)>0)=suborthvec(abs(suborthvec)>0);
            end
        end
        for v=1:nAx2
            if abs(OrthAxis2(v,d))>0 & isempty(dclassinds)
                OrthAxis_vec2(v,:)=NaN;
            else
                suborthvec2=ind2vec(nROI,dclassinds,OrthAxis2(v,d)/length(dclassinds));
                OrthAxis_vec2(v,abs(suborthvec2)>0)=suborthvec2(abs(suborthvec2)>0);
            end
        end
    end

    % project the (robust-dF/F) subthreshold onto the orthogonal axes
    filteredSubthres_orthAxis  = OrthAxis_vec * Sub_c;                       % [nAx1 x nTime]
    filteredSubthres_orthAxis(:,blue_f>0)=NaN;
    filteredSubthres_orthAxis2 = cat(3, OrthAxis_vec2*Sub_ch1, OrthAxis_vec2*Sub_ch2); % [nAx2 x nTime x 2ch]
    filteredSubthres_orthAxis2(:,blue_f>0,:)=NaN;

    for dclass=1:4
        [~, STmat_SubthOrth{f,dclass}]=get_STA(filteredSubthres_orthAxis, spkClass(dclass,:).*(blue_f==0),30,30);
    end

    silenttime=setdiff([1:nTime],unique(find(max([somaSpike; spkClass; CS_f; blue_f])>0)'+[-5:30]));
    silenttime_vec=ind2vec(nTime,silenttime,1);
    filteredSubthres_orthAxis_silent2=filteredSubthres_orthAxis2;
    filteredSubthres_orthAxis_silent2(:,~silenttime_vec,:)=NaN;
    ph_src=interpolateNaN(filteredSubthres_orthAxis_silent2(:,:,1));   % axes x time, channel 1

    for oi=1:nAx2
        for oj=oi:nAx2
            OrthaxisXcorr_silent(oi,oj,:,f)=nanXCorr(filteredSubthres_orthAxis_silent2(oi,:,1),filteredSubthres_orthAxis_silent2(oj,:,2),max_xcorrT,1);
            OrthaxisXcorr(oi,oj,:,f)=nanXCorr(filteredSubthres_orthAxis2(oi,:,1),filteredSubthres_orthAxis2(oj,:,2),max_xcorrT,1);
        end
        if any(~isnan(ph_src(oi,:)))
            [PhaseTr{f}(oi,:),~,ThetaPow{f}(oi,:)] = get_phase(ph_src(oi,:), 1000, FilterFreq);
            [frq, pw] = nanPSD(filteredSubthres_orthAxis_silent2(oi,:,1), 1000, 3000);
            psd_norm{f}=[frq pw];
        end
    end
end

% Basal & apical & soma autocorr and crosscorr
figure(127); clf; t_lag=[-max_xcorrT:max_xcorrT]; ax2=[];
subplot_loc=[1 3 2 9 6 5]; g=1; ROI_str={'Basal','Soma','Trunk','Distal'};
cmap_label=jet(9);
for oi=1:size(OrthAxis2,1)
    for oj=oi:size(OrthAxis2,1)
        nexttile([1 1])
        plot(t_lag,squeeze(OrthaxisXcorr_silent(oi,oj,:,foi)),'color',[0.8 0.8 0.8]); hold all
        errorbar_shade(t_lag,mean(squeeze(OrthaxisXcorr_silent(oi,oj,:,foi))',1,'omitnan'),std(squeeze(OrthaxisXcorr_silent(oi,oj,:,foi))',0,1,'omitnan'),cmap_label(g,:));
        title([ROI_str{oi} ' & ' ROI_str{oj}])
        g=g+1;
        xlabel('\tau (ms)')
        ylabel('R(\tau)')
    end
end

%% CS vs SS
figure(300); clf;
stype_str={'SS','CS','BS'};
OrthAxis=[1 0 0 0 0; 0 0 1 0 1; 0 0 0 0 1; 0 1 0 0 0]; %rows: basal, apical(trunk+distal), distal, soma
cmaps=gen_colormap(Plasma,size(OrthAxis,1));
ax1=[]; Prespike_t=5; nTau2show=[150 50];
foi2show=[18 23 25 26];
cax=[-0.02 0.1]; cax_scatter=[-0.3 0.5]; presubtime2ave=[10 3];
branch_str={'Basal','Apical','Distal','Soma'};   % rows of OrthAxis
branchPlotOrder=[1 4 2 3];                        % trace display order: Basal, Soma, Apical, Distal
cmap_SSCS=[1 0 1; 0 0.7 0.7];                     % SS (magenta), CS (teal)
for f=foi2show
    % --- source everything from robust-dF/F PCresult{f} (adapter for legacy names) ---
    Sub_f         = PCresult{f}.Subthreshold;          % robust dF/F subthreshold
    NT_dirt       = PCresult{f}.NormalizedTrace_dirt;  % robust dF/F, distance-sorted
    SpikeMat_f    = PCresult{f}.SpikeMat;              % row 1 = somatic spike train
    SpikeClass_f  = PCresult{f}.SpikeClassVecMat;      % 4 x nTime (1 SS, 2 CS, 3 BS, 4 dSP)
    CS_f          = PCresult{f}.CStrace;
    blue_f        = PCresult{f}.Blue;
    roisD_order_f = PCresult{f}.roisD_order;
    dax_f         = PCresult{f}.dendaxis;

    [nROI nTime]=size(Sub_f); perispike_time=[-2:20];
    [SubV SubD subTrace]=get_eigvector(Sub_f(:,sum(isnan(Sub_f),1)==0),nROI);
    subTrace_onFrame=NaN(nROI,nTime);
    subTrace_onFrame(:,sum(isnan(Sub_f),1)==0)=subTrace';

    perispike_frame=ind2vec(nTime,unique([tovec(find(double(SpikeMat_f(1,:)==1))'+perispike_time); find(CS_f)']),1);
    Blue_on_frame=imdilate(blue_f>0, [ones(1, 1), 1, ones(1, 200)])>0;
    Blue_off_frame=imdilate(blue_f==0, [1, ones(1, 50)])>0;

    roisD_order_ind=cellfun(@find,roisD_order_f,'UniformOutput',false);
    OrthAxis_vec=zeros(size(OrthAxis,1),nROI);
    for d=1:5 %dClass
        dclassinds=cell2mat(roisD_order_ind(d,:)');
        for v=1:size(OrthAxis,1)
            if abs(OrthAxis(v,d))>0 & isempty(dclassinds)
                OrthAxis_vec(v,:)=NaN;
            else
                suborthvec=ind2vec(nROI,dclassinds,OrthAxis(v,d)/length(dclassinds));
                OrthAxis_vec(v,abs(suborthvec)>0)=suborthvec(abs(suborthvec)>0);
            end
        end
    end
    %OrthAxis_vec=OrthAxis_vec(:,cell2mat(roisD_order_ind(:)));

    SubBasalApical=(Sub_f'*OrthAxis_vec')';
    %NormTr_filtered=pcafilterTrace(NT_dirt,[1:5]);
    NormBasalApical=(NT_dirt'*OrthAxis_vec')';
    % ---- STA kymographs, branch traces, and isolated pre-spike scatter (SS, CS) ----
    kymoStore=cell(1,2); kymoVecStore=cell(1,2); nSPstore=zeros(1,2); scatterStore=cell(1,2);
    preCols=nTau2show(1)+1+[-presubtime2ave(1):-presubtime2ave(2)];   % -10..-3 ms pre-spike window
    for stype=1:2
        [~, STAmat, spframes]=get_STA(NT_dirt,max(SpikeClass_f(stype,:),[],1).*Blue_off_frame,nTau2show(1),nTau2show(2));
        [~, SpClassMat]=get_STA(SpikeClass_f+repmat(CS_f,4,1),spframes,nTau2show(1),nTau2show(2));
        SPisolated=sum(double(SpClassMat(:,:,1:nTau2show(1))),[1 3])==0; % no other event in the 150 ms pre-window
        nSPstore(stype)=sum(SPisolated);
        kymo2show=permute(mean(STAmat(:,SPisolated,:),2,'omitnan'),[1 3 2]);
        kymo2show=kymo2show-median(kymo2show(:,[1:nTau2show(1)]),2,'omitnan');
        kymoStore{stype}=kymo2show;
        kymoVecStore{stype}=(kymo2show'*OrthAxis_vec')';
        % pre-spike subthreshold of the SAME isolated events (for the figure-300 scatter)
        [~, subBAwin]=get_STA(SubBasalApical,spframes,nTau2show(1),nTau2show(2));
        subPre=mean(subBAwin(:,:,preCols),3,'omitnan');
        scatterStore{stype}=subPre(:,SPisolated);
    end

    % ---- figure 300: basal-vs-apical density with the isolated SS / CS points ----
    figure(300);
    nexttile([1 1])
    try
        [hist_subV_A_B]=scatter_heatmap(SubBasalApical(1,Blue_off_frame),SubBasalApical(2,Blue_off_frame),50,50);
        hold all
        hSS=scatter(scatterStore{1}(1,:),scatterStore{1}(2,:),10,cmap_SSCS(1,:),'filled');
        hCS=scatter(scatterStore{2}(1,:),scatterStore{2}(2,:),10,cmap_SSCS(2,:),'filled');
        xlabel('Basal'); ylabel('Apical');
        colormap('bone'); title(['ID# ' num2str(f) '  (SS ' num2str(nSPstore(1)) ', CS ' num2str(nSPstore(2)) ')']);
        %xlim(cax_scatter); ylim(cax_scatter);
        axis tight;
        if f==foi2show(1), legend([hSS hCS],{'SS','CS'},'Box','on'); end
    end

    % ---- one figure per neuron: kymographs (top) + branch traces (bottom, incl. soma) ----
    figure(3030+f); clf; tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
    sgtitle(['Neuron ID# ' num2str(f)]);
    t_axis=[-nTau2show(1):nTau2show(2)];
    for stype=1:2
        nexttile([1 2]);
        imagesc(kymoStore{stype},cax)
        title([stype_str{stype} ', N=' num2str(nSPstore(stype))])
        if stype==1, set_kymoYtick(dax_f); else, set(gca,'ytick',[]); end
        set(gca,'XTick',[1 nTau2show(1)+1 sum(nTau2show)+1],'XTickLabel',num2str([-nTau2show(1) 0 nTau2show(end)]',3))
    end
    for v=branchPlotOrder
        nexttile([1 1]); hold on; box off;
        plot(t_axis,kymoVecStore{1}(v,:),'color',cmap_SSCS(1,:),'linewidth',1.5);
        plot(t_axis,kymoVecStore{2}(v,:),'color',cmap_SSCS(2,:),'linewidth',1.5);
        ylim([-0.2 1.5]*PCresult{f}.SpikeHeight); xlim([-150 50]);
        title(branch_str{v}); xlabel('Time (ms)');
        if v==branchPlotOrder(1), ylabel('\DeltaF/F'); legend({'SS','CS'},'Box','off','location','northwest'); end
    end
    set_fontsize(12); colormap(turbo(256));
end
figure(300); set_fontsize(12);

MeanPrespikeSubthreshold=[];
nIsolatedSp=zeros(max(foi),3);   % # isolated spikes averaged per neuron per class (1 SS, 2 CS, 3 BS)
for f=foi
    % --- source everything from robust-dF/F PCresult{f} (adapter for legacy names) ---
    Sub_f         = PCresult{f}.Subthreshold;          % robust dF/F subthreshold
    NT_dirt       = PCresult{f}.NormalizedTrace_dirt;  % robust dF/F, distance-sorted
    SpikeMat_f    = PCresult{f}.SpikeMat;              % row 1 = somatic spike train
    SpikeClass_f  = PCresult{f}.SpikeClassVecMat;      % 4 x nTime (1 SS, 2 CS, 3 BS, 4 dSP)
    CS_f          = PCresult{f}.CStrace;
    blue_f        = PCresult{f}.BlueStim;
    roisD_order_f = PCresult{f}.roisD_order;

    [nROI nTime]=size(Sub_f); perispike_time=[-2:20];
    [SubV SubD subTrace]=get_eigvector(Sub_f(:,sum(isnan(Sub_f),1)==0),nROI);
    subTrace_onFrame=NaN(nROI,nTime);
    subTrace_onFrame(:,sum(isnan(Sub_f),1)==0)=subTrace';

    perispike_frame=ind2vec(nTime,unique([tovec(find(double(SpikeMat_f(1,:)==1))'+perispike_time); find(CS_f)']),1);
    Blue_on_frame=imdilate(blue_f>0, [ones(1, 1), 1, ones(1, 200)])>0;
    Blue_off_frame=imdilate(blue_f==0, [1, ones(1, 50)])>0;

    roisD_order_ind=cellfun(@find,roisD_order_f,'UniformOutput',false);
    OrthAxis_vec=zeros(size(OrthAxis,1),nROI);
    for d=1:5 %dClass
        dclassinds=cell2mat(roisD_order_ind(d,:)');
        for v=1:size(OrthAxis,1)
            if abs(OrthAxis(v,d))>0 & isempty(dclassinds)
                OrthAxis_vec(v,:)=NaN;
            else
                suborthvec=ind2vec(nROI,dclassinds,OrthAxis(v,d)/length(dclassinds));
                OrthAxis_vec(v,abs(suborthvec)>0)=suborthvec(abs(suborthvec)>0);
            end
        end
    end

    SubBasalApical=(Sub_f'*OrthAxis_vec')';
    NormBasalApical=(NT_dirt'*OrthAxis_vec')';
    for stype=1:3
        [~, STAmat, spframes]=get_STA(NT_dirt,max(SpikeClass_f(stype,:),[],1).*Blue_off_frame,nTau2show(1),nTau2show(2));
        [~, SpClassMat]=get_STA(SpikeClass_f+repmat(CS_f,4,1),spframes,nTau2show(1),nTau2show(2));
        if length(spframes)>3
        SPisolated=sum(double(SpClassMat(:,:,1:nTau2show(1))),[1 3])==0; % no SS/CS/BS/dSP in the pre-spike window
        nIsolatedSp(f,stype)=sum(SPisolated);
        kymo2show=permute(mean(STAmat(:,SPisolated,:),2,'omitnan'),[1 3 2]);
        kymo2show_std=permute(std(STAmat(:,SPisolated,:),0,2,'omitnan')./sqrt(sum(SPisolated)),[1 3 2]);

        kymo2show=kymo2show-median(kymo2show(:,[1:nTau2show(1)]),2,'omitnan');        
        kymo2show_vec=(kymo2show'*OrthAxis_vec')'; kymo2show_std_vec=(kymo2show_std'*OrthAxis_vec')';

        MeanPrespikeSubthreshold(f,:,stype)=mean(kymo2show_vec(:,nTau2show(1)+[-presubtime2ave(1):-presubtime2ave(2)]),2,'omitnan');
        else
        MeanPrespikeSubthreshold(f,:,stype)=NaN(size(OrthAxis_vec,1),1);
        end
    end
end

figure(302); clf; cmap_stype=[1 0 1; 0 0.7 0.7];
title_str={'Basal','Apical','Distal'}; g=1;
ylim_302={[-0.03 0.05],[-0.03 0.05],[-0.01 0.1]};   % per-panel y-limits (Basal, Apical, Distal)
SpikeHeight_foi=arrayfun(@(ff) PCresult{ff}.SpikeHeight, foi)'; % robust dF/F spike height per neuron
for dclass=[1 3 2]
    nexttile([1 1]);
    MeanPrespikeSubthreshold_Sp=permute(MeanPrespikeSubthreshold(foi,dclass,[1 2]),[1 3 2]);
    MeanPrespikeSubthreshold_Sp(nIsolatedSp(foi,2)<10,:)=NaN;   % drop neurons with <10 isolated CS
    p=Boxplot_wPoints2(MeanPrespikeSubthreshold_Sp,cmap_stype);
    drawPValueLines(p,0);
    box off;
    % # neurons (finite points) and # isolated spikes averaged per group
    validSS=isfinite(MeanPrespikeSubthreshold_Sp(:,1)); validCS=isfinite(MeanPrespikeSubthreshold_Sp(:,2));
    nSpSS=sum(nIsolatedSp(foi(validSS),1)); nSpCS=sum(nIsolatedSp(foi(validCS),2));
    set(gca,'xtick',[1 2],'xticklabel',{sprintf('SS'), ...
                                        sprintf('CS')});
    ylabel('Pre-spike subthreshold (∆F/F)');
    ylim(ylim_302{g});
    title(title_str{g}); g=g+1;
end
set_figsize(200,90)

%% [3-setup] Compute pre-spike theta phase / amplitude / dV/dt / PLV / STA / ETA  (feeds 3a-3e)
% Multi-compartment (basal / apical / distal / soma) PRE-spike theta phase at the
% onset of each event, estimated from pre-onset dynamics ONLY (get_prespike_phase),
% plus the causal distal-basal theta phase-locking (coherence) in the pre-onset
% window (get_prespike_coherence). We compare where SS, CS and (silent-period)
% dSpikes sit on each compartment's theta, and whether SS/CS occur when distal
% and basal theta are PHASE-MATCHED (signed index <cos(dphi)>: +1 in-phase, -1
% anti-phase). Run this cell once; then run any of [3a]-[3e].
%
% Isolation (req 2): an SS/CS onset is kept only if NO CS falls in its pre-window
% (a preceding SS is allowed). dSpikes are taken from SILENT periods (no SS/CS/BS
% in the pre-window). "Silent" = random frames whose 600 ms pre-window has no
% somatic event, used as the coherence baseline.
%
% NOTE: pre-spike V / dV/dt are a straight-line fit over preWin on the RAW traces
% (not the subthreshold) -- raw is not smoothed, so it has no post-spike leakage and
% the window can run right up to -2 ms. The spike upstroke is later than -2 ms.
if isempty(which('circ_rtest'))     % make CircStat2012a available if it isn't
    addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'..','CircStat2012a')));
end
rng(0);
OrthAxis=[1 0 0 0 0; 0 0 1 0 1; 0 0 0 0 1; 0 1 0 0 0]; % rows: basal, apical, distal, soma
compName={'Basal','Apical','Distal','Soma'}; basalAx=1; apicalAx=2; distalAx=3; somaAx=4;
thetaBand=[5 12]; phaseW_ms=600; phaseGuard_ms=40;     % Hz; causal pre-window / edge guard (ms)
etypeName={'SS','CS','BS','dSpike','Silent'};              % real spike types 1-4, silent baseline = 5
cmap_et=[1 0 1; 0 0.7 0.7; 0.5 0.2 0.75; 0.9 0.6 0; 0.6 0.6 0.6];
nET=numel(etypeName); nRealET=4; nComp=numel(compName);   % nRealET = SS,CS,BS,dSpike
preWin=-20:-2;                                      % pre-spike dV/dt / V window (ms), fit on RAW traces

staTau=[50 20];                                     % STA window (ms pre, ms post) for fig 308
csAmpWin=0:80; ssAUCwin=0:50;                       % soma post-onset windows for fig-309 colors
Phase_all  = repmat({nan(0,nComp)},max(foi),nET);   % {f,et}: onsets x compartment phase (rad)
Mag_all    = repmat({nan(0,nComp)},max(foi),nET);   % {f,et}: onsets x compartment theta amplitude
dVdtPre_all= repmat({nan(0,nComp)},max(foi),nET);   % {f,et}: onsets x compartment pre-spike dV/dt
Vpre_all   = repmat({nan(0,nComp)},max(foi),nET);   % {f,et}: onsets x compartment pre-spike V (mean of preWin)
PeakSlow_all=repmat({nan(0,nComp)},max(foi),nET);   % {f,et}: onsets x compartment PEAK slow component (fig 310)
peakSlowWin=0:150;                                  % post-onset window (ms) for the peak slow component
ColorVal_all=repmat({nan(0,1)},   max(foi),nET);    % {f,et}: onsets x soma color value (fig 309)
AUCpost_all =repmat({nan(0,1)},   max(foi),nET);    % {f,et}: onsets x soma post-spike AUC (fig 309 phase-diff)
Coh_all    = repmat({nan(0,1)},   max(foi),nET);    % {f,et}: onsets x distal-basal PLV
Dphi_all   = repmat({nan(0,1)},   max(foi),nET);    % {f,et}: onsets x distal-basal phase offset
STA_comp   = repmat({nan(nComp,sum(staTau)+1)},max(foi),2);  % {f,stype}: subthreshold compartment STA (SS,CS)
STA_comp_raw=repmat({nan(nComp,sum(staTau)+1)},max(foi),2);  % {f,stype}: RAW compartment STA (fig 307 dV/dt)
phaseEdges=-pi:pi/5:pi; nPB=numel(phaseEdges)-1;             % 10 phase bins for the binned STA (fig 309)
phaseCtr=(phaseEdges(1:end-1)+phaseEdges(2:end))/2;
staTau9=[150 150]; nT9=sum(staTau9)+1;                        % fig-309 STA window (-150..50 ms, >=1 theta cycle)
verifTau=[150 50]; nVT=sum(verifTau)+1;                      % window for the phase-detection check (fig 303)
PV_sum=zeros(nComp,nPB,nVT); PV_cnt=zeros(nComp,nPB);        % filtered trace binned by estimated phase
STAbin_sum=zeros(nComp,nPB,nT9,2); STAbin_cnt=zeros(nComp,nPB,2);   % phase-binned STA (comp x bin x t x SS/CS)
STAdiff_sum=zeros(nPB,nT9,2); STAdiff_cnt=zeros(nPB,2);             % distal-basal-phase-diff-binned soma STA
% fig-311 (SS) / fig-312 (CS): dendritic STA kymograph, binned into a 2x2 grid of
% (basal-phase sign) x (distal-phase sign) -- the four quadrants of the fig-306 basal-vs-distal
% phase scatter. Pooled across neurons by dendritic distance (RAW trace, baseline-subtr, /spike ht).
kymoTau=[100 250]; nKT=sum(kymoTau)+1; kymoT=-kymoTau(1):kymoTau(2);      % kymograph window (ms)
kymoEdges=-300:40:460; kymoCtr=(kymoEdges(1:end-1)+kymoEdges(2:end))/2; nKB=numel(kymoCtr);
KymoSum=zeros(nKB,nKT,4,2); KymoN=zeros(nKB,nKT,4,2); nEvQuad=zeros(4,2); % 4 quadrants x (SS,CS)
for f=foi
    Sub_f        = PCresult{f}.Subthreshold;
    SpikeClass_f = PCresult{f}.SpikeClassVecMat;
    CS_f         = PCresult{f}.CStrace;
    blue_f       = PCresult{f}.BlueStim;
    [nROI, nTime]= size(Sub_f);
    W = round(phaseW_ms);

    OrthAxis_vec = build_orthaxis_vec(OrthAxis, PCresult{f}.roisD_order, nROI);
    SubBA = OrthAxis_vec * Sub_f;                 % 4 x nTime: basal, apical, distal, soma (subthreshold)
    rawComp = OrthAxis_vec * PCresult{f}.NormalizedTrace_dirt;  % 4 x nTime RAW voltage (fig 309 STA)
                                                  % (missing compartments -> NaN row -> dropped later)
    rawFull = PCresult{f}.NormalizedTrace_dirt;   % nROI x nTime RAW voltage (fig 311 kymograph)
    dax_f   = PCresult{f}.dendaxis(:)';           % dendritic distance per ROI (sorted order)
    kbin    = discretize(dax_f, kymoEdges);       % distance-bin index per ROI (NaN outside range)
    Mbin    = zeros(nKB, nROI);                   % distance-bin membership (bins x ROI), fig 311
    for r=1:nROI, if ~isnan(kbin(r)), Mbin(kbin(r),r)=1; end, end
    shF     = PCresult{f}.SpikeHeight;
    % req 4: if this neuron has no basal ROI, use soma as the basal proxy for the
    % distal-basal phase relationships (coherence / phase match / phase difference).
    basalTr = SubBA(basalAx,:); basalProxy = all(isnan(basalTr));
    if basalProxy, basalTr = SubBA(somaAx,:); end
    % event traces
    CS_any  = (SpikeClass_f(2,:)>0) | (CS_f>0);
    Somatic = (max(SpikeClass_f(1:3,:),[],1)>0) | (CS_f>0);   % SS | CS | BS
    dSP     = SpikeClass_f(4,:)>0;

    % mask peri-event & blue so transients don't leak into the theta band
    exclMask = imdilate(max(SpikeClass_f,[],1)>0 | CS_f>0, ones(1,11))>0 | blue_f>0;

    % causal isolation tests (strictly BEFORE the onset, within the pre-window)
    hasCSpre  = @(t) any(CS_any(max(1,t-W):t-1));
    hasSomPre = @(t) any(Somatic(max(1,t-W):t-1));

    % ---- onsets per event type ----
    Et  = EventTable(EventTable.Neuron==f,:);
    frSS = Et.Frame(Et.EventClass=="SS" & Et.IsValid)';
    frCS = Et.Frame(Et.EventClass=="CS" & Et.IsValid)';
    frBS = Et.Frame(Et.EventClass=="BS" & Et.IsValid)';
    frDS = find(diff([0 dSP])==1);                  % dSpike onsets (rising edges)
    frDS = frDS(~(blue_f(min(max(frDS,1),nTime))>0));
    silent = ~imdilate(Somatic, ones(1,201)) & ~(blue_f>0);   % >=100 ms from any somatic event
    silCand = find(silent); silCand = silCand(silCand>W);
    frRND = silCand(randperm(numel(silCand), min(200,numel(silCand))));

    frSS = frSS(arrayfun(@(t) t>W && ~hasCSpre(t),  frSS));    % req 2: reject if preceding CS
    frCS = frCS(arrayfun(@(t) t>W && ~hasCSpre(t),  frCS));
    frBS = frBS(arrayfun(@(t) t>W && ~hasCSpre(t),  frBS));
    frDS = frDS(arrayfun(@(t) t>W && ~hasSomPre(t), frDS));    % silent-period dSpikes
    frRND= frRND(arrayfun(@(t) t>W && ~hasSomPre(t), frRND));  % silent baseline
    onsetSet = {frSS, frCS, frBS, frDS, frRND};               % 1 SS, 2 CS, 3 BS, 4 dSpike, 5 Silent

    for et=1:nET
        ons = onsetSet{et}; if isempty(ons), continue; end
        % per-compartment pre-spike theta phase + amplitude + pre-spike dV/dt (-20..0 ms)
        ph  = nan(numel(ons),nComp); mg = nan(numel(ons),nComp); dvd = nan(numel(ons),nComp);
        for c=1:nComp
            [phc,mgc,~,okc] = get_prespike_phase(SubBA(c,:),1000,thetaBand,ons,exclMask,phaseW_ms,phaseGuard_ms);
            phc(~okc)=NaN; ph(:,c)=phc(:);
            mgc(~okc)=NaN; mg(:,c)=mgc(:);
        end
        vpre = nan(numel(ons),nComp);
        for k=1:numel(ons)
            t0=ons(k); if t0+preWin(1)<1 || t0+preWin(end)>nTime, continue; end
            seg=rawComp(:,t0+preWin);                 % RAW traces (linear fit over preWin)
            for c=1:nComp
                v=seg(c,:); gd=isfinite(v);
                if nnz(gd)>=5, pp=polyfit(preWin(gd),v(gd),1); dvd(k,c)=pp(1); vpre(k,c)=mean(v(gd)); end
            end
        end
        Phase_all{f,et}=ph; Mag_all{f,et}=mg; dVdtPre_all{f,et}=dvd; Vpre_all{f,et}=vpre;
        % per-onset PEAK slow component per compartment (fig 310): max of the smoothed
        % subthreshold (slow) trace over peakSlowWin (0..150 ms), minus the pre-onset baseline.
        pk=nan(numel(ons),nComp); SubBAsm=movmean(SubBA,5,2,'omitnan');
        for k=1:numel(ons)
            t0=ons(k); w=t0+peakSlowWin; w=w(w>=1 & w<=nTime); b0=max(1,t0-20):t0-1;
            if isempty(w) || isempty(b0), continue; end
            base=mean(SubBAsm(:,b0),2,'omitnan');
            pk(k,:)=max(SubBAsm(:,w),[],2,'omitnan')' - base';
        end
        PeakSlow_all{f,et}=pk./shF;                  % /spike ht so neurons pool
        % ---- fig 311 (SS) / fig 312 (CS): bin the basal-vs-distal phase scatter into a 2x2
        % grid and accumulate the dendritic STA kymograph of each quadrant (pooled by distance) ----
        if et<=2
            Pb=ph(:,basalAx); if all(isnan(Pb)), Pb=ph(:,somaAx); end   % scatter x (soma proxy)
            Pd=ph(:,distalAx);                                          % scatter y
            for k=1:numel(ons)
                t0=ons(k);
                if t0-kymoTau(1)<1 || t0+kymoTau(2)>nTime, continue; end
                if ~isfinite(Pb(k)) || ~isfinite(Pd(k)), continue; end
                q = 1 + (Pb(k)>0) + 2*(Pd(k)<=0);   % 1 TL,2 TR,3 BL,4 BR (top row = distal>0)
                W = rawFull(:, t0+kymoT);                            % nROI x nKT raw window
                W = (W - mean(W(:,1:10),2,'omitnan'))./shF;          % baseline-subtract, /spike ht
                fin=isfinite(W); Wz=W; Wz(~fin)=0;
                meanB=(Mbin*Wz)./(Mbin*fin);                         % nKB x nKT distance-binned (NaN where empty)
                good=isfinite(meanB); meanB(~good)=0;
                KymoSum(:,:,q,et)=KymoSum(:,:,q,et)+meanB; KymoN(:,:,q,et)=KymoN(:,:,q,et)+good;
                nEvQuad(q,et)=nEvQuad(q,et)+1;
            end
        end
        % per-event soma color value for fig 309:
        %   CS -> peak somatic slow-depol amplitude; SS -> somatic post-spike AUC (0..50 ms)
        if et<=2
            cv=nan(numel(ons),1); aucp=nan(numel(ons),1); somaTr=SubBA(somaAx,:);
            for k=1:numel(ons)
                t0=ons(k);
                base=mean(somaTr(max(1,t0-20):t0-1),'omitnan');
                wA=t0+ssAUCwin; wA=wA(wA<=nTime); aucp(k)=sum(somaTr(wA)-base,'omitnan');   % post-spike soma AUC
                if et==2
                    w=t0+csAmpWin; w=w(w<=nTime); cv(k)=max(somaTr(w),[],'omitnan')-base;    % CS soma amplitude
                else
                    cv(k)=aucp(k);                                                            % SS post-spike AUC
                end
            end
            ColorVal_all{f,et}=cv; AUCpost_all{f,et}=aucp;
        end
        % compartment STA aligned to onset (SS / CS), baseline-subtracted, /SpikeHeight
        if et<=2 && ~isempty(ons)
            sta = get_STA(SubBA, ons, staTau(1), staTau(2));      % nComp x nT (subthreshold)
            if ~isempty(sta)
                sta = sta - mean(sta(:,1:10),2,'omitnan');        % baseline = first 10 ms of window
                STA_comp{f,et} = sta ./ PCresult{f}.SpikeHeight;  % normalize across neurons
            end
            staR = get_STA(rawComp, ons, staTau(1), staTau(2));   % nComp x nT (RAW, for fig-307 dV/dt)
            if ~isempty(staR)
                staR = staR - mean(staR(:,1:10),2,'omitnan');
                STA_comp_raw{f,et} = staR ./ PCresult{f}.SpikeHeight;
            end
            % phase-binned STA accumulation (fig 309 STA columns): per-event windows of the
            % SOMA raw trace, grouped by EACH compartment's pre-spike theta phase (row c),
            % and by the distal-basal phase difference (bottom row). /SpikeHeight so neurons pool.
            okWin = ons+(-staTau9(1))>=1 & ons+staTau9(2)<=nTime;   % fig-309 STA window (-150..50)
            onsK=ons(okWin); phK=ph(okWin,:);
            if ~isempty(onsK)
                [~, staMat] = get_STA(rawComp, onsK, staTau9(1), staTau9(2)); % nComp x nK x nT (RAW voltage)
                staMat = staMat ./ PCresult{f}.SpikeHeight;
                for c=1:nComp
                    bidx=discretize(phK(:,c),phaseEdges);               % bins from compartment c's phase
                    for bb=1:nPB
                        sel=bidx==bb; if ~any(sel), continue; end
                        wseg=reshape(staMat(somaAx,sel,:),sum(sel),[]);  % SOMA trace (nSel x nT)
                        wseg=wseg-mean(wseg(:,1:10),2,'omitnan');
                        STAbin_sum(c,bb,:,et)=STAbin_sum(c,bb,:,et)+reshape(sum(wseg,1,'omitnan'),1,1,[]);
                        STAbin_cnt(c,bb,et)=STAbin_cnt(c,bb,et)+sum(sel);
                    end
                end
                Pb=phK(:,basalAx); if all(isnan(Pb)), Pb=phK(:,somaAx); end   % req 4 proxy
                dpd=mod(phK(:,distalAx)-Pb+pi,2*pi)-pi;
                bidxD=discretize(dpd,phaseEdges);
                for bb=1:nPB
                    sel=bidxD==bb; if ~any(sel), continue; end
                    wseg=reshape(staMat(somaAx,sel,:),sum(sel),[]);
                    wseg=wseg-mean(wseg(:,1:10),2,'omitnan');
                    STAdiff_sum(bb,:,et)=STAdiff_sum(bb,:,et)+sum(wseg,1,'omitnan');
                    STAdiff_cnt(bb,et)=STAdiff_cnt(bb,et)+sum(sel);
                end
            end
        end
        % causal distal-basal theta phase match over the pre-window (soma proxies basal if absent)
        [plv,dph,okco]=get_prespike_coherence(SubBA(distalAx,:),basalTr,1000,thetaBand,ons,exclMask,phaseW_ms,phaseGuard_ms);
        plv(~okco)=NaN; dph(~okco)=NaN;
        Coh_all{f,et}=plv(:); Dphi_all{f,et}=dph(:);
    end

    % ---- phase-detection check: theta-FILTERED trace binned by the estimated phase (fig 303) ----
    % Pool all onset types; within each phase bin the filtered trace should stay a clean
    % oscillation (peaks shift diagonally with bin) and be at cos(phi) at t=0 if the phase
    % estimate is correct. Normalized to unit theta amplitude so neurons pool.
    for c=1:nComp
        [~, filtC] = get_phase(SubBA(c,:), 1000, thetaBand, exclMask);   % continuous band-passed
        sfc=std(filtC,0,'omitnan'); if ~(sfc>0), continue; end
        filtC=filtC/sfc;
        for et=1:nET
            ons=onsetSet{et}; if isempty(ons), continue; end
            phc=Phase_all{f,et}(:,c);
            okW = ons+(-verifTau(1))>=1 & ons+verifTau(2)<=nTime;
            ons2=ons(okW); phc2=phc(okW);
            if isempty(ons2), continue; end
            bidx=discretize(phc2,phaseEdges);
            for bb=1:nPB
                sel=bidx==bb; if ~any(sel), continue; end
                frB=ons2(sel); frB=frB(:);
                seg=filtC(frB+(-verifTau(1):verifTau(2)));      % nSel x nVT
                PV_sum(c,bb,:)=PV_sum(c,bb,:)+reshape(sum(seg,1,'omitnan'),1,1,[]);
                PV_cnt(c,bb)=PV_cnt(c,bb)+sum(sel);
            end
        end
    end
end

%% [3-verify] Theta phase-detection check: filtered traces binned by the estimated phase
% If the causal phase estimate is correct, then binning events by their estimated
% pre-spike phase and averaging the theta-FILTERED trace should give clean oscillations
% whose peak shifts diagonally across bins, and whose value at t=0 tracks cos(phi).
figure(303); clf; tiledlayout(2,nComp,'TileSpacing','compact','Padding','compact');
tV=-verifTau(1):verifTau(2); cmapPhase=hsv(nPB);
for c=1:nComp                                        % row 1: binned averaged filtered traces
    nexttile; hold on; box off;
    for bb=1:nPB
        if PV_cnt(c,bb)<5, continue; end
        plot(tV, squeeze(PV_sum(c,bb,:))/PV_cnt(c,bb), 'color',cmapPhase(bb,:),'LineWidth',1);
    end
    xline(0,'k:'); title(compName{c}); xlabel('Time from onset (ms)'); xlim([-150 50]);
    colormap(gca,hsv(nPB)); clim([-pi pi]);
    if c==1, ylabel('theta-filtered (norm.)'); end
    if c==nComp, cb=colorbar; cb.Label.String='estimated phase @ onset'; cb.Ticks=[-pi 0 pi]; cb.TickLabels={'-\pi','0','\pi'}; end
end
for c=1:nComp                                        % row 2: value at t=0 vs bin (should track cos)
    nexttile; hold on; box off;
    val0=nan(1,nPB);
    for bb=1:nPB
        if PV_cnt(c,bb)>=5, val0(bb)=PV_sum(c,bb,verifTau(1)+1)/PV_cnt(c,bb); end
    end
    hd=plot(phaseCtr, val0, 'ko-','MarkerFaceColor','k');
    A=max(abs(val0)); if ~(A>0), A=1; end
    hc=plot(phaseCtr, A*cos(phaseCtr), 'r--','LineWidth',1);
    xline(0,'k:'); yline(0,'k:'); xlim([-pi pi]);
    set(gca,'xtick',[-pi 0 pi],'xticklabel',{'-\pi','0','\pi'}); xlabel('estimated phase bin');
    if c==1, ylabel('filtered value @ t=0'); legend([hd hc],{'data','cos(\phi)'},'Box','off','FontSize',7); end
end
set_fontsize(9);
sgtitle('[3-verify] Phase-detection check: filtered trace binned by causal phase estimate (value@0 should track cos \phi)');

%% [3a] Figure 304 - pre-spike theta phase, compartment (rows) x event (cols)
figure(304); clf; tiledlayout(nComp,nRealET,'TileSpacing','compact','Padding','compact');
phEdges=linspace(-pi,pi,19);
for c=1:nComp
    for et=1:nRealET
        allPhi = cell2mat(cellfun(@(M) M(:,c), Phase_all(foi,et),'UniformOutput',false));
        allPhi = allPhi(~isnan(allPhi));
        nexttile;
        polarhistogram(allPhi, phEdges,'Normalization','probability','FaceColor',cmap_et(et,:),'EdgeColor','none'); hold on;
        if numel(allPhi)>=3
            R=circ_r(allPhi); mu=circ_mean(allPhi); p=circ_rtest(allPhi);
            polarplot([mu mu],[0 max(0.001,R)],'k','LineWidth',1.5);
            title(sprintf('%s / %s  R=%.2f p=%.1g', compName{c}, etypeName{et}, R, p),'FontSize',8);
        else
            title(sprintf('%s / %s  (n<3)', compName{c}, etypeName{et}),'FontSize',8);
        end
    end
end
set_fontsize(9);
sgtitle('[3a] Pre-spike theta phase: compartment (rows) x event (cols)');

%% [3b] Figure 306 - distal-basal theta phase match & basal-vs-distal phase per spike type
figure(306); clf; tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on; box off;                    % (a) phase-match per spike type -- SAME data as panel (c)
% cos(phi_distal - phi_basal) from the causal per-onset phase estimates (Phase_all), on the
% isolated subset, so this histogram is the 1-D marginal of the panel-(c) scatter.
pmEdges=linspace(-1,1,21); pmCtr=(pmEdges(1:end-1)+pmEdges(2:end))/2;
for et=1:nRealET
    pm=[];
    for f=foi
        P=Phase_all{f,et}; if isempty(P), continue; end
        Pb=P(:,basalAx); if all(isnan(Pb)), Pb=P(:,somaAx); end   % req 4: soma proxy
        pm=[pm; cos(P(:,distalAx)-Pb)];   %#ok<AGROW>
    end
    pm=pm(isfinite(pm));
    if isempty(pm), continue; end
    hc=histcounts(pm,pmEdges,'Normalization','probability');
    plot(pmCtr, hc, '-', 'color',cmap_et(et,:), 'LineWidth',1.5);
end
xlabel('Distal-basal phase match cos(\Delta\phi) at onset'); ylabel('P(event)');
xlim([-1 1]); legend(etypeName(1:nRealET),'Box','off','Location','best'); title('Phase-match at onset per spike type (isolated)');
nexttile; hold on; box off;                    % (b) pre-onset phase match per type vs silent
PLVneuron=nan(max(foi),nET);
for f=foi
    for et=1:nET
        if ~isempty(Coh_all{f,et}), PLVneuron(f,et)=mean(Coh_all{f,et},'omitnan'); end
    end
end
p=Boxplot_wPoints2(PLVneuron(foi,:), cmap_et); drawPValueLines(p,0); box off;
set(gca,'xtick',1:nET,'xticklabel',etypeName); ylabel('Pre-onset phase match'); title('Phase match at onset');
% (c) basal (or soma proxy) phase vs distal phase -- ONE PANEL PER SPIKE TYPE (SS,CS,BS,dSpike)
for et=1:nRealET
    nexttile; hold on; box off;
    bx=[]; dy=[];
    for f=foi
        P=Phase_all{f,et}; if isempty(P), continue; end
        Pb=P(:,basalAx); if all(isnan(Pb)), Pb=P(:,somaAx); end   % req 4: soma proxy
        bx=[bx; Pb]; dy=[dy; P(:,distalAx)];   %#ok<AGROW>
    end
    v=isfinite(bx)&isfinite(dy);
    if nnz(v)>=10
        scatter_density(bx(v), dy(v), 10); colorbar(gca,'off');   % color = local point density
    elseif any(v)
        scatter(bx(v), dy(v), 8, cmap_et(et,:),'filled','MarkerFaceAlpha',0.5);
    end
    plot([-pi pi],[-pi pi],'k:');               % diagonal = in-phase (matched)
    axis square; xlim([-pi pi]); ylim([-pi pi]);
    set(gca,'xtick',[-pi 0 pi],'xticklabel',{'-\pi','0','\pi'},'ytick',[-pi 0 pi],'yticklabel',{'-\pi','0','\pi'});
    xlabel('basal phase'); ylabel('distal phase');
    title(sprintf('%s: basal vs distal phase (n=%d, color=density)', etypeName{et}, sum(v)));
end
set_fontsize(9);
sgtitle('[3b] Distal-basal phase match (<cos \Delta\phi>) & joint basal-distal phase per spike type');

%% [3c] Figure 307 - pre-spike dV/dt by compartment (slope of the RAW STA), SS vs CS
figure(307); clf; tiledlayout(1,nComp,'TileSpacing','compact','Padding','compact');
preLbl=sprintf('%d..%d ms', preWin(1), preWin(end));
staCols=preWin+staTau(1)+1;                          % STA columns spanning preWin (tSTA(col)=time)
for c=1:nComp
    M=nan(max(foi),2);
    for f=foi
        for et=1:2                                   % dV/dt = slope of the raw STA over preWin
            S=STA_comp_raw{f,et};
            if isempty(S), continue; end
            seg=S(c,staCols); gd=isfinite(seg);
            if nnz(gd)>=3, pp=polyfit(preWin(gd),seg(gd),1); M(f,et)=pp(1); end
        end
    end
    nexttile;
    p=Boxplot_wPoints2(M(foi,:), cmap_et(1:2,:)); drawPValueLines(p,0); box off;
    set(gca,'xtick',[1 2],'xticklabel',{'SS','CS'}); title(compName{c});
    if c==1, ylabel(sprintf('STA dV/dt, %s (\\DeltaF/F/spike ht per ms)',preLbl)); end
end
set_fontsize(11);
sgtitle(sprintf('[3c] Pre-spike (%s) dV/dt from the RAW STA, by compartment: SS vs CS', preLbl));

%% [3d] Figure 308 - compartment STA + pre-spike V-vs-dV/dt phase plane
% Top STA panels are the SUBTHRESHOLD compartment STA. The bottom phase plane uses the
% RAW traces (per event): dV/dt is the linear-fit slope over preWin, V is the mean raw
% over the same window; both are divided by each neuron's spike height so events pool.
figure(308); clf; tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
tSTA=-staTau(1):staTau(2); pw0=preWin(1); pw1=preWin(end);
dispOrd=[1 4 2 3];                                  % display order: Basal, Soma, Apical, Distal
cmap_comp=[0.2 0.4 1; 0 0 0; 0.2 0.7 0.2; 0.85 0.3 0.1];
compLbl=compName(dispOrd);
stypeOrder=[2 1];                                   % top panels: CS then SS
% --- top row: compartment STA per event (each spans 2 tiles) ---
for pp=1:2
    stype=stypeOrder(pp);
    nexttile([1 2]); hold on; box off;
    h=gobjects(1,numel(dispOrd));
    for j=1:numel(dispOrd)
        c=dispOrd(j);
        M = cell2mat(cellfun(@(S) S(c,:), STA_comp(foi,stype),'UniformOutput',false));  % neurons x nT
        M = M(all(isfinite(M),2),:);
        if isempty(M), continue; end
        h(j)=errorbar_shade(tSTA, mean(M,1,'omitnan'), std(M,0,1,'omitnan')./sqrt(size(M,1)), cmap_comp(j,:));
    end
    yl=ylim; hp=patch([pw0 pw1 pw1 pw0],[yl(1) yl(1) yl(2) yl(2)],[0.9 0.9 0.6], ...
        'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off');   % pre-spike dV/dt window
    uistack(hp,'bottom'); xline(0,'k:'); title([etypeName{stype} ' STA']); xlabel('Time from onset (ms)');
    if pp==1
        ylabel('Subthreshold (\DeltaF/F / spike ht)');
        valid=isgraphics(h); legend(h(valid), compLbl(valid),'Box','off','Location','northwest');
    end
end
% --- bottom row: phase plane V vs dV/dt (one point per event), SS & CS, per compartment ---
for j=1:numel(dispOrd)
    c=dispOrd(j);
    nexttile; hold on; box off;
    hS=gobjects(1,2);
    for stype=1:2
        Vv=[]; dV=[];
        for f=foi
            if isempty(Vpre_all{f,stype}), continue; end
            sh=PCresult{f}.SpikeHeight;
            Vv=[Vv; Vpre_all{f,stype}(:,c)./sh]; dV=[dV; dVdtPre_all{f,stype}(:,c)./sh];   %#ok<AGROW>
        end
        v=isfinite(Vv)&isfinite(dV); Vv=Vv(v); dV=dV(v);
        if ~isempty(Vv)
            hS(stype)=scatter(Vv, dV, 8, cmap_et(stype,:),'filled','MarkerFaceAlpha',0.25);
            % 2D mean +- std errorbar (both axes linear)
            errorbar(mean(Vv), mean(dV), std(dV), std(dV), std(Vv), std(Vv), 'o', ...
                'Color','k','MarkerFaceColor',cmap_et(stype,:),'MarkerSize',7,'LineWidth',1.5,'CapSize',6);
        end
    end
    yline(0,'k:'); xline(0,'k:'); title(compLbl{j});
    if j==1
        xlabel('V (\DeltaF/F / spike ht)'); ylabel('dV/dt (per ms)');
        lbl={'SS','CS'}; valid=isgraphics(hS); legend(hS(valid),lbl(valid),'Box','off','Location','best');
    end
end
set_fontsize(10);
sgtitle(sprintf('[3d] Top: subthreshold STA (CS, SS).  Bottom: pre-spike (%d..%d ms) V vs dV/dt (RAW)', pw0, pw1));

%% [3e] Figure 309 - phase x amplitude scatter (cols 1-2) + SOMA STA binned by phase (cols 3-4)
% Columns: 1=SS scatter, 2=CS scatter, 3=SS soma STA, 4=CS soma STA.
% Rows 1-4 = compartments. Scatter: phase x, theta amplitude y, color=post-spike soma AUC.
% STA columns: the SOMA raw-voltage STA, grouped by bins of THAT ROW's compartment phase
% (so you see how the somatic waveform depends on each compartment's theta phase). Row 5 =
% distal-basal phase DIFFERENCE (scatter: x=phase diff, y=joint distal*basal amplitude;
% STA: soma STA grouped by the phase difference). Soma proxies basal where absent (req 4).
figure(309); clf; tiledlayout(nComp+1,4,'TileSpacing','compact','Padding','compact');
ampScale=nan(max(foi),1);                            % per-neuron amp scale (shared SS & CS)
for f=foi, ampScale(f)=median([Mag_all{f,1}(:); Mag_all{f,2}(:)],'omitnan'); end
phaseOffsets=[-2*pi 0 2*pi];                         % tile phase so -3pi/2..3pi/2 is continuous
wrappi=@(x) mod(x+pi,2*pi)-pi;
tSTA=-staTau9(1):staTau9(2); cmapPhase=hsv(nPB);     % fig-309 STA time axis (-150..50 ms)
for c=1:nComp
    % --- cols 1-2: SS / CS scatter, phase(x, 3 cycles) x amplitude(y), color=post AUC ---
    for stype=1:2
        ph=[]; am=[]; cv=[];
        for f=foi
            P=Phase_all{f,stype}; A=Mag_all{f,stype}; C=AUCpost_all{f,stype};
            if isempty(P) || ~isfinite(ampScale(f)) || ampScale(f)<=0, continue; end
            ph=[ph; P(:,c)]; am=[am; A(:,c)./ampScale(f)]; cv=[cv; C(:)];   %#ok<AGROW>
        end
        v=isfinite(ph)&isfinite(am)&isfinite(cv); ph=ph(v); am=am(v); cv=cv(v);
        nexttile; hold on; box off;
        yyaxis left;                                       % left axis: phase x amplitude scatter
        for oo=phaseOffsets, scatter(ph+oo, am, 8, cv, 'filled','MarkerFaceAlpha',0.5); end
        colormap(gca, parula(256));
        if ~isempty(cv), cr=prctile(cv,[5 95]); if diff(cr)>0, clim(cr); end, end
        if stype==1, ylabel(sprintf('%s\ntheta amp', compName{c})); end
        set(gca,'YColor','k');
        yyaxis right;                                      % right axis: # events vs phase (marginal)
        if ~isempty(ph)
            phEdgesM=linspace(-pi,pi,19); phCtrM=(phEdgesM(1:end-1)+phEdgesM(2:end))/2;
            cntM=histcounts(ph, phEdgesM);
            for oo=phaseOffsets, plot(phCtrM+oo, cntM, '-','color',[0 0 0],'LineWidth',1.3); end
            ylim([0 max(cntM)*3.2]);                        % keep the count line in the lower third
        end
        set(gca,'YColor',[0.3 0.3 0.3]); if stype==2, ylabel('# events'); end
        yyaxis left;
        cb=colorbar; cb.Label.String='post AUC';
        xlim([-1.5*pi 1.5*pi]);
        set(gca,'xtick',(-1:1)*pi,'xticklabel',{'-\pi','0','\pi'});
        if c==nComp, xlabel('theta phase'); end
        title(sprintf('%s  n=%d (all neurons)', etypeName{stype}, numel(ph)), 'FontSize',7);
    end
    % --- cols 3-4: SS / CS SOMA raw-voltage STA, binned by THIS ROW's compartment phase ---
    for stype=1:2
        nexttile; hold on; box off;
        for bb=1:nPB
            if STAbin_cnt(c,bb,stype)<3, continue; end
            m=squeeze(STAbin_sum(c,bb,:,stype))/STAbin_cnt(c,bb,stype);
            plot(tSTA, m, 'color',cmapPhase(bb,:),'LineWidth',1);
        end
        xline(0,'k:'); xlim([tSTA(1) tSTA(end)]); colormap(gca,hsv(nPB)); clim([-pi pi]);
        if c==1, cb=colorbar; cb.Label.String='pre-spike phase'; cb.Ticks=[-pi 0 pi]; cb.TickLabels={'-\pi','0','\pi'}; end
        if c==nComp, xlabel('Time from onset (ms)'); end
        ylabel('soma raw V (norm.)');
        title(sprintf('%s soma STA / %s phase', etypeName{stype}, compName{c}));
    end
end
% --- row 5: distal-basal phase DIFFERENCE (x) x offset-strength (y); soma STA by phase diff ---
for stype=1:2
    dpb=[]; amO=[]; cv=[];
    for f=foi
        P=Phase_all{f,stype}; A=Mag_all{f,stype}; C=AUCpost_all{f,stype};
        if isempty(P) || ~isfinite(ampScale(f)) || ampScale(f)<=0, continue; end
        Pb=P(:,basalAx); Ab=A(:,basalAx);
        if all(isnan(Pb)), Pb=P(:,somaAx); Ab=A(:,somaAx); end     % req 4 proxy
        dpb=[dpb; wrappi(P(:,distalAx)-Pb)]; amO=[amO; sqrt(A(:,distalAx).*Ab)./ampScale(f)]; cv=[cv; C(:)];   %#ok<AGROW>
    end
    v=isfinite(dpb)&isfinite(amO)&isfinite(cv); dpb=dpb(v); amO=amO(v); cv=cv(v);
    nexttile; hold on; box off;
    yyaxis left;                                       % left axis: phase-diff x offset-strength scatter
    if ~isempty(cv)
        scatter(dpb, amO, 8, cv, 'filled','MarkerFaceAlpha',0.5); colormap(gca,parula(256));
        cr=prctile(cv,[5 95]); if diff(cr)>0, clim(cr); end
    end
    ylabel('offset strength (rel. amp)'); set(gca,'YColor','k');
    yyaxis right;                                      % right axis: # events vs phase difference
    if ~isempty(dpb)
        phEdgesD=linspace(-pi,pi,19); phCtrD=(phEdgesD(1:end-1)+phEdgesD(2:end))/2;
        cntD=histcounts(dpb, phEdgesD);
        plot(phCtrD, cntD, '-','color',[0 0 0],'LineWidth',1.3);
        ylim([0 max(cntD)*3.2]);                        % keep the count line in the lower third
    end
    set(gca,'YColor',[0.3 0.3 0.3]); if stype==2, ylabel('# events'); end
    yyaxis left;
    cb=colorbar; cb.Label.String='post AUC'; xline(0,'k:'); xlim([-pi pi]);
    set(gca,'xtick',[-pi 0 pi],'xticklabel',{'-\pi','0','\pi'});
    xlabel('\phi_{distal}-\phi_{basal}');
    title(sprintf('%s offset  n=%d', etypeName{stype}, numel(dpb)));
end
for stype=1:2                                         % soma STA binned by distal-basal phase diff
    nexttile; hold on; box off;
    for bb=1:nPB
        if STAdiff_cnt(bb,stype)<3, continue; end
        m=squeeze(STAdiff_sum(bb,:,stype))/STAdiff_cnt(bb,stype);
        plot(tSTA, m, 'color',cmapPhase(bb,:),'LineWidth',1);
    end
    xline(0,'k:'); xlim([tSTA(1) tSTA(end)]); colormap(gca,hsv(nPB)); clim([-pi pi]);
    cb=colorbar; cb.Label.String='\phi_{distal}-\phi_{basal}'; cb.Ticks=[-pi 0 pi]; cb.TickLabels={'-\pi','0','\pi'};
    xlabel('Time from onset (ms)'); title([etypeName{stype} ' soma raw STA / \Delta\phi']);
end
set_fontsize(8);
sgtitle('[3e] Cols 1-2: phase x amplitude scatter (n=all-neuron events; color=post AUC).  Cols 3-4: RAW-voltage STA binned by phase / phase diff');

%% [3e+] Figures 311 (SS) & 312 (CS) - dendritic STA kymograph, binned 2x2 by (basal phase) x (distal phase)
% Companion to fig 309/306: the isolated basal-vs-distal phase scatter is split at phase 0 into
% four quadrants; within each quadrant the dendritic STA (RAW, baseline-subtracted per event,
% /spike ht) is pooled across neurons by dendritic distance and shown as a kymograph.
% Tiles match the scatter layout: columns = basal phase (<0 | >0), rows = distal phase (>0 | <0).
qTitle={'basal<0, distal>0','basal>0, distal>0','basal<0, distal<0','basal>0, distal<0'};
kymoFig=[311 312]; kymoName={'SS','CS'};             % last dim of KymoSum: 1=SS, 2=CS
for et=1:2
    figure(kymoFig(et)); clf; tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    Kymo=KymoSum(:,:,:,et)./KymoN(:,:,:,et); Kymo(KymoN(:,:,:,et)==0)=NaN;   % mean kymograph per quadrant
    finK=Kymo(isfinite(Kymo));
    if isempty(finK), cax=[-1 1]; else, cax=[-1 1]*max(prctile(abs(finK),99), eps); end   % symmetric color scale
    for q=1:4
        nexttile;
        imagesc(kymoT, kymoCtr, Kymo(:,:,q), cax); axis xy; hold on;
        xline(0,'k:'); yline(0,'k:');                % t=0 (onset) & soma (distance 0)
        title(sprintf('%s (n=%d)', qTitle{q}, nEvQuad(q,et)));
        if q>=3, xlabel(sprintf('Time from %s onset (ms)', kymoName{et})); end
        if mod(q,2)==1, ylabel('dendritic distance (\mum)'); end
    end
    colormap(turbo); cb=colorbar; cb.Label.String='\DeltaF/F / spike ht'; cb.Layout.Tile='east';
    set_fontsize(10);
    sgtitle(sprintf('[3e+] %s dendritic STA kymograph, binned 2x2 by (basal \\phi) x (distal \\phi) at onset', kymoName{et}));
end

%% [3g] Figure 310 - CS peak slow component as a function of theta phase
% CS events only. For each CS, take the PEAK of the slow (subthreshold) depolarization in the
% 0..150 ms window after the first spike (PeakSlow_all, /spike ht) and relate it to the pre-spike
% theta phase. y = peak slow component of the DISTAL dendrite and the SOMA; x = theta phase,
% shown against five phase references (basal / apical / soma / distal, and the distal-basal
% phase difference). Peak and phase are matched per event (same onsets => same rows).
figure(310); clf; tiledlayout(1,5,'TileSpacing','compact','Padding','compact');
wrappi=@(x) mod(x+pi,2*pi)-pi; etCS=2;         % CS only
PH=[]; PK=[];
for f=foi
    P=Phase_all{f,etCS}; K=PeakSlow_all{f,etCS};
    if isempty(P) || isempty(K), continue; end
    Pb=P(:,basalAx); if all(isnan(Pb)), Pb=P(:,somaAx); end          % soma proxy for the diff
    PH=[PH; P(:,basalAx) P(:,apicalAx) P(:,somaAx) P(:,distalAx) wrappi(P(:,distalAx)-Pb)]; %#ok<AGROW>
    PK=[PK; K(:,distalAx) K(:,somaAx)];        %#ok<AGROW>  distal, soma peak slow component
end
refName={'basal \phi','apical \phi','soma \phi','distal \phi','\phi_{distal}-\phi_{basal}'};
phB=linspace(-pi,pi,11); phC=(phB(1:end-1)+phB(2:end))/2;      % 10 phase bins
cmapDS=[0.85 0.3 0.1; 0 0 0];                                  % distal (red), soma (black)
for r=1:5
    nexttile; hold on; box off; hl=gobjects(1,2);
    for j=1:2
        [mu,se]=phase_tuning(PH(:,r), PK(:,j), phB);
        hl(j)=errorbar_shade(phC, mu, se, cmapDS(j,:));
    end
    xline(0,'k:'); xlim([-pi pi]); set(gca,'xtick',[-pi 0 pi],'xticklabel',{'-\pi','0','\pi'});
    xlabel(refName{r}); title(refName{r});
    if r==1
        ylabel('CS peak slow component (/spike ht, 0..150 ms)');
        vv=isgraphics(hl); legend(hl(vv),{'distal','soma'},'Box','off','Location','best');
    end
end
set_fontsize(9);
sgtitle(sprintf('[3g] CS peak slow component (distal & soma, %d..%d ms) vs pre-spike theta phase', peakSlowWin(1), peakSlowWin(end)));

%% [3f] Figure 305 - CS slow depolarization: soma-vs-dendrite timing & dV/dt vs bAP order  [Task 2]
% (1) Cross-correlate the soma and distal SLOW (spike-removed, subthreshold)
%     depolarization of each CS to measure their timing difference (lag of the
%     xcorr peak). Positive lag => distal LAGS soma (i.e. soma leads distal).
% (2) Relate the dendritic slow-depolarization rate (dV/dt of apical / distal)
%     to the order of the fast somatic bAP within the CS burst.
if ~exist('nTau_EV','var'), nTau_EV=[-200 500]; end
OrthAxis=[1 0 0 0 0; 0 0 1 0 1; 0 0 0 0 1; 0 1 0 0 0]; % rows: basal, apical, distal, soma
apicalAxis=2; distalAxis=3; somaAxis=4;
tWin_ms=[-50 300]; maxLagCS=60;      % xcorr analysis window & half-window (ms)
lagStepFine=0.05;                    % sub-ms lag resolution via spline interp (req 5)
onsetCol=-nTau_EV(1)+1;              % column of the event onset in AlignedEvntall
winCols=onsetCol+[tWin_ms(1):tWin_ms(2)];
baseCols=1:20;                       % baseline = first 20 ms of the analysis window
maxOrder=8; smoothMs=5; dVdtWin=2;   % soma-spike order cap; dV/dt smoothing & +-window (ms)
maxSearch=0:200;                     % post-onset window (ms) to locate the peak distal dV/dt & V
spkPre=1; spkPost=3; slowSm=5;       % distal SLOW build: blank spike frames [-spkPre..+spkPost], interp, movmean(slowSm)

CSlag=nan(max(foi),1); CSlag_perEv=cell(max(foi),1);
dVdt_distal=nan(max(foi),maxOrder);     dVdt_apical=nan(max(foi),maxOrder);       % at spike
CSsomaGrand=[]; CSdistGrand=[];      % pooled (peak-normalized) mean CS waveforms for display
% per-CS-event distal-dendrite peak timing (relative to first spike / CS onset)
CSt_maxdVdt=[]; CSt_maxV=[]; CSt_rise=[]; CSdVdt_val=[]; CSV_val=[]; CSorder_at=[];
% per-event feature records for the interactive browser (fig 3052)
exWin=-50:250;                       % kymograph display window (ms rel. first spike)
CSmeta=struct('f',{},'ev',{},'feat',{},'scalar',{},'dist',{},'somaSlow',{},'spk',{},'label',{});
for f=foi
    [nROI, nTime]=size(PCresult{f}.Subthreshold);
    OrthAxis_vec=build_orthaxis_vec(OrthAxis, PCresult{f}.roisD_order, nROI);
    if any(isnan(OrthAxis_vec([somaAxis distalAxis apicalAxis],1))), continue; end

    % ---- (1) per-CS soma vs distal slow-component cross-correlation ----
    A3 = AlignedEvntall{f,3};        % ROI x time x event, subthreshold (fast spikes removed)
    A1 = AlignedEvntall{f,1};        % ROI x time x event, RAW (for the spike-excised slow build)
    Et = EventTable(EventTable.Neuron==f,:);
    csIdx = Et.LocalIdx(Et.EventClass=="CS" & Et.IsValid);
    somaEv=[]; distEv=[]; lags_ev=[];
    for ii=1:numel(csIdx)
        slow = OrthAxis_vec * A3(:,:,csIdx(ii));      % 4 x time
        s  = slow(somaAxis, winCols);  dd = slow(distalAxis, winCols);
        s  = s  - mean(s(baseCols),'omitnan');
        dd = dd - mean(dd(baseCols),'omitnan');
        if all(isnan(s)) || all(isnan(dd)), continue; end
        [cc,lg]=nanXCorr(s, dd, maxLagCS, true);
        lagThis=subms_peaklag(cc, lg, lagStepFine);   % sub-ms soma->distal lag for THIS CS (req 5)
        lags_ev(end+1,1)=lagThis;   %#ok<SAGROW>
        somaEv=[somaEv; s]; distEv=[distEv; dd];        %#ok<AGROW>

        % ---- peak distal dV/dt & V timing relative to CS onset (first spike) ----
        % Build DISTAL and SOMA slow traces from the RAW compartments: blank frames
        % [-spkPre..+spkPost] around every somatic spike, interpolate, then movmean(slowSm) --
        % tracks the depolarization rise without the heavy subthreshold smoothing.
        onF  = Et.Frame(Et.LocalIdx==csIdx(ii));
        stev = SpikeTable(SpikeTable.Neuron==f & SpikeTable.EventLocalIdx==csIdx(ii),:);
        nAl  = size(A1,2); spCols = onsetCol + (stev.Frame - onF)';   % somatic spike columns
        bad  = false(1,nAl);
        for sc = spCols(:)', bad(max(1,sc-spkPre):min(nAl,sc+spkPost))=true; end
        rawD = OrthAxis_vec(distalAxis,:) * A1(:,:,csIdx(ii)); rawD(bad)=NaN;   % RAW distal
        rawS = OrthAxis_vec(somaAxis,:)   * A1(:,:,csIdx(ii)); rawS(bad)=NaN;   % RAW soma
        ddF = movmean(interpolateNaN(rawD), slowSm, 'omitnan');   % cleaned distal slow
        sdF = movmean(interpolateNaN(rawS), slowSm, 'omitnan');   % cleaned soma slow
        % restrict the peak search to WITHIN this CS: stop before the next somatic event onset
        maxT = maxSearch(end); nextOn = min(Et.Frame(Et.IsValid & Et.Frame>onF));
        if ~isempty(nextOn), maxT = min(maxT, (nextOn-onF)-1); end
        srch = 0:max(0,maxT); sCols = onsetCol+srch;
        keep = sCols>=2 & sCols<=nAl; sCols=sCols(keep); srch=srch(keep);
        if numel(sCols)>10 && any(isfinite(ddF(sCols)))
            base = mean(ddF(max(1,onsetCol-20):onsetCol-1),'omitnan');
            baseS= mean(sdF(max(1,onsetCol-20):onsetCol-1),'omitnan');
            Vseg = ddF(sCols)-base; dvseg = ddF(sCols)-ddF(sCols-1);   % V and dV/dt over search window
            [~,iV]=max(Vseg); tVv=srch(iV);                            % peak distal V (within this CS)
            [~,iR]=min(Vseg(1:iV)); tRise=srch(iR);                    % rise onset = trough before the peak
            seg=iR:iV; [mdv,rel]=max(dvseg(seg)); iD=seg(rel); tDd=srch(iD);   % max dV/dt LEADING to the peak
            ordAt = sum((stev.Frame-onF)>=0 & (stev.Frame-onF)<=tDd);  % # bAPs by the peak-dV/dt time
            CSt_maxdVdt(end+1,1)=tDd; CSt_maxV(end+1,1)=tVv; CSt_rise(end+1,1)=tRise;   %#ok<AGROW>
            CSdVdt_val(end+1,1)=mdv; CSV_val(end+1,1)=Vseg(iD); CSorder_at(end+1,1)=ordAt;   %#ok<AGROW>
            % record this event for the interactive browser (fig 3052)
            ec=onsetCol+exWin;
            if all(ec>=1 & ec<=nAl)
                spk   =spikes_in_window(PCresult{f}.SpikeClassVecMat, onF, exWin, nTime);   % SS/CS/BS/dSpike times
                feat  =struct('riseOnset',tRise,'maxdVdt',tDd,'maxV',tVv);   % TIME markers (drawn on panels)
                scalar=struct('riseOnset',tRise,'maxdVdt',tDd,'maxV',tVv,'lag',lagThis);   % scatter-axis scalars
                CSmeta(end+1)=struct('f',f,'ev',csIdx(ii),'feat',feat,'scalar',scalar, ...
                    'dist',ddF(ec)-base,'somaSlow',sdF(ec)-baseS,'spk',spk, ...
                    'label',sprintf('f%d ev%d',f,csIdx(ii)));   %#ok<SAGROW>
            end
        end
    end
    CSlag_perEv{f}=lags_ev;
    if size(somaEv,1)>=3
        sM=mean(somaEv,1,'omitnan'); dM=mean(distEv,1,'omitnan');
        [cc,lg]=nanXCorr(sM, dM, maxLagCS, true);
        CSlag(f)=subms_peaklag(cc, lg, lagStepFine);           % sub-ms lag on mean CS waveform
        if max(sM)>0 && max(dM)>0                        % peak-normalize for shape overlay
            CSsomaGrand=[CSsomaGrand; sM./max(sM)];      %#ok<AGROW>
            CSdistGrand=[CSdistGrand; dM./max(dM)];      %#ok<AGROW>
        end
    end

    % ---- (2) dendritic slow dV/dt vs fast somatic bAP order within the CS ----
    SubBA  = (PCresult{f}.Subthreshold' * OrthAxis_vec')';    % 4 x nTime slow components
    distTr = movmean(SubBA(distalAxis,:), smoothMs, 'omitnan');
    apicTr = movmean(SubBA(apicalAxis,:), smoothMs, 'omitnan');
    dV_dist=[NaN diff(distTr)];  dV_apic=[NaN diff(apicTr)];  % dF/F per ms (fs=1000 Hz)
    St = SpikeTable(SpikeTable.Neuron==f & SpikeTable.EventClass=="CS" & SpikeTable.IsValid,:);
    for o=1:maxOrder
        fr = St.Frame(St.SpikeOrder==o);
        fr = fr(fr>dVdtWin & fr<=nTime-dVdtWin);
        if isempty(fr), continue; end
        % dV/dt at the somatic spike (+-dVdtWin ms)
        dVdt_distal(f,o)=mean(arrayfun(@(x) mean(dV_dist(x+[-dVdtWin:dVdtWin]),'omitnan'), fr),'omitnan');
        dVdt_apical(f,o)=mean(arrayfun(@(x) mean(dV_apic(x+[-dVdtWin:dVdtWin]),'omitnan'), fr),'omitnan');
    end
end

% ---- figure ----
figure(305); clf; tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
tAx=tWin_ms(1):tWin_ms(2);
nexttile; hold on; box off;                              % (a) grand-avg CS slow soma vs distal
errorbar_shade(tAx, mean(CSsomaGrand,1,'omitnan'), std(CSsomaGrand,0,1,'omitnan')./sqrt(max(1,size(CSsomaGrand,1))), [0 0 0]);
errorbar_shade(tAx, mean(CSdistGrand,1,'omitnan'), std(CSdistGrand,0,1,'omitnan')./sqrt(max(1,size(CSdistGrand,1))), [0.85 0.3 0.1]);
xline(0,'k:'); xlabel('Time from CS onset (ms)'); ylabel('Slow \DeltaF/F (peak-norm.)');
legend({'Soma','Distal'},'Box','off','Location','southeast'); title('CS slow depolarization');
nexttile; hold on; box off;                              % (b) soma-vs-distal lag: EACH event
allLags=cell2mat(CSlag_perEv(foi));
if ~isempty(allLags)
    histogram(allLags, -maxLagCS:2:maxLagCS, 'FaceColor',[0.2 0.4 0.8],'EdgeColor','none','Normalization','probability');
    xline(median(allLags,'omitnan'),'r-','LineWidth',1.5); xline(0,'k:');
end
xlabel('xcorr lag soma\rightarrowdistal (ms)'); ylabel('P(CS event)');
title(sprintf('Each CS (n=%d); - = distal precedes soma (median %.1f)', numel(allLags), median(allLags,'omitnan')));
nexttile; hold on; box off;                              % (c) dendritic dV/dt vs bAP order (@spike)
ord=1:maxOrder; sem=@(M) std(M,0,1,'omitnan')./sqrt(max(1,sum(~isnan(M),1)));
errorbar_shade(ord, mean(dVdt_distal(foi,:),1,'omitnan'), sem(dVdt_distal(foi,:)), [0.85 0.3 0.1]);
errorbar_shade(ord, mean(dVdt_apical(foi,:),1,'omitnan'), sem(dVdt_apical(foi,:)), [0.2 0.5 0.2]);
xlabel('Fast somatic bAP order within CS'); ylabel('Slow dV/dt (\DeltaF/F per ms)');
legend({'Distal','Apical'},'Box','off','Location','best'); title('Dendritic dV/dt vs bAP order'); xlim([0.5 maxOrder+0.5]);
nexttile; hold on; box off;                              % (d) timing of rise onset / peak dV/dt / peak V vs 1st spike
tEdges=maxSearch(1):10:maxSearch(end);
histogram(CSt_rise,    tEdges,'FaceColor',[0.2 0.6 0.2], 'EdgeColor','none','FaceAlpha',0.5,'Normalization','probability');
histogram(CSt_maxdVdt, tEdges,'FaceColor',[0.85 0.3 0.1],'EdgeColor','none','FaceAlpha',0.5,'Normalization','probability');
histogram(CSt_maxV,    tEdges,'FaceColor',[0.2 0.3 0.8], 'EdgeColor','none','FaceAlpha',0.5,'Normalization','probability');
xlabel('Time from first spike (ms)'); ylabel('P(CS event)');
legend({'rise onset','peak distal dV/dt','peak distal V'},'Box','off','Location','best'); title('Peak dendritic timing');
nexttile; hold on; box off;                              % (e) at peak dV/dt: dV/dt(x) vs V(y), color=bAP order
if ~isempty(CSdVdt_val)
    scatter(CSdVdt_val, CSV_val, 14, CSorder_at, 'filled','MarkerFaceAlpha',0.6);
    colormap(gca, turbo(256)); cb=colorbar; cb.Label.String='# bAPs by peak dV/dt';
    if max(CSorder_at)>min(CSorder_at), clim([min(CSorder_at) max(CSorder_at)]); end
end
xlabel('peak distal dV/dt (\DeltaF/F per ms)'); ylabel('distal V at that time (\DeltaF/F)');
title('Distal phase-plane at peak dV/dt');
set_fontsize(11);
sgtitle('[3f] CS soma-vs-dendrite timing (each event, sub-ms) & distal peak dV/dt / V vs first spike');

%% [3f++] Figure 3052 - INTERACTIVE feature browser: click a point -> that event's raw kymograph + traces
% LEFT: scatter, x = peak dV/dt time, y = soma->distal slow-component lag (ms; + = distal lags
% soma). RIGHT (top->bottom): RAW dendritic kymograph (soma top -> distal bottom), the somatic
% slow-component trace, and the dendritic slow-component trace -- time-aligned. Every panel marks
% the TIME features (green=rise onset, red=peak dV/dt, blue=peak V) as dotted lines and all spikes
% in the window (SS=magenta, CS=teal, BS=purple, dSpike=orange). Generalizable: add fields to
% CSmeta(i).feat (time markers) and/or .scalar (scatter axes); pass any two scalar names as x/y.
if isempty(CSmeta)
    warning('fig 3052: CSmeta is empty; run the fig-305 loop first.');
else
    featColors=struct('riseOnset',[0.2 0.6 0.2], 'maxdVdt',[0.85 0.1 0.1], 'maxV',[0.1 0.3 0.85]);
    spkColors =struct('SS',[1 0 1], 'CS',[0 0.7 0.7], 'BS',[0.5 0.2 0.75], 'dSpike',[0.9 0.6 0]);
    getKymo=@(i) fetch_raw_kymo(CSmeta(i), AlignedEvntall, PCresult, onsetCol, exWin);
    interactive_feature_kymo(3052, CSmeta, featColors, 'maxdVdt', 'lag', getKymo, exWin, spkColors, ...
        '[3f++] Click a CS: kymo + soma/distal slow traces.  features green=rise red=dV/dt blue=V.  spikes SS=magenta CS=teal BS=purple dSpike=orange');
end

%% [3h] Figures 313 (pooled) & 314 (per neuron) - CS first-spike STA + spike-order timing
% Somatic (black) and distal-dendrite (red) STA (soma/distal-class ROIs of NormalizedTrace_dirt,
% /spike ht, pre-onset baseline-subtracted) aligned via get_STA to the FIRST spike of each complex
% spike (SpikeOrder==1 within a CS), EXCLUDING CS onsets during blue. Overlaid (right y-axis) are the timing histograms
% of the 1st, 2nd, 3rd ... somatic spikes within the CS relative to the first spike -- obtained by
% aligning each per-order spike train with the same get_STA call -- each order a different colour.
tWinSTA=-50:250; tauPre=50; tauPost=250; maxOrd=6; cmapOrd=turbo(maxOrd); binE=tWinSTA(1):1:tWinSTA(end);
cSoma=[0 0 0]; cDist=[0.85 0.3 0.1];                                                 % soma / distal STA colours
somaEvN=cell(1,numel(foi)); distEvN=cell(1,numel(foi)); ordCntN=cell(1,numel(foi)); nEvN=zeros(1,numel(foi));
somaEvAll=[]; distEvAll=[]; ordCntAll=zeros(maxOrd,numel(tWinSTA)); nEvAll=0;         % pooled
for i=1:numel(foi)
    f=foi(i);
    somaROI=ismember(PCresult{f}.roi_dClass,2);
    distROI=ismember(PCresult{f}.roi_dClass,5);                                 % distal dendrite class
    somaTr = mean(PCresult{f}.NormalizedTrace_dirt(somaROI,:),1,'omitnan') / PCresult{f}.SpikeHeight;
    distTr = mean(PCresult{f}.NormalizedTrace_dirt(distROI,:),1,'omitnan') / PCresult{f}.SpikeHeight;
    csSpk  = (PCresult{f}.SpikeMat(1,:)>0) & (PCresult{f}.CStrace>0);            % somatic spikes in a CS
    trig   = csSpk & (PCresult{f}.SpikeOrder==1) & ~(PCresult{f}.BlueStim>0);    % CS first spikes, no blue
    ordCntN{i}=zeros(maxOrd,numel(tWinSTA));
    if ~any(somaROI) || ~any(trig), continue; end
    [~, evSoma] = get_STA(somaTr, trig, tauPre, tauPost);        % 1 x nEv x W soma windows
    evSoma = squeeze(evSoma); if size(evSoma,2)~=numel(tWinSTA), evSoma=evSoma(:)'; end   % nEv x W
    evSoma = evSoma - mean(evSoma(:, tWinSTA<-5), 2, 'omitnan');  % per-event pre-onset baseline
    somaEvN{i}=evSoma; nEvN(i)=size(evSoma,1);
    if any(distROI)                                              % distal dendrite STA (same triggers)
        [~, evDist] = get_STA(distTr, trig, tauPre, tauPost);
        evDist = squeeze(evDist); if size(evDist,2)~=numel(tWinSTA), evDist=evDist(:)'; end
        evDist = evDist - mean(evDist(:, tWinSTA<-5), 2, 'omitnan');
        distEvN{i}=evDist; distEvAll=[distEvAll; evDist];        %#ok<AGROW>
    end
    for o=1:maxOrd                                                % per-order spike-train alignment
        [~, evO] = get_STA(double(csSpk & (PCresult{f}.SpikeOrder==o)), trig, tauPre, tauPost);
        if ~isempty(evO), ordCntN{i}(o,:)=reshape(sum(evO,2),1,[]); end   % spike counts per ms offset
    end
    somaEvAll=[somaEvAll; evSoma]; ordCntAll=ordCntAll+ordCntN{i}; nEvAll=nEvAll+nEvN(i);   %#ok<AGROW>
end
% ---- pooled across neurons (fig 313) ----
figure(313); clf;
yyaxis left; hold on; set(gca,'YColor','k'); hSTA=gobjects(1,2);
if ~isempty(somaEvAll)
    hSTA(1)=errorbar_shade(tWinSTA, mean(somaEvAll,1,'omitnan'), std(somaEvAll,0,1,'omitnan')./sqrt(max(1,size(somaEvAll,1))), cSoma);
end
if ~isempty(distEvAll)
    hSTA(2)=errorbar_shade(tWinSTA, mean(distEvAll,1,'omitnan'), std(distEvAll,0,1,'omitnan')./sqrt(max(1,size(distEvAll,1))), cDist);
end
ylabel('STA (/spike ht)'); xline(0,'k:');
yyaxis right; hold on; set(gca,'YColor',[0.3 0.3 0.3]); hh=gobjects(1,maxOrd);
for o=1:maxOrd
    ts=repelem(tWinSTA, ordCntAll(o,:)); if isempty(ts), continue; end   % reconstruct spike-offset times
    hh(o)=histogram(ts, binE, 'DisplayStyle','stairs', 'EdgeColor',cmapOrd(o,:), 'LineWidth',1.3);
end
ylabel('# spikes'); xlim([tWinSTA(1) tWinSTA(end)]);
hLeg=[hSTA hh]; lblLeg=[{'soma STA','distal STA'}, arrayfun(@(o)sprintf('spike %d',o),1:maxOrd,'UniformOutput',false)];
vv=isgraphics(hLeg); legend(hLeg(vv), lblLeg(vv), 'Box','off','Location','northeast');
xlabel('Time from first CS spike (ms)'); box off;
title(sprintf('CS first-spike STA (soma & distal) + spike-order timing (pooled: %d CS, %d neurons, no blue)', nEvAll, numel(foi)));
set_fontsize(11);

% ---- per neuron (fig 314) ----
figure(314); clf; nN=numel(foi); nc=ceil(sqrt(nN)); nr=ceil(nN/nc);
tiledlayout(nr,nc,'TileSpacing','compact','Padding','compact');
for i=1:nN
    nexttile; yyaxis left; hold on; set(gca,'YColor','k'); hSTA=gobjects(1,2);
    if ~isempty(somaEvN{i})
        hSTA(1)=errorbar_shade(tWinSTA, mean(somaEvN{i},1,'omitnan'), std(somaEvN{i},0,1,'omitnan')./sqrt(max(1,size(somaEvN{i},1))), cSoma);
    end
    if ~isempty(distEvN{i})
        hSTA(2)=errorbar_shade(tWinSTA, mean(distEvN{i},1,'omitnan'), std(distEvN{i},0,1,'omitnan')./sqrt(max(1,size(distEvN{i},1))), cDist);
    end
    if mod(i-1,nc)==0, ylabel('STA'); end
    xline(0,'k:');
    yyaxis right; hold on; set(gca,'YColor',[0.3 0.3 0.3]); hh=gobjects(1,maxOrd);
    for o=1:maxOrd
        ts=repelem(tWinSTA, ordCntN{i}(o,:)); if isempty(ts), continue; end
        hh(o)=histogram(ts, binE, 'DisplayStyle','stairs', 'EdgeColor',cmapOrd(o,:), 'LineWidth',1);
    end
    xlim([tWinSTA(1) tWinSTA(end)]); box off;
    title(sprintf('f%d (n=%d CS)',foi(i),nEvN(i)),'FontSize',8);
    if i>nN-nc, xlabel('t from 1st CS spike (ms)'); end
    if i==1
        hLeg=[hSTA hh]; lblLeg=[{'soma','distal'}, arrayfun(@(o)sprintf('spk %d',o),1:maxOrd,'UniformOutput',false)];
        vv=isgraphics(hLeg); if any(vv), legend(hLeg(vv), lblLeg(vv), 'Box','off','FontSize',6,'Location','northeast'); end
    end
end
sgtitle('[3h] CS first-spike STA (soma black, distal red) + spike-order timing, per neuron (no blue)');
set_fontsize(8);

%% ===== Local helper functions (shared by the cells above) =====

function K = fetch_raw_kymo(m, AlignedEvntall, PCresult, onsetCol, exWin)
% On-demand RAW dendritic kymograph for one event record m (fields .f, .ev, .dist).
% Returns K.img (dist x time), K.y (dendritic distance), K.t (ms rel. onset), K.overlay (a
% 1-D trace over K.t to draw on top, e.g. the detector's slow trace). Kept light so the
% interactive browser can fetch per click instead of caching every event's image.
ec = onsetCol + exWin;
img = AlignedEvntall{m.f,1}(:, ec, m.ev);                 % RAW aligned window (nROI x nT)
img = img - mean(img(:, exWin>=-40 & exWin<0), 2, 'omitnan');   % per-ROI baseline
K = struct('img', img, 'y', PCresult{m.f}.dendaxis(:), 't', exWin, 'overlay', m.dist(:)');
end

function interactive_feature_kymo(figNum, EV, featColors, xName, yName, getKymo, tOv, spkColors, ttl)
% Generalizable interactive event browser.
%   LEFT (tall): scatter of scalar xName vs yName, one point per event in EV(i).scalar.
%   RIGHT top->bottom: RAW kymograph via getKymo(i) (soma top -> distal bottom); the somatic
%   slow trace EV(i).somaSlow; the dendritic slow trace EV(i).dist -- all sampled on tOv and
%   sharing the time axis. Every panel marks each field of EV(i).feat as a dotted line coloured
%   by featColors, and the spikes in EV(i).spk (struct of per-class time vectors) coloured by
%   spkColors (somatic classes on the soma trace, dSpike on the distal trace, all on the kymo).
% To visualize other features: add fields to EV.feat / EV.scalar and pass their names.
nE = numel(EV); if nE==0, return; end
X = arrayfun(@(e) e.scalar.(xName), EV);
Y = arrayfun(@(e) e.scalar.(yName), EV);
fnames  = fieldnames(EV(1).feat);              % time-marker features (drawn on every panel)
somaCls = {'SS','CS','BS'};                    % somatic classes -> soma trace; dSpike -> distal
fig = figure(figNum); clf(fig); set(fig,'Color','w');
ax1 = subplot(3,2,[1 3 5]); hold(ax1,'on'); box(ax1,'on');
scatter(ax1, X, Y, 18, 'filled', 'MarkerFaceAlpha',0.45, 'MarkerFaceColor',[0.35 0.35 0.85], ...
    'HitTest','off','PickableParts','none');
xline(0,'k:'); yline(0,'k:');                   % reference lines (0 = first spike / zero lag)
hSel = plot(ax1, NaN, NaN, 'ko', 'MarkerSize',11, 'LineWidth',1.8);            % selected-point ring
xlabel(ax1, sprintf('%s (ms)', xName)); ylabel(ax1, sprintf('%s (ms)', yName));
title(ax1, sprintf('click a point  (n=%d events)', nE));
ax1.ButtonDownFcn = @onClick;
ax2 = subplot(3,2,2); ax3 = subplot(3,2,4); ax4 = subplot(3,2,6);
sgtitle(ttl);
showEvent(1);                                             % show one by default
    function onClick(src, ~)
        cp = src.CurrentPoint(1,1:2);
        rx = max(max(X)-min(X),eps); ry = max(max(Y)-min(Y),eps);   % normalize so "nearest" is visual
        [~, sel] = min(((X-cp(1))/rx).^2 + ((Y-cp(2))/ry).^2);
        showEvent(sel);
    end
    function showEvent(idx)
        E = EV(idx); K = getKymo(idx);
        % --- kymograph (soma top -> distal bottom) ---
        cla(ax2); imagesc(ax2, K.t, K.y, K.img); set(ax2,'YDir','reverse'); hold(ax2,'on');
        vv = K.img(isfinite(K.img));
        if ~isempty(vv), c = max(abs(prctile(vv,[2 98]))); if c>0, clim(ax2,[-1 1]*c); end, end
        colormap(ax2, turbo);
        mark_panel(ax2, E, fnames, featColors, spkColors, fieldnames(spkColors), true);
        ylabel(ax2, 'dist. from soma (\mum)'); title(ax2, E.label);
        % --- somatic slow-component trace ---
        cla(ax3); hold(ax3,'on'); plot(ax3, tOv, E.somaSlow, 'k-', 'LineWidth',1.2);
        mark_panel(ax3, E, fnames, featColors, spkColors, somaCls, false);
        ylabel(ax3, 'soma slow (\DeltaF/F)');
        % --- dendritic slow-component trace ---
        cla(ax4); hold(ax4,'on'); plot(ax4, tOv, E.dist, '-', 'Color',[0.85 0.3 0.1], 'LineWidth',1.2);
        mark_panel(ax4, E, fnames, featColors, spkColors, {'dSpike'}, false);
        ylabel(ax4, 'distal slow (\DeltaF/F)'); xlabel(ax4, 't from 1st spike (ms)');
        linkaxes([ax2 ax3 ax4], 'x'); xlim(ax4, [tOv(1) tOv(end)]);
        set(hSel, 'XData', X(idx), 'YData', Y(idx));       % ring the selected point
    end
end

function mark_panel(ax, E, fnames, featColors, spkColors, spkTypes, topTicks)
% Overlay the first-spike line, the feature dotted lines, and the spike markers on one panel.
yl = ylim(ax); if diff(yl)<=0, yl=[-1 1]; end
plot(ax, [0 0], yl, 'k:');                                  % first spike (t=0)
for q = 1:numel(fnames)                                     % feature time markers
    fn = fnames{q}; tv = E.feat.(fn);
    if isfield(featColors,fn), col = featColors.(fn); else, col = [0.4 0.4 0.4]; end
    plot(ax, [tv tv], yl, ':', 'Color',col, 'LineWidth',1.8);
end
for q = 1:numel(spkTypes)                                   % spikes in the window
    tp = spkTypes{q}; if ~isfield(E.spk,tp) || isempty(E.spk.(tp)), continue; end
    col = spkColors.(tp); ts = E.spk.(tp);
    if topTicks                                             % short ticks along the top edge (kymo)
        for tt = ts(:)', plot(ax, [tt tt], yl(1)+[0 0.06]*(yl(2)-yl(1)), '-', 'Color',col, 'LineWidth',1.5); end
    else                                                    % full-height thin lines (trace panels)
        for tt = ts(:)', plot(ax, [tt tt], yl, '-', 'Color',col, 'LineWidth',0.8); end
    end
end
ylim(ax, yl);
end

function spk = spikes_in_window(scv, onF, tRel, nTime)
% Times (rel. to onset, ms) of each spike class within a display window.
% scv rows: 1 SS, 2 CS, 3 BS, 4 dSpike (PCresult.SpikeClassVecMat). onF = onset frame.
names = {'SS','CS','BS','dSpike'}; spk = struct();
fr = onF + tRel; ok = fr>=1 & fr<=nTime;
for c = 1:numel(names)
    hit = false(1, numel(tRel));
    if size(scv,1) >= c, hit(ok) = scv(c, fr(ok)) > 0; end
    spk.(names{c}) = tRel(hit);
end
end

function [mu, se] = phase_tuning(x, y, edges)
% Mean +- s.e.m. of y within bins of x (edges). Bins with <3 finite points -> NaN.
b = discretize(x, edges); nb = numel(edges)-1; mu = nan(1,nb); se = nan(1,nb);
for k = 1:nb
    yy = y(b==k & isfinite(y));
    if numel(yy) >= 3, mu(k) = mean(yy); se(k) = std(yy)/sqrt(numel(yy)); end
end
end

function lag = subms_peaklag(cc, lg, step)
% Sub-millisecond peak lag of a cross-correlation: spline-interpolate cc (given
% at integer-sample lags lg) onto a fine lag grid and return the argmax lag.
if nargin<3 || isempty(step), step = 0.05; end
cc = cc(:)'; lg = lg(:)';
good = isfinite(cc);
if nnz(good) < 4, [~,mi]=max(cc); lag = lg(mi); return; end
lagFine = lg(1):step:lg(end);
ccFine  = interp1(lg(good), cc(good), lagFine, 'spline');
[~,mi]  = max(ccFine);
lag = lagFine(mi);
end

function OrthAxis_vec = build_orthaxis_vec(OrthAxis, roisD_order_f, nROI)
% Map the dendrite-class weighting OrthAxis (rows = projection axes, columns =
% dClass 1..5) onto per-ROI weights, matching the inline construction used by
% the CS-vs-SS cells: each axis distributes its per-class weight equally over the
% ROIs of that class; a requested-but-absent class marks the whole axis NaN.
roisD_order_ind = cellfun(@find, roisD_order_f, 'UniformOutput', false);
OrthAxis_vec = zeros(size(OrthAxis,1), nROI);
for d = 1:5
    dclassinds = cell2mat(roisD_order_ind(d,:)');
    for v = 1:size(OrthAxis,1)
        if abs(OrthAxis(v,d))>0 && isempty(dclassinds)
            OrthAxis_vec(v,:) = NaN;
        else
            suborthvec = ind2vec(nROI, dclassinds, OrthAxis(v,d)/length(dclassinds));
            OrthAxis_vec(v,abs(suborthvec)>0) = suborthvec(abs(suborthvec)>0);
        end
    end
end
end

function showScaleScatter_rot(vals, str, ftprnt, cmapName, cmapRange, Rv)
% Like showScaleScatter, but ROI membership/coloring is computed in the ORIGINAL
% coordinates while the points are plotted at coordinates rotated by Rv (a 2x2
% screen-space rotation; build it as  T=viewmtx(az,90); Rv=T(1:2,1:2) ). This
% rotates the morphology WITHOUT view(), so a vertical scalebar stays vertical.
if nargin<4 || isempty(cmapName),  cmapName='parula'; end
if nargin<5 || isempty(cmapRange), cmapRange=[min(vals) max(vals)]; end
if nargin<6 || isempty(Rv),        Rv=eye(2); end
if size(str,2)<4, str(:,4)=1; end
str(str(:,4)==0,4)=1;

N  = size(ftprnt,3);
sz = [size(ftprnt,1), size(ftprnt,2)];
rot = @(xy) xy*Rv';                       % xy = [X Y] = [str(:,2) str(:,1)]

Pall = rot([str(:,2) str(:,1)]);
scatter(Pall(:,1), Pall(:,2), str(:,4), [0.5 0.5 0.5], 'filled');
axis equal tight off; hold on;

colors = vec2cmap(vals, cmapName, cmapRange);
valid  = str(:,1)>0.5 & str(:,1)<(sz(2)-0.5) & str(:,2)>0.5 & str(:,2)<(sz(1)-0.5);
for r = 1:N
    DMDmask = ftprnt(:,:,r) > 0;
    pts = str(valid,:);
    ind = sub2ind(sz, round(pts(:,2)), round(pts(:,1)));
    inside = DMDmask(ind);
    P = rot([pts(inside,2) pts(inside,1)]);
    scatter(P(:,1), P(:,2), pts(inside,4), colors(r,:), 'filled');
end
colormap(cmapName); caxis(cmapRange);
end

function p = pair_p(a, b)
% Pairwise (unpaired) p-value between two groups via get_pValue, with the
% function's console printout suppressed. NaNs are dropped inside get_pValue.
[~, pp] = evalc('get_pValue({a(:), b(:)}, 0)');
p = pp(1,2);
end

function s = sig_star(p)
% Significance string: '*' <0.05, '**' <0.01, '***' <0.001, else 'n.s.'
if p < 0.001, s = '***'; elseif p < 0.01, s = '**'; elseif p < 0.05, s = '*'; else, s = 'n.s.'; end
end

function SI = skaggs_info(rate, occ)
% Skaggs spatial information (bits/spike). rate: 1xB firing rate (Hz) per bin,
% occ: 1xB dwell time (frames) per bin. Non-finite bins are treated as empty.
occ(~isfinite(occ))=0; rate(~isfinite(rate))=0;
if sum(occ)<=0, SI=0; return; end
p    = occ/sum(occ);            % occupancy probability per bin
rbar = sum(p.*rate);           % occupancy-weighted mean rate
if rbar<=0, SI=0; return; end
idx  = rate>0 & occ>0;
SI   = sum(p(idx).*(rate(idx)/rbar).*log2(rate(idx)/rbar));
end

function [SpikeTable, EventTable, NeuronTable] = build_spike_event_tables(PCresult, bAPPropsMat, EventPropsMat, foi)
%BUILD_SPIKE_EVENT_TABLES Consolidate the per-spike / per-event descriptors.
%
% Replaces the scattered bAPPropsMat{f,1..2} / EventPropsMat{f,1..2} variables
% with three tidy tables that span every neuron and are sorted by (Neuron,Frame):
%
%   SpikeTable  one row per somatic spike (finest granularity). Carries the bAP
%               propagation features, the behavioural state at the spike, and the
%               attributes of the event the spike belongs to (denormalized join).
%   EventTable  one row per spike event (SS / CS / BS), linked back to its spikes
%               through SpikeUIDs / FirstSpikeUID.
%   NeuronTable one row per neuron: dendritic axis, ROI dendrite classes, geometry.
%
% The full ROI x time aligned traces are deliberately NOT copied into the tables
% (they would duplicate several GB). Slice them with (Neuron, LocalIdx), e.g.
%   r = SpikeTable(k,:);  tr = AlignedbAPall{r.Neuron,2}(:,:,r.LocalIdx);
% ROI-resolved AUC vectors and the distance-binned traces (small) ARE included.

spikeClassNames = ["SS" "CS" "BS" "dSP"];   % PCresult.SpikeClassVecMat row order
eventClassNames = ["SS" "CS" "BS"];         % EventPropsMat.Spike_type

Scell = cell(numel(foi),1); Ecell = cell(numel(foi),1); Ncell = cell(numel(foi),1);
Dcell = cell(numel(foi),1);   % dSpike (dendritic-spike) event rows, appended to EventTable
spikeOffset = 0; eventOffset = 0;

for k = 1:numel(foi)
    f  = foi(k);
    Bt = bAPPropsMat{f,1};  B2 = bAPPropsMat{f,2};
    Et = EventPropsMat{f,1}; E2 = EventPropsMat{f,2};
    ns = height(Bt); nev = height(Et);

    dax = PCresult{f}.dendaxis(:);       % nROI x 1, signed distance from soma (um)
    VR  = PCresult{f}.VR;
    nROI = numel(dax);
    % roi_dClass is filled by logical assignment, so it can be short if the last
    % ROIs belong to no dendrite class; pad with 0 (= unclassified).
    dc = PCresult{f}.roi_dClass(:);      % 1 basal 2 soma 3 trunk 4 oblique 5 distal
    if numel(dc) < nROI, dc(numel(dc)+1:nROI,1) = 0; end

    AUCraw   = B2.AUCbAPraw{1};          % nROI x ns, bAP AUC per ROI
    AUCshort = B2.AUCbAPrawShort{1};     % nROI x ns, short-window (peak) AUC

    clsMean = @(M,c) mean(M(dc==c,:),1,'omitnan')';   % mean over ROIs of one dendrite class

    % ---- spatial description of each bAP ----
    W    = AUCraw; W(W<0) = 0;                        % rectify -> well-defined centroid
    Wsum = sum(W,1,'omitnan');
    AUC_centroid_pos = (sum(dax.*W,1,'omitnan')./Wsum)';
    AUC_spread = sqrt(sum(W.*(dax - AUC_centroid_pos').^2,1,'omitnan')./Wsum)';
    [AUC_peak, imax] = max(AUCraw,[],1);
    AUC_peak = AUC_peak(:);  AUC_peakDist = dax(imax(:));
    nROI_active = sum(AUCraw>0,1)';

    AUCshort_soma = clsMean(AUCshort,2);

    % ---- behavioural state at each spike ----
    frames  = Bt.SpikeFrame;
    speedTr = movmean(VR(end,:),200);
    VRpos = VR(5,frames)';  Lap = VR(8,frames)';  Speed = speedTr(frames)';

    % ---- which event does each spike belong to? ----
    evStart = Et.SpikeFrame;  evEnd = evStart + Et.Length;
    EventLocalIdx = NaN(ns,1);  PosInEvent = NaN(ns,1);
    for ev = 1:nev
        m = find(frames>=evStart(ev) & frames<=evEnd(ev));
        if isempty(m), continue; end
        EventLocalIdx(m) = ev;
        [~,ord] = sort(frames(m));
        PosInEvent(m(ord)) = (1:numel(m))';
    end
    inEv = ~isnan(EventLocalIdx);

    % ---- spike table for this neuron ----
    S = table();
    S.SpikeUID = (1:ns)' + spikeOffset;
    S.Neuron   = repmat(f,ns,1);
    S.LocalIdx = (1:ns)';                 % slices AlignedbAPall{f,:}(:,:,LocalIdx)
    S.Frame    = frames;
    S.Time_s   = frames/1000;

    S.SpikeClass = categorical(spikeClassNames(Bt.Spike_Type)', spikeClassNames);
    S.SpikeOrder = Bt.SpikeOrder;
    S.IsFirstInEvent = PosInEvent == 1;
    S.IsCS       = logical(Bt.IsCS);
    S.IsIsolated = logical(Bt.IsIsolated);

    S.IsBlue  = logical(Bt.IsBlue);
    S.IsNA    = logical(Bt.IsNA);
    S.IsPF    = logical(Bt.IsPF);
    S.IsRun   = logical(Bt.IsRun);
    S.IsValid = ~logical(Bt.IsBlue) & ~logical(Bt.IsNA);   % usable for analysis

    S.ISI    = Bt.ISI;                                   % ms since previous spike
    S.logISI = log10(Bt.ISI);
    S.ISI_next = [Bt.ISI(2:end); NaN];

    % bAP propagation (existing)
    S.AUC_soma        = Bt.AUC_soma;
    S.AUC_apical      = Bt.AUC_apical;
    S.AUC_apical_norm = Bt.AUC_apical_norm;
    S.AUC_centroid    = Bt.AUC_centroid;
    S.TransmissionRate       = Bt.TransmissionRate;
    S.TransmissionRate_equal = Bt.TransmissionRate_equal;
    S.TransmissionRate_norm  = Bt.TransmissionRate_norm;

    % bAP propagation (new, per dendrite class)
    S.AUC_basal   = clsMean(AUCraw,1);
    S.AUC_somaROI = clsMean(AUCraw,2);
    S.AUC_trunk   = clsMean(AUCraw,3);
    S.AUC_oblique = clsMean(AUCraw,4);
    S.AUC_distal  = clsMean(AUCraw,5);
    S.AUC_somaShort = AUCshort_soma;
    S.Attenuation   = clsMean(AUCraw,5) ./ AUCshort_soma;  % distal / somatic peak
    S.AUC_peak      = AUC_peak;
    S.AUC_peakDist  = AUC_peakDist;       % dendritic distance of the largest AUC (um)
    S.AUC_centroidDist = AUC_centroid_pos;% AUC-weighted mean distance (um)
    S.AUC_spread    = AUC_spread;         % AUC-weighted s.d. of distance (um)
    S.nROI_active   = nROI_active;

    % behaviour
    S.VRpos = VRpos; S.Lap = Lap; S.Speed = Speed;

    % event linkage (denormalized)
    S.EventUID = NaN(ns,1);  S.EventUID(inEv) = EventLocalIdx(inEv) + eventOffset;
    S.EventLocalIdx = EventLocalIdx;
    S.PosInEvent    = PosInEvent;
    evClassStr = strings(ns,1);  evClassStr(:) = missing;
    evClassStr(inEv) = eventClassNames(Et.Spike_type(EventLocalIdx(inEv)));
    S.EventClass  = categorical(evClassStr, eventClassNames);
    S.EventLength = NaN(ns,1); S.EventLength(inEv) = Et.Length(EventLocalIdx(inEv));
    S.EventNspike = NaN(ns,1); S.EventNspike(inEv) = Et.Nspike(EventLocalIdx(inEv));

    % ROI-resolved / binned quantities (small enough to carry)
    S.AUC_perROI      = num2cell(AUCraw,1)';
    S.AUCshort_perROI = num2cell(AUCshort,1)';
    S.AUC_perROI_ch1  = num2cell(B2.AUCbAPraw_ch1{1},1)';
    S.AUC_perROI_ch2  = num2cell(B2.AUCbAPraw_ch2{1},1)';
    S.AUCshort_perROI_ch1 = num2cell(B2.AUCbAPrawShort_ch1{1},1)';
    S.AUCshort_perROI_ch2 = num2cell(B2.AUCbAPrawShort_ch2{1},1)';
    S.AUC_binned      = B2.AUC_binned{1}(:);
    S.AUC_binned_fine = B2.AUC_binned_fine{1}(:);
    S.Amp_binned      = B2.Amp_binned{1}(:);
    S.AlignedbAP_binned = squeeze(num2cell(B2.AlignedbAP_binned{1},[1 2]));

    % ---- event table for this neuron ----
    AUCev = E2.AUC{1};                      % nROI x nev
    Wev = AUCev; Wev(Wev<0) = 0; WevSum = sum(Wev,1,'omitnan');
    EV_centroidDist = (sum(dax.*Wev,1,'omitnan')./WevSum)';

    E = table();
    E.EventUID = (1:nev)' + eventOffset;
    E.Neuron   = repmat(f,nev,1);
    E.LocalIdx = (1:nev)';                  % slices AlignedEvntall{f,:}(:,:,LocalIdx)
    E.Frame    = Et.SpikeFrame;
    E.Time_s   = Et.SpikeFrame/1000;

    E.EventClass = categorical(eventClassNames(Et.Spike_type)', eventClassNames);
    E.Length = Et.Length;                   % event duration (ms)
    E.Nspike = Et.Nspike;

    E.IsBlue  = logical(Et.IsBlue);
    E.IsNA    = logical(Et.IsNA);
    E.IsPF    = logical(Et.IsPF);
    E.IsRun   = logical(Et.IsRun);
    E.IsValid = ~logical(Et.IsBlue) & ~logical(Et.IsNA);

    E.ISI = Et.ISI;  E.logISI = log10(Et.ISI);

    E.AUC_soma = Et.AUC_soma;  E.AUC_apical = Et.AUC_apical;
    E.AUC_centroid = Et.AUC_centroid;  E.Amp_apical = Et.Amp_apical;
    E.TransmissionRate = Et.TransmissionRate;
    E.TRfrst_bAP = Et.TRfrst_bAP;  E.TRmax_bAP = Et.TRmax_bAP;
    E.TRfrst_bAP_norm = Et.TRfrst_bAP_norm;  E.TRmax_bAP_norm = Et.TRmax_bAP_norm;
    E.AUCapical_frst = Et.AUCapical_frst;  E.AUCapical_frst_norm = Et.AUCapical_frst_norm;

    E.PreSub_Soma = Et.PreSub_Soma;  E.PreSub_Basal = Et.PreSub_Basal;
    E.PreSub_Distal = Et.PreSub_Distal;  E.PreSub_Apical = Et.PreSub_Apical;
    E.PostSub_Basal = Et.PostSub_Basal;  E.PostSub_Distal = Et.PostSub_Distal;

    % new features: the basal/distal "seesaw" and global drive used by the figures
    E.SeesawAmp = Et.PreSub_Basal - Et.PreSub_Distal;
    E.GlobalAmp = Et.PreSub_Basal + Et.PreSub_Distal;
    E.AUC_basal  = clsMean(AUCev,1);
    E.AUC_distal = clsMean(AUCev,5);
    E.AUC_centroidDist = EV_centroidDist;

    E.AUC_perROI = num2cell(AUCev,1)';
    E.AUC_binned = E2.AUC_binned{1}(:);
    E.AlignedEv_binned = squeeze(num2cell(E2.AlignedEv_binned{1},[1 2]));
    E.bAP_TransmissionRates = E2.bAPtransmissionRate{1}(:);

    % link events -> spikes
    E.SpikeUIDs = arrayfun(@(ev) S.SpikeUID(EventLocalIdx==ev), (1:nev)', 'UniformOutput', false);
    E.FirstSpikeUID = NaN(nev,1);
    for ev = 1:nev
        if ~isempty(E.SpikeUIDs{ev}), E.FirstSpikeUID(ev) = E.SpikeUIDs{ev}(1); end
    end

    % ---- neuron table ----
    N = table();
    N.Neuron = f;  N.nROI = nROI;  N.nSpike = ns;  N.nEvent = nev;
    N.Pixelsize = PCresult{f}.pixelsize;
    N.SpikeHeight = PCresult{f}.SpikeHeight;
    N.RecDuration_s = size(VR,2)/1000;
    N.dendaxis = {dax'};  N.roi_dClass = {dc'};

    % ---- dSpike (dendritic-spike, class 4) events for this neuron ----
    % dSpikes are NOT in AlignedEvntall, so LocalIdx / bAP / AUC fields are NaN; only
    % Frame, class, validity, ISI and behavioural placeholders are filled. Built from a
    % 1-row template of E so the schema matches exactly for vertcat.
    dsp = PCresult{f}.SpikeClassVecMat(4,:) > 0;
    dspFrames = find(diff([0 dsp])==1)';        % rising edges = dSpike onsets (column)
    if ~isempty(dspFrames) && height(E)>=1      % need >=1 somatic event as a schema template
        nds = numel(dspFrames);
        Dtab = repmat(E(1,:), nds, 1);          % template row -> exact schema
        vn = E.Properties.VariableNames;
        for kk = 1:numel(vn)                    % clear all inherited values to NaN/false/[]
            col = E.(vn{kk});
            if     isnumeric(col), Dtab.(vn{kk}) = nan(nds, size(col,2));
            elseif islogical(col), Dtab.(vn{kk}) = false(nds,1);
            elseif iscell(col),    Dtab.(vn{kk}) = repmat({[]}, nds,1);
            end
        end
        blue = PCresult{f}.BlueStim;
        Dtab.EventUID   = nan(nds,1);           % renumbered after assembly
        Dtab.Neuron     = repmat(f,nds,1);
        Dtab.LocalIdx   = nan(nds,1);           % no AlignedEvntall slice for dSpikes
        Dtab.Frame      = dspFrames;
        Dtab.Time_s     = dspFrames/1000;
        Dtab.EventClass = repmat(categorical("dSpike"), nds, 1);
        Dtab.Length     = ones(nds,1);
        Dtab.Nspike     = ones(nds,1);
        Dtab.IsBlue     = blue(dspFrames)' > 0;    % blue is a row vector -> transpose to nds x 1
        Dtab.IsNA       = false(nds,1);
        Dtab.IsPF       = false(nds,1);
        Dtab.IsRun      = false(nds,1);
        Dtab.IsValid    = ~Dtab.IsBlue;
        Dtab.ISI        = [NaN; diff(dspFrames)];
        Dtab.logISI     = log10(Dtab.ISI);
        Dcell{k} = Dtab;
    end

    Scell{k} = S; Ecell{k} = E; Ncell{k} = N;
    spikeOffset = spikeOffset + ns;  eventOffset = eventOffset + nev;
end

SpikeTable  = sortrows(vertcat(Scell{:}), {'Neuron','Frame'});
Dall = Dcell(~cellfun(@isempty, Dcell));
if isempty(Dall)
    EventTable = sortrows(vertcat(Ecell{:}), {'Neuron','Frame'});
else
    EventTable = sortrows([vertcat(Ecell{:}); vertcat(Dall{:})], {'Neuron','Frame'});
end
NeuronTable = vertcat(Ncell{:});
% assign fresh unique EventUIDs to the appended dSpike rows (above the somatic ones)
isD = EventTable.EventClass=="dSpike";
if any(isD)
    base = max(EventTable.EventUID(~isD)); if isempty(base)||isnan(base), base=0; end
    EventTable.EventUID(isD) = base + (1:sum(isD))';
end

SpikeTable.Properties.Description = 'One row per somatic spike (bAP), all neurons, sorted by neuron then frame.';
EventTable.Properties.Description = 'One row per spike event (SS/CS/BS/dSpike), all neurons, sorted by neuron then frame.';
NeuronTable.Properties.Description = 'One row per neuron: dendritic axis, ROI dendrite classes, geometry.';
SpikeTable.Properties.VariableUnits = table_units(SpikeTable);
EventTable.Properties.VariableUnits = table_units(EventTable);
end

function u = table_units(T)
% Units for the columns we can name; '' for everything else.
u = repmat({''},1,width(T));
ms  = {'ISI','ISI_next','Length','EventLength'};
um  = {'AUC_peakDist','AUC_centroidDist','AUC_spread'};
sec = {'Time_s','RecDuration_s'};
for n = ms,  u(strcmp(T.Properties.VariableNames,n{1})) = {'ms'}; end
for n = um,  u(strcmp(T.Properties.VariableNames,n{1})) = {'um'}; end
for n = sec, u(strcmp(T.Properties.VariableNames,n{1})) = {'s'};  end
u(strcmp(T.Properties.VariableNames,'Frame')) = {'frame'};
u(strcmp(T.Properties.VariableNames,'VRpos')) = {'cm'};
end

function [bAPmat, NonBlueInd] = build_branchNorm_bAPmat(SpikeTable, PCresult, f)
% Per-spike bAP AUC matrix (ROI x spike), soma rows replaced by the short-window
% AUC, normalized to the mean somatic short AUC of the same spike, then centered
% in time. Only non-blue-stim spikes are kept.
St = SpikeTable(SpikeTable.Neuron==f,:);
somaROI    = ismember(PCresult{f}.roi_dClass,2);
NonBlueInd = ~St.IsBlue;
AUCraw   = cell2mat(St.AUC_perROI');        % nROI x nSpike
AUCshort = cell2mat(St.AUCshort_perROI');   % nROI x nSpike
bAPmat = AUCraw(:,NonBlueInd);
bAPmat(somaROI,:) = AUCshort(somaROI,NonBlueInd);
bAPmat = bAPmat ./ mean(AUCshort(somaROI,NonBlueInd),1,'omitnan');
bAPmat = bAPmat - mean(bAPmat,2,'omitnan');
end

function [corrTXnPreSub, corrTXnPreSub_timelapse, t2avg] = get_presub_apicalAUC_corr(EventTable, AlignedEvntall, PCresult, foi, nTau_EV)
% For each session, correlate the (normalized) first-bAP apical AUC of every
% event with the pre-spike subthreshold in distal (col 1/3) and basal (col 2/4)
% dendrite. corrTXnPreSub holds [rho_distal rho_basal p_distal p_basal];
% corrTXnPreSub_timelapse resolves the same correlations over pre-spike time
% windows t2avg.
t2avg = [-155:20:-5]; Presubmat = [];
corrTXnPreSub = NaN(max(foi),4);
corrTXnPreSub_timelapse = NaN(max(foi),4,length(t2avg)-1);
for f = foi
    Et = EventTable(EventTable.Neuron==f,:);
    Et = Et(~isnan(Et.LocalIdx),:);            % keep only AlignedEvntall-backed events (drop dSpike)
    apAUC = Et.AUCapical_frst_norm;
    validInd = ~isnan(apAUC) & ~isnan(Et.PreSub_Distal);
    [corrTXnPreSub(f,1), corrTXnPreSub(f,3)] = corr(apAUC(validInd), Et.PreSub_Distal(validInd), 'Type', 'Pearson');
    [corrTXnPreSub(f,2), corrTXnPreSub(f,4)] = corr(apAUC(validInd), Et.PreSub_Basal(validInd), 'Type', 'Pearson');
    for t = 1:length(t2avg)-1
        Presubtmp = squeeze(mean(AlignedEvntall{f}(:,-nTau_EV(1)+[t2avg(t):t2avg(t+1)],:),2,'omitnan'));
        Presubmat{f,1}(t,:) = mean(Presubtmp(PCresult{f}.roi_dClass==5,:),1,'omitnan');
        Presubmat{f,2}(t,:) = mean(Presubtmp(PCresult{f}.roi_dClass==1,:),1,'omitnan');
        validIndt = ~isnan(apAUC) & ~isnan(Presubmat{f,1}(t,:)');
        if any(validIndt)
            [corrTXnPreSub_timelapse(f,1,t) corrTXnPreSub_timelapse(f,3,t)] = corr(apAUC(validIndt), Presubmat{f,1}(t,validIndt)', 'Type', 'Pearson');
        end
        validIndt2 = ~isnan(apAUC) & ~isnan(Presubmat{f,2}(t,:)');
        if any(validIndt2)
            [corrTXnPreSub_timelapse(f,2,t) corrTXnPreSub_timelapse(f,4,t)] = corr(apAUC(validIndt2), Presubmat{f,2}(t,validIndt2)', 'Type', 'Pearson');
        end
    end
end
end

function [dAPCandidates_threshold, dt2show, dst2show, bAPthreshold, catCandidates] = classify_dAP_candidates(props)
% Separate back-propagating APs (bAP) from dendritically-initiated spikes using
% a conduction-speed threshold on (delay-to-soma, distance-from-soma): a spike
% whose distance exceeds conductionspeed*(dt+dt_error) is too far to be a bAP.
conductionspeed = 230; % um/ms
dt_error        = 1;   % ms
bAPthreshold    = @(x) conductionspeed*(x+dt_error);
dAPCandidates  = find(~props.IsbAP & ~props.IsCS & props.dClass~=2 & props.dt_post_somAP<2.5);
bAPsCandidates = find( props.IsbAP & ~props.IsCS & props.dClass~=2 & props.dt_pre_somAP <4);
catCandidates  = [bAPsCandidates; dAPCandidates];
dt2show  = [props.dt_pre_somAP_int(bAPsCandidates); -props.dt_post_somAP_int(dAPCandidates)];
dst2show = [props.Distance_from_soma([bAPsCandidates; dAPCandidates])];
dAPCandidates_threshold = catCandidates(find((dst2show-bAPthreshold(dt2show))>0));
end
