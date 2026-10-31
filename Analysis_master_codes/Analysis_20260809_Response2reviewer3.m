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
mouse=cell2mat(raw(:,2))';
StructureData=raw(:,8);
BadROI=cellfun(@(x) (str2num(num2str(x))),raw(:,17),'UniformOutput',false);
EndFrame=cell2mat(raw(:,15));
ifmotionReject=cell2mat(raw(:,16));
ifdirtRemov=cell2mat(raw(:,18));
Pixelsize=cell2mat(raw(:,6));
NSeesawComponent=cell2mat(raw(:,25));
save_figto='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism';
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
nTau = [-200:50]; % SS, CS, Brst
clear PCresult

for f = 20%foi
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
    PCresult{f}.Ftprnts = Result.ftprnt;

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

    ftprintscoord=get_coord(PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist));
    SomaApicalcoord=ftprintscoord(ismember(PCresult{f}.branch_dClass,[2 3 5]),:);
    Somacoord=get_coord(PCresult{f}.Ftprnts(:,:,1));
    [~, principal_axis]= projectTrunkaxis(SomaApicalcoord);
    PCresult{f}.dendaxis1d = (ftprintscoord - Somacoord)* principal_axis;
    PCresult{f}.principal_axis=principal_axis;

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
    PCresult{f}.Subthreshold = get_subthreshold(NormalizedTrace_dirt, max(Result.spike(1,:), [], 1) > 0, 7, 20);
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

    roisD_order_ind=cellfun(@find,PCresult{f}.roisD_order,'UniformOutput',false);
    PCresult{f}.labelvec=NaN(1,nROI);
    for dClass=1:5
        PCresult{f}.labelvec(cell2mat(roisD_order_ind(dClass,:)'))=dClass;
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

        PCresult{f}.Subthreshold_ch{ch}= get_subthreshold(PCresult{f}.NormalizedTrace_ch{ch}, max(Result.spike(1,:), [], 1) > 0, 7, 17);
        PCresult{f}.Subthreshold_ch{ch}(isnan(NormalizedTrace_dirt)) = NaN;
    end

    PCresult{f}.BlueStim = Result.Blue;
    PCresult{f}.CStrace = CS_trace;
    PCresult{f}.BStrace = BS_trace;
    PCresult{f}.avgImg = Result.ref_im;
    PCresult{f}.VR = Result.VR;

    % % Detect peaks
    % StopFreq=[12 200]; %filter high frequency
    % [nROI nTime]=size(PCresult{f}.Subthreshold);
    % Subthreshold_int=interpolateNaN(PCresult{f}.Subthreshold);
    % FilterTr=[];
    % PCresult{f}.peakvec=zeros(nROI,nTime); PCresult{f}.troughvec=zeros(nROI,nTime);
    % Npeak=[]; Ntrough=[];
    % % figure(f+400); clf;
    % for b=unique(Result.BranchLabel)
    %     branchLabel=Result.BranchLabel(PCresult{f}.sortDist);
    %     branchInd=find(branchLabel==b);
    %     %[PhaseTrace BasalSubFilt BasalthetaPower] = get_phase(Subthreshold_int(n,:), 1000, FilterFreq);
    %     perispikefrm=unique(find(PCresult{f}.SpikeMat(1,:))'+[-5:50]);
    %     perispikefrm(perispikefrm<1 | perispikefrm>nTime)=[];
    %     FilterTr(b,:)=get_bandstop(mean(Subthreshold_int(branchInd,:),1,'omitnan')./PCresult{f}.SpikeHeight,1000,StopFreq);
    %     FilterTr(b,:)=FilterTr(b,:)-movmedian(FilterTr(b,:),300,2);
    %     [pks, locs] = findpeaks(FilterTr(b,:),'MinPeakHeight', 0.5,'MinPeakDistance',100, ...
    %         'MinPeakProminence', 0.2);  % Prominence can also be tuned
    %     [trough, loc_trgh] = findpeaks(-FilterTr(b,:),'MinPeakHeight', 0.5,'MinPeakDistance',100, ...
    %         'MinPeakProminence', 0.2);  % Prominence can also be tuned
    %     % Turn to 0 during blue Stim and peri-spike frame
    %     PCresult{f}.peakvec(branchInd,:)=repmat(ind2vec(nTime,locs,1),length(branchInd),1);
    %     PCresult{f}.peakvec(branchInd,PCresult{f}.BlueStim>0 | ind2vec(nTime,perispikefrm,1))=0;
    %     PCresult{f}.troughvec(branchInd,:)=repmat(ind2vec(nTime,loc_trgh,1),length(branchInd),1);
    %     PCresult{f}.troughvec(branchInd,PCresult{f}.BlueStim>0 | ind2vec(nTime,perispikefrm,1))=0;
    % end

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

%% ===== Correlation between ROIs — calculation (spatial & temporal) =====
% For each neuron and each frame condition (condName: Good, Silent, Run, PF) a
% single spatiotemporal correlation tensor Cmat = N x N x (2*maxLag+1) is built
% with get_corrTensor (raw traces, nanXCorrFFT on a continuous time axis).
%   * spatial correlation matrix = zero-lag slice  Cmat(:,:,maxLag+1)
%   * temporal correlation       = the lag profile  Cmat(i,j,:)
% Diagonal (same ROI) crosses the two optical channels (NormalizedTrace_ch) so
% independent shot noise does not inflate the zero-lag value; off-diagonal uses
% the combined trace (NormalizedTrace_dirt). Spatial decay (corr vs distance)
% and temporal decay (corr vs lag) are both extracted from this tensor. Auto
% class pairs (dClass2temporal{dc,1}==dClass2temporal{dc,2}) use diagonal pairs
% only; dClass2temporal row 5 ({1..5} auto) is the all-ROI pooled autocorr.
corrDataFile = fullfile(save_figto, 'Correlation_Fig2_data.mat');   % tensors + fit inputs
if ~isfile(corrDataFile)
perispike_time  = -3:10;                  % peri-spike frames to exclude
dilate_sub      = 7;  avgwnd_sub = 20;    % get_subthreshold parameters
maxLag          = 500;                    % temporal correlogram max lag (frames)
lagAxis         = (0:maxLag)';            % expfit sampling points
minPairs_xcorr  = 100;                    % min valid frame-pairs per lag for nanXCorr

dClass2temporal = [{1} {1}; {[3 5]} {[3 5]}; {[1]} {[5]}; {[5]} {[5]}; {[1 2 3 4 5]} {[1 2 3 4 5]}]; % class pairs
classNames      = {'Basal','Soma','Trunk','Oblique','Apical (distal)'};
title_str       = {'Basal, Basal (auto)', 'Apical, Apical (auto)', ...
                   'Basal, Distal (cross)', 'Distal, Distal (auto)',...
                   'All, All (auto)'};

nSess = max(foi);  nPair = size(dClass2temporal,1);
condName = {'Good','Silent','Run','PF'};   % frame conditions; one tensor per condition
nCond    = numel(condName);
heatCond = 2;                              % condition kept as full tensor (Silent) for Fig 127

corrMat         = cell(nSess, nCond);        % zero-lag spatial matrix per condition
temporalCorr    = cell(nSess, nCond, nPair); % class-pair correlograms
temporalCorrMat = cell(nSess, 1);            % full N x N x T tensor (heatCond only, single)
lagsCorr        = cell(nSess, 1);

% Accumulators for the fit-INPUT table (one row per would-be fit; no fitting here)
nD = 0;
CD_f = []; CD_type = {}; CD_class = {}; CD_cond = {}; CD_dA = {}; CD_dB = {};
CD_x = {}; CD_y = {};

for f = foi
    fprintf('Calculating #neuron %d ...\n', f);
    [nROI, nTime] = size(PCresult{f}.NormalizedTrace_dirt);
    nLag = 2*maxLag + 1;

    % raw traces: per channel (diagonal/auto) + combined (off-diagonal/cross)
    tr1 = PCresult{f}.NormalizedTrace_ch{1};
    tr2 = PCresult{f}.NormalizedTrace_ch{2};
    trC = PCresult{f}.NormalizedTrace_dirt;

    % --- frame selections ---
    spk = PCresult{f}.SpikeMat(1,:) > 0;
    perispike_frame = unique([tovec(find(spk)' + perispike_time); find(PCresult{f}.CStrace)']);
    perispike_frame(perispike_frame <= 0 | perispike_frame > nTime) = [];
    silenttime_vec  = ind2vec(nTime, perispike_frame, 0, 1);
    nonvalid_frame  = find(sum(isnan(PCresult{f}.Subthreshold),1) > 0);
    Blue_on_frame   = find(imdilate(PCresult{f}.BlueStim > 0, [ones(1,1), 1, ones(1,200)]));
    Goodframe       = setdiff(1:nTime, unique([Blue_on_frame nonvalid_frame]));
    goodvec         = ind2vec(nTime, Goodframe, 1) > 0;
    runVec = false(1,nTime);  if isfield(PCresult{f},'runVec'), runVec = PCresult{f}.runVec > 0; end
    PFvec  = false(1,nTime);  if isfield(PCresult{f},'PFvec'),  PFvec  = PCresult{f}.PFvec  > 0; end
    condMask = { goodvec, goodvec & silenttime_vec, goodvec & runVec, goodvec & PFvec };

    % pairwise dendritic distance (symmetric), upper triangle incl. diagonal
    Dpair   = PCresult{f}.interDendDist(PCresult{f}.sortDist, PCresult{f}.sortDist);
    Dpair   = max(cat(3,Dpair,Dpair'),[],3);
    triuIdx = triu(true(nROI));
    xd_all  = Dpair(triuIdx);
    allClasses = unique(PCresult{f}.labelvec(~isnan(PCresult{f}.labelvec)));

    for c = 1:nCond
        if nnz(condMask{c}) <= 2*maxLag, continue; end   % not enough frames

        % --- one spatiotemporal correlation tensor for this condition (raw) ---
        [Cmat, lags] = get_corrTensor(tr1, tr2, trC, condMask{c}, maxLag, minPairs_xcorr);
        lagsCorr{f}  = lags(:)';
        corrMat{f,c} = Cmat(:, :, maxLag+1);            % zero-lag = spatial matrix
        if c == heatCond, temporalCorrMat{f} = single(Cmat); end
        Cr = reshape(Cmat, nROI*nROI, nLag);            % rows indexed by (i,j) linear index

        % --- spatial fit input: zero-lag correlation vs distance, pooled over all pairs ---
        yc = corrMat{f,c}(triuIdx);
        ok = ~isnan(xd_all) & ~isnan(yc);
        if nnz(ok) > 5
            nD = nD + 1;
            CD_f(nD) = f;               CD_type{nD} = 'spatial';   CD_class{nD} = 'all pairs';
            CD_cond{nD} = condName{c};  CD_dA{nD} = allClasses;    CD_dB{nD} = allClasses;
            CD_x{nD} = xd_all(ok);      CD_y{nD} = yc(ok);
        end

        % --- per class pair: spatial (zero-lag) + temporal correlogram fit inputs ---
        % auto pair (classA==classB): diagonal only (true autocorrelation);
        % cross pair: A x B off-diagonal block.
        for dc = 1:nPair
            idxA = find(ismember(PCresult{f}.labelvec, dClass2temporal{dc,1}));
            idxB = find(ismember(PCresult{f}.labelvec, dClass2temporal{dc,2}));
            if isempty(idxA) || isempty(idxB), continue; end

            if isequal(sort(dClass2temporal{dc,1}), sort(dClass2temporal{dc,2}))
                linT = sub2ind([nROI nROI], idxA(:), idxA(:));       % temporal: diagonal
                mS   = triu(true(numel(idxA)));                      % spatial: unique+self
                Ds = Dpair(idxA, idxA);  Cs = corrMat{f,c}(idxA, idxA);
                xdc = Ds(mS);  ycc = Cs(mS);
            else
                [II, JJ] = ndgrid(idxA, idxB);
                linT = sub2ind([nROI nROI], II(:), JJ(:));           % temporal: A x B block
                Ds = Dpair(idxA, idxB);  Cs = corrMat{f,c}(idxA, idxB);
                xdc = Ds(:);  ycc = Cs(:);                           % spatial: A x B block
            end

            % spatial fit input for this class pair
            okc = ~isnan(xdc) & ~isnan(ycc);
            if nnz(okc) > 5
                nD = nD + 1;
                CD_f(nD) = f;               CD_type{nD} = 'spatial';   CD_class{nD} = title_str{dc};
                CD_cond{nD} = condName{c};  CD_dA{nD} = dClass2temporal{dc,1};  CD_dB{nD} = dClass2temporal{dc,2};
                CD_x{nD} = xdc(okc);        CD_y{nD} = ycc(okc);
            end

            % temporal fit input for this class pair (from the tensor)
            meanCorr = mean(Cr(linT, :), 1, 'omitnan');
            mc  = meanCorr(:);  okL = ~isnan(mc);
            if nnz(okL) > 5
                temporalCorr{f,c,dc} = meanCorr;         % full-length (NaN preserved)
                nD = nD + 1;
                CD_f(nD) = f;               CD_type{nD} = 'temporal';  CD_class{nD} = title_str{dc};
                CD_cond{nD} = condName{c};  CD_dA{nD} = dClass2temporal{dc,1};  CD_dB{nD} = dClass2temporal{dc,2};
                CD_x{nD} = abs(lags(okL));  CD_y{nD} = mc(okL);
            end
        end
    end
    fprintf('  neuron %d done.\n', f);
end

% ===== Fit-input table (calculation output; NO fitting yet) =====
% One row per would-be fit: spatial (correlation vs dendritic distance) and
% temporal (correlogram vs lag) for each condition / dendrite-class pair.
% Xdata/Ydata are the values to be fed to the fit in the next cell.
CorrelationData = table( ...
    CD_f(:), string(CD_type(:)), string(CD_class(:)), string(CD_cond(:)), CD_dA(:), CD_dB(:), ...
    CD_x(:), CD_y(:), ...
    'VariableNames', {'FileInd','FitType','ClassPair','Condition','DendClassA','DendClassB', ...
                      'Xdata','Ydata'});

save(corrDataFile, ...
    'corrMat', 'condName', 'heatCond', 'temporalCorr', 'temporalCorrMat', 'lagsCorr', ...
    'CorrelationData', 'dClass2temporal', 'title_str', 'classNames', ...
    'nPair', 'nCond', 'maxLag', 'foi', '-v7.3');
disp('CorrMat is saved');
else
load(corrDataFile);
end

%% ################################################################
%  REVIEWER 3: short-timescale subthreshold dynamics (single neuron f=20)
%  Q: does 20-ms averaging hide faster dynamics, and is the extra
%     high-frequency content blood-vessel contamination?
%
%  This block (f=20 only) produces:
%   1. Subthreshold correlation matrices at 5 / 10 / 20 ms  -> are they different?
%   2. RMS of each smoothing vs the raw trace              -> how much is smoothed away
%   3. Cumulative PCA variance vs component, per window    -> shorter window = more
%                                                            variance in later PCs?
%   4. Movie-level proof it is vascular: high-frequency residual spatial-variance
%      map + neuron-masked spatial PCs vs Result.bvMask / bv-regression component
%   5. Correlation matrix from the blood-vessel-corrected traces
%      (Analysis_20260709_BloodVesselTest) vs the 20-ms matrix on the same frames
% ################################################################
f  = 20;
Pr = PCresult{f};
ch1 = Pr.NormalizedTrace_ch{1};          % per-channel traces (sorted, NaN at motion/dirt)
ch2 = Pr.NormalizedTrace_ch{2};
trC = Pr.NormalizedTrace_dirt;           % combined trace (sorted)
[nROI, nTime] = size(trC);
spk = Pr.SpikeMat(1,:) > 0;              % soma spike train (NaN->false)
dax = Pr.dendaxis;                       % dendritic position (um), sorted order
Dpair = Pr.interDendDist(Pr.sortDist, Pr.sortDist);   % pairwise dendritic distance (um)
Dpair = max(cat(3,Dpair,Dpair'),[],3);

minPairs = 100;
windows  = [5 10 20];                    % smoothing windows (ms == frames @ 1 kHz)
win_all  = [1 windows];                  % include raw (=1) as reference
nW = numel(win_all);
wlab = arrayfun(@(w) sprintf('%d ms',w), win_all, 'uni', 0); wlab{1}='raw';
cmapW = lines(nW);

% good frames = blue-stim off (motion/dirt NaNs are handled by the nan-aware corr)
blueOn  = imdilate(Pr.BlueStim>0, ones(1,201));
goodvec = ~blueOn;

% subthreshold at each window (per channel + combined); keep motion/dirt as NaN
sub1 = cell(1,nW); sub2 = cell(1,nW); subC = cell(1,nW);
for wi = 1:nW
    w = win_all(wi);
    sub1{wi} = get_subthreshold(ch1, spk, 7, w); sub1{wi}(isnan(ch1)) = NaN;
    sub2{wi} = get_subthreshold(ch2, spk, 7, w); sub2{wi}(isnan(ch2)) = NaN;
    subC{wi} = get_subthreshold(trC, spk, 7, w); subC{wi}(isnan(trC)) = NaN;
end

%% (1) Correlation matrices: non-filtered (raw) / 5 / 10 / 20 ms  (cross-channel diagonal, zero lag)
Cmat_w = cell(1,nW);
for k = 1:nW
    [Cm,~] = get_corrTensor(sub1{k}, sub2{k}, subC{k}, goodvec, 0, minPairs);
    Cmat_w{k} = Cm(:,:,1);
end
i20 = find(win_all==20);

% ----- leaf order: hierarchical clustering of the 20-ms correlation matrix -----
% (same leaf order applied to every window so blocks line up; missing pairs -> max distance)
C20  = Cmat_w{i20};
Dd   = 1 - C20;  Dd = (Dd+Dd')/2;  Dd(1:nROI+1:end) = 0;  Dd(isnan(Dd)) = 1;
dvec = squareform(Dd);
Zc   = linkage(dvec, 'average');
try
    ord = optimalleaforder(Zc, dvec);            % Stats toolbox
catch
    fT = figure('Visible','off'); [~,~,ord] = dendrogram(Zc, 0); close(fT);
end

cmapCorr = gen_colormap([0 0.5 1; 1 1 1; 1 0 0], 256);   % blue - white - red
figure('Name','R3-corrMats','NumberTitle','off'); clf
tl = tiledlayout(1,nW,'TileSpacing','compact','Padding','compact');
for k = 1:nW
    nexttile(k);
    imagesc(Cmat_w{k}(ord,ord), [-0.25 1]); axis image; colormap(gca,cmapCorr);
    title(wlab{k}); if k==1, ylabel('ROI (leaf order)'); end
end
cb=colorbar; cb.Layout.Tile='east';
title(tl, sprintf('f=%d  subthreshold correlation (leaf order): non-filtered vs smoothing window', f));

% quantify similarity + correlation-vs-distance
triu1 = triu(true(nROI),1);
fprintf('\n[R3] Correlation-matrix similarity (off-diagonal Pearson r vs 20 ms):\n');
for k = 1:nW
    if k==i20, continue; end
    a = Cmat_w{k}(triu1); b = Cmat_w{i20}(triu1); ok = ~isnan(a)&~isnan(b);
    fprintf('   %-5s vs 20 ms : r = %.4f , mean|Δ| = %.4f\n', ...
        wlab{k}, corr(a(ok),b(ok)), mean(abs(a(ok)-b(ok))));
end

figure('Name','R3-corrVsDist','NumberTitle','off'); clf; hold on
xd = Dpair(triu1);
for k = 1:nW
    yc = Cmat_w{k}(triu1); ok = ~isnan(xd)&~isnan(yc);
    [xs,si] = sort(xd(ok)); ys = yc(ok); ys = ys(si);
    plot(xs, movmean(ys,15), '-', 'LineWidth',1.5, 'DisplayName',wlab{k});
end
xlabel('dendritic distance (\mum)'); ylabel('pairwise correlation'); legend; grid on
title(sprintf('f=%d  correlation vs distance', f));

%% (1b) 20-ms subthreshold correlation matrix over time (thirds, f=20)
% Reviewer: does bleaching / drift over the recording change the subthreshold
% dynamics? Split the recording into three equal time-thirds and recompute the
% 20-ms cross-channel correlation matrix in each (same leaf order + good frames).
% If the dynamics are stable, the three matrices (and correlation-vs-distance)
% overlap.
edgesT     = round(linspace(0, nTime, 4));
Cmat_third = cell(1,3);   thirdMask = false(3, nTime);
for t = 1:3
    thirdMask(t, edgesT(t)+1:edgesT(t+1)) = true;
    cm = get_corrTensor(sub1{i20}, sub2{i20}, subC{i20}, goodvec & thirdMask(t,:), 0, minPairs);
    Cmat_third{t} = cm(:,:,1);
end
tlab3 = {'1st third','2nd third','3rd third'};

figure('Name','R3-corrMats-thirds','NumberTitle','off'); clf
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for t = 1:3
    nexttile(t);
    imagesc(Cmat_third{t}(ord,ord), [-0.25 1]); axis image; colormap(gca,cmapCorr);
    title(sprintf('%s (%d frm)', tlab3{t}, nnz(goodvec & thirdMask(t,:))));
    if t==1, ylabel('ROI (leaf order)'); end
end
cb = colorbar; cb.Layout.Tile = 'east';
title(tl, sprintf('f=%d  20-ms subthreshold correlation over time (thirds)', f));

triu1 = triu(true(nROI),1);
fprintf('\n[R3] 20-ms correlation-matrix similarity across time-thirds (off-diagonal Pearson r):\n');
prs = [1 2; 2 3; 1 3];
for p = 1:size(prs,1)
    a = Cmat_third{prs(p,1)}(triu1);  b = Cmat_third{prs(p,2)}(triu1);
    ok = ~isnan(a) & ~isnan(b);
    fprintf('   third %d vs %d : r = %.4f , mean|Δ| = %.4f\n', ...
        prs(p,1), prs(p,2), corr(a(ok),b(ok)), mean(abs(a(ok)-b(ok))));
end

figure('Name','R3-corrVsDist-thirds','NumberTitle','off'); clf; hold on
xd = Dpair(triu1);
for t = 1:3
    yc = Cmat_third{t}(triu1);  ok = ~isnan(xd) & ~isnan(yc);
    [xs,si] = sort(xd(ok));  ys = yc(ok);  ys = ys(si);
    plot(xs, movmean(ys,15), '-', 'LineWidth',1.5, 'DisplayName',tlab3{t});
end
xlabel('dendritic distance (\mum)'); ylabel('pairwise correlation'); legend; grid on
title(sprintf('f=%d  20-ms correlation vs distance over time', f));

%% (1c) Peak-triggered average over time (thirds, f=23)
% Reviewer (bleaching/drift): recompute the subthreshold peak-triggered average in
% each time-third. Peaks are detected as in Fig. 2 (band-stop-filtered, per-ROID
% subthreshold, excluding peri-spike and blue-on frames); for each ROI we take the
% auto (same-ROI) peak-triggered subthreshold waveform, then average across ROIs.
% Stable amplitude/shape across thirds => no bleaching-driven change in dynamics.
fp   = 23;
Prp  = PCresult{fp};
[nRp, nTp] = size(Prp.Subthreshold);
SubN = interpolateNaN(Prp.Subthreshold ./ Prp.SpikeHeight);        % spike-height units
StopFreq = [12 200];

perispk = unique([tovec(find(Prp.SpikeMat(1,:))' + (-10:40)); find(Prp.CStrace)']);
perispk(perispk < 1 | perispk > nTp) = [];
badp = ind2vec(nTp, perispk, 1) > 0 | imdilate(Prp.BlueStim > 0, ones(1,201));

peakvec = zeros(nRp, nTp);
for n = 1:nRp
    ft = get_bandstop(SubN(n,:), 1000, StopFreq);
    ft = ft - movmedian(ft, 300, 2);
    [~, locs] = findpeaks(ft, 'MinPeakHeight',0.2, 'MinPeakDistance',50, 'MinPeakProminence',0.1);
    pv = ind2vec(nTp, locs, 1);  pv(badp) = 0;
    peakvec(n,:) = pv;
end

edgesTp = round(linspace(0, nTp, 4));
preP = 300;  postP = 300;  tauP = (-preP:postP);
labs3 = {'1st third','2nd third','3rd third'};
staThird = cell(1,3);  nPk = zeros(1,3);
for t = 1:3
    inThird = false(1,nTp);  inThird(edgesTp(t)+1:edgesTp(t+1)) = true;
    acc = [];
    for n = 1:nRp
        pv = peakvec(n,:) > 0 & inThird;
        if nnz(pv) < 3, continue; end
        sta = get_STA(SubN(n,:), pv, preP, postP);            % 1 x L, auto (same ROI)
        sta = sta - mean(sta(tauP < -100));                   % pre-peak baseline
        acc = [acc; sta];   %#ok<AGROW>
        nPk(t) = nPk(t) + nnz(pv);
    end
    staThird{t} = acc;
end

figure('Name','R3-PTA-thirds','NumberTitle','off'); clf; hold on
cmap3 = lines(3);  h3 = gobjects(1,3);
for t = 1:3
    S = staThird{t};
    m = mean(S,1,'omitnan');  s = std(S,0,1,'omitnan') ./ sqrt(size(S,1));
    h3(t) = errorbar_shade(tauP, m, s, cmap3(t,:));  hold on;
end
xlabel('time from subthreshold peak (ms)'); ylabel('subthreshold (spike height)');
legend(h3, arrayfun(@(t) sprintf('%s (%d peaks)', labs3{t}, nPk(t)), 1:3, 'uni',0), 'box','off');
title(sprintf('f=%d  peak-triggered average over time (thirds)', fp));
xlim([-preP postP]); grid on;

fprintf('\n[R3] f=%d peak-triggered average amplitude per time-third:\n', fp);
for t = 1:3
    m = mean(staThird{t},1,'omitnan');
    fprintf('   %s: peak amp = %.3f spike-height  (%d peaks, %d ROIs)\n', ...
        labs3{t}, max(m), nPk(t), size(staThird{t},1));
end

%% (2) How much of the raw trace does smoothing change? (% of variance)
% Per ROI and window w, compare the raw (non-filtered) subthreshold trace with its
% w-ms moving average, over analysed (blue-off) frames, as a fraction of the raw
% trace's variance (power):
%     retained(w) = var(smoothed_w) / var(raw)     % what the moving average keeps
%     removed(w)  = 100% - retained(w)             % the fast part it discards
% (raw = window 1, so removed(1) = 0 and grows with the window). "removed" is the
% percentage of the original subthreshold trace that w-ms smoothing takes out --
% the short-timescale content the reviewer asks about. Sections (3)-(4) test how
% much of it is vascular. Both curves are plotted (they sum to 100%).
iRaw   = find(win_all==1);
Vraw   = var(subC{iRaw}(:,goodvec), 0, 2, 'omitnan');    % nROI x 1, raw variance per ROI
retMat = nan(nROI, nW);
for wi = 1:nW
    Vw = var(subC{wi}(:,goodvec), 0, 2, 'omitnan');       % variance of the smoothed trace
    retMat(:,wi) = 100 * Vw ./ Vraw;                      % % retained
end
remMat = 100 - retMat;                                    % % removed

figure('Name','R3-RMS','NumberTitle','off'); clf; hold on
errorbar(win_all, mean(remMat,1,'omitnan'), std(remMat,0,1,'omitnan')./sqrt(nROI), '-o','LineWidth',1.5,'DisplayName','removed by smoothing');
errorbar(win_all, mean(retMat,1,'omitnan'), std(retMat,0,1,'omitnan')./sqrt(nROI), '-o','LineWidth',1.5,'DisplayName','retained (smoothed)');
xlabel('moving-average window (ms)   [1 = non-filtered]');
ylabel('% of raw subthreshold variance');
grid on; legend('Location','east'); ylim([0 100]);
title(sprintf('f=%d  fraction of the subthreshold trace changed by smoothing', f));

mRem = mean(remMat,1,'omitnan');
fprintf('\n[R3] %% of raw subthreshold variance REMOVED by smoothing (mean over ROIs):\n');
for wi = 1:nW, fprintf('   %-5s : %5.1f%% removed  (%5.1f%% retained)\n', wlab{wi}, mRem(wi), 100-mRem(wi)); end

%% (3) Cumulative PCA variance vs component, per window
figure('Name','R3-PCA','NumberTitle','off'); clf
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
axC = nexttile; hold(axC,'on'); title(axC,'cumulative variance');
axV = nexttile; hold(axV,'on'); title(axV,'per-PC variance'); set(axV,'YScale','log','Xscale','log');
for wi = 1:nW
    X = subC{wi}; X(:,~goodvec) = NaN;
    vc = all(~isnan(X),1);                 % complete columns only
    Xv = X(:,vc); Xv = Xv - mean(Xv,2);
    s  = svd(Xv,'econ');
    ve = s.^2 / sum(s.^2);
    plot(axC, cumsum(ve), '-o','Color',cmapW(wi,:),'DisplayName',wlab{wi});
    plot(axV, ve,         '-o','Color',cmapW(wi,:),'DisplayName',wlab{wi});
end
xlabel(axC,'PC #'); ylabel(axC,'cumulative variance explained'); grid(axC,'on'); legend(axC,'Location','se');
xlabel(axV,'PC #'); ylabel(axV,'variance explained'); grid(axV,'on'); xlim(axC,[1 min(nROI,20)]); xlim(axV,[1 min(nROI,20)]);
sgtitle(sprintf('f=%d  subthreshold PCA vs smoothing window (shorter window -> variance spreads to later PCs = blood)', f));

%% (4) Movie-level: is the high-frequency content vascular?  (f=20 segment)
fpath2readmov = '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20240807/201405BHLm141_N1_VR_LowStim';
Resultsmall=importdata("/Volumes/BHL18TB_D2/20240807/201405BHLm141_N1_VR_LowStim/PC_Result_smallerROI.mat");
dMask=point2img(Result.SWC,3,size(Result.ref_im));
segFrames = 78000:88000;                       % == actualframe(1:2000)
nrow = size(Pr.ref_im,1);  ncol = size(Pr.ref_im,2);   % (180, 552)

%mov_mc = double(readBinMov(fullfile(fpath2readmov,'mov_mc_segment.bin'), nrow, ncol));
mov_mc=readBinMov_BHL_multiple(fpath{f},3,segFrames ,15000,'mc_ShutterReg');
mov_mc = mov_mc(:,:,1:numel(segFrames));
%S = load(fullfile(fpath2readmov,'mcTrace_seg.mat'),'mcTrace'); mcTrace = S.mcTrace(1:numel(segFrames),:);
mcTrace=Result.mcTrace(segFrames ,:);

ftprnt = Pr.Ftprnts;      % (180,552,nROIall)
bvMask = Pr.bvMask;       % (180,552,nBV)

% blood-vessel temporal trace (regress out the neuron footprint), as in BloodVesselTest
bv_trace = tovec(mov_mc)' * tovec(bvMask);
V_trace  = tovec(mov_mc)' * tovec(ftprnt);
bv_trace = squeeze(SeeResiduals(permute(bv_trace,[2 3 1]), V_trace))';

% motion residual movie (spike-positive convention handled later by sign)
mov_res  = mov_mc - mean(mov_mc,3);
mov_res  = SeeResiduals(mov_res, mcTrace);
mov_res  = SeeResiduals(mov_res, mcTrace.^2);
mov_res  = SeeResiduals(mov_res, mcTrace(:,1).*mcTrace(:,2));
mov_res2 = SeeResiduals(mov_res, bv_trace);        % blood-vessel-regressed movie

% high-frequency residual (what 20-ms smoothing throws away) and its spatial variance
hf        = mov_res - movmean(mov_res, 20, 3);
hfvar     = var(hf, 0, 3);
bloodComp = mov_res - mov_res2;                     % component removed by bv-regression
bloodVar  = var(bloodComp, 0, 3);
bvMap     = max(bvMask, [], 3);
neuronMap = max(ftprnt, [], 3);

% blood-vessel mask outline
bvBW = bvMap > 0.2*max(bvMap(:));

figure('Name','R3-spatialVar','NumberTitle','off'); clf
nexttile([1 1]);
imagesc(hfvar); axis image off; colormap(turbo); colorbar; hold on;
contour(bvBW, [0.5 0.5], 'g', 'LineWidth', 1);
title(sprintf('f=%d  HF residual variance (raw - 20ms) + bvMask contour', f));

nexttile([1 1]);
imagesc(bloodVar); axis image off; colormap(turbo); colorbar; hold on;
contour(bvBW, [0.5 0.5], 'g', 'LineWidth', 1);
title(sprintf('f=%d  bv component + bvMask contour', f));

bvpix = bvMap(:) > 0;
fprintf('\n[R3] HF-variance spatial map vs vascular reference:\n');
fprintf('   corr2(HF var, bv-regression var) = %.3f\n', corr2(hfvar, bloodVar));
fprintf('   HF var in bvMask pixels / outside = %.2f\n', mean(hfvar(bvpix))/mean(hfvar(~bvpix)));

ftprintsmallest=divide_ROIs(Resultsmall.ftprint_smaller.*dMask,30,dMask);

%% (4c) Small-footprint fast-dynamics: is what the 20-ms average removes vascular?
% ftprintsmallest: [H x W x M] NON-overlapping small sub-ROIs (~30 px, a few um) tiled along
% the dendrite. Per sub-ROI we dot the motion-residual movie mov_res with (all) the
% footprint and its (odd)/(even) checkerboard pixel halves. odd/even share the same signal
% but have INDEPENDENT shot noise.
%   Q1 fast component = raw - movmean(raw,20) = EXACTLY what the paper's 20-ms subthreshold discards.
%   Q2 distance = geodesic (contour) length along the skeleton (geodesic_distance), not Euclidean.
%   Q3 spatiotemporal cross-correlation C(i,j,tau): AUTO (i==j) from odd x even (shot-noise-free);
%      CROSS (i!=j) from the raw trace.
%   Q4 fast variance vs spatial-length bins uses odd/even (shot-noise-free).
% Hypothesis: the <20 ms content is blood -> large variance but coherent only over a few um.
fs   = 1000;                                   % sampling rate (Hz)
Msub = size(ftprintsmallest,3);
dMaskSkel=point2img(Result.SWC,3,size(Result.ref_im));

% --- footprints all / odd / even (checkerboard), unit-sum normalized (spatial mean) ---
[ccg,rrg] = meshgrid(1:ncol,1:nrow);  oddPix = mod(rrg+ccg,2)==1;
Fall = double(ftprintsmallest);
Fodd = Fall .*  oddPix;   Feven = Fall .* ~oddPix;
nrm  = @(F) F ./ max(sum(sum(F,1),2), eps);
Fall = nrm(Fall);  Fodd = nrm(Fodd);  Feven = nrm(Feven);

% --- extract traces from mov_res  (M x T), spike-positive dF ---
% mov_res is the MOTION-RESIDUAL movie (mean image subtracted, mcTrace / mcTrace^2 /
% cross term regressed out), so these traces are dF in camera counts with the motion
% artifact already removed -- the same movie the blood analysis below operates on.
% Sign is flipped (-mov_res) to the usual spike-positive voltage convention, as in
% loadVoltChunk. Every downstream quantity here is a correlation / covariance between
% pairs of traces, so the flip does not change any of them.
Xres  = -tovec(mov_res);                        % pixels x T
Tall  = (Xres' * tovec(Fall ))';                % M x T
Todd  = (Xres' * tovec(Fodd ))';
Teven = (Xres' * tovec(Feven))';

% --- frames to drop: peri-spike (+/-Wsp, avoids filter ringing) and blue-on (global) ---
Wsp   = 25;
spk   = full(Result.spike(1, segFrames)) > 0;
badfr = imdilate(spk, ones(1,2*Wsp+1)) > 0;
blue  = Result.Blue(:)';  if numel(blue)>=max(segFrames), badfr = badfr | blue(segFrames)>0; end

% --- Geodesic (contour) distance along the skeleton between sub-ROIs (um)  [Q2] ---
% Path is constrained to dMaskSkel, so distance follows the neurite, not a straight line.
% Computed once and cached (keyed by sub-ROI count + centroids); reloaded next time.
coordS  = get_coord(ftprintsmallest);           % M x 2 [x y]=[col row] px
triuM   = triu(true(Msub),1);
geoFile = fullfile(fpath{f}, sprintf('geodist_smallROI_M%d.mat', Msub));

Dgeo_px = [];
if isfile(geoFile)
    Sg = load(geoFile, 'Dgeo_px', 'coordS');
    if isequal(size(Sg.Dgeo_px),[Msub Msub]) && isequal(size(Sg.coordS),size(coordS)) ...
            && max(abs(Sg.coordS(:)-coordS(:))) < 1e-6
        Dgeo_px = Sg.Dgeo_px;
        fprintf('[R3] loaded cached geodesic distances (%s)\n', geoFile);
    else
        fprintf('[R3] cache %s does not match current sub-ROIs -> recomputing\n', geoFile);
    end
end
if isempty(Dgeo_px)
    Dgeo_px = zeros(Msub);
    fprintf('[R3] computing geodesic distances for %d sub-ROIs (one time, ~a minute) ...\n', Msub);
    for ii = 1:Msub
        for jj = ii+1:Msub
            Dgeo_px(ii,jj) = geodesic_distance(dMaskSkel, coordS(ii,:), coordS(jj,:));
        end
    end
    Dgeo_px = Dgeo_px + Dgeo_px';                % symmetric, pixels
    save(geoFile, 'Dgeo_px', 'coordS');
    fprintf('[R3] saved geodesic distances -> %s\n', geoFile);
end
Dgeo = Dgeo_px * Pixelsize(f);                  % um

% --- FAST (<20 ms) component = exactly what the paper's 20-ms average removes  [Q1] ---
% subthreshold = movmean(trace,20); the discarded fast part = raw - movmean(raw,20).
% We analyse THIS component, so the statistics describe precisely what the reviewer
% claims the 20-ms average throws away (a boxcar high-pass, matched to the paper).
Wsm    = 20;                                     % ms == frames @ 1 kHz  (paper's window)
fTall  = Tall  - movmean(Tall , Wsm, 2);
fTodd  = Todd  - movmean(Todd , Wsm, 2);
fTeven = Teven - movmean(Teven, Wsm, 2);

% --- spatiotemporal ODD/EVEN cross-correlation tensor  C(i,j,tau)  [Q3] ---
% EVERY pair (incl. the diagonal) is the standard Pearson correlation of odd_i and even_j:
% C(i,j,tau) = corr(odd_i(t), even_j(t+tau)), normalized by each trace's own variance
% (nanXCorrFFT, normalized=true). Same normalization for auto and cross -> distance 0 and
% distance>0 are on the same footing (fixes the earlier auto-vs-raw mismatch).
maxLag = 30;  lags = (-maxLag:maxLag)';  lagms = lags(:)';  nLag = numel(lags);  tau0 = maxLag+1;
minPairsX = 100;
Om = fTodd;  Om(:,badfr) = NaN;      % odd,  bad frames -> NaN (nanXCorrFFT handles masking + centering)
Em = fTeven; Em(:,badfr) = NaN;      % even
Cxt = nan(Msub, Msub, nLag);
for i = 1:Msub
    for j = i:Msub
        xc = nanXCorrFFT(Om(i,:), Em(j,:), maxLag, true, minPairsX);   % corr(odd_i, even_j)
        Cxt(i,j,:) = xc(:);
        if i ~= j, Cxt(j,i,:) = flip(xc(:)); end                       % reciprocity C(j,i,tau)=C(i,j,-tau)
    end
end

% --- average correlogram binned by geodesic distance -> (distance x tau)  [Q4 fig] ---
distEdges = 0:2:60;  distCtr = (distEdges(1:end-1)+distEdges(2:end))/2;  nDist = numel(distCtr);
diagC = zeros(Msub,nLag);
for i = 1:Msub, diagC(i,:) = squeeze(Cxt(i,i,:))'; end   % odd/even autocorrelograms (distance 0)
Ci   = reshape(Cxt, Msub*Msub, nLag);
dv   = Dgeo(triuM);   Cpair = Ci(triuM(:), :);           % off-diagonal (cross) pairs
heatDT = nan(nDist+1, nLag);
heatDT(1,:) = mean(diagC, 1, 'omitnan');                 % row 1 = distance 0 (auto)
for bdi = 1:nDist
    sel = dv>=distEdges(bdi) & dv<distEdges(bdi+1);
    if any(sel), heatDT(bdi+1,:) = mean(Cpair(sel,:), 1, 'omitnan'); end
end
distAxis = [0 distCtr];

% ---- correlation colormap (blue-white-red) + clim for the fast maps ----
ccFast = gen_colormap([0 0.5 1; 1 1 1; 1 0 0], 256);   ccLim = [-0.1 0.3];
% (heatmap / tau-profile / distance-profile are drawn in the combined figure R3-fast-summary below)

% ---- tau=0 correlation matrix ordered by contour (geodesic) distance from soma  [Q1] ----
% Reveals where the distance-profile peak comes from (which sub-ROI pairs are correlated).
R0 = Cxt(:,:,tau0);
somaXY = get_coord(Result.ftprnt(:,:,1));            % soma centroid [x y] (Result ROI #1)
[~, somaSub] = min(sum((coordS - somaXY).^2, 2));    % sub-ROI nearest the soma
[dsoma, ordSoma] = sort(Dgeo(somaSub,:));            % order sub-ROIs by geodesic distance from soma
figure('Name','R3-xcorr-tau0-matrix','NumberTitle','off'); clf
imagesc(R0(ordSoma,ordSoma)); axis image; colormap(ccFast); caxis(ccLim); colorbar;
ti = unique(round(linspace(1,Msub,8)));
set(gca,'XTick',ti,'XTickLabel',round(dsoma(ti)),'YTick',ti,'YTickLabel',round(dsoma(ti)));
xlabel('contour distance from soma (\mum)'); ylabel('contour distance from soma (\mum)');
title(sprintf('f=%d  \\tau=0 fast correlation, ordered by contour distance from soma (sub-ROI %d)', f, somaSub));

% ---- tau=0 seed-correlation MOVIE projected onto the morphology (+ bvMask contour) ----
% Frame s paints every sub-ROI footprint with its tau=0 fast correlation to seed s, and
% sweeps the seed across all sub-ROIs. Blood -> a small blob that stays near the seed.
R0  = Cxt(:,:,tau0);                             % M x M  tau=0 correlation (auto=odd/even, cross=raw)
movieFile = fullfile(fpath{f}, sprintf('f%d_tau0_fastcorr_seedmovie.mp4', f));
if isfile(movieFile)
    fprintf('[R3] seed-correlation movie already exists, skipping (%s)\n', movieFile);
else
    ftL = double(ftprintsmallest);  den = sum(ftL,3);
    [~, labelImg] = max(ftL, [], 3);  inMask = den > eps;  labelImg(~inMask) = 0;   % pixel -> sub-ROI (non-overlapping)
    vw = VideoWriter(movieFile, 'MPEG-4');  vw.FrameRate = 10;  vw.Quality = 100;  open(vw);
    figM = figure('Color','w','WindowStyle','normal','Position',[100 100 950 340]);
    for s = 1:Msub
        val  = R0(:,s);
        cmap = nan(nrow,ncol);  cmap(inMask) = val(labelImg(inMask));
        clf(figM);  ax = axes(figM);
        hh = imagesc(ax, cmap);  set(hh,'AlphaData',~isnan(cmap));  axis(ax,'image','off');
        colormap(ax, ccFast);  caxis(ax, ccLim);  colorbar(ax);  hold(ax,'on');
        contour(ax, bvBW, [0.5 0.5], 'k', 'LineWidth', 1);
        plot(ax, coordS(s,1), coordS(s,2), 'o','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',[0 1 0]);
        title(ax, sprintf('f=%d   seed %d/%d   \\tau=0 fast correlation', f, s, Msub));
        drawnow;  writeVideo(vw, getframe(figM));
    end
    close(vw);  close(figM);
    fprintf('[R3] saved tau=0 seed-correlation movie -> %s\n', movieFile);
end

% --- Q4: shared variance vs spatial-length bins -- RAW vs HIGH-PASS, violin distributions ---
% Per ROI pair, the tau=0 shot-noise-free cross-covariance cov(odd_i,even_j) (variance units).
%   raw       : full-band mov_res trace (residual bleach removed by a linear detrend) = all timescales
%   high-pass : raw - movmean(raw,20) = the fast (<20 ms) part the paper's average removes
% For each geodesic-length bin we show the DISTRIBUTION over pairs as a violin + dots.
lenEdges = [0 5 10 20 50 100];  lenLab = {'<5','5-10','10-20','20-50','50-100'};  nLen = numel(lenLab);

% good-frame odd/even traces (mean-subtracted). raw is also linearly detrended (removes bleach).
gO_raw  = detrend(Todd(:,~badfr)','linear')';   gE_raw  = detrend(Teven(:,~badfr)','linear')';
gO_fast = fTodd(:,~badfr) - mean(fTodd(:,~badfr),2);   gE_fast = fTeven(:,~badfr) - mean(fTeven(:,~badfr),2);
M2 = size(gO_raw,2);
Kraw  = (gO_raw *gE_raw' )/(M2-1);   Kraw  = (Kraw +Kraw' )/2;   % raw       cross-cov (odd,even)
Kfast = (gO_fast*gE_fast')/(M2-1);   Kfast = (Kfast+Kfast')/2;   % high-pass cross-cov (odd,even)
covRaw  = Kraw(triuM);   covFast = Kfast(triuM);                 % per-pair, aligned with dv

% build interleaved violin groups: [bin1 raw][bin1 HP][bin2 raw][bin2 HP]...
Xcell = {};  vcolors = [];  vlab = {};
rawC = [0.60 0.60 0.60];  fastC = [0.85 0.20 0.20];
for L = 1:nLen
    sel = dv>=lenEdges(L) & dv<lenEdges(L+1);
    Xcell{end+1} = covRaw(sel);   vcolors(end+1,:) = rawC;   vlab{end+1} = [lenLab{L} ' raw'];   %#ok<SAGROW>
    Xcell{end+1} = covFast(sel);  vcolors(end+1,:) = fastC;  vlab{end+1} = [lenLab{L} ' HP'];    %#ok<SAGROW>
end

figure('Name','R3-fastVar-vs-length','NumberTitle','off'); clf
pV = Violin_wPoints(Xcell, vcolors);            % violins + dots (+ ranksum p-values, printed)
set(gca,'XTick',1:numel(Xcell),'XTickLabel',vlab); xtickangle(45);
ylabel('cross-covariance @ \tau=0  (odd/even, variance units)'); grid on; yline(0,'k:');
title(sprintf('f=%d  shared variance vs geodesic length:  raw (gray) vs high-pass <20 ms (red)', f));

fprintf('\n[R3] Shared variance (tau=0 cross-cov) per geodesic-length bin  [median | raw vs HP ranksum p]:\n');
for L = 1:nLen
    sel = dv>=lenEdges(L) & dv<lenEdges(L+1);
    fprintf('   %-7s um (n=%5d): raw=%.3g  HP=%.3g   p=%.2g\n', lenLab{L}, nnz(sel), ...
        median(covRaw(sel),'omitnan'), median(covFast(sel),'omitnan'), pV(2*L-1,2*L));
end

% --- HP/raw variance ratio for coarse spatial ranges (raw = denominator)  [Q2] ---
% mean cross-covariance (variance units) in each range; ratio = HP / raw.
ratEdges = {[0 20],[20 100],[100 inf]};  ratLab = {'<20','20-100','>100'};  nR = numel(ratLab);
mRaw = nan(1,nR); mHP = nan(1,nR);
for r = 1:nR
    sel = dv>=ratEdges{r}(1) & dv<ratEdges{r}(2);
    mRaw(r) = mean(covRaw(sel),'omitnan');  mHP(r) = mean(covFast(sel),'omitnan');
end
ratio_HPraw = mHP ./ mRaw;      % fraction of shared variance that is fast (<20 ms)

fprintf('\n[R3] HP/raw shared-variance ratio (mean cross-cov, raw = denominator):\n');
for r = 1:nR
    fprintf('   %-7s um : raw=%.3g  HP=%.3g   HP/raw=%.2f\n', ratLab{r}, mRaw(r), mHP(r), ratio_HPraw(r));
end

% --- fast (<20 ms) variance decomposition: blood (local) vs synaptic (30-100 um) ---
% tau=0 shot-noise-free cross-covariance C(d) (odd/even, variance units). Nested layers
% telescope to the total fast variance Vtot = C(0):
%   independent = Vtot - C(<10)      not shared even at the nearest ROI
%   BLOOD       = C(<10) - C(30-100) coherent within ~10 um, gone by 30 um (vascular scale)
%   SYNAPTIC    = C(30-100) - C(>120) coherent at synaptic scale
%   global      = C(>120)            long-range baseline (~0)
autoVarFast = sum(gO_fast.*gE_fast,2)/(M2-1);            % per-ROI fast variance (odd/even)
Vtot = mean(autoVarFast,'omitnan');
C10  = mean(covFast(dv<20),          'omitnan');
C30  = mean(covFast(dv>=20 & dv<100),'omitnan');
C120 = mean(covFast(dv>=100),        'omitnan');
% blood/local folds the sub-ROI (independent) variance in, since a sub-ROI is itself < ~10 um.
comp    = [Vtot-C30, C30-C120 C120];                   % [blood/local  synaptic  global]
compLab = {'blood/local (<~20 µm, incl. sub-ROI)','synaptic (20-100 µm)','global (>100 µm)'};

fb = 0:3:150; fc = (fb(1:end-1)+fb(2:end))/2;           % C(d) curve (cross-cov vs distance)
Cd = arrayfun(@(k) mean(covFast(dv>=fb(k)&dv<fb(k+1)),'omitnan'), 1:numel(fc));

compPos = max(comp,0);   pctBlood = 100*compPos(1)/sum(compPos);   pctSyn = 100*compPos(2)/sum(compPos);

%% ---- Combined fast (<20 ms) dynamics summary figure (6 panels) ----
figure('Name','R3-fast-summary','NumberTitle','off'); clf
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% (i) heatmap  distance x tau
% ax=nexttile; imagesc(lagms, distAxis, heatDT); set(ax,'YDir','normal'); colormap(ax,ccFast); caxis(ax,ccLim); colorbar(ax);
% xlabel('\tau (ms)'); ylabel('Contour distance (µm)'); %title('(i) cross-corr  (distance \times \tau)');

heatDTnorm=heatDT./heatDT(:,lagms==0);
ax=nexttile; imagesc(lagms, distAxis, heatDTnorm); set(ax,'YDir','normal'); colormap(ax,ccFast); caxis(ax,ccLim); colorbar(ax);
xlabel('\tau (ms)'); ylabel('Contour distance (µm)'); %title('(i) cross-corr  (distance \times \tau)');

% (ii) tau profile at distance 0
nexttile; plot(lagms, heatDT(1,:), '-','LineWidth',1.5); grid on;
xlabel('\tau (ms)'); ylabel('Autocorrelation'); %title('(ii) \tau profile @ distance 0'); box off;

% (iii) distance profile at tau=0
nexttile; plot(distAxis, heatDT(:,tau0), '-o','LineWidth',1.5); grid on;
xlabel('Distance (µm)'); ylabel('Correlation'); %title('(iii) distance profile @ \tau=0'); box off;
set(gca,'xscale','log');

% (iv) C(d) with spatial bands
nexttile; hold on
hD = plot([0 fc], [Vtot Cd], '-o','Color',[.15 .15 .15],'LineWidth',1.3);
yl = ylim;  xr = [0.1 20; 20 100; 100 400];  ccb = num2cell(lines(3),2);
for r = 1:3, patch([xr(r,1) xr(r,2) xr(r,2) xr(r,1)],[yl(1) yl(1) yl(2) yl(2)], ccb{r},'EdgeColor','none','FaceAlpha',.4); end
uistack(hD,'top'); ylim(yl); yline(0,'k:');
xlabel('Distance (µm)'); ylabel('Covariance (A.U.)'); box off;
% y-limits follow the data: the covariance scale depends on which movie the traces
% were extracted from (mov_res vs mov_mc), so it must not be hardcoded.
visD = fc >= 10.^1 & fc <= 10.^2.5;
set(gca,'xscale','log'); xlim(10.^[1 2.5]);
ylim([min(0, min(Cd(visD))) max(Cd(visD))*1.1]);
%title('(iv) C(d): blood | synaptic | baseline');

% (v) blood / synaptic / global pie
nexttile; pie(compPos,[1 0 0]); colormap(gca,lines(3));
legend(compLab,'Location','southoutside');
%title(sprintf('(v) blood %.0f%% / syn %.0f%%', pctBlood, pctSyn));

% (vi) HP/raw ratio (fast fraction of shared variance)
nexttile; bar(categorical(ratLab,ratLab), ratio_HPraw); box off;
ylabel('Covariance ratio (high/raw)'); xlabel('Contour distance (µm)'); %title('(vi) fast fraction  (HP/raw)');

%sgtitle(sprintf('f=%d  fast (<20 ms = raw-movmean20) dynamics summary, %d sub-ROIs', f, Msub));
set_font('Arial'); set_fontsize(14);
set_figsize(280,160);

fprintf('\n[R3] Fast (<20 ms) variance decomposition (odd/even):  Vtot=%.3g\n', Vtot);
for c = 1:3, fprintf('   %-32s = %.3g  (%.0f%% of total)\n', compLab{c}, comp(c), 100*comp(c)/Vtot); end

% ---- ROI map colored by the tau=0 correlation at the ~6 um distance tick ----
% Spatially decomposes the (iii) distance-profile point at ~6 um: each sub-ROI is painted
% with its mean tau=0 correlation to the sub-ROIs that are ~targetD um away (geodesic).
% The average of this map over ROIs = the distance-profile value at that distance.
targetD = 6;  tolD = 2;                                  % um (2nd distance tick, adjust as needed)
perROI6 = nan(Msub,1);
for i = 1:Msub
    nb = find(abs(Dgeo(i,:) - targetD) <= tolD & (1:Msub) ~= i);
    if ~isempty(nb), perROI6(i) = mean(R0(i,nb), 'omitnan'); end
end
ftL6 = double(ftprintsmallest);  den6 = sum(ftL6,3);
[~, labelImg6] = max(ftL6,[],3);  inMask6 = den6 > eps;  labelImg6(~inMask6) = 0;
map6 = nan(nrow,ncol);  map6(inMask6) = perROI6(labelImg6(inMask6));
figure;
figure('Name','R3-corr6um-map','NumberTitle','off'); clf
hh = imagesc(map6);  set(hh,'AlphaData',~isnan(map6));  axis image off;
colormap(ccFast);  caxis(ccLim);  colorbar;  hold on;
contour(bvBW, [0.5 0.5], 'k', 'LineWidth', 1);
title(sprintf('f=%d  per-ROI mean \\tau=0 correlation to %g\\pm%g \\mum neighbors', f, targetD, tolD));

%% (4b) ROI-based variance map: raw vs subthreshold, colored footprints + bvMask contour
% Per ROI, temporal variance (over analysed frames) of:
%   (i)  the raw (non-filtered) subthreshold trace,
%   (ii) the 20-ms subthreshold trace,
%   (iii) their difference (raw - 20 ms) = the fast content removed by 20-ms smoothing.
% Each ROI footprint is painted with its value and the bvMask contour is overlaid.
% If the (raw - 20 ms) map is high on ROIs lying over vessels, the fast content is vascular.
sortDist = Pr.sortDist;
iRaw = find(win_all==1);  i20 = find(win_all==20);
Vraw = var(subC{iRaw}(:,goodvec), 0, 2, 'omitnan');
Vsub = var(subC{i20 }(:,goodvec), 0, 2, 'omitnan');
Vhf  = var(subC{iRaw}(:,goodvec) - subC{i20}(:,goodvec), 0, 2, 'omitnan');
roiSets = {Vraw,'var(raw ROI trace)'; Vsub,'var(20 ms subthreshold)'; Vhf,'var(raw - 20 ms) = removed'};

ftUse = ftprnt(:,:,sortDist);  ftUse(ftUse<0) = 0;  den = sum(ftUse,3);   % footprint weights
cmax  = prctile([Vraw;Vsub], 99);

figure('Name','R3-ROIvar','NumberTitle','off'); clf
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for p = 1:3
    val = roiSets{p,1};
    roiMap = sum(ftUse .* reshape(val,1,1,[]), 3) ./ den;   % footprint-weighted per-ROI value
    roiMap(den < eps) = NaN;
    ax = nexttile;  hh = imagesc(roiMap);  axis image off;  set(hh,'AlphaData',~isnan(roiMap));
    colormap(ax,turbo);  colorbar(ax);  hold(ax,'on');
    contour(ax, bvBW, [0.5 0.5], 'm', 'LineWidth', 1);
    title(roiSets{p,2});
    if p < 3, caxis(ax,[0 cmax]); end
end
sgtitle(sprintf('f=%d  per-ROI variance (colored footprints) + bvMask contour', f));

% neuron-masked movie PCA: (a) full variance spectrum (ALL PCs) for raw/5/20 ms,
% (b) leading spatial-PC images. Shorter window -> variance spreads to later PCs and
% vascular structures appear among the low modes.
maskIdx = find(imdilate(neuronMap>0, strel('disk',3)));
nPCshow = 12;
win_mov = [1 5 20];  wl_mov = {'raw','5 ms','20 ms'};  cmapM = lines(numel(win_mov));

figVE = figure('Name','R3-PCA-movie','NumberTitle','off'); clf
tlM  = tiledlayout(figVE,1,2,'TileSpacing','compact','Padding','compact');
axMc = nexttile(tlM); hold(axMc,'on'); title(axMc,'cumulative variance');
axMv = nexttile(tlM); hold(axMv,'on'); title(axMv,'per-PC variance'); set(axMv,'YScale','log','xscale','log');

for wi = 1:numel(win_mov)
    w = win_mov(wi);
    movsm = movmean(mov_res, w, 3);
    M = tovec(movsm); M = M(maskIdx,:); M = M - mean(M,2);   % masked pixels x time, centered

    % ALL-component variance from the time x time Gram eigenvalues (= singular values^2)
    ev = sort(eig(M'*M), 'descend'); ev = max(ev,0);
    ve = ev / sum(ev);
    plot(axMc, cumsum(ve), '-', 'Color',cmapM(wi,:), 'LineWidth',1.2, 'DisplayName',wl_mov{wi});
    plot(axMv, ve,          '-', 'Color',cmapM(wi,:), 'LineWidth',1.2, 'DisplayName',wl_mov{wi});

    % leading spatial-PC images for this window
    [u,~,~] = svds(M, nPCshow);
    uimg = zeros(nrow*ncol, nPCshow); uimg(maskIdx,:) = u;
    uimg = reshape(uimg, nrow, ncol, nPCshow);
    figure('Name',sprintf('R3-spatialPCA-%s',strrep(wl_mov{wi},' ','')),'NumberTitle','off'); clf
    imshow2_patch(uimg);
    title(sprintf('f=%d  neuron-masked spatial PCs, %s movie', f, wl_mov{wi}));
end
xlabel(axMc,'PC #'); ylabel(axMc,'cumulative variance'); grid(axMc,'on'); legend(axMc,'Location','se'); axis equal;
xlabel(axMv,'PC #'); ylabel(axMv,'variance explained'); grid(axMv,'on'); legend(axMv,'Location','ne'); axis equal;
sgtitle(figVE, sprintf('f=%d  movie-level PCA (neuron-masked) vs smoothing window', f));

%% (5) 20-ms subthreshold matrix  vs  RAW correlation from the FULL bleach-dr trace
% A = standard 20-ms subthreshold correlation matrix (uncorrected; cross-channel
%     diagonal, combined off-diagonal) -- reused from cell (1), Cmat_w{20 ms}.
% B = RAW (unsmoothed) correlation of the FULL bleach-dr-corrected trace, i.e. the
%     whole-movie output BVcorr_traces.mat (tr_corrected) written by
%     Analysis_20260806_BloodVesselUtkuAlgorithm.
%
% Channels: A's DIAGONAL uses the two optical channels (odd/even pixels,
% NormalizedTrace_ch{1} x {2}) so shot noise does not inflate the self term;
% A's OFF-diagonal (different ROIs) uses the combined trace. B now matches this:
% the whole-movie driver projected the corrected movie onto odd/even-pixel
% footprints (tr_corrected_ch1/ch2), so B also gets a noise-free cross-channel
% diagonal and a combined off-diagonal, via the same get_corrTensor.
%
% RAW here = NO smoothing and NO spike interpolation: only spike (peri-spike)
% frames are set to NaN, then a zero-lag correlation is taken over what remains.
sortDist = Pr.sortDist;
Cmat20   = Cmat_w{find(win_all==20)};        % A: 20-ms subthreshold matrix (uncorrected)

bvtr_file = fullfile(fpath{f}, 'BVcorr_traces.mat');   % whole-movie bleach-dr traces
if isfile(bvtr_file)
    Sbv = load(bvtr_file, 'tr_corrected', 'tr_corrected_ch1', 'tr_corrected_ch2');
    tbC = Sbv.tr_corrected(:,     sortDist)';   % combined      (nROI x Nframes), dend order
    tb1 = Sbv.tr_corrected_ch1(:, sortDist)';   % odd  pixels
    tb2 = Sbv.tr_corrected_ch2(:, sortDist)';   % even pixels
    assert(size(tbC,2)==nTime, ...
        'BVcorr_traces length (%d) ~= nTime (%d); re-run the whole-movie driver on the same frames.', ...
        size(tbC,2), nTime);

    % B: RAW -> no smoothing, no interpolation; only spike (peri-spike) frames NaN.
    % NB: use the FULL-recording spike train here -- cell (4c) rebinds `spk` to the
    % segFrames segment, which would truncate the mask to the first 5001 frames.
    spkDil = imdilate(Pr.SpikeMat(1,:) > 0, ones(1,7)) > 0;   % same +/-3 window as get_subthreshold(dilate=7)
    tb1(:, spkDil) = NaN;  tb1(isnan(trC)) = NaN;   % remove spikes + same motion/dirt mask as A
    tb2(:, spkDil) = NaN;  tb2(isnan(trC)) = NaN;
    tbC(:, spkDil) = NaN;  tbC(isnan(trC)) = NaN;
    % zero-lag correlation with cross-channel diagonal (odd x even), combined off-diagonal
    [Cb,~] = get_corrTensor(tb1, tb2, tbC, goodvec, 0, minPairs);
    Cbleach_raw = Cb(:,:,1);

    if ~exist('ord','var') || numel(ord) ~= nROI, ord = 1:nROI; end
    figure('Name','R3-20ms-vs-bleachRaw','NumberTitle','off'); clf
    tl5 = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    if ~exist('cmapCorr','var'), cmapCorr = gen_colormap([0 0.5 1; 1 1 1; 1 0 0], 256); end
    nexttile; imagesc(Cmat20(ord,ord),[-0.25 1]);      axis image; colormap(gca,cmapCorr); title('20 ms subthreshold (uncorrected)');
    nexttile; imagesc(Cbleach_raw(ord,ord),[-0.25 1]); axis image; colormap(gca,cmapCorr); title('bleach-dr, raw (no smoothing)');
    nexttile; imagesc(Cbleach_raw(ord,ord)-Cmat20(ord,ord),[-.3 .3]); axis image; colormap(gca,cmapCorr); title('bleach-raw  -  20 ms');
    cb=colorbar; cb.Layout.Tile='east';
    title(tl5, sprintf('f=%d  (leaf order)', f));

    % off-diagonal scatter + similarity
    a = Cmat20(triu1); b = Cbleach_raw(triu1); ok = ~isnan(a)&~isnan(b);
    figure('Name','R3-20ms-vs-bleachRaw-scatter','NumberTitle','off'); clf; hold on
    plot(a(ok), b(ok), '.'); plot([0 1],[0 1],'k--'); axis square; xlim([0 1]); ylim([0 1]);
    xlabel('20 ms subthreshold corr (uncorrected)'); ylabel('bleach-dr raw corr');
    title(sprintf('off-diagonal pairs (r = %.3f)', corr(a(ok),b(ok)))); grid on

    fprintf('\n[R3] 20-ms (uncorrected) vs bleach-dr-raw, off-diagonal correlation:\n');
    fprintf('   Pearson r between the two matrices = %.3f\n', corr(a(ok),b(ok)));
    fprintf('   mean 20-ms   off-diag = %.4f\n', mean(a(ok)));
    fprintf('   mean bleach  off-diag = %.4f\n', mean(b(ok)));
    fprintf('   -> if the raw bleach-dr matrix matches the 20-ms matrix, the 20-ms\n');
    fprintf('      averaging was not hiding real structure; it was suppressing blood.\n');
else
    warning(['BVcorr_traces.mat not found in %s.\n' ...
        'Run the whole-movie bleach-dr (Analysis_20260806_BloodVesselUtkuAlgorithm) first.'], fpath{f});
end

%% ################################################################
%  REVIEWER 3 -- ADDENDUM.  Everything below runs on the f=20 segment
%  (segFrames) and reuses the variables built in cells (4)/(4c):
%     mov_mc, mov_res, mov_res2, ftprintsmallest, Fall/Fodd/Feven,
%     Tall/Todd/Teven, fTall/fTodd/fTeven, Dgeo, coordS, Cxt/R0, badfr, bvBW.
%
%   (6)  two example small-footprint traces: raw / 20 ms / raw-20 ms
%   (7)  neuron image + contours of those two sub-ROIs
%   (8)  cluster the fast cross-correlation matrix -> blood regressors, and map
%        where each regressor lives in the movie (covariance / amplitude / phase)
%  (10)  PCA on the ODD/EVEN cross-covariance (shot-noise-free) vs the
%        ordinary covariance PCA  -> get_eigvector_oddeven
% ################################################################

% free the big cell-(4) intermediates that are no longer needed
% (their summaries hfvar / bloodVar / Tall.. are already computed)
clear hf bloodComp Xres

nTseg  = numel(segFrames);
gdfr   = ~badfr;                             % frames used for every fit below (no spike, no blue)
exColr = [0.85 0.33 0.10; 0.00 0.45 0.74];   % colors of the two example sub-ROIs
FvecS  = single(tovec(Fall));                % pixels x Msub, unit-sum small footprints
refImg = mean(mov_mc, 3);                    % fluorescence image for the ROI overlay

%% (6) Example small-footprint traces: raw / 20 ms average / (raw - 20 ms)   [item 1]
% Two sub-ROIs (~30 px, a few um) plotted as dF (camera counts) in three versions:
%   raw            - the full-band mov_res trace projected on the small footprint
%   20 ms average  - movmean(raw,20) = what the paper's subthreshold analysis keeps
%   raw - 20 ms    - exactly what the 20-ms average discards (the "fast" component)
% The pair is picked automatically as two sub-ROIs exSep um apart along the
% dendrite with high fast correlation and (if a vessel is present) high bvMask
% overlap -- the situation the reviewer is asking about.  Set exROI by hand to
% plot any other pair.
exROI = [46 114];                 % <- set two sub-ROI indices to override the automatic pick
exSep = [4 12];             % otherwise pick a pair this far apart (geodesic, um)

if isempty(exROI)
    ftBin  = ftprintsmallest > 0;
    npixS  = squeeze(sum(sum(ftBin, 1), 2));
    nbvS   = squeeze(sum(sum(ftBin & bvBW, 1), 2));
    bvFrac = nbvS ./ max(npixS, 1);                       % vessel overlap per sub-ROI
    pairOK = triuM & Dgeo >= exSep(1) & Dgeo <= exSep(2);
    [pi_, pj_] = find(pairOK);
    if isempty(pi_)                                       % fall back to the closest pair
        Dtmp = Dgeo;  Dtmp(~triuM) = inf;
        [~, li] = min(Dtmp(:));  [pi_, pj_] = ind2sub([Msub Msub], li);
    end
    scorePair = R0(sub2ind([Msub Msub], pi_, pj_)) + bvFrac(pi_) + bvFrac(pj_);
    [~, bp]   = max(scorePair);
    exROI     = [pi_(bp) pj_(bp)];
end
fprintf('\n[R3] example sub-ROIs: #%d and #%d  (geodesic separation %.1f um, tau=0 fast corr %.3f)\n', ...
    exROI(1), exROI(2), Dgeo(exROI(1),exROI(2)), R0(exROI(1),exROI(2)));

% dF in camera counts. Tall already comes from -mov_res, so it is a spike-positive
% residual; only a slow running median is removed so that residual bleach drift
% does not dominate the display (baseWin is a display parameter, not a filter used
% by any statistic below).
baseWin = 1001;                                 % ms, display baseline window
dFS   = Tall - movmedian(Tall, baseWin, 2, 'omitnan');   % dF (counts), spike-positive
dFS20 = movmean(dFS, Wsm, 2);                   % 20 ms moving average
dFSf  = dFS - dFS20;                            % fast component (what 20 ms removes)

% display window = longest spike-free / blue-off stretch (up to 2 s)
bwGood = bwlabel(gdfr);
runLen = accumarray(bwGood(bwGood>0)', 1);
[~, bestRun] = max(runLen);
tShow = find(bwGood == bestRun);
tShow = tShow(1:min(2000, numel(tShow)));
tms   = tShow - tShow(1);

% Row 1 overlays the raw trace (pale, thin) with its 20 ms average (full color,
% thick) so the reader sees directly what the average keeps; row 2 is the
% difference, i.e. what it throws away.
exPale = 1 - 0.45*(1 - exColr);                          % lighter tint for the raw trace
rowDat = {{dFS, dFS20}, sprintf('raw + %d ms average', Wsm); ...
          {dFSf},        sprintf('raw - %d ms', Wsm)};

get_figure('R3-exampleTraces');
tlE = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
for r = 1:2
    ax = nexttile; hold(ax,'on');
    Y   = rowDat{r,1}{1}(exROI, tShow);                  % first entry sets the scale/offset
    amp = max(max(Y,[],2) - min(Y,[],2));
    off = [0 -1.2*amp];
    for k = 1:2
        if numel(rowDat{r,1}) > 1                        % raw underneath, pale and thin
            plot(ax, tms, Y(k,:) + off(k), '-', 'Color', exPale(k,:), 'LineWidth', 0.5);
            Ys = rowDat{r,1}{2}(exROI(k), tShow);
            plot(ax, tms, Ys + off(k), '-', 'Color', exColr(k,:), 'LineWidth', 1.4);
        else
            plot(ax, tms, Y(k,:) + off(k), '-', 'Color', exColr(k,:), 'LineWidth', 0.8);
        end
    end
    ylabel(ax, rowDat{r,2}); box(ax,'off'); set(ax,'YTick',[]); xlim(ax,[0 tms(end)]);
    sb = round(amp/4, 1, 'significant');                 % \DeltaF scale bar
    y0 = off(2) + min(Y(2,:));
    plot(ax, tms(end)*[1 1]-25, y0 + [0 sb], 'k-', 'LineWidth', 2);
    text(ax, tms(end)-35, y0 + sb/2, sprintf('%g counts', sb), 'HorizontalAlignment','right','FontSize',9);
    if r == 2, xlabel(ax,'time (ms)'); end
end
title(tlE, sprintf('f=%d  small footprints #%d (orange) / #%d (blue), %.1f \\mum apart -- \\DeltaF (counts)', ...
    f, exROI(1), exROI(2), Dgeo(exROI(1),exROI(2))));

fprintf('[R3] example sub-ROI decomposition (dF variance over analysed frames):\n');
for k = 1:2
    vR = var(dFS  (exROI(k), gdfr), 0, 2, 'omitnan');
    vS = var(dFS20(exROI(k), gdfr), 0, 2, 'omitnan');
    vF = var(dFSf (exROI(k), gdfr), 0, 2, 'omitnan');
    fprintf('   sub-ROI %3d : raw=%.3g  20ms=%.3g (%.0f%%)  fast=%.3g (%.0f%%)\n', ...
        exROI(k), vR, vS, 100*vS/vR, vF, 100*vF/vR);
end
fprintf('   between the two ROIs: corr(20 ms)=%.3f   corr(fast)=%.3f\n', ...
    corr(dFS20(exROI(1),gdfr)', dFS20(exROI(2),gdfr)'), ...
    corr(dFSf (exROI(1),gdfr)', dFSf (exROI(2),gdfr)'));

% (7) Neuron image + contours of the two example small footprints          [item 2]
% Top: whole FOV, all sub-ROI outlines in grey, the two example sub-ROIs in their
% trace colors, bvMask contour in cyan.  Bottom: zoom around the pair.
climRef = prctile(refImg(:), [20 99.8]);

get_figure('R3-exampleROImap');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
for p = 1:2
    ax = nexttile; imagesc(ax, refImg, climRef); colormap(ax, gray); axis(ax,'image','off'); hold(ax,'on');
    for m = 1:Msub                                    % all sub-ROIs, faint
        B = bwboundaries(ftprintsmallest(:,:,m) > 0);
        for b = 1:numel(B), plot(ax, B{b}(:,2), B{b}(:,1), '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.4); end
    end
    contour(ax, bvBW, [0.5 0.5], 'c', 'LineWidth', 1);      % blood-vessel mask
    for k = 1:2                                       % the two example sub-ROIs
        B = bwboundaries(ftprintsmallest(:,:,exROI(k)) > 0);
        for b = 1:numel(B), plot(ax, B{b}(:,2), B{b}(:,1), '-', 'Color', exColr(k,:), 'LineWidth', 2); end
        text(ax, coordS(exROI(k),1)+4, coordS(exROI(k),2)-5, sprintf('#%d', exROI(k)), ...
            'Color', exColr(k,:), 'FontSize', 11, 'FontWeight','bold');
    end
    if p == 1
        sbpx = 20 / Pixelsize(f);                     % 20 um scale bar
        plot(ax, ncol-10-[sbpx 0], (nrow-8)*[1 1], 'w-', 'LineWidth', 3);
        text(ax, ncol-10-sbpx/2, nrow-16, '20 \mum', 'Color','w','HorizontalAlignment','center');
        title(ax, sprintf('f=%d  full FOV (grey = all %d sub-ROIs, cyan = bvMask)', f, Msub));
    else                                              % zoom on the example pair
        cxy = mean(coordS(exROI,:), 1);  hw = 60;
        xlim(ax, cxy(1) + [-hw hw]); ylim(ax, cxy(2) + [-hw/1.8 hw/1.8]);
        title(ax, 'zoom on the two example sub-ROIs');
    end
end

%% (7b) Coherence between the blood-vessel mask ROIs and every sub-ROI
% Magnitude-squared coherence
%       C_xy(f) = |P_xy(f)|^2 / ( P_xx(f) * P_yy(f) )
% is the frequency-resolved, PHASE-INVARIANT version of correlation: it asks how much
% of a sub-ROI's power at frequency f is linearly predictable from the vessel trace,
% regardless of the delay between them. That makes it the natural tool here -- a
% vessel imposes the same fluctuation on nearby dendrite with an arbitrary phase, which
% a zero-lag correlation underestimates and a max-over-lag correlation only patches up.
%
% Two things this gives that the correlation matrix cannot:
%   * WHICH BAND the vessel and the dendrite share, instead of one lumped number.
%   * The cross-spectral PHASE at that band -> a delay in ms, per sub-ROI, per vessel.
%
% Coherence is biased upward when few segments are averaged, so the null level is
% computed and drawn: with n independent segments, P(C > c) = alpha for
%       c_crit = 1 - alpha^(1/(n-1))
% Anything below that line is indistinguishable from two unrelated traces.

% ------------------------------- parameters -------------------------------
cohWin       = 128;       % samples per Welch segment. This is the critical knob: the
                          % coherence null scales as 1-alpha^(1/(n-1)) with n = T/cohWin
                          % independent segments, so a long window leaves too few
                          % segments and NOTHING can clear the null. On T=5001 with a
                          % synthetic 30 Hz vessel: win 1024 -> null 0.63, 0/20 coupled
                          % ROIs detected; win 128 -> null 0.076, 20/20 detected and
                          % 0/40 false positives. Resolution here is fs/128 = 7.8 Hz.
cohOverlap   = 0.5;       % fractional overlap between segments
cohBand      = [10 100];  % Hz, band averaged to make the spatial map
cohAlpha     = 0.05;      % significance level for the coherence null
cohDetrend   = 200;       % ms, movmedian detrend applied to every trace
cohRegNeuron = true;      % regress the neuron footprint traces out of the vessel trace
                          % first, so "vessel" is not just neuron signal bleeding in
% --------------------------------------------------------------------------

nBV = size(bvMask, 3);
fprintf('\n[R3] coherence with %d blood-vessel mask ROI(s), %d sub-ROIs\n', nBV, Msub);

% ---- vessel reference traces from mov_res (same movie the sub-ROIs come from) ----
bvW  = tovec(bvMask);  bvW = bvW ./ max(sum(bvW, 1), eps);     % unit-sum masks
bvTr = -(bvW' * tovec(mov_res));                                % nBV x T, spike-positive
if cohRegNeuron
    ftW  = tovec(ftprnt);  ftW = ftW ./ max(sum(ftW, 1), eps);
    Vtr  = -(ftW' * tovec(mov_res));                            % nROIfull x T
    Vaug = [Vtr; ones(1, nTseg)];
    bvTr = bvTr - (bvTr / Vaug) * Vaug;                         % project the neuron out
    clear ftW Vtr Vaug
end

% ---- detrend everything the same way; interpolate ONLY so Welch has continuous data --
dt = @(A) A - movmedian(A, cohDetrend, 2, 'omitnan');
Sub = dt(Tall);   Sub(:, badfr) = NaN;
Bv  = dt(bvTr);   Bv(:,  badfr) = NaN;
SubI = interpolateNaN(Sub);  SubI(isnan(SubI)) = 0;
BvI  = interpolateNaN(Bv);   BvI(isnan(BvI))  = 0;

nOv    = round(cohOverlap * cohWin);
nSegEf = floor(nTseg / cohWin);                       % independent (non-overlapping) count
cCrit  = 1 - cohAlpha^(1/max(nSegEf-1, 1));
fprintf('   %d independent segments of %d samples -> coherence null at p=%.2g is %.3f\n', ...
    nSegEf, cohWin, cohAlpha, cCrit);

% ---- coherence and cross-spectral phase, every sub-ROI x every vessel ROI ----
[~, fCoh] = mscohere(SubI(1,:), BvI(1,:), hann(cohWin), nOv, [], fs);
fCoh = fCoh(:);  nF = numel(fCoh);
Coh   = nan(Msub, nF, nBV);
PhCoh = nan(Msub, nF, nBV);
for v = 1:nBV
    for m = 1:Msub
        Coh(m,:,v) = mscohere(SubI(m,:), BvI(v,:), hann(cohWin), nOv, [], fs);
        PhCoh(m,:,v) = angle(cpsd(SubI(m,:), BvI(v,:), hann(cohWin), nOv, [], fs));
    end
end

inBand  = fCoh >= cohBand(1) & fCoh <= cohBand(2);
cohMapV = squeeze(mean(Coh(:, inBand, :), 2));        % Msub x nBV, band-averaged
if nBV == 1, cohMapV = cohMapV(:); end

% delay per sub-ROI: fit phase = 2*pi*f*tau across the band, weighted by coherence,
% using only the bins that clear the null (a phase estimate at C ~ 0 is meaningless)
delMs = nan(Msub, nBV);
for v = 1:nBV
    for m = 1:Msub
        ok = inBand & squeeze(Coh(m,:,v))' > cCrit;
        if nnz(ok) < 4, continue; end
        ph = unwrap(squeeze(PhCoh(m,ok,v)));
        pf = polyfit(fCoh(ok), ph(:), 1);
        % cpsd(x,y) has phase -2*pi*f*d when x lags y by d, so the delay is MINUS the
        % slope: positive delMs = this sub-ROI lags the vessel.
        delMs(m,v) = -1000 * pf(1) / (2*pi);
    end
end

% ---- is a sub-ROI ON a vessel?  (used only to split the curves below) ----
ftBinC  = ftprintsmallest > 0;
onVesss = squeeze(sum(sum(ftBinC & (max(bvMask,[],3) > 0), 1), 2)) ./ ...
          max(squeeze(sum(sum(ftBinC,1),2)), 1) > 0.2;
fprintf('   %d of %d sub-ROIs overlap a vessel mask by >20%% of their area\n', ...
    nnz(onVesss), Msub);

% ---- (1) coherence spectra: on-vessel vs off-vessel sub-ROIs ----
get_figure('R3-bvCoherence-spectra');
tlH = tiledlayout(1, nBV, 'TileSpacing','compact','Padding','compact');
for v = 1:nBV
    axh = nexttile; hold(axh,'on');
    grp = {onVesss, ~onVesss};  gcol = [0.85 0.2 0.2; 0.3 0.3 0.3];
    gnm = {'sub-ROIs on a vessel', 'sub-ROIs off vessel'};
    for gg = 1:2
        if ~any(grp{gg}), continue; end
        Cg = Coh(grp{gg}, :, v);
        mu = mean(Cg, 1, 'omitnan');
        lo = prctile(Cg, 25, 1);  hi = prctile(Cg, 75, 1);
        patch(axh, [fCoh; flipud(fCoh)], [lo(:); flipud(hi(:))], gcol(gg,:), ...
            'EdgeColor','none', 'FaceAlpha', 0.2, 'HandleVisibility','off');
        plot(axh, fCoh, mu, '-', 'Color', gcol(gg,:), 'LineWidth', 1.6, 'DisplayName', gnm{gg});
    end
    yline(axh, cCrit, 'k--', sprintf('p=%.2g null', cohAlpha), 'DisplayName','null');
    xline(axh, cohBand(1), ':'); xline(axh, cohBand(2), ':');
    set(axh,'XScale','log'); xlim(axh, [max(fCoh(2),1) fs/2]); ylim(axh, [0 1]);
    xlabel(axh,'frequency (Hz)'); ylabel(axh,'magnitude-squared coherence');
    grid(axh,'on'); box(axh,'off'); if v == 1, legend(axh,'Location','northeast'); end
    title(axh, sprintf('vessel ROI %d', v));
end
title(tlH, sprintf(['f=%d  coherence with the blood-vessel mask (shaded = IQR over ' ...
    'sub-ROIs, dashed = null level)'], f));

% ---- (2) coherence and delay painted on the morphology ----
ftLc7 = double(ftprintsmallest);  denc7 = sum(ftLc7, 3);
[~, labImg7] = max(ftLc7, [], 3);  inMask7 = denc7 > eps;  labImg7(~inMask7) = 0;
clear ftLc7 denc7

get_figure('R3-bvCoherence-map');
tlG = tiledlayout(ceil(nBV/2), 4, 'TileSpacing','compact','Padding','compact');
dLim = max(1, prctile(abs(delMs(:)), 95));
for v = 1:nBV
    for pnl = 1:2
        if pnl == 1
            val = cohMapV(:,v);  cm = turbo;      cl = [0 max(prctile(val,99), cCrit*1.5)];
            ttl = sprintf('vessel %d: coherence, %g-%g Hz', v, cohBand(1), cohBand(2));
        else
            val = delMs(:,v);    cm = cmapCorr;   cl = dLim*[-1 1];
            ttl = sprintf('vessel %d: delay from phase slope (ms)', v);
        end
        axg = nexttile;
        vmap = nan(nrow, ncol);  vmap(inMask7) = val(labImg7(inMask7));
        hG = imagesc(axg, vmap, cl);  set(hG, 'AlphaData', ~isnan(vmap));
        axis(axg,'image','off'); colormap(axg, cm); colorbar(axg); hold(axg,'on');
        contour(axg, bvMask(:,:,v) > 0.2*max(tovec(bvMask(:,:,v))), [0.5 0.5], ...
            'w', 'LineWidth', 1);
        title(axg, ttl);
    end
end
title(tlG, sprintf(['f=%d  per-sub-ROI coherence with each vessel mask (white = that ' ...
    'vessel). Coherence null = %.3f'], f, cCrit));

fprintf('   band-averaged coherence (%g-%g Hz), median [on vessel | off vessel]:\n', ...
    cohBand(1), cohBand(2));
for v = 1:nBV
    fprintf('      vessel %d : %.3f | %.3f   (null %.3f, %.0f%% of sub-ROIs above null)\n', ...
        v, median(cohMapV(onVesss,v),'omitnan'), median(cohMapV(~onVesss,v),'omitnan'), ...
        cCrit, 100*mean(cohMapV(:,v) > cCrit));
    dv = delMs(:,v);  dv = dv(isfinite(dv));
    if ~isempty(dv)
        fprintf('                 delay median %+.1f ms, IQR %+.1f..%+.1f ms (n=%d ROIs)\n', ...
            median(dv), prctile(dv,25), prctile(dv,75), numel(dv));
    end
end


%% (7c0) Sweep the spatial passband before committing to one: var(S)/var(mov_res)
% Which spatial band actually isolates the vessel signal is an empirical question, so
% scan a grid of (lo, hi) and score each one. Three pixel populations are compared:
%       vessel   = bvMask
%       dendrite = neuron footprint but NOT vessel
%       background = neither
% and the number that matters is the SELECTIVITY
%       [var(S) on vessel / var(S) on dendrite]  /  [same ratio for mov_res]
% i.e. how much the band enriches vessel signal RELATIVE to dendritic signal, over and
% above what the unfiltered movie already gives. A band that just attenuates everything
% equally scores 1 no matter how much total variance it removes.
%
% The FFT of the movie is taken ONCE and reused for every (lo,hi), so the cost is one
% ifft2 per grid point. Only sbSwFrames frames are used -- a variance estimate converges
% long before 5000 frames.

% ------------------------------- parameters -------------------------------
sbSwLo     = [2 3 4 5 6 8];        % um, fine edges to try
sbSwHi     = [6 8 10 14 20 30];    % um, coarse edges to try
sbSwOrder  = 4;                    % Butterworth order (same as (7c))
sbSwFrames = 1000;                 % frames used for the variance estimate
voltRegress = true;                % regress PCresult{f}.normTraces out of mov_res first
% --------------------------------------------------------------------------
% ---- remove the neuron's VOLTAGE signal before anything spatial ------------------
% The dendritic Vm is spatially smooth over the length constant, so a spatial bandpass
% suppresses it -- but only partly, because the neurite is a thin bright structure and
% the filter responds to its EDGES. Projecting the measured voltage traces out of the
% movie first removes that signal at its source, wherever it sits spatially, so what is
% left is much closer to "blood only".
%   normTraces is nROI x nTimeFull for the whole recording; only segFrames belong to
%   mov_res. Rows that are all-NaN or constant on this segment are dropped, and a
%   constant term is included so SeeResiduals also takes out the mean image.
% Set voltRegress = false to work on the raw mov_res instead.
if ~exist('movNvFlag','var'), movNvFlag = []; end
if ~exist('mov_nv','var') || size(mov_nv,3) ~= nTseg || ~isequal(voltRegress, movNvFlag)
    if voltRegress
        Vtr = PCresult{f}.normTraces;
        assert(size(Vtr,2) >= segFrames(end), ...
            'normTraces has %d frames, segFrames needs %d.', size(Vtr,2), segFrames(end));
        Vtr = Vtr(:, segFrames);
        keepV = all(isfinite(Vtr), 2) & std(Vtr, 0, 2) > 0;
        Vtr = Vtr(keepV, :);
        Vtr = (Vtr - mean(Vtr,2)) ./ std(Vtr,0,2);          % scale-free, better conditioned
        fprintf('\n[R3] regressing %d of %d voltage traces out of mov_res ...\n', ...
            nnz(keepV), numel(keepV));
        mov_nv = single(SeeResiduals(mov_res, Vtr));
        vr = var(double(mov_nv), 0, 3);  v0 = var(mov_res, 0, 3);
        fprintf(['     variance left after voltage regression: %.1f%% overall, ' ...
            '%.1f%% on the neuron, %.1f%% on bvMask\n'], ...
            100*mean(vr(:))/mean(v0(:)), ...
            100*mean(vr(neuronMap>0))/mean(v0(neuronMap>0)), ...
            100*mean(vr(bvBW))/mean(v0(bvBW)));
        clear Vtr vr v0
    else
        mov_nv = single(mov_res);
    end
    movNvFlag = voltRegress;
end


pxUm = Pixelsize(f);
frSw = 1:min(sbSwFrames, nTseg);

% pixel populations
pxVess = bvBW & true(nrow,ncol);
pxDend = (neuronMap > 0) & ~pxVess;
pxBg   = ~pxVess & ~pxDend;
fprintf('\n[R3] passband sweep on %d frames: %d vessel / %d dendrite / %d background pixels\n', ...
    numel(frSw), nnz(pxVess), nnz(pxDend), nnz(pxBg));
assert(nnz(pxVess) > 20 && nnz(pxDend) > 20, 'Need both vessel and dendrite pixels.');

% reference variance of the unfiltered movie on the same frames
varRef  = var(double(mov_nv(:,:,frSw)), 0, 3);
refRatio = mean(varRef(pxVess)) / mean(varRef(pxDend));
fprintf('     input vessel/dendrite variance ratio = %.3f (the baseline to beat)\n', refRatio);

% radial frequency grid, and the movie's spectrum computed once
[fxG, fyG] = meshgrid(((0:ncol-1)-floor(ncol/2))/ncol, ((0:nrow-1)-floor(nrow/2))/nrow);
frG  = sqrt(fxG.^2 + fyG.^2) / pxUm;                  % cycles per um
Fmov = fft2(mov_nv(:,:,frSw));            % mov_nv is already single

nLo = numel(sbSwLo);  nHi = numel(sbSwHi);
swKeep = nan(nLo, nHi);   % fraction of total variance kept
swVsD  = nan(nLo, nHi);   % var(S) vessel / dendrite
swSel  = nan(nLo, nHi);   % that ratio, relative to mov_res  <- the score
swBand = nan(nLo, nHi, 2);% achieved half-power band
for a = 1:nLo
    for b = 1:nHi
        if sbSwHi(b) <= sbSwLo(a)*1.2, continue; end          % need a real band
        Hb = 1 ./ (1 + (frG*sbSwLo(a)).^(2*sbSwOrder)) .* ...
             1 ./ (1 + (1./max(frG*sbSwHi(b),eps)).^(2*sbSwOrder));
        Ss = real(ifft2(Fmov .* ifftshift(Hb)));
        vs = var(double(Ss), 0, 3);
        swKeep(a,b) = mean(vs(:)) / mean(varRef(:));
        swVsD(a,b)  = mean(vs(pxVess)) / mean(vs(pxDend));
        swSel(a,b)  = swVsD(a,b) / refRatio;
        % achieved half-power band for this grid point
        eE = linspace(0, max(frG(:)), 200);  eC = (eE(1:end-1)+eE(2:end))/2;
        Hr = arrayfun(@(q) mean(Hb(frG>=eE(q) & frG<eE(q+1)),'omitnan'), 1:numel(eC));
        ok = isfinite(Hr) & eC > 0;  pu = 1./eC(ok);  Ho = Hr(ok);
        hh = Ho >= max(Ho)/2;
        swBand(a,b,:) = [min(pu(hh)) max(pu(hh))];
        clear Ss vs
    end
end
clear Fmov

% ---- report ----
[~, iBest] = max(swSel(:));
[aB, bB]   = ind2sub([nLo nHi], iBest);
fprintf('     best selectivity %.2fx at lo=%g hi=%g um (achieved %.1f-%.1f um), keeps %.1f%% of variance\n', ...
    swSel(aB,bB), sbSwLo(aB), sbSwHi(bB), swBand(aB,bB,1), swBand(aB,bB,2), 100*swKeep(aB,bB));
fprintf('     selectivity grid (rows = lo um, cols = hi um):\n         ');
fprintf('%7g', sbSwHi); fprintf('\n');
for a = 1:nLo
    fprintf('   %5g ', sbSwLo(a));  fprintf('%7.2f', swSel(a,:));  fprintf('\n');
end
if swSel(aB,bB) < 1.1
    warning(['No passband enriches vessel over dendrite by more than 10%%. ' ...
             'The 5-10 um band may simply not separate the two here.']);
end

% ---- heatmaps ----
get_figure('R3-spatialBP-sweep');
tlSw = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
swSets = {swSel,  'selectivity: (vessel/dendrite) / mov_res', [], cmapCorr; ...
          swVsD,  'var(S) vessel / dendrite',                 [], turbo; ...
          100*swKeep, '% of mov\_res variance kept',          [], turbo};
for q = 1:3
    axq = nexttile;  M = swSets{q,1};
    if q == 1
        cl = [2-max(M(:),[],'omitnan') max(M(:),[],'omitnan')];       % 1 = white
        if ~all(isfinite(cl)) || cl(2) <= cl(1), cl = [0 2]; end
    else
        cl = [0 max(M(:),[],'omitnan')];
    end
    hI = imagesc(axq, M, cl);  set(hI, 'AlphaData', ~isnan(M));
    colormap(axq, swSets{q,4}); colorbar(axq);
    set(axq, 'XTick', 1:nHi, 'XTickLabel', sbSwHi, 'YTick', 1:nLo, 'YTickLabel', sbSwLo);
    xlabel(axq, 'hi (\mum)'); ylabel(axq, 'lo (\mum)');
    title(axq, swSets{q,2});  axis(axq,'square');
    hold(axq,'on'); plot(axq, bB, aB, 'ks', 'MarkerSize', 14, 'LineWidth', 2);
    for a = 1:nLo
        for b = 1:nHi
            if isfinite(M(a,b))
                text(axq, b, a, sprintf('%.2f', M(a,b)), 'HorizontalAlignment','center', ...
                    'FontSize', 7, 'Color', [0 0 0]);
            end
        end
    end
end
title(tlSw, sprintf(['f=%d  spatial passband sweep (black square = best selectivity). ' ...
    'Panel 1 is the one to read: >1 means the band favours vessel over dendrite'], f));


%% (7c) Spatial-bandpass movie S: keep only structure 5-10 um wide
% The three things in mov_res live at very different SPATIAL scales:
%     shot noise      - independent per pixel            (< 1 um)
%     blood           - coherent over ~10 um, gone by 30 um  (measured in cell (4c))
%     dendritic Vm    - coherent over the length constant  (>> 50 um)
% So a purely SPATIAL bandpass is a way to isolate the vascular component that uses no
% temporal assumption at all -- no filter band, no lag search, no clustering. Whatever
% survives a 5-10 um passband is, by construction, the part of the movie that varies on
% the vessel's spatial scale.
%
% Implementation: a radial band applied frame by frame, either in the 2D Fourier domain
% (Butterworth) or as a Difference-of-Gaussians kernel. The passband that actually comes
% out is MEASURED from the filter's radial transfer function, not assumed -- see the
% note under the parameters for why that matters.
%
% Caveat worth keeping in mind: the dendrite is a thin bright structure on a dark
% background, so the DoG also responds to its EDGES. Read the maps inside the neuron
% mask, and treat a rim of signal around the dendrite as an artefact of the filter.

% ------------------------------- parameters -------------------------------
sbLoUm    = 6;      % um, fine end of the spatial passband
sbHiUm    = 30;     % um, coarse end
sbMethod  = 'dog';  % 'fourier' = radial Butterworth band in the 2D FFT (sharp, hits
                        %             the requested band; mild ringing)
                        % 'dog'     = Difference of Gaussians (no ringing, but see below)
sbOrder   = 4;      % Butterworth order for 'fourier'
sbChunk   = 200;    % frames filtered per pass (memory)
sbMovFrm  = 200;    % frames written to the comparison mp4
voltRegress = true;         % regress PCresult{f}.normTraces out of mov_res first
% -- tracked-movie overlay (the [input | S | S+tracks] mp4 at the end) -----------
sbMovRate   = 20;           % playback fps
sbTailFrm   = 8;            % frames a trajectory point stays painted (the tail)
sbMarkPx    = 1;            % half-width of the trajectory marker, px (1 -> 3x3)
sbColBright = [1 0 0];      % RED  : bright-polarity tracklets
sbColDark   = [0 0.45 1];   % BLUE : dark-polarity tracklets
% --------------------------------------------------------------------------
% ---- remove the neuron's VOLTAGE signal before anything spatial ------------------
% The dendritic Vm is spatially smooth over the length constant, so a spatial bandpass
% suppresses it -- but only partly, because the neurite is a thin bright structure and
% the filter responds to its EDGES. Projecting the measured voltage traces out of the
% movie first removes that signal at its source, wherever it sits spatially, so what is
% left is much closer to "blood only".
%   normTraces is nROI x nTimeFull for the whole recording; only segFrames belong to
%   mov_res. Rows that are all-NaN or constant on this segment are dropped, and a
%   constant term is included so SeeResiduals also takes out the mean image.
% Set voltRegress = false to work on the raw mov_res instead.
if ~exist('movNvFlag','var'), movNvFlag = []; end
if ~exist('mov_nv','var') || size(mov_nv,3) ~= nTseg || ~isequal(voltRegress, movNvFlag)
    if voltRegress
        Vtr = PCresult{f}.normTraces;
        assert(size(Vtr,2) >= segFrames(end), ...
            'normTraces has %d frames, segFrames needs %d.', size(Vtr,2), segFrames(end));
        Vtr = Vtr(:, segFrames);
        keepV = all(isfinite(Vtr), 2) & std(Vtr, 0, 2) > 0;
        Vtr = Vtr(keepV, :);
        Vtr = (Vtr - mean(Vtr,2)) ./ std(Vtr,0,2);          % scale-free, better conditioned
        fprintf('\n[R3] regressing %d of %d voltage traces out of mov_res ...\n', ...
            nnz(keepV), numel(keepV));
        mov_nv = single(SeeResiduals(mov_res, Vtr));
        vr = var(double(mov_nv), 0, 3);  v0 = var(mov_res, 0, 3);
        fprintf(['     variance left after voltage regression: %.1f%% overall, ' ...
            '%.1f%% on the neuron, %.1f%% on bvMask\n'], ...
            100*mean(vr(:))/mean(v0(:)), ...
            100*mean(vr(neuronMap>0))/mean(v0(neuronMap>0)), ...
            100*mean(vr(bvBW))/mean(v0(bvBW)));
        clear Vtr vr v0
    else
        mov_nv = single(mov_res);
    end
    movNvFlag = voltRegress;
end

% WHY NOT A PLAIN DoG BY DEFAULT: a Difference of Gaussians cannot make a band this
% narrow. Measured on this geometry, the narrowest half-power band a DoG reaches is a
% ratio of ~3.4:1 whatever sigmas are used, so a request for 5-10 um (ratio 2) comes
% out as roughly 4.4-14.7 um -- it leaks in structure well above the vessel scale.
% A 4th-order radial Butterworth band gives 4.9-10.1 um, i.e. the band actually asked
% for, at the cost of a -0.24 side lobe in its impulse response. Both are implemented;
% the achieved band is MEASURED below either way rather than assumed.

pxUm = Pixelsize(f);
fLoC = 1/sbHiUm;  fHiC = 1/sbLoUm;                  % passband in cycles per um

switch sbMethod
    case 'fourier'
        [fxG, fyG] = meshgrid(((0:ncol-1)-floor(ncol/2))/ncol, ...
                              ((0:nrow-1)-floor(nrow/2))/nrow);
        frG   = sqrt(fxG.^2 + fyG.^2) / pxUm;       % cycles per um
        Hband = 1 ./ (1 + (frG/fHiC).^(2*sbOrder)) .* ...
                1 ./ (1 + (fLoC./max(frG,eps)).^(2*sbOrder));
        Hmask = ifftshift(Hband);
        frMeas = frG;  Hmeas = Hband;  kDoG = [];
    case 'dog'
        % calibrate the two sigmas numerically so the MEASURED band is centred on the
        % requested one (the analytic sigma = d/2 rule is wrong by ~2.7x for a Gaussian,
        % whose half-power period is ~10.7*sigma, not 2*sigma)
        rat = 1.2;  a = 0.2;  b = 20;  tgt = sqrt(sbLoUm*sbHiUm);
        for it = 1:40
            s1 = sqrt(a*b);
            kk = fspecial('gaussian', 2*ceil(4*max(s1*rat,1))+1, s1) - ...
                 fspecial('gaussian', 2*ceil(4*max(s1*rat,1))+1, s1*rat);
            KK = abs(fftshift(fft2(kk, 256, 256)));
            [gx, gy] = meshgrid((-128:127)/256);  gr = sqrt(gx.^2+gy.^2);
            eE = linspace(0,0.5,200);  eC = (eE(1:end-1)+eE(2:end))/2;
            Kr = arrayfun(@(q) mean(KK(gr>=eE(q) & gr<eE(q+1)),'omitnan'), 1:numel(eC));
            [~, ipk] = max(Kr);  pk = 1/(eC(ipk)/pxUm);
            if pk > tgt, b = s1; else, a = s1; end
        end
        sLo = s1;  sHi = s1*rat;
        kSz = 2*ceil(4*sHi)+1;
        kDoG = fspecial('gaussian', kSz, sLo) - fspecial('gaussian', kSz, sHi);
        [fxG, fyG] = meshgrid((-128:127)/256);
        frMeas = sqrt(fxG.^2+fyG.^2) / pxUm;
        Hmeas  = abs(fftshift(fft2(kDoG, 256, 256)));
        Hmask  = [];
        fprintf('\n[R3] DoG calibrated: sigma %.2f / %.2f px\n', sLo, sHi);
    otherwise
        error('sbMethod must be ''fourier'' or ''dog''.');
end

% ---- what band does this filter ACTUALLY pass? (radial transfer function) ----
eE = linspace(0, max(frMeas(:)), 300);  eC = (eE(1:end-1)+eE(2:end))/2;
Hrad = arrayfun(@(q) mean(Hmeas(frMeas>=eE(q) & frMeas<eE(q+1)), 'omitnan'), 1:numel(eC));
okH  = isfinite(Hrad) & eC > 0;
perUm = 1 ./ eC(okH);  Hok = Hrad(okH);
[hPk, ipk] = max(Hok);   half = Hok >= hPk/2;
bandUm = [min(perUm(half)) max(perUm(half))];
fprintf('\n[R3] spatial bandpass (%s, %.2f um/px): requested %g-%g um\n', ...
    sbMethod, pxUm, sbLoUm, sbHiUm);
fprintf('     ACHIEVED half-power band %.1f-%.1f um (peak %.1f um, ratio %.2f vs %.2f asked)\n', ...
    bandUm(1), bandUm(2), perUm(ipk), bandUm(2)/bandUm(1), sbHiUm/sbLoUm);
if bandUm(1) < 0.7*sbLoUm || bandUm(2) > 1.4*sbHiUm
    warning(['The achieved band is much wider than requested -- interpret S as the ' ...
             '%.1f-%.1f um band, not %g-%g um.'], bandUm(1), bandUm(2), sbLoUm, sbHiUm);
end

get_figure('R3-spatialBP-transfer');
tlT = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;
if strcmp(sbMethod,'dog')
    imagesc(kDoG); axis image; colormap(gca, cmapCorr); caxis(max(abs(kDoG(:)))*[-1 1]);
    title(sprintf('DoG kernel (\\sigma %.2f / %.2f px)', sLo, sHi));
else
    imagesc(Hband); axis image; colormap(gca, turbo);
    title(sprintf('2D Butterworth band (order %d)', sbOrder));
end
colorbar;
nexttile; semilogx(perUm, Hok/hPk, '-', 'LineWidth', 1.5); hold on
xline(sbLoUm, 'r--', 'requested lo'); xline(sbHiUm, 'r--', 'requested hi');
yline(0.5, 'k:', 'half power'); grid on; box off;
set(gca,'XDir','reverse'); xlim([1 200]);
xlabel('spatial period (\mum)'); ylabel('normalised gain');
title(sprintf('achieved %.1f-%.1f \\mum', bandUm(1), bandUm(2)));
title(tlT, sprintf('f=%d  spatial bandpass: filter and its measured radial transfer', f));

% ---- build S ----
S = zeros(nrow, ncol, nTseg, 'single');
for c0 = 1:sbChunk:nTseg
    ci = c0:min(c0+sbChunk-1, nTseg);
    if strcmp(sbMethod,'fourier')
        S(:,:,ci) = single(real(ifft2(fft2(mov_nv(:,:,ci)) .* Hmask)));
    else
        S(:,:,ci) = imfilter(mov_nv(:,:,ci), kDoG, 'replicate');
    end
end
fprintf('   S built: %s single (%.1f GB)\n', mat2str(size(S)), numel(S)*4/2^30);

% ---- does S concentrate on the vessels? ----
varRes = var(double(mov_nv), 0, 3);
varS   = var(double(S), 0, 3);
bvpixS = bvBW(:) > 0;
neurS  = neuronMap(:) > 0;
fprintf('   variance in bvMask / outside:  input %.2f  ->  S %.2f\n', ...
    mean(varRes(bvpixS))/mean(varRes(~bvpixS)), mean(varS(bvpixS))/mean(varS(~bvpixS)));
fprintf('   variance on neuron / off:      input %.2f  ->  S %.2f\n', ...
    mean(varRes(neurS))/mean(varRes(~neurS)), mean(varS(neurS))/mean(varS(~neurS)));
fprintf('   S keeps %.1f%% of the input movie variance\n', 100*mean(varS(:))/mean(varRes(:)));

get_figure('R3-spatialBP-var');
tlS = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
vSets = {varRes, 'var(input movie)'; varS, sprintf('var(S), %g-%g \\mum', sbLoUm, sbHiUm); ...
         varS./max(varRes,eps), 'var(S) / var(input)'};
for q = 1:3
    axv = nexttile;  vv = vSets{q,1};
    imagesc(axv, vv, [0 prctile(vv(:), 99)]); axis(axv,'image','off'); hold(axv,'on');
    colormap(axv, turbo); colorbar(axv);
    contour(axv, bvBW, [0.5 0.5], 'w', 'LineWidth', 0.8);
    title(axv, vSets{q,2});
end
title(tlS, sprintf(['f=%d  where the %.1f-%.1f \\mum spatial band lives (white = bvMask). ' ...
    'If it is vascular, panel 3 is bright on vessels and dark on the dendrite'], ...
    f, bandUm(1), bandUm(2)));


ptDir = fullfile(fileparts(which('divide_ROIs')), 'Image Processing', 'Particle Tracking');
if isempty(which('pkfnd')) && isfolder(ptDir), addpath(ptDir); end
assert(~isempty(which('pkfnd')) && ~isempty(which('track')), ...
    'pkfnd/track not found; expected them in %s', ptDir);

%% (7d) ---- RBC particle tracking in S: detect per frame, link into tracklets -----------
% A red blood cell moving through a capillary is a small spot that is displaced frame to
% frame. S already isolates that spatial scale, so the pipeline is the classic
% Crocker-Grier one, using the repo's implementation:
%       pkfnd  -> local maxima above threshold, separated by trkFeatPx
%       cntrd  -> sub-pixel centroid + brightness + radius of gyration
%       track  -> link detections across frames into trajectories
% RBCs can appear either BRIGHT or DARK depending on whether they scatter or absorb, so
% both polarities are detected -- dark ones by running the same detector on -S.
%
% BRIGHT AND DARK ARE DETECTED AND LINKED SEPARATELY, never linked to each other. That
% is deliberate: a physical object keeps its polarity as it moves, so a bright-to-dark
% link would almost always be a spurious jump between an RBC and the plasma gap next to
% it. The cost is that a cell which flips apparent polarity (e.g. drifting through
% focus) breaks into two tracklets rather than one. trkPol records which population each
% tracklet came from, and the two are reported and drawn separately.
%
% The output that matters is the FLOW MAP: every tracklet contributes its displacement
% to a spatial bin, and the binned mean velocity is drawn as a quiver field. Real flow
% gives long tracklets whose directions agree within a vessel; noise gives short
% tracklets pointing every which way, which is exactly what the shuffled control shows.
% =========================== PARAMETERS (all of them) ============================
% THE FULL SET OF THRESHOLDS, in the order they are applied:
%
%   DETECT   trkThrSD    a pixel must exceed this many SD of S to be a peak
%            trkFeatPx   peaks closer than this are merged (pkfnd separation)
%   LINK     trkMaxDisp  a particle may move at most this far between frames
%            trkMem      frames it may vanish and still be linked (0 = none)
%            trkMinLen   a linked trajectory shorter than this is discarded here
%   KEEP     trkMinLenKeep  the STRONGEST filter: minimum trajectory length. A long
%                       trajectory is hard to build by chance and gives a far better
%                       velocity (the velocity uses only the two endpoints), so this is
%                       the threshold to tune first.
%            trkMinAmpPrc   drop the weakest tracklets by |amplitude| percentile
%            trkMinPersist  travel consistency, mean cos(angle between steps).
%                       0 = random walk regardless of length, 1 = straight travel.
%                       Secondary to length: it is left permissive by default and is
%                       mainly there to be READ, not to do the filtering.
trkFrames   = 1:min(2000, nTseg);  % frames tracked (the whole segment is slow)
trkThrSD    = 4;                   % detection threshold, in SD of S over trkFrames
trkFeatPx   = 2*round(sbLoUm/pxUm/2)+1;  % odd feature size, from the passband
trkMaxDisp  = 5;                   % px, largest displacement allowed between frames
trkMinLen   = 5;                   % frames, shortest tracklet kept at the linking stage
trkMem      = 0;                   % frames a particle may vanish and still be linked.
                                   % NOTE: the Crocker-Grier track() in this repo aborts
                                   % ("You are not going any further") for mem > 0 on
                                   % dense data, so 0 is the safe default; the code
                                   % falls back to 0 automatically if a larger value
                                   % fails.
trkPolarity = 'both';              % 'bright' | 'dark' | 'both'
trkBinPx    = 12;                  % px, spatial bin for the flow map
% -- keep filter (applied after linking; see the note above) ---------------------
trkMinLenKeep = 4;    % frames. THE PRIMARY THRESHOLD -- raise this first.
trkMinAmpPrc  = 0;    % drop tracklets below this percentile of |amplitude|
trkMinPersist = -1;     % mean cos(step angle); -1 keeps everything, ~0.3 is strict
% =================================================================================

vSm   = imgaussfilt(varS, 1);          % smoothed variance map, used as the backdrop
Strk  = S(:,:,trkFrames);
sdS   = std(double(Strk(:)));
thrTrk = trkThrSD * sdS;
switch trkPolarity
    case 'bright', pols = 1;    polNm = {'bright'};
    case 'dark',   pols = -1;   polNm = {'dark'};
    case 'both',   pols = [1 -1]; polNm = {'bright','dark'};
    otherwise, error('trkPolarity must be bright | dark | both');
end
fprintf(['\n[R3] RBC tracking on %d frames: threshold %.3g (=%gSD), feature %d px, ' ...
    'max displacement %d px\n'], numel(trkFrames), thrTrk, trkThrSD, trkFeatPx, trkMaxDisp);

trkAll = {};   trkPol = [];
for pp = 1:numel(pols)
    posT = cell(1, numel(trkFrames));
    nDet = 0;
    for t = 1:numel(trkFrames)
        im = pols(pp) * double(Strk(:,:,t));
        pk = pkfnd(im, thrTrk, trkFeatPx);
        if isempty(pk), posT{t} = zeros(0,4); continue; end
        ct = cntrd(im, pk, trkFeatPx+2);
        if isempty(ct), posT{t} = zeros(0,4); continue; end
        posT{t} = ct(:,1:4);   % [x  y  brightness  Rg^2] -- cntrd returns all four
        nDet = nDet + size(ct,1);
    end
    fprintf('   %-6s : %d detections (%.1f per frame)\n', polNm{pp}, nDet, nDet/numel(trkFrames));

    % link -> trajectories. track() takes [x y <extra...> frame] with param.dim = 2, and
    % carries the extra columns through, so brightness and Rg^2 stay attached to every
    % detection. Columns of each tracklet: [x y brightness Rg2 frame].
    ppos = [];
    for t = 1:numel(posT)
        if ~isempty(posT{t}), ppos = [ppos; posT{t}, repmat(t, size(posT{t},1), 1)]; end %#ok<AGROW>
    end
    if size(ppos,1) < 10, warning('Too few %s detections to link.', polNm{pp}); continue; end
    lk = [];
    for memTry = unique([trkMem 0], 'stable')
        try
            lk = track(ppos, trkMaxDisp, ...
                struct('mem',memTry,'good',0,'dim',2,'quiet',1));
            if memTry ~= trkMem
                fprintf('   (track() failed at mem=%d, succeeded at mem=%d)\n', trkMem, memTry);
            end
            break
        catch ME
            fprintf('   track() failed at mem=%d: %s\n', memTry, ME.message);
        end
    end
    if isempty(lk)
        warning(['track() could not link the %s detections. Reduce trkMaxDisp, raise ' ...
                 'trkThrSD so there are fewer detections, or shorten trkFrames.'], polNm{pp});
        continue
    end
    for idp = unique(lk(:,6))'
        r = lk(lk(:,6)==idp, 1:5);   % [x y bright Rg2 frame]
        if size(r,1) >= trkMinLen, trkAll{end+1} = r; trkPol(end+1) = pols(pp); end %#ok<AGROW>
    end
end
nTrk = numel(trkAll);
fprintf('   linked %d tracklets of >= %d frames  (%d bright, %d dark; linked separately)\n', ...
    nTrk, trkMinLen, nnz(trkPol > 0), nnz(trkPol < 0));

if nTrk > 0
    % per-tracklet displacement, speed and direction (total displacement / duration)
    trkVel = nan(nTrk, 2);  trkSpd = nan(nTrk,1);  trkLen = nan(nTrk,1);  trkCen = nan(nTrk,2);
    trkAmp = nan(nTrk,1);   trkSig = nan(nTrk,1);   % brightness and Gaussian width (um)
    for q = 1:nTrk
        r = trkAll{q};
        trkAmp(q) = median(r(:,3));
        trkSig(q) = median(sqrt(max(r(:,4),0)/2)) * pxUm;   % Rg^2 = 2*sigma^2
        dxy = (r(end,1:2) - r(1,1:2)) * pxUm;          % um
        dtS = (r(end,5) - r(1,5)) / fs;                % s
        trkVel(q,:) = dxy / max(dtS, eps);             % um/s
        trkSpd(q)   = norm(trkVel(q,:));
        trkLen(q)   = size(r,1);
        trkCen(q,:) = mean(r(:,1:2), 1);
    end
    fprintf(['   speed: median %.0f um/s (IQR %.0f-%.0f), tracklet length ' ...
        'median %d frames\n'], ...
        median(trkSpd), prctile(trkSpd,25), prctile(trkSpd,75), median(trkLen));
    for pv = [1 -1]
        sel = trkPol(:) == pv;
        if ~any(sel), continue; end
        if pv > 0, nmP = 'bright'; else, nmP = 'dark  '; end
        fprintf(['     %s: n=%3d, speed %.0f um/s, direction %+.0f deg, ' ...
                 'width sigma %.2f um (FWHM %.2f um)\n'], ...
            nmP, nnz(sel), median(trkSpd(sel)), ...
            median(mod(atan2d(trkVel(sel,2), trkVel(sel,1)), 180)), ...
            median(trkSig(sel),'omitnan'), 2.355*median(trkSig(sel),'omitnan'));
    end
    fprintf(['   particle width from cntrd''s radius of gyration: median sigma %.2f um ' ...
             '(IQR %.2f-%.2f), i.e. FWHM %.2f um\n'], ...
        median(trkSig,'omitnan'), prctile(trkSig,25), prctile(trkSig,75), ...
        2.355*median(trkSig,'omitnan'));

    % ---- filter the tracklets, then measure how CONSISTENTLY they travel -------------
    % What separates real flow from chance links is not which way a tracklet points, but
    % whether it keeps going that way. Two measures, both per tracklet:
    %
    %   straightness = |net displacement| / path length.  1 = a straight line, small =
    %       a wandering path. Its chance level DEPENDS ON LENGTH: an n-step random walk
    %       gives about sqrt(pi/(4n)), so a 6-step tracklet scores ~0.36 by luck alone.
    %       That chance level is computed and drawn, so the axis is interpretable.
    %
    %   persistence  = mean cos(angle between consecutive steps).  0 for a random walk
    %       REGARDLESS of length, 1 for perfectly straight travel. This is the fairer of
    %       the two when tracklets differ in length, so it is what the filter uses.
    %
    % Filters on amplitude and length come first: weak and short tracklets are the ones
    % most easily produced by chance, and the velocity estimate uses only the endpoints.
    trkStr = nan(nTrk,1);  trkPer = nan(nTrk,1);
    for q = 1:nTrk
        r  = trkAll{q};
        st = diff(r(:,1:2), 1, 1);                       % steps
        sl = vecnorm(st, 2, 2);
        trkStr(q) = norm(r(end,1:2)-r(1,1:2)) / max(sum(sl), eps);
        ok = sl > eps;
        if nnz(ok) >= 2
            u = st(ok,:) ./ sl(ok);
            trkPer(q) = mean(sum(u(1:end-1,:) .* u(2:end,:), 2));   % mean cos(delta theta)
        end
    end
    chanceStr = sqrt(pi ./ (4*max(trkLen-1,1)));         % random-walk straightness

    cutLen = trkLen >= trkMinLenKeep;                      % the primary filter
    cutAmp = trkAmp >= prctile(trkAmp, trkMinAmpPrc);
    cutPer = trkPer >= trkMinPersist | isnan(trkPer);
    keepTrk = cutLen & cutAmp & cutPer;
    fprintf('   keep filter, each cut applied on its own (of %d tracklets):\n', nTrk);
    fprintf('      length    >= %2d frames : %4d pass  <- primary\n', trkMinLenKeep, nnz(cutLen));
    fprintf('      |amp|     >= p%-2g       : %4d pass\n', trkMinAmpPrc, nnz(cutAmp));
    fprintf('      persistence >= %.2f     : %4d pass\n', trkMinPersist, nnz(cutPer));
    fprintf('      ALL THREE              : %4d kept (%d bright, %d dark)\n', ...
        nnz(keepTrk), nnz(keepTrk & trkPol(:)>0), nnz(keepTrk & trkPol(:)<0));
    fprintf(['   persistence: median %+.3f over all tracklets (0 = random walk); ' ...
             'straightness %.3f vs chance %.3f\n'], median(trkPer,'omitnan'), ...
        median(trkStr,'omitnan'), median(chanceStr));

    if nnz(keepTrk) >= 4
        kIdx = find(keepTrk);

        % Note on colour: an image and a scatter in the SAME axes share one colormap, so
        % painting the points from a scalar makes them take the image's gray map. Every
        % marker colour below is therefore an explicit RGB triplet, and the background
        % image is converted to truecolor, so the two never interfere.
        % Colour encodes WHY a tracklet was rejected (first failing cut, in the order
        % length -> amplitude -> persistence), so a grey cloud can never hide which
        % threshold is actually doing the work.
        cKeep = [0.85 0.15 0.15];                      % red   = kept
        cRejL = [0.45 0.45 0.45];                      % dark grey = too short
        cRejA = [0.20 0.45 0.85];                      % blue  = too weak
        cRejP = [0.95 0.70 0.15];                      % amber = not persistent enough
        why   = zeros(nTrk,1);                         % 0 kept, 1 len, 2 amp, 3 persist
        why(~cutLen) = 1;
        why(cutLen & ~cutAmp) = 2;
        why(cutLen & cutAmp & ~cutPer) = 3;
        colAll = repmat(cKeep, nTrk, 1);
        colAll(why==1,:) = repmat(cRejL, nnz(why==1), 1);
        colAll(why==2,:) = repmat(cRejA, nnz(why==2), 1);
        colAll(why==3,:) = repmat(cRejP, nnz(why==3), 1);
        whyNm = {sprintf('kept (%d)', nnz(why==0)), ...
                 sprintf('too short, <%d (%d)', trkMinLenKeep, nnz(why==1)), ...
                 sprintf('too weak, <p%g (%d)', trkMinAmpPrc, nnz(why==2)), ...
                 sprintf('persistence <%.2f (%d)', trkMinPersist, nnz(why==3))};
        whyCol = [cKeep; cRejL; cRejA; cRejP];
        fprintf('   rejection reasons: %s\n', strjoin(whyNm, ' | '));
        pathUm = nan(nTrk,1);
        for q = 1:nTrk
            r = trkAll{q};
            pathUm(q) = sum(vecnorm(diff(r(:,1:2),1,1), 2, 2)) * pxUm;
        end

        get_figure('R3-trackScatter');
        tlSc = tiledlayout(2, 2, 'TileSpacing','compact','Padding','compact');

        % (1) where the kept tracklets are -- truecolor backdrop, RGB markers
        axs = nexttile;
        gBg = mat2gray(vSm, [0 prctile(vSm(:),99.5)]);
        image(axs, repmat(gBg,1,1,3)); axis(axs,'image','off'); hold(axs,'on');
        contour(axs, bvBW, [0.5 0.5], 'c', 'LineWidth', 0.8);
        for w = [1 2 3 0]                              % kept drawn last, on top
            sel = why == w;
            if ~any(sel), continue; end
            scatter(axs, trkCen(sel,1), trkCen(sel,2), 14 + 14*(w==0), whyCol(w+1,:), ...
                'filled', 'MarkerFaceAlpha', 0.4 + 0.5*(w==0), 'DisplayName', whyNm{w+1});
        end
        legend(axs, 'Location','eastoutside');
        title(axs, 'tracklets coloured by why they were rejected');

        % (2) length vs straightness, against the random-walk chance level
        axs2 = nexttile; hold(axs2,'on');
        scatter(axs2, trkLen, trkStr, 22, colAll, 'filled', 'MarkerFaceAlpha', 0.7);
        [ls, si2] = sort(trkLen);
        plot(axs2, ls, chanceStr(si2), 'k-', 'LineWidth', 1.5);
        xline(axs2, trkMinLenKeep, 'r--', sprintf('length >= %d', trkMinLenKeep));
        xlabel(axs2,'tracklet length (frames)'); ylabel(axs2,'straightness  |net| / path');
        ylim(axs2,[0 1]); grid(axs2,'on'); box(axs2,'off');
        title(axs2, 'length vs straightness (black = random-walk chance)');

        % (3) distance travelled vs speed
        axs3 = nexttile; hold(axs3,'on');
        scatter(axs3, pathUm, trkSpd, 22, colAll, 'filled', 'MarkerFaceAlpha', 0.7);
        xlabel(axs3,'path length travelled (\mum)'); ylabel(axs3,'speed (\mum/s)');
        grid(axs3,'on'); box(axs3,'off');
        title(axs3, 'distance travelled vs speed');

        % (4) width vs amplitude
        axs4 = nexttile; hold(axs4,'on');
        scatter(axs4, trkSig, trkAmp, 22, colAll, 'filled', 'MarkerFaceAlpha', 0.7);
        if trkMinAmpPrc > 0
            yline(axs4, prctile(trkAmp, trkMinAmpPrc), 'r--', sprintf('p%g', trkMinAmpPrc));
        end
        set(axs4,'YScale','log'); xlabel(axs4,'\sigma (\mum)'); ylabel(axs4,'|amplitude|');
        grid(axs4,'on'); box(axs4,'off');
        title(axs4, 'width vs amplitude');

        title(tlSc, sprintf(['f=%d  %d of %d tracklets kept (red): length >= %d, ' ...
            '|amp| >= p%g, persistence >= %.2f'], f, nnz(keepTrk), nTrk, ...
            trkMinLenKeep, trkMinAmpPrc, trkMinPersist));

        % shuffled control: reorder each tracklet's points in time and re-measure
        pShuf = nan(nTrk,1);
        for q = 1:nTrk
            r = trkAll{q};  rp = r(randperm(size(r,1)), 1:2);
            st = diff(rp,1,1);  sl = vecnorm(st,2,2);  ok = sl > eps;
            if nnz(ok) >= 2
                u = st(ok,:)./sl(ok);
                pShuf(q) = mean(sum(u(1:end-1,:).*u(2:end,:), 2));
            end
        end
        fprintf(['   persistence vs a frame-shuffled control: real %+.3f, shuffled %+.3f ' ...
                 '-> %s\n'], median(trkPer,'omitnan'), median(pShuf,'omitnan'), ...
            'a real gap means the tracklets are following something');
    else
        warning(['Only %d tracklets passed. Lower trkMinPersist (now %.2f) or ' ...
            'trkMinAmpPrc/trkMinLenKeep.'], nnz(keepTrk), trkMinPersist);
    end

    % --- directional-consistency control: is this flow, or linked noise? ---
    % Shuffle the frame labels of the detections and re-link. Real flow survives the
    % length/speed test only when positions are ordered in time; noise does not.
    ppS = [];
    for q = 1:nTrk
        r = trkAll{q};
        ppS = [ppS; r(:,1:2), r(randperm(size(r,1)),5)]; %#ok<AGROW>
    end
    ppS = sortrows(ppS, 3);
    nTrkS = NaN;
    try
        lkS = track(ppS, trkMaxDisp, struct('mem',0,'good',trkMinLen,'dim',2,'quiet',1));
        nTrkS = numel(unique(lkS(:,4)));
    catch
        fprintf('   (shuffled control could not be linked; skipped)\n');
    end
    fprintf('   frame-shuffled control: %d tracklets vs %d real (ratio %.2f)\n', ...
        nTrkS, nTrk, nTrkS/max(nTrk,1));

    % --- flow map: binned mean velocity, from the KEPT tracklets only ---
    xb = 1:trkBinPx:ncol;  yb = 1:trkBinPx:nrow;
    [Xb, Yb] = meshgrid(xb + trkBinPx/2, yb + trkBinPx/2);
    Ub = zeros(size(Xb)); Vb = zeros(size(Xb)); Nb = zeros(size(Xb));
    for q = find(keepTrk)'
        bi = max(1, min(numel(yb), ceil(trkCen(q,2)/trkBinPx)));
        bj = max(1, min(numel(xb), ceil(trkCen(q,1)/trkBinPx)));
        Ub(bi,bj) = Ub(bi,bj) + trkVel(q,1);
        Vb(bi,bj) = Vb(bi,bj) + trkVel(q,2);
        Nb(bi,bj) = Nb(bi,bj) + 1;
    end
    Ub = Ub ./ max(Nb,1);  Vb = Vb ./ max(Nb,1);
    okB = Nb >= 2;

    get_figure('R3-spatialBP-tracks');
    tlTk = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
    angC = hsv(180);
    kSet = {find(keepTrk),  sprintf('KEPT (length>=%d, |amp|>=p%g, persist>=%.2f)', ...
                trkMinLenKeep, trkMinAmpPrc, trkMinPersist); ...
            find(~keepTrk), 'REJECTED by the same filter'};

    % (i) kept and (ii) rejected trajectories, same axes and same colour code, so the
    % two panels can be compared directly: real flow should leave the top panel with
    % long, ordered, spatially clustered tracks and the bottom one with short scatter.
    for pnl = 1:2
        axk = nexttile;
        imagesc(axk, vSm, [0 prctile(vSm(:), 99.5)]); axis(axk,'image','off'); hold(axk,'on');
        colormap(axk, gray); contour(axk, bvBW, [0.5 0.5], 'c', 'LineWidth', 0.8);
        idxSet = kSet{pnl,1};
        for q = idxSet(:)'
            r = trkAll{q};
            aa = mod(round(atan2d(trkVel(q,2), trkVel(q,1))), 180) + 1;
            plot(axk, r(:,1), r(:,2), '-', 'Color', angC(aa,:), 'LineWidth', 1.1);
        end
        if isempty(idxSet)
            medLen = NaN;
        else
            medLen = median(trkLen(idxSet));
        end
        title(axk, sprintf('%s: n=%d, median length %.0f frames  (colour = direction)', ...
            kSet{pnl,2}, numel(idxSet), medLen));
    end

    % (iii) binned flow field from the kept tracklets
    axf = nexttile;
    imagesc(axf, vSm, [0 prctile(vSm(:), 99.5)]); axis(axf,'image','off'); hold(axf,'on');
    colormap(axf, gray); contour(axf, bvBW, [0.5 0.5], 'c', 'LineWidth', 0.8);
    spdK = trkSpd(keepTrk);
    qs = trkBinPx / max(prctile(spdK, 90), eps);
    quiver(axf, Xb(okB), Yb(okB), Ub(okB)*qs, Vb(okB)*qs, 0, 'r', 'LineWidth', 1.2, ...
        'MaxHeadSize', 2);
    title(axf, sprintf(['binned mean velocity of the kept tracklets (%d px bins, ' ...
        '>=2 tracklets), arrow = %.0f um/s'], trkBinPx, prctile(spdK, 90)));

    title(tlTk, sprintf(['f=%d  RBC tracklets in S (%.1f-%.1f \\mum band, %s spots, ' ...
        '%d frames)'], f, bandUm(1), bandUm(2), trkPolarity, numel(trkFrames)));
else
    warning(['No tracklets found. Lower trkThrSD, raise trkMaxDisp, or check that S ' ...
             'actually contains moving spots (watch the mp4 first).']);
end
clear Strk


% ---- three-panel COLOUR movie: input | S | S with the tracklets overlaid ----------
% Panel 3 blends the linked trajectories into S as coloured overlays with a fading tail,
% so direction of travel is visible while it plays:
%       RED  = bright-polarity tracklets
%       BLUE = dark-polarity tracklets
% Grayscale would force the two populations to the two ends of the display range, which
% makes them fight with the image itself; colour keeps the underlying S readable. The
% frames are therefore RGB and written straight through VideoWriter rather than writeMov
% (which renders a grayscale figure per frame).
sbMovFile   = fullfile(fpath{f}, sprintf('f%d_spatialBP_%g-%gum_tracks.mp4', f, sbLoUm, sbHiUm));

if isfile(sbMovFile)
    fprintf('[R3] tracked spatial-bandpass movie already exists, skipping (%s)\n', sbMovFile);
else
    nMovF = min([sbMovFrm, nTseg, numel(trkFrames)]);   % stay inside the tracked window
    pA = -double(mov_nv(:,:,1:nMovF));                  % spike-positive input
    pB = -double(S(:,:,1:nMovF));
    pB = pB * (std(pA(:)) / max(std(pB(:)), eps));      % same display scale as panel 1
    caxM = prctile(pA(:), [1 99]);

    % trajectory overlay weights (0..1, fading with tail age), one per polarity
    wR = zeros(nrow, ncol, nMovF, 'single');
    wB = zeros(nrow, ncol, nMovF, 'single');
    nPaint = 0;
    if exist('trkAll','var') && ~isempty(trkAll)
        for q = 1:numel(trkAll)
            r = trkAll{q};
            isBright = trkPol(q) > 0;
            for u = 1:size(r,1)
                rr = round(r(u,2));  cc = round(r(u,1));
                if rr < 1 || rr > nrow || cc < 1 || cc > ncol, continue; end
                ri = max(1,rr-sbMarkPx):min(nrow,rr+sbMarkPx);
                ci = max(1,cc-sbMarkPx):min(ncol,cc+sbMarkPx);
                for w = 0:sbTailFrm
                    ft = r(u,5) + w;
                    if ft < 1 || ft > nMovF, continue; end
                    aw = single(1 - w/(sbTailFrm+1));       % fades with age
                    if isBright, wR(ri,ci,ft) = max(wR(ri,ci,ft), aw);
                    else,        wB(ri,ci,ft) = max(wB(ri,ci,ft), aw);  end
                end
                nPaint = nPaint + 1;
            end
        end
        fprintf('[R3] overlaying %d trajectory points from %d tracklets (%d red, %d blue)\n', ...
            nPaint, numel(trkAll), nnz(trkPol>0), nnz(trkPol<0));
    else
        warning('No tracklets available; panel 3 is a plain copy of S.');
    end

    vw = VideoWriter(sbMovFile, 'MPEG-4');
    vw.FrameRate = sbMovRate;  vw.Quality = 100;  open(vw);
    nrm3 = @(X) min(max((X - caxM(1)) / max(diff(caxM), eps), 0), 1);
    for t = 1:nMovF
        gA = repmat(nrm3(pA(:,:,t)), 1, 1, 3);
        gB = repmat(nrm3(pB(:,:,t)), 1, 1, 3);
        g3 = repmat(nrm3(pB(:,:,t)), 1, 1, 3);
        aR = double(wR(:,:,t));  g3 = g3 .* (1-aR) + reshape(sbColBright,1,1,3) .* aR;
        aB = double(wB(:,:,t));  g3 = g3 .* (1-aB) + reshape(sbColDark  ,1,1,3) .* aB;
        fr = cat(1, gA, gB, g3);
        fr = fr(1:2*floor(size(fr,1)/2), 1:2*floor(size(fr,2)/2), :);   % MPEG-4 wants even dims
        writeVideo(vw, fr);
    end
    close(vw);
    fprintf(['[R3] saved [input | S | S+tracks] -> %s\n' ...
             '     red = bright-polarity tracklets, blue = dark-polarity\n'], sbMovFile);
    clear pA pB wR wB
end

%% (7d2) Image correlation: which pixels carry a FLOWING component, without tracking
% Particle tracking has to RESOLVE one red cell from the next -- pkfnd cannot separate
% two peaks closer than its feature size -- so it thins out exactly where cells are
% densest. Spatiotemporal image correlation (STICS) never identifies a particle at all,
% so it does not care.
%
% ================================ HOW IT WORKS ===================================
% Remove each pixel's time-average and call what is left dI(x,y,t). Then form
%
%       G(xi, eta, tau) = < dI(x, y, t) * dI(x+xi, y+eta, t+tau) >
%
% averaged over every pixel and every start time. Read it as: "if something is HERE
% now, how much is there something THERE, tau frames later?"
%
%   tau = 0   G is a bump at the origin. Its WIDTH is the size of one particle, and
%             shot noise -- white in space and time -- adds a single spike at exactly
%             (0,0) and nowhere else. That is why noise cannot fake flow.
%   tau > 0   a STATIONARY pattern keeps its bump at the origin and merely fades.
%             A FLOWING pattern carries the bump with it: the centre sits at
%             (vx*tau, vy*tau).
%
% So VELOCITY IS THE SLOPE of bump position against lag. Panel 2 below plots exactly
% that, and the number you read off it is the flow speed.
%
% The fit is done in VELOCITY SPACE rather than by following the bump lag by lag,
% because a quasi-regular train of cells makes G periodic in space: on a synthetic
% train 8 px apart, the true peak at +6 px reappears at -2 px and a per-lag maximum
% picks the wrong branch. Scoring each candidate velocity by how much correlation
% lies along the line it predicts avoids that, and how far the winner stands above the
% rest of velocity space is a real confidence -- on synthetic data an empty region
% scored 0 of 15 tiles while the vessels scored every tile they covered.
% =================================================================================
%
% Cost: this is the slow cell. Start with a few hundred frames and a coarse tileStep.

sticsFrames = 1:min(100, nTseg);   % frames used (raise once the settings look right)
sticsOpt = struct( ...
    'tilePx',   16, ...   % tile side. Bigger = faster flow measurable, coarser map
    'tileStep', 4, ...   % 50% overlap
    'maxTau',   32,  ...   % largest lag, frames
    'maxShift', [], ...   % [] = tilePx/2 - 1
    'dv',       0.25, ... % velocity grid, px/frame
    'minZ',     4,  ...   % confidence floor
    'fs',       fs, 'pxUm', pxUm, 'verbose', true);

stFlow = flow_stics(S(:,:,sticsFrames), sticsOpt);

% ---- 1. the maps: where is there flow, how fast, which way ---------------------------
get_figure('R3-stics-maps');  flow_theme(gcf);
tlS2 = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

axA = nexttile;
imagesc(axA, var(double(S(:,:,sticsFrames)),0,3)); axis(axA,'image','off');
colormap(axA, gray); colorbar(axA); hold(axA,'on');
contour(axA, bvBW, [0.5 0.5], 'c', 'LineWidth', 0.8);
title(axA, 'var(S)  -- for reference');

axB = nexttile;
h = imagesc(axB, stFlow.score); set(h,'AlphaData',~isnan(stFlow.score));
axis(axB,'image','off'); colormap(axB, turbo); colorbar(axB); hold(axB,'on');
contour(axB, bvBW, [0.5 0.5], 'w', 'LineWidth', 0.8);
title(axB, sprintf('FLOW SCORE (confidence z, threshold %g)', sticsOpt.minZ));

axC = nexttile;
h = imagesc(axC, stFlow.speed); set(h,'AlphaData',~isnan(stFlow.speed));
axis(axC,'image','off'); colormap(axC, turbo); cb = colorbar(axC);
cb.Label.String = '\mum/s'; hold(axC,'on');
contour(axC, bvBW, [0.5 0.5], 'w', 'LineWidth', 0.8);
ok = isfinite(stFlow.tileSpeed);
qx = stFlow.tileXY(ok,1); qy = stFlow.tileXY(ok,2);
qu = stFlow.tileVx(ok);   qv = stFlow.tileVy(ok);
sc = sticsOpt.tilePx / max(prctile(hypot(qu,qv),90), eps);
quiver(axC, qx, qy, qu*sc, qv*sc, 0, 'r', 'LineWidth', 1.1);
title(axC, sprintf('speed + direction (%d of %d tiles)', nnz(ok), numel(ok)));

axD = nexttile;
if exist('trkSpd','var') && exist('trkCen','var') && any(keepTrk)
    imagesc(axD, var(double(S(:,:,sticsFrames)),0,3)); axis(axD,'image','off');
    colormap(axD, gray); hold(axD,'on');
    scatter(axD, trkCen(keepTrk,1), trkCen(keepTrk,2), 14, trkSpd(keepTrk), 'filled');
    colormap(axD, turbo); cb2 = colorbar(axD); cb2.Label.String = '\mum/s';
    title(axD, sprintf('(7d) tracking, for comparison: %d tracklets', nnz(keepTrk)));
    fprintf(['\n[R3] speed, the two methods on the same field:\n' ...
             '     tracking (7d) : median %.0f um/s from %d tracklets\n' ...
             '     correlation   : median %.0f um/s from %d tiles\n'], ...
        median(trkSpd(keepTrk)), nnz(keepTrk), ...
        median(stFlow.tileSpeed(ok)), nnz(ok));
else
    axis(axD,'off'); text(axD,0.5,0.5,'run (7d) to compare','HorizontalAlignment','center');
end
title(tlS2, sprintf(['f=%d  spatiotemporal image correlation on %d frames ' ...
    '(%d px tiles, lags 0-%d)'], f, numel(sticsFrames), sticsOpt.tilePx, sticsOpt.maxTau));

% ---- 2. THE TEACHING PANEL: watch the correlation bump move --------------------------
% Pick the tile with the strongest flow and show G(xi,eta,tau) lag by lag. If the method
% is working you SEE the bump walk away from the origin, and the walk is a straight line
% whose slope is the speed. This is the whole method in one figure.
[~, ibest] = max(stFlow.tileR2(:));
[iy, ix]   = ind2sub(size(stFlow.tileR2), ibest);
Cb = stFlow.tileC{iy,ix};
if ~isempty(Cb)
    tp2 = sticsOpt.tilePx;  lagAx = -tp2:tp2-1;
    nSh = min(6, numel(stFlow.tau));
    get_figure('R3-stics-howItWorks');  flow_theme(gcf);
    tlS3 = tiledlayout(2, nSh, 'TileSpacing','compact','Padding','compact');
    for k = 1:nSh
        axk = nexttile;
        imagesc(axk, lagAx, lagAx, Cb(:,:,k), [-0.4 1]);
        axis(axk,'image'); colormap(axk, cmapCorr); hold(axk,'on');
        plot(axk, 0, 0, 'k+', 'MarkerSize', 8, 'LineWidth', 1);
        xlim(axk, stFlow.opt.maxShift*[-1 1]); ylim(axk, stFlow.opt.maxShift*[-1 1]);
        if k == 1, ylabel(axk, '\eta lag (px)'); end
        xlabel(axk, '\xi lag (px)');
        title(axk, sprintf('\\tau = %d', stFlow.tau(k)));
    end
    axL = nexttile([1 nSh]); hold(axL,'on');
    vxp = stFlow.tileVx(iy,ix)/(pxUm*fs);   % back to px per frame
    vyp = stFlow.tileVy(iy,ix)/(pxUm*fs);
    pk = nan(numel(stFlow.tau),2);
    for k = 1:numel(stFlow.tau)
        Ck = Cb(:,:,k);
        msk = abs(lagAx) <= stFlow.opt.maxShift;
        Cm = Ck(msk,msk);  ax2 = lagAx(msk);
        [~, im] = max(Cm(:));  [ri,ci] = ind2sub(size(Cm), im);
        pk(k,:) = [ax2(ci) ax2(ri)];
    end
    plot(axL, stFlow.tau, pk(:,1), 'o', 'Color',[0.2 0.4 0.9], 'DisplayName','\xi peak (measured)');
    plot(axL, stFlow.tau, pk(:,2), 's', 'Color',[0.9 0.4 0.2], 'DisplayName','\eta peak (measured)');
    plot(axL, stFlow.tau, vxp*stFlow.tau, '-', 'Color',[0.2 0.4 0.9], 'DisplayName','\xi fitted');
    plot(axL, stFlow.tau, vyp*stFlow.tau, '-', 'Color',[0.9 0.4 0.2], 'DisplayName','\eta fitted');
    yline(axL,  stFlow.opt.maxShift, 'k:', 'search limit');
    yline(axL, -stFlow.opt.maxShift, 'k:');
    box(axL,'off'); legend(axL,'Location','best'); xlabel(axL, '\tau (frames)');
    ylabel(axL, 'peak position (px)');
    title(axL, sprintf(['SLOPE = VELOCITY: %.2f, %.2f px/frame = %.0f um/s ' ...
        '(confidence z = %.1f). Circles that leave the dotted lines are past the ' ...
        'measurable range.'], vxp, vyp, stFlow.tileSpeed(iy,ix), stFlow.tileR2(iy,ix)));
    title(tlS3, sprintf(['HOW IT WORKS: G(\\xi,\\eta,\\tau) for the strongest tile ' ...
        '(centre %.0f, %.0f). The bump WALKS -- its slope is the flow.'], ...
        stFlow.tileXY(ibest,1), stFlow.tileXY(ibest,2)));
end

fprintf(['\n[R3] STICS done. The flow score map answers "which pixel has a flowing ' ...
         'component";\n     the speed map answers "how fast, and which way".\n' ...
         '     A suspiciously exact 0 um/s on a busy vessel means the cell train is too\n' ...
         '     REGULAR for correlation to tell v = 0 from v = spacing/interval -- see the\n' ...
         '     periodicity note in flow_stics.\n']);


%% (7e) Rebuild the vascular movie from the tracks: does the tracking explain S?
% Every linked detection carries a position, an amplitude and a WIDTH, so the tracking
% output is a generative model of the vascular movie:
%
%       P(x,y,t) = SUM_k  a_k(t) * G_sigma_k( x - x_k(t), y - y_k(t) )
%
% Fitting nothing and simply RENDERING that sum is a strong test: if P reproduces S,
% then "moving Gaussian blobs" is a sufficient description of what S contains, and the
% tracking has captured it. If P looks nothing like S, the detections are noise.
%
%   a_k(t)  cntrd's integrated brightness, converted to a PEAK amplitude
%           (for a 2-D Gaussian, integral = 2*pi*sigma^2 * peak), signed by polarity
%   sigma_k cntrd's radius of gyration:  Rg^2 = 2*sigma^2  ->  sigma = sqrt(Rg^2/2)
%           (this answers "is there width information?" -- yes, per detection, and it
%            is what sets the footprint here rather than a single global sigma)
%
% (7e) produces no regressor. It exists to answer two questions that (7f) depends on:
%   1. IS THE TRACKING REAL?  If P reproduces S, then "moving Gaussian blobs" is a
%      sufficient description of what S contains and the detections are not noise. If it
%      does not, nothing downstream is worth building.
%   2. WHERE ARE THE VESSELS?  The detection-density mask -- pixels where particles were
%      repeatedly found -- is the data-driven vessel map that (7f) cuts into segments.
%      It replaces the hand-drawn bvMask with something the data itself defines.

% ------------------------------- parameters -------------------------------
recSigMode  = 'perParticle';  % 'perParticle' = sigma from each detection's Rg
                              % 'global'      = one median sigma for all
recSigGlob  = [];             % um, used when recSigMode = 'global' ([] = median of data)
recPatchSig = 3;              % render each Gaussian out to this many sigma
recDensPrc  = 85;             % percentile of the detection-density map used as the
                              % data-driven "where the particles are" mask
recDensSm   = 3;              % px, smoothing of the density map before thresholding
recUseKeptOnly = true;        % use only the tracklets that passed the (7d) filter.
                              % Without this the rejected ones would still land in P, in
                              % the density mask, and so in the (7f) segments -- the
                              % filter would be cosmetic.
% -- display / output -----------------------------------------------------------
recMovFrm   = 600;            % frames in the 2x2 movie (4 panels x nrow x ncol x this)
recMovSep   = 3;              % px, divider width in that movie
% --------------------------------------------------------------------------

assert(exist('trkAll','var') && ~isempty(trkAll), 'Run the tracking in (7c) first.');
nRecF = numel(trkFrames);

% Which tracklets feed everything downstream. (7d) filters them; without this the
% rejected ones would still land in P, in the detection-density mask, and therefore in
% the (7f) segments -- the filter would be cosmetic.
if recUseKeptOnly && exist('keepTrk','var') && numel(keepTrk) == numel(trkAll)
    trkUse = find(keepTrk(:))';
else
    trkUse = 1:numel(trkAll);
    if recUseKeptOnly
        warning('keepTrk not available; falling back to ALL %d tracklets.', numel(trkAll));
    end
end
if recUseKeptOnly, useTag = 'kept only'; else, useTag = 'all'; end
fprintf('\n[R3] reconstruction uses %d of %d tracklets (%s)\n', ...
    numel(trkUse), numel(trkAll), useTag);
assert(~isempty(trkUse), 'No tracklets left after the (7d) filter.');

% ---- flatten every detection that belongs to a used tracklet ----
detXY = [];  detA = [];  detS = [];  detT = [];
for q = trkUse
    r = trkAll{q};
    sg = sqrt(max(r(:,4),0)/2);                 % sigma in PIXELS
    detXY = [detXY; r(:,1:2)];                  %#ok<AGROW>
    detA  = [detA;  trkPol(q) * r(:,3)];        %#ok<AGROW>  signed integrated brightness
    detS  = [detS;  sg];                        %#ok<AGROW>
    detT  = [detT;  r(:,5)];                    %#ok<AGROW>
end
if strcmp(recSigMode,'global')
    sg0 = recSigGlob;  if isempty(sg0), sg0 = median(detS)*pxUm; end
    detS(:) = sg0 / pxUm;
end
% sigma = 0 means cntrd returned Rg^2 <= 0 for that detection: the intensity-weighted
% second moment inside the window came out non-positive because the noise there beat the
% peak. It is a FIT FAILURE, not a real zero-width particle, and it is purely an SNR
% effect -- measured on a synthetic sigma = 1.8 px Gaussian in unit noise:
%       peak SNR   0.5   1.0   1.5   2.0   3.0   5.0   10
%       Rg^2<=0    22%   24%   18%    8%    3%    0%    0%
%       sigma     2.68  2.22  1.98  1.92  1.89  1.84  1.83
% (signed input is fine, incidentally: subtracting the window minimum biases sigma to
% 2.75 and clipping at 3 SD biases it to 1.06, both worse than feeding the residual as is)
% So these detections are DROPPED rather than clamped -- clamping would have rendered
% them as bright 0.5 px spikes in P. A large failure fraction means the detection
% threshold is too low.
badSig = ~isfinite(detS) | detS <= 0;
if any(badSig)
    fprintf(['     dropping %d of %d detections (%.1f%%) where cntrd returned Rg^2 <= 0 ' ...
             '(low-SNR fit failures)\n'], nnz(badSig), numel(badSig), 100*mean(badSig));
    if mean(badSig) > 0.15
        warning(['%.0f%% of detections failed the width fit -- the detection threshold ' ...
            'trkThrSD = %g is letting through peaks that are barely above noise. ' ...
            'Raise it.'], 100*mean(badSig), trkThrSD);
    end
    detXY(badSig,:) = [];  detA(badSig) = [];  detS(badSig) = [];  detT(badSig) = [];
end
assert(~isempty(detT), 'Every detection failed the width fit; raise trkThrSD.');
fprintf(['\n[R3] rebuilding P from %d detections in %d tracklets; sigma median %.2f um ' ...
         '(FWHM %.2f um), mode = %s\n'], numel(detT), numel(trkAll), ...
    median(detS)*pxUm, 2.355*median(detS)*pxUm, recSigMode);

% ---- how trustworthy are the width and amplitude that P is built from? -----------
% Both come from cntrd, which integrates inside a circular window of diameter
% (trkFeatPx+2). Calibrated on clean synthetic Gaussians, sqrt(Rg^2/2) is UNBIASED only
% when sigma is about 0.18 x the window; a smaller particle picks up background and Rg
% is overestimated, a larger one is truncated and Rg is underestimated:
%       window  7  9  11  15  19  23 px
%       best  1.45 1.75 2.0 2.55 3.05 3.50 px sigma
% At the matched width the peak amplitude recovered as brightness/(2*pi*sigma^2) is
% accurate to 1% over a 100x amplitude range, so the conversion itself is sound -- it is
% the WIDTH that carries the window-dependence, and the amplitude inherits it.
% A constant amplitude error is harmless (the global scale alpha below absorbs it); a
% width-dependent one is not, which is why the check is printed rather than assumed.
cntWin  = trkFeatPx + 2;
sigBest = 0.18 * cntWin;
sigMedPx = median(detS);
fprintf(['     cntrd window %d px -> widths are unbiased near sigma %.2f px; ' ...
         'measured median %.2f px\n'], cntWin, sigBest, sigMedPx);
if sigMedPx < 0.6*sigBest || sigMedPx > 1.5*sigBest
    warning(['Particle width (%.2f px) is far from the %.2f px this cntrd window ' ...
        'measures without bias. Set trkFeatPx to about %d px so the window matches, ' ...
        'or use recSigMode = ''global'' so at least the bias is constant.'], ...
        sigMedPx, sigBest, max(3, 2*round((sigMedPx/0.18 - 2)/2)+1));
end

% ---- render P(x,y,t) ----
P = zeros(nrow, ncol, nRecF, 'single');
for d = 1:numel(detT)
    t = detT(d);
    if t < 1 || t > nRecF, continue; end
    sg = detS(d);   hw = ceil(recPatchSig*sg);
    r0 = max(1, round(detXY(d,2))-hw) : min(nrow, round(detXY(d,2))+hw);
    c0 = max(1, round(detXY(d,1))-hw) : min(ncol, round(detXY(d,1))+hw);
    if isempty(r0) || isempty(c0), continue; end
    [cc, rr] = meshgrid(c0, r0);
    Gk = exp(-((cc-detXY(d,1)).^2 + (rr-detXY(d,2)).^2) / (2*sg^2));
    pk = detA(d) / (2*pi*sg^2);                 % integrated brightness -> peak amplitude
    P(r0, c0, t) = P(r0, c0, t) + single(pk * Gk);
end

% ---- how well does P reproduce S? ----
Sref = -S(:,:,trkFrames);                       % same sign convention as the detections
Pv   = reshape(P,    [], nRecF);
Sv   = reshape(double(Sref), [], nRecF);
Pv   = Pv - mean(Pv, 2);   Sv = Sv - mean(Sv, 2);
alpha = sum(Sv(:).*double(Pv(:))) / max(sum(double(Pv(:)).^2), eps);   % one global scale
Rv   = Sv - alpha*double(Pv);
r2map = reshape(1 - var(Rv,0,2)./max(var(Sv,0,2),eps), nrow, ncol);
% NOTE on how to read these numbers: the OVERALL figure is diluted by the many pixels
% that never contain a particle and are pure noise -- no model can explain those. The
% number to judge the reconstruction on is the one restricted to the pixels where
% particles were actually detected. Validated on synthetic data (5 Gaussian particles,
% sigma 1.8 px, noise added): P recovered sigma to within 8%, the global scale came out
% at 0.99 (so the integrated-brightness -> peak conversion is right), and P explained
% 8.6% of var(S) against a noise-free ceiling of 10.2% -- i.e. 84% of what was there.
fprintf('   P explains %.1f%% of the variance of S overall (global scale %.3g)\n', ...
    100*(1 - var(Rv(:))/var(Sv(:))), alpha);
fprintf('   on bvMask pixels %.1f%%, on the dendrite %.1f%%, on background %.1f%%\n', ...
    100*mean(r2map(bvBW)), 100*mean(r2map(neuronMap>0 & ~bvBW)), ...
    100*mean(r2map(~bvBW & ~(neuronMap>0))));
clear Pv Sv Rv

% ---- detection-density mask: where were particles actually found? ----
densMap = zeros(nrow, ncol);
li = sub2ind([nrow ncol], min(max(round(detXY(:,2)),1),nrow), ...
                          min(max(round(detXY(:,1)),1),ncol));
for d = 1:numel(li), densMap(li(d)) = densMap(li(d)) + 1; end
densSm  = imgaussfilt(densMap, recDensSm);
densBW  = densSm > prctile(densSm(densSm>0), recDensPrc);
fprintf('   detection-density mask: %d px (p%g of the non-empty density map), %.0f%% inside bvMask\n', ...
    nnz(densBW), recDensPrc, 100*mean(bvBW(densBW)));
fprintf('   >> variance of S explained by P INSIDE that mask: %.1f%% (the fair number)\n', ...
    100*mean(r2map(densBW)));

get_figure('R3-trackRecon-maps');
tlR = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
rSets = {var(double(Sref),0,3), 'var(S)', turbo, []; ...
         var(double(P),0,3),    'var(P) rebuilt from tracks', turbo, []; ...
         r2map, 'variance of S explained by P', turbo, [0 max(prctile(r2map(:),99),eps)]; ...
         densSm, sprintf('detection density (mask = p%g)', recDensPrc), turbo, []};
for q = 1:4
    axr = nexttile;  vv = rSets{q,1};
    cl = rSets{q,4};  if isempty(cl), cl = [0 max(prctile(vv(:),99),eps)]; end
    imagesc(axr, vv, cl); axis(axr,'image','off'); hold(axr,'on');
    colormap(axr, rSets{q,3}); colorbar(axr);
    contour(axr, bvBW, [0.5 0.5], 'w', 'LineWidth', 0.8);
    if q == 4, contour(axr, densBW, [0.5 0.5], 'g', 'LineWidth', 1); end
    title(axr, rSets{q,2});
end
title(tlR, sprintf(['f=%d  vascular movie rebuilt from %d tracklets (white = bvMask, ' ...
    'green = density mask)'], f, numel(trkAll)));

% NOTE: P was also tried as a REGRESSOR (project it onto bvMask / the density mask and
% regress that out of the movie). It did not work well -- P is zero wherever detection
% failed, so the regressor has holes exactly where the vascular signal continues -- and
% that block has been removed. P is kept here only as a VALIDATION of the tracking: if
% it reproduces S, then 'moving Gaussian blobs' is a sufficient description of S. The
% regressor that is actually used comes from (7f).

% ---- (3) particle width distribution ----
get_figure('R3-trackRecon-width');
tlW2 = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; histogram(detS*pxUm, 40); box off; grid on;
xlabel('\sigma per detection (\mum)'); ylabel('detections');
xline(median(detS)*pxUm, 'r--', sprintf('median %.2f \\mum', median(detS)*pxUm));
title(sprintf('particle width (FWHM median %.2f \\mum)', 2.355*median(detS)*pxUm));
nexttile; scatter(detS*pxUm, abs(detA), 6, 'filled', 'MarkerFaceAlpha', 0.15);
set(gca,'YScale','log'); box off; grid on;
xlabel('\sigma (\mum)'); ylabel('|integrated brightness|');
title('width vs amplitude');
title(tlW2, sprintf(['f=%d  particle width from cntrd radius of gyration ' ...
    '(Rg^2 = 2\\sigma^2), %d detections'], f, numel(detS)));

clear Pvec

% ---- 2x2 comparison movie, played with moviefixsc --------------------------------
%   [ input (mov_res)     | S                 ]
%   [ P (rebuilt)         | S - P (residual)  ]
% This is a RECONSTRUCTION check, not a regression: it asks whether the tracked
% particles account for what is in S. Panel 4 should look like unstructured noise; any
% moving spot still visible there is a particle the detection missed.
% P carries the global scale alpha so panels 2, 3 and 4 are all in S's units; the input
% panel is rescaled to S's standard deviation for display only (it is a much larger
% signal, which is what the spatial bandpass removed). A bright divider separates the
% quadrants.

nQ  = min(recMovFrm, nRecF);
qS  = single(Sref(:,:,1:nQ));                        % S, spike-positive
qP  = single(alpha) * P(:,:,1:nQ);                   % rebuilt, on S's scale (panel 3)
qI  = single(-mov_res(:,:,1:nQ));                    % raw motion residual
qR  = qS - qP;                                       % S - P (panel 4)
fprintf('     residual RMS / S RMS = %.3f (lower means P captured more of S)\n', ...
    std(qR(:)) / max(std(qS(:)), eps));

sI  = std(qS(:)) / max(std(qI(:)), eps);             % display gain for the input panel
qI  = qI * sI;

climQ = double(prctile(qS(:), [1 99]));
sepV  = single(climQ(2));
rowSep = sepV * ones(recMovSep, 2*ncol + recMovSep, nQ, 'single');
colSep = sepV * ones(nrow, recMovSep, nQ, 'single');
movQ = cat(1, cat(2, qI, colSep, qS), rowSep, cat(2, qP, colSep, qR));
fprintf(['[R3] 2x2 movie: [input | S ; P | S-P], %d frames, %dx%d px. ' ...
         'Scale %.3g..%.3g (input panel gain %.3g)\n'], nQ, size(movQ,1), size(movQ,2), ...
    climQ(1), climQ(2), sI);
moviefixsc(movQ, climQ);
title(sprintf(['f=%d  [input | S ; P | S-P]   P from %d tracklets, explains %.1f%% ' ...
    'of var(S) inside the density mask'], f, numel(trkUse), 100*mean(r2map(densBW))));

clear qS qP qR qI rowSep colSep


%% (7f) Phase-continuous vascular regressor from the flow, via segment + Hilbert
% WHY THE (7e) REGRESSOR IS NOT THE BEST USE OF THE TRACKING
% P in (7e) is a RENDERING of the detections, so it is zero wherever detection failed --
% a missed cell, a frame under threshold, a cell crossing focus. Those dropouts punch
% holes in the regressor exactly where the vascular signal is still present.
%
% RBCs inside one capillary pass a point quasi-periodically, so the local signal is
%       s(t) = A(t) * cos(phi(t)),     dphi/dt = 2*pi*f(t)
% and a regressor written that way runs continuously through the dropouts.
%
% HOW *NOT* TO GET f(t) -- measured, not guessed:
% The obvious estimate is the particle PASSAGE RATE (detections per second). It fails,
% because every missed detection lowers the count, and phi is the running INTEGRAL of f,
% so the error accumulates until the oscillator is out of phase with the signal. On a
% synthetic segment with 34% dropout the rate came out at 20 Hz for a true 40 Hz, and
% that oscillator explained R^2 = 0.17 -- WORSE than (7e)'s rendering at 0.30.
%
% WHAT WORKS: take the phase and amplitude from the MEASURED segment signal (bandpass +
% Hilbert), and use the tracking for what it is reliable at -- telling you WHERE the
% vessel is, and WHICH FREQUENCY BAND to look in (from v/d, which does not depend on how
% many other cells were detected). Same synthetic test: R^2 = 0.82 vs 0.30 for the
% rendering. The band matters: centred right with a +/-60% width it works, and a narrow
% band centred on the biased rate gives 0.06, so oscBandFrac is deliberately generous.
%
% Both cos and sin are emitted. Being in quadrature, a per-pixel fit on the pair matches
% ANY phase offset between the segment and that pixel.

% ------------------------------- parameters -------------------------------
segMinTrk   = 8;     % minimum tracklets for a segment to be used
segMinPx    = 10;    % minimum pixels in a segment
segSmoothMs = 50;    % ms, smoothing of the DIAGNOSTIC-ONLY detection-count trace. It
                     % feeds fCnt (the count-based frequency estimate) and nothing else;
                     % fCnt is printed purely as a warning reference because detection
                     % dropout biases it low, so this value does not affect any result.
oscBandFrac = 0.5;   % bandpass half-width as a fraction of the centre frequency
oscFmin     = 40;    % Hz, HARD FLOOR on every frequency this cell uses. Dendritic Vm
                     % lives below this: bAPs are ~1-2 ms events but their trains, the
                     % subthreshold depolarisations and the 20-ms-average signal the paper
                     % reports all have their power under ~40 Hz, so a regressor centred
                     % (or merely REACHING) below the floor can fit real voltage. RBC
                     % passage is much faster -- v/d with v ~ 0.5-2 mm/s and d ~ 5-15 um
                     % puts it at 40-400 Hz -- so a floor here costs no real vascular
                     % signal. It clamps THREE things, and all three matter:
                     %   1. the search band, so the PSD peak cannot be found below it;
                     %   2. the PASSBAND after the +/-oscBandFrac widening -- without this
                     %      a legal fPeak = 50 Hz still opens the band down to 20 Hz;
                     %   3. the fallback band used when v/d comes out unusable.
                     % Set it to 0 to recover the old behaviour (floor = 1 Hz).
oscDropPinned = true;% drop a segment whose PSD peak lands ON the floor: its real peak is
                     % below oscFmin, so the band was forced onto it rather than measured,
                     % and whatever it fits is as likely to be voltage as blood.
oscSearch   = [1 3]; % the spectral peak is searched in fTrk .* [1/oscSearch(2), oscSearch(2)]
oscScoreMin = 0.01;  % get_oscScore threshold: a PC must look like an oscillation.
                     % Calibrated on synthetic traces -- oscillatory ones score >0.2,
                     % drift and noise <0.06, so 0.10-0.15 separates them.
oscMaxPC    = 10;     % at most this many PCs kept per segment
segBinPix   = 12;    % px per PCA mode SCANNED: a segment is searched with
                     % ceil(nPixels/segBinPix) modes, so a bigger segment gets more.
                     % Set it to roughly one resolution element along the vessel --
                     % below that the extra modes are just noise.
oscCPCA     = true;  % COMPLEX PCA instead of ordinary PCA: band-limit each pixel to
                     % srchBand, take its analytic signal, then SVD the complex matrix.
                     % A travelling wave is rank 1 in that domain instead of rank 2, so
                     % one mode holds the whole wave -- |U| is the envelope along the
                     % vessel and angle(U) is the propagation. See get_cpca.
segMinPxPCA = 20;    % px. Segments at least this big get PCA (they can hold a phase
                     % gradient); smaller ones use their own trace, because their leading
                     % PC would just be the mean with extra noise.
% -- design matrix --------------------------------------------------------------
flowUseDeriv = 0;    % 0 = cos/sin only, 1 = also d/dt, 2 = also d2/dt2.
                     % Mostly redundant: d/dt(A cos) is the sin term weighted by the
                     % instantaneous frequency, so corr(d/dt cos, sin) measured -1.000 at
                     % constant frequency, -0.995 at 15% FM, -0.964 at 40% FM. It only
                     % adds a direction when omega swings a lot. The code measures this
                     % and warns if the design turns near-singular.
% -- display / output -----------------------------------------------------------
segMaxShow  = 4;     % segments plotted
regShowMax  = Inf;   % regressors drawn in R3-flowRegressors (Inf = all of them)
regDropBelowScore = true;  % drop design rows whose re-scored oscScore is below
                     % oscScoreMin, instead of regressing them out anyway
flowMovFrm  = 600;   % frames shown in the 2x2 residual movie
flowMovSep  = 3;     % px, divider width in that movie
% --------------------------------------------------------------------------

assert(exist('detT','var') && ~isempty(detT), 'Run (7e) first (it builds detXY/detA/detT).');
if oscCPCA, pcaNm = 'complex PCA'; else, pcaNm = 'PCA'; end
smK    = max(3, round(segSmoothMs));
nRecF2 = numel(trkFrames);

% ---- vessel segments = connected components of the detection-density mask ----
CCs   = bwconncomp(densBW);
segPx = CCs.PixelIdxList;
segPx = segPx(cellfun(@numel, segPx) >= segMinPx);
nSeg  = numel(segPx);
assert(nSeg > 0, 'No segment reached segMinPx=%d; lower recDensPrc in (7e).', segMinPx);

segOf = zeros(nrow, ncol);
for sgi = 1:nSeg, segOf(segPx{sgi}) = sgi; end
detLin = sub2ind([nrow ncol], min(max(round(detXY(:,2)),1),nrow), ...
                              min(max(round(detXY(:,1)),1),ncol));
detSeg = segOf(detLin);
trkSeg = zeros(numel(trkAll),1);          % 0 = not used / not in any segment
for q = trkUse                             % same tracklet set as (7e)
    cl = sub2ind([nrow ncol], min(max(round(trkCen(q,2)),1),nrow), ...
                              min(max(round(trkCen(q,1)),1),ncol));
    trkSeg(q) = segOf(cl);
end
useSeg = find(arrayfun(@(sgi) nnz(trkSeg==sgi), 1:nSeg) >= segMinTrk);
fprintf('\n[R3] %d density segments, %d with >= %d tracklets\n', nSeg, numel(useSeg), segMinTrk);
assert(~isempty(useSeg), 'No segment reached segMinTrk=%d.', segMinTrk);

% ---- per segment: measured signal, tracking-predicted band, Hilbert phase ----
Svv2   = reshape(-double(S(:,:,trkFrames)), [], nRecF2);
% mov_nv, un-bandpassed, for the browser's source panel only (reshape is a view, so only
% each segment's own pixels are ever copied out of it)
haveMnv = exist('mov_nv','var') && size(mov_nv,3) == nTseg;
if haveMnv, Mnv = reshape(mov_nv, [], nTseg); else, Mnv = []; end
if ~haveMnv
    warning(['mov_nv is not in the workspace, so the regressor browser will show the ' ...
        'S-derived (spatially band-passed) source trace instead. Re-run (7c) to get it.']);
end
oscCos = zeros(numel(useSeg), nRecF2);  oscSin = oscCos;  oscY = oscCos;
segPCall = cell(1,numel(useSeg));
segU     = cell(1,numel(useSeg));   % spatial eigenvectors of the kept PCs
segPCidx = cell(1,numel(useSeg));   % their index in the SVD (not their rank after sorting)
segYsrc  = zeros(numel(useSeg), nRecF2);   % the trace BEFORE bandpass+Hilbert
segNPC   = zeros(1,numel(useSeg));         % how many PCs were scored in that segment
segUsc   = cell(1,numel(useSeg));   % their oscillation scores
segMode  = zeros(1,numel(useSeg));  % 1 = PCA, 2 = PCA fell back to mean, 3 = own trace
segStat = nan(numel(useSeg), 6);   % [nTrk speed spacing fTrk fCnt fPeak]
segPinned = false(1,numel(useSeg));% PSD peak landed on the oscFmin floor
segSrchB = nan(numel(useSeg), 2);  % search band actually used, per segment
segBnd   = nan(numel(useSeg), 2);  % Hilbert passband actually used, per segment
fFloor  = max(oscFmin, 1);         % never let any band reach below this
for k = 1:numel(useSeg)
    sgi  = useSeg(k);
    dsel = detSeg == sgi;
    tq   = detT(dsel);

    % speed of the tracklets in this segment (robust to dropout)
    vMed = median(trkSpd(trkSeg == sgi));

    % inter-particle spacing: nearest neighbour among SIMULTANEOUS detections
    dd = [];
    for t = unique(round(tq))'
        pxy = detXY(dsel & round(detT)==t, :);
        if size(pxy,1) < 2, continue; end
        Dp = squareform(pdist(pxy)) + diag(inf(size(pxy,1),1));
        dd = [dd; min(Dp,[],2)]; %#ok<AGROW>
    end
    spacingNN = median(dd) * pxUm;
    fTrk = vMed / max(spacingNN, eps);                      % Hz, from the flow

    % fTrk = speed / spacing can come out physically impossible: if simultaneous
    % detections sit ~1 px apart (the same cell found twice, or noise passing the size
    % gate) the spacing is far too small and fTrk lands above Nyquist. Clamp it into a
    % usable range and say so, rather than handing an invalid band to get_bandpass.
    fLim = [fFloor, 0.45*fs];
    if ~isfinite(fTrk) || fTrk < fLim(1) || fTrk > fLim(2)
        fprintf(['   seg %3d: v/d = %.3g Hz is outside %g-%g Hz (spacing %.2f um, ' ...
                 'speed %.0f um/s) -> searching the full band instead\n'], ...
            sgi, fTrk, fLim(1), fLim(2), spacingNN, vMed);
        srchBand = [fLim(1), fLim(2)];
    else
        srchBand = [fTrk/oscSearch(2), fTrk*oscSearch(2)];
    end
    srchBand = [max(srchBand(1), fFloor), min(srchBand(2), fs/2 - 1)];
    if srchBand(2) <= srchBand(1)*1.05, srchBand = [fFloor, fs/2 - 1]; end
    segSrchB(k,:) = srchBand;
    cnt  = accumarray(max(min(round(tq),nRecF2),1), 1, [nRecF2 1]);
    fCnt = median(movmean(cnt, smK)) * fs;                  % the BIASED estimate, for reference

    % ---- big segment: PCA.  small segment: its own trace ---------------------------
    % A segment longer than the inter-particle spacing holds an ALTERNATING dark/bright
    % train, so its spatial MEAN cancels: on a synthetic vessel 6x longer than the
    % spacing the mean suppressed the envelope 68-fold and its correlation with the true
    % A(t) fell to 0.05. PCA fixes that without guessing where to cut -- a travelling
    % wave puts its two phases into two leading PCs about 90 deg apart, so the quadrature
    % pair falls out of the decomposition (verified: PC1/PC2 phase difference 90 deg,
    % oscScore 0.99/0.97, while the segment mean scored 0.03).
    %
    % But that only applies to a segment big enough to CONTAIN a phase gradient. Below
    % segMinPxPCA there is no room for one, the leading PC is just the mean with extra
    % noise, and PCA buys nothing -- so those segments use their own filtered trace.
    %
    % Which PCs to keep is decided by get_oscScore -- spectral concentration inside the
    % band, times the autocorrelation side lobe of the band-passed trace. Both are needed:
    % concentration alone scores white noise 0.29 and the side lobe alone scores pure
    % drift 0.999, while the product separates cleanly (sine 0.97, narrowband 0.89,
    % drift 0.04, white 0.02). Anything above oscScoreMin is treated as an oscillation.
    pxIdx = segPx{sgi};
    Yseg  = Svv2(pxIdx, :);
    Yseg  = Yseg - mean(Yseg, 2);
    gS    = ~badfr(trkFrames);
    if numel(pxIdx) >= segMinPxPCA
        % How many modes to LOOK AT scales with the segment: a longer segment spans
        % more independent resolution elements and so can hold more distinct spatial
        % modes, while a short one cannot. Scanning a fixed number would under-search a
        % big segment and waste time on a small one.
        nPCscan = max(2, ceil(numel(pxIdx) / segBinPix));
        nPCscan = min(nPCscan, min(size(Yseg,1), nnz(gS)) - 1);
        if oscCPCA
            % COMPLEX PCA. The block comment above says a travelling wave lands in two
            % real PCs 90 deg apart -- that is a degeneracy, not a feature: the wave's
            % variance is split, neither mode alone is the wave, and which of the two
            % leads depends on noise. On the analytic signal the same wave is rank ONE,
            % with the propagation carried in angle(U) and the envelope in |U|. Measured
            % on a synthetic 60 Hz wave down a 60-px vessel: real PCA 49%/49% across two
            % modes, complex PCA 96% in one, with the phase gradient recovered at
            % corr 1.000. real(pc) and imag(pc) then come out already in quadrature, so
            % the oscillator no longer depends on a second Hilbert transform finding it.
            [Us, pcTr, varPC] = get_cpca(Yseg, fs, srchBand, nPCscan, gS);
            scPC = get_oscScore(real(pcTr), fs, srchBand);
        else
            [Us, Ss, ~] = svds(Yseg(:, gS), nPCscan);
            pcTr  = Us' * Yseg;                          % PC time courses over ALL frames
            varPC = diag(Ss).^2;  varPC = varPC / max(sum(varPC), eps);
            scPC  = get_oscScore(pcTr, fs, srchBand);
        end
        segNPC(k) = numel(scPC);
        keepPC = find(scPC >= oscScoreMin);
        [~, ordPC] = sort(scPC(keepPC) .* varPC(keepPC), 'descend');  % oscillatory AND strong
        keepPC = keepPC(ordPC);
        if numel(keepPC) > oscMaxPC, keepPC = keepPC(1:oscMaxPC); end
        fprintf(['   seg %3d: %4d px  %s over %d modes (px/%d) -> %d of %d kept ' ...
                 '(score %s | var %s%%)\n'], ...
            sgi, numel(pxIdx), pcaNm, nPCscan, segBinPix, numel(keepPC), numel(scPC), ...
            mat2str(round(scPC(keepPC)',2)), mat2str(round(100*varPC(keepPC)')));
        if isempty(keepPC)
            warning(['Segment %d: no PC reached oscScoreMin = %.2f (best %.2f). Falling ' ...
                'back to the segment mean, which may be phase-cancelled.'], ...
                sgi, oscScoreMin, max(scPC));
            y = mean(Yseg, 1);  segMode(k) = 2;
        else
            y = real(pcTr(keepPC(1), :));                % the most oscillatory PC
            segMode(k) = 1;
            segU{k}    = Us(:, keepPC);                  % spatial modes, for the map below
            segUsc{k}  = scPC(keepPC);
            segPCidx{k} = keepPC;                        % TRUE svd index of each kept PC
            segPCall{k} = pcTr(keepPC, :);               % all kept PCs -> extra regressors
        end
    else
        % too small to hold a phase gradient: use the segment's own trace
        y = mean(Yseg, 1);  segMode(k) = 3;
        fprintf('   seg %3d: %4d px  own trace (< segMinPxPCA = %d)\n', ...
            sgi, numel(pxIdx), segMinPxPCA);
    end
    % ---- the browser's "before any filtering" trace, taken from mov_nv --------------
    % NOT from S. y above is the PC time course of S, which has already been through the
    % 6-30 um spatial bandpass, so showing it as the source hides what the bandpass did.
    % The point of the top panel is to see what the regressor was pulled OUT of, so the
    % SAME spatial weights are applied to mov_nv instead -- same pixels, same eigenvector,
    % no spatial filtering. Only the display changes: y still drives the PSD, the
    % bandpass and the Hilbert below, because choosing the vessel's spatial scale is
    % exactly what S is for.
    if haveMnv
        Ynv = -double(Mnv(pxIdx, trkFrames));
        Ynv = Ynv - mean(Ynv, 2);
        if segMode(k) == 1 && ~isempty(segU{k})
            segYsrc(k,:) = real(segU{k}(:,1)' * Ynv);   % complex weights -> in-phase part
        else
            segYsrc(k,:) = mean(Ynv, 1);
        end
    else
        segYsrc(k,:) = y;                          % mov_nv gone: fall back to S
    end
    yg = y;  yg(badfr(trkFrames)) = NaN;
    [Py, Fy] = nanPSD(yg, fs, min(1024, 2^floor(log2(nRecF2/2))));
    Py = Py(:);  Fy = Fy(:);
    inS = Fy >= srchBand(1) & Fy <= srchBand(2);
    if ~any(inS & isfinite(Py)), inS = Fy > fFloor & Fy < fs/2-1; end
    Pm = Py;  Pm(~inS | ~isfinite(Pm)) = -inf;
    [~, ipk] = max(Pm);  fPeak = Fy(ipk);

    % Did the peak land ON the floor? Then the true maximum is below oscFmin, the search
    % just could not go there -- the band was imposed, not measured. Flag it.
    segPinned(k) = oscFmin > 0 && fPeak <= fFloor * 1.05;
    if segPinned(k)
        fprintf(['   seg %3d: PSD peak pinned at the %.0f Hz floor -- its real peak is ' ...
                 'below oscFmin%s\n'], sgi, fFloor, ...
            repmat(' -> DROPPED', 1, oscDropPinned));
    end

    % phase and amplitude from the measured signal, band centred on that peak.
    % The floor is applied AFTER the widening: fPeak = 50 Hz with oscBandFrac = 0.6 would
    % otherwise open the passband down to 20 Hz, straight into the voltage band.
    bnd = fPeak * [1-oscBandFrac, 1+oscBandFrac];
    bnd = [max(bnd(1), fFloor), min(bnd(2), fs/2-1)];
    if ~all(isfinite(bnd)) || bnd(2) <= bnd(1)*1.05
        warning(['Segment %d: degenerate passband [%.3g %.3g] from fPeak = %.3g Hz; ' ...
                 'falling back to the search band.'], sgi, bnd(1), bnd(2), fPeak);
        bnd = srchBand;
    end
    segBnd(k,:) = bnd;
    yb  = get_bandpass(y, fs, bnd);
    yb(~isfinite(yb)) = 0;
    z   = hilbert(yb);
    oscCos(k,:) = abs(z).*cos(angle(z));
    oscSin(k,:) = abs(z).*sin(angle(z));
    oscY(k,:)   = y;
    segStat(k,:) = [nnz(trkSeg==sgi) vMed spacingNN fTrk fCnt fPeak];
end
clear Svv2 Mnv Ynv

fprintf('   seg  nTrk  speed(um/s)  NNspacing(um)  f=v/d(Hz)  f=count(Hz)  f=PSD peak(Hz)\n');
for k = 1:numel(useSeg)
    fprintf('   %3d  %4d  %11.0f  %13.1f  %9.1f  %11.1f  %14.1f\n', ...
        useSeg(k), segStat(k,1), segStat(k,2), segStat(k,3), segStat(k,4), ...
        segStat(k,5), segStat(k,6));
end
fprintf(['   (f=count is shown only as a warning: detection dropout biases it LOW and it ' ...
         'must not drive the phase)\n']);
if oscFmin > 0
    fprintf('   passbands used (floor oscFmin = %g Hz): %s\n', oscFmin, ...
        strjoin(arrayfun(@(k) sprintf('s%d %.0f-%.0f', useSeg(k), segBnd(k,1), segBnd(k,2)), ...
        1:numel(useSeg), 'uni', 0), ', '));
end

% ---- drop the segments whose peak was pinned on the floor --------------------------
% A pinned segment has its true spectral maximum BELOW oscFmin. The oscillator built on
% it is therefore not measuring the RBC train -- it is fitting whatever the lowest legal
% band contains, which on a dendrite is voltage. Keeping it would put a low-frequency
% regressor into the design under a vascular label, which is exactly the failure the
% floor exists to prevent.
if oscDropPinned && any(segPinned)
    kp = ~segPinned;
    fprintf('[R3] dropping %d of %d segments pinned at the %g Hz floor: %s\n', ...
        nnz(segPinned), numel(useSeg), fFloor, mat2str(useSeg(segPinned)));
    useSeg  = useSeg(kp);       segMode  = segMode(kp);     segStat  = segStat(kp,:);
    oscCos  = oscCos(kp,:);     oscSin   = oscSin(kp,:);    oscY     = oscY(kp,:);
    segYsrc = segYsrc(kp,:);    segNPC   = segNPC(kp);      segPinned = segPinned(kp);
    segU    = segU(kp);         segUsc   = segUsc(kp);      segPCidx = segPCidx(kp);
    segPCall = segPCall(kp);    segSrchB = segSrchB(kp,:);  segBnd   = segBnd(kp,:);
    assert(~isempty(useSeg), ['Every segment was pinned at oscFmin = %g Hz: no vascular ' ...
        'oscillation above the floor was found anywhere. Lower oscFmin (and accept the ' ...
        'voltage-contamination risk) or revisit the tracking.'], oscFmin);
elseif any(segPinned)
    warning(['%d segment(s) have their PSD peak pinned at the %g Hz floor and ' ...
        'oscDropPinned = false, so they ARE in the design -- their band was imposed, ' ...
        'not measured, and they may be fitting voltage.'], nnz(segPinned), fFloor);
end

% ---- how much of the measured segment signal does the oscillator explain? ----
gRec2 = ~badfr(trkFrames);
r2osc = nan(1,numel(useSeg));
for k = 1:numel(useSeg)
    yv = oscY(k, gRec2)';  yv = yv - mean(yv);
    X1 = [oscCos(k,gRec2); oscSin(k,gRec2); ones(1,nnz(gRec2))]';
    r2osc(k) = 1 - var(yv - X1*(X1\yv)) / max(var(yv), eps);
end
fprintf('\n[R3] variance of the measured segment signal explained by its oscillator:\n');
for k = 1:numel(useSeg)
    fprintf('   segment %3d : R^2 = %.3f\n', useSeg(k), r2osc(k));
end
fprintf('   median R^2 = %.3f\n', median(r2osc));

% ---- figures ----
% Top row: which segments are used and how each was treated. Rows 2..: the SPATIAL
% EIGENVECTORS of the PCs that passed oscScoreMin, painted on the dendrite. Segments are
% disjoint, so every segment's rank-j kept PC can share one map. A travelling wave shows
% a sign flip across the segment in PC1 and a shifted pattern in PC2 -- that spatial
% structure is exactly what the segment mean would have cancelled.
nPCmap = min(oscMaxPC, max([0, cellfun(@(u) size(u,2), segU)]));
get_figure('R3-flowOscillator');
nShow = min(segMaxShow, numel(useSeg));
tlO = tiledlayout(nShow + 1 + nPCmap, 1, 'TileSpacing','compact','Padding','compact');

axm = nexttile;
segMap = nan(nrow, ncol);
for k = 1:numel(useSeg), segMap(segPx{useSeg(k)}) = segMode(k); end
hS2 = imagesc(axm, segMap, [0.5 3.5]); set(hS2,'AlphaData',~isnan(segMap));
axis(axm,'image','off'); colormap(axm, [0.10 0.55 0.25; 0.85 0.55 0.10; 0.30 0.45 0.85]);
cbS = colorbar(axm, 'Ticks', 1:3, 'TickLabels', ...
    {sprintf('PCA (>=%d px)', segMinPxPCA), 'PCA -> mean fallback', 'own trace (small)'});
hold(axm,'on'); contour(axm, bvBW, [0.5 0.5], 'k', 'LineWidth', 1);
title(axm, sprintf('%d of %d segments used: %d by PCA, %d fallback, %d own trace', ...
    numel(useSeg), nSeg, nnz(segMode==1), nnz(segMode==2), nnz(segMode==3)));

% Under COMPLEX PCA the eigenvector is complex and a signed map would throw away the
% half that matters: angle(U) is the PROPAGATION, and a travelling wave shows up as a
% smooth phase ramp along the vessel. So the phase is drawn on a cyclic colormap with
% |U| as opacity -- bright where the mode actually lives, and the colour sweeping
% through the spectrum along a vessel that carries a wave. A mode that is a standing
% pattern rather than a wave shows two flat phase plateaus 180 deg apart instead.
for j = 1:nPCmap
    axe = nexttile;
    eMap = nan(nrow, ncol);  aMap = zeros(nrow, ncol);  lbl = {};
    for k = 1:numel(useSeg)
        if isempty(segU{k}) || size(segU{k},2) < j, continue; end
        u = segU{k}(:,j);
        u = u / max(abs(u));                             % per-segment normalisation
        if oscCPCA
            eMap(segPx{useSeg(k)}) = angle(u);           % get_cpca already fixed the
            aMap(segPx{useSeg(k)}) = abs(u);             % global rotation
        else
            if sum(u) < 0, u = -u; end                   % fix the arbitrary SVD sign
            eMap(segPx{useSeg(k)}) = u;
            aMap(segPx{useSeg(k)}) = 1;
        end
        lbl{end+1} = sprintf('s%d:%.2f', useSeg(k), segUsc{k}(j));   %#ok<AGROW>
    end
    if oscCPCA
        hE = imagesc(axe, eMap, [-pi pi]);
        set(hE, 'AlphaData', aMap .* ~isnan(eMap));
        colormap(axe, hsv); cbE = colorbar(axe, 'Ticks', [-pi 0 pi], ...
            'TickLabels', {'-\pi','0','\pi'});  cbE.Label.String = 'phase';
        ttl = sprintf('kept cPC%d: phase (colour) x |U| (opacity)', j);
    else
        hE = imagesc(axe, eMap, [-1 1]);  set(hE,'AlphaData',~isnan(eMap));
        colormap(axe, cmapCorr); colorbar(axe);
        ttl = sprintf('kept PC%d eigenvector', j);
    end
    axis(axe,'image','off');
    hold(axe,'on'); contour(axe, bvBW, [0.5 0.5], 'k', 'LineWidth', 0.8);
    title(axe, sprintf('%s  (oscScore per segment: %s)', ttl, strjoin(lbl, ', ')));
end
tO = 1:min(1200, nRecF2);
for k = 1:nShow
    axo = nexttile; hold(axo,'on');
    yy = oscY(k,tO);   yy = (yy-mean(yy))/max(std(yy),eps);
    cc = oscCos(k,tO); cc = (cc-mean(cc))/max(std(cc),eps);
    plot(axo, tO, yy, 'k-', 'LineWidth', 0.8, 'DisplayName','measured S in segment');
    plot(axo, tO, cc, 'r-', 'LineWidth', 1.0, 'DisplayName','A(t)cos\phi(t)');
    box(axo,'off'); xlim(axo,[tO(1) tO(end)]);
    if k == 1, legend(axo,'Location','northeast'); end
    ylabel(axo, sprintf('seg %d', useSeg(k)));
    title(axo, sprintf('R^2 = %.3f  |  %.0f um/s, v/d %.0f Hz, PSD peak %.0f Hz', ...
        r2osc(k), segStat(k,2), segStat(k,4), segStat(k,6)));
    if k == nShow, xlabel(axo,'frame'); end
end
title(tlO, sprintf(['f=%d  flow regressor: phase from the measured segment signal, ' ...
    'band anchored by the tracked flow (v/d)'], f));

% ---- optional: add the time derivative of each oscillator --------------------------
% BE AWARE THIS IS MOSTLY REDUNDANT. For s = A(t)cos(phi),
%       ds/dt = A' cos(phi) - A*omega(t) sin(phi)
% so the derivative is the SIN component already in the design, weighted by the
% instantaneous frequency. Measured on synthetic oscillators:
%       frequency modulation   corr(d/dt cos, sin)   cond(X X')
%             0 %                    -1.000            2.6e5
%            15 %                    -0.995            372
%            40 %                    -0.964             54
% It only becomes a genuinely new direction when omega swings a lot. Set
% flowUseDeriv = 0 to leave it out. When it IS included and the design turns out
% near-singular, the FITTED COMPONENT (and so the residual movie) is still well defined
% -- MATLAB's mldivide handles the rank deficiency -- but the individual coefficients
% stop being interpretable, so do not read the per-regressor betas in that case.

zrow = @(V) (V - mean(V,2)) ./ max(std(V,0,2), eps);
flowReg   = [zrow(oscCos); zrow(oscSin)];
srcNm = cell(1, numel(useSeg));
for k = 1:numel(useSeg)
    if segMode(k) == 1
        if oscCPCA, pcTag = 'cpcaPC'; else, pcTag = 'svdPC'; end
        srcNm{k} = sprintf('seg%d %s%d (%d kept of %d)', useSeg(k), pcTag, ...
            segPCidx{k}(1), numel(segPCidx{k}), segNPC(k));
    elseif segMode(k) == 2
        srcNm{k} = sprintf('seg%d mean (%s found none of %d)', useSeg(k), pcaNm, segNPC(k));
    else
        srcNm{k} = sprintf('seg%d own trace (small)', useSeg(k));
    end
end
flowRegNm = [cellfun(@(v) [v ' | cos'], srcNm, 'uni',0), ...
             cellfun(@(v) [v ' | sin'], srcNm, 'uni',0)];
regSeg    = [1:numel(useSeg), 1:numel(useSeg)];   % which segment each row came from
% The SECOND and later oscillatory modes of a segment are genuinely different spatial
% patterns -- a second vessel in the same component, or a different propagation direction
% -- so they are added directly.
%
% Under complex PCA each of them is a complex time course, i.e. its own quadrature pair,
% so real AND imag go in. Under real PCA there is only the one real trace to add. This is
% the one place the two paths differ in the number of rows they produce.
for k = 1:numel(useSeg)
    if isempty(segPCall{k}) || size(segPCall{k},1) < 2, continue; end
    extra = segPCall{k}(2:end, :);
    qIdx  = 2:size(segPCall{k},1);
    if oscCPCA
        flowReg   = [flowReg; zrow(real(extra)); zrow(imag(extra))];            %#ok<AGROW>
        flowRegNm = [flowRegNm, ...
            arrayfun(@(q) sprintf('seg%d cpcaPC%d (%d kept of %d) | extra re', ...
                useSeg(k), segPCidx{k}(q), numel(segPCidx{k}), segNPC(k)), qIdx, 'uni',0), ...
            arrayfun(@(q) sprintf('seg%d cpcaPC%d (%d kept of %d) | extra im', ...
                useSeg(k), segPCidx{k}(q), numel(segPCidx{k}), segNPC(k)), qIdx, 'uni',0)]; %#ok<AGROW>
        regSeg    = [regSeg, repmat(k, 1, 2*size(extra,1))];                    %#ok<AGROW>
    else
        flowReg   = [flowReg; zrow(extra)];                                     %#ok<AGROW>
        flowRegNm = [flowRegNm, arrayfun(@(q) ...
            sprintf('seg%d svdPC%d (%d kept of %d) | extra PC', useSeg(k), ...
            segPCidx{k}(q), numel(segPCidx{k}), segNPC(k)), qIdx, 'uni',0)];    %#ok<AGROW>
        regSeg    = [regSeg, repmat(k, 1, size(extra,1))];                      %#ok<AGROW>
    end
end
dC = oscCos;  dS = oscSin;
for dv = 1:flowUseDeriv
    dC = gradient(dC')';  dS = gradient(dS')';
    flowReg   = [flowReg; zrow(dC); zrow(dS)];                                  %#ok<AGROW>
    flowRegNm = [flowRegNm, ...
        cellfun(@(v) sprintf('%s | d%d cos', v, dv), srcNm, 'uni',0), ...
        cellfun(@(v) sprintf('%s | d%d sin', v, dv), srcNm, 'uni',0)];         %#ok<AGROW>
    regSeg    = [regSeg, 1:numel(useSeg), 1:numel(useSeg)];                    %#ok<AGROW>
end

if flowUseDeriv > 0
    gD  = ~badfr(trkFrames);
    nSg = numel(useSeg);
    rDS = nan(1,nSg);
    for k = 1:nSg
        rDS(k) = corr(flowReg(2*nSg+k, gD)', flowReg(nSg+k, gD)');   % d1(cos) vs sin
    end
    cnd = cond(flowReg(:,gD) * flowReg(:,gD)');
    fprintf(['\n[R3] derivative regressors added: |corr(d/dt cos, sin)| median %.3f, ' ...
             'design cond = %.3g\n'], median(abs(rDS)), cnd);
    if median(abs(rDS)) > 0.98 || cnd > 1e6
        warning(['The derivative is almost exactly the sin term already present ' ...
            '(|r| = %.3f, cond = %.3g): the flow frequency is nearly constant here, so ' ...
            'it adds no information and only destabilises the coefficients. Consider ' ...
            'flowUseDeriv = 0.'], median(abs(rDS)), cnd);
    end
end
fprintf(['\n[R3] flowReg is %d x %d: %d segments x %d components ' ...
         '(cos, sin%s) -- ready for SeeResiduals\n'], ...
    size(flowReg,1), size(flowReg,2), numel(useSeg), 2*(flowUseDeriv+1), ...
    repmat(', d/dt', 1, flowUseDeriv > 0));

% ---- every regressor: its trace, where it came from, and its spatial footprint -----
% The footprint is the COVARIANCE of each regressor with the movie, computed on the good
% frames only. Covariance is used rather than the fitted beta because the regressors are
% correlated with each other (the derivatives especially), and a beta map splits shared
% structure between them in a way that is not interpretable per row. Covariance is
% marginal: it answers "where in the movie does THIS waveform appear", independently of
% the rest of the design.
%
% Each map is outlined in green at the segment the regressor was derived from, so a
% regressor whose footprint sits far from its own segment is a warning sign -- it is
% picking up something other than that vessel.
gF  = ~badfr(trkFrames);               % frames used for every fit / statistic below
Mvr = single(tovec(mov_res(:,:,trkFrames)));   % same movie the regression below uses
Mvr = Mvr(:, gF);  Mvr = Mvr - mean(Mvr, 2);
Rz  = flowReg(:, gF) - mean(flowReg(:, gF), 2);
covReg7 = double(Mvr * single(Rz)') / max(nnz(gF)-1, 1);      % pixels x nReg
clear Mvr Rz

% Score EVERY row of the final design the same way the PCs were scored. This is the
% direct answer to "why is that regressor in here": a row well above oscScoreMin is a
% real oscillation, a row near or below it got in on the borderline and is mostly noise.
% Note the ranking that admitted it was score x VARIANCE, so a high-variance PC scoring
% just over the threshold can outrank a cleaner but weaker one -- if a regressor looks
% like noise, this table is where that shows up.
[regScore, regF0] = get_oscScore(flowReg, fs, [1 fs/2-1]);
regPk = max(abs(covReg7), [], 1)';

% ---- SECOND CRITERION: does the regressor's spatial map look like the NEURON? -------
% oscScore only asks whether a row looks like an oscillation in time. A row can pass that
% and still be voltage -- and a periodic component sitting exactly on the dendrite is far
% more likely to be voltage the design has latched onto than blood. So correlate each
% regressor's spatial map (its covariance with the movie) against the voltage footprint
% and crush the score when they match: |r| <= oscVoltLo costs nothing, |r| >= oscVoltHi
% multiplies the score by oscVoltMin, which puts it under any sensible oscScoreMin.
%
% The correlation is taken over the neuron and vessel pixels only. Over the whole frame
% both maps are ~0 across the background, which alone would push r towards 1 for
% everything and make the criterion meaningless.
%
% CAVEAT worth checking when the penalty fires on a regressor you believe: a vessel
% running ALONG the dendrite overlaps the footprint, so it can score a high r honestly.
% Both reference maps are computed and printed for exactly that comparison -- blood is
% not spike-locked, so the spike-triggered image separates them where the anatomical
% footprint cannot.
oscVoltPen = true;    % apply the penalty (false = compute and print it only)
oscVoltRef = 'auto';  % WHICH MAP COUNTS AS "VOLTAGE".
                      % 'auto'   = 'sta' when there are more than oscVoltMinSpk spikes,
                      %            'ftprnt' otherwise. This is the right default: the STA
                      %            is the DIRECT spatial signature of voltage and, unlike
                      %            the anatomical footprint, blood cannot contaminate it
                      %            (the RBC train is not spike-locked). The footprint is
                      %            only a proxy -- it flags anything shaped like the
                      %            neuron, including a vessel running along the dendrite,
                      %            so it over-penalises real vascular regressors.
                      %            But an STA over too few spikes is mostly noise, and a
                      %            noisy reference correlates with nothing, which silently
                      %            switches the whole criterion off. Hence the fallback.
                      % 'sta' / 'ftprnt' force one regardless.
oscVoltMinSpk = 20;   % more than this many spikes required before the STA is trusted
oscVoltLo  = 0.3;     % |r| below this: no penalty
oscVoltHi  = 0.6;     % |r| at or above this: full penalty
oscVoltMin = 0.02;    % the multiplier applied at full penalty
% ------------------------------------------------------------------------------------
voltMask = (neuronMap > 0) | bvBW | densBW;
voltFp   = neuronMap ./ max(neuronMap(:), [], 'omitnan');       % 'ftprnt' reference

% The STA is built over the WHOLE segment, not just trkFrames: it is a reference map, so
% nothing ties it to the tracking window, and trkFrames is only 2 s -- rarely enough
% spikes to clear oscVoltMinSpk. segFrames gives 5x more.
voltSTA = [];
spkV    = find(full(Result.spike(1, segFrames)) > 0);
spkV    = spkV(spkV > 3 & spkV <= nTseg - 3);
nSpkV   = numel(spkV);
if nSpkV > oscVoltMinSpk
    accV = zeros(nrow, ncol);
    for q = spkV(:)'
        accV = accV + mean(-double(mov_res(:,:,q:q+2)), 3) ...
                    - mean(-double(mov_res(:,:,q-3:q-1)), 3);
    end
    voltSTA = accV / nSpkV;
    voltSTA = voltSTA / max(abs(voltSTA(:)));
end

switch oscVoltRef
    case 'auto'
        if ~isempty(voltSTA)
            oscVoltRefUse = 'sta';
        else
            oscVoltRefUse = 'ftprnt';
            fprintf(['[R3] voltage reference: only %d spikes in segFrames (<= ' ...
                     'oscVoltMinSpk = %d), so the STA would be noise -> using ' ...
                     '''ftprnt''\n'], nSpkV, oscVoltMinSpk);
        end
    case 'sta'
        if isempty(voltSTA)
            warning(['oscVoltRef = ''sta'' but only %d spikes are available (need more ' ...
                'than oscVoltMinSpk = %d); falling back to ''ftprnt''.'], nSpkV, oscVoltMinSpk);
            oscVoltRefUse = 'ftprnt';
        else
            oscVoltRefUse = 'sta';
        end
    case 'ftprnt'
        oscVoltRefUse = 'ftprnt';
    otherwise
        error('oscVoltRef must be ''auto'', ''sta'' or ''ftprnt''.');
end
if strcmp(oscVoltRefUse,'sta')
    voltRef = voltSTA;
    fprintf('[R3] voltage reference: spike-triggered image over %d spikes\n', nSpkV);
else
    voltRef = voltFp;
end

[wVolt, rVolt] = get_voltPenalty(covReg7, voltRef, voltMask, oscVoltLo, oscVoltHi, oscVoltMin);
rVoltAlt = nan(size(rVolt));
if ~isempty(voltSTA)
    altRef = voltFp;  if strcmp(oscVoltRefUse,'ftprnt'), altRef = voltSTA; end
    [~, rVoltAlt] = get_voltPenalty(covReg7, altRef, voltMask, oscVoltLo, oscVoltHi, oscVoltMin);
end
regScoreOsc = regScore;                       % keep the un-penalised score for the table
if oscVoltPen, regScore = regScore .* wVolt; end

altNm = 'r(sta)';  if strcmp(oscVoltRefUse,'sta'), altNm = 'r(ftprnt)'; end
fprintf(['[R3] every regressor in the design (oscScoreMin %.2f; voltage penalty ref ' ...
         '''%s'', |r| %.2f-%.2f -> x%.2f):\n'], oscScoreMin, oscVoltRefUse, oscVoltLo, ...
         oscVoltHi, oscVoltMin);
fprintf('     %-18s %-8s %-8s %-8s %-9s %-8s %-10s %s\n', 'name', 'oscScore', ...
    'r(volt)', altNm, 'penalty', 'score', 'peak|cov|', 'segment');
[~, ordScore] = sort(regScore, 'ascend');
for r = ordScore'
    flag = '';
    if regScore(r) < oscScoreMin
        flag = '   <- BELOW oscScoreMin';
        if wVolt(r) < 0.5, flag = '   <- KILLED by the voltage map'; end
    end
    fprintf('     %-18s %8.3f %8.3f %8.3f %9.3f %8.3f %10.4g %7d%s\n', flowRegNm{r}, ...
        regScoreOsc(r), rVolt(r), rVoltAlt(r), wVolt(r), regScore(r), regPk(r), ...
        useSeg(regSeg(r)), flag);
end
if ~isempty(voltSTA)
    fprintf('     (the two r columns disagreeing means the vessel overlaps the dendrite -- ');
    fprintf('trust r(sta))\n');
end
fprintf(['     (rows below the threshold are derivatives or re-scored PCs; a PC scored ' ...
         'on its own\n      can fall under it once z-scored and differentiated)\n']);

% ---- DROP the rows that do not clear oscScoreMin ---------------------------------
% Until now oscScoreMin only gated which PCs were SELECTED; every row that reached
% flowReg was regressed out regardless of how it scored afterwards. That is not what the
% threshold is for -- a row scoring 0.03 is noise, and regressing it out of the movie can
% only remove real signal. cos and sin are one quadrature UNIT, so they are judged
% together (kept if either passes); extra PCs and derivatives are judged on their own.
regKeep = regScore >= oscScoreMin;
nSegF   = numel(useSeg);
for k = 1:nSegF                                    % tie each cos to its sin partner
    regKeep([k, k+nSegF]) = any(regKeep([k, k+nSegF]));
end
if regDropBelowScore && any(~regKeep)
    fprintf('[R3] dropping %d of %d regressors below oscScoreMin = %.2f:\n', ...
        nnz(~regKeep), numel(regKeep), oscScoreMin);
    for r = find(~regKeep)'
        fprintf('       %-46s score %.3f\n', flowRegNm{r}, regScore(r));
    end
    flowReg = flowReg(regKeep,:);   flowRegNm = flowRegNm(regKeep);
    regSeg  = regSeg(regKeep);      regScore  = regScore(regKeep);
    regF0   = regF0(regKeep);       covReg7   = covReg7(:,regKeep);
    rVolt   = rVolt(regKeep);       rVoltAlt  = rVoltAlt(regKeep);
    wVolt   = wVolt(regKeep);       regScoreOsc = regScoreOsc(regKeep);
    regPk   = regPk(regKeep);
    assert(~isempty(flowReg), 'Every regressor fell below oscScoreMin = %.2f.', oscScoreMin);
elseif any(~regKeep)
    fprintf(['[R3] %d regressors are below oscScoreMin but regDropBelowScore = false, ' ...
             'so they ARE still regressed out of the movie\n'], nnz(~regKeep));
end

% ---- order by segment, so the browser walks through the segmentation --------------
[~, ordSeg] = sortrows([regSeg(:), (1:numel(regSeg))']);
flowReg = flowReg(ordSeg,:);  flowRegNm = flowRegNm(ordSeg);
regSeg  = regSeg(ordSeg);     regScore  = regScore(ordSeg);
regF0   = regF0(ordSeg);      covReg7   = covReg7(:,ordSeg);
rVolt   = rVolt(ordSeg);      rVoltAlt  = rVoltAlt(ordSeg);
wVolt   = wVolt(ordSeg);      regScoreOsc = regScoreOsc(ordSeg);
regPk   = regPk(ordSeg);
segYsrcReg = segYsrc(regSeg,:);            % the un-filtered source of each regressor
fprintf('[R3] %d regressors, ordered by segment: %s\n', size(flowReg,1), ...
    mat2str(useSeg(regSeg)));


% ---- interactive: one regressor at a time, Prev/Next to step through them ----------
% A stacked figure makes each regressor too small to judge; browse_regressors shows one
% at full size with its trace above and its spatial map below. Navigate with the
% buttons, the index box, or the left/right arrow keys.
segIdxForReg = arrayfun(@(r) regSeg(r), 1:size(flowReg,1));

% ---- per-regressor statistics, for the scatter panel and the readout ---------------
% Everything already computed above, gathered in one place so a regressor can be judged
% RELATIVE to the others rather than on its own: an outlier in frequency, an amplitude
% far above the rest, or the one map that correlates with the neuron.
%   f0Hz / bandLo / bandHi   where in frequency it lives, and the passband it came from
%   oscScore / score         the oscillation score before and after the voltage penalty
%   rVolt / rVoltAlt         correlation of its SPATIAL MAP with the voltage reference
%                            (rVoltAlt is the other reference -- STA vs ftprnt)
%   voltPen                  the multiplier that correlation cost it
%   ampCov                   peak |covariance| with the movie: how strongly it loads
%   ampEnv                   the oscillator's own envelope in the segment, in -S units
%   varPct                   share of the movie variance its own segment's fit removes
segAmp = median(sqrt(oscCos.^2 + oscSin.^2), 2);        % envelope per segment
statFlow = struct( ...
    'f0Hz',     regF0(:), ...
    'oscScore', regScoreOsc(:), ...
    'score',    regScore(:), ...
    'rVolt',    rVolt(:), ...
    'rVoltAlt', rVoltAlt(:), ...
    'voltPen',  wVolt(:), ...
    'ampCov',   regPk(:), ...
    'ampEnv',   segAmp(regSeg(:)), ...
    'bandLo',   segBnd(regSeg(:), 1), ...
    'bandHi',   segBnd(regSeg(:), 2), ...
    'segment',  useSeg(regSeg(:))');

browse_regressors(flowReg, covReg7, [nrow ncol], flowRegNm, struct( ...
    'preTrace', segYsrcReg, ...
    'score',    regScore, ...
    'scoreMin', oscScoreMin, ...
    'f0',       regF0, ...
    'seg',      segIdxForReg, ...
    'segMask',  {cellfun(@(k) segPx{k}, num2cell(useSeg), 'uni', 0)}, ...
    'refMask',  bvBW, ...
    'cmap',     cmapCorr, ...
    'tWin',     1:min(1000, nRecF2), ...
    'stats',    statFlow, ...
    'statXY',   {{'f0Hz','score'}}, ...     % the two that decide whether a row survives
    'statC',    'rVolt', ...                % ...coloured by what the penalty saw
    'statLog',  {{'ampCov','ampEnv'}}, ...
    'figName',  'R3-flowRegressors'));
fprintf(['[R3] R3-flowRegressors is interactive: Prev/Next, the index box, or the ' ...
         'arrow keys step through all %d regressors\n'], size(flowReg,1));

% ---- what does the flow regressor actually take out of the movie? -----------------
%   [ input (mov_res)  | S                       ]
%   [ fitted flow      | input - fitted          ]
% The fit and the display both use mov_res -- the RAW motion residual, with the neuron's
% voltage still in it -- not the voltage-regressed mov_nv. That is the honest final test:
% what does this regressor actually do to the data you started from?
%
% The cost is that the voltage signal is now present for the regressor to eat if the two
% are correlated, which is exactly what mov_nv was protecting against. So the number to
% watch below is the variance removed ON THE DENDRITE AWAY FROM VESSELS: it should stay
% small while the bvMask / density-mask figures are large. If it does not, the regressor
% is taking neural signal and the band or the segment selection needs tightening.
%
% Coefficients are fitted on the good tracked frames only (so spikes cannot bias them)
% and then applied everywhere, the same pattern used elsewhere in this script. The
% variance the regressor removes at each pixel is b'*Sf*b with Sf the regressor
% covariance -- exact, and it avoids materialising a residual movie over all frames.
%
% Note this is a fairer test than the per-segment R^2 above: that one is scored on the
% very segment signal the oscillator was derived from, while this is the whole movie.

[~, Cf] = SeeResiduals(mov_res(:,:,trkFrames(gF)), flowReg(:,gF));
Cfv  = tovec(Cf);                                   % pixels x (1 + nReg), col 1 = DC
Bf   = Cfv(:, 2:end);                               % regressor coefficients

varIn   = var(mov_res(:,:,trkFrames(gF)), 0, 3);
Sf      = cov(flowReg(:, gF)');
varFit  = reshape(sum((Bf * Sf) .* Bf, 2), nrow, ncol);
dropF   = varFit ./ max(varIn, eps);
fprintf(['\n[R3] flow regressor removes, from mov_res (voltage still in it): ' ...
         '%.1f%% overall, %.1f%% in the density mask, %.1f%% on bvMask\n'], ...
    100*mean(dropF(:)), 100*mean(dropF(densBW)), 100*mean(dropF(bvBW)));
fprintf('     %.1f%% on the DENDRITE away from vessels  <- keep this one small\n', ...
    100*mean(dropF(neuronMap>0 & ~bvBW)));
if mean(dropF(neuronMap>0 & ~bvBW)) > 0.5*mean(dropF(bvBW))
    warning(['The regressor removes %.1f%% on the dendrite vs %.1f%% on bvMask -- that ' ...
        'ratio is high enough that it is probably taking neural signal. Narrow ' ...
        'oscBandFrac or raise oscScoreMin.'], ...
        100*mean(dropF(neuronMap>0 & ~bvBW)), 100*mean(dropF(bvBW)));
end


nQF = min(flowMovFrm, nRecF2);
IallF = [ones(1,nQF); flowReg(:,1:nQF) - mean(flowReg(:,gF), 2)];
fitQ  = reshape(single(Cfv) * single(IallF), nrow, ncol, nQF);
fI    = single(-mov_res(:,:,trkFrames(1:nQF)));
fS    = single(-S(:,:,trkFrames(1:nQF)));
fF    = -fitQ;                                       % fitted flow, spike-positive
fR    = fI - fF;                                     % input minus fitted

sF    = std(fS(:)) / max(std(fI(:)), eps);           % one display gain for panels 1,3,4
fI = fI*sF;  fF = fF*sF;  fR = fR*sF;
climF = double(prctile(fS(:), [1 99]));
sepF  = single(climF(2));
rowSepF = sepF * ones(flowMovSep, 2*ncol + flowMovSep, nQF, 'single');
colSepF = sepF * ones(nrow, flowMovSep, nQF, 'single');
movF = cat(1, cat(2, fI, colSepF, fS), rowSepF, cat(2, fF, colSepF, fR));
fprintf('[R3] 2x2 movie: [input | S ; fitted flow | input-fitted], %d frames\n', nQF);
get_figure('[R3] 2x2 movie: [input | S ; fitted flow | input-fitted]');
moviefixsc(movF, climF);
title(sprintf(['f=%d  [input | S ; fitted flow | input-fitted]   %d segments x %d ' ...
    'components, removes %.1f%% of the input variance in the density mask'], ...
    f, numel(useSeg), 2*(flowUseDeriv+1), 100*mean(dropF(densBW))));
clear fI fS fF fR fitQ rowSepF colSepF Cfv IallF



%% (7f-ALT) Flow regressor v2: ring voltage reference, then ICA / disk tiling per segment
% An alternate to (7f). Four things change; (7g) and (7h) run on the result unchanged.
%
% -------------------------------------------------------------------------------------
% 1. THE VOLTAGE REFERENCE IS TAKEN FROM BESIDE EACH VESSEL, NOT FROM THE WHOLE ROI.
%
% (7c) built mov_nv by regressing out PCresult{f}.normTraces -- the per-ROI voltage
% traces. Those are read off the WHOLE footprint, and the footprint overlaps the vessels,
% so the "voltage" being removed already has blood in it. Regressing it out takes some of
% the vascular signal away before the tracking ever sees it, and leaves a voltage estimate
% that is itself part vascular.
%
% Here the reference is built where blood cannot be: a RING just outside each density
% segment, intersected with the neuron's footprint. Those pixels are dendrite adjacent to
% the vessel -- close enough to carry the same membrane potential, far enough to carry
% none of that vessel's RBC train. One trace per (segment, ROI) pair that overlaps enough.
%
% WHY THE RING EXCLUDES EVERY VESSEL, not just its own segment: a ring drawn around
% segment A can easily land on segment B, or on a bvMask vessel that was never dense
% enough to become a segment. Either would put blood back into the reference. So densBW
% and bvBW are BOTH removed from every ring.
%
% 2. THE SMOOTHING IS PART OF THE ARGUMENT, NOT INCIDENTAL. A ring trace contains
% dendritic Vm (slow) plus shot noise at every frequency, and regressing raw shot noise
% out of the movie removes variance everywhere for nothing. nvSmoothMs keeps the part
% that is actually voltage. Read the note at nvSmoothMs before changing the ring width.
%
% 3. S IS REBUILT FROM THE NEW mov_nv. S is what the decomposition reads, and (7c) built
% it from the OLD mov_nv -- leaving it alone would apply the new voltage removal to
% nothing. The same spatial filter (kDoG / Hmask, still in the workspace from (7c)) is
% reused, so the voltage step is the ONLY difference between S and S_alt.
%
% 4. THE DECOMPOSITION. Three methods via get_segModes; altCompare scores all three on
% every segment so they can be read side by side.
%
%   'cpca'  what (7f) does -- splits the segment by VARIANCE.
%   'ica'   splits by INDEPENDENCE instead. Variance is the wrong criterion when one
%           connected component holds TWO capillaries of similar amplitude: the leading
%           PCs are then mixtures of both and neither is a vessel. Two RBC trains at
%           different speeds are close to independent, so ICA can separate what PCA
%           blends. See the collinearity check below for the specific way it fails.
%   'disk'  no decomposition at all: overlapping disks, one local trace each. Assumes
%           nothing about how many sources the segment holds.
%
% ONE DIFFERENCE THAT APPLIES TO ALL THREE: in (7f) only the leading mode became a cos/sin
% pair while the extra complex modes contributed real/imag. Here EVERY kept mode gets its
% own Hilbert quadrature pair -- the ICA and disk modes are real and have no imaginary
% half to donate. Row bookkeeping is therefore uniform: one pair per kept mode.
%
% THE VOLTAGE CRITERION IS A HARD VETO HERE, not the graded penalty (7f) applies. A
% regressor whose spatial map looks like the neuron is not a weak vascular regressor to be
% down-weighted; it is the wrong thing, and subtracting it removes signal.
% -------------------------------------------------------------------------------------

% ------------------------------- parameters -------------------------------
altDecomp   = 'disk'; % 'local' | 'cpca' | 'ica' | 'disk'  -- what goes into the design
                       %
                       % 'local' IS THE ONE THAT CANNOT PRODUCE THE BACKGROUND BAND. The
                       % other three build a design and fit it to EVERY pixel, so a
                       % segment's regressor is free to remove variance anywhere in the
                       % frame -- which is what the background-vs-chance row measured, and
                       % why the % removed map has a smooth arc under it. 'local' takes the
                       % same complex-PCA modes but fits and subtracts them ONLY inside the
                       % segment they came from. A pixel outside every segment is not
                       % touched at all, so background removal is exactly zero by
                       % construction rather than by luck, and the dendrite away from
                       % vessels is untouched for the same reason.
                       % The rows are still band-limited to the segment's own passband, so
                       % the correction is confined in frequency as well as in space: it
                       % cannot reach sub-oscFmin voltage even inside the segment.
altCompare  = true;    % also score the other methods on every segment and print the table
altExport   = true;    % publish under the (7f) names so (7g)/(7h) use THIS design.
                       % false = leave (7f)'s variables alone; only the alt* copies exist.
% -- 1/2: the ring voltage reference ----------------------------------------------
nvRingIn    = 1;     % px. Gap between the density mask and the inner edge of the ring.
                     % densBW is a threshold on a SMOOTHED density map, so its boundary is
                     % soft: the pixel just outside it still carries vessel signal.
nvRingOut   = 8;     % px, outer edge. Wide enough to catch footprint, narrow enough that
                     % the dendrite it lands on is electrotonically the same place as the
                     % dendrite beside the vessel.
nvFtprntThr = 0.2;   % fraction of an ROI's own peak weight below which its footprint is
                     % treated as absent. Same convention as bvBW in the loading cell.
nvMinPx     = 50;    % px. Minimum size of a CONNECTED PATCH of ring-and-footprint overlap.
                     % Each ROI's overlap with each ring is split into its connected
                     % components and every component this big becomes its own reference
                     % trace -- see the block comment at the split. Below it the patch is
                     % dropped, not merged into a neighbour: a small patch is mostly shot
                     % noise, and merging two distant patches averages two different
                     % points on the neurite into a potential that is neither.
nvSmoothMs  = 5;     % ms, moving average on each reference trace before regression.
                     % NOTE WHAT THIS DOES AND DOES NOT REMOVE: a 5 ms boxcar is -3 dB near
                     % 90 Hz and first-nulls at 200 Hz, so it does NOT confine the
                     % reference to the sub-40 Hz voltage band -- it still passes power
                     % through the lower half of the RBC range (oscFmin = 40 Hz up). That
                     % is acceptable ONLY because the ring excludes every vessel pixel, so
                     % there is no RBC train in the trace for it to pass. If you widen the
                     % ring or drop the densBW/bvBW exclusion, raise this a lot.
nvDedupR    = 0.95;  % |corr| above which two reference traces count as the same trace.
                     % Rings around neighbouring segments overlap the same piece of
                     % dendrite and would otherwise contribute it several times over.
nvMaxTrace  = 40;    % hard cap. Every reference trace is a degree of freedom taken out of
                     % the movie; past a few dozen the removal starts to be measurable on
                     % the dendrite for its own sake.
% -- 4: the decompositions ---------------------------------------------------------
icaMaxIC    = 8;     % at most this many ICs per segment
icaSeed     = 0;     % FastICA starts from a random guess; fix it so a re-run reproduces
icaPCfloor  = 3;     % skip ICA on a segment that cannot support this many PCs
diskSpacing = 6;     % px, MINIMUM separation between disk centres. Centres are placed by
                     % farthest-point sampling over the segment's OWN pixels, so this is
                     % what actually decides how many disks a segment gets: a long vessel
                     % gets many, a small blob gets few, automatically and in proportion.
diskR       = 4;     % px, disk radius. diskR > diskSpacing/2 means neighbours OVERLAP,
                     % deliberately: a vessel crossing the gap between two centres would
                     % otherwise fall between them and be measured by neither.
% -- 4: the S_alt subtraction ---------------------------------------------------------
% NOT A REGRESSION AND NOT A NOTCH. For each segment: take S_alt's OWN trace at each
% pixel, band-pass it to that segment's measured band, scale it, and subtract it from
% mov_res. Only inside the segment, only in its band.
%
% What makes it different from the other methods, and why it is worth trying:
%
%   vs 'local' (and every row-based method) -- those subtract ONE waveform per segment,
%     A(t)cos(phi(t)), fitted per pixel in amplitude and phase. That assumes every pixel
%     in the segment shares an envelope. It does not: A(t) varies pixel to pixel, and the
%     part that does not match the shared waveform is exactly what those methods leave
%     behind. Here every pixel gets its OWN band-limited trace, so a per-pixel envelope
%     costs nothing to represent.
%
%   vs a notch on mov_res -- a notch removes ALL power in the band, voltage included; a
%     spike is a 1-2 ms event and has real power from 40 Hz up, so a segment crossing the
%     dendrite would lose some of it. Here what is subtracted has already been through
%     the 6-30 um SPATIAL bandpass, which suppresses the dendritic Vm because it is smooth
%     over the length constant. The protection against eating voltage is SPATIAL, and the
%     temporal band only decides which part of the vascular signal to take.
%
% THE CAVEAT THAT DECIDES WHETHER IT WORKS: (7c) warns that the DoG also responds to the
% dendrite's EDGES, so S_alt is not voltage-free there. If this method eats voltage, that
% is where it will happen, and the dendrite column of the removal table is the test.
altSfilt   = true;   % compute it and show it beside the others
sfiltGain  = 'fit';  % 'fit' = one scalar per pixel, least squares against mov_res. NEEDED,
                     % because S_alt is a DoG-filtered movie and the DoG's gain at the
                     % vessel's spatial scale is not 1 -- subtracting raw would be off by
                     % that factor. One degree of freedom per pixel, against 2 per mode for
                     % the row-based methods, so this is the most conservative fit here.
                     % 'one' = subtract as-is; only sensible if you have a reason to
                     % believe the filter gain is unity, which you generally do not.
% -- 4: the 'local' method -----------------------------------------------------------
localDilate = 0;     % px the segment is grown by before its pixels are fitted. 0 keeps
                     % the correction strictly inside densBW, which is what makes the
                     % background removal identically zero. A small value (1-2) is
                     % defensible because densBW is a threshold on a smoothed map and the
                     % vessel bleeds a little past it -- but every pixel added is a pixel
                     % where the guarantee no longer holds, so raise it deliberately.
                     % Where two grown segments overlap, the pixel is fitted by whichever
                     % segment reaches it first and the run reports how many that was.
diskMinPx   = 7;     % px of segment inside a disk for it to be used
diskMaxPer  = 24;    % CEILING on disks per segment, NOT a target. diskSpacing is meant to
                     % be what stops the tiling; if this cap binds instead, one end of the
                     % segment is left unsampled -- which is exactly how the disks came
                     % out bunched in the fat part before. The run prints when it binds.
% -- scoring and the voltage veto ---------------------------------------------------
altScoreMin = oscScoreMin;   % inherited from (7f) so the two are judged the same way
altVoltCut  = 0.5;   % |r| between a row's spatial map and the voltage reference at or
                     % above which the row is DROPPED, whatever it scored. A hard cut, not
                     % a ramp: there is no |r| at which eating some voltage becomes an
                     % acceptable trade, so there is nothing for a multiplier to express.
altDropBelowScore = true;
% -- display -------------------------------------------------------------------------
altSegShow  = 4;     % segments plotted in the trace figure
altModeMap  = 4;     % kept modes drawn as spatial maps
altMovFrm   = 600;   % frames in the uncorrected-vs-corrected movie
altMovSep   = 3;     % px, divider width in that movie
altMovMark  = true;  % outline densBW on every panel of both movies. Worth having on: the
                     % whole point of watching them is whether the flicker stops WHERE THE
                     % VESSELS ARE, and without the outline you are judging that from
                     % memory of a different figure.
% --------------------------------------------------------------------------

assert(exist('densBW','var') && ~isempty(densBW), 'Run (7e) first (it builds densBW).');
assert(exist('detT','var') && ~isempty(detT), 'Run (7e) first (it builds detXY/detT).');
assert(exist('sbMethod','var') == 1, 'Run (7c) first: sbMethod / kDoG / Hmask are needed.');
assert(ismember(altDecomp, {'local','ica','disk','cpca'}), ...
    'altDecomp must be ''local'', ''cpca'', ''ica'' or ''disk''.');
% 'local' differs from the others in HOW THE DESIGN IS FITTED, not in how the modes are
% found -- its modes ARE complex PCA's. So it maps onto 'cpca' everywhere a decomposition
% is asked for, and is distinguished only where the fit happens.
decompOf = containers.Map({'local','cpca','ica','disk'}, {'cpca','cpca','ica','disk'});
if strcmp(altDecomp,'ica')
    assert(~isempty(which('fastica')) || ~isempty(which('rica')), ...
        'altDecomp = ''ica'' needs fastica (FastICA_25) or rica (Statistics Toolbox).');
end
npxA   = nrow * ncol;
nRecF2 = numel(trkFrames);
gF     = ~badfr(trkFrames);
smKa   = max(3, round(segSmoothMs));

% =====================================================================================
% 1. VESSEL SEGMENTS FROM THE DENSITY MASK  (same as (7f), repeated so this cell stands
%    on its own and does not depend on (7f) having been run)
% =====================================================================================
CCa    = bwconncomp(densBW);
segPxA = CCa.PixelIdxList;
segPxA = segPxA(cellfun(@numel, segPxA) >= segMinPx);
nSegA  = numel(segPxA);
assert(nSegA > 0, 'No segment reached segMinPx = %d; lower recDensPrc in (7e).', segMinPx);

segOfA = zeros(nrow, ncol);
for sgi = 1:nSegA, segOfA(segPxA{sgi}) = sgi; end
detLinA = sub2ind([nrow ncol], min(max(round(detXY(:,2)),1),nrow), ...
                               min(max(round(detXY(:,1)),1),ncol));
detSegA = segOfA(detLinA);
trkSegA = zeros(numel(trkAll),1);
for q = trkUse
    cl = sub2ind([nrow ncol], min(max(round(trkCen(q,2)),1),nrow), ...
                              min(max(round(trkCen(q,1)),1),ncol));
    trkSegA(q) = segOfA(cl);
end
useSegA = find(arrayfun(@(sgi) nnz(trkSegA==sgi), 1:nSegA) >= segMinTrk);
assert(~isempty(useSegA), 'No segment reached segMinTrk = %d.', segMinTrk);
nUseA = numel(useSegA);
fprintf('\n[R3-alt] %d density segments, %d with >= %d tracklets\n', nSegA, nUseA, segMinTrk);

% =====================================================================================
% 2. THE RING REFERENCE: ftprnt just outside each segment, with every vessel removed
% =====================================================================================
vesselAll = densBW | bvBW;                 % nothing from here may enter a ring
nROIa     = size(ftprnt, 3);
seOut = strel('disk', nvRingOut);
if nvRingIn > 0, seIn = strel('disk', nvRingIn); else, seIn = []; end

nvW     = zeros(npxA, 0);     % the non-vessel footprints, one column each
nvNm    = {};
nvOfSeg = [];
ringAll = false(nrow, ncol);  % union of all rings, for the figure
for k = 1:nUseA
    sgi  = useSegA(k);
    mSeg = false(nrow, ncol);  mSeg(segPxA{sgi}) = true;
    if isempty(seIn), inner = mSeg; else, inner = imdilate(mSeg, seIn); end
    ring = imdilate(mSeg, seOut) & ~inner & ~vesselAll;
    if ~any(ring(:)), continue, end
    ringAll = ringAll | ring;
    for j = 1:nROIa
        fj = double(ftprnt(:,:,j));
        mx = max(fj(:));
        if ~isfinite(mx) || mx <= 0, continue, end
        ov = (fj >= nvFtprntThr * mx) & ring;
        if ~any(ov(:)), continue, end
        % ONE ROI'S OVERLAP WITH ONE RING IS USUALLY SEVERAL SEPARATE PATCHES, and they
        % must not be averaged into a single trace. A density segment can be far longer
        % than the stretch of dendrite that carries a given ROI's signal, so the ring
        % around it clips that footprint in two, three, four disconnected places -- on
        % opposite sides of the vessel, or at either end of it. Those patches sit at
        % different points along the neurite, so they are not at the same potential, and
        % one mean over all of them is a reference that matches none of them.
        %
        % So the overlap is split into connected components and each one above nvMinPx
        % becomes its OWN footprint and its own trace. Small crumbs are dropped rather
        % than merged: a handful of pixels is shot noise, and regressing shot noise out
        % of the movie removes variance everywhere for nothing.
        CCo = bwconncomp(ov);
        for p = 1:CCo.NumObjects
            pxo = CCo.PixelIdxList{p};
            if numel(pxo) < nvMinPx, continue, end
            w = zeros(nrow, ncol);
            w(pxo) = fj(pxo);                       % keep the footprint's own weighting
            nvW(:, end+1)  = w(:) / sum(w(:));                                %#ok<SAGROW>
            nvNm{end+1}    = sprintf('seg%d/ROI%d p%d (%d px)', sgi, j, p, numel(pxo)); %#ok<SAGROW>
            nvOfSeg(end+1) = sgi;                                             %#ok<SAGROW>
        end
    end
end
nNv0 = size(nvW, 2);
assert(nNv0 > 0, ['No ring-and-footprint overlap held a connected patch of at least ' ...
    'nvMinPx = %d px. Widen nvRingOut, lower nvFtprntThr or nvMinPx, or check that ' ...
    'ftprnt actually reaches the vessels in this field.'], nvMinPx);
fprintf('[R3-alt] ring %d-%d px outside each segment, minus densBW and bvBW: %d px\n', ...
    nvRingIn, nvRingOut, nnz(ringAll));
fprintf('         %d footprint PATCHES (connected, >= %d px) over %d segments\n', ...
    nNv0, nvMinPx, numel(unique(nvOfSeg)));

% ---- the traces: one chunked pass over mov_res, all footprints at once -------------
% Chunked because a full [npx x nTseg] copy of mov_res is several GB. The weights are
% sparse, so the product is cheap once a chunk is in memory.
nvTr = zeros(nNv0, nTseg);
nvChunk = 500;
for c0 = 1:nvChunk:nTseg
    ci = c0:min(c0+nvChunk-1, nTseg);
    nvTr(:, ci) = nvW' * double(reshape(mov_res(:,:,ci), npxA, numel(ci)));
end

% ---- smooth, drop duplicates, cap --------------------------------------------------
smNv = max(1, round(nvSmoothMs * fs / 1000));
nvTr = movmean(nvTr, smNv, 2);
% nvW is subset alongside the traces at every step below, so that after the cap it still
% holds EXACTLY the footprints that end up being regressed out -- that is what the first
% panel of R3-alt-ringReference draws, and a picture of the candidates rather than the
% survivors would be a picture of something the movie never saw.
okNv = all(isfinite(nvTr), 2) & std(nvTr, 0, 2) > 0;
nvTr = nvTr(okNv, :);  nvNm = nvNm(okNv);  nvOfSeg = nvOfSeg(okNv);
nvW  = nvW(:, okNv);
nvZ  = (nvTr - mean(nvTr,2)) ./ std(nvTr,0,2);

keepNv = true(1, size(nvZ,1));
for a = 2:size(nvZ,1)
    prev = find(keepNv(1:a-1));
    if isempty(prev), continue, end
    if any(abs(corr(nvZ(a,:)', nvZ(prev,:)')) > nvDedupR), keepNv(a) = false; end
end
nDup = nnz(~keepNv);
nvZ = nvZ(keepNv,:);  nvNm = nvNm(keepNv);  nvOfSeg = nvOfSeg(keepNv);
nvW = nvW(:, keepNv);
if nDup > 0
    fprintf('         %d dropped as duplicates (|r| > %.2f with one already kept)\n', ...
        nDup, nvDedupR);
end
if size(nvZ,1) > nvMaxTrace
    [~, ordNv] = sort(var(nvZ,0,2), 'descend');     % the ones carrying signal
    ordNv = sort(ordNv(1:nvMaxTrace));
    fprintf('         capped at nvMaxTrace = %d (of %d)\n', nvMaxTrace, size(nvZ,1));
    nvZ = nvZ(ordNv,:);  nvNm = nvNm(ordNv);  nvOfSeg = nvOfSeg(ordNv);
    nvW = nvW(:, ordNv);
end
nNv = size(nvZ, 1);
fprintf('         -> %d reference traces, smoothed %g ms (%d frames)\n', nNv, nvSmoothMs, smNv);

% =====================================================================================
% 3. mov_nv_alt, AND S_alt REBUILT WITH (7c)'S OWN SPATIAL FILTER
% =====================================================================================
mov_nv_alt = single(SeeResiduals(mov_res, nvZ));
vrA = var(double(mov_nv_alt), 0, 3);  v0A = var(double(mov_res), 0, 3);
onDendA = neuronMap > 0 & ~bvBW;
fprintf(['[R3-alt] variance left after the ring regression: %.1f%% overall, %.1f%% on the ' ...
         'neuron,\n         %.1f%% on bvMask, %.1f%% in the density mask  <- these last ' ...
         'two should stay HIGH\n'], ...
    100*mean(vrA(:))/mean(v0A(:)), ...
    100*mean(vrA(onDendA))/mean(v0A(onDendA)), ...
    100*mean(vrA(bvBW))/mean(v0A(bvBW)), ...
    100*mean(vrA(densBW))/mean(v0A(densBW)));
if mean(vrA(densBW))/mean(v0A(densBW)) < 0.8
    warning(['The ring reference removed %.0f%% of the variance INSIDE the density mask. ' ...
        'That is the vascular signal this cell exists to measure -- the ring is probably ' ...
        'touching vessel pixels. Raise nvRingIn.'], ...
        100*(1 - mean(vrA(densBW))/mean(v0A(densBW))));
end

S_alt = zeros(nrow, ncol, nTseg, 'single');
for c0 = 1:sbChunk:nTseg
    ci = c0:min(c0+sbChunk-1, nTseg);
    if strcmp(sbMethod,'fourier')
        S_alt(:,:,ci) = single(real(ifft2(fft2(mov_nv_alt(:,:,ci)) .* Hmask)));
    else
        S_alt(:,:,ci) = imfilter(mov_nv_alt(:,:,ci), kDoG, 'replicate');
    end
end
fprintf('[R3-alt] S_alt rebuilt with the (7c) ''%s'' filter\n', sbMethod);

% ---- figure: what the ring is, and what removing it did ---------------------------
get_figure('R3-alt-ringReference');
tlR = tiledlayout(1, 3, 'TileSpacing','compact','Padding','compact');
% Panel 1 shows ONLY the two things that matter: the density mask, and the voltage ROIs
% that are actually regressed out. Not the candidate rings, and not the whole footprint --
% those are intermediates, and drawing them makes the panel look like far more of the
% frame is being touched than really is. What is painted here is nvW AFTER the dedup and
% the cap, so it is exactly the set of weights that went into SeeResiduals.
axr = nexttile;
nvUsed = reshape(sum(nvW, 2), nrow, ncol);          % the union of the kept footprints
nvUsed = nvUsed / max(nvUsed(:), [], 'omitnan');
rgbR = zeros(nrow, ncol, 3);
rgbR(:,:,1) = 0.95*double(densBW);                  % red   = densBW
rgbR(:,:,2) = 0.95*min(nvUsed*3, 1);                % green = the voltage ROIs used
rgbR(:,:,3) = 0.95*min(nvUsed*3, 1);                % (cyan, so it cannot be read as red)
image(axr, rgbR); axis(axr,'image','off');
hold(axr,'on');
for a = 1:nNv                                       % one outline per ROI regressed out
    contour(axr, reshape(nvW(:,a) > 0, nrow, ncol), [0.5 0.5], 'w', 'LineWidth', 0.5);
end
title(axr, sprintf('red = densBW (%d px) | cyan = the %d voltage ROIs regressed out (%d px)', ...
    nnz(densBW), nNv, nnz(nvUsed > 0)));
axr2 = nexttile;
imagesc(axr2, vrA ./ max(v0A, eps), [0 1]); axis(axr2,'image','off');
colormap(axr2, turbo); colorbar(axr2);
hold(axr2,'on'); contour(axr2, densBW, [0.5 0.5], 'w', 'LineWidth', 0.8);
title(axr2, 'var(mov\_nv\_alt) / var(mov\_res)   (1 = untouched)');
axr3 = nexttile; hold(axr3,'on');
tR = 1:min(1500, nTseg);
for a = 1:min(6, nNv), plot(axr3, tR, nvZ(a,tR) + 4*(a-1), 'LineWidth', 0.7); end
box(axr3,'off'); xlim(axr3, [tR(1) tR(end)]); xlabel(axr3,'frame');
title(axr3, sprintf('first %d reference traces (%g ms smoothed)', min(6,nNv), nvSmoothMs));
title(tlR, sprintf(['f=%d  the voltage reference comes from the ring OUTSIDE each vessel: ' ...
    '%d traces kept of %d connected footprint patches (>= %d px each)'], ...
    f, nNv, nNv0, nvMinPx));

% =====================================================================================
% 4. PER SEGMENT: BAND FROM THE FLOW, THEN cpca / ica / disk
% =====================================================================================
SvvA    = reshape(-double(S_alt(:,:,trkFrames)), [], nRecF2);
fFloorA = max(oscFmin, 1);
gS      = ~badfr(trkFrames);

aSrch   = nan(nUseA, 2);   aBnd  = nan(nUseA, 2);
aStat   = nan(nUseA, 6);   % [nTrk speed spacing fTrk fCnt fPeak] -- (7g) reads 4 and 6
aPinned = false(1, nUseA);
aModeTr = cell(1, nUseA);  aModeU  = cell(1, nUseA);
aModeSc = cell(1, nUseA);  aModeNm = cell(1, nUseA);
aMode   = zeros(1, nUseA); % 1 = modes found, 2 = fell back to the segment mean
aY      = zeros(nUseA, nRecF2);
% Every method's modes are kept, not just the chosen one's, so that altCompare can build
% and FIT a full design per method rather than only scoring their leading modes.
aAll    = cell(1, nUseA);  % aAll{k}.(method) = struct('U',..,'tr',..,'sc',..,'nm',..)
cmpNM   = {'cpca','ica','disk','local'};
nCmp    = numel(cmpNM);
cmpR2   = nan(nUseA, nCmp);

rng(icaSeed);
for k = 1:nUseA
    sgi  = useSegA(k);
    dsel = detSegA == sgi;
    tq   = detT(dsel);

    % ---- the band, from the tracked flow (v/d) -- same logic as (7f) ---------------
    vMed = median(trkSpd(trkSegA == sgi));
    dd = [];
    for t = unique(round(tq))'
        pxy = detXY(dsel & round(detT)==t, :);
        if size(pxy,1) < 2, continue, end
        Dp = squareform(pdist(pxy)) + diag(inf(size(pxy,1),1));
        dd = [dd; min(Dp,[],2)]; %#ok<AGROW>
    end
    spacingNN = median(dd) * pxUm;
    fTrk = vMed / max(spacingNN, eps);
    fLim = [fFloorA, 0.45*fs];
    if ~isfinite(fTrk) || fTrk < fLim(1) || fTrk > fLim(2)
        fprintf(['   seg %3d: v/d = %.3g Hz outside %g-%g Hz (spacing %.2f um, %.0f um/s)' ...
                 ' -> full band\n'], sgi, fTrk, fLim(1), fLim(2), spacingNN, vMed);
        srchBand = fLim;
    else
        srchBand = [fTrk/oscSearch(2), fTrk*oscSearch(2)];
    end
    srchBand = [max(srchBand(1), fFloorA), min(srchBand(2), fs/2 - 1)];
    if srchBand(2) <= srchBand(1)*1.05, srchBand = [fFloorA, fs/2 - 1]; end
    aSrch(k,:) = srchBand;
    cnt  = accumarray(max(min(round(tq),nRecF2),1), 1, [nRecF2 1]);
    fCnt = median(movmean(cnt, smKa)) * fs;         % the BIASED estimate, for reference

    pxIdx = segPxA{sgi};
    Yseg  = SvvA(pxIdx, :);
    Yseg  = Yseg - mean(Yseg, 2);

    optM = struct('scoreMin', altScoreMin, 'goodMask', gS, 'segId', sgi, ...
        'binPix', segBinPix, 'minPxPCA', segMinPxPCA, 'pcFloor', icaPCfloor, ...
        'pxIdx', pxIdx, 'sz', [nrow ncol], 'spacing', diskSpacing, 'radius', diskR, ...
        'minPxDisk', diskMinPx);

    % ---- run the requested method, and the others too when comparing --------------
    if altCompare, want = cmpNM; else, want = {altDecomp}; end
    trAll = struct();  Uall = struct();  scAll = struct();  nmAll = struct();
    for w = 1:numel(want)
        dcp = decompOf(want{w});
        if strcmp(want{w}, 'local') && isfield(trAll, 'cpca')
            % identical modes to 'cpca' -- reuse rather than run the SVD twice
            Uw = Uall.cpca;  tw = trAll.cpca;  sw = scAll.cpca;
            nw = strrep(nmAll.cpca, 'cpca', 'local');
        else
            o = optM;
            switch dcp
                case 'cpca', o.maxMode = oscMaxPC;
                case 'ica',  o.maxMode = icaMaxIC;
                case 'disk', o.maxMode = diskMaxPer;
            end
            [Uw, tw, sw, nw] = get_segModes(Yseg, fs, srchBand, dcp, o);
            if isempty(tw)  % same mean fallback every method gets, so the comparison is fair
                tw = mean(Yseg, 1);  Uw = ones(numel(pxIdx),1)/numel(pxIdx);
                sw = get_oscScore(tw, fs, srchBand);
                nw = {sprintf('seg%d mean (%s found none)', sgi, want{w})};
            end
            if strcmp(want{w}, 'local'), nw = strrep(nw, 'cpca', 'local'); end
        end
        Uall.(want{w}) = Uw;  trAll.(want{w}) = tw;
        scAll.(want{w}) = sw; nmAll.(want{w}) = nw;
        aAll{k}.(want{w}) = struct('U', Uw, 'tr', tw, 'sc', sw(:), 'nm', {nw});
    end

    tr = trAll.(altDecomp);  Uk = Uall.(altDecomp);
    sc = scAll.(altDecomp);  nmk = nmAll.(altDecomp);
    aMode(k) = 1 + contains(nmk{1}, 'mean (');   % 2 = the fallback fired
    if aMode(k) == 2
        fprintf('   seg %3d: %4d px  %s found no mode above %.2f -> segment mean\n', ...
            sgi, numel(pxIdx), altDecomp, altScoreMin);
    end
    y = tr(1,:);

    % ---- the passband, from the leading mode's PSD peak ---------------------------
    yg = y;  yg(~gS) = NaN;
    [Py, Fy] = nanPSD(yg, fs, min(1024, 2^floor(log2(nRecF2/2))));
    Py = Py(:);  Fy = Fy(:);
    inS = Fy >= srchBand(1) & Fy <= srchBand(2);
    if ~any(inS & isfinite(Py)), inS = Fy > fFloorA & Fy < fs/2-1; end
    Pm = Py;  Pm(~inS | ~isfinite(Pm)) = -inf;
    [~, ipk] = max(Pm);  fPeak = Fy(ipk);
    aPinned(k) = oscFmin > 0 && fPeak <= fFloorA * 1.05;
    bnd = fPeak * [1-oscBandFrac, 1+oscBandFrac];
    bnd = [max(bnd(1), fFloorA), min(bnd(2), fs/2-1)];
    if ~all(isfinite(bnd)) || bnd(2) <= bnd(1)*1.05, bnd = srchBand; end
    aBnd(k,:)  = bnd;
    aStat(k,:) = [nnz(trkSegA==sgi) vMed spacingNN fTrk fCnt fPeak];

    aModeTr{k} = tr;  aModeU{k} = Uk;  aModeSc{k} = sc(:);  aModeNm{k} = nmk;
    aY(k,:) = y;

    fprintf('   seg %3d: %4d px  %-4s -> %d modes (score %s)  band %.0f-%.0f Hz%s\n', ...
        sgi, numel(pxIdx), altDecomp, size(tr,1), mat2str(round(sc(:)',2)), ...
        bnd(1), bnd(2), repmat('  <- PINNED', 1, aPinned(k)));

    % ---- the comparison: each method's leading mode against a common yardstick ----
    if altCompare
        ymeas = mean(Yseg, 1);  ymeas = ymeas - mean(ymeas(gS));
        for w = 1:nCmp
            tw = trAll.(cmpNM{w});
            if isempty(tw), continue, end
            yb_ = get_bandpass(tw(1,:), fs, bnd);  yb_(~isfinite(yb_)) = 0;
            z_  = hilbert(yb_);
            X1  = [real(z_(gS)); imag(z_(gS)); ones(1,nnz(gS))]';
            yv  = ymeas(gS)';
            cmpR2(k,w) = 1 - var(yv - X1*(X1\yv)) / max(var(yv), eps);
        end
    end
end
clear SvvA

% ---- drop the segments whose peak was pinned on the floor -------------------------
if oscDropPinned && any(aPinned)
    kp = ~aPinned;
    fprintf('[R3-alt] dropping %d segments pinned at the %g Hz floor: %s\n', ...
        nnz(aPinned), fFloorA, mat2str(useSegA(aPinned)));
    useSegA = useSegA(kp);  aModeTr = aModeTr(kp);  aModeU  = aModeU(kp);
    aModeSc = aModeSc(kp);  aModeNm = aModeNm(kp);  aBnd    = aBnd(kp,:);
    aSrch   = aSrch(kp,:);  aStat   = aStat(kp,:);  aY      = aY(kp,:);
    cmpR2   = cmpR2(kp,:);  aMode   = aMode(kp);    aPinned = aPinned(kp);
    aAll    = aAll(kp);
    nUseA   = numel(useSegA);
    assert(nUseA > 0, 'Every segment was pinned at oscFmin = %g Hz.', oscFmin);
end

if altCompare
    fprintf('\n[R3-alt] leading-mode R^2 against the segment MEAN signal, by method:\n');
    fprintf('     %6s%s   %s\n', 'seg', sprintf('%9s', cmpNM{:}), 'best');
    for k = 1:nUseA
        [~, bw] = max(cmpR2(k,:));
        fprintf('     %6d%s   %s\n', useSegA(k), sprintf('%9.3f', cmpR2(k,:)), cmpNM{bw});
    end
    fprintf('     %6s%s   (median)\n', 'all', ...
        sprintf('%9.3f', median(cmpR2, 1, 'omitnan')));
    fprintf(['     READ THIS WITH CARE: the yardstick is the segment MEAN, which is the one\n' ...
             '     signal none of these optimises -- a common reference, not the quantity\n' ...
             '     any of them maximises. And the mean is exactly what phase cancellation\n' ...
             '     destroys in a long segment, so on those it is a weak target for all of\n' ...
             '     them. A method can win here and still be the wrong design; the number\n' ...
             '     that decides is what the design REMOVES on the dendrite versus in the\n' ...
             '     density mask, printed below.\n']);
    fprintf(['     ''local'' scores the SAME as ''cpca'' here by construction -- identical\n' ...
             '     modes, and this table only looks at the mode. They differ in the fit,\n' ...
             '     which is what the removal table and the movies show.\n']);
end

% =====================================================================================
% 5. THE DESIGN: every kept mode gets its own quadrature pair
% =====================================================================================
% Built for EVERY method that was run, not only the chosen one. With altCompare on that
% is what lets each method be fitted and its corrected movie shown beside the others --
% comparing leading-mode R^2 alone would compare the modes, not the designs, and the
% design is what actually touches the movie.
zrowA = @(V) (V - mean(V,2)) ./ max(std(V,0,2), eps);
if altCompare, methRun = cmpNM; else, methRun = {altDecomp}; end
dsn = struct();
for w = 1:numel(methRun)
    mth = methRun{w};
    R_ = [];  N_ = {};  S_ = [];  M_ = [];  P_ = [];  np_ = 0;
    for k = 1:nUseA
        tr = aAll{k}.(mth).tr;  nmk = aAll{k}.(mth).nm;
        for m = 1:size(tr,1)
            yb = get_bandpass(tr(m,:), fs, aBnd(k,:));
            yb(~isfinite(yb)) = 0;
            z  = hilbert(yb);
            np_ = np_ + 1;
            R_ = [R_; zrowA(real(z)); zrowA(imag(z))];              %#ok<AGROW>
            N_ = [N_, {[nmk{m} ' | cos'], [nmk{m} ' | sin']}];      %#ok<AGROW>
            S_ = [S_, k, k];                                        %#ok<AGROW>
            M_ = [M_, m, m];                                        %#ok<AGROW>
            P_ = [P_, np_, np_];                                    %#ok<AGROW>
        end
    end
    dsn.(mth) = struct('reg',R_, 'nm',{N_}, 'segOf',S_, 'modeOf',M_, ...
                       'pairOf',P_, 'nPair',np_);
end
altReg = dsn.(altDecomp).reg;      altRegNm  = dsn.(altDecomp).nm;
altSegOf = dsn.(altDecomp).segOf;  altModeOf = dsn.(altDecomp).modeOf;
altPairOf = dsn.(altDecomp).pairOf;  nPair   = dsn.(altDecomp).nPair;
nRegA = size(altReg, 1);
assert(nRegA > 0, 'The alternate design came out empty.');
fprintf('\n[R3-alt] design: %d rows = %d quadrature pairs over %d segments (%s)\n', ...
    nRegA, nPair, nUseA, altDecomp);
if altCompare
    fprintf('         rows per method: %s\n', strjoin(arrayfun(@(w) ...
        sprintf('%s %d', cmpNM{w}, size(dsn.(cmpNM{w}).reg,1)), 1:nCmp, 'uni', 0), ' | '));
end

% ---- collinearity: did the method split one source across several modes? ----------
% This is the specific way ICA fails here, so it is MEASURED rather than assumed away.
% ICA has no variance ordering and no orthogonality constraint on the time courses, so a
% single travelling wave can come back as two ICs 90 deg apart -- four rows for one
% physical source. The fit stays well defined (mldivide handles the rank deficiency) but
% the individual coefficients stop meaning anything.
cndA = cond(altReg(:,gF) * altReg(:,gF)');
Rall = corr(altReg(:,gF)');
Rall(logical(eye(nRegA))) = 0;
[mxR, mxI] = max(abs(Rall(:)));
[r1, r2]   = ind2sub([nRegA nRegA], mxI);
fprintf('         design cond = %.3g | largest |corr| between two rows = %.3f\n', cndA, mxR);
if mxR > 0.9
    fprintf('           %s\n           %s\n', altRegNm{r1}, altRegNm{r2});
    warning(['Two rows are %.0f%% correlated -- with altDecomp = ''%s'' that usually ' ...
        'means one source was split across two modes. The fit is still well defined, ' ...
        'but do not read the individual coefficients.'], 100*mxR, altDecomp);
end

% =====================================================================================
% 6. SPATIAL MAPS, THE VOLTAGE VETO, AND WHAT THE DESIGN REMOVES
% =====================================================================================
% Every method's covariance map is taken in the SAME pass. MvA is [nPix x nGoodFrames]
% and costs the better part of a gigabyte, so it is materialised once, all the products
% that need it are taken, and it is dropped -- rather than rebuilding it per method later
% when the side-by-side movie needs the other two.
MvA     = single(tovec(mov_res(:,:,trkFrames)));
MvA     = MvA(:, gF);
pxMeanA = mean(MvA, 2);          % kept: it is column 1 of every coefficient matrix below
MvA     = MvA - pxMeanA;
covByM = struct();
for w = 1:numel(methRun)
    Zw = dsn.(methRun{w}).reg(:, gF);
    Zw = Zw - mean(Zw, 2);
    covByM.(methRun{w}) = double(MvA * single(Zw)') / max(nnz(gF)-1, 1);
end
clear Zw
covA = covByM.(altDecomp);

% Which pixels each segment is allowed to correct, for the 'local' fit. Built once and
% reused by the chosen design and by every method in the comparison.
% Segments are connected components and so disjoint; localDilate can make them overlap,
% and a pixel claimed twice is given to the segment that reached it first -- the
% alternative, fitting it with both segments' rows at once, is a global fit in miniature
% and would reintroduce exactly the reach this method exists to remove.
locPx = cell(1, nUseA);  claimed = false(npxA, 1);  nOverlap = 0;
for k = 1:nUseA
    pk = segPxA{useSegA(k)};
    if localDilate > 0
        mk = false(nrow, ncol);  mk(pk) = true;
        pk = find(imdilate(mk, strel('disk', localDilate)));
    end
    dup = claimed(pk);
    nOverlap = nOverlap + nnz(dup);
    pk = pk(~dup);
    claimed(pk) = true;
    locPx{k} = pk;
end
fprintf('[R3-alt] ''local'' correction covers %d px (%.1f%% of the frame), dilate %d px%s\n', ...
    nnz(claimed), 100*nnz(claimed)/npxA, localDilate, ...
    repmat(sprintf(', %d px claimed twice', nOverlap), 1, nOverlap > 0));

[scoreA, f0A] = get_oscScore(altReg, fs, [1 fs/2-1]);

% the voltage reference: reuse whichever (7f) settled on, else the footprint
if exist('voltRef','var') && ~isempty(voltRef) && isequal(size(voltRef), [nrow ncol])
    vRefA = voltRef;  vRefNmA = oscVoltRefUse;
else
    vRefA = neuronMap ./ max(neuronMap(:), [], 'omitnan');  vRefNmA = 'ftprnt';
end
vMaskA = (neuronMap > 0) | bvBW | densBW;
rVoltA = get_voltCorr(covA, vRefA, vMaskA);

% HARD VETO, propagated across a quadrature pair. The pair spans every phase between the
% two rows, so if either half's map looks like the neuron the pair can reconstruct it --
% keeping the partner would let the voltage back in through the side its own map does not
% show. The score rule moves the OTHER way: either half passing carries the pair.
voltKillA = abs(rVoltA) >= altVoltCut;
scoreOKA  = scoreA >= altScoreMin;
for p = 1:nPair
    ii = altPairOf == p;
    scoreOKA(ii)  = any(scoreOKA(ii));
    voltKillA(ii) = any(voltKillA(ii));
end
if altDropBelowScore, keepA = scoreOKA & ~voltKillA; else, keepA = ~voltKillA; end

fprintf('\n[R3-alt] every row (scoreMin %.2f, voltage veto |r| >= %.2f vs ''%s''):\n', ...
    altScoreMin, altVoltCut, vRefNmA);
fprintf('     %-32s %8s %8s %10s  %s\n', 'name', 'oscScore', 'r(volt)', 'peak|cov|', 'seg');
[~, ordA] = sort(scoreA, 'ascend');
for r = ordA'
    flg = '';
    if voltKillA(r),  flg = sprintf('   <- VETOED (|r| >= %.2f)', altVoltCut);
    elseif ~keepA(r), flg = '   <- dropped';
    end
    fprintf('     %-32s %8.3f %+8.3f %10.4g  %3d%s\n', altRegNm{r}, scoreA(r), ...
        rVoltA(r), max(abs(covA(:,r))), useSegA(altSegOf(r)), flg);
end
altVetoed = struct('name', {altRegNm(voltKillA)}, 'rVolt', rVoltA(voltKillA), ...
    'oscScore', scoreA(voltKillA), 'cut', altVoltCut);
if any(voltKillA)
    fprintf('[R3-alt] %d rows vetoed by the voltage map, %d more below scoreMin\n', ...
        nnz(voltKillA), nnz(~keepA & ~voltKillA));
end
altKeptPair = unique(altPairOf(keepA));
if any(~keepA)
    altReg    = altReg(keepA,:);     altRegNm  = altRegNm(keepA);
    altSegOf  = altSegOf(keepA);     altModeOf = altModeOf(keepA);
    altPairOf = altPairOf(keepA);    scoreA    = scoreA(keepA);
    f0A       = f0A(keepA);          rVoltA    = rVoltA(keepA);
    covA      = covA(:, keepA);
    assert(~isempty(altReg), ['Every row was rejected: %d vetoed by the voltage map ' ...
        '(altVoltCut = %.2f), the rest below %.2f.'], nnz(voltKillA), altVoltCut, altScoreMin);
    nRegA = size(altReg,1);
end
fprintf('[R3-alt] %d rows survive\n', nRegA);

% ---- what does it take out of mov_res? --------------------------------------------
% Fitted on mov_res -- the RAW motion residual, voltage still in it -- because that is the
% honest final test. The number to watch is the DENDRITE row: it should stay small while
% the density-mask and bvMask rows are large.
%
% GLOBAL vs LOCAL. The global fit gives every pixel every row, so a segment's oscillator
% can remove variance on the far side of the frame if the noise there happens to correlate
% with it -- harmless per pixel, but it is what puts the smooth arc under the % removed
% map. The local fit gives pixel p only the rows of the segment that contains it, and
% gives a pixel outside every segment nothing at all.
%
% The local result is expressed as the SAME [nPix x (1+nReg)] coefficient matrix, just
% with zeros outside each pixel's own block. That is not a convenience: it means every
% downstream formula -- varFit = b'*Sf*b, the fitted movie Cfv*[1; design], the export to
% (7g) -- keeps working with no special case, and the two fits stay directly comparable
% because they are the same object. Column 1 is the per-pixel mean, matching what
% SeeResiduals returns, so `Cfv * [ones; centred design]` reconstructs either way.
if strcmp(altDecomp, 'local')
    CfvA = zeros(npxA, 1 + nRegA);
    CfvA(:, 1) = pxMeanA;                       % untouched pixels keep only their mean
    ZcA = altReg(:, gF) - mean(altReg(:, gF), 2);
    nFitPx = 0;
    for k = 1:nUseA
        rws = find(altSegOf == k);
        pk  = locPx{k};
        if isempty(rws) || isempty(pk), continue, end
        % centred data on centred regressors: the intercept is then exactly the pixel
        % mean, which is already in column 1, so no ones-row is needed here
        CfvA(pk, 1 + rws) = double(MvA(pk, :)) / ZcA(rws, :);
        nFitPx = nFitPx + numel(pk);
    end
    fprintf('[R3-alt] LOCAL fit: %d px fitted inside their own segment, %d px untouched\n', ...
        nFitPx, npxA - nFitPx);
    clear ZcA
else
    [~, CfA] = SeeResiduals(mov_res(:,:,trkFrames(gF)), altReg(:,gF));
    CfvA     = tovec(CfA);
end
BfA     = CfvA(:, 2:end);
varInA  = var(mov_res(:,:,trkFrames(gF)), 0, 3);
SfA     = cov(altReg(:, gF)');
varFitA = reshape(sum((BfA * SfA) .* BfA, 2), nrow, ncol);
dropA   = varFitA ./ max(varInA, eps);

% THE CHANCE FLOOR. Any design of nReg rows fitted to nGood frames removes about
% nReg/nGood of the variance of PURE NOISE, with no signal involved -- least squares
% always finds some projection. So no percentage below this line means anything, and a
% map that looks uniformly "a bit corrected" everywhere is usually just showing it.
% Printed first because every number under it has to be read against it.
% Rank, not row count: collinear rows (a quadrature pair split across two modes, two
% overlapping disks) do not each buy a degree of freedom, so nReg would overstate the
% floor for exactly the designs most likely to have one.
rankA   = rank(altReg(:, gF));
chanceA = rankA / max(nnz(gF), 1);
bgA     = ~(neuronMap > 0) & ~bvBW & ~densBW;      % off the cell and off every vessel
bgMean  = mean(dropA(bgA));
fprintf('\n[R3-alt] chance floor: rank %d of %d rows / %d good frames = %.2f%% from noise alone\n', ...
    rankA, nRegA, nnz(gF), 100*chanceA);
fprintf('     %-26s %8s %9s %s\n', 'region', 'removed%', 'x chance', 'x background');
rgA = {densBW, 'density mask'; bvBW, 'bvMask'; onDendA, 'DENDRITE off-vessel'; ...
       bgA, 'background (no cell)'};
for q = 1:4
    mq = rgA{q,1};
    if ~any(mq(:)), continue, end
    fprintf('     %-26s %8.2f %8.1fx %8.2fx\n', rgA{q,2}, 100*mean(dropA(mq)), ...
        mean(dropA(mq))/max(chanceA, eps), mean(dropA(mq))/max(bgMean, eps));
end
% READING THE TABLE. The two ratio columns answer two different questions and the second
% is the one that matters.
%
%   x chance   compares against white noise. It is an idealisation: the rows are NARROW-
%              BAND, so a pixel whose noise happens to carry power in the design's
%              passbands loses more than rank/nFrames even with no vascular signal in it.
%              Background above 1x is therefore normal, not an alarm.
%
%   x background  compares against a pixel with no cell and no vessel in it -- an
%              EMPIRICAL null that already contains whatever global component the movie
%              has: scattered and out-of-focus light from the sample, illumination drift,
%              residual motion. This is the selectivity that matters.
%
% And the comparison that settles whether voltage is being eaten is DENDRITE vs
% BACKGROUND, not dendrite vs chance. Voltage lives on the cell and not in the
% background, so a design eating voltage puts the dendrite ABOVE the background. A
% dendrite that merely tracks the background is being touched by the global component
% and nothing else -- whatever that costs, it costs the empty frame equally.
selD = mean(dropA(densBW)) / max(bgMean, eps);
selN = mean(dropA(onDendA)) / max(bgMean, eps);
fprintf('     -> vessel selectivity %.2fx background | dendrite %.2fx background\n', selD, selN);
if selN > 1.25
    warning(['The dendrite sits %.2fx the background while the density mask sits %.2fx. ' ...
        'The design is reaching the cell for its own sake, not just through the global ' ...
        'component.'], selN, selD);
elseif selD > 1.5
    fprintf(['     dendrite ~= background, so no voltage is being eaten; the excess over ' ...
             'chance is\n     global (scatter / out-of-focus signal) and costs the cell ' ...
             'and the empty frame alike.\n']);
end
if mean(dropA(onDendA)) > 0.5*mean(dropA(bvBW))
    warning(['%.1f%% removed on the dendrite vs %.1f%% on bvMask: the design is probably ' ...
        'taking neural signal.'], 100*mean(dropA(onDendA)), 100*mean(dropA(bvBW)));
end

% =====================================================================================
% 6b. THE S_alt SUBTRACTION: band-pass S_alt inside each segment, subtract from mov_res
% =====================================================================================
% No design matrix and no shared waveform. Each covered pixel gets its own band-limited
% trace out of S_alt, one scalar gain, and that is what comes off mov_res. See the block
% at altSfilt for why this is neither a regression nor a notch.
if altSfilt
    sPx    = vertcat(locPx{:});                 % the same pixels 'local' is allowed to touch
    sSegOf = zeros(numel(sPx), 1);
    p0 = 0;
    for k = 1:nUseA
        sSegOf(p0 + (1:numel(locPx{k}))) = k;
        p0 = p0 + numel(locPx{k});
    end
    nCovS = numel(sPx);
    assert(nCovS > 0, 'No segment pixels to subtract from.');

    % pull S_alt at those pixels in one chunked pass -- a full [npx x nRecF2] copy of
    % S_alt would be most of a gigabyte, and only the segment rows are ever needed
    Sseg = zeros(nCovS, nRecF2, 'single');
    for c0 = 1:500:nRecF2
        ci = c0:min(c0+499, nRecF2);
        Bc = reshape(S_alt(:,:,trkFrames(ci)), npxA, numel(ci));
        Sseg(:, ci) = Bc(sPx, :);
    end
    clear Bc

    % each pixel band-passed to ITS OWN segment's band
    for k = 1:nUseA
        ii = sSegOf == k;
        if ~any(ii), continue, end
        tmp = get_bandpass(double(Sseg(ii,:)), fs, aBnd(k,:));
        tmp(~isfinite(tmp)) = 0;
        Sseg(ii,:) = single(tmp);
    end

    % one scalar per pixel. Centred on the good frames, so the gain is a pure projection
    % and carries no offset; the DC stays in pxMeanA exactly as for the other methods.
    SgS = double(Sseg(:, gF));  SgS = SgS - mean(SgS, 2);
    McS = double(MvA(sPx, :));                  % mov_res, already mean-removed
    if strcmp(sfiltGain, 'fit')
        gainS = sum(McS .* SgS, 2) ./ max(sum(SgS.^2, 2), eps);
    else
        gainS = ones(nCovS, 1);
    end
    varFitS = var(gainS .* SgS, 0, 2);
    dropS   = zeros(nrow, ncol);
    dropS(sPx) = varFitS ./ max(varInA(sPx), eps);

    fprintf('\n[R3-alt] S_alt subtraction over %d px, gain ''%s'' (median %.3f, IQR %.3f-%.3f)\n', ...
        nCovS, sfiltGain, median(gainS), prctile(gainS,25), prctile(gainS,75));
    fprintf('     %-26s %8s %s\n', 'region', 'removed%', 'x background');
    for q = 1:4
        mq = rgA{q,1};
        if ~any(mq(:)), continue, end
        fprintf('     %-26s %8.2f %10.2fx\n', rgA{q,2}, 100*mean(dropS(mq)), ...
            mean(dropS(mq))/max(mean(dropS(bgA)), eps));
    end
    fprintf('     background is 0 by construction here too -- nothing outside a segment is touched\n');
    % HOW MUCH OF S_alt AT A PIXEL IS THAT PIXEL'S OWN NOISE. The DoG is a local average,
    % so S_alt(p) contains mov_res(p) with the kernel's centre weight. Fitting a gain
    % between them therefore fits a little noise to itself and biases the gain up. The
    % weight is small for a 6-30 um DoG, but it is printed rather than assumed, because it
    % is the one way this method can inflate its own removal number.
    if strcmp(sbMethod, 'dog') && exist('kDoG','var')
        w0 = kDoG(ceil(size(kDoG,1)/2), ceil(size(kDoG,2)/2)) / max(sum(abs(kDoG(:))), eps);
        fprintf(['     DoG centre weight %.3f: that fraction of S_alt at a pixel is the ' ...
                 'pixel itself,\n     so the fitted gain is biased up by roughly that ' ...
                 'much through shared noise\n'], w0);
    end
    if mean(dropS(onDendA)) > 0.5*mean(dropS(bvBW))
        warning(['S_alt subtraction removes %.1f%% on the dendrite vs %.1f%% on bvMask. ' ...
            'The DoG responds to the dendrite EDGES ((7c) says so), and this is what that ' ...
            'costs -- it is the failure mode to check for this method.'], ...
            100*mean(dropS(onDendA)), 100*mean(dropS(bvBW)));
    end
    clear McS
end

% ---- the disk tiling, when that is what is being used -----------------------------
% Worth drawing because the disk method has no spatial model at all: the ROIs ARE the
% result, and whether they cover the vessel sensibly is something only a picture answers.
% Overlap between neighbours is expected (diskR > diskSpacing/2 is deliberate) -- what to
% look for instead is a segment whose disks miss its middle, or one disk stretched across
% two vessels.
if strcmp(altDecomp, 'disk')
    get_figure('R3-alt-diskROIs');
    axd = axes; hold(axd, 'on');
    bg = 0.85*ones(nrow, ncol, 3);
    bg(:,:,1) = bg(:,:,1) - 0.35*double(densBW);      % densBW as a light grey-blue bed
    bg(:,:,2) = bg(:,:,2) - 0.15*double(densBW);
    image(axd, bg); axis(axd, 'image', 'off'); set(axd, 'YDir', 'reverse');
    nDk = sum(cellfun(@(u) size(u,2), aModeU));
    ck  = turbo(max(nDk, 2));
    q = 0;  lblD = {};
    for k = 1:nUseA
        Uk = aModeU{k};
        for m = 1:size(Uk,2)
            q = q + 1;
            dm = false(nrow, ncol);
            dm(segPxA{useSegA(k)}) = Uk(:,m) > 0;
            contour(axd, dm, [0.5 0.5], 'LineWidth', 1.4, 'LineColor', ck(q,:));
            [rr_, cc_] = find(dm);
            text(axd, mean(cc_), mean(rr_), sprintf('%d', m), 'Color', ck(q,:), ...
                'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        end
        lblD{end+1} = sprintf('s%d:%d', useSegA(k), size(Uk,2)); %#ok<AGROW>
    end
    contour(axd, bvBW, [0.5 0.5], 'k', 'LineWidth', 0.8);
    title(axd, sprintf(['f=%d  (7f-ALT) the %d disk ROIs actually used (r = %g px, ' ...
        'spacing %g px, >= %d px inside)\ndisks per segment: %s   |   grey = densBW, ' ...
        'black = bvMask'], f, nDk, diskR, diskSpacing, diskMinPx, strjoin(lblD, ', ')));
end

% ---- figures -----------------------------------------------------------------------
get_figure('R3-alt-modes');
nShowA = min(altSegShow, nUseA);
nMapA  = min(altModeMap, max(cellfun(@(u) size(u,2), aModeU)));
tlA = tiledlayout(nShowA + nMapA, 1, 'TileSpacing','compact','Padding','compact');
for j = 1:nMapA
    axj = nexttile;  eM = nan(nrow, ncol);  aM = zeros(nrow, ncol);  lb = {};
    for k = 1:nUseA
        Uk = aModeU{k};
        if size(Uk,2) < j, continue, end
        u = Uk(:,j);  u = u / max(abs(u(:)));
        if ~isreal(u)
            eM(segPxA{useSegA(k)}) = angle(u);   aM(segPxA{useSegA(k)}) = abs(u);
        else
            if sum(u) < 0, u = -u; end           % the sign of a real mode is arbitrary
            eM(segPxA{useSegA(k)}) = u;          aM(segPxA{useSegA(k)}) = abs(u);
        end
        lb{end+1} = sprintf('s%d:%.2f', useSegA(k), aModeSc{k}(j)); %#ok<AGROW>
    end
    hM = imagesc(axj, eM);
    set(hM, 'AlphaData', (aM ./ max(aM(:))) .* ~isnan(eM));
    if ismember(altDecomp, {'cpca','local'})   % 'local' modes ARE complex PCA's
        clim(axj, [-pi pi]);  colormap(axj, hsv);
    else
        clim(axj, [-1 1]);    colormap(axj, cmapCorr);
    end
    colorbar(axj); axis(axj,'image','off');
    hold(axj,'on'); contour(axj, bvBW, [0.5 0.5], 'k', 'LineWidth', 0.8);
    title(axj, sprintf('%s mode %d  (score: %s)', altDecomp, j, strjoin(lb, ', ')));
end
tA = 1:min(1200, nRecF2);
for k = 1:nShowA
    axk = nexttile; hold(axk,'on');
    yy = aY(k,tA);  yy = (yy - mean(yy))/max(std(yy),eps);
    yb = get_bandpass(aY(k,:), fs, aBnd(k,:));  yb(~isfinite(yb)) = 0;
    zz = hilbert(yb);
    cc = real(zz(tA));  cc = (cc - mean(cc))/max(std(cc),eps);
    plot(axk, tA, yy, 'k-', 'LineWidth', 0.8, 'DisplayName','leading mode');
    plot(axk, tA, cc, 'r-', 'LineWidth', 1.0, 'DisplayName','A(t)cos\phi(t)');
    box(axk,'off'); xlim(axk, [tA(1) tA(end)]);
    if k == 1, legend(axk, 'Location','northeast'); end
    ylabel(axk, sprintf('seg %d', useSegA(k)));
    title(axk, sprintf('%.0f um/s, v/d %.0f Hz, PSD peak %.0f Hz, band %.0f-%.0f Hz', ...
        aStat(k,2), aStat(k,4), aStat(k,6), aBnd(k,1), aBnd(k,2)));
    if k == nShowA, xlabel(axk, 'frame'); end
end
title(tlA, sprintf(['f=%d  (7f-ALT) %s on the ring-corrected movie: %d rows, removes ' ...
    '%.1f%% in the density mask / %.1f%% on the dendrite'], f, altDecomp, nRegA, ...
    100*mean(dropA(densBW)), 100*mean(dropA(onDendA))));

statA = struct('f0Hz', f0A(:), 'oscScore', scoreA(:), 'score', scoreA(:), ...
    'rVolt', rVoltA(:), 'ampCov', max(abs(covA),[],1)', ...
    'bandLo', aBnd(altSegOf(:),1), 'bandHi', aBnd(altSegOf(:),2), ...
    'segment', useSegA(altSegOf(:))');
browse_regressors(altReg, covA, [nrow ncol], altRegNm, struct( ...
    'score', scoreA, 'scoreMin', altScoreMin, 'f0', f0A, 'seg', altSegOf, ...
    'segMask', {cellfun(@(kk) segPxA{kk}, num2cell(useSegA), 'uni', 0)}, ...
    'refMask', bvBW, 'cmap', cmapCorr, 'tWin', 1:min(1000, nRecF2), ...
    'stats', statA, 'statXY', {{'f0Hz','score'}}, 'statC', 'rVolt', ...
    'figName', 'R3-alt-flowRegressors'));

% ---- the corrected movie, against the movie it was corrected FROM ------------------
%   [ mov_res (uncorrected) | mov_res - fitted (corrected) ]
%   [ fitted flow           | S_alt                        ]
%
% The top row is the comparison asked for, side by side on the SAME colour scale -- a
% correction shown on its own scale can look like anything. The bottom row is what makes
% the top row readable: the fitted component is what was subtracted, and S_alt is where
% the design came from, so a moving cell visible in S_alt should be visible in the fitted
% panel and gone from the corrected one.
%
% WHAT TO LOOK FOR, in this order:
%   1. the vessels stop flickering between the two top panels;
%   2. the DENDRITE looks the same in both -- if it dims or its transients soften, the
%      design is eating voltage, and the dendrite number printed above already said so;
%   3. the fitted panel is dark everywhere except on vessels. Structure in it away from
%      the density mask is the design reaching somewhere it should not.
%
% Fitted on the good tracked frames, applied to the first altMovFrm of them. Spike-
% positive convention (-mov_res) throughout, matching every other movie in this script.
%
% densBW is outlined on EVERY panel, tiled to match the layout. It is drawn as an overlay
% after moviefixsc rather than burnt into the pixels: moviesc clears the figure on entry
% and then animates by writing imh.CData only, so anything added afterwards survives every
% frame -- and the pixels underneath stay their real values instead of being overwritten
% by a marker. Its colormap is gray(2^11), so a coloured line is guaranteed to read
% against both the dark and the bright ends of the movie.
nQA   = min(altMovFrm, nRecF2);
IallA = [ones(1,nQA); altReg(:,1:nQA) - mean(altReg(:,gF), 2)];
fitA  = reshape(single(CfvA) * single(IallA), nrow, ncol, nQA);
aI = single(-mov_res(:,:,trkFrames(1:nQA)));    % uncorrected
aF = -fitA;                                     % the fitted vascular component
aC = aI - aF;                                   % corrected
aS = single(-S_alt(:,:,trkFrames(1:nQA)));      % what the design was built on

sA = std(aS(:)) / max(std(aI(:)), eps);         % one display gain for the -mov_res panels
aI = aI*sA;  aF = aF*sA;  aC = aC*sA;
climA = double(prctile(aS(:), [1 99]));
sepA  = single(climA(2));
rowSepA = sepA * ones(altMovSep, 2*ncol + altMovSep, nQA, 'single');
colSepA = sepA * ones(nrow, altMovSep, nQA, 'single');
movA = cat(1, cat(2, aI, colSepA, aC), rowSepA, cat(2, aF, colSepA, aS));

% densBW replicated into all four quadrants of the tile, in the movie's own coordinates
r2A = nrow + altMovSep;   c2A = ncol + altMovSep;
densTile = false(2*nrow + altMovSep, 2*ncol + altMovSep);
densTile(1:nrow,      1:ncol)      = densBW;
densTile(1:nrow,      c2A+(1:ncol)) = densBW;
densTile(r2A+(1:nrow), 1:ncol)      = densBW;
densTile(r2A+(1:nrow), c2A+(1:ncol)) = densBW;

get_figure('R3-alt-correctedMovie');
moviefixsc(movA, climA);
if altMovMark
    hold(gca, 'on');
    contour(gca, double(densTile), [0.5 0.5], 'LineWidth', 0.75, 'LineColor', [0 1 1]);
end
title(sprintf(['f=%d  (7f-ALT/%s)  [uncorrected | CORRECTED ; fitted flow | S\\_alt]   ' ...
    'cyan = densBW\n%d rows remove %.1f%% in the density mask, %.1f%% on bvMask, ' ...
    '%.1f%% on the dendrite'], ...
    f, altDecomp, nRegA, 100*mean(dropA(densBW)), 100*mean(dropA(bvBW)), ...
    100*mean(dropA(onDendA))));
fprintf(['[R3-alt] R3-alt-correctedMovie: [uncorrected | corrected ; fitted | S_alt], ' ...
         '%d frames, shared colour scale\n'], nQA);

% ---- altCompare: the corrected movie from EVERY method, side by side --------------
%   [ uncorrected | cpca  | ica   ]
%   [ disk        | local | S_alt ]
%
% Each method's FULL design is fitted here, not just its leading mode. That distinction
% is the point: the R^2 table above compares the modes, and a method can win there while
% its design does something quite different to the movie -- more rows, more collinearity,
% more reach onto the dendrite. This compares what actually gets subtracted.
%
% All panels share ONE colour scale, taken from the uncorrected movie. A per-panel scale
% would rescale each correction to fill the range and make every method look equally
% effective.
%
% The same score threshold and voltage veto are applied to each design, so the comparison
% is between the METHODS and not between four different admission policies. 'local' also
% gets the same rows as 'cpca' -- only its FIT is confined -- so the pair of panels
% isolates the effect of the confinement on its own, with the decomposition held fixed.
if altCompare
    fprintf('\n[R3-alt] fitting all %d designs for the side-by-side movie:\n', nCmp);
    fprintf('     %-6s %6s %8s %10s %10s %10s %12s\n', 'method', 'rows', 'vetoed', ...
        'densBW %', 'bvMask %', 'dendrite %', 'background %');
    cmpPanel = cell(1,nCmp);  cmpDrop = nan(nCmp,4);  cmpRows = zeros(1,nCmp);
    for w = 1:nCmp
        mth = cmpNM{w};
        Rw  = dsn.(mth).reg;   Pw = dsn.(mth).pairOf;   Sw = dsn.(mth).segOf;
        % same admission rules as the chosen design, applied to this one
        scW  = get_oscScore(Rw, fs, [1 fs/2-1]);
        rvW  = get_voltCorr(covByM.(mth), vRefA, vMaskA);
        kilW = abs(rvW) >= altVoltCut;
        okW  = scW >= altScoreMin;
        for p = unique(Pw)
            ii = Pw == p;
            okW(ii)  = any(okW(ii));
            kilW(ii) = any(kilW(ii));
        end
        if altDropBelowScore, kpW = okW & ~kilW; else, kpW = ~kilW; end
        if ~any(kpW)
            fprintf('     %-6s  every row rejected -- panel left as the uncorrected movie\n', mth);
            cmpPanel{w} = [];  cmpRows(w) = 0;  continue
        end
        Rw = Rw(kpW, :);  Sw = Sw(kpW);
        if strcmp(mth, 'local')
            CfvW = zeros(npxA, 1 + size(Rw,1));
            CfvW(:, 1) = pxMeanA;
            ZcW = Rw(:, gF) - mean(Rw(:, gF), 2);
            for k = 1:nUseA
                rws = find(Sw == k);  pk = locPx{k};
                if isempty(rws) || isempty(pk), continue, end
                CfvW(pk, 1 + rws) = double(MvA(pk, :)) / ZcW(rws, :);
            end
            clear ZcW
        else
            [~, CfW] = SeeResiduals(mov_res(:,:,trkFrames(gF)), Rw(:,gF));
            CfvW = tovec(CfW);
        end
        BfW  = CfvW(:, 2:end);
        SfW  = cov(Rw(:, gF)');
        vFitW = reshape(sum((BfW * SfW) .* BfW, 2), nrow, ncol);
        dW    = vFitW ./ max(varInA, eps);
        cmpDrop(w,:) = [mean(dW(densBW)), mean(dW(bvBW)), mean(dW(onDendA)), mean(dW(bgA))];
        cmpRows(w)   = size(Rw,1);
        IallW = [ones(1,nQA); Rw(:,1:nQA) - mean(Rw(:,gF), 2)];
        fitW  = reshape(single(CfvW) * single(IallW), nrow, ncol, nQA);
        cmpPanel{w} = single(-mov_res(:,:,trkFrames(1:nQA))) + fitW;   % corrected, spike-positive
        fprintf('     %-6s %6d %8d %10.2f %10.2f %10.2f %12.2f\n', mth, size(Rw,1), ...
            nnz(kilW), 100*cmpDrop(w,1), 100*cmpDrop(w,2), 100*cmpDrop(w,3), ...
            100*cmpDrop(w,4));
        clear CfW CfvW BfW SfW vFitW IallW fitW
    end
    % ---- the S_alt subtraction as one more panel, built the same way ---------------
    % It has no design and no rows, so it cannot go through the loop above; the corrected
    % movie is formed directly from the per-pixel band-limited traces and their gains.
    sfPanel = [];
    if altSfilt
        remS = zeros(nrow*ncol, nQA, 'single');
        remS(sPx, :) = single(gainS .* double(Sseg(:, 1:nQA) - mean(Sseg(:, gF), 2)));
        sfPanel = single(-mov_res(:,:,trkFrames(1:nQA))) + reshape(remS, nrow, ncol, nQA);
        fprintf('     %-6s %6s %8s %10.2f %10.2f %10.2f %12.2f\n', 'sfilt', ...
            sprintf('%dpx', nCovS), '-', 100*mean(dropS(densBW)), 100*mean(dropS(bvBW)), ...
            100*mean(dropS(onDendA)), 100*mean(dropS(bgA)));
        clear remS
    end
    fprintf(['     the columns that decide are DENDRITE and BACKGROUND. More in densBW ' ...
             'bought with\n     more on the dendrite is not a better method, only a ' ...
             'harder fit; and background\n     above ~0 is reach the correction should ' ...
             'not have had. ''local'' and ''sfilt''\n     are 0.00 there by construction ' ...
             '-- for them it is a guarantee, not a result.\n']);

    % ---- tile: [uncorrected | cpca | ica ; disk | local | sfilt] -------------------
    uncM  = single(-mov_res(:,:,trkFrames(1:nQA)));
    gainC = std(aS(:)) / max(std(uncM(:)), eps);       % same gain the panel above used
    for w = 1:nCmp
        if isempty(cmpPanel{w}), cmpPanel{w} = uncM; end
        cmpPanel{w} = cmpPanel{w} * gainC;
    end
    uncM = uncM * gainC;
    if ~isempty(sfPanel)
        panC = [{uncM}, cmpPanel, {sfPanel * gainC}];
        labC = [{'uncorrected'}, cmpNM, {'sfilt'}];
    else
        panC = [{uncM}, cmpPanel, {aS}];
        labC = [{'uncorrected'}, cmpNM, {'S\_alt'}];
    end
    tRw = 2;  tCl = ceil(numel(panC) / tRw);
    sepC = single(climA(2));
    colSepC = sepC * ones(nrow, altMovSep, nQA, 'single');
    rowSepC = sepC * ones(altMovSep, tCl*ncol + (tCl-1)*altMovSep, nQA, 'single');
    movC = [];
    for rr = 1:tRw
        blk = [];
        for cc = 1:tCl
            iq = (rr-1)*tCl + cc;
            if iq <= numel(panC), P = panC{iq}; else, P = zeros(nrow,ncol,nQA,'single'); end
            if cc == 1, blk = P; else, blk = cat(2, blk, colSepC, P); end
        end
        if rr == 1, movC = blk; else, movC = cat(1, movC, rowSepC, blk); end
    end
    % densBW outline, retiled for THIS grid (the 2x2 densTile no longer fits)
    densTileC = false(tRw*nrow + (tRw-1)*altMovSep, tCl*ncol + (tCl-1)*altMovSep);
    for rr = 1:tRw
        for cc = 1:tCl
            ro = (rr-1)*(nrow+altMovSep);  co = (cc-1)*(ncol+altMovSep);
            densTileC(ro+(1:nrow), co+(1:ncol)) = densBW;
        end
    end
    get_figure('R3-alt-methodMovies');
    moviefixsc(movC, climA);
    if altMovMark
        hold(gca, 'on');
        contour(gca, double(densTileC), [0.5 0.5], 'LineWidth', 0.75, 'LineColor', [0 1 1]);
    end
    allDrop = cmpDrop;                       % rows follow cmpNM, then sfilt if present
    if ~isempty(sfPanel)
        allDrop(end+1,:) = [mean(dropS(densBW)), mean(dropS(bvBW)), ...
                            mean(dropS(onDendA)), mean(dropS(bgA))];
    end
    title(sprintf(['f=%d  (7f-ALT) corrected movie by method: [%s ; %s]   cyan = densBW\n' ...
        'densBW removed %s %%   DENDRITE %s %%   BACKGROUND %s %%   (shared colour scale)'], ...
        f, strjoin(labC(1:tCl), ' | '), strjoin(labC(tCl+1:end), ' | '), ...
        strjoin(compose('%.1f', 100*allDrop(:,1)'), '/'), ...
        strjoin(compose('%.1f', 100*allDrop(:,3)'), '/'), ...
        strjoin(compose('%.2f', 100*allDrop(:,4)'), '/')));
    fprintf('[R3-alt] R3-alt-methodMovies: [%s ; %s], %d frames\n', ...
        strjoin(labC(1:tCl), ' | '), strjoin(labC(tCl+1:end), ' | '), nQA);
    clear cmpPanel panC uncM movC rowSepC colSepC densTileC sfPanel
end

% ---- the same comparison as a still: variance before and after --------------------
% The movie shows single frames; this shows the whole tracked window at once, which is
% where a small persistent bias is visible and a frame-by-frame look is not.
varUnc = varInA;                          % already computed above; variance ignores the sign
varCor = max(varUnc - varFitA, 0);        % var(raw) - var(fitted) IS the corrected variance
get_figure('R3-alt-correctedVariance');
tlC = tiledlayout(2, 3, 'TileSpacing','compact','Padding','compact');
climV = [0 prctile(varUnc(:), 99)];
axc1 = nexttile; imagesc(axc1, varUnc, climV); axis(axc1,'image','off');
colormap(axc1, turbo); colorbar(axc1); hold(axc1,'on');
contour(axc1, densBW, [0.5 0.5], 'w', 'LineWidth', 0.8);
title(axc1, 'var(uncorrected)');
axc2 = nexttile; imagesc(axc2, varCor, climV); axis(axc2,'image','off');
colormap(axc2, turbo); colorbar(axc2); hold(axc2,'on');
contour(axc2, densBW, [0.5 0.5], 'w', 'LineWidth', 0.8);
title(axc2, 'var(corrected)   -- same scale');
axc3 = nexttile; imagesc(axc3, 100*dropA, [0 100*prctile(dropA(:), 99)]);
axis(axc3,'image','off'); colormap(axc3, turbo); colorbar(axc3); hold(axc3,'on');
contour(axc3, densBW, [0.5 0.5], 'w', 'LineWidth', 0.8);
contour(axc3, onDendA, [0.5 0.5], 'g', 'LineWidth', 0.8);
title(axc3, sprintf('%% removed  (white = densBW, green = dendrite)   band = global, %.0f-%.0f%%', ...
    100*prctile(dropA(bgA), 5), 100*prctile(dropA(bgA), 95)));

% ---- the same map with the global band divided out --------------------------------
% The raw map above is hard to read because a smooth global component -- scatter and
% out-of-focus light from the sample, strongest where the tissue is thickest -- sits
% under everything and forms a horizontal band. It is only a few percentage points wide,
% but the colour scale is set by the 99th percentile of the whole map, so those few
% points occupy a large part of the ramp and read as structure.
%
% Dividing each row by its OWN background mean removes that baseline and leaves
% selectivity: 1 means "this pixel lost exactly what an empty pixel in the same row
% lost", and only values above 1 are the design finding something. The band flattens,
% and what survives is the vessels -- or, if nothing survives, there was never anything
% but the band.
%
% Rows with too little background to average are filled from their neighbours rather
% than left as NaN, so a row crossing the middle of the cell does not punch a hole.
bgRow = arrayfun(@(r) mean(dropA(r, bgA(r,:)), 'omitnan'), 1:nrow)';
bgCnt = sum(bgA, 2);
bad   = ~isfinite(bgRow) | bgCnt < 10;
if any(~bad)
    bgRow(bad) = interp1(find(~bad), bgRow(~bad), find(bad), 'linear', 'extrap');
    bgRow = movmean(bgRow, 5);
    dropRel = dropA ./ max(bgRow, eps);
    axc5 = nexttile;
    imagesc(axc5, dropRel, [0 prctile(dropRel(:), 99)]);
    axis(axc5,'image','off'); colormap(axc5, turbo); colorbar(axc5); hold(axc5,'on');
    contour(axc5, densBW, [0.5 0.5], 'w', 'LineWidth', 0.8);
    title(axc5, sprintf('removed / row background   (1 = no selectivity, densBW %.2fx)', selD));
else
    dropRel = dropA;
    nexttile; axis off
    title('too little background to normalise by row');
end

% ---- the row profile: is the removal a BAND, or is it on the vessels? --------------
% A broad horizontal band in the map above is ambiguous by eye, because the tissue is a
% horizontal band too -- removal that merely follows the cell looks the same as removal
% that follows the vessels. This collapses each image row to its mean and puts the three
% populations side by side against the chance floor, which separates the two readings:
%   removal tracks the VESSELS  -> the densBW curve stands well above the other two, and
%                                  background sits flat on the chance line;
%   removal tracks the TISSUE   -> dendrite rises with densBW across the whole cell;
%   removal is GLOBAL           -> background rises too, and the band is an artefact
%                                  (residual motion, illumination, an unfiltered row),
%                                  not vascular signal at all.
axc4 = nexttile; hold(axc4, 'on');
rowMean = @(M, msk) arrayfun(@(r) mean(M(r, msk(r,:)), 'omitnan'), 1:nrow);
pDens = rowMean(100*dropA, densBW);
pDend = rowMean(100*dropA, onDendA);
pBg   = rowMean(100*dropA, bgA);
plot(axc4, pDens, 1:nrow, 'r-', 'LineWidth', 1.4, 'DisplayName','density mask');
plot(axc4, pDend, 1:nrow, 'g-', 'LineWidth', 1.2, 'DisplayName','dendrite off-vessel');
plot(axc4, pBg,   1:nrow, 'k-', 'LineWidth', 1.0, 'DisplayName','background');
% HandleVisibility off, or the reference lines join the legend as "data1"/"data2"
xline(axc4, 100*chanceA, 'b--', 'chance', 'LineWidth', 1.2, 'HandleVisibility', 'off');
xline(axc4, 100*bgMean, 'm:', 'background', 'LineWidth', 1.2, 'HandleVisibility', 'off');
set(axc4, 'YDir', 'reverse'); ylim(axc4, [1 nrow]); box(axc4, 'off');
xlabel(axc4, '% removed'); ylabel(axc4, 'image row'); legend(axc4, 'Location','southeast');
% The question the panel answers is which curve separates from which. Black and green
% together, red apart = vessel-selective removal on a global baseline. Green above black
% = the design is on the cell for its own sake.
title(axc4, sprintf('per-row mean  (dendrite %.2fx background)', selN));

title(tlC, sprintf(['f=%d  (7f-ALT/%s) what the correction removed: %.2f%% density mask / ' ...
    '%.2f%% dendrite / %.2f%% background, chance %.2f%%  <- background at chance is the ' ...
    'test that the band is real'], f, altDecomp, 100*mean(dropA(densBW)), ...
    100*mean(dropA(onDendA)), 100*mean(dropA(bgA)), 100*chanceA));
% MvA is held all the way to here because the 'local' fits read it; it is the largest
% thing this cell allocates, so it goes as soon as the last fit is done.
clear aI aF aC aS fitA rowSepA colSepA IallA movA MvA Sseg SgS

% ---- hand the design to (7g) / (7h) under the names they expect -------------------
% (7g) does NOT reuse flowReg: it re-projects the SPATIAL modes onto the full record, so
% what it actually needs is segU / segPx / useSeg / segStat (columns 4 and 6) and the
% per-segment bookkeeping. Those are what is published here.
if altExport
    useSeg   = useSegA;
    segPx    = segPxA;
    segStat  = aStat;                        % [nTrk speed spacing fTrk fCnt fPeak]
    segBnd   = aBnd;
    segU     = aModeU;
    segUsc   = aModeSc;
    segNPC   = cellfun(@(t) size(t,1), aModeTr);
    segMode  = aMode;                        % 1 = modes found, 2 = mean fallback
    segPCidx = cellfun(@(t) (1:size(t,1))', aModeTr, 'uni', 0);
    flowReg    = altReg;     flowRegNm = altRegNm;   regSeg  = altSegOf;
    regScore   = scoreA;     regF0     = f0A;        covReg7 = covA;
    rVolt      = rVoltA;     rVoltAlt  = nan(size(rVoltA));
    wVolt      = double(~voltKillA(keepA));  % 1 = survived the veto; no ramp in this cell
    regScoreOsc = scoreA;    regPk     = max(abs(covA),[],1)';
    oscCos = altReg(1:2:end, :);   oscSin = altReg(2:2:end, :);
    oscCPCA = ismember(altDecomp, {'cpca','local'});  % (7g) infers its own path from this
    if strcmp(altDecomp, 'local')
        % (7g) fits its rebuilt design to the WHOLE frame. It has no notion of confining a
        % segment's rows to that segment, so carrying 'local' through it silently turns it
        % back into 'cpca'. Said here rather than discovered from the numbers later.
        fprintf(['         WARNING (7g) fits globally and has no local mode, so running it ' ...
                 'on this\n         export reverts the confinement -- the (7f-ALT) figures ' ...
                 'are where ''local'' is\n         actually evaluated.\n']);
    end
    flowRegSrc = sprintf('7f-ALT/%s', altDecomp);
    fprintf(['\n[R3-alt] EXPORTED under the (7f) names: useSeg / segPx / segU / segStat / ' ...
             'flowReg.\n         (7g) and (7h) will now run on this design. ' ...
             'flowRegSrc = ''%s''.\n'], flowRegSrc);
    fprintf(['         NOTE (7g) defaults to g7RefitPCA = true, which RE-DERIVES the ' ...
             'modes with\n         complex PCA on the full record and so discards the ' ...
             '%s decomposition.\n         Set g7RefitPCA = false there to carry this ' ...
             'cell''s modes through.\n'], altDecomp);
else
    flowRegSrc = '7f';
    fprintf('\n[R3-alt] altExport = false: (7f)''s variables are untouched; alt* only.\n');
end

%% (7g) The flow regressor over the WHOLE segment, and the corrected movie
% (7f) built the oscillator on trkFrames only -- the first min(2000, nTseg) frames, i.e.
% the window the particle tracking ran on. This cell extends it to every frame of
% segFrames and shows what it does to the movie.
%
% DO THE VESSEL FOOTPRINTS HAVE TO BE RE-COMPUTED?  No -- and they should NOT be.
% Two different objects came out of (7d)-(7f) and they behave differently in time:
%
%   SPATIAL (reusable):  densBW -> segPx (which pixels belong to which vessel segment)
%                        and segU{k} (the spatial eigenvectors inside each segment).
%     These describe WHERE the vessel is. The movie is motion-corrected and the mean
%     image is subtracted, so the vessel does not move between frame 2000 and frame
%     10000; the geometry that 2000 frames measured is the geometry for all of them.
%     Re-running the detection on the rest of the movie would only re-estimate the same
%     mask from a different noise realisation. It is reused verbatim below.
%
%   TEMPORAL (must be recomputed):  oscCos / oscSin.
%     These are A(t)cos(phi(t)) -- defined only on the frames they were built from. They
%     are rebuilt here by projecting the SAME spatial modes onto S over all nTseg frames
%     and redoing the bandpass + Hilbert on the full-length trace. No new tracking is
%     needed: the tracking only ever supplied (i) the segments and (ii) the search band
%     via v/d, and both are already in hand (segPx, segStat(:,4)).
%
%   The one thing that CAN go stale is the passband. fPeak was the PSD peak measured on
%   2000 frames; if the flow speeds up or slows down over the 10 s, the peak moves. That
%   is measured explicitly below (per-block PSD peak and per-block R^2). If the peak
%   drifts by less than oscBandFrac the single generous band still covers it -- if it
%   drifts more, the warning fires and the segments should be re-tracked in blocks.
%
% ------------------------------- parameters -------------------------------
g7BandMode  = 'reuse';   % 'refit' = re-find the PSD peak on the FULL-length trace
                         % 'reuse' = keep fPeak from (7f) (segStat(:,6))
g7Fmin      = [];        % Hz, frequency floor for the full-record refit. [] = inherit
                         % oscFmin from (7f). Same three jobs as there: it bounds the
                         % search, it bounds the passband AFTER the +/-oscBandFrac
                         % widening, and it bounds the fallback. Refitting on 5x more
                         % data can move the peak, so the floor has to be re-applied here
                         % -- a segment that was legal on the tracked window can pin on
                         % the full record, and that is a finding, not a nuisance.
g7RefitPCA  = true;      % false = reuse the (7f) spatial eigenvectors segU{k} AND the set
                         %         of PCs it selected; only their TIME COURSES are redone
                         % true  = re-run the per-segment PCA on the full record and
                         %         re-select the oscillatory PCs from scratch. Reuses only
                         %         segPx (the density-mask segments) and the v/d search
                         %         band. Forces g7KeepMode = 'rescore', because the PC
                         %         indices and counts -- and so the regressor names --
                         %         change and can no longer be matched to (7f).
g7KeepMode  = 'rescore'; % 'match7f' = use exactly the regressor set (7f) settled on,
                         %             matched by name and in the same order
                         % 'rescore' = re-score every row on the full record and re-apply
                         %             oscScoreMin (more data, so a better-founded cut).
                         %             This is a FILTER ON ROWS THAT ALREADY EXIST -- it
                         %             cannot add a PC (7f) rejected and it cannot change
                         %             an eigenvector. Re-deriving the footprint is
                         %             g7RefitPCA, a separate switch.
% -- kymograph of the dendritic ROI traces (Result.ftprnt) ------------------------
% -- how a dendritic ROI trace is extracted (used by BOTH 7g kymographs) ----------
g7FpMode    = 'mean';    % 'mean'     plain average over every pixel of the ROI footprint
                         %            (mask/nPix). All pixels count equally.
                         % 'weighted' unit-sum footprint weights, i.e. a weighted average.
                         %            Result.ftprnt spans 0..8.3e6, so this is dominated
                         %            by the few brightest centre pixels.
g7FpThr     = 0;         % a pixel joins the mask if ftprnt > g7FpThr * max(ftprnt) for
                         % that ROI. 0 = every non-zero pixel (~685 px per ROI here).
g7Kymo      = true;      % draw [mov_res | fitted vascular | corrected] as kymographs
g7KymoBands = {'full','fast'};   % rows of that figure. 'full' = dF, 'fast' = dF minus
                         % its 20-ms average (where the blood is), 'sub' = the 20-ms
                         % average (the trace the paper reports). Any subset, any order.
g7KymoFrames = [];       % [] = the whole segment, else a frame range
g7KymoBase  = 1001;      % ms, running-median baseline for the displayed dF
g7BlockS    = 1;         % s, block length for the drift / stability diagnostic
g7FitChunk  = 500;       % frames per pass when accumulating the normal equations
g7Browse    = true;      % open the interactive regressor browser on the full record
g7MovFrm    = 2000;      % frames in the 1x3 movie. Set to Inf for the whole segment --
                         % cost is 3*nrow*ncol*nFrames*4 bytes (printed before it builds)
g7MovBin    = 1;         % temporal binning for that movie (>1 lets the whole record fit,
                         % at the price of smearing the oscillation you are looking at)
g7MovSep    = 3;         % px, divider width
% --------------------------------------------------------------------------

assert(exist('segPx','var') && exist('useSeg','var') && exist('segU','var') && ...
       exist('flowRegNm','var') && exist('segStat','var'), ...
       'Run (7f) first (it builds segPx/useSeg/segU/segStat/flowRegNm).');
assert(size(S,3) == nTseg, 'S has %d frames but nTseg = %d; re-run (7c).', size(S,3), nTseg);

zrow  = @(V) (V - mean(V,2)) ./ max(std(V,0,2), eps);
if ~exist('oscCPCA','var'), oscCPCA = ~isreal(segU{find(~cellfun(@isempty,segU),1)}); end
if oscCPCA, pcaNm = 'complex PCA'; else, pcaNm = 'PCA'; end
nSegU = numel(useSeg);
gAll  = ~badfr(:)';                      % good frames over the WHOLE segment
nG    = nnz(gAll);
fprintf(['\n[R3] (7g) extending the flow regressor from %d tracked frames to all %d ' ...
         'frames (%.1f s, %d good)\n'], numel(trkFrames), nTseg, nTseg/fs, nG);
fprintf('     segments REUSED from (7f): %d, %d px total (densBW is never recomputed here)\n', ...
    nSegU, sum(cellfun(@numel, segPx(useSeg))));

% ---- 1-2. per segment: search band, spatial modes, source trace, phase --------------
% WHAT IS REUSED AND WHAT IS RE-DERIVED -- there are FOUR separate decisions stacked on
% top of each other, and each has its own switch. Going down the list, each one is
% cheaper to redo and less likely to be wrong:
%
%   (a) densBW -> segPx        WHICH PIXELS ARE A VESSEL SEGMENT.
%                              ALWAYS reused. Changing it means re-running the detection
%                              and tracking of (7d)/(7e); geometry does not change with
%                              time in a motion-corrected movie, so there is no reason to.
%   (b) segU{k}, keepPC        THE PCA FOOTPRINT INSIDE THE SEGMENT and which of its modes
%                              count as oscillatory.  <- g7RefitPCA
%                              (7f) estimated these from 2000 frames. They are the only
%                              reused quantity that is a STATISTICAL estimate rather than
%                              a geometric fact, so they are the ones worth redoing on 5x
%                              more data.
%   (c) fPeak, passband        WHERE IN FREQUENCY to look.  <- g7BandMode
%   (d) which design rows survive oscScoreMin.  <- g7KeepMode
%
% NOTE g7KeepMode = 'rescore' ONLY touches (d). It re-scores rows that already exist and
% drops the weak ones; it cannot add a PC that (7f) rejected, and it cannot change a
% single spatial eigenvector. If you want the footprint itself re-derived, that is
% g7RefitPCA -- 'rescore' is not a substitute for it.
%
% reshape(S,...) is a view, so only each segment's own pixels are ever copied.
if isempty(g7Fmin)
    if exist('oscFmin','var'), g7Fmin = oscFmin; else, g7Fmin = 0; end
end
fFloor7 = max(g7Fmin, 1);
if g7RefitPCA && strcmp(g7KeepMode, 'match7f')
    fprintf(['[R3] g7RefitPCA = true re-derives the PCs, so their names and count no ' ...
             'longer match (7f) -- switching g7KeepMode to ''rescore''\n']);
    g7KeepMode = 'rescore';
end
if g7RefitPCA, pcaWhat = 're-derived on the full record'; else, pcaWhat = 'reused from (7f)'; end
fprintf('     PCA footprint %s | band %s | row set %s\n', pcaWhat, g7BandMode, g7KeepMode);

% search band per segment: the (7f) one if it was stored, else rebuilt from
% segStat(:,4) = fTrk = v/d with the identical clamp logic. Either way no re-tracking.
haveSrch = exist('segSrchB','var') && isequal(size(segSrchB), [nSegU 2]);
fLim     = [fFloor7, 0.45*fs];
segSrch  = nan(nSegU, 2);
for k = 1:nSegU
    if haveSrch
        srchBand = segSrchB(k,:);
    else
        fTrk = segStat(k,4);
        if ~isfinite(fTrk) || fTrk < fLim(1) || fTrk > fLim(2)
            srchBand = fLim;
        else
            srchBand = [fTrk/oscSearch(2), fTrk*oscSearch(2)];
        end
    end
    srchBand = [max(srchBand(1), fFloor7), min(srchBand(2), fs/2 - 1)];
    if srchBand(2) <= srchBand(1)*1.05, srchBand = [fFloor7, fs/2 - 1]; end
    segSrch(k,:) = srchBand;
end

Sv7       = reshape(S, [], nTseg);
% mov_nv, un-bandpassed, for the browser's source panel ONLY -- segYfull (from S) still
% drives the PSD, the bandpass and the Hilbert, because picking the vessel's spatial
% scale is what S is for. See the same split in (7f).
haveMnv7 = exist('mov_nv','var') && size(mov_nv,3) == nTseg;
if haveMnv7, Mnv7 = reshape(mov_nv, [], nTseg); else, Mnv7 = []; end
segYnv    = zeros(nSegU, nTseg);          % display source, from mov_nv
segYfull  = zeros(nSegU, nTseg);          % derivation source, from S
segPCfull = cell(1, nSegU);
segU7     = cell(1, nSegU);   segPCidx7 = cell(1, nSegU);   segUsc7 = cell(1, nSegU);
segNPC7   = zeros(1, nSegU);  segMode7  = zeros(1, nSegU);  segAlign = nan(1, nSegU);
oscCosF   = zeros(nSegU, nTseg);  oscSinF = oscCosF;
segBandF  = nan(nSegU, 2);  segFPeakF = nan(1, nSegU);  segPinF = false(1, nSegU);
srcNm7    = cell(1, nSegU);
if g7RefitPCA
    fprintf('[R3] re-deriving the PCA footprint of each segment on all %d good frames:\n', nG);
end
for k = 1:nSegU
    sgi   = useSeg(k);
    pxIdx = segPx{sgi};
    Yseg  = -double(Sv7(pxIdx, :));          % same sign convention as (7f)
    Yseg  = Yseg - mean(Yseg, 2);            % mean over the FULL record now

    if ~g7RefitPCA
        % ---- reuse (7f)'s modes: same eigenvectors, new time courses ----------------
        segU7{k} = segU{k};  segPCidx7{k} = segPCidx{k};  segUsc7{k} = segUsc{k};
        segNPC7(k) = segNPC(k);  segMode7(k) = segMode(k);
        if segMode(k) == 1 && ~isempty(segU{k})
            % Complex modes must be projected onto the ANALYTIC signal, not onto Yseg:
            % they were derived from hilbert(bandpass(Yseg)), so projecting the raw real
            % traces would mix the two quadratures and lose the phase the mode encodes.
            if oscCPCA
                Yb7 = get_bandpass(Yseg, fs, segSrch(k,:));  Yb7(~isfinite(Yb7)) = 0;
                Zseg = hilbert(Yb7.').';  Zseg = Zseg - mean(Zseg, 2);
                segPCfull{k}  = segU{k}' * Zseg;
                segYfull(k,:) = real(segPCfull{k}(1,:));
            else
                segPCfull{k}  = segU{k}' * Yseg;
                segYfull(k,:) = segPCfull{k}(1,:);
            end
        else
            segYfull(k,:) = mean(Yseg, 1);
        end
    else
        % ---- re-derive the footprint: the (7f) PCA block, on the full record --------
        % Same rules (segMinPxPCA, segBinPix, oscScoreMin, oscMaxPC, score x variance
        % ranking) applied to 5x more data, so the modes are better determined and a PC
        % that was ambiguous over 2 s can now be resolved either way.
        if numel(pxIdx) >= segMinPxPCA
            nPCscan = max(2, ceil(numel(pxIdx) / segBinPix));
            nPCscan = min(nPCscan, min(size(Yseg,1), nG) - 1);
            if oscCPCA
                [Us, pcTr, varPC] = get_cpca(Yseg, fs, segSrch(k,:), nPCscan, gAll);
                scPC = get_oscScore(real(pcTr), fs, segSrch(k,:));
            else
                [Us, Ss, ~] = svds(Yseg(:, gAll), nPCscan);
                pcTr  = Us' * Yseg;
                varPC = diag(Ss).^2;  varPC = varPC / max(sum(varPC), eps);
                scPC  = get_oscScore(pcTr, fs, segSrch(k,:));
            end
            segNPC7(k) = numel(scPC);
            keepPC = find(scPC >= oscScoreMin);
            [~, ordPC] = sort(scPC(keepPC) .* varPC(keepPC), 'descend');
            keepPC = keepPC(ordPC);
            if numel(keepPC) > oscMaxPC, keepPC = keepPC(1:oscMaxPC); end
            if isempty(keepPC)
                warning(['Segment %d: no PC reached oscScoreMin = %.2f on the full ' ...
                    'record (best %.2f); falling back to the segment mean.'], ...
                    sgi, oscScoreMin, max(scPC));
                segYfull(k,:) = mean(Yseg, 1);  segMode7(k) = 2;
            else
                segU7{k} = Us(:, keepPC);  segPCidx7{k} = keepPC;  segUsc7{k} = scPC(keepPC);
                segPCfull{k}  = pcTr(keepPC, :);
                segYfull(k,:) = real(pcTr(keepPC(1), :));
                segMode7(k)   = 1;
                % how much did re-deriving actually move the footprint? |u_new . u_old|
                % is 1 if the 2000-frame estimate was already right, and well below 1 if
                % the leading mode genuinely changed.
                if ~isempty(segU{k})
                    segAlign(k) = abs(segU7{k}(:,1)' * segU{k}(:,1));   % ' = conj'
                end
            end
            fprintf(['   seg %3d: %4d px  PCA over %d modes -> %d of %d kept ' ...
                     '(was %d of %d in (7f)) | leading-mode overlap %.3f\n'], ...
                sgi, numel(pxIdx), nPCscan, numel(keepPC), numel(scPC), ...
                numel(segPCidx{k}), segNPC(k), segAlign(k));
        else
            segYfull(k,:) = mean(Yseg, 1);  segMode7(k) = 3;
            fprintf('   seg %3d: %4d px  own trace (< segMinPxPCA = %d)\n', ...
                sgi, numel(pxIdx), segMinPxPCA);
        end
    end

    % ---- the browser's source trace: same spatial weights, applied to mov_nv -----
    if haveMnv7
        Ynv7 = -double(Mnv7(pxIdx, :));
        Ynv7 = Ynv7 - mean(Ynv7, 2);
        if segMode7(k) == 1 && ~isempty(segU7{k})
            segYnv(k,:) = real(segU7{k}(:,1)' * Ynv7);
        else
            segYnv(k,:) = mean(Ynv7, 1);
        end
    else
        segYnv(k,:) = segYfull(k,:);      % mov_nv gone: fall back to S
    end

    % ---- band + Hilbert on the full-length trace --------------------------------
    y = segYfull(k,:);
    if strcmp(g7BandMode, 'reuse')
        fPeak = max(segStat(k,6), fFloor7);
    else
        yg = y;  yg(~gAll) = NaN;
        [Py, Fy] = nanPSD(yg, fs, min(1024, 2^floor(log2(nTseg/2))));
        Py = Py(:);  Fy = Fy(:);
        inS = Fy >= segSrch(k,1) & Fy <= segSrch(k,2);
        if ~any(inS & isfinite(Py)), inS = Fy > fFloor7 & Fy < fs/2-1; end
        Pm = Py;  Pm(~inS | ~isfinite(Pm)) = -inf;
        [~, ipk] = max(Pm);  fPeak = Fy(ipk);
    end
    segFPeakF(k) = fPeak;
    segPinF(k)   = g7Fmin > 0 && fPeak <= fFloor7 * 1.05;

    bnd = fPeak * [1-oscBandFrac, 1+oscBandFrac];
    bnd = [max(bnd(1), fFloor7), min(bnd(2), fs/2-1)];   % floor AFTER the widening
    if ~all(isfinite(bnd)) || bnd(2) <= bnd(1)*1.05, bnd = segSrch(k,:); end
    segBandF(k,:) = bnd;

    yb = get_bandpass(y, fs, bnd);  yb(~isfinite(yb)) = 0;
    z  = hilbert(yb);
    oscCosF(k,:) = abs(z).*cos(angle(z));
    oscSinF(k,:) = abs(z).*sin(angle(z));

    % name the source, exactly as (7f) does, but from the modes actually in use here
    if segMode7(k) == 1
        if oscCPCA, pcTag7 = 'cpcaPC'; else, pcTag7 = 'svdPC'; end
        srcNm7{k} = sprintf('seg%d %s%d (%d kept of %d)', sgi, pcTag7, ...
            segPCidx7{k}(1), numel(segPCidx7{k}), segNPC7(k));
    elseif segMode7(k) == 2
        srcNm7{k} = sprintf('seg%d mean (%s found none of %d)', sgi, pcaNm, segNPC7(k));
    else
        srcNm7{k} = sprintf('seg%d own trace (small)', sgi);
    end
end
clear Sv7 Mnv7 Ynv7 Yseg pcTr Us Ss Yb7 Zseg

if g7RefitPCA
    okA = isfinite(segAlign);
    if ~any(okA)
        fprintf(['[R3] footprint re-derivation: no segment has a leading PC in both ' ...
                 '(7f) and (7g), so there is nothing to compare\n']);
    else
    fprintf(['[R3] footprint re-derivation: leading-mode overlap with (7f) median %.3f ' ...
             '(min %.3f over %d segments)\n'], median(segAlign(okA)), min(segAlign(okA)), nnz(okA));
    if any(segAlign(okA) < 0.8)
        fprintf(['     %d segment(s) below 0.8 -- their 2000-frame footprint was NOT ' ...
                 'the full-record one, so refitting mattered: %s\n'], ...
            nnz(segAlign(okA) < 0.8), mat2str(useSeg(okA & segAlign < 0.8)));
    end
    end
end
if g7Fmin > 0
    fprintf('     frequency floor %g Hz; passbands: %s\n', g7Fmin, ...
        strjoin(arrayfun(@(k) sprintf('s%d %.0f-%.0f', useSeg(k), segBandF(k,1), ...
        segBandF(k,2)), 1:nSegU, 'uni', 0), ', '));
end
if any(segPinF)
    fprintf(['     %d segment(s) pin at the %g Hz floor on the full record: %s ' ...
             '(their peak sits below the floor over 10 s even if it did not over 2 s)\n'], ...
        nnz(segPinF), fFloor7, mat2str(useSeg(segPinF)));
end

% ---- 3. does the 2000-frame band still hold at frame 10000? -------------------------
% The direct test of "can I reuse this": per block, (a) where the spectral peak actually
% sits and (b) how much of the segment signal the FULL-RECORD oscillator explains there.
% If the footprint and band had gone stale, R^2 would fall off in the later blocks.
blkLen = max(round(g7BlockS*fs), 64);
nBlk   = max(1, floor(nTseg/blkLen));
fPkBlk = nan(nSegU, nBlk);  r2Blk = nan(nSegU, nBlk);
for k = 1:nSegU
    for b = 1:nBlk
        ib = (b-1)*blkLen + (1:blkLen);
        yb0 = segYfull(k, ib);  yb0 = yb0 - mean(yb0);
        yb0(~gAll(ib)) = 0;                          % gate, do not shift the spectrum
        w   = 0.5*(1 - cos(2*pi*(0:blkLen-1)/(blkLen-1)));
        Pw  = abs(fft(yb0 .* w)).^2;
        fv  = (0:blkLen-1) * fs / blkLen;
        inB = fv >= segSrch(k,1) & fv <= segSrch(k,2) & fv < fs/2;
        if any(inB)
            Pw(~inB) = -inf;  [~, ip] = max(Pw);  fPkBlk(k,b) = fv(ip);
        end
        gb = gAll(ib);
        if nnz(gb) > 20
            yv = segYfull(k, ib(gb))';  yv = yv - mean(yv);
            X1 = [oscCosF(k, ib(gb)); oscSinF(k, ib(gb)); ones(1, nnz(gb))]';
            r2Blk(k,b) = 1 - var(yv - X1*(X1\yv)) / max(var(yv), eps);
        end
    end
end
drift = (max(fPkBlk,[],2,'omitnan') - min(fPkBlk,[],2,'omitnan')) ./ segFPeakF(:);

fprintf('\n[R3] is the (7f) band still valid over the whole record?  (%d blocks of %.2f s)\n', ...
    nBlk, blkLen/fs);
fprintf('     %-5s %-10s %-10s %-12s %-12s %s\n', 'seg', 'fPeak 7f', 'fPeak full', ...
    'peak drift', 'R^2 trk win', 'R^2 later blocks');
for k = 1:nSegU
    nTrkBlk = max(1, floor(numel(trkFrames)/blkLen));
    r2Early = mean(r2Blk(k, 1:min(nTrkBlk,nBlk)), 'omitnan');
    r2Late  = mean(r2Blk(k, min(nTrkBlk,nBlk)+1:end), 'omitnan');
    fl = '';
    if drift(k) > oscBandFrac, fl = '   <- DRIFTS OUT OF BAND'; end
    fprintf('     %-5d %-10.1f %-10.1f %-12.2f %-12.3f %.3f%s\n', useSeg(k), ...
        segStat(k,6), segFPeakF(k), drift(k), r2Early, r2Late, fl);
end
if any(drift > oscBandFrac)
    warning(['%d segment(s) move their spectral peak by more than oscBandFrac = %.2f ' ...
        'across the record: one band cannot cover the whole movie there. Either raise ' ...
        'oscBandFrac, or re-run (7d)-(7f) on blocks and concatenate the regressors.'], ...
        nnz(drift > oscBandFrac), oscBandFrac);
end

% ---- 4. rebuild the design in the (7f) canonical order --------------------------------
% srcNm7 / segPCidx7 / segNPC7 are the modes ACTUALLY in use: they are (7f)'s when
% g7RefitPCA is false, and the newly derived ones when it is true.
flowRegF   = [zrow(oscCosF); zrow(oscSinF)];
flowRegNmF = [cellfun(@(v) [v ' | cos'], srcNm7, 'uni',0), ...
              cellfun(@(v) [v ' | sin'], srcNm7, 'uni',0)];
regSegF    = [1:nSegU, 1:nSegU];
for k = 1:nSegU
    if isempty(segPCfull{k}) || size(segPCfull{k},1) < 2, continue; end
    extra = segPCfull{k}(2:end, :);
    qIdx  = 2:size(segPCfull{k},1);
    if oscCPCA                                   % complex mode -> its own quadrature pair
        flowRegF   = [flowRegF; zrow(real(extra)); zrow(imag(extra))];          %#ok<AGROW>
        flowRegNmF = [flowRegNmF, ...
            arrayfun(@(q) sprintf('seg%d cpcaPC%d (%d kept of %d) | extra re', ...
                useSeg(k), segPCidx7{k}(q), numel(segPCidx7{k}), segNPC7(k)), qIdx, 'uni',0), ...
            arrayfun(@(q) sprintf('seg%d cpcaPC%d (%d kept of %d) | extra im', ...
                useSeg(k), segPCidx7{k}(q), numel(segPCidx7{k}), segNPC7(k)), qIdx, 'uni',0)]; %#ok<AGROW>
        regSegF    = [regSegF, repmat(k, 1, 2*size(extra,1))];                  %#ok<AGROW>
    else
        flowRegF   = [flowRegF; zrow(extra)];                                   %#ok<AGROW>
        flowRegNmF = [flowRegNmF, arrayfun(@(q) ...
            sprintf('seg%d svdPC%d (%d kept of %d) | extra PC', useSeg(k), ...
            segPCidx7{k}(q), numel(segPCidx7{k}), segNPC7(k)), qIdx, 'uni',0)];  %#ok<AGROW>
        regSegF    = [regSegF, repmat(k, 1, size(extra,1))];                    %#ok<AGROW>
    end
end
dC = oscCosF;  dS = oscSinF;
for dv = 1:flowUseDeriv
    dC = gradient(dC')';  dS = gradient(dS')';
    flowRegF   = [flowRegF; zrow(dC); zrow(dS)];                               %#ok<AGROW>
    flowRegNmF = [flowRegNmF, ...
        cellfun(@(v) sprintf('%s | d%d cos', v, dv), srcNm7, 'uni',0), ...
        cellfun(@(v) sprintf('%s | d%d sin', v, dv), srcNm7, 'uni',0)];        %#ok<AGROW>
    regSegF    = [regSegF, 1:nSegU, 1:nSegU];                                  %#ok<AGROW>
end
clear dC dS

[regScoreF, regF0F] = get_oscScore(flowRegF, fs, [1 fs/2-1]);

% ---- the voltage-map criterion, same as (7f) ---------------------------------------
% The keep decision happens HERE, before the fit that produces covRegF, so the spatial
% maps have to be measured first. They are estimated on a SUBSAMPLE of the good frames:
% a covariance map is a spatial pattern and a few thousand frames pin it down, so there
% is no reason to pay for a second full pass over an 8 GB movie.
% If (7g) is run without (7f) in the workspace, rebuild the references it needs.
if ~exist('oscVoltPen','var') || ~exist('oscVoltRefUse','var')
    oscVoltPen = true;  oscVoltRefUse = 'ftprnt';
    oscVoltLo  = 0.3;   oscVoltHi  = 0.6;  oscVoltMin = 0.02;
end
if ~exist('voltMask','var'), voltMask = (neuronMap > 0) | bvBW | densBW; end
if ~exist('voltFp','var'),   voltFp   = neuronMap ./ max(neuronMap(:), [], 'omitnan'); end
if ~exist('voltSTA','var'),  voltSTA  = []; end

wVoltF = ones(size(regScoreF));  rVoltF = zeros(size(regScoreF));
if true                                   % always measured; oscVoltPen only gates the use
    subF = find(gAll);  subF = subF(1:g7VoltSub:end);
    Rz2  = flowRegF(:, subF) - mean(flowRegF(:, subF), 2);
    Mc2  = zeros(nrow*ncol, size(flowRegF,1));
    for c0 = 1:g7FitChunk:numel(subF)
        jj = c0:min(c0+g7FitChunk-1, numel(subF));
        Mc2 = Mc2 + tovec(mov_res(:,:,subF(jj))) * Rz2(:,jj)';
    end
    covCand = Mc2 / max(numel(subF)-1, 1);
    if strcmp(oscVoltRefUse,'sta') && ~isempty(voltSTA), vRefF = voltSTA; else, vRefF = voltFp; end
    [wVoltF, rVoltF] = get_voltPenalty(covCand, vRefF, voltMask, ...
        oscVoltLo, oscVoltHi, oscVoltMin);
    fprintf(['[R3] voltage-map criterion on %d of %d frames (every %dth): |r| median ' ...
             '%.2f, max %.2f\n'], numel(subF), nG, g7VoltSub, median(abs(rVoltF)), ...
             max(abs(rVoltF)));
    if oscVoltPen
        nKill = nnz(regScoreF >= oscScoreMin & regScoreF .* wVoltF < oscScoreMin);
        regScoreF = regScoreF .* wVoltF;
        if nKill > 0
            fprintf('     %d row(s) dropped below oscScoreMin by the penalty alone:\n', nKill);
            for rr = find(wVoltF < 0.5)'
                fprintf('       %-46s r = %+.2f, x%.2f\n', flowRegNmF{rr}, rVoltF(rr), wVoltF(rr));
            end
        end
    end
    clear Mc2 Rz2 covCand
end

switch g7KeepMode
    case 'match7f'
        [tf, loc] = ismember(flowRegNm, flowRegNmF);
        assert(all(tf), ['%d regressor name(s) from (7f) are not in the rebuilt design; ' ...
            '(7f) and (7g) are out of sync -- re-run (7f).'], nnz(~tf));
        selF = loc(:);
    case 'rescore'
        keepF = regScoreF >= oscScoreMin;
        for k = 1:nSegU, keepF([k, k+nSegU]) = any(keepF([k, k+nSegU])); end
        if regDropBelowScore
            fprintf('[R3] re-scored on %d frames: dropping %d of %d rows below %.2f\n', ...
                nTseg, nnz(~keepF), numel(keepF), oscScoreMin);
        else
            keepF(:) = true;
        end
        selF = find(keepF);
        [~, oF] = sortrows([regSegF(selF)', (1:numel(selF))']);
        selF = selF(oF);
    otherwise
        error('g7KeepMode must be ''match7f'' or ''rescore''.');
end
flowRegF = flowRegF(selF,:);   flowRegNmF = flowRegNmF(selF);
regSegF  = regSegF(selF);      regScoreF  = regScoreF(selF);
regF0F   = regF0F(selF);       wVoltF     = wVoltF(selF);
rVoltF   = rVoltF(selF);

% A segment that pins at the floor over the full record has no measured peak above
% g7Fmin, so its band was imposed rather than found -- the same reason (7f) drops pinned
% segments. Drop its rows here too, whatever keep mode selected them.
if exist('oscDropPinned','var') && oscDropPinned && any(segPinF)
    kpR = ~segPinF(regSegF);
    fprintf('[R3] dropping %d rows from %d segment(s) pinned at the %g Hz floor: %s\n', ...
        nnz(~kpR), nnz(segPinF), fFloor7, mat2str(useSeg(segPinF)));
    flowRegF = flowRegF(kpR,:);  flowRegNmF = flowRegNmF(kpR);
    regSegF  = regSegF(kpR);     regScoreF  = regScoreF(kpR);
    regF0F   = regF0F(kpR);      wVoltF     = wVoltF(kpR);
    rVoltF   = rVoltF(kpR);
    assert(~isempty(flowRegF), ['Every segment pins at g7Fmin = %g Hz over the full ' ...
        'record: nothing above the floor is left to regress out.'], g7Fmin);
end
nRegF    = size(flowRegF,1);
fprintf('[R3] full-record design: %d x %d (%s), segments %s\n', nRegF, nTseg, ...
    g7KeepMode, mat2str(useSeg(regSegF)));

% ---- 5. fit on every good frame, without copying the movie ---------------------------
% Same normal equations SeeResiduals solves -- C = (Y*I')/(I*I') with I = [1; centred
% regressors] -- accumulated in time chunks so no second copy of mov_res is made. The
% per-pixel mean and variance of the input are accumulated in the same pass.
Ig    = [ones(1,nG); flowRegF(:,gAll) - mean(flowRegF(:,gAll), 2)];
AtA   = Ig * Ig';
gIdx  = find(gAll);
npx   = nrow*ncol;
MtI   = zeros(npx, nRegF+1);
s1    = zeros(npx,1);  s2 = zeros(npx,1);
for c0 = 1:g7FitChunk:nG
    jj = c0:min(c0+g7FitChunk-1, nG);
    Yc = tovec(mov_res(:,:,gIdx(jj)));           % px x nt
    MtI = MtI + Yc * Ig(:,jj)';
    s1  = s1 + sum(Yc,2);   s2 = s2 + sum(Yc.^2,2);
end
clear Yc
Cfv7 = MtI / AtA;                                 % px x (1+nReg), col 1 = DC
Bf7  = Cfv7(:, 2:end);

varIn7  = reshape((s2 - s1.^2/nG) / max(nG-1,1), nrow, ncol);
Sf7     = cov(flowRegF(:, gAll)');
varFit7 = reshape(sum((Bf7 * Sf7) .* Bf7, 2), nrow, ncol);
dropF7  = varFit7 ./ max(varIn7, eps);
covRegF = MtI(:, 2:end) / max(nG-1, 1);           % Ig rows 2: are centred -> this is cov

fprintf(['\n[R3] over ALL %d good frames the flow regressor removes from mov_res: ' ...
         '%.1f%% overall, %.1f%% in the density mask, %.1f%% on bvMask\n'], ...
    nG, 100*mean(dropF7(:)), 100*mean(dropF7(densBW)), 100*mean(dropF7(bvBW)));
fprintf('     %.1f%% on the DENDRITE away from vessels  <- keep this one small\n', ...
    100*mean(dropF7(neuronMap>0 & ~bvBW)));
if exist('dropF','var')
    fprintf('     (7f), tracked window only:  %.1f%% on bvMask, %.1f%% on the dendrite\n', ...
        100*mean(dropF(bvBW)), 100*mean(dropF(neuronMap>0 & ~bvBW)));
end
if mean(dropF7(neuronMap>0 & ~bvBW)) > 0.5*mean(dropF7(bvBW))
    warning(['Over the full record the regressor removes %.1f%% on the dendrite vs ' ...
        '%.1f%% on bvMask -- probably taking neural signal.'], ...
        100*mean(dropF7(neuronMap>0 & ~bvBW)), 100*mean(dropF7(bvBW)));
end

% ---- 5c. does a TIME-VARYING coefficient fit the vessel better? ----------------------
% The regressor's waveform can be right while its per-pixel coefficient is wrong, because
% one global beta assumes the vessel's contribution to each pixel is constant for the
% whole 10 s. It is not obviously constant: flow speed drifts, the RBC train thins and
% thickens, and the vessel's optical contribution moves with it. So refit beta in
% overlapping windows and see.
%
% THE COMPARISON HAS TO BE OUT-OF-SAMPLE, or it is not a comparison at all. A windowed
% fit has ~nWin times as many free parameters as the global one, so it removes more
% variance in-sample ALWAYS -- including the neural signal it should not touch. Both fits
% therefore train on the same frames and are scored on the same held-out frames: whole
% g7CVChunk-frame blocks are withheld (not single frames, which the neighbouring frames
% would predict anyway), and the winner is whichever removes more variance ON HELD-OUT
% DATA over the vessels while removing no more on the dendrite.
%
% Windows are blended with a linear cross-fade over the overlap, weights normalised to
% sum to 1 per frame, so the fitted movie has no step at a window edge.
%
% ------------------------------- parameters -------------------------------
g7Split     = true;      % run the split-fit comparison
g7SplitWin  = 1000;      % frames per window
g7SplitOvl  = 200;       % frames of overlap between neighbouring windows
g7CVChunk   = 50;        % frames per held-out block
g7CVFrac    = 0.2;       % fraction of frames held out for scoring
g7SplitMinR = 5;         % a window needs this many training frames per regressor, else
                         % it falls back to the global coefficients
g7SplitKymo = true;      % draw R3-flowSplitKymo: raw / fitted / corrected for both fits,
                         % as linked kymographs over Result.ftprnt
g7VoltSub   = 5;         % use every g7VoltSub-th good frame to estimate the regressor
                         % spatial maps for the voltage-map criterion. A covariance MAP is
                         % a spatial pattern -- a few thousand frames pin it down, so this
                         % keeps the criterion from costing a second full pass.
% --------------------------------------------------------------------------
if g7Split
    stepW = max(1, g7SplitWin - g7SplitOvl);
    wS    = 1:stepW:max(1, nTseg - g7SplitWin + 1);
    if wS(end) + g7SplitWin - 1 < nTseg, wS(end+1) = nTseg - g7SplitWin + 1; end
    wS    = unique(max(wS, 1));
    nWin  = numel(wS);
    wE    = min(wS + g7SplitWin - 1, nTseg);

    % cross-fade weights, then normalise so every frame's weights sum to 1
    Wgt = zeros(nWin, nTseg);
    for k = 1:nWin
        L = wE(k) - wS(k) + 1;
        tp = ones(1, L);
        nr = min(g7SplitOvl, floor(L/2));
        if nr > 0
            ramp = (1:nr) / (nr + 1);
            if k > 1,     tp(1:nr)         = ramp;        end   % first window keeps its head
            if k < nWin,  tp(end-nr+1:end) = fliplr(ramp); end   % last keeps its tail
        end
        Wgt(k, wS(k):wE(k)) = tp;
    end
    Wsum = sum(Wgt, 1);
    assert(all(Wsum > 0), 'Some frames are covered by no window; check g7SplitWin/g7SplitOvl.');
    Wgt = Wgt ./ Wsum;

    % held-out blocks, shared by both fits
    chunkId  = ceil((1:nTseg) / g7CVChunk);
    nHold    = max(2, round(1/max(g7CVFrac, eps)));
    testMask = mod(chunkId, nHold) == 0;
    trainF   = gAll & ~testMask;
    testF    = gAll & testMask;
    nTr = nnz(trainF);  nTe = nnz(testF);
    fprintf(['\n[R3] split fit: %d windows of %d frames, %d overlap | train %d frames, ' ...
             'held out %d (%d-frame blocks)\n'], nWin, g7SplitWin, g7SplitOvl, nTr, nTe, g7CVChunk);

    Iall7 = [ones(1,nTseg); flowRegF - mean(flowRegF(:,trainF), 2)];   % centred on TRAIN

    % ---- global fit, trained on trainF only (so the two are comparable) -------------
    % Not the step-5 Cfv7: that one saw every good frame, including the held-out ones,
    % and scoring it here would flatter it.
    gI2 = find(trainF);
    A_g = zeros(nRegF+1);  M_g = zeros(npx, nRegF+1);
    for c0 = 1:g7FitChunk:numel(gI2)
        jj = c0:min(c0+g7FitChunk-1, numel(gI2));  ci = gI2(jj);
        Yc = tovec(mov_res(:,:,ci));
        M_g = M_g + Yc * Iall7(:,ci)';
        A_g = A_g + Iall7(:,ci) * Iall7(:,ci)';
    end
    Cg7 = single(M_g / A_g);
    clear Yc M_g

    % ---- per-window fits ----------------------------------------------------------
    Cw7   = cell(1, nWin);
    nFall = 0;
    for k = 1:nWin
        idxW = find(trainF & (1:nTseg) >= wS(k) & (1:nTseg) <= wE(k));
        if numel(idxW) < g7SplitMinR * (nRegF+1)
            Cw7{k} = Cg7;  nFall = nFall + 1;  continue        % too few frames: use global
        end
        A_w = zeros(nRegF+1);  M_w = zeros(npx, nRegF+1);
        for c0 = 1:g7FitChunk:numel(idxW)
            jj = c0:min(c0+g7FitChunk-1, numel(idxW));  ci = idxW(jj);
            Yc = tovec(mov_res(:,:,ci));
            M_w = M_w + Yc * Iall7(:,ci)';
            A_w = A_w + Iall7(:,ci) * Iall7(:,ci)';
        end
        if rcond(A_w) < 1e-12
            Cw7{k} = Cg7;  nFall = nFall + 1;
        else
            Cw7{k} = single(M_w / A_w);
        end
        clear Yc M_w
    end
    if nFall > 0
        fprintf(['     %d of %d windows fell back to the global coefficients (too few ' ...
                 'training frames or singular design)\n'], nFall, nWin);
    end
    fprintf('     training frames per regressor: %.0f global, %.0f per window\n', ...
        nTr/(nRegF+1), median(arrayfun(@(k) nnz(trainF & (1:nTseg)>=wS(k) & (1:nTseg)<=wE(k)), ...
        1:nWin))/(nRegF+1));

    % ---- score both on the held-out frames -----------------------------------------
    teI = find(testF);
    s1R = zeros(npx,1); s2R = zeros(npx,1);
    s1G = zeros(npx,1); s2G = zeros(npx,1);
    s1S = zeros(npx,1); s2S = zeros(npx,1);
    for c0 = 1:g7FitChunk:numel(teI)
        jj = c0:min(c0+g7FitChunk-1, numel(teI));  ci = teI(jj);
        Yc = single(tovec(mov_res(:,:,ci)));
        Ic = single(Iall7(:,ci));
        Rg = Yc - Cg7 * Ic;
        Rs = Yc;
        for k = 1:nWin
            act = find(Wgt(k, ci) > 0);
            if isempty(act), continue; end
            Rs(:,act) = Rs(:,act) - Cw7{k} * (Ic(:,act) .* single(Wgt(k, ci(act))));
        end
        s1R = s1R + sum(double(Yc),2);  s2R = s2R + sum(double(Yc).^2,2);
        s1G = s1G + sum(double(Rg),2);  s2G = s2G + sum(double(Rg).^2,2);
        s1S = s1S + sum(double(Rs),2);  s2S = s2S + sum(double(Rs).^2,2);
    end
    clear Yc Ic Rg Rs
    vTe = @(s1,s2) (s2 - s1.^2/nTe) / max(nTe-1,1);
    vR7 = reshape(vTe(s1R,s2R), nrow, ncol);
    vG7 = reshape(vTe(s1G,s2G), nrow, ncol);
    vS7 = reshape(vTe(s1S,s2S), nrow, ncol);
    dropG = 1 - vG7 ./ max(vR7, eps);          % held-out variance removed, global
    dropS = 1 - vS7 ./ max(vR7, eps);          % ...and split

    mskSet = {bvBW, 'bvMask'; densBW, 'density mask'; neuronMap>0 & ~bvBW, 'dendrite (no vessel)'};
    fprintf('     HELD-OUT variance removed (%%):   global   split    split - global\n');
    for q = 1:3
        mq = mskSet{q,1};
        fprintf('       %-22s %7.2f %7.2f %+12.2f\n', mskSet{q,2}, ...
            100*mean(dropG(mq)), 100*mean(dropS(mq)), 100*mean(dropS(mq)-dropG(mq)));
    end
    gainVes = 100*mean(dropS(bvBW) - dropG(bvBW));
    gainDen = 100*mean(dropS(neuronMap>0 & ~bvBW) - dropG(neuronMap>0 & ~bvBW));
    if gainVes > 0.5 && gainDen < 0.5*gainVes
        fprintf(['     -> SPLIT WINS: +%.2f%% more removed on the vessels out of sample ' ...
                 'while the dendrite only\n        moves +%.2f%%. The coefficient really ' ...
                 'is time-varying.\n'], gainVes, gainDen);
    elseif gainVes <= 0
        fprintf(['     -> GLOBAL WINS: the windowed fit removes no more on the vessels ' ...
                 'out of sample (%+.2f%%),\n        so the extra parameters are only ' ...
                 'fitting noise. Keep the single fit.\n'], gainVes);
    else
        fprintf(['     -> AMBIGUOUS: +%.2f%% on the vessels but +%.2f%% on the dendrite. ' ...
                 'The windowed fit is\n        buying vessel variance partly with neural ' ...
                 'signal -- widen g7SplitWin.\n'], gainVes, gainDen);
    end

    % ---- how much does beta actually move between windows? -------------------------
    % Mean coefficient over the vessel pixels, per regressor, per window. A flat line is
    % a constant coefficient (and then the global fit was right all along).
    bvIdx = find(bvBW);
    betaW = zeros(nWin, nRegF);
    for k = 1:nWin, betaW(k,:) = mean(double(Cw7{k}(bvIdx, 2:end)), 1); end
    betaG = mean(double(Cg7(bvIdx, 2:end)), 1);
    cvBeta = std(betaW, 0, 1) ./ max(abs(mean(betaW,1)), eps);
    fprintf('     beta drift over windows (SD / |mean|, on bvMask): median %.2f, max %.2f\n', ...
        median(cvBeta), max(cvBeta));

    % ---- kymographs: the two fits and what each one leaves behind -------------------
    % The maps below say WHERE the two fits differ but not WHEN, which is the whole point
    % of a time-varying coefficient. These six panels put both fits on the same time axis:
    %
    %   raw                 the dendritic ROI traces going in
    %   fitted GLOBAL       one coefficient for the whole record
    %   fitted SPLIT        the blended per-window coefficients
    %   corrected GLOBAL    raw - fitted, the trace you would keep
    %   corrected SPLIT     the same with the windowed fit
    %   split - global      where and when the two fits actually disagree
    %
    % Read the two CORRECTED rows against each other: if the split fit is genuinely
    % better, its residual vessel texture is fainter over the bvMask rows and no fainter
    % on the dendrite rows. If the last row lights up mainly at window boundaries, or
    % looks like noise rather than vessel structure, the extra parameters are fitting
    % noise -- which is what the held-out numbers above are testing formally.
    %
    % The held-out blocks are marked in black along the top panel. ONLY those frames
    % count for the verdict; everywhere else both fits saw the data during training and
    % the split one is expected to look better whether or not it really is.
    if g7SplitKymo && exist('ftprnt','var') && exist('Pr','var')
        sdS = Pr.sortDist;
        [FvS, nPxS] = get_roiOp(ftprnt, sdS, g7FpMode, g7FpThr);   % same extraction as 6c
        nRs = size(FvS,2);
        fprintf('     building split-fit kymographs (%d ROIs, ''%s'', %d px each) ...\n', ...
            nRs, g7FpMode, round(median(nPxS)));

        TrawS = -(FvS' * tovec(mov_res));                  % spike-positive
        TgS   = -(FvS' * double(Cg7)) * Iall7;             % fitted, global
        TsS   = zeros(nRs, nTseg);                         % fitted, split (blended)
        for k = 1:nWin
            sp = wS(k):wE(k);
            Bk = -(FvS' * double(Cw7{k}));
            TsS(:,sp) = TsS(:,sp) + (Bk * Iall7(:,sp)) .* Wgt(k,sp);
        end

        dfS = @(T) T - movmedian(T, g7KymoBase, 2, 'omitnan');
        dR  = dfS(TrawS);
        scS = 1 ./ max(std(dR(:,gAll) - movmean(dR(:,gAll), Wsm, 2), 0, 2), eps);
        pS  = {dR .* scS,                       'raw'; ...
               dfS(TgS) .* scS,                 'fitted GLOBAL'; ...
               dfS(TsS) .* scS,                 'fitted SPLIT'; ...
               dfS(TrawS - TgS) .* scS,         'corrected GLOBAL'; ...
               dfS(TrawS - TsS) .* scS,         'corrected SPLIT'; ...
               dfS(TsS - TgS) .* scS,           'split - global (fits)'};
        clS = prctile(abs(pS{1,1}(:)), 99.5);
        if ~isfinite(clS) || clS <= 0, clS = 1; end

        if isfield(Pr,'dendaxis') && numel(Pr.dendaxis) == nRs
            daxS = Pr.dendaxis(:);
        else
            daxS = (1:nRs)';
        end
        tickS = unique(round(linspace(1, nRs, min(8, nRs))));
        tS    = (1:nTseg) / fs;

        get_figure('R3-flowSplitKymo');
        tlSk = tiledlayout(6, 1, 'TileSpacing','compact','Padding','compact');
        axS  = gobjects(1,6);
        for q = 1:6
            axS(q) = nexttile;
            imagesc(axS(q), tS, 1:nRs, pS{q,1}, [-clS clS]);
            colormap(axS(q), cmapCorr); hold(axS(q),'on');
            set(axS(q), 'YTick', tickS, 'YTickLabel', round(daxS(tickS)));
            ylabel(axS(q), {pS{q,2}, '\mum'});
            if q == 6, xlabel(axS(q), 'time (s)'); end
        end
        te0 = find(diff([false testMask false]) ==  1);      % held-out blocks, on top
        te1 = find(diff([false testMask false]) == -1) - 1;
        for b = 1:numel(te0)
            plot(axS(1), [te0(b) te1(b)]/fs, [0.9 0.9], 'k-', 'LineWidth', 4);
        end
        cbS = colorbar(axS(end));  cbS.Layout.Tile = 'east';
        cbS.Label.String = 'raw fast-noise SD';
        linkaxes(axS, 'x');  xlim(axS(1), [tS(1) tS(end)]);
        title(tlSk, sprintf(['f=%d  global vs %d-window split fit, same colour scale ' ...
            '(black bars on top = held-out blocks, the only frames the verdict uses)'], ...
            f, nWin));
        clear TrawS TgS TsS dR pS FvS FkS
    end

    get_figure('R3-flowSplitFit');
    tlSp = tiledlayout(2, 3, 'TileSpacing','compact','Padding','compact');
    mSp = {dropG, 'held-out removed, GLOBAL', [0 prctile(dropG(:),99)]; ...
           dropS, 'held-out removed, SPLIT',  [0 prctile(dropG(:),99)]; ...
           dropS-dropG, 'split - global', prctile(abs(dropS(:)-dropG(:)),99)*[-1 1]};
    for q = 1:3
        axs = nexttile;  cl = mSp{q,3};
        if ~all(isfinite(cl)) || cl(2) <= cl(1), cl = [-1 1]; end
        imagesc(axs, mSp{q,1}, cl); axis(axs,'image','off'); hold(axs,'on');
        if q < 3, colormap(axs, turbo); else, colormap(axs, cmapCorr); end
        colorbar(axs);
        contour(axs, bvBW,   [0.5 0.5], 'w', 'LineWidth', 0.8);
        contour(axs, densBW, [0.5 0.5], 'g', 'LineWidth', 0.8);
        title(axs, mSp{q,2});
    end

    axbw = nexttile([1 2]); hold(axbw,'on');
    tW8c = (wS + wE)/2 / fs;
    for rr = 1:nRegF
        plot(axbw, tW8c, betaW(:,rr), '-o', 'MarkerSize', 3, 'LineWidth', 0.9);
        yline(axbw, betaG(rr), ':', 'Color', [0.6 0.6 0.6], 'HandleVisibility','off');
    end
    box(axbw,'off'); xlabel(axbw, 'window centre (s)');
    ylabel(axbw, 'mean \beta over bvMask');
    title(axbw, sprintf(['per-window coefficient (dotted = the single global fit). ' ...
        'drift SD/|mean| median %.2f'], median(cvBeta)));

    axbr = nexttile; hold(axbr,'on');
    bar(axbr, [100*mean(dropG(bvBW)) 100*mean(dropS(bvBW)); ...
               100*mean(dropG(densBW)) 100*mean(dropS(densBW)); ...
               100*mean(dropG(neuronMap>0 & ~bvBW)) 100*mean(dropS(neuronMap>0 & ~bvBW))]);
    set(axbr, 'XTick', 1:3, 'XTickLabel', {'bvMask','density','dendrite'});
    ylabel(axbr, '% held-out variance removed'); legend(axbr, {'global','split'}, ...
        'Location','northeast'); box(axbr,'off');
    title(axbr, 'out-of-sample, the only fair comparison');
    title(tlSp, sprintf(['f=%d  one global coefficient vs %d windows of %d frames ' ...
        '(%d overlap), scored on %d held-out frames'], f, nWin, g7SplitWin, g7SplitOvl, nTe));

    clear s1R s2R s1G s2G s1S s2S A_g A_w
end

% ---- 5b. did re-deriving the footprint change it? -----------------------------------
% Leading eigenvector painted on the dendrite, (7f)'s from 2000 frames next to the
% full-record one, with their overlap |u_new . u_old| per segment. Two modes that differ
% only by the arbitrary SVD sign are the same footprint, so the sign is normalised first
% and the overlap is taken in absolute value.
if g7RefitPCA
    get_figure('R3-flowFootprintRefit');
    tlF = tiledlayout(3, 1, 'TileSpacing','compact','Padding','compact');
    uSet = {segU, sprintf('(7f) leading PC, %d frames', numel(trkFrames)); ...
            segU7, sprintf('(7g) leading PC, %d frames', nTseg)};
    for q = 1:2
        axu = nexttile;  eMap = nan(nrow, ncol);
        for k = 1:nSegU
            if isempty(uSet{q,1}{k}), continue; end
            u = uSet{q,1}{k}(:,1);  u = u / max(abs(u));
            if sum(u) < 0, u = -u; end
            eMap(segPx{useSeg(k)}) = u;
        end
        hU = imagesc(axu, eMap, [-1 1]);  set(hU, 'AlphaData', ~isnan(eMap));
        axis(axu,'image','off'); colormap(axu, cmapCorr); colorbar(axu); hold(axu,'on');
        contour(axu, bvBW, [0.5 0.5], 'k', 'LineWidth', 0.8);
        title(axu, uSet{q,2});
    end
    axb = nexttile;
    bar(axb, 1:nSegU, segAlign, 'FaceColor', [0.4 0.6 0.85]); hold(axb,'on');
    yline(axb, 0.8, 'r--', 'refitting changed the footprint below here');
    set(axb, 'XTick', 1:nSegU, 'XTickLabel', arrayfun(@(s) sprintf('s%d', s), ...
        useSeg, 'uni', 0)); ylim(axb, [0 1]); box(axb,'off');
    ylabel(axb, '|u_{new} \cdot u_{old}|');
    title(axb, 'overlap of the re-derived leading mode with the (7f) one');
    title(tlF, sprintf(['f=%d  the PCA footprint re-derived on the full record ' ...
        '(segments themselves unchanged)'], f));
end

% ---- 6. the regressor over the whole movie -------------------------------------------
get_figure('R3-flowRegressorFull');
tlG = tiledlayout(3 + min(segMaxShow, nSegU), 1, 'TileSpacing','compact','Padding','compact');

axg = nexttile;
imagesc(axg, (1:nTseg)/fs, 1:nRegF, flowRegF, [-3 3]);
colormap(axg, cmapCorr); colorbar(axg);
shortNm = arrayfun(@(r) sprintf('s%d', useSeg(regSegF(r))), 1:nRegF, 'uni', 0);
set(axg, 'YTick', 1:nRegF, 'YTickLabel', shortNm, 'FontSize', 7);
xlabel(axg, 'time (s)');
title(axg, sprintf('all %d regressors over the whole %.1f s segment (z-scored)', ...
    nRegF, nTseg/fs));

axd = nexttile; hold(axd,'on');
for k = 1:nSegU
    plot(axd, ((1:nBlk)-0.5)*blkLen/fs, fPkBlk(k,:), '-o', 'MarkerSize', 3, ...
        'DisplayName', sprintf('seg %d', useSeg(k)));
end
xline(axd, numel(trkFrames)/fs, 'k--', 'end of tracked window');
box(axd,'off'); ylabel(axd, 'PSD peak (Hz)'); legend(axd,'Location','eastoutside');
title(axd, 'per-block spectral peak: does the band the tracking set still fit later?');

axr = nexttile; hold(axr,'on');
for k = 1:nSegU
    plot(axr, ((1:nBlk)-0.5)*blkLen/fs, r2Blk(k,:), '-o', 'MarkerSize', 3);
end
xline(axr, numel(trkFrames)/fs, 'k--');
if exist('r2osc','var'), yline(axr, median(r2osc), 'r:', 'median R^2 in (7f)'); end
box(axr,'off'); ylabel(axr, 'R^2 per block'); xlabel(axr, 'time (s)');
title(axr, 'variance of the segment signal the full-record oscillator explains, per block');

tG = 1:min(1200, nTseg);
for k = 1:min(segMaxShow, nSegU)
    axo = nexttile; hold(axo,'on');
    yy = segYfull(k,tG);   yy = (yy-mean(yy))/max(std(yy),eps);
    cc = oscCosF(k,tG);    cc = (cc-mean(cc))/max(std(cc),eps);
    plot(axo, tG/fs, yy, 'k-', 'LineWidth', 0.8, 'DisplayName','measured S in segment');
    plot(axo, tG/fs, cc, 'r-', 'LineWidth', 1.0, 'DisplayName','A(t)cos\phi(t), full record');
    box(axo,'off'); xlim(axo, [tG(1) tG(end)]/fs);
    if k == 1, legend(axo,'Location','northeast'); end
    ylabel(axo, sprintf('seg %d', useSeg(k)));
    title(axo, sprintf('band %.0f-%.0f Hz (fPeak %.0f Hz, was %.0f Hz in (7f))', ...
        segBandF(k,1), segBandF(k,2), segFPeakF(k), segStat(k,6)));
    if k == min(segMaxShow, nSegU), xlabel(axo, 'time (s)'); end
end
title(tlG, sprintf(['f=%d  flow regressor on the WHOLE segment (%.1f s): same segment ' ...
    'footprints and spatial modes as (7f), time courses re-derived'], f, nTseg/fs));

get_figure('R3-flowFullVarianceMaps');
tlV = tiledlayout(1, 3, 'TileSpacing','compact','Padding','compact');
vS7 = {varIn7, 'var(mov_res), all good frames'; ...
       varFit7, 'var explained by the flow regressor'; ...
       dropF7, 'fraction removed'};
for q = 1:3
    axv = nexttile;  vv = vS7{q,1};
    imagesc(axv, vv, [0 prctile(vv(:), 99)]); axis(axv,'image','off'); hold(axv,'on');
    colormap(axv, turbo); colorbar(axv);
    contour(axv, bvBW,  [0.5 0.5], 'w', 'LineWidth', 0.8);
    contour(axv, densBW,[0.5 0.5], 'g', 'LineWidth', 0.8);
    title(axv, vS7{q,2});
end
title(tlV, sprintf(['f=%d  full-record fit (white = bvMask, green = density mask): ' ...
    '%.1f%% removed on bvMask, %.1f%% on the dendrite'], f, ...
    100*mean(dropF7(bvBW)), 100*mean(dropF7(neuronMap>0 & ~bvBW))));

% ---- 6c. kymograph of the dendritic ROI traces: before | blood | after --------------
% The traces the paper actually reports, laid out as space x time: one row per ROI of
% Result.ftprnt, ordered along the dendrite by Pr.sortDist and labelled in um by
% Pr.dendaxis. Three columns -- mov_res, the fitted vascular component, and the
% difference -- so the middle panel is exactly what leaves the left one.
%
% As in (7h) no corrected movie is needed: the correction is linear, so
%       T_blood = -(F' * Cfv7) * I     and     T_corr = T_raw - T_blood.
% Footprints are clipped at zero and unit-sum normalised (the convention used for the
% weighted maps elsewhere in this script) so ROIs of different size are comparable.
%
% Every panel in a row shares one colour scale, taken from the RAW panel -- a per-panel
% scale would renormalise the corrected trace and hide exactly the change being looked
% for. Traces are divided by each ROI's own fast-band noise measured on the RAW trace,
% so the units are "SD of that ROI's raw fast noise" and one clim suits every ROI.
if g7Kymo && exist('ftprnt','var') && exist('Pr','var')
    kFr = g7KymoFrames;  if isempty(kFr), kFr = 1:nTseg; end
    sd7 = Pr.sortDist;
    [Fkd, nPxK] = get_roiOp(ftprnt, sd7, g7FpMode, g7FpThr);   % px x nROI
    nRk = size(Fkd,2);
    fprintf(['[R3] kymograph: %d dendritic ROIs, ''%s'' extraction over %d px per ROI ' ...
             '(median), %d frames ...\n'], nRk, g7FpMode, round(median(nPxK)), nTseg);
    Ifull7 = [ones(1,nTseg); flowRegF - mean(flowRegF(:,gAll), 2)];
    Traw7  = -(Fkd' * tovec(mov_res));                % nROI x nTseg, spike-positive
    Tbld7  = -(Fkd' * Cfv7) * Ifull7;                 % the fitted vascular component
    Tcor7  = Traw7 - Tbld7;

    dfk  = @(T) T - movmedian(T, g7KymoBase, 2, 'omitnan');
    dRaw = dfk(Traw7);  dBld = dfk(Tbld7);  dCor = dfk(Tcor7);
    fRawK = dRaw - movmean(dRaw, Wsm, 2);            % fast component of the RAW trace
    scK   = 1 ./ max(std(fRawK(:,gAll), 0, 2), eps); % ...its noise sets the display scale

    % Rows are ROI INDEX in dendritic order, not micrometres: ROI positions along the
    % dendrite are not uniformly spaced, and imagesc would stretch them onto a uniform
    % axis using only the first and last value. The position is put on the tick labels
    % instead, so the numbers stay true.
    hasDax = isfield(Pr,'dendaxis') && numel(Pr.dendaxis) == nRk;
    if hasDax
        daxK = Pr.dendaxis(:);  yLab = 'dendritic position (\mum)';
    else
        daxK = (1:nRk)';        yLab = 'ROI (dendritic order)';
    end
    tickK = unique(round(linspace(1, nRk, min(8, nRk))));
    spkK  = find(full(Result.spike(1, segFrames(kFr))) > 0);

    get_figure('R3-flowKymograph');
    nBk = numel(g7KymoBands);
    tlK = tiledlayout(nBk, 3, 'TileSpacing','compact','Padding','compact');
    for b = 1:nBk
        switch g7KymoBands{b}
            case 'full', pick = @(D) D;                            bNm = 'full band';
            case 'fast', pick = @(D) D - movmean(D, Wsm, 2);        bNm = sprintf('fast (<%d ms)', Wsm);
            case 'sub',  pick = @(D) movmean(D, Wsm, 2);            bNm = sprintf('%d ms subthreshold', Wsm);
            otherwise,   error('g7KymoBands entries must be ''full'', ''fast'' or ''sub''.');
        end
        Xr = pick(dRaw(:,kFr)) .* scK;
        Xb = pick(dBld(:,kFr)) .* scK;
        Xc = pick(dCor(:,kFr)) .* scK;
        clK = prctile(abs(Xr(:)), 99.5);   if ~isfinite(clK) || clK <= 0, clK = 1; end
        pSet = {Xr, 'mov\_res'; Xb, 'fitted vascular'; Xc, 'corrected (mov\_res - fitted)'};
        for q = 1:3
            axk = nexttile;
            imagesc(axk, kFr/fs, 1:nRk, pSet{q,1}, [-clK clK]);
            colormap(axk, cmapCorr);  hold(axk,'on');
            set(axk, 'YTick', tickK, 'YTickLabel', round(daxK(tickK)));
            if q == 1 && ~isempty(spkK)                 % soma spikes, as a rug on top
                plot(axk, kFr(spkK)/fs, ones(1,numel(spkK)), 'k|', ...
                    'MarkerSize', 4, 'LineWidth', 0.5);
            end
            if q == 3, colorbar(axk); end
            if b == nBk, xlabel(axk, 'time (s)'); end
            if q == 1,   ylabel(axk, {bNm, yLab}); end
            title(axk, pSet{q,2});
        end
    end
    title(tlK, sprintf(['f=%d  Result.ftprnt traces before and after the flow correction ' ...
        '(%d ROIs, colour = SD of each ROI''s raw fast noise, shared per row)'], f, nRk));

    remK  = 100*(1 - var(Tcor7(:,gAll),0,2) ./ max(var(Traw7(:,gAll),0,2), eps));
    fCorK = dCor - movmean(dCor, Wsm, 2);
    remKf = 100*(1 - var(fCorK(:,gAll),0,2) ./ max(var(fRawK(:,gAll),0,2), eps));
    fprintf(['     variance removed from the dendritic ROIs: full band median %.2f%% ' ...
             '(max %.2f%%), fast band median %.2f%% (max %.2f%%)\n'], ...
        median(remK), max(remK), median(remKf), max(remKf));
    clear Fk Fkd Ifull7 dRaw dBld dCor Xr Xb Xc fRawK fCorK
end

if g7Browse
    % Same statistics as (7f), measured on the full record. varPct is extra here: it is
    % the share of the movie variance THIS row's segment removes, which (7f) cannot know
    % because its fit had not run yet at the point the browser opens.
    segAmpF = median(sqrt(oscCosF.^2 + oscSinF.^2), 2);
    varSegF = nan(1, nSegU);
    for k = 1:nSegU
        rk = find(regSegF == k);
        if isempty(rk), continue; end
        varSegF(k) = 100 * mean(reshape(sum((Bf7(:,rk) * Sf7(rk,rk)) .* Bf7(:,rk), 2), ...
            nrow, ncol) ./ max(varIn7, eps), 'all');
    end
    statFlowF = struct( ...
        'f0Hz',     regF0F(:), ...
        'score',    regScoreF(:), ...
        'rVolt',    rVoltF(:), ...
        'voltPen',  wVoltF(:), ...
        'ampCov',   max(abs(covRegF), [], 1)', ...
        'ampEnv',   segAmpF(regSegF(:)), ...
        'bandLo',   segBandF(regSegF(:), 1), ...
        'bandHi',   segBandF(regSegF(:), 2), ...
        'fPeakDrift', drift(regSegF(:)), ...
        'varPctSeg', varSegF(regSegF(:))', ...
        'segment',  useSeg(regSegF(:))');
    browse_regressors(flowRegF, covRegF, [nrow ncol], flowRegNmF, struct( ...
        'preTrace', segYnv(regSegF,:), ...
        'score',    regScoreF, ...
        'scoreMin', oscScoreMin, ...
        'f0',       regF0F, ...
        'seg',      regSegF, ...
        'segMask',  {cellfun(@(k) segPx{k}, num2cell(useSeg), 'uni', 0)}, ...
        'refMask',  bvBW, ...
        'cmap',     cmapCorr, ...
        'tWin',     1:min(2000, nTseg), ...
        'stats',    statFlowF, ...
        'statXY',   {{'f0Hz','score'}}, ...
        'statC',    'rVolt', ...
        'statLog',  {{'ampCov','ampEnv'}}, ...
        'figName',  'R3-flowRegressorsFull'));
end

% ---- 7. side by side: mov_res | fitted blood | mov_res - fitted -----------------------
%   [ mov_res | fitted vascular component | corrected = mov_res - fitted ]
% All three on ONE colour scale (the input's 1-99th percentile), sign flipped to the
% spike-positive convention, so the panels are directly comparable: whatever leaves panel
% 1 for panel 2 is gone from panel 3.
nQ7 = min(g7MovFrm, floor(nTseg/g7MovBin));
q7  = 1:(nQ7*g7MovBin);
fprintf('\n[R3] 1x3 movie: %d frames (bin %d) -> %.2f GB displayed, ~%.2f GB peak while it builds\n', ...
    nQ7, g7MovBin, 3*npx*nQ7*4/2^30, 6*npx*nQ7*4/2^30);

IallQ = [ones(1,numel(q7)); flowRegF(:,q7) - mean(flowRegF(:,gAll), 2)];
qFit  = reshape(single(Cfv7) * single(IallQ), nrow, ncol, numel(q7));
qIn   = single(-mov_res(:,:,q7));
qFit  = -qFit;                                   % fitted vascular, spike-positive
qRes  = qIn - qFit;                              % corrected movie
if g7MovBin > 1
    rb   = @(A) squeeze(mean(reshape(A, nrow, ncol, g7MovBin, nQ7), 3));
    qIn  = rb(qIn);  qFit = rb(qFit);  qRes = rb(qRes);
end

clim7 = double(prctile(qIn(:), [1 99]));
sep7  = single(clim7(2)) * ones(nrow, g7MovSep, nQ7, 'single');
mov7  = cat(2, qIn, sep7, qFit, sep7, qRes);
get_figure('[R3] mov_res | fitted blood | corrected');
moviefixsc(mov7, clim7);
title(sprintf(['f=%d  [mov\\_res | fitted vascular | mov\\_res - fitted]   %d regressors ' ...
    'from %d segments, %d/%d frames, removes %.1f%% in the density mask / %.1f%% on the dendrite'], ...
    f, nRegF, nSegU, nQ7*g7MovBin, nTseg, ...
    100*mean(dropF7(densBW)), 100*mean(dropF7(neuronMap>0 & ~bvBW))));
clear qIn qFit qRes sep7 IallQ MtI Ig s1 s2

% ---- 8. the corrected movie, for whatever comes next ---------------------------------
% Materialised only on request: mov_res is already ~8 GB here, so a second full-length
% copy is a deliberate choice, not a default. Build it in chunks when you want it:
%
%   mov_flow = zeros(nrow, ncol, nTseg, 'single');
%   for c0 = 1:500:nTseg
%       ci = c0:min(c0+499, nTseg);
%       Iq = [ones(1,numel(ci)); flowRegF(:,ci) - mean(flowRegF(:,gAll),2)];
%       mov_flow(:,:,ci) = single(mov_res(:,:,ci)) - ...
%           reshape(single(Cfv7)*single(Iq), nrow, ncol, numel(ci));
%   end


%% (7h) Traces before and after the flow correction: does it remove blood or signal?
% (7g) showed the correction as a movie and as a per-pixel variance map. The question
% that actually decides whether it is usable is what it does to the TRACES -- the
% quantities the paper reports -- so this cell extracts the same sub-ROI traces from
% mov_res and from the corrected movie and puts them side by side.
%
% NO CORRECTED MOVIE IS BUILT. The correction is a linear projection, so for any
% footprint F the corrected trace follows in closed form:
%       T_corr = -F' * (mov_res - C*I) = T_raw + (F'*C) * I
% with C = Cfv7 (pixels x 1+nReg) and I = [1; centred regressors]. F'*C is Msub x (1+nReg),
% so the whole thing costs one small matrix product instead of a second 8 GB movie.
% The odd and even checkerboard halves are projected the same way, which keeps the
% shot-noise-free odd/even correlation of cell (4c) available after correction.
%
% FOUR THINGS ARE CHECKED, in increasing order of how much they matter:
%   1. Where the removed variance sits.  It must concentrate on sub-ROIs that overlap
%      bvMask. A correction that removes as much off the vessels as on them is removing
%      something that is not blood.
%   2. WHICH FREQUENCIES it comes out of.  The regressor is a narrowband oscillator, so
%      the before/after PSDs should differ inside the segment bands and be on top of each
%      other elsewhere. A broadband drop means the design is soaking up general variance.
%   3. The spike-triggered average.  Spike frames were excluded from the fit, so the
%      BLOOD component's own STA should be flat. If it has a deflection at t=0 the
%      regressor has learned the AP and the corrected traces are biased.
%   4. The fast (<20 ms) correlation vs geodesic distance -- the reviewer's actual claim.
%      If the short-range fast correlation is vascular it should collapse after the
%      correction; if it survives, it is not blood.
%
% ------------------------------- parameters -------------------------------
g8BvFrac    = 0.2;      % a sub-ROI counts as "on a vessel" if this fraction of its
                        % pixels fall inside bvBW
g8BaseWin   = 1001;     % ms, display baseline (running median), as in cell (6)
g8PsdWin    = 1024;     % samples per PSD window (50% overlap, Hann)
g8StaWin    = [-30 30]; % ms around each spike for the STA
g8StaBase   = [-30 -10];% ms, STA baseline window
g8StaPeak   = [0 5];    % ms, where the spike amplitude is measured
g8BigROI    = true;     % also summarise the paper's dendritic ROIs (Result.ftprnt), and
                        % compare against the independent bleach-dr correction if its
                        % BVcorr_traces.mat is on disk
g8BleachBP  = [200 220];% Hz, the bp band the bleach-dr run used (from the driver's `bp`).
                        % Only shades the PSD panel -- set [] to leave it unshaded.
g8NormMode  = 'bap';    % display units for R3-traceCompare-bleach:
                        %  'bap'      ./ the bAP measured on each pipeline's own STA
                        %             -> 1.0 = one spike  (default)
                        %  'dff'      that, times RobustdFFfit.dFFrobustbAP -> genuine dF/F
                        %  'bap_file' the literal ./ScaleFactor ./dFFrobustbAP. WRONG for
                        %             these traces (see the note in the cell) -- kept only
                        %             to reproduce it
                        %  'own'      each trace by its own fast-band noise SD
g8ScaleFile = '';       % RobustdFFfit.mat. '' = <fpath{f}>/RobustdFFfit.mat
g8ExBigN    = 2;        % dendritic ROIs drawn in the bleach-dr comparison
g8ExShowMs  = 2000;     % ms plotted in the example-trace figure
% --------------------------------------------------------------------------

assert(exist('Cfv7','var') && exist('flowRegF','var'), 'Run (7g) first (it builds Cfv7/flowRegF).');
assert(exist('Tall','var') && exist('Fall','var') && exist('Fodd','var'), ...
       'Run (4c) first (it builds Fall/Fodd/Feven and Tall/Todd/Teven).');

% ---- 1. corrected traces, in closed form --------------------------------------------
Ifull = [ones(1,nTseg); flowRegF - mean(flowRegF(:,gAll), 2)];
prj8  = @(Fv) double(single(tovec(Fv))' * single(Cfv7)) * Ifull;   % ROIs x nTseg, mov_res units
bldAll  = -prj8(Fall);          % spike-positive vascular component of each trace
bldOdd  = -prj8(Fodd);
bldEven = -prj8(Feven);

TallC  = Tall  - bldAll;        % corrected = raw - blood, same sign convention as Tall
ToddC  = Todd  - bldOdd;
TevenC = Teven - bldEven;

Wsm8    = Wsm;                                   % 20 ms, the paper's window
fTallC  = TallC  - movmean(TallC , Wsm8, 2);     % fast (<20 ms) after correction
fToddC  = ToddC  - movmean(ToddC , Wsm8, 2);
fTevenC = TevenC - movmean(TevenC, Wsm8, 2);
sTall   = movmean(Tall , Wsm8, 2);               % 20 ms subthreshold, before
sTallC  = movmean(TallC, Wsm8, 2);               % and after

% ---- 2. which sub-ROIs sit on a vessel ----------------------------------------------
ftBin8   = ftprintsmallest > 0;
npx8     = squeeze(sum(sum(ftBin8, 1), 2));
bvFrac8  = squeeze(sum(sum(ftBin8 & bvBW,   1), 2)) ./ max(npx8, 1);
dsFrac8  = squeeze(sum(sum(ftBin8 & densBW, 1), 2)) ./ max(npx8, 1);
isVess8  = bvFrac8 >= g8BvFrac;
fprintf('\n[R3] (7h) %d of %d sub-ROIs have >= %.0f%% of their pixels on bvMask\n', ...
    nnz(isVess8), Msub, 100*g8BvFrac);

vR8 = var(Tall (:,gdfr), 0, 2);   vC8 = var(TallC (:,gdfr), 0, 2);
vFR = var(fTall(:,gdfr), 0, 2);   vFC = var(fTallC(:,gdfr), 0, 2);
vSR = var(sTall(:,gdfr), 0, 2);   vSC = var(sTallC(:,gdfr), 0, 2);
remAll  = 100*(1 - vC8./max(vR8,eps));
remFast = 100*(1 - vFC./max(vFR,eps));
remSub  = 100*(1 - vSC./max(vSR,eps));

fprintf('     variance removed from the sub-ROI traces (median over ROIs):\n');
fprintf('     %-22s %-12s %-12s %s\n', 'band', 'on vessel', 'off vessel', 'ratio');
prt8 = @(nm, v) fprintf('     %-22s %-12.1f %-12.1f %.1f\n', nm, ...
    median(v(isVess8)), median(v(~isVess8)), median(v(isVess8))/max(median(v(~isVess8)),eps));
prt8('full band  (%)',      remAll);
prt8('fast <20 ms  (%)',    remFast);
prt8('20 ms subthr  (%)',   remSub);
if median(remAll(~isVess8)) > 0.5*median(remAll(isVess8))
    warning(['The correction removes nearly as much off the vessels as on them ' ...
        '(%.1f%% vs %.1f%%) -- it is not behaving like a vascular regressor.'], ...
        median(remAll(~isVess8)), median(remAll(isVess8)));
end

% ---- 3. what the blood component looks like at the spikes ---------------------------
% The fit excluded peri-spike frames, so this STA is a genuine out-of-sample test.
spkSeg = full(Result.spike(1, segFrames)) > 0;
wSta   = g8StaWin(1):g8StaWin(2);
sIdx   = find(spkSeg);
sIdx   = sIdx(sIdx + wSta(1) >= 1 & sIdx + wSta(end) <= nTseg);
staRaw = nan(Msub, numel(wSta));  staCor = staRaw;  staBld = staRaw;
if numel(sIdx) >= 5
    idxM  = sIdx(:)' + wSta(:);                        % nW x nSpk
    idxV  = idxM(:);                                   % A(r, idxM) would flatten to a
    nWs   = numel(wSta);                               % row, so reshape it back by hand
    dFraw = Tall  - movmedian(Tall , g8BaseWin, 2, 'omitnan');
    dFcor = TallC - movmedian(TallC, g8BaseWin, 2, 'omitnan');
    bBase = wSta >= g8StaBase(1) & wSta <= g8StaBase(2);
    bPeak = wSta >= g8StaPeak(1) & wSta <= g8StaPeak(2);
    sta1  = @(X, r) mean(reshape(X(r, idxV), nWs, []), 2)';    % 1 x nW spike average
    for r = 1:Msub
        a = sta1(dFraw,  r);   staRaw(r,:) = a - mean(a(bBase));
        a = sta1(dFcor,  r);   staCor(r,:) = a - mean(a(bBase));
        a = sta1(bldAll, r);   staBld(r,:) = a - mean(a(bBase));
    end
    ampRaw = max(staRaw(:,bPeak), [], 2);   ampCor = max(staCor(:,bPeak), [], 2);
    ampBld = max(abs(staBld(:,bPeak)), [], 2);
    nseRaw = std(fTall (:,gdfr), 0, 2);     nseCor = std(fTallC(:,gdfr), 0, 2);
    fprintf('\n[R3] spike-triggered average over %d spikes (median over sub-ROIs):\n', numel(sIdx));
    fprintf('     spike amplitude   raw %.2f -> corrected %.2f counts (%.1f%% change)\n', ...
        median(ampRaw), median(ampCor), 100*(median(ampCor)/max(median(ampRaw),eps) - 1));
    fprintf('     fast-band noise   raw %.2f -> corrected %.2f counts (%.1f%% change)\n', ...
        median(nseRaw), median(nseCor), 100*(median(nseCor)/max(median(nseRaw),eps) - 1));
    fprintf('     spike SNR         raw %.2f -> corrected %.2f  (%+.1f%%)\n', ...
        median(ampRaw./nseRaw), median(ampCor./nseCor), ...
        100*(median(ampCor./nseCor)/max(median(ampRaw./nseRaw),eps) - 1));
    fprintf('     STA of the REMOVED component: %.2f counts at the spike vs %.2f amplitude\n', ...
        median(ampBld), median(ampRaw));
    if median(ampBld) > 0.15*median(ampRaw)
        warning(['The removed component itself has a %.0f%%-of-a-spike deflection in its ' ...
            'STA: the regressor is partly tracking the AP, not only the vessel.'], ...
            100*median(ampBld)/max(median(ampRaw),eps));
    end
else
    fprintf('\n[R3] only %d usable spikes in this segment -- STA check skipped\n', numel(sIdx));
end

% ---- 4. where in frequency did the variance go? -------------------------------------
% Block-averaged Hann periodogram. Bad frames are zeroed AFTER removing the window mean,
% and raw and corrected get identical treatment, so the comparison is fair even though
% the gating itself spreads a little power.
nW8   = min(g8PsdWin, 2^floor(log2(nTseg/2)));
hopW  = round(nW8/2);
wHan  = 0.5*(1 - cos(2*pi*(0:nW8-1)/(nW8-1)));
st8   = 1:hopW:(nTseg - nW8 + 1);
fvec8 = (0:floor(nW8/2)) * fs / nW8;
Pr8 = zeros(Msub, numel(fvec8));  Pc8 = Pr8;
for b = 1:numel(st8)
    ib = st8(b) + (0:nW8-1);
    for q = 1:2
        X = Tall(:,ib);  if q == 2, X = TallC(:,ib); end
        X = X - mean(X, 2);
        X(:, ~gdfr(ib)) = 0;
        Z = abs(fft(X .* wHan, [], 2)).^2;
        Z = Z(:, 1:numel(fvec8));
        if q == 1, Pr8 = Pr8 + Z; else, Pc8 = Pc8 + Z; end
    end
end
Pr8 = Pr8 / numel(st8);  Pc8 = Pc8 / numel(st8);

% ---- 5. fast correlation vs geodesic distance, before and after ---------------------
% tau = 0 slice of the (4c) tensor, recomputed the same way: corr(odd_i, even_j), which
% is shot-noise-free on the diagonal because odd and even share signal but not noise.
symc8 = @(O,E) 0.5*(corr(O(:,gdfr)', E(:,gdfr)') + corr(E(:,gdfr)', O(:,gdfr)'));
R0raw8 = symc8(fTodd,  fTeven);
R0cor8 = symc8(fToddC, fTevenC);
if exist('R0','var')
    fprintf('\n[R3] sanity: recomputed tau=0 matrix vs the (4c) one, max |difference| = %.3g\n', ...
        max(abs(R0raw8(triuM) - R0(triuM))));
end
dvec8  = Dgeo(triuM);
prRaw8 = R0raw8(triuM);  prCor8 = R0cor8(triuM);
if ~exist('distEdges','var'), distEdges = 0:2:60; end
distCtr8 = (distEdges(1:end-1) + distEdges(2:end)) / 2;
cRawB = nan(1,numel(distCtr8));  cCorB = cRawB;  nB8 = zeros(1,numel(distCtr8));
for b = 1:numel(distCtr8)
    sel = dvec8 >= distEdges(b) & dvec8 < distEdges(b+1);
    nB8(b) = nnz(sel);
    if any(sel), cRawB(b) = median(prRaw8(sel)); cCorB(b) = median(prCor8(sel)); end
end
autoRaw8 = mean(diag(R0raw8), 'omitnan');   autoCor8 = mean(diag(R0cor8), 'omitnan');
shortSel = dvec8 < 10;
fprintf(['\n[R3] fast (<%d ms) tau=0 correlation, raw -> corrected:\n' ...
         '     within a sub-ROI (odd/even auto) : %.3f -> %.3f\n' ...
         '     pairs < 10 um apart (median)     : %.3f -> %.3f\n' ...
         '     pairs > 30 um apart (median)     : %.3f -> %.3f\n'], Wsm8, ...
    autoRaw8, autoCor8, ...
    median(prRaw8(shortSel)), median(prCor8(shortSel)), ...
    median(prRaw8(dvec8>30)), median(prCor8(dvec8>30)));

% ---- 6. figures ---------------------------------------------------------------------
ftL8 = double(ftprintsmallest);  den8 = sum(ftL8,3);
[~, labImg8] = max(ftL8, [], 3);  inMask8 = den8 > eps;  labImg8(~inMask8) = 0;
mapAll8 = nan(nrow,ncol);  mapAll8(inMask8) = remAll (labImg8(inMask8));
mapFst8 = nan(nrow,ncol);  mapFst8(inMask8) = remFast(labImg8(inMask8));

% -- 6a. what the correction takes out of the example sub-ROIs --
% The RESIDUAL only: bldAll = mov_res - corrected, i.e. the fitted vascular component of
% each sub-ROI trace, full band. This is the quantity that has to look like blood -- an
% oscillation while the vessel is flowing and nothing at the spikes. Its size relative to
% the trace it came out of is in the panel title rather than as an extra curve.
tW8 = tShow(1:min(g8ExShowMs, numel(tShow)));
get_figure('R3-traceCompare-examples');
axe = axes; hold(axe,'on');  %#ok<LAXES>
Yd  = bldAll(exROI, tW8);
off = 1.2 * max(max(Yd,[],2) - min(Yd,[],2));
if ~isfinite(off) || off <= 0, off = 1; end
for k = 1:2
    yo = (2-k)*off;
    plot(axe, tW8-tW8(1), Yd(k,:) + yo, '-', 'Color', exColr(k,:), 'LineWidth', 0.9);
    text(axe, 0, yo + 0.45*off, sprintf('#%d (bvMask %.0f%%, %.1f%% of its variance)', ...
        exROI(k), 100*bvFrac8(exROI(k)), remAll(exROI(k))), ...
        'Color', exColr(k,:), 'FontSize', 8);
end
% (no spike markers here: tShow is by construction the longest spike-free, blue-off
%  stretch, so this window never contains one. The STA in the summary figure is where
%  the residual's behaviour at spikes is checked.)
box(axe,'off'); xlim(axe, [0 numel(tW8)-1]);
xlabel(axe, 'time (ms)');  ylabel(axe, 'residual: mov\_res - corrected (counts)');
title(axe, sprintf(['f=%d  what the flow correction removes from two sub-ROIs %.1f \\mum ' ...
    'apart'], f, Dgeo(exROI(1),exROI(2))));

% -- 6b. summary panels --
get_figure('R3-traceCompare-summary');
tlS8 = tiledlayout(2, 3, 'TileSpacing','compact','Padding','compact');

for q = 1:2
    axm = nexttile;  vv = remAll;  mp = mapAll8;  nm = 'full band';
    if q == 2, vv = remFast;  mp = mapFst8;  nm = sprintf('fast <%d ms', Wsm8); end
    hI = imagesc(axm, mp, [0 prctile(vv, 98)]); set(hI, 'AlphaData', inMask8);
    axis(axm,'image','off'); colormap(axm, turbo); colorbar(axm); hold(axm,'on');
    contour(axm, bvBW,   [0.5 0.5], 'w', 'LineWidth', 0.9);
    contour(axm, densBW, [0.5 0.5], 'g', 'LineWidth', 0.9);
    title(axm, sprintf('%% variance removed, %s', nm));
end

axsc = nexttile; hold(axsc,'on');
scatter(axsc, 100*bvFrac8, remFast, 18, dsFrac8, 'filled');
colormap(axsc, turbo); cb8 = colorbar(axsc); cb8.Label.String = 'density-mask overlap';
xlabel(axsc, '% of sub-ROI pixels on bvMask'); ylabel(axsc, '% fast variance removed');
box(axsc,'off');
rSc = corr(bvFrac8(:), remFast(:), 'rows','complete');
title(axsc, sprintf('removal follows the vessel mask (r = %.2f)', rSc));

axp = nexttile; hold(axp,'on');
for k = unique(regSegF(:))'                         % bands of the SURVIVING regressors
    patch(axp, segBandF(k,[1 2 2 1]), [1e-6 1e-6 1e12 1e12], [1 0.9 0.6], ...
        'EdgeColor','none', 'FaceAlpha', 0.45, 'HandleVisibility','off');
end
grp8 = {isVess8, 'on vessel'; ~isVess8, 'off vessel'};
col8 = [0.85 0.33 0.10; 0 0.45 0.74];
for g = 1:2
    if ~any(grp8{g,1}), continue; end
    plot(axp, fvec8, median(Pr8(grp8{g,1},:),1), '-',  'Color', col8(g,:), ...
        'LineWidth', 1.6, 'DisplayName', [grp8{g,2} ', mov\_res']);
    plot(axp, fvec8, median(Pc8(grp8{g,1},:),1), '--', 'Color', col8(g,:), ...
        'LineWidth', 1.2, 'DisplayName', [grp8{g,2} ', corrected']);
end
set(axp, 'XScale','log', 'YScale','log'); xlim(axp, [1 fs/2]); box(axp,'off');
xlabel(axp, 'frequency (Hz)'); ylabel(axp, 'power'); legend(axp, 'Location','southwest');
title(axp, 'PSD before / after (shaded = regressor bands)');

axd = nexttile; hold(axd,'on');
plot(axd, 0, autoRaw8, 'ko', 'MarkerFaceColor','k', 'DisplayName','within ROI, mov\_res');
plot(axd, 0, autoCor8, 'ks', 'MarkerFaceColor','w', 'DisplayName','within ROI, corrected');
plot(axd, distCtr8, cRawB, '-o', 'Color', [0.3 0.3 0.3], 'MarkerSize', 4, ...
    'LineWidth', 1.4, 'DisplayName','mov\_res');
plot(axd, distCtr8, cCorB, '-o', 'Color', [0.85 0.2 0.2], 'MarkerSize', 4, ...
    'LineWidth', 1.4, 'DisplayName','corrected');
yline(axd, 0, 'k:'); box(axd,'off'); legend(axd, 'Location','northeast');
xlabel(axd, 'geodesic distance (\mum)'); ylabel(axd, sprintf('\\tau=0 fast corr'));
title(axd, sprintf('fast (<%d ms) correlation vs distance', Wsm8));

axa = nexttile; hold(axa,'on');
if numel(sIdx) >= 5
    plot(axa, wSta, median(staRaw,1), '-',  'Color', [0.3 0.3 0.3], 'LineWidth', 1.6, ...
        'DisplayName','mov\_res');
    plot(axa, wSta, median(staCor,1), '-',  'Color', [0.85 0.2 0.2], 'LineWidth', 1.4, ...
        'DisplayName','corrected');
    plot(axa, wSta, median(staBld,1), '-',  'Color', [0 0.45 0.74], 'LineWidth', 1.4, ...
        'DisplayName','removed component');
    xline(axa, 0, 'k:'); yline(axa, 0, 'k:'); box(axa,'off');
    legend(axa, 'Location','northwest');
    title(axa, sprintf('STA over %d spikes: the removed part must be flat', numel(sIdx)));
else
    text(axa, 0.5, 0.5, 'too few spikes for an STA', 'HorizontalAlignment','center');
    axis(axa,'off');
end
xlabel(axa, 'time from spike (ms)'); ylabel(axa, 'dF (counts)');
title(tlS8, sprintf(['f=%d  what the flow correction does to the sub-ROI traces ' ...
    '(%d ROIs, %d on bvMask)'], f, Msub, nnz(isVess8)));

% (the full tau=0 correlation MATRIX before/after is not drawn -- its distance profile is
%  the panel in R3-traceCompare-summary, and the numbers are in the printout above.)

% ---- 7. the paper's dendritic ROIs: mov_res / flow-corrected / bleach-dr -------------
% The flow correction is not the only blood correction in this project. The whole-movie
% bleach-dr run (Analysis_20260806_BloodVesselUtkuAlgorithm) wrote BVcorr_traces.mat,
% which is an INDEPENDENT correction of the same data -- different algorithm, different
% band (its bp parameter, ~200-220 Hz), applied over the whole recording. Putting the
% three traces side by side is the strongest available check: two methods that both
% remove blood should agree about what is left, even though nothing about them is shared.
%
% THREE CAVEATS, all of which the code respects:
%   1. bleach-dr exists ONLY for the dendritic ROIs (Result.ftprnt), not for the ~30 px
%      sub-ROIs, so this comparison lives here and not in the panels above.
%   2. Its input movie is loadVoltChunk's photobleach-normalised residual, NOT this
%      script's mov_res (mean image out, mcTrace regressed). The amplitude scale is
%      therefore different and absolute counts are NOT comparable. Everything reported
%      below is scale-free -- variance ratios, correlations, PSD shape -- and each trace
%      is displayed after division by its OWN fast-band noise.
%   3. The driver projects onto the UNNORMALISED tovec(Result.ftprnt). The same
%      unnormalised footprints are used here, so the only difference between the traces
%      is the movie behind them, not the weighting.
if g8BigROI && exist('ftprnt','var')
    Fbg   = tovec(double(ftprnt));                 % driver convention: no normalisation
    Tbig  = -(Fbg' * tovec(mov_res));              % nROI x nTseg, spike-positive
    bldB  = -(Fbg' * Cfv7) * Ifull;
    TbigC = Tbig - bldB;

    bl8   = struct('have', false);
    bvtrF = fullfile(fpath{f}, 'BVcorr_traces.mat');
    if isfile(bvtrF)
        % the file covers the whole 358000-frame recording; read only segFrames out of it
        try
            mf8 = matfile(bvtrF);   sz8 = size(mf8, 'tr_corrected');
            getBl = @() mf8.tr_corrected(segFrames, :)';
        catch
            Sb8 = load(bvtrF, 'tr_corrected');   sz8 = size(Sb8.tr_corrected);
            getBl = @() Sb8.tr_corrected(segFrames, :)';
        end
        if sz8(1) >= max(segFrames) && sz8(2) == size(Tbig,1)
            bl8.T = getBl();            % nROI x nTseg, already spike-positive
            bl8.have = true;
        else
            warning(['BVcorr_traces.mat is %d x %d but segFrames needs %d frames and ' ...
                '%d ROIs -- skipping the bleach-dr comparison.'], ...
                sz8(1), sz8(2), max(segFrames), size(Tbig,1));
        end
        clear Sb8 mf8
    else
        fprintf(['\n[R3] BVcorr_traces.mat not found in %s -- showing mov_res vs the flow ' ...
                 'correction only\n'], fpath{f});
    end

    % display order along the dendrite, and per-method band decomposition
    sd8  = Pr.sortDist;
    if isfield(Pr,'dendaxis') && numel(Pr.dendaxis) == numel(sd8)
        dax8 = Pr.dendaxis(:);
    else
        dax8 = (1:numel(sd8))';        % no position available: fall back to ROI index
    end
    dfb  = @(T) T - movmedian(T, g8BaseWin, 2, 'omitnan');
    mset = {dfb(Tbig(sd8,:)), 'mov\_res', [0.35 0.35 0.35]; ...
            dfb(TbigC(sd8,:)), 'flow-corrected', [0.85 0.20 0.20]};
    if bl8.have
        mset(end+1,:) = {dfb(bl8.T(sd8,:)), 'bleach-dr', [0 0.45 0.74]};
    end
    nM8  = size(mset,1);
    nRb8 = size(mset{1,1}, 1);
    fst8 = cell(1,nM8);  sub8 = cell(1,nM8);
    for m = 1:nM8
        fst8{m} = mset{m,1} - movmean(mset{m,1}, Wsm8, 2);
        sub8{m} = movmean(mset{m,1}, Wsm8, 2);
    end

    % ---- display normalisation: put every method in units of ONE bAP ---------------
    % RobustdFFfit defines  ScaleFactor = (raw STA bAP peak) / (robust dF/F of the bAP),
    % applied as normTraces ./ ScaleFactor = genuine dF/F.
    %
    % WHY ScaleFactor IS NOT APPLIED TO THESE TRACES DIRECTLY.
    % It was fitted on Result.normTraces, which are ALREADY normalised: on this dataset
    % their median |value| is 0.02 and their measured raw STA bAP peak is 0.1066, against
    % ScaleFactor.*dFFrobustbAP = 0.1042 -- ratio 1.013, i.e. the file is exact for that
    % trace and that trace only. The traces here are projections onto the UNNORMALISED
    % tovec(Result.ftprnt), kept because that is what the bleach-dr driver used and the
    % comparison has to be like-for-like. Its per-ROI weight sum is ~3.6e8, so a 100-count
    % dF reads ~3.6e10; dividing THAT by 0.104 gives ~1e11, a number with no meaning.
    % Cross-pipeline the mismatch is worse still: bleach-dr ran on loadVoltChunk's
    % photobleach-normalised movie, so it carries yet another gain.
    %
    % The bridge that does not depend on the extraction is the bAP itself. Its peak is
    % measured on each trace's own spike-triggered average, and then
    %       T ./ bapMeas                  -> units of one bAP   (1.0 = a spike)
    %       T ./ bapMeas .* dFFrobustbAP  -> genuine dF/F, calibrated by the robust dF/F
    %                                        of the bAP from the same file
    % Both use the file for the quantity it can legitimately supply (dFFrobustbAP, a
    % physical dF/F independent of how the trace was extracted); ScaleFactor is still
    % loaded, but only to VALIDATE the scales, printed below.
    %
    % mov_res and flow-corrected SHARE one divisor (mov_res's own bAP) because they come
    % from the same movie -- giving each its own would rescale the corrected trace and
    % hide what the correction did. bleach-dr, a different pipeline, gets its own.
    %   g8NormMode = 'bap'      1.0 = one bAP  (default)
    %                'dff'      genuine dF/F via dFFrobustbAP
    %                'bap_file' the literal ./ScaleFactor ./dFFrobustbAP -- WRONG for
    %                           these traces, kept only so the 1e11 can be reproduced
    %                'own'      each trace by its own fast-band noise SD
    sfFile = g8ScaleFile;
    if isempty(sfFile), sfFile = fullfile(fpath{f}, 'RobustdFFfit.mat'); end
    sfS = [];  bapS = [];
    needSF = ismember(g8NormMode, {'dff','bap_file'});   % 'bap' only uses it to validate
    if needSF || isfile(sfFile)
        assert(isfile(sfFile), ['RobustdFFfit.mat not found at %s, needed for ' ...
            'g8NormMode = ''%s''. Set g8ScaleFile or use ''bap''.'], sfFile, g8NormMode);
        Rdff = importdata(sfFile);
        % ScaleFactor / dFFrobustbAP are indexed in Result.ftprnt order (ALL ROIs), while
        % sortDist is the dendritic order of the GOOD ROIs only -- a subset, since
        % Pr.sortDist drops BadROI. So the length check is against the full ROI count and
        % the subsetting happens in the same step as the reordering.
        nROIall8 = size(ftprnt, 3);
        assert(numel(Rdff.ScaleFactor) == nROIall8 && numel(Rdff.dFFrobustbAP) == nROIall8, ...
            ['RobustdFFfit has %d ScaleFactor / %d dFFrobustbAP entries but Result.ftprnt ' ...
             'has %d ROIs -- it was fitted on a different ROI set.'], ...
            numel(Rdff.ScaleFactor), numel(Rdff.dFFrobustbAP), nROIall8);
        assert(max(sd8) <= nROIall8, 'Pr.sortDist indexes ROI %d but there are only %d.', ...
            max(sd8), nROIall8);
        sfS  = Rdff.ScaleFactor(sd8);  sfS  = sfS(:);     % Result order -> dendritic order
        bapS = Rdff.dFFrobustbAP(sd8); bapS = bapS(:);
    end

    % bAP peak measured on each method's own STA, for the consistency check
    bapMeas = nan(nRb8, nM8);
    if numel(sIdx) >= 5
        idxV8 = reshape(sIdx(:)' + wSta(:), [], 1);   nW8s = numel(wSta);
        bB8   = wSta >= g8StaBase(1) & wSta <= g8StaBase(2);
        bP8   = wSta >= g8StaPeak(1) & wSta <= g8StaPeak(2);
        for m = 1:nM8
            for r = 1:nRb8
                a = mean(reshape(mset{m,1}(r, idxV8), nW8s, []), 2)';
                a = a - mean(a(bB8));
                bapMeas(r,m) = max(a(bP8));
            end
        end
    end

    % divisor per method, with a floor so an ROI whose bAP is ~0 (or negative, far out on
    % a dendrite) cannot blow the display up
    divM = nan(nRb8, nM8);  nFloor = 0;
    if numel(sIdx) >= 5
        bFloor = 0.05 * median(bapMeas(:,1), 'omitnan');
        dRawM  = max(bapMeas(:,1), bFloor);            % mov_res AND flow-corrected
        nFloor = nnz(bapMeas(:,1) < bFloor);
        divM(:,1:min(2,nM8)) = repmat(dRawM, 1, min(2,nM8));
        if bl8.have
            divM(:,3) = max(bapMeas(:,3), 0.05*median(bapMeas(:,3), 'omitnan'));
        end
    end

    scM8 = cell(1,nM8);
    switch g8NormMode
        case 'own'
            for m = 1:nM8, scM8{m} = 1 ./ max(std(fst8{m}(:,gdfr),0,2), eps); end
            normLbl = 'own fast-noise SD';
        case 'bap'
            assert(numel(sIdx) >= 5, ['g8NormMode = ''bap'' needs spikes to measure the ' ...
                'bAP; only %d usable in this segment. Use ''own''.'], numel(sIdx));
            for m = 1:nM8, scM8{m} = 1 ./ divM(:,m); end
            normLbl = 'bAP peak (1.0 = one spike)';
        case 'dff'
            assert(numel(sIdx) >= 5 && ~isempty(bapS), ...
                'g8NormMode = ''dff'' needs spikes and RobustdFFfit.');
            for m = 1:nM8, scM8{m} = bapS ./ divM(:,m); end
            normLbl = 'dF/F';
        case 'bap_file'
            warning(['g8NormMode = ''bap_file'' divides by ScaleFactor.*dFFrobustbAP, ' ...
                'which was fitted on Result.normTraces. These traces use the ' ...
                'UNNORMALISED footprints, so the result is off by ~%.3g and is not ' ...
                'interpretable.'], median(bapMeas(:,1)) / max(median(sfS.*bapS), eps));
            for m = 1:nM8, scM8{m} = 1 ./ max(sfS .* bapS, eps); end
            normLbl = 'ScaleFactor x dFFrobustbAP (NOT valid here)';
        otherwise
            error('g8NormMode must be ''bap'', ''dff'', ''bap_file'' or ''own''.');
    end

    if numel(sIdx) >= 5
        fprintf('\n[R3] bAP peak per trace, measured on %d spikes (median over ROIs):\n', numel(sIdx));
        for m = 1:nM8
            fprintf('     %-16s raw trace p99 %9.4g   measured bAP %9.4g\n', mset{m,2}, ...
                prctile(abs(mset{m,1}(:)), 99), median(bapMeas(:,m), 'omitnan'));
        end
        if ~isempty(sfS)
            fprintf('     RobustdFFfit (on normTraces): ScaleFactor x dFFrobustbAP = %.4g\n', ...
                median(sfS .* bapS));
            fprintf(['     -> mov_res sits %.3g x that scale (unnormalised footprints, ' ...
                     'weight sum ~%.3g),\n        so ScaleFactor is used only as this ' ...
                     'check; the display divides by the MEASURED bAP.\n'], ...
                median(bapMeas(:,1)) / max(median(sfS.*bapS), eps), ...
                median(sum(tovec(double(ftprnt)), 1)));
        end
        if nFloor > 0
            fprintf('     %d ROI(s) had a bAP below 5%% of the median and were floored\n', nFloor);
        end
    end

    remB  = 100*(1 - var(mset{2,1}(:,gdfr),0,2) ./ max(var(mset{1,1}(:,gdfr),0,2), eps));
    remBf = 100*(1 - var(fst8{2}(:,gdfr),0,2)   ./ max(var(fst8{1}(:,gdfr),0,2), eps));
    remBs = 100*(1 - var(sub8{2}(:,gdfr),0,2)   ./ max(var(sub8{1}(:,gdfr),0,2), eps));
    fprintf(['\n[R3] the paper''s dendritic ROIs (%d of them), variance the FLOW ' ...
             'correction removes [median (min-max) %%]:\n'], size(Tbig,1));
    fprintf('     full band     %.2f (%.2f - %.2f)\n', median(remB),  min(remB),  max(remB));
    fprintf('     fast <%d ms    %.2f (%.2f - %.2f)\n', Wsm8, median(remBf), min(remBf), max(remBf));
    fprintf('     %d ms subthr   %.2f (%.2f - %.2f)   <- the trace the paper reports\n', ...
        Wsm8, median(remBs), min(remBs), max(remBs));
    fprintf(['     (large footprints average over many pixels, so a vessel inside one ' ...
             'contributes only a\n      fraction of its area -- expect much smaller ' ...
             'numbers than the sub-ROIs above)\n']);

    if bl8.have
        % Scale-free agreement. If the two corrections remove the same thing, what they
        % LEAVE should agree better than either agrees with the uncorrected trace.
        cc8 = @(A,B) mean(arrayfun(@(r) corr(A(r,gdfr)', B(r,gdfr)'), 1:size(A,1)), 'omitnan');
        fprintf('\n[R3] flow correction vs bleach-dr, per-ROI correlation (mean over ROIs):\n');
        fprintf('     fast <%d ms : raw-flow %.3f | raw-bleach %.3f | flow-bleach %.3f\n', ...
            Wsm8, cc8(fst8{1},fst8{2}), cc8(fst8{1},fst8{3}), cc8(fst8{2},fst8{3}));
        fprintf('     %d ms subthr: raw-flow %.3f | raw-bleach %.3f | flow-bleach %.3f\n', ...
            Wsm8, cc8(sub8{1},sub8{2}), cc8(sub8{1},sub8{3}), cc8(sub8{2},sub8{3}));
        rfb = 100*(1 - var(fst8{3}(:,gdfr),0,2) ./ max(var(fst8{1}(:,gdfr),0,2), eps));
        fprintf(['     fast-band variance removed: flow %.1f%% vs bleach-dr %.1f%% ' ...
                 '(median over ROIs)\n'], median(remBf), median(rfb));
        fprintf(['     NOTE bleach-dr ran on a different preprocessing, so its variance ' ...
                 'ratio also carries\n          that difference -- read it as a rough ' ...
                 'comparison, not a like-for-like one.\n']);
    end

    % ---- figure 1: raw traces, then a kymograph per method AND per residual ---------
    % Everything is on the same time axis in seconds and linkaxes ties them together, so
    % zooming or panning any panel moves all of them -- pick a feature in a kymograph
    % and the traces above follow it.
    %
    % Every panel shares ONE colour scale, taken from the raw kymograph. Scaling each to
    % its own range would renormalise the corrected panels and, worse, blow the residuals
    % up to full contrast however small they are -- the whole point is to see the residual
    % AGAINST the raw.
    %
    % THE RESIDUALS ARE FORMED AFTER NORMALISATION, and that ordering matters:
    %   raw - flow      both share one divisor (same movie), so this is exactly the fitted
    %                   vascular component, in units of one bAP.
    %   raw - bleach    the two have different divisors because they come from different
    %                   pipelines. After normalisation both peak at 1.0 on a spike, so the
    %                   bAP cancels and the difference is meaningful -- but it carries
    %                   EVERY difference between the pipelines (photobleach normalisation,
    %                   motion regression), not only the blood bleach-dr removed. Read it
    %                   as "where the two traces disagree", not as "what bleach-dr took".
    get_figure('R3-traceCompare-bleach');
    [~, ordRem] = sort(remBf, 'descend');
    exBig = ordRem(1:min(g8ExBigN, numel(ordRem)));      % ROIs the flow fit touches most
    tAll8 = (1:nTseg) / fs;
    Xk8   = cellfun(@(B,s) B .* s, mset(:,1)', scM8, 'uni', 0);   % all in bAP units

    kImg = Xk8;  kNm = mset(:,2)';  kCol = mset(:,3)';
    kImg{end+1} = Xk8{1} - Xk8{2};
    kNm{end+1}  = 'raw - flow';         kCol{end+1} = [0.93 0.69 0.13];
    if bl8.have
        kImg{end+1} = Xk8{1} - Xk8{3};
        kNm{end+1}  = 'raw - bleach-dr'; kCol{end+1} = [0.49 0.18 0.56];
    end
    nP8  = numel(kImg);
    tlB8 = tiledlayout(1 + nP8, 1, 'TileSpacing','compact','Padding','compact');

    axT8 = nexttile; hold(axT8,'on');
    offB = 0;  hB = gobjects(1,nP8);
    for e = 1:numel(exBig)
        r  = exBig(e);
        yy = cellfun(@(X) X(r,:), kImg, 'uni', 0);
        stepB = 1.3 * max(cellfun(@(y) prctile(y,99.5) - prctile(y,0.5), yy));
        if ~isfinite(stepB) || stepB <= 0, stepB = 1; end
        for m = 1:nP8
            lw = 0.7;  if m > nM8, lw = 0.9; end       % residuals a touch heavier
            h = plot(axT8, tAll8, yy{m} - offB, '-', 'Color', kCol{m}, 'LineWidth', lw);
            if e == 1, hB(m) = h; end
        end
        text(axT8, tAll8(1), -offB + 0.45*stepB, ...
            sprintf('ROI %d (%.0f \\mum), flow removes %.1f%% of the fast band', ...
            r, dax8(r), remBf(r)), 'FontSize', 8);
        offB = offB + stepB;
    end
    box(axT8,'off'); legend(axT8, hB, kNm, 'Location','northeast');
    ylabel(axT8, sprintf('trace  [%s]', normLbl));
    title(axT8, sprintf('traces and residuals, %d ROI(s) the flow fit touches most', ...
        numel(exBig)));

    tickB8 = unique(round(linspace(1, nRb8, min(8, nRb8))));
    clB8   = prctile(abs(Xk8{1}(:)), 99.5);
    if ~isfinite(clB8) || clB8 <= 0, clB8 = 1; end
    axK8 = gobjects(1, nP8);
    for m = 1:nP8
        axK8(m) = nexttile;
        imagesc(axK8(m), tAll8, 1:nRb8, kImg{m}, [-clB8 clB8]);
        colormap(axK8(m), cmapCorr);
        set(axK8(m), 'YTick', tickB8, 'YTickLabel', round(dax8(tickB8)));
        ylabel(axK8(m), {kNm{m}, 'position (\mum)'});
        if m == nP8, xlabel(axK8(m), 'time (s)'); end
    end
    cbB8 = colorbar(axK8(end));  cbB8.Layout.Tile = 'east';   % one scale for every panel
    cbB8.Label.String = normLbl;
    linkaxes([axT8 axK8], 'x');
    xlim(axT8, [tAll8(1) tAll8(end)]);
    title(tlB8, sprintf(['f=%d  dendritic ROI traces and residuals: mov\\_res, ' ...
        'flow-corrected, bleach-dr (x-axes linked, shared colour scale)'], f));

    % how big is each residual next to the trace it came out of?
    fprintf('\n[R3] residual size in %s (RMS over good frames, median over ROIs):\n', normLbl);
    fprintf('     %-18s %.4g\n', 'raw', median(rms(Xk8{1}(:,gdfr), 2)));
    for m = nM8+1:nP8
        fprintf('     %-18s %.4g   (%.1f%% of raw)\n', kNm{m}, ...
            median(rms(kImg{m}(:,gdfr), 2)), ...
            100*median(rms(kImg{m}(:,gdfr), 2)) / max(median(rms(Xk8{1}(:,gdfr), 2)), eps));
    end

    % ---- figure 2: where in frequency each method acts, and whether they agree ------
    get_figure('R3-traceCompare-bleachSpectra');
    tlP8 = tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
    axpb = nexttile; hold(axpb,'on');
    for k = unique(regSegF(:))'
        patch(axpb, segBandF(k,[1 2 2 1]), [1e-6 1e-6 1e12 1e12], [1 0.9 0.6], ...
            'EdgeColor','none', 'FaceAlpha', 0.45, 'HandleVisibility','off');
    end
    if bl8.have && ~isempty(g8BleachBP)
        patch(axpb, g8BleachBP([1 2 2 1]), [1e-6 1e-6 1e12 1e12], [0.7 0.85 1], ...
            'EdgeColor','none', 'FaceAlpha', 0.5, 'HandleVisibility','off');
    end
    for m = 1:nM8
        Pm8 = zeros(size(mset{m,1},1), numel(fvec8));
        for b = 1:numel(st8)
            ib = st8(b) + (0:nW8-1);
            Xw = mset{m,1}(:,ib) .* scM8{m};  Xw = Xw - mean(Xw,2);  Xw(:,~gdfr(ib)) = 0;
            Zw = abs(fft(Xw .* wHan, [], 2)).^2;
            Pm8 = Pm8 + Zw(:, 1:numel(fvec8));
        end
        plot(axpb, fvec8, median(Pm8/numel(st8), 1), '-', 'Color', mset{m,3}, ...
            'LineWidth', 1.5, 'DisplayName', mset{m,2});
    end
    set(axpb, 'XScale','log', 'YScale','log'); xlim(axpb, [1 fs/2]); box(axpb,'off');
    xlabel(axpb, 'frequency (Hz)'); ylabel(axpb, 'power'); legend(axpb, 'Location','southwest');
    title(axpb, 'PSD (orange = flow regressor bands, blue = bleach-dr bp)');

    axsb = nexttile;
    if bl8.have
        hold(axsb,'on');
        scatter(axsb, remBf, rfb, 22, dax8, 'filled');
        colormap(axsb, turbo); cbB = colorbar(axsb); cbB.Label.String = 'dendritic position (\mum)';
        lm = [min([remBf; rfb]) max([remBf; rfb])];
        plot(axsb, lm, lm, 'k--'); axis(axsb,'square'); box(axsb,'off');
        xlabel(axsb, '% fast variance removed, flow');
        ylabel(axsb, '% fast variance removed, bleach-dr');
        title(axsb, sprintf('per-ROI agreement (r = %.2f)', corr(remBf, rfb, 'rows','complete')));
    else
        text(axsb, 0.5, 0.5, 'BVcorr\_traces.mat not loaded', 'HorizontalAlignment','center');
        axis(axsb, 'off');
    end
    title(tlP8, sprintf(['f=%d  the two corrections act in different bands -- if they ' ...
        'agree anyway, that agreement means something'], f));

    clear Fbg fst8 sub8 mset Xk8 kImg Pm8 Xw Zw
end

clear Ifull Pr8 Pc8 Z X


%% (7i) DEMO: the same pipeline as staged functions (Cohen_lab_BHL_Code/flow)
% Run this AFTER (7d) -- it reuses the tracklets already in the workspace rather than
% re-tracking. Nothing here modifies any variable the cells above depend on; the results
% land in `Rflow` and `cfgFlow`, so the script can be re-run over it freely.
%
% WHAT IS PORTED SO FAR. The refactor is partial and this demo is honest about it:
%
%   flow_defaults   ready   all 78 knobs, grouped by the stage that owns them
%   flow_stamp      ready   config hash -> what actually needs recomputing
%   flow_segments   ready   the ~15 lines of (7e) that anything downstream used
%   flow_pipeline   ready   the driver, but it calls stages that do not exist yet
%   flow_bandpass   TODO    (7c)      flow_regressors  TODO   (7f)
%   flow_track      TODO    (7d)      flow_apply       TODO   (7g)
%                                     flow_compare     TODO   (7h)
%
% So this demo does three things that ARE runnable today:
%   (a) builds cfg and shows which stage owns which parameter;
%   (b) runs flow_segments on the real tracklets and checks it reproduces the script's
%       own densBW / segPx / useSeg -- an equivalence test, not a demonstration;
%   (c) exercises the cache contract, i.e. which stages a given parameter change
%       invalidates. That is the part the script version could not do at all.

flowDir = fullfile(fileparts(which('divide_ROIs')), 'flow');
if ~isfolder(flowDir), flowDir = '/Users/bhlee1117/Documents/GitHub/Cohen_lab_BHL_Code/flow'; end
assert(isfolder(flowDir), 'flow/ not found; expected it next to the other helpers.');
addpath(flowDir);

% ---- (a) the config ------------------------------------------------------------------
cfgFlow = flow_defaults();
fprintf('\n=== flow_defaults: %d knobs across %d stages ===\n', ...
    sum(structfun(@(s) numel(fieldnames(s)), cfgFlow)), numel(fieldnames(cfgFlow)));
fs_ = fieldnames(cfgFlow);
for q = 1:numel(fs_)
    fprintf('   cfg.%-11s %2d  %s\n', fs_{q}, numel(fieldnames(cfgFlow.(fs_{q}))), ...
        strjoin(fieldnames(cfgFlow.(fs_{q}))', ' '));
end

% Match the values this script is actually running with, so (b) compares like for like.
cfgFlow.common.fs        = fs;
cfgFlow.common.smoothMs  = Wsm;
cfgFlow.bandpass.loUm    = sbLoUm;      cfgFlow.bandpass.hiUm   = sbHiUm;
cfgFlow.bandpass.method  = sbMethod;    cfgFlow.bandpass.voltRegress = voltRegress;
cfgFlow.segments.densSm  = recDensSm;   cfgFlow.segments.densPrc = recDensPrc;
cfgFlow.segments.useKeptOnly = recUseKeptOnly;
cfgFlow.segments.minPx   = segMinPx;    cfgFlow.segments.minTrk = segMinTrk;
cfgFlow.regressors.cpca  = oscCPCA;     cfgFlow.regressors.fMin = oscFmin;
cfgFlow.regressors.scoreMin = oscScoreMin;

% ---- (b) flow_segments on the real tracklets, against the script's own answer ---------
assert(exist('trkAll','var') && ~isempty(trkAll), 'Run (7d) first: trkAll is needed.');
trkS = struct('all', {trkAll}, 'pol', trkPol(:), 'cen', trkCen, 'spd', trkSpd, ...
              'keep', keepTrk(:));
segS = flow_segments(trkS, [nrow ncol], cfgFlow, ...
                     struct('bvBW', bvBW, 'pxUm', pxUm));

fprintf('\n=== flow_segments vs the script''s (7e)/(7f) ===\n');
if exist('densBW','var')
    fprintf('   densBW identical            : %d  (%d px vs %d px)\n', ...
        isequal(segS.densBW, densBW), nnz(segS.densBW), nnz(densBW));
end
if exist('segPx','var') && exist('useSeg','var')
    okN  = numel(segS.px) == numel(useSeg);
    okPx = okN && all(arrayfun(@(q) isequal(sort(segS.px{q}), sort(segPx{useSeg(q)})), ...
                               1:numel(useSeg)));
    fprintf('   same number of used segments: %d  (%d vs %d)\n', okN, numel(segS.px), numel(useSeg));
    fprintf('   same pixels in every segment: %d\n', okPx);
    fprintf('   tracklets per segment        : %s   vs   %s\n', ...
        mat2str(segS.nTrkSeg), mat2str(arrayfun(@(q) nnz(trkSeg==q), useSeg)));
    if ~okPx
        fprintf(['   -> a mismatch here is a REAL difference, not a rounding artefact. ' ...
                 'Check that\n      the cfg overrides above match the values (7e)/(7f) ' ...
                 'actually ran with.\n']);
    end
end

% ---- (c) what a parameter change invalidates -----------------------------------------
% This is the thing the six-cell version could not answer. Each stamp folds in the
% stage's own knobs AND its upstream stamp, so the dependency is checked rather than
% remembered.
stg   = {'bandpass','track','segments','regressors','apply','compare'};
probe = {'(nothing)',          '',            '',         0
         'bandpass.loUm',      'bandpass',    'loUm',     7
         'segments.densPrc',   'segments',    'densPrc',  75
         'regressors.fMin',    'regressors',  'fMin',     60
         'apply.splitWin',     'apply',       'splitWin', 800};
stampTab = cell(size(probe,1), numel(stg));
for q = 1:size(probe,1)
    cq = cfgFlow;
    if ~isempty(probe{q,2}), cq.(probe{q,2}).(probe{q,3}) = probe{q,4}; end
    prev = flow_stamp(cq.bandpass, cq.common);      % chained exactly as flow_pipeline does
    stampTab{q,1} = prev;
    for j = 2:numel(stg)
        prev = flow_stamp(cq.(stg{j}), prev);
        stampTab{q,j} = prev;
    end
end
fprintf('\n=== changing one knob re-runs these stages ===\n');
hdr = sprintf('%-11s', stg{:});
fprintf('   %-20s %s\n', 'knob changed', hdr);
for q = 2:size(probe,1)
    row = '';
    for j = 1:numel(stg)
        if strcmp(stampTab{1,j}, stampTab{q,j}), row = [row sprintf('%-11s','.')];
        else,                                    row = [row sprintf('%-11s','RERUN')]; end %#ok<AGROW>
    end
    fprintf('   %-20s %s\n', probe{q,1}, row);
end
fprintf('   (RERUN = invalidated; . = reused from cache)\n');

% ---- (d) the GUI ---------------------------------------------------------------------
% This is the intended way to drive the pipeline: one window, the six steps across the
% top, the selected step's parameters in an editable table, its result in the preview.
% It blocks until Done and hands back the parameters and the results, so it drops into a
% script without changing how anything downstream is written.
%
% The step buttons are colour-coded -- green up to date, AMBER stale (a parameter it
% depends on changed since it ran), grey not run, red not ported yet. Amber is the part
% the cell-based version could not do: it tells you what a change actually invalidated
% instead of leaving you to re-run everything or re-run one cell and hope.
dataFlow = struct( ...
    'mov_res',  mov_res, ...
    'fs',       fs, ...
    'pxUm',     pxUm, ...
    'badfr',    badfr, ...
    'spike',    full(Result.spike(1, segFrames)) > 0, ...
    'ftprnt',   ftprnt, ...
    'bvBW',     bvBW, ...
    'neuronMap',neuronMap, ...
    'voltTr',   PCresult{f}.normTraces(:, segFrames), ...
    'tag',      fpath{f});

fprintf(['\n=== the GUI ===\n' ...
         '   [cfgFlow, Rflow] = flow_gui(dataFlow, cfgFlow);\n' ...
         '   [cfgFlow, Rflow] = flow_gui(dataFlow, cfgFlow, Rflow);   %% resume\n' ...
         '   (bandpass/track/regressors/apply/compare show red until they are ported;\n' ...
         '    segments runs once its upstream two exist.)\n' ...
         '   Headless equivalent, same cache rules:\n' ...
         '   Rflow = flow_pipeline(dataFlow, cfgFlow, Rflow, ''segments'');\n']);

Rflow = struct('segments', segS, 'tag', 'partial: flow_segments only');
% [cfgFlow, Rflow] = flow_gui(dataFlow, cfgFlow, Rflow);   %% <- uncomment to open it
clear fs_ okN okPx probe stampTab stg prev cq hdr row
set_fig         
%%
Volttr_tmp=tovec(mov_res)'*tovec(Result.ftprnt);
STAmovie=importdata(fullfile(fpath{f},'SS_STAmovie.mat'));
STAmap = zeros(size(mov_res(:,:,1)));
STAmap(6:end-6,6:end-6) = max(-STAmovie,[],3);
dMask=point2img(Result.SWC, 5, size(Result.ref_im)) > 0;
nUse = 1:500;
dataFlow = struct( ...
    'mov_res', mov_res(:,:,nUse), ... 
    'fs',      1000, ...
    'pxUm',    Pixelsize(f), ...
    'bvBW',    bvBW, ...
    'neuronMap', STAmap, ...
    'anaMask',dMask,...
    'ftprnt',Result.ftprnt,...
    'tag',     'demo-1000');

[cfgFlow, Rflow] = flow_gui(dataFlow,cfgFlow,Rflow);
%[cfgFlow, Rflow] = flow_gui(dataFlow);
%% (9) Blood estimate from clusters of the fast cross-correlation matrix    [item 3]
%  Cxt(i,j,tau) = corr(odd_i(t), even_j(t+tau)) of the FAST (<20 ms) component --
%  built for EVERY pair including the diagonal from the odd/even checkerboard pixel
%  halves, so it is shot-noise-free throughout (see cell (4c)).  Blood fluctuation is
%  spatially coherent over ~10 um, so it appears as compact, strongly correlated
%  BLOCKS of neighbouring sub-ROIs.  Here we
%     (a) collapse Cxt to a LAG-TOLERANT matrix Rmax (see below) and cluster it
%         (average linkage on 1-Rmax), showing the dendrogram, the cut threshold
%         and the correlation matrix reordered by cluster,
%     (b) keep the BIG clusters (>= minClustSize sub-ROIs) and average their
%         z-scored fast traces -> one blood regressor per cluster,
%     (c) map every pixel's covariance with those regressors -- no fit, no
%         interpolation, bad frames simply excluded -- and split the bandpassed
%         component into an amplitude and a phase map via its bp/d1 quadrature pair.
%     Regressing them out of the movie comes later, once (a)-(c) are settled.
% ============================ PARAMETERS (all of them) ============================
% -- correlation & clustering --------------------------------------------------
lagBlood       = 30;   % ms, max |lag| searched for the lag-tolerant correlation Rmax
nClustFast     = 25;   % clusters requested from the fast correlation matrix
% -- cluster selection: a blood cluster is BIG, SPATIALLY DISTRIBUTED, above background
minClustSize   = 15;   % minimum number of sub-ROIs
minSpreadUm    = 50;   % minimum median pairwise centroid distance (um)
minWithinRatio = 1.1;    % within-corr floor, as a multiple of the background mean
nRegMax        = 20;   % at most this many clusters become blood regressors
% -- trace preprocessing --------------------------------------------------------
detrendWin     = 200;      % ms, movmedian window used to detrend each member trace
                           % (the ONLY temporal processing applied)
psdWinP        = 2048;     % samples per Welch / Lomb-Scargle window in (b)
cwtMinCycles   = 5;        % the 1/f-corrected wavelet is only trusted above
                           % cwtMinCycles/(median gap-free segment length)
% -- display / output ----------------------------------------------------------
kSweep          = 2:60;  % cluster counts scanned for the coherence-vs-k diagnostic panel
nMapShow        = 12;    % total panels in R3-bloodClusterMaps: every KEPT cluster plus
                         % the best-scoring rejected ones until this many are drawn
clusterMapCols  = 2;     % columns in that figure (the FOV is wide, so keep it small)
% =================================================================================

% --- (a1) lag-tolerant correlation: blood flow has a phase delay between locations ---
% Two sub-ROIs on the same vessel see the same fluctuation shifted in time, so the
% tau=0 slice underestimates their coherence.  Take the MAXIMUM of the correlogram
% over |tau| <= lagBlood instead.  The search window is symmetric and the tensor obeys
% Cxt(j,i,tau) = Cxt(i,j,-tau), so Rmax comes out symmetric by construction; the lag at
% which the maximum sits is antisymmetric and is the propagation delay between the pair.
% Caveat: max-over-lag is an upward-biased estimate of correlation (it picks the largest
% of 2*lagBlood+1 noisy values), but the bias is the same for every pair, so it shifts
% the whole matrix rather than distorting the cluster structure.
lagSel   = abs(lags) <= lagBlood;
lagVals  = lags(lagSel);
[Rmax, iMax] = max(Cxt(:,:,lagSel), [], 3);            % max ignores NaN unless all-NaN
lagAtMax = reshape(lagVals(iMax), Msub, Msub);         % ms, delay of the maximum
allNaN   = all(isnan(Cxt(:,:,lagSel)), 3);
Rmax(allNaN) = NaN;  lagAtMax(allNaN) = NaN;
Rmax = (Rmax + Rmax')/2;                               % kill residual asymmetry from ties

offD = triu(true(Msub),1);
fprintf('\n[R3] lag tolerance (|tau| <= %d ms) on the fast odd/even correlation:\n', lagBlood);
fprintf('   mean off-diagonal  tau=0 : %+.4f   max-over-lag : %+.4f\n', ...
    mean(R0(offD),'omitnan'), mean(Rmax(offD),'omitnan'));
lagOff = lagAtMax(offD);  lagOff = lagOff(~isnan(lagOff));
fprintf('   pairs whose maximum is NOT at tau=0 : %.0f%%   median |delay| = %.1f ms\n', ...
    100*mean(lagOff ~= 0), median(abs(lagOff)));

% --- (a2) hierarchical clustering of the lag-tolerant matrix ---
% Clustering directly on D = 1 - Rmax with average linkage CHAINS: every off-diagonal
% correlation here is small (background ~0.04, blocks ~0.1), so all distances sit in
% a narrow band around 1, average linkage cannot separate them, and 'maxclust' returns
% one giant cluster plus a few peeled-off singletons. Two changes fix it:
%   (i)  cluster on the correlation PROFILE, not on the raw correlation. Row i of Rmax
%        is sub-ROI i's correlation with everything else; two sub-ROIs on the same
%        vessel have SIMILAR PROFILES even when their mutual r is only 0.1. Rows are
%        z-scored, so Euclidean distance between them is monotone in 1 - corr(profiles)
%        and spans the full range instead of a narrow band near 1.
%   (ii) 'ward' linkage (legitimate here: the input really is Euclidean), which
%        minimises within-cluster scatter and therefore refuses to chain.
% The DIAGONAL must not be forced to 1: off-diagonal values here are ~0.02-0.05, so a
% single 1.0 per row inflates that row's std by ~2x and adds a large constant to every
% pairwise distance. Setting it to the row's own off-diagonal mean makes it neutral --
% it contributes nothing to the row's deviation profile. (Measured on synthetic data
% matched to this dataset: cluster purity 0.67 with diag=1 vs 0.75 with diag=row mean.)
Rc = Rmax;  Rc(isnan(Rc)) = 0;  Rc = (Rc + Rc')/2;
dIdx = 1:Msub+1:Msub^2;
Rc(dIdx) = 0;
rowMuOff = sum(Rc,2) ./ (Msub-1);                      % mean of the off-diagonal entries
Rc(dIdx) = rowMuOff;                                   % neutral diagonal
Rprof = (Rc - mean(Rc,2)) ./ max(std(Rc,0,2), eps);    % z-scored correlation profile
Dvec  = pdist(Rprof, 'euclidean');
Zf    = linkage(Dvec, 'ward');
cluFast = cluster(Zf, 'maxclust', nClustFast);
% height that produces exactly nClustFast clusters (midpoint of the two merge heights
% that straddle the cut), drawn on the dendrogram below
hSorted  = sort(Zf(:,3), 'descend');
cutH     = mean(hSorted(nClustFast-1:nClustFast));
cluSize  = accumarray(cluFast, 1);

% --- keep clusters that are LARGE, SPATIALLY DISTRIBUTED and above background -----
% A vessel runs across the field, so its cluster is spread over tens of um and has
% many members, but only a MODERATE pairwise correlation. The small, very-high-r
% blocks hugging the diagonal are the opposite: a handful of ADJACENT sub-ROIs
% sharing PSF / dendritic signal over a few um. A coherence-only criterion therefore
% keeps exactly the wrong clusters. Require a minimum spatial extent and size, then
% rank by the TOTAL excess shared correlation
%       score = (#pairs) x (within-corr - background)
% which puts a large, moderately correlated vessel cluster ahead of a tiny, tightly
% correlated neighbour pair.
Deuc = squareform(pdist(coordS)) * Pixelsize(f);       % um between sub-ROI centroids
bgMu = mean(Rmax(offD), 'omitnan');
nClu = max(cluFast);
cluWithin = nan(nClu,1);  cluSpread = nan(nClu,1);  cluGeo = nan(nClu,1);  cluScore = zeros(nClu,1);
for c = 1:nClu
    mk = cluFast == c;  nc = nnz(mk);
    if nc < 2, continue; end
    mIn = triu(true(nc),1);
    w  = Rmax(mk,mk);   cluWithin(c) = mean(w(mIn), 'omitnan');
    de = Deuc(mk,mk);   cluSpread(c) = median(de(mIn));
    dg = Dgeo(mk,mk);   cluGeo(c)    = median(dg(mIn), 'omitnan');
    cluScore(c)  = nnz(mIn) * max(cluWithin(c) - bgMu, 0);
end
keepClu = find(cluSize >= minClustSize & cluSpread >= minSpreadUm & ...
               cluWithin >= bgMu*minWithinRatio);
[~, si] = sort(cluScore(keepClu), 'descend');          % largest shared variance first
bigClu  = keepClu(si);
if numel(bigClu) > nRegMax, bigClu = bigClu(1:nRegMax); end
nBig    = numel(bigClu);
assert(nBig > 0, ['No cluster passed size>=%d, spread>=%g um and within-corr>=%.4f. ' ...
    'Relax minClustSize/minSpreadUm/minWithinRatio, or change nClustFast.'], ...
    minClustSize, minSpreadUm, bgMu*minWithinRatio);
fprintf('\n[R3] clustering (ward on correlation profiles): %d clusters, sizes %d-%d\n', ...
    nClu, min(cluSize), max(cluSize));
fprintf('   background off-diagonal corr = %+.4f  ->  coherence floor %+.4f (%.2fx)\n', ...
    bgMu, bgMu*minWithinRatio, minWithinRatio);
fprintf('   kept %d clusters (size>=%d, spread>=%g um, ranked by #pairs x excess corr)\n', ...
    nBig, minClustSize, minSpreadUm);

% --- (a3) dendrogram + cut threshold + correlation matrices reordered by cluster ---
% ordD is the dendrogram leaf order, so members of a cluster are contiguous; the same
% order is applied to the tau=0 matrix and to the delay map for direct comparison.
get_figure('R3-bloodDendrogram');
tlD = tiledlayout(2, 5, 'TileSpacing','compact','Padding','compact');

axDen = nexttile([1 5]);
[~, ~, ordD] = dendrogram(Zf, 0, 'ColorThreshold', cutH, 'Parent', axDen);
set(axDen, 'XTick', []); box(axDen,'off'); hold(axDen,'on');
yline(axDen, cutH, 'r--', sprintf('cut = %.3g  ->  %d clusters', cutH, nClustFast), ...
    'LineWidth', 1.2, 'LabelHorizontalAlignment','left');
ylabel(axDen, 'ward distance (correlation profile)');
title(axDen, sprintf('f=%d  ward linkage on the fast odd/even correlation profiles (%d sub-ROIs)', f, Msub));

% cluster block boundaries in leaf order
cluLeaf  = cluFast(ordD);
bnd      = find(diff(cluLeaf(:)) ~= 0) + 0.5;

% shared color limits for the two correlation panels (max-over-lag is biased upward,
% so a fixed ccLim would saturate it and hide the comparison with tau=0)
climC = [min(0, prctile(Rmax(offD), 1)) prctile(Rmax(offD), 99.5)];
matSets = {Rmax(ordD,ordD),     sprintf('max over |\\tau| \\leq %d ms', lagBlood), climC, ccFast; ...
           R0(ordD,ordD),       '\tau = 0',                                       climC, ccFast; ...
           lagAtMax(ordD,ordD), 'lag at the maximum (ms)',    [-lagBlood lagBlood], cmapCorr};
for s = 1:3
    axm = nexttile;
    imagesc(axm, matSets{s,1}, matSets{s,3}); axis(axm,'image'); colormap(axm, matSets{s,4});
    colorbar(axm); hold(axm,'on');
    for b = bnd(:)'                                     % cluster boundaries
        plot(axm, [0.5 Msub+0.5], [b b], 'k-', 'LineWidth', 0.6);
        plot(axm, [b b], [0.5 Msub+0.5], 'k-', 'LineWidth', 0.6);
    end
    set(axm, 'XTick', [], 'YTick', []);
    title(axm, matSets{s,2});
    if s == 1
        ylabel(axm, 'sub-ROI (dendrogram leaf order)');
        for q = 1:nBig                                  % outline the big clusters
            r = find(cluLeaf == bigClu(q));
            plot(axm, [r(1)-0.5 r(end)+0.5 r(end)+0.5 r(1)-0.5 r(1)-0.5], ...
                      [r(1)-0.5 r(1)-0.5 r(end)+0.5 r(end)+0.5 r(1)-0.5], ...
                 'g-', 'LineWidth', 1.5);
        end
    end
end

% (4th panel) how coherent are the clusters as a function of how many you ask for?
% For each k: the size-weighted mean within-cluster correlation over clusters that pass
% minClustSize.  Too few clusters -> one blob at background level; too many -> tiny
% clusters.  Pick nClustFast where the curve is comfortably above the coherence floor.
swWithin = nan(size(kSweep));  swKept = nan(size(kSweep));
for kk = 1:numel(kSweep)
    ct = cluster(Zf, 'maxclust', kSweep(kk));
    ws = []; ns = [];
    for c = 1:max(ct)
        mk = ct == c;
        if nnz(mk) < minClustSize, continue; end
        w = Rmax(mk,mk);
        ws(end+1) = mean(w(triu(true(nnz(mk)),1)), 'omitnan');   %#ok<SAGROW>
        ns(end+1) = nnz(mk);                                     %#ok<SAGROW>
    end
    if ~isempty(ws), swWithin(kk) = sum(ws.*ns)/sum(ns);  swKept(kk) = numel(ws); end
end
axSw = nexttile; yyaxis(axSw,'left'); hold(axSw,'on');
plot(axSw, kSweep, swWithin, '-', 'LineWidth', 1.4);
yline(axSw, bgMu*minWithinRatio, 'k--');
yline(axSw, bgMu, ':', 'background');
xline(axSw, nClustFast, 'r-', 'LineWidth', 1.2);
ylabel(axSw, 'within-cluster corr (size-weighted)');
yyaxis(axSw,'right'); plot(axSw, kSweep, swKept, '-');
ylabel(axSw, sprintf('# clusters with n \\geq %d', minClustSize));
xlabel(axSw, 'number of clusters'); grid(axSw,'on'); box(axSw,'off');
title(axSw, 'coherence vs cluster count');

% (5th panel) the selection itself: size vs coherence, colored by spatial spread.
% Blood clusters live to the RIGHT (many members) and are WARM (spread over tens of
% um); the tiny adjacent-ROI blocks sit top-left and cold, and are rejected.
vc = find(cluSize >= 2);
axSel = nexttile; hold(axSel,'on');
scatter(axSel, cluSize(vc), cluWithin(vc), 45, cluSpread(vc), 'filled');
plot(axSel, cluSize(bigClu), cluWithin(bigClu), 'ko', 'MarkerSize', 12, 'LineWidth', 1.3);
for q = 1:nBig
    text(axSel, cluSize(bigClu(q))*1.08, cluWithin(bigClu(q)), sprintf('c%d', bigClu(q)), 'FontSize', 8);
end
xline(axSel, minClustSize, 'k--'); yline(axSel, bgMu*minWithinRatio, 'k--');
yline(axSel, bgMu, ':');
set(axSel, 'XScale','log'); colormap(axSel, turbo); caxis(axSel, [0 max(cluSpread(vc))]);
cbS = colorbar(axSel); cbS.Label.String = 'spatial spread (\mum)';
xlabel(axSel, 'cluster size (sub-ROIs)'); ylabel(axSel, 'within-cluster corr');
grid(axSel,'on'); box(axSel,'off');
title(axSel, sprintf('selection: circled = kept (spread \\geq %g \\mum)', minSpreadUm));

title(tlD, sprintf(['f=%d  fast (<%d ms) odd/even correlation, ordered by cluster ' ...
    '(green = the %d clusters kept as blood regressors)'], f, Wsm, nBig));

% kept clusters, ranked by score; the rejected ones are listed underneath so it is
% obvious WHY they were dropped (too small / too local / not above background)
fprintf('[R3] kept clusters (within-corr max-over-lag | tau=0):\n');
for k = 1:nBig
    c = bigClu(k);  mk = cluFast == c;  mIn = triu(true(nnz(mk)),1);
    wIn = Rmax(mk,mk);  wIn0 = R0(mk,mk);  wLag = abs(lagAtMax(mk,mk));  wOut = Rmax(mk,~mk);
    fprintf(['   c%-3d n=%3d  spread %5.1f um (geo %5.1f)  within %+.3f | %+.3f  ' ...
             'between %+.3f  |delay| %.1f ms  score %.1f\n'], ...
        c, nnz(mk), cluSpread(c), cluGeo(c), mean(wIn(mIn),'omitnan'), mean(wIn0(mIn),'omitnan'), ...
        mean(wOut(:),'omitnan'), median(wLag(mIn),'omitnan'), cluScore(c));
end
rej = setdiff(find(cluSize >= 2), bigClu);
[~, rs] = sort(cluScore(rej), 'descend');  rej = rej(rs);
fprintf('[R3] top rejected clusters (reason):\n');
for k = 1:min(6, numel(rej))
    c = rej(k);  why = {};
    if cluSize(c)   <  minClustSize,          why{end+1} = 'too small';   end %#ok<SAGROW>
    if cluSpread(c) <  minSpreadUm,           why{end+1} = 'too local';   end %#ok<SAGROW>
    if cluWithin(c) <  bgMu*minWithinRatio,   why{end+1} = 'at background'; end %#ok<SAGROW>
    if isempty(why),                          why = {sprintf('below top %d', nRegMax)}; end
    fprintf('   c%-3d n=%3d  spread %5.1f um  within %+.3f  score %5.1f   -> %s\n', ...
        c, cluSize(c), cluSpread(c), cluWithin(c), cluScore(c), strjoin(why, ', '));
end

% --- (a4) each cluster drawn on the dendrite map, with its statistics in the title ---
% One panel per cluster: the members are painted on the morphology so you can see
% whether the cluster really follows a vessel (extended, crossing the dendrite) or is
% just a local patch. Kept clusters have a black title, the closest REJECTED ones a
% red title with the reason -- use these to decide whether minClustSize / minSpreadUm /
% minWithinRatio / nRegMax need moving.
ftLc = double(ftprintsmallest);  denc = sum(ftLc,3);
[~, labImg] = max(ftLc, [], 3);  inMaskC = denc > eps;  labImg(~inMaskC) = 0;

% Panel count = nMapShow. Every KEPT cluster is always drawn (never hide a regressor);
% the remaining slots go to the best-scoring REJECTED clusters, so raising nMapShow
% simply shows more near-miss candidates.
nRejUse  = max(0, min(nMapShow - nBig, numel(rej)));
showClu  = [bigClu(:); rej(1:nRejUse)];
showKept = [true(nBig,1); false(nRejUse,1)];
if nBig <= 7, cluCmap = lines(max(nBig,2)); else, cluCmap = turbo(nBig); end
rejColr  = [1 0 1];        % magenta: a saturated hue stays visible on the gray image
                           % (a gray overlay disappears into the bright dendrite pixels)
fprintf(['[R3] cluster maps: %d panels = %d kept + %d of %d rejected ' ...
         '(%d clusters total, %d with n>=2)\n'], ...
    numel(showClu), nBig, nRejUse, numel(rej), nClu, numel(rej)+nBig);

get_figure('R3-bloodClusterMaps');
tlM = tiledlayout(ceil(numel(showClu)/clusterMapCols), clusterMapCols, ...
                  'TileSpacing','compact','Padding','compact');
for q = 1:numel(showClu)
    c   = showClu(q);
    mk  = cluFast == c;
    inC = inMaskC & ismember(labImg, find(mk));            % pixels belonging to this cluster
    if showKept(q), colq = cluCmap(q,:); else, colq = rejColr; end

    axq = nexttile;
    imagesc(axq, refImg, climRef); colormap(axq, gray); axis(axq,'image','off'); hold(axq,'on');
    ovl = cat(3, colq(1)*ones(nrow,ncol), colq(2)*ones(nrow,ncol), colq(3)*ones(nrow,ncol));
    hOv = imagesc(axq, ovl);  set(hOv, 'AlphaData', 0.85*double(inC));
    % outline every member as well, so small clusters stay visible whatever the background
    Bc = bwboundaries(inC);
    for b = 1:numel(Bc), plot(axq, Bc{b}(:,2), Bc{b}(:,1), '-', 'Color', colq, 'LineWidth', 1); end
    contour(axq, bvBW, [0.5 0.5], 'c', 'LineWidth', 0.8);
    plot(axq, coordS(mk,1), coordS(mk,2), '.', 'Color', colq, 'MarkerSize', 6);

    mIn  = triu(true(nnz(mk)),1);
    wIn  = Rmax(mk,mk);  wOut = Rmax(mk,~mk);  wLag = abs(lagAtMax(mk,mk));
    ttl = sprintf('c%d   n=%d   spread %.0f\\mum (geo %.0f)   score %.0f', ...
        c, nnz(mk), cluSpread(c), cluGeo(c), cluScore(c));
    ttl2 = sprintf('r_{in}=%.3f  r_{out}=%.3f  (bg %.3f, floor %.3f)  |lag| %.1f ms', ...
        mean(wIn(mIn),'omitnan'), mean(wOut(:),'omitnan'), bgMu, bgMu*minWithinRatio, ...
        median(wLag(mIn),'omitnan'));
    if showKept(q)
        title(axq, {ttl, ttl2}, 'Color', 'k');
    else
        why = {};
        if cluSize(c)   <  minClustSize,        why{end+1} = 'too small';      end %#ok<SAGROW>
        if cluSpread(c) <  minSpreadUm,         why{end+1} = 'too local';      end %#ok<SAGROW>
        if cluWithin(c) <  bgMu*minWithinRatio, why{end+1} = 'at background';  end %#ok<SAGROW>
        if isempty(why), why = {sprintf('below top %d', nRegMax)}; end
        title(axq, {['REJECTED (' strjoin(why,', ') ')   ' ttl], ttl2}, 'Color', [0.75 0 0]);
    end
end
title(tlM, sprintf(['f=%d  clusters on the morphology (cyan = bvMask). Kept = black title, ' ...
    'rejected = red. Thresholds: n\\geq%d, spread\\geq%g \\mum, r\\geq%.2fx background, top %d by score'], ...
    f, minClustSize, minSpreadUm, minWithinRatio, nRegMax));

% --- (b) PSD of the cluster sub-ROI traces, measured three independent ways --------
% Estimators disagree for well-understood reasons, so all three are drawn instead of
% trusting one:
%   Welch         pwelch on an interpolated copy   -- needs continuous data
%   Lomb-Scargle  nanPSD on the gapped trace       -- handles gaps natively, more leakage
%   wavelet       time-average of |CWT|^2, DIVIDED BY f
%
% The 1/f correction on the wavelet is not optional. A Morlet CWT is constant-Q: its
% analysis band at frequency f is proportional to f, so the raw scalogram is
% scalogram(f) ~ f * PSD(f) -- tilted by exactly one decade per decade. Measured on
% synthetic noise: pwelch log-log slope +0.04 vs scalogram +1.01 for white noise, and
% -1.95 vs -0.94 for 1/f^2 noise; the offset is +1.00 in both cases. An uncorrected
% scalogram therefore moves the apparent peak of a falling spectrum to a HIGHER
% frequency than the PSD shows, which is exactly the kind of mismatch to watch for.
% The correction is only applied above fMinCWT (see below): dividing by f where the
% wavelet has too few cycles to work with would hand the low end a spurious peak.
%
% The spectra below are the MEAN OF THE MEMBERS' PSDs, not the PSD of their average, so
% the averaging step cannot manufacture or hide a peak. A size-matched control drawn
% from the sub-ROIs outside every kept cluster is shown as a dashed line.
psdWin  = psdWinP;
psdBand = [2 fs/2-1];                 % Hz displayed

% detrended, z-scored member traces; bad frames NaN and never interpolated for a statistic
Fdet = Tall - movmedian(Tall, detrendWin, 2, 'omitnan');
Fdet(:, badfr) = NaN;
Fdet = (Fdet - mean(Fdet, 2, 'omitnan')) ./ std(Fdet, 0, 2, 'omitnan');
Fint = interpolateNaN(Fdet);  Fint(isnan(Fint)) = 0;      % Welch copy only

% lowest frequency the wavelet can be trusted at: it needs several cycles inside a
% gap-free segment, and the 1/f correction amplifies whatever it gets wrong down there.
segA    = regionprops(bwlabel(gdfr), 'Area');
segLenM = median([segA.Area]);
fMinCWT = cwtMinCycles * fs / max(segLenM, 1);
fprintf('\n[R3] median gap-free segment = %.0f ms -> wavelet trusted above %.1f Hz\n', ...
    segLenM, fMinCWT);

outPool = find(~ismember(cluFast, bigClu));               % everything outside the clusters
psdGrp  = cell(1, nBig+1);  psdLab = cell(1, nBig+1);
for k = 1:nBig
    psdGrp{k} = find(cluFast == bigClu(k));
    psdLab{k} = sprintf('c%d (n=%d)', bigClu(k), cluSize(bigClu(k)));
end
nCtl = min(max(cluSize(bigClu)), numel(outPool));
psdGrp{end} = outPool(unique(round(linspace(1, numel(outPool), max(nCtl,1)))));
psdLab{end} = sprintf('outside clusters (n=%d)', numel(psdGrp{end}));

estName = {'Welch (interpolated)', 'Lomb-Scargle (gapped)', 'wavelet / f'};
Pest = cell(3, nBig+1);  Fest = cell(1,3);
for g = 1:nBig+1
    idx = psdGrp{g};
    % (1) Welch on the interpolated copy
    [Pw, fw] = pwelch(Fint(idx,:)', hann(psdWin), psdWin/2, [], fs);
    Pest{1,g} = mean(Pw, 2);  Fest{1} = fw(:);
    % (2) Lomb-Scargle straight on the gapped traces
    [Pl, fl] = nanPSD(Fdet(idx,:), fs, psdWin);
    Pest{2,g} = mean(Pl, 1)';  Fest{2} = fl(:);
    % (3) wavelet time-average, corrected for the constant-Q bandwidth
    acc = [];
    for q = 1:numel(idx)
        [wt, fc] = get_waveletTransform(Fdet(idx(q),:), fs, psdBand);
        sc = mean(abs(wt).^2, 2, 'omitnan') ./ fc(:);     % <-- the 1/f correction
        if isempty(acc), acc = sc; else, acc = acc + sc; end
    end
    sc = acc / max(numel(idx),1);
    sc(fc(:) < fMinCWT) = NaN;      % below this the wavelet has too few cycles per
                                    % gap-free segment and 1/f blows the estimate up
    Pest{3,g} = sc;  Fest{3} = fc(:);
end

get_figure('R3-clusterPSD');
tlP = tiledlayout(1, 3, 'TileSpacing','compact','Padding','compact');
for e = 1:3
    axe = nexttile; hold(axe,'on'); set(axe,'XScale','log','YScale','log');
    for g = 1:nBig+1
        if g <= nBig, st = '-';  lw = 1.3;  co = cluCmap(g,:);
        else,         st = '--'; lw = 1.1;  co = [0 0 0];  end
        plot(axe, Fest{e}, Pest{e,g}, st, 'Color', co, 'LineWidth', lw, 'DisplayName', psdLab{g});
    end
    xlim(axe, psdBand); grid(axe,'on'); box(axe,'off');
    xlabel(axe,'frequency (Hz)'); ylabel(axe,'PSD');
    title(axe, estName{e});
    if e == 1, legend(axe, 'Location','southwest'); end
end
title(tlP, sprintf(['f=%d  mean PSD of the member sub-ROIs (dashed = size-matched ' ...
    'non-cluster control). The wavelet is 1/f-corrected so all three are comparable.'], f));

fprintf('\n[R3] spectral peak of each cluster, by estimator (search %g-%g Hz):\n', psdBand(1), psdBand(2));
for g = 1:nBig+1
    pk = nan(1,3);
    for e = 1:3
        okf = Fest{e} >= psdBand(1) & Fest{e} <= psdBand(2) & isfinite(Pest{e,g});
        if any(okf)
            ff = Fest{e}(okf);  pp = Pest{e,g}(okf);
            [~, ip] = max(pp);  pk(e) = ff(ip);
        end
    end
    fprintf('   %-28s Welch %6.1f Hz | Lomb-Scargle %6.1f Hz | wavelet/f %6.1f Hz\n', ...
        psdLab{g}, pk(1), pk(2), pk(3));
end
fprintf(['   (large disagreement between the three means there is no real peak -- the\n' ...
         '    spectrum is smooth and "max of the PSD" is picking noise.)\n']);

% --- (c) do the members actually share a signal?  overlay + pairwise correlation ----
% If the member traces are only weakly correlated there is nothing to average: the mean
% would be mostly noise, and the honest move is to hand ALL member traces to
% SeeResiduals as separate regressors instead of collapsing them into one.
%
% coherenceGain = std(mean of members) / (mean std of members / sqrt(n))
%   ~ 1     members are independent -> averaging buys nothing beyond sqrt(n) noise
%           averaging, i.e. there is no shared component
%   ~ sqrt(n) members are perfectly coherent -> the average IS the shared signal
ovWin = min(600, numel(tShow));                 % ms shown in the overlay
tOv   = tShow(1:ovWin);  tOvMs = tOv - tOv(1);

get_figure('R3-clusterTraceOverlay');
tlO = tiledlayout(nBig, 2, 'TileSpacing','compact','Padding','compact');
fprintf('\n[R3] within-cluster trace similarity (are these worth averaging?):\n');
for k = 1:nBig
    mk = find(cluFast == bigClu(k));
    A  = Fdet(mk, :);

    % pairwise correlation over the frames where both members are defined
    C  = corr(A(:, gdfr)', 'rows', 'pairwise');
    mIn = triu(true(numel(mk)),1);  cv = C(mIn);  cv = cv(isfinite(cv));

    mAll = mean(A, 1, 'omitnan');
    gk   = gdfr & sum(~isnan(A),1) >= max(2, ceil(0.5*numel(mk)));
    gain = std(mAll(gk), 0, 2, 'omitnan') / ...
           max(mean(std(A(:,gk), 0, 2, 'omitnan'), 'omitnan') / sqrt(numel(mk)), eps);
    fprintf(['   c%-3d n=%3d : pairwise corr median %+.3f (IQR %+.3f..%+.3f) | ' ...
             'coherenceGain %.2f of sqrt(n)=%.2f\n'], ...
        bigClu(k), numel(mk), median(cv), prctile(cv,25), prctile(cv,75), ...
        gain, sqrt(numel(mk)));

    % left: every member trace overlaid, cluster mean on top
    axo = nexttile; hold(axo,'on');
    plot(axo, tOvMs, A(:, tOv)', '-', 'Color', [cluCmap(k,:) 0.25], 'LineWidth', 0.4);
    plot(axo, tOvMs, mAll(tOv), 'k-', 'LineWidth', 1.4);
    xlim(axo, [0 tOvMs(end)]); box(axo,'off'); ylabel(axo, sprintf('c%d  (z)', bigClu(k)));
    if k == nBig, xlabel(axo, 'time (ms)'); end
    title(axo, sprintf('c%d  %d members overlaid (black = mean)', bigClu(k), numel(mk)));

    % right: distribution of pairwise correlations, against the non-cluster baseline
    axh = nexttile; hold(axh,'on');
    histogram(axh, cv, 40, 'Normalization','probability', 'FaceColor', cluCmap(k,:), ...
        'EdgeColor','none', 'DisplayName', 'within cluster');
    Cout = corr(Fdet(psdGrp{end}(1:min(30,numel(psdGrp{end}))), gdfr)', 'rows','pairwise');
    ov = triu(true(size(Cout,1)),1);  cvo = Cout(ov);  cvo = cvo(isfinite(cvo));
    histogram(axh, cvo, 40, 'Normalization','probability', 'FaceColor',[.5 .5 .5], ...
        'EdgeColor','none', 'FaceAlpha',0.5, 'DisplayName','outside clusters');
    xline(axh, median(cv), '-', 'Color', cluCmap(k,:), 'LineWidth', 1.5, 'HandleVisibility','off');
    box(axh,'off'); xlabel(axh, 'pairwise correlation');
    if k == 1, legend(axh,'Location','northeast'); end
    title(axh, sprintf('median %+.3f  |  coherenceGain %.2f / %.2f', ...
        median(cv), gain, sqrt(numel(mk))));
end
title(tlO, sprintf(['f=%d  cluster member traces (detrended movmedian %d ms, z-scored). ' ...
    'Weak correlation => do not average; regress the members out individually.'], f, detrendWin));


%% (10) PCA on the ODD/EVEN cross-covariance vs the ordinary covariance     [item 5]
% get_eigvector         : eig of  <x x'>   = S + N   (N = shot noise -> the
%                         spectrum flattens into a noise floor; late PCs are noise)
% get_eigvector_oddeven : eig of  <x1 x2'> = S       (odd and even pixels see the
%                         same signal with independent noise, so N cancels)
% The odd/even spectrum keeps falling and its tail goes NEGATIVE; the magnitude of
% that negative tail is the estimator's own sampling floor, so an eigenvalue clearly
% above it is real signal variance.  Run on (a) small ROIs raw, (b) small ROIs fast
% (<20 ms), (c) the whole-ROI traces used in cells (1)-(5).
spkFull  = Pr.SpikeMat(1,:) > 0;                 % NaN (motion-reject) -> false
goodFull = goodvec & ~imdilate(spkFull, ones(1,51));

oeSets = { Todd,  Teven,  Tall,  'small ROI, raw'; ...
           fTodd, fTeven, fTall, sprintf('small ROI, fast (<%d ms)', Wsm); ...
           Pr.NormalizedTrace_ch{1}, Pr.NormalizedTrace_ch{2}, Pr.NormalizedTrace_dirt, ...
                                  'whole ROI, raw (full recording)' };
oeMask = { gdfr, gdfr, goodFull };

nOE = size(oeSets,1);
Voe = cell(nOE,1);  Doe = cell(nOE,1);  Vst = cell(nOE,1);  Dst = cell(nOE,1);

get_figure('R3-oddevenPCA-spectra');
tiledlayout(1, nOE, 'TileSpacing','compact','Padding','compact');
for s = 1:nOE
    X1 = oeSets{s,1};  X2 = oeSets{s,2};  Xc = oeSets{s,3};
    % identical valid frames for both methods so the eigenvalues are comparable
    vcom = oeMask{s} & all(~isnan(X1),1) & all(~isnan(X2),1) & all(~isnan(Xc),1);
    X1(:,~vcom) = NaN;  X2(:,~vcom) = NaN;  Xc(:,~vcom) = NaN;

    [Voe{s}, Doe{s}, Soe] = get_eigvector_oddeven(X1, X2);       % shot-noise-free
    [Vst{s}, Dst{s}, Sst] = get_eigvector(Xc);                   % original method
    Dst{s} = Dst{s} / (nnz(vcom) - 1);                           % X*X' -> covariance

    ax = nexttile; hold(ax,'on'); set(ax,'XScale','log','YScale','log');
    dS = Dst{s};  dO = Doe{s};
    plot(ax, find(dS>0), dS(dS>0), '-o','Color',[.3 .3 .3],'MarkerSize',3,'LineWidth',1.2, ...
        'DisplayName','ordinary cov (get\_eigvector)');
    plot(ax, find(dO>0), dO(dO>0), '-o','Color',[.85 .2 .2],'MarkerSize',3,'LineWidth',1.2, ...
        'DisplayName','odd/even cov');
    plot(ax, find(dO<0), -dO(dO<0), ':s','Color',[.85 .2 .2],'MarkerSize',3,'LineWidth',1.0, ...
        'DisplayName','odd/even, |negative|');
    nfloor = mean(abs(dO(dO<0)));
    if ~isnan(nfloor), yline(ax, nfloor, 'k--', 'DisplayName','odd/even sampling floor'); end
    grid(ax,'on'); xlabel(ax,'PC #'); ylabel(ax,'eigenvalue (variance)'); title(ax, oeSets{s,4});
    if s == 1, legend(ax,'Location','sw'); end

    nSig = nnz(dO > max(nfloor, 0));
    fprintf('\n[R3] odd/even PCA -- %s  (%d ROIs, %d frames)\n', oeSets{s,4}, numel(dO), nnz(vcom));
    fprintf('   ordinary : %d PCs for 90%% of the variance, last/first eigenvalue = %.2g\n', ...
        find(cumsum(dS)/sum(dS) >= 0.9, 1), dS(end)/dS(1));
    fprintf('   odd/even : %d eigenvalues above the sampling floor (%.3g)\n', nSig, nfloor);
    fprintf('   shot-noise variance removed = %.3g  (%.1f%% of the ordinary total)\n', ...
        sum(dS) - sum(dO), 100*(sum(dS)-sum(dO))/sum(dS));
    ncmp = min(6, min(size(Soe,2), size(Sst,2)));
    fprintf('   |corr| of the two methods'' temporal scores, PC1..PC%d: %s\n', ncmp, ...
        sprintf('%.3f ', abs(diag(corr(Soe(:,1:ncmp), Sst(:,1:ncmp))))));
end
sgtitle(sprintf('f=%d  ordinary vs odd/even-covariance PCA (log-log eigenvalue spectra)', f));

% leading spatial modes of both methods painted on the small-footprint morphology
nModeShow = 4;  Vpair = {Vst{1}, Voe{1}};  mlab = {'ordinary cov','odd/even cov'};
get_figure('R3-oddevenPCA-modes');
tiledlayout(2, nModeShow, 'TileSpacing','compact','Padding','compact');
for meth = 1:2
    Vshow = Vpair{meth};
    for m = 1:nModeShow
        axm = nexttile;
        vmap = nan(nrow, ncol);
        vmap(inMaskC) = Vshow(labImg(inMaskC), m);
        hhm = imagesc(axm, vmap); set(hhm,'AlphaData',~isnan(vmap)); axis(axm,'image','off'); hold(axm,'on');
        colormap(axm, ccFast); caxis(axm, max(abs(Vshow(:,m)))*[-1 1]); colorbar(axm);
        contour(axm, bvBW, [0.5 0.5], 'k', 'LineWidth', 0.8);
        title(axm, sprintf('%s  PC%d', mlab{meth}, m));
    end
end
sgtitle(sprintf('f=%d  leading spatial modes, small ROIs (raw): ordinary vs odd/even covariance', f));
