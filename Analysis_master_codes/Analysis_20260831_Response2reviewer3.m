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

%% (1a) How much of the inter-ROI correlation survives smoothing?
% Reviewer 3: the moving average removes the fast subthreshold content -- does it
% also destroy the inter-ROI correlation structure the paper rests on? Three
% numbers per window, over the off-diagonal pairs of Cmat_w:
%   slope : least-squares slope of C_w vs C_raw THROUGH THE ORIGIN = the fraction
%           of each pair's non-filtered correlation that survives. This is the
%           headline "% preserved"; a per-pair ratio C_w/C_raw is useless because
%           it explodes for pairs whose raw correlation sits near zero.
%   mean  : mean(C_w)/mean(C_raw) over the same pairs, as a simple cross-check.
%   r     : Pearson r between the two matrices = how much of the STRUCTURE
%           (which pairs are correlated) survives, independent of scale.
% Two opposing effects are folded into these numbers: smoothing discards genuine
% fast shared signal (lowers r_ij) but also discards independent noise, which
% un-attenuates the correlation coefficient (raises r_ij). A value above 100%
% means the second effect won -- not that correlation was created.
iRawC = find(win_all==1);
triuC = triu(true(nROI),1);
araw  = Cmat_w{iRawC}(triuC);

preserveSlope = nan(1,nW); preserveMean = nan(1,nW); preserveR = nan(1,nW);
fprintf('\n[R3] Inter-ROI correlation preserved by smoothing (vs non-filtered):\n');
fprintf('   %-6s %9s %9s %9s %8s\n','window','slope','mean','shape r','nPairs');
for k = 1:nW
    b  = Cmat_w{k}(triuC);
    ok = ~isnan(araw) & ~isnan(b);
    preserveSlope(k) = 100 * (araw(ok)'*b(ok)) / (araw(ok)'*araw(ok));
    preserveMean(k)  = 100 * mean(b(ok)) / mean(araw(ok));
    preserveR(k)     = corr(araw(ok), b(ok));
    fprintf('   %-6s %8.1f%% %8.1f%% %9.4f %8d\n', ...
        wlab{k}, preserveSlope(k), preserveMean(k), preserveR(k), nnz(ok));
end

figure('Name','R3-corrPreserved','NumberTitle','off'); clf
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on
plot(win_all, preserveSlope, '-o','LineWidth',1.5,'DisplayName','slope through origin');
plot(win_all, preserveMean,  '-s','LineWidth',1.5,'DisplayName','ratio of means');
yline(100,'k:','HandleVisibility','off');
xlabel('moving-average window (ms)   [1 = non-filtered]');
ylabel('% of non-filtered inter-ROI correlation');
grid on; legend('Location','best'); title(sprintf('f=%d  correlation preserved', f));

nexttile; hold on
for k = 1:nW
    if k==iRawC, continue; end
    b = Cmat_w{k}(triuC); ok = ~isnan(araw) & ~isnan(b);
    scatter(araw(ok), b(ok), 6, cmapW(k,:), 'filled', ...
        'MarkerFaceAlpha',0.25, 'DisplayName',wlab{k});
end
xl = xlim; plot(xl, xl, 'k:','HandleVisibility','off'); axis square
xlabel('pairwise correlation, non-filtered'); ylabel('pairwise correlation, smoothed');
grid on; legend('Location','se'); title('per-pair comparison');

% ---- the same numbers, corrected for noise attenuation ------------------
% Values above 100% are real but reflect an estimation artefact, not smoothing
% creating correlation. Each ROI's combined trace (trC) carries its own
% independent noise, which inflates the denominator of r_ij and ATTENUATES it:
%     r_obs(i,j) = r_true(i,j) * sqrt(rho_i * rho_j)
% Smoothing removes that noise, so r_obs rises even while genuine fast shared
% signal is lost. Disattenuating both matrices before comparing isolates the
% loss of real shared signal (Spearman's correction for attenuation).
% rho is free here: diag(Cmat_w) is the odd/even split-half correlation of each
% ROI, i.e. the reliability of a HALF-pixel trace, and the Spearman-Brown
% prophecy formula converts it to the full trace that the off-diagonal uses:
%     rho_full = 2*rho_half / (1 + rho_half)
Cdis = cell(1,nW); relFull = nan(nROI,nW);
for k = 1:nW
    rh = diag(Cmat_w{k});
    rh(rh <= 0 | ~isfinite(rh)) = NaN;          % ROI too noisy to correct
    rf = 2*rh ./ (1 + rh);                      % half-pixel trace -> full trace
    relFull(:,k) = rf;
    Cdis{k} = Cmat_w{k} ./ sqrt(rf*rf');
    Cdis{k}(1:nROI+1:end) = 1;
end

adis = Cdis{iRawC}(triuC);
preserveSlopeC = nan(1,nW);
fprintf('\n[R3] Same, after correcting both matrices for noise attenuation:\n');
fprintf('   %-6s %9s %9s %9s %9s %8s\n','window','slope','mean','shape r','med rho','nPairs');
for k = 1:nW
    b  = Cdis{k}(triuC);
    ok = ~isnan(adis) & ~isnan(b);
    preserveSlopeC(k) = 100 * (adis(ok)'*b(ok)) / (adis(ok)'*adis(ok));
    fprintf('   %-6s %8.1f%% %8.1f%% %9.4f %9.3f %8d\n', wlab{k}, preserveSlopeC(k), ...
        100*mean(b(ok))/mean(adis(ok)), corr(adis(ok),b(ok)), ...
        median(relFull(:,k),'omitnan'), nnz(ok));
end
fprintf('   (disattenuated r can exceed 1 for low-reliability ROIs; median rho shows how\n');
fprintf('    much correction was applied -- rho near 1 means the correction is negligible)\n');

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

mov_sub=get_subthreshold(tovec(mov_res),max(Pr.SpikeClassMat(:,segFrames),[],1)>0,9,20);
mov_sub=toimg(mov_sub,nrow,ncol);

ftprintsmallest=divide_ROIs(Resultsmall.ftprint_smaller.*dMask,50,dMask);
dMask2=point2img(Result.SWC,7,size(Result.ref_im));
%% Show variance image
% Four panels:
%   (1) lag-1 pixel variance maps (total / subthreshold / difference) + scale bar
%   (2) per-sub-ROI map of the variance fraction retained by smoothing
%   (3) scatter: bvMask overlap of each sub-ROI  vs  that retained fraction
%       -> if the fast subthreshold content is vascular, ROIs sitting on a
%          vessel should RETAIN less (lose more) when smoothed = negative slope
%   (4) the clicky traces as dF/F0
spike_frame=find(max(Pr.SpikeClassMat(:,segFrames),[],1)>0);
silent_frame=setdiff([1:length(segFrames)],spike_frame'+[-3:3]);
ROI_fpath=fullfile(fpath{f},'bvROI_clicky.mat');

% 'Clip' is deliberately NOT used here: max(V,0) biases low-variance pixels
% upward (E[max(V,0)] > 0 where the true value is 0), which would inflate the
% VarFraction ratio below. Keep the signed maps for the statistics and clip
% only at display time.
[V_mov] = get_lag1pixelVar(mov_res(:,:,silent_frame),'DetrendWindow',500,'Clip',1);
[V_sub] = get_lag1pixelVar(mov_sub(:,:,silent_frame),'DetrendWindow',500,'Clip',1);
V_mov=maskEdge(V_mov,5,0);
V_sub=maskEdge(V_sub,5,0);
if isfile(ROI_fpath)
    load(ROI_fpath);
else
    [roi, tr]=clicky(mov_res,max(V_mov,0));
    save(ROI_fpath,'roi','tr','-v7.3');
end

edgeW  = 5;                                   % must match the maskEdge calls above
edgeOK = true(nrow,ncol);
edgeOK([1:edgeW end-edgeW+1:end],:) = false;
edgeOK(:,[1:edgeW end-edgeW+1:end]) = false;

% ---- per-sub-ROI variance ratio and blood-vessel overlap -----------------
% (the previous one-liner was  tovec(ftprintsmallest>0.*dMask) : MATLAB reads
%  "0." as the number 0, so it evaluated ftprintsmallest>(0*dMask) and dMask
%  never entered. Harmless -- dMask is already applied at the divide_ROIs call
%  above -- but the loop states the intent, uses means instead of sums, and
%  drops the maskEdge border whose zeros would otherwise dilute the ROI mean.)
bvMapL = max(Pr.bvMask,[],3);
bvBWL  = bvMapL > 0.2*max(bvMapL(:));         % same threshold as bvBW above
nSub   = size(ftprintsmallest,3);
Vm_mov = nan(nSub,1); Vm_sub = nan(nSub,1);
bvOvl  = nan(nSub,1); bvOvlW = nan(nSub,1); fracOK = zeros(nSub,1);
for k = 1:nSub
    mk  = ftprintsmallest(:,:,k) > 0;
    mok = mk & edgeOK;
    fracOK(k) = sum(mok(:))/max(sum(mk(:)),1);
    if ~any(mok(:)), continue; end
    Vm_mov(k) = mean(V_mov(mok));             % mean variance over the sub-ROI
    Vm_sub(k) = mean(V_sub(mok));
    bvOvl(k)  = mean(bvBWL(mok));             % fraction of ROI area on a vessel
    bvOvlW(k) = mean(bvMapL(mok))/max(bvMapL(:));   % intensity-weighted version
end
VarFraction = Vm_sub ./ Vm_mov;               % retained after smoothing
keep = fracOK > 0.9 & Vm_mov > 0 & isfinite(VarFraction);
fprintf('\n[R3] sub-ROIs: %d total, %d kept (>=90%% inside the field, Vm_mov>0)\n', nSub, sum(keep));

% ---- dF/F0 for the clicky ROIs ------------------------------------------
% mov_res is mean-subtracted AND spike-negative (the rest of the script uses
% -mov_res as spike-positive), so tr has ~zero mean: F0 must come from the RAW
% movie over the same polygons, and the sign is flipped here.
if exist('mov_mc','var')
    F0img = mean(mov_mc,3);
else
    warning('mov_mc not in memory - falling back to Pr.ref_im for F0 (check its scaling).');
    F0img = double(Pr.ref_im);
end
nROIc = numel(roi);
F0roi = nan(1,nROIc);
[xg,yg] = meshgrid(1:ncol,1:nrow);
for k = 1:nROIc
    inp = inpolygon(xg,yg,roi{k}(:,1),roi{k}(:,2));   % identical mask to clicky's
    F0roi(k) = mean(F0img(inp));
end
dFF = -100 * tr ./ F0roi;                     % [nFrames x nROIc], %, spike-positive

figure(256); clf;
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% ---- (1) variance maps + blood-vessel contour + scale bar ---------------
ax1 = nexttile(1,[2 1]);
imshow2([max(V_mov,0); max(V_sub,0); max(V_mov-V_sub,0)],[0 1.5]*10^4);
colormap(ax1,'gray'); hold(ax1,'on');
rowLab = {'total','subthreshold','difference'};
Bbv = bwboundaries(bvBWL);                      % blood vessels (green)
if exist('dMask2','var') && isequal(size(dMask2),[nrow ncol])
    BdM = bwboundaries(dMask2 > 0);              % SWC dendrite mask (magenta)
else
    BdM = {};
    warning('dMask missing or not %dx%d - dendrite contour skipped.', nrow, ncol);
end
for r = 1:3
    yoff = (r-1)*nrow;
    % for b = 1:numel(Bbv)
    %     plot(ax1, Bbv{b}(:,2), Bbv{b}(:,1)+yoff, 'g-','LineWidth',0.8);
    % end
    for b = 1:numel(BdM)
        plot(ax1, BdM{b}(:,2), BdM{b}(:,1)+yoff, '-','Color',[1 0 1],'LineWidth',0.8);
    end
    text(ax1, 6, yoff+12, rowLab{r}, 'Color','y','FontSize',9,'FontWeight','bold');
end
if isfield(Pr,'pixelsize') && ~isempty(Pr.pixelsize)
    pxs = Pr.pixelsize;
elseif exist('Pixelsize','var')
    pxs = Pixelsize(f);
else
    pxs = [];
end
if ~isempty(pxs)
    barUm = 100; barPx = barUm/pxs;
    if barPx > 0.4*ncol, barUm = 50; barPx = barUm/pxs; end   % keep it inside the frame
    drawScaleBar(barPx,'horizontal','Position',[ncol-12, 3*nrow-12],'Color',[1 1 1],'LineWidth',3);
    text(ax1, ncol-12-barPx/2, 3*nrow-24, sprintf('%d µm',barUm), ...
        'Color','w','HorizontalAlignment','center','FontSize',9);
else
    warning('No pixelsize found for f=%d - scale bar skipped.', f);
end
title(ax1, sprintf('f=%d  lag-1 pixel variance (clipped for display),  magenta = dMask', f),'FontSize',9);
cb = colorbar(ax1); cb.Label.String = 'Variance (count^2)';
% ---- (2) retained fraction per sub-ROI ----------------------------------
% excluded sub-ROIs are dropped from the stack, so they render as background
vfLim = [0 1];                                  % colour bounds, shared with the colorbar
show_footprint_heatmap(ftprintsmallest(:,:,keep), min(max(VarFraction(keep),vfLim(1)),vfLim(2)), vfLim(1), vfLim(2));
% show_footprint_heatmap draws an RGB image (no CData scale of its own), so the
% colorbar is attached by hand: same turbo map and the same bounds it used.
ax2 = gca;
colormap(ax2,'turbo'); clim(ax2, vfLim);
cb2 = colorbar(ax2); cb2.Label.String = 'retained';   % short: avoids colliding with panel 3

% ---- (3) bvMask overlap vs retained variance fraction -------------------
ax3 = nexttile;
xo = bvOvl(keep); yv = VarFraction(keep);
scatter(ax3, xo, yv, 26, log10(Vm_mov(keep)), 'filled','MarkerFaceAlpha',0.75);
hold(ax3,'on'); colormap(ax3,'parula');
cb = colorbar(ax3); cb.Label.String = 'log_{10} mean V\_mov';
if numel(xo) > 3
    pf = polyfit(xo,yv,1); xx = linspace(min(xo),max(xo),50);
    plot(ax3, xx, polyval(pf,xx), 'k-','LineWidth',1.5);
    [rP,pP] = corr(xo,yv);
    [rS,pS] = corr(xo,yv,'Type','Spearman');
    [rW,pW] = corr(bvOvlW(keep),yv);
    title(ax3, sprintf('r=%.2f (p=%.1e),  \\rho=%.2f (p=%.1e),  n=%d', rP,pP,rS,pS,numel(xo)),'FontSize',9);
    fprintf('[R3] bv overlap vs retained variance fraction (n=%d):\n', numel(xo));
    fprintf('   binary overlap    : Pearson r=%+.3f (p=%.2e),  Spearman rho=%+.3f (p=%.2e)\n', rP,pP,rS,pS);
    fprintf('   weighted overlap  : Pearson r=%+.3f (p=%.2e)\n', rW,pW);
    fprintf('   slope = %+.3f retained-fraction per unit overlap\n', pf(1));
end
xlabel(ax3,'sub-ROI overlap with bvMask   (fraction of ROI pixels)');
ylabel(ax3,'variance retained by 20-ms smoothing');
grid(ax3,'on'); ylim(ax3,[0 1]);

% ---- (4) clicky traces as dF/F0 -----------------------------------------
ax4 = nexttile(5,[1 2]); hold(ax4,'on');
showMs = min(3000, size(dFF,1));              % 1 frame = 1 ms (frameRate 1000)
tt   = 0:showMs-1;
dOff = 1.3*max(range(dFF(1:showMs,:),1));
for k = 1:nROIc
    yk = dFF(1:showMs,k) - (k-1)*dOff;
    plot(ax4, tt, yk, '-','Color',[0.55 0.55 0.55],'LineWidth',0.5);
    plot(ax4, tt, movmean(yk,20), '-','Color',[0.85 0.33 0.10],'LineWidth',1.2);
    text(ax4, -0.012*showMs, -(k-1)*dOff, sprintf('#%d  F_0=%.0f',k,F0roi(k)), ...
        'HorizontalAlignment','right','FontSize',8);
end
sf = spike_frame(spike_frame<=showMs);
plot(ax4, sf, 0.55*dOff*ones(size(sf)), 'k|','MarkerSize',4);
sbAmp = max(1, round(0.25*dOff));             % vertical scale bar, % dF/F
ybase = -(nROIc-1)*dOff - 0.55*dOff;
plot(ax4, 0.99*showMs*[1 1], ybase+[0 sbAmp], 'k-','LineWidth',3);
text(ax4, 0.975*showMs, ybase+sbAmp/2, sprintf('%d%%',sbAmp), ...
    'HorizontalAlignment','right','FontSize',8);
xlabel(ax4,'time (ms)'); set(ax4,'YTick',[]); box(ax4,'off');
xlim(ax4,[-0.06*showMs showMs]);
title(ax4,'clicky ROIs as dF/F_0  (grey = raw, orange = 20 ms moving average, | = spikes)','FontSize',9);
colormap(ax1,'gray');                         % re-assert: the heatmap sets turbo

fprintf('[R3] clicky ROI F0 (raw counts) and dF/F0 SD:\n');
for k = 1:nROIc
    fprintf('   #%d  F0=%8.1f   SD(dF/F0)=%5.2f%%\n', k, F0roi(k), std(dFF(:,k),'omitnan'));
end



%% ################################################################
%  (5) BLOOD-VESSEL CORRECTION (bleach-dr): PARAMETER SWEEP  (f = 20)
%
%  Reviewer 3, final step: if the short-timescale subthreshold content that
%  20-ms smoothing removes is vascular, then removing the blood signal at the
%  MOVIE level should (a) shrink the high-frequency residual, (b) leave the
%  spikes untouched, and (c) leave the subthreshold correlation matrix
%  essentially equal to the 20-ms matrix on the same frames.
%
%  What the bleach-dr pipeline does (bleach_dr.CorrectionStrategy.run_correction):
%    0. input movie = spike-POSITIVE motion residual from loadVoltChunk
%       (photobleach-normalised, mcTrace / mcTrace^2 / cross-term regressed out,
%        then sign-flipped so spikes go up).
%    1. rank1 (neural) handling -- rank1_image is the spike STA image; with
%       rank1_method = 'ols'/'kalman' the neuron's rank-1 contribution is
%       estimated and SUBTRACTED before the blood estimation, then added back at
%       the end; 'skip' leaves the movie untouched (the rank1 image is still used
%       as the projection target inside blood_correction_projection).
%    2. PMD compression of the (rank1-subtracted) movie -- block_sizes,
%       max_components, spatial/temporal_avg_factor.
%    3. ITERATIVE blood estimation: for each block size in the strategy list, a
%       band-passed (bandpass_range Hz) block-wise decomposition of the current
%       residual is projected away from the rank1 image -> one BloodArray. The
%       residual is updated and the next (larger) block size runs.
%          strategy_preset 0 -> {B}          (1 iteration)
%          strategy_preset 1 -> {B, 2B}      (2 iterations)
%          strategy_preset 2 -> {B, 2B, 3B}  (3 iterations)
%       So the preset controls how aggressively/multi-scale the blood is peeled.
%    4. corrected = PMD movie - sum(blood arrays) (+ rank1 added back).
%
%  This cell sweeps the parameters on ONE short segment and scores each setting
%  with (i) blood removal, (ii) spike preservation, (iii) vascular specificity of
%  what was removed, and (iv) the subthreshold correlation matrix.
% ################################################################
f  = 20;
Pr = PCresult{f};
if ~exist('Result','var') || ~isfield(Result,'spike')
    load(fullfile(fpath{f},'PC_Result.mat'),'Result');
end

% ---------------------------------------------------------------------------
% 5a. ALL bleach-dr PARAMETERS IN ONE PLACE
%     Each field maps 1:1 onto a bleach_dr.CorrectionStrategy constructor
%     argument (see bleach_dr/correction_strategy.py) via run_correction.py.
% ---------------------------------------------------------------------------
bdr = struct();

% -- runner / IO -------------------------------------------------------------
bdr.venvpy   = '/Users/bhlee1117/Documents/GitHub/bleach-dr/.venv/bin/python';
bdr.script   = '/Users/bhlee1117/Documents/GitHub/bleach-dr/run_correction.py';
bdr.dataset           = 'mov';           % --dataset            : movie dataset in the input .h5, (T,Ly,Lx)
bdr.rank1ImageDataset = 'rank1_image';   % --rank1-image-dataset: (Ly,Lx) spike image; '' forces bleach-dr's own spike detection
bdr.noBlood           = false;           % --no-blood           : skip writing /blood (faster, half the disk). Keep FALSE here: the blood movie is a metric.
bdr.chunkSize         = 0;               % --chunk-size         : 0 = whole segment in one pass. Each processed window must be >= 1024 frames.
bdr.chunkOverlap      = 0;               % --chunk-overlap      : frames read on each side of a chunk and trimmed after correction (edge artifacts)

% -- PMD / compression (masknmf CompressStrategy) ---------------------------
% NOTE: every bdr.* field holds a SINGLE value -- the default for one run. To
%       sweep a parameter, put the list in sweep.<field> below, not here.
bdr.blockSize     = 20;                  % --block-size     : PMD block edge in pixels; ALSO the base block of the iterative blood strategy
bdr.maxComponents = 40;                  % --max-components : max SVD rank per block
bdr.device        = 'cpu';               % --device         : 'cpu' | 'auto' | 'cuda'  (cpu on Apple Silicon)

% -- blood band / iteration --------------------------------------------------
bdr.frameRate      = 1000;               % --frame-rate      : Hz
bdr.bandpass       = [150 220];          % --bandpass-low/high: Hz band used to isolate the blood fluctuation
bdr.strategyPreset = 2;                  % --strategy-preset  : 0 = {B}, 1 = {B,2B}, 2 = {B,2B,3B} block sizes -> #iterations

% -- rank1 (neural) handling -------------------------------------------------
bdr.rank1Method = 'skip';                % --rank1-method : 'skip' | 'ols' | 'kalman'
bdr.rank1Source = 'STAmap';              % where the /rank1_image comes from:
                                         %  'STAmap'   : SS_STAmovie.mat -> max(-STAmovie,[],3), 6-px border
                                         %               (the clean, whole-recording spike template)
                                         %  'spikeSTA' : STA of the residual movie over its own spike frames
                                         %  'custom'   : use bdr.rank1Custom (nRow x nCol) as-is
bdr.rank1Custom = [];                    % only used when rank1Source = 'custom'

% -- CorrectionStrategy knobs that run_correction.py does NOT expose ---------
%    (listed here for completeness; sweeping them needs 6 extra argparse flags
%     in run_correction.py + the matching kwargs in build_strategy())
bdr.unexposed = struct( ...
    'frame_range',                [], ...    % frames used for the PMD fit (default: all)
    'frame_batch_size',           1024, ...  % PMD / rank1 batching
    'spatial_avg_factor',         1, ...     % spatial downsample before PMD
    'temporal_avg_factor',        1, ...     % temporal downsample before PMD
    'bandpass_batch_size',        1024, ...  % batching of the band-pass filter
    'manual_strategy',            [], ...    % explicit list of blood block sizes (overrides strategy_preset)
    'compute_rank1_image_on_pmd', false, ... % detect the spike image on the PMD movie instead of the raw movie
    'rank_prune_prct',            0.95);     % hard-coded in correction_iteration()

% -- test segment ------------------------------------------------------------
bdr.segSource = 'mov_res';               % 'mov_res'  : reuse the movie already in memory from section (4)
                                         %              (mov_res covers global frames segFrames)
                                         % 'chunk'    : reload one mc chunk with loadVoltChunk (matches the
                                         %              whole-movie production run in Analysis_20260806)
bdr.nSweepFrames = 1100;                  % frames of the test segment, counted from the start of the source
                                         % NOTE: bleach-dr needs >= 1024 frames per processed window (PMD +
                                         % bandpass). 500 frames WILL warn and most likely fail -- raise to
                                         % >= 1024 (e.g. 2000) if run_correction.py errors out.
bdr.fpath_mc  = fpath{f};                % folder holding mc_ShutterReg##.bin + mcTrace##.mat ('chunk' source)
bdr.mc_prefix = 'mc_ShutterReg';
bdr.seg_len   = 15000;                   % frames per mc chunk file
bdr.chunkStart= 78000;                   % first global frame of the test window ('chunk' source)
bdr.scratch   = bdr.fpath_mc;            % where the shared INPUT HDF5 lives (use a fast/local disk)

% -- corrected-movie cache ---------------------------------------------------
% One .h5 per parameter combination (/corrected + /blood for the whole test
% segment), named f<f>_<tag>_seg<digest>.h5. The digest is an MD5 of the input
% movie + rank1 image, so changing the segment or the rank1 source gives new
% files instead of silently reusing stale ones. On a rerun, a file whose stored
% provenance attributes match the requested parameters is REUSED and the python
% run is skipped -- so widening the sweep only computes the new combinations.
bdr.cacheDir  = fullfile(fpath{f}, 'BVsweep_corrected');
bdr.useCache  = true;                    % false = always recompute (existing files are overwritten)
bdr.reuseMetrics = true;                 % if the sweep .mat already holds results for THIS
                                         % exact grid and input digest, load them and skip the
                                         % loop entirely. Without this the metrics were
                                         % recomputed every run even when every movie was
                                         % cached -- which means re-reading corrected+blood
                                         % (~0.9 GB per combination) and redoing every
                                         % subthreshold / correlation / STA. Set false to
                                         % force the metrics to be recomputed.
bdr.keepCache = true;                    % false = delete each .h5 after its metrics are computed
bdr.nPreviewFrames = 200;                % frames of corrected+blood kept per combo for the section (6) movie
                                         % memory: nPreview x 180 x 552 x 4 B x 2 x nCombo (200 & 9 combos ~ 1.4 GB)
bdr.savePreview    = true;               % also write those previews into the sweep .mat (large but makes (6) standalone)

% -- SCORING PARAMETERS (how each corrected movie is judged) -----------------
%    Applied identically to the uncorrected reference and to every combination.
scr = struct();
scr.sub_win     = 20;    % ms. Smoothing window defining "fast": the correlation matrix
                         % uses the sub_win subthreshold, and the fast component is
                         % (raw - sub_win moving average). Must match sections (1)-(3).
scr.sub_dilate  = 7;     % frames dilated around each spike inside get_subthreshold
scr.minPairs    = 100;   % min valid frame-pairs per lag for get_corrTensor / nanXCorrFFT
scr.bAP_pre     = 5;     % frames before the spike in the bAP STA window
scr.bAP_post    = 15;    % frames after  the spike in the bAP STA window
scr.bAP_baseprc = 30;    % percentile of the STA used as its baseline before taking the peak
scr.bvThresh    = 0.2;   % bvMask contour level, as a fraction of max(bvMask) (vascular reference)
scr.blueDilate  = 201;   % frames dilated around blue-stim ON to define bad frames

% ---------------------------------------------------------------------------
% 5a'. THE SWEEP GRID -- add values to any line to widen the sweep
%
% GOAL: remove blood WITHOUT losing the small, fast, branch-local voltage
% signal. Section (6b) scores exactly that (bAP_min = worst branch's bAP kept).
% Sweep ONE axis at a time, in this order -- each stage fixes the winner of the
% previous one, so a 3-value stage costs 3 runs instead of 3^5 combinations:
%
%   Stage 0  no run. Just look at the rank1 image figure from (5b-ii). It is the
%            ONLY thing protecting the neuron (every blood component is
%            orthogonalised against it). If distal branches are dim in it, they
%            are the least protected -- expect them to be the worst class in 6b.
%   Stage 1  bandpass, preset=0, blockSize=20.  Higher band = less neural
%            content in the spatial template = safer. Pick the HIGHEST band that
%            still removes blood.        {[150 170],[200 220],[250 300]}
%   Stage 2  blockSize, preset=0, best band.  THE MAIN RISK AXIS: 20 px = 23 um
%            is the length of a dendritic segment, so a small block can absorb a
%            branch's own signal. Larger = safer.      {20, 30, 40, 60}
%   Stage 3  preset, best band + block.  Removal accumulates over iterations, so
%            watch bAP_min fall as you go 0 -> 1 -> 2.        {0, 1, 2}
%   Stage 4  maxComp {15,30,60} and rank1 {'skip','kalman'} only if stage 1-3
%            left no setting with bAP_min > 0.9.
%
% NOISE FLOOR FIRST: truncated_random_svd has no seed, so two identical runs
% differ. Put the SAME combination in twice and treat any metric difference
% smaller than that pair's spread as meaningless, e.g.
%   sweep.bandpass = {[200 220], [200 220], [150 170]};
% ---------------------------------------------------------------------------
sweep.bandpass  = {[150 220],[200 220]};   % --bandpass-low/high
sweep.preset    = {0, 1, 2};                           % --strategy-preset
sweep.rank1     = {'skip'};                            % --rank1-method   e.g. {'skip','kalman','ols'}
sweep.blockSize = [20 30 40];                     % --block-size     e.g. {10,20,40}
sweep.maxComp   = [20 30 40 50];                 % --max-components e.g. {15,30,60}

% -- normalise + validate the grid ------------------------------------------
% Each sweep.<field> is a CELL with one entry per value to try. A plain numeric
% vector is accepted too and expanded (bandpass: an N x 2 matrix = N bands).
% The common mistake is {[10 20 40]} -- ONE combination whose blockSize is a
% vector -- instead of {10, 20, 40}. That used to fail later with a cryptic
% matrix-multiply error, so it is caught here.
if isnumeric(sweep.bandpass)
    assert(size(sweep.bandpass,2) == 2, '[R3-5] numeric sweep.bandpass must be N x 2 ([low high] per row)');
    sweep.bandpass = num2cell(sweep.bandpass, 2)';
end
for q = {'preset','blockSize','maxComp'}
    if isnumeric(sweep.(q{1})), sweep.(q{1}) = num2cell(sweep.(q{1})); end
end
if ischar(sweep.rank1), sweep.rank1 = {sweep.rank1}; end
for j = 1:numel(sweep.bandpass)
    assert(isnumeric(sweep.bandpass{j}) && numel(sweep.bandpass{j}) == 2, ...
        '[R3-5] sweep.bandpass{%d} must be [low high]', j);
end
for q = {'preset','blockSize','maxComp'}
    for j = 1:numel(sweep.(q{1}))
        assert(isnumeric(sweep.(q{1}){j}) && isscalar(sweep.(q{1}){j}), ...
            ['[R3-5] sweep.%s{%d} has %d elements -- each cell entry must be ONE value. ' ...
             'Write {10, 20, 40} (three settings), not {[10 20 40]} (one setting with a vector).'], ...
            q{1}, j, numel(sweep.(q{1}){j}));
    end
end

fnS  = fieldnames(sweep);
idxS = cellfun(@(x) 1:numel(sweep.(x)), fnS, 'UniformOutput', false);
G    = cell(1, numel(fnS));  [G{:}] = ndgrid(idxS{:});
nCombo = numel(G{1});
P = struct('bandpass',{},'preset',{},'rank1',{},'blockSize',{},'maxComp',{},'tag',{});
for k = 1:nCombo
    for q = 1:numel(fnS), P(k).(fnS{q}) = sweep.(fnS{q}){G{q}(k)}; end
    P(k).tag = sprintf('bp%g-%g_p%d_%s_b%d_c%d', P(k).bandpass(1), P(k).bandpass(2), ...
                       P(k).preset, P(k).rank1, P(k).blockSize, P(k).maxComp);
end
% -- pre-flight: the LAST blood iteration uses (preset+1)*blockSize as its block,
%    and blood_correction_projection crashes (IndexError on an empty block grid)
%    as soon as that exceeds the smaller FOV dimension. Catch it here instead of
%    losing a python run per bad combination.
fovMin = min(size(Pr.ref_im));
for k = 1:nCombo
    Bmax = (P(k).preset + 1) * P(k).blockSize;
    if Bmax > fovMin
        error(['[R3-5] %s: preset %d needs blocks up to %d px but the FOV is only %d px ' ...
               '(blood_correction_projection will raise IndexError). ' ...
               'Use blockSize <= %d with preset %d.'], ...
               P(k).tag, P(k).preset, Bmax, fovMin, floor(fovMin/(P(k).preset+1)), P(k).preset);
    end
    % The blood model on the FIRST (smallest) block has maxComp components for
    % blockSize^2 pixels. As that ratio approaches 1 the model can represent
    % almost anything inside the block -- including a dendritic branch's own
    % fast signal -- and the only thing still protecting the neuron is the
    % single rank1-image direction that gets orthogonalised out.
    expressiveness = P(k).maxComp / P(k).blockSize^2;
    if expressiveness > 0.25
        warning(['[R3-5] %s: %d components for %d px in a %dx%d block (%.0f%% of the block''s ' ...
                 'pixel space). The blood model is nearly unconstrained at this scale and is ' ...
                 'likely to absorb branch-local fast signal -- watch bAP_min in section (6b).'], ...
                 P(k).tag, P(k).maxComp, P(k).blockSize^2, P(k).blockSize, P(k).blockSize, ...
                 100*expressiveness);
    end
end
fprintf('\n[R3-5] %d parameter combination(s) to run (FOV %d px -> max block %d px)\n', nCombo, fovMin, fovMin);
for k = 1:nCombo
    fprintf('        %-34s blood blocks: %s\n', P(k).tag, ...
            mat2str(P(k).blockSize * (1:P(k).preset+1)));
end

% ---------------------------------------------------------------------------
% 5b-i. TEST SEGMENT
%     bleach-dr expects a SPIKE-POSITIVE motion residual (that is what
%     loadVoltChunk returns as movR = -mov_res). Section (4)'s mov_res is the
%     same movie with the raw, negative-going voltage sign, so it is negated
%     here; it also skips the photobleach normalisation loadVoltChunk applies,
%     which is fine for a short test but is a difference from the production run.
% ---------------------------------------------------------------------------
switch bdr.segSource
    case 'mov_res'      % ---- reuse section (4): first nSweepFrames of mov_res ----
        assert(exist('mov_res','var')==1, 'Run section (4) first (mov_res / segFrames must exist).');
        nTS         = min(bdr.nSweepFrames, size(mov_res,3));
        sweepFrames = segFrames(1:nTS);                         % global frame indices
        movS        = single(-mov_res(:,:,1:nTS));              % spike-positive

    case 'chunk'        % ---- reload one mc chunk (matches the production run) ----
        load(fullfile(fpath{f},'output_data.mat'),'Device_Data');
        szS  = double(Device_Data{1,3}.ROI([2 4]));             % szS(1)=nCol, szS(2)=nRow
        sweepFrames = bdr.chunkStart + (0:bdr.nSweepFrames-1);
        jSeg = ceil(sweepFrames(1)/bdr.seg_len);
        assert(ceil(sweepFrames(end)/bdr.seg_len) == jSeg, ...
               'test window must stay inside one %d-frame mc chunk', bdr.seg_len);
        locS = sweepFrames - (jSeg-1)*bdr.seg_len;              % local index inside the chunk
        fprintf('[R3-5] loading mc chunk %d (%s) ...\n', jSeg, bdr.mc_prefix);
        movChunk  = loadVoltChunk(bdr.fpath_mc, jSeg, szS, bdr.seg_len, bdr.mc_prefix);  % spike-positive
        movS      = single(movChunk(:,:,locS));
        movChunkSpk = find(Result.spike(1, (jSeg-1)*bdr.seg_len + (1:size(movChunk,3))) > 0);
        nTS = size(movS,3);

    otherwise
        error('bdr.segSource must be ''mov_res'' or ''chunk''');
end
[nRowS, nColS, ~] = size(movS);
if nTS < 1024
    warning(['[R3-5] test segment has %d frames; bleach-dr requires >= 1024 per processed ' ...
             'window (PMD + bandpass) and run_correction.py will likely fail. ' ...
             'Set bdr.nSweepFrames >= 1024 if it does.'], nTS);
end

% ---------------------------------------------------------------------------
% 5b-ii. rank1 (spike) IMAGE  -> written to /rank1_image in the input .h5
%     This is the spatial template bleach-dr projects the blood estimate AWAY
%     from (and, for rank1_method ~= 'skip', the image whose rank-1 weight is
%     subtracted before / added back after the blood estimation). It must be
%     spike-POSITIVE and the same (nRow, nCol) as the movie.
%     'STAmap' is the preferred one: the SS spike-triggered-average movie
%     (SS_STAmovie.mat), peak over time, sign-flipped -- far cleaner than an STA
%     built from the few spikes inside a short test window. The 6-pixel border
%     is zero because the STA movie is cropped by `bound`.
% ---------------------------------------------------------------------------
switch bdr.rank1Source
    case 'STAmap'       % ---- SS_STAmovie.mat: max over time of -STAmovie ----
        STAmovie  = importdata(fullfile(fpath{f},'SS_STAmovie.mat'));   % SS_STA, 169 x 541 x 51 for f=20
        STAmap    = max(-STAmovie, [], 3);                              % peak of the spike-positive STA
        assert(size(STAmap,1) <= nRowS && size(STAmap,2) <= nColS, ...
               '[R3-5] SS_STAmovie (%dx%d) is larger than the movie frame (%dx%d)', ...
               size(STAmap,1), size(STAmap,2), nRowS, nColS);
        % centre it in the full frame; the STA movie is cropped by `bound`, so the
        % offset comes out as 5 -> rows 6:end-6, cols 6:end-6 (180x552 <- 169x541)
        r0 = floor((nRowS - size(STAmap,1))/2);
        c0 = floor((nColS - size(STAmap,2))/2);
        spikeImgS = zeros(nRowS, nColS, 'single');
        spikeImgS(r0 + (1:size(STAmap,1)), c0 + (1:size(STAmap,2))) = STAmap;
        rank1Info = sprintf('SS_STAmovie.mat, max(-STA), %d/%d-px border', r0, c0);

    case 'spikeSTA'     % ---- STA of the residual movie over its own spike frames ----
        if strcmp(bdr.segSource,'mov_res')
            spkSTA = find(Result.spike(1, segFrames) > 0);       % wider than the test window
            assert(numel(spkSTA) >= 5, '[R3-5] only %d spikes in segFrames', numel(spkSTA));
            spikeImgS = single(mean(-mov_res(:,:,spkSTA), 3));
        else
            assert(numel(movChunkSpk) > 0, '[R3-5] no soma spikes inside mc chunk %d', jSeg);
            spikeImgS = single(mean(movS(:,:,movChunkSpk(movChunkSpk<=nTS)), 3));
            spkSTA    = movChunkSpk;
        end
        rank1Info = sprintf('STA over %d spike frames', numel(spkSTA));

    case 'custom'       % ---- whatever you put in bdr.rank1Custom ----
        assert(isequal(size(bdr.rank1Custom), [nRowS nColS]), ...
               '[R3-5] bdr.rank1Custom must be %d x %d', nRowS, nColS);
        spikeImgS = single(bdr.rank1Custom);
        rank1Info = 'bdr.rank1Custom';

    otherwise
        error('bdr.rank1Source must be ''STAmap'', ''spikeSTA'' or ''custom''');
end
if strcmp(bdr.segSource,'chunk'), clear movChunk; end

fprintf('[R3-5] segment %d x %d x %d (global frames %d:%d), rank1 image = %s\n', ...
        nRowS, nColS, nTS, sweepFrames(1), sweepFrames(end), rank1Info);
figure('Name','R3-BVsweep-rank1','NumberTitle','off'); clf
imagesc(spikeImgS); axis image off; colormap(gray); colorbar
title(sprintf('f=%d  rank1 / spike image (%s)', f, rank1Info));

% footprints in the SORTED ROI order used everywhere else + odd/even channels
oddPixS   = mod((1:nRowS)' + (1:nColS), 2) == 1;                 % checkerboard
ftS       = Pr.Ftprnts(:,:,Pr.sortDist);
ft_vecS   = tovec(ftS);                                          % (pix x nROIs)
ft_oddS   = tovec(ftS .*  oddPixS);
ft_evenS  = tovec(ftS .* ~oddPixS);
ft_somaS  = tovec(Pr.Ftprnts(:,:,1));                            % soma footprint (metrics)
nROIs     = size(ft_vecS, 2);

% frame selections inside the segment
spkS     = Result.spike(1, sweepFrames) > 0;
goodvecS = ~imdilate(Pr.BlueStim(sweepFrames) > 0, ones(1,scr.blueDilate));

% vascular reference maps
bvMapS   = max(Pr.bvMask, [], 3);
bvBWS    = bvMapS > scr.bvThresh*max(bvMapS(:));
bvpixS   = bvMapS(:) > 0;

% ---------------- reference = UNCORRECTED segment ---------------------------
vecS   = tovec(movS)';                                            % (nT x pix)
trC_in = (vecS * ft_vecS)';   tr1_in = (vecS * ft_oddS)';   tr2_in = (vecS * ft_evenS)';
soma_in = vecS * ft_somaS;    soma_in = soma_in - median(soma_in);

% two correlation matrices for the reference: NO temporal smoothing (window 1)
% and the sub_win-smoothed one. The 20-ms matrix is the paper's; the unsmoothed
% one still contains the blood, and section (6b') asks whether blood-correcting
% the movie makes the UNSMOOTHED matrix look like the smoothed original.
subRaw_in  = get_subthreshold(trC_in, spkS, scr.sub_dilate, 1);
sub1Raw_in = get_subthreshold(tr1_in, spkS, scr.sub_dilate, 1);
sub2Raw_in = get_subthreshold(tr2_in, spkS, scr.sub_dilate, 1);
Cm_rin      = get_corrTensor(sub1Raw_in, sub2Raw_in, subRaw_in, goodvecS, 0, scr.minPairs);
Cmat_raw_in = Cm_rin(:,:,1);                       % uncorrected, NO smoothing
sub1_in   = get_subthreshold(tr1_in, spkS, scr.sub_dilate, scr.sub_win);
sub2_in   = get_subthreshold(tr2_in, spkS, scr.sub_dilate, scr.sub_win);
subC_in   = get_subthreshold(trC_in, spkS, scr.sub_dilate, scr.sub_win);
Cm_in     = get_corrTensor(sub1_in, sub2_in, subC_in, goodvecS, 0, scr.minPairs);
Cmat_in   = Cm_in(:,:,1);                          % uncorrected, sub_win smoothed <- the reference
% % of the raw subthreshold variance that sub_win smoothing removes (as in section 2)
rem_in = 100 * (1 - var(subC_in(:,goodvecS),0,2,'omitnan') ./ var(subRaw_in(:,goodvecS),0,2,'omitnan'));
% high-frequency (raw - sub_win) spatial variance and its vascular specificity
hfvar_in = var(movS - movmean(movS, scr.sub_win, 3), 0, 3);
hfratio_in = mean(hfvar_in(bvpixS)) / mean(hfvar_in(~bvpixS));

% --- BRANCH-LEVEL fast signal (the thing the correction could eat) ----------
% Two per-ROI readouts, both computed identically before and after correction:
%
%  (a) bAP amplitude -- a KNOWN real, fast (~2 ms), branch-local voltage event.
%      Positive control: if the bAP survives on a branch, fast branch-local
%      signal survives there. This is the hard requirement.
%  (b) fast cross-channel POWER -- variance of the (raw - 20 ms) component
%      measured as the odd x even channel COVARIANCE, so independent shot noise
%      contributes zero. It is the total real fast power per ROI: blood +
%      genuine fast dendritic signal. Correction should take out the blood part
%      of it while (a) stays at 1.
bAP_in     = get_STA(trC_in, spkS, scr.bAP_pre, scr.bAP_post);
bAP_in     = bAP_in - prctile(bAP_in, scr.bAP_baseprc, 2);
bAPamp_in  = max(bAP_in, [], 2);                                  % nROIs x 1
fast1_in   = tr1_in - movmean(tr1_in, scr.sub_win, 2);
fast2_in   = tr2_in - movmean(tr2_in, scr.sub_win, 2);
fastPow_in = mean(fast1_in(:,goodvecS).*fast2_in(:,goodvecS), 2, 'omitnan');   % nROIs x 1
labelvecS  = Pr.labelvec(:);            % 1 basal, 2 soma, 3 trunk, 4 oblique, 5 distal
classNm    = {'Basal','Soma','Trunk','Oblique','Distal'};
dendaxisS  = Pr.dendaxis(:);            % dendritic position (um), sorted order

fprintf(['[R3-5] uncorrected reference: %d good frames, %d soma spikes, ' ...
         '%%var removed by %d-ms smoothing = %.1f%%, HF var in/out bvMask = %.2f\n'], ...
        nnz(goodvecS), nnz(spkS), scr.sub_win, mean(rem_in,'omitnan'), hfratio_in);
if nnz(goodvecS) < 200
    warning('[R3-5] only %d blue-off frames in the test window -- the correlation matrix will be noisy', nnz(goodvecS));
end
if nnz(spkS) < 3
    warning('[R3-5] only %d soma spikes in the test window -- spike_slope will be NaN', nnz(spkS));
end

% ---------------------------------------------------------------------------
% 5c. SWEEP
%     Each combination's corrected+blood movie is written to its OWN file under
%     bdr.cacheDir, named by the parameter tag plus a digest of the input movie
%     and rank1 image. On a rerun, a file whose stored provenance attributes
%     match the requested parameters is reused and the python run is skipped.
%     Delete a file (or set bdr.useCache=false) to force recomputation.
% ---------------------------------------------------------------------------
in_h5S    = fullfile(bdr.scratch, 'bv_sweep_in.h5');
sweepFile = fullfile(save_figto, sprintf('R3_bleachdr_sweep_f%d.mat', f));
if ~isfolder(bdr.cacheDir), mkdir(bdr.cacheDir); end

% input digest: any change to the segment or the rank1 image invalidates the cache
stepSig  = max(1, floor(numel(movS)/1e5));          % subsample: ~1e5 values is plenty
mdS      = java.security.MessageDigest.getInstance('MD5');
sigVecS  = [reshape(movS(1:stepSig:end),[],1); spikeImgS(:)];   % both as columns
mdS.update(typecast(single(sigVecS), 'int8'));                  % int8 == Java byte
segSig   = sprintf('%02x', mod(double(mdS.digest()), 256));
segSig   = segSig(1:10);
cacheFn  = @(tag) fullfile(bdr.cacheDir, sprintf('f%d_%s_seg%s.h5', f, tag, segSig));
fprintf('[R3-5] cache dir: %s   (input digest %s)\n', bdr.cacheDir, segSig);
fprintf('[R3-5] ~%.1f GB per combination (%d frames, corrected+blood, single)\n', ...
        2*nTS*nRowS*nColS*4/2^30, nTS);

% ---- level 1 reuse: are the METRICS for this exact grid already on disk? ---
% The movie cache only ever skipped the python call; everything after it -- the
% ~0.9 GB read-back per combination and every subthreshold / correlation / STA --
% still ran. This check skips all of that when nothing has changed.
metricsReused = false;
if bdr.reuseMetrics && isfile(sweepFile)
    try
        Sm = load(sweepFile, 'SW','segSig','sweepFrames');
        sameGrid = isfield(Sm,'SW') && numel(Sm.SW) == nCombo && ...
                   isequal({Sm.SW.tag}, {P.tag});
        sameSeg  = isfield(Sm,'segSig') && strcmp(Sm.segSig, segSig) && ...
                   isfield(Sm,'sweepFrames') && isequal(Sm.sweepFrames(:)', sweepFrames(:)');
        if sameGrid && sameSeg
            load(sweepFile);                       % SW, Cmat_*, prev_*, everything
            metricsReused = true;
            fprintf(['[R3-5] REUSING the saved metrics: same %d-combination grid and the ' ...
                     'same input digest (%s).\n       Nothing is recomputed. ' ...
                     'Set bdr.reuseMetrics = false to force it.\n'], nCombo, segSig);
        elseif ~sameGrid
            fprintf('[R3-5] sweep .mat holds a different grid -> recomputing\n');
        else
            fprintf('[R3-5] sweep .mat was made from a different input segment -> recomputing\n');
        end
    catch ME
        fprintf('[R3-5] could not reuse the sweep .mat (%s) -> recomputing\n', ME.message);
    end
end

if ~metricsReused
% ---- level 2 reuse: which combinations still need a python run? ------------
% Checked BEFORE anything is written, so the ~%.1f GB input file is only created
% when at least one combination actually needs it.
cacheOK = false(1, nCombo);
for k = 1:nCombo
    fnk = cacheFn(P(k).tag);
    if ~(bdr.useCache && isfile(fnk)), continue; end
    try
        att = @(n) h5readatt(fnk, '/', n);
        okk = att('bandpass_low')  == P(k).bandpass(1) && ...
              att('bandpass_high') == P(k).bandpass(2) && ...
              att('block_size')    == P(k).blockSize   && ...
              att('max_components')== P(k).maxComp     && ...
              att('strategy_preset')== P(k).preset     && ...
              att('frame_rate')    == bdr.frameRate    && ...
              strcmp(strtrim(char(att('rank1_method'))), P(k).rank1);
        if okk && ~bdr.noBlood
            okk = any(strcmp({h5info(fnk).Datasets.Name}, 'blood'));
        end
        cacheOK(k) = okk;
    catch
        cacheOK(k) = false;
    end
end
fprintf('[R3-5] %d/%d corrected movies already cached; %d need a python run\n', ...
        nnz(cacheOK), nCombo, nnz(~cacheOK));

% write the shared input ONLY if something has to be computed
if any(~cacheOK)
    if exist(in_h5S,'file'), delete(in_h5S); end
    fprintf('[R3-5] writing the shared input h5 (%.2f GB) ...\n', nTS*nRowS*nColS*4/2^30);
    mov_outS   = permute(single(movS),      [2 1 3]);      % -> h5 reads (nT,nRow,nCol)
    spike_outS = permute(single(spikeImgS), [2 1]);        % -> h5 reads (nRow,nCol)
    h5create(in_h5S,'/mov',         size(mov_outS),  'Datatype','single'); h5write(in_h5S,'/mov',         mov_outS);
    h5create(in_h5S,'/rank1_image', size(spike_outS),'Datatype','single'); h5write(in_h5S,'/rank1_image', spike_outS);
    clear mov_outS spike_outS
else
    fprintf('[R3-5] every movie is cached -- the input h5 is not written at all\n');
end

% one row per combination (index-aligned with P / Cmat_sw); a failed run stays NaN
SW = repmat(struct('tag','','bandpass',[NaN NaN],'preset',NaN,'rank1','','blockSize',NaN,'maxComp',NaN, ...
            'blood_removed',NaN,'spike_slope',NaN,'spike_corr',NaN, ...
            'hf_removed',NaN,'hfratio_corr',NaN,'bloodratio',NaN,'hf_vs_blood',NaN, ...
            'rem_smooth',NaN,'corr_vs_in',NaN,'seconds',NaN, ...
            'bAP_ratio',[],'fast_keep',[],'bAP_min',NaN,'bAP_median',NaN, ...
            'bAP_worstclass_val',NaN,'bAP_worstclass','', ...
            'corrRaw_vs_in20',NaN,'corrRaw_vs_inRaw',NaN,'corrRaw_absdiff',NaN, ...
            'cacheFile','','fromCache',false), 1, nCombo);
for k = 1:nCombo
    SW(k).tag = P(k).tag;          SW(k).bandpass  = P(k).bandpass;  SW(k).preset  = P(k).preset;
    SW(k).rank1 = P(k).rank1;      SW(k).blockSize = P(k).blockSize; SW(k).maxComp = P(k).maxComp;
end
Cmat_sw     = cell(1,nCombo);   % corrected, sub_win smoothed
Cmat_raw_sw = cell(1,nCombo);   % corrected, NO temporal smoothing
blood_tr   = cell(1,nCombo);
prev_corr  = cell(1,nCombo);    % first nPreviewFrames of the corrected movie, per combo
prev_blood = cell(1,nCombo);    % ... and of the blood movie (section 6 viewer)

for k = 1:nCombo
    p       = P(k);
    out_h5S = cacheFn(p.tag);
    fprintf('\n[R3-5] (%d/%d) %s ...\n', k, nCombo, p.tag);

    useCached = cacheOK(k);          % decided in the pre-pass above

    if useCached
        d = dir(out_h5S);
        fprintf('   REUSING %s (%.1f GB, written %s)\n', out_h5S, d.bytes/2^30, ...
                string(datetime(d.datenum,'ConvertFrom','datenum','Format','yyyy-MM-dd HH:mm')));
        elapsedS = 0;
    else
        if isfile(out_h5S), delete(out_h5S); end
        if bdr.noBlood, bloodflagS = '--no-blood'; else, bloodflagS = ''; end
        cmdS = sprintf(['"%s" "%s" --input "%s" --output "%s" ' ...
            '--dataset %s --rank1-image-dataset "%s" ' ...
            '--frame-rate %g --bandpass-low %g --bandpass-high %g ' ...
            '--block-size %d --max-components %d --strategy-preset %d ' ...
            '--rank1-method %s --device %s --chunk-size %d --chunk-overlap %d %s'], ...
            bdr.venvpy, bdr.script, in_h5S, out_h5S, ...
            bdr.dataset, bdr.rank1ImageDataset, ...
            bdr.frameRate, p.bandpass(1), p.bandpass(2), ...
            p.blockSize, p.maxComp, p.preset, ...
            p.rank1, bdr.device, bdr.chunkSize, bdr.chunkOverlap, bloodflagS);

        tk = tic;
        [statusS, outS] = system(cmdS, '-echo');
        elapsedS = toc(tk);
        if statusS ~= 0
            warning('[R3-5] python failed on %s:\n%s', p.tag, outS);
            if isfile(out_h5S), delete(out_h5S); end     % never cache a partial run
            continue
        end
    end
    SW(k).cacheFile = out_h5S;
    SW(k).fromCache = useCached;

    corrS = permute(h5read(out_h5S,'/corrected'), [2 1 3]);        % (nRow,nCol,nT) single
    if ~bdr.noBlood
        bloodS = permute(h5read(out_h5S,'/blood'), [2 1 3]);
    else
        bloodS = movS - corrS;                                     % blood == what was taken out
    end

    % ---- traces on the corrected movie (combined + odd/even channels) ------
    vecC    = tovec(corrS)';
    trC_c   = (vecC * ft_vecS)';   tr1_c = (vecC * ft_oddS)';   tr2_c = (vecC * ft_evenS)';
    soma_c  = vecC * ft_somaS;     soma_c = soma_c - median(soma_c);
    blood_tr{k} = (tovec(bloodS)' * ft_somaS)';                    % soma-projected blood trace

    % ---- headline metrics (same definitions as sweep_parameters.py) --------
    bloodRemoved = 1 - var(soma_c) / var(soma_in);                 % >0 = variance taken out
    spikeCorr    = corr(soma_in, soma_c);
    spL          = find(spkS);
    if numel(spL) >= 3                                             % slope through origin at spike peaks
        x = soma_in(spL);  y = soma_c(spL);
        spikeSlope = (x'*y) / (x'*x);                              % ~1 = spikes preserved
    else
        spikeSlope = NaN;
    end

    % ---- is what was removed VASCULAR? (movie level) ----------------------
    hfvar_c   = var(corrS  - movmean(corrS,  scr.sub_win, 3), 0, 3);   % HF left after correction
    bloodVarS = var(bloodS - mean(bloodS,3), 0, 3);                % spatial map of what was removed
    hfRemoved   = 1 - mean(hfvar_c(:)) / mean(hfvar_in(:));        % fraction of HF power removed
    hfratio_c   = mean(hfvar_c(bvpixS))   / mean(hfvar_c(~bvpixS));% HF still concentrated on vessels?
    bloodRatio  = mean(bloodVarS(bvpixS)) / mean(bloodVarS(~bvpixS)); % removed power in/out bvMask
    hfVsBlood   = corr2(hfvar_in, bloodVarS);                      % does it match the HF map?

    % ---- subthreshold correlation matrix on the corrected traces ----------
    subRaw_c  = get_subthreshold(trC_c, spkS, scr.sub_dilate, 1);
    sub1Raw_c = get_subthreshold(tr1_c, spkS, scr.sub_dilate, 1);
    sub2Raw_c = get_subthreshold(tr2_c, spkS, scr.sub_dilate, 1);
    Cm_rc     = get_corrTensor(sub1Raw_c, sub2Raw_c, subRaw_c, goodvecS, 0, scr.minPairs);
    Cmat_raw_sw{k} = Cm_rc(:,:,1);                 % corrected, NO temporal smoothing
    sub1_c   = get_subthreshold(tr1_c, spkS, scr.sub_dilate, scr.sub_win);
    sub2_c   = get_subthreshold(tr2_c, spkS, scr.sub_dilate, scr.sub_win);
    subC_c   = get_subthreshold(trC_c, spkS, scr.sub_dilate, scr.sub_win);
    Cm_c     = get_corrTensor(sub1_c, sub2_c, subC_c, goodvecS, 0, scr.minPairs);
    Cmat_sw{k} = Cm_c(:,:,1);                      % corrected, sub_win smoothed
    rem_c    = 100 * (1 - var(subC_c(:,goodvecS),0,2,'omitnan') ./ var(subRaw_c(:,goodvecS),0,2,'omitnan'));

    tri1S = triu(true(nROIs),1);
    rOff  = @(A,B) corr(A(tri1S & ~isnan(A) & ~isnan(B)), B(tri1S & ~isnan(A) & ~isnan(B)));
    corrVsIn = rOff(Cmat_sw{k}, Cmat_in);            % smoothed corrected vs smoothed original
    % THE KEY COMPARISON: does blood-correcting the movie make the UNSMOOTHED
    % correlation matrix equal to the 20-ms-smoothed ORIGINAL? If yes, the 20-ms
    % smoothing in the paper was doing nothing but removing blood.
    SW(k).corrRaw_vs_in20  = rOff(Cmat_raw_sw{k}, Cmat_in);
    SW(k).corrRaw_vs_inRaw = rOff(Cmat_raw_sw{k}, Cmat_raw_in);
    dRaw = Cmat_raw_sw{k} - Cmat_in;
    SW(k).corrRaw_absdiff  = mean(abs(dRaw(tri1S)), 'omitnan');

    % ---- BRANCH-LEVEL survival: did the fast signal on each branch survive? --
    bAP_c     = get_STA(trC_c, spkS, scr.bAP_pre, scr.bAP_post);
    bAP_c     = bAP_c - prctile(bAP_c, scr.bAP_baseprc, 2);
    bAPamp_c  = max(bAP_c, [], 2);
    fast1_c   = tr1_c - movmean(tr1_c, scr.sub_win, 2);
    fast2_c   = tr2_c - movmean(tr2_c, scr.sub_win, 2);
    fastPow_c = mean(fast1_c(:,goodvecS).*fast2_c(:,goodvecS), 2, 'omitnan');

    SW(k).bAP_ratio  = bAPamp_c  ./ bAPamp_in;      % nROIs x 1, want ~1 EVERYWHERE
    SW(k).fast_keep  = fastPow_c ./ fastPow_in;     % nROIs x 1, blood removal -> < 1
    SW(k).bAP_min    = min(SW(k).bAP_ratio);        % worst branch
    SW(k).bAP_median = median(SW(k).bAP_ratio, 'omitnan');
    % worst class = the dendrite class whose bAP suffered most (usually distal)
    cls = unique(labelvecS(~isnan(labelvecS)))';
    cmed = arrayfun(@(c) median(SW(k).bAP_ratio(labelvecS==c),'omitnan'), cls);
    [SW(k).bAP_worstclass_val, wi] = min(cmed);
    SW(k).bAP_worstclass = classNm{cls(wi)};

    SW(k).blood_removed = bloodRemoved;  SW(k).spike_slope = spikeSlope;  SW(k).spike_corr = spikeCorr;
    SW(k).hf_removed    = hfRemoved;     SW(k).hfratio_corr = hfratio_c;  SW(k).bloodratio = bloodRatio;
    SW(k).hf_vs_blood   = hfVsBlood;     SW(k).rem_smooth = mean(rem_c,'omitnan');
    SW(k).corr_vs_in    = corrVsIn;      SW(k).seconds    = elapsedS;

    % ---- keep the first nPreviewFrames for the side-by-side movie (section 6) --
    nPv = min(bdr.nPreviewFrames, nTS);
    prev_corr{k}  = corrS(:,:,1:nPv);
    prev_blood{k} = bloodS(:,:,1:nPv);

    fprintf(['        blood_rm=%.3f  spk_slope=%.3f  spk_r=%.3f | HF_rm=%.3f  ' ...
             'blood_in/out=%.2f  r(HF,blood)=%.2f | r(Cmat,Cmat_in)=%.3f  (%.0fs)\n'], ...
        bloodRemoved, spikeSlope, spikeCorr, hfRemoved, ...
        bloodRatio, hfVsBlood, corrVsIn, elapsedS);
    fprintf(['        >> BRANCH bAP kept: median=%.3f  min=%.3f  worst class=%s (%.3f) | ' ...
             'fast power kept: median=%.3f\n'], ...
        SW(k).bAP_median, SW(k).bAP_min, SW(k).bAP_worstclass, SW(k).bAP_worstclass_val, ...
        median(SW(k).fast_keep,'omitnan'));

    if ~bdr.keepCache && isfile(out_h5S), delete(out_h5S); end
    clear corrS bloodS vecC
end
if exist(in_h5S,'file'), delete(in_h5S); end   % the shared input is cheap to rebuild
end                                            % <- end of "if ~metricsReused"

if bdr.keepCache
    dCache = dir(fullfile(bdr.cacheDir, sprintf('f%d_*_seg%s.h5', f, segSig)));
    fprintf('\n[R3-5] corrected movies kept in %s\n', bdr.cacheDir);
    fprintf('        %d file(s), %.1f GB total. Reload one with:\n', numel(dCache), sum([dCache.bytes])/2^30);
    fprintf('        mov = permute(h5read(SW(k).cacheFile,''/corrected''), [2 1 3]);\n');
end

% ---------------------------------------------------------------------------
% 5d. WHAT WAS SWEPT AND SAVED
%     One row per combination: the parameters that went to run_correction.py,
%     the blood block ladder they imply, where the corrected movie is on disk,
%     and whether this run computed it or reused a cached file.
% ---------------------------------------------------------------------------
okRun       = ~isnan([SW.blood_removed]);                       % combinations that completed
sweptFields = fnS(cellfun(@(x) numel(sweep.(x)) > 1,  fnS));     % axes that actually vary
fixedFields = fnS(cellfun(@(x) isscalar(sweep.(x)),   fnS));     % axes held constant
fprintf('\n[R3-5] ===== SWEPT PARAMETERS =====\n');
fprintf('   segment : %s frames %d:%d (%d frames), rank1 = %s, digest %s\n', ...
        bdr.segSource, sweepFrames(1), sweepFrames(end), nTS, bdr.rank1Source, segSig);
if isempty(sweptFields)
    fprintf('   varying : (none -- single combination)\n');
else
    fprintf('   varying : %s\n', strjoin(sweptFields', ', '));
end
fixedStr = cell(1, numel(fixedFields));
for q = 1:numel(fixedFields)
    v = sweep.(fixedFields{q}){1};
    if ischar(v), vs = v; else, vs = mat2str(v); end
    fixedStr{q} = sprintf('%s=%s', fixedFields{q}, vs);
end
fprintf('   fixed   : %s\n', strjoin(fixedStr, '  '));

status  = repmat({'FAILED'}, 1, nCombo);
status(okRun &  [SW.fromCache]) = {'cached'};
status(okRun & ~[SW.fromCache]) = {'computed'};
gb = nan(1,nCombo);
for k = 1:nCombo
    if ~isempty(SW(k).cacheFile) && isfile(SW(k).cacheFile)
        dd = dir(SW(k).cacheFile);  gb(k) = dd.bytes/2^30;
    end
end
sweptTable = table( (1:nCombo)', string({SW.tag})', ...
    vertcat(SW.bandpass), [SW.preset]', string({SW.rank1})', [SW.blockSize]', [SW.maxComp]', ...
    arrayfun(@(k) string(mat2str(P(k).blockSize*(1:P(k).preset+1))), (1:nCombo)'), ...
    string(status)', [SW.seconds]', gb(:), string({SW.cacheFile})', ...
    'VariableNames', {'k','tag','bandpass','preset','rank1','blockSize','maxComp', ...
                      'bloodBlocks','status','seconds','GB','file'});
disp(sweptTable(:, {'k','tag','bandpass','preset','rank1','blockSize','maxComp', ...
                    'bloodBlocks','status','seconds','GB'}));
fprintf('   corrected movies -> %s\n', bdr.cacheDir);

% ---------------------------------------------------------------------------
% 5e. SAVE — everything section 6 needs, so the viewer can run after a restart
% ---------------------------------------------------------------------------
if ~metricsReused
prev_in  = movS(:,:,1:min(bdr.nPreviewFrames,nTS));   % uncorrected reference preview
saveVars = {'SW','P','bdr','scr','sweep','sweptTable','segSig','sweepFrames','okRun','nROIs', ...
            'Cmat_in','Cmat_sw','Cmat_raw_in','Cmat_raw_sw','rem_in','hfratio_in','blood_tr','soma_in', ...
            'spikeImgS','bvBWS','bvpixS','Dpair','nTS', ...
            'bAPamp_in','fastPow_in','labelvecS','classNm','bAP_in','dendaxisS'};
if bdr.savePreview, saveVars = [saveVars, {'prev_in','prev_corr','prev_blood'}]; end
save(sweepFile, saveVars{:}, '-v7.3');
else
    fprintf('[R3-5] metrics reused -- the sweep .mat is left untouched\n');
end

if ~any(okRun)
    error('[R3-5] every parameter combination failed -- see the python output above');
end
fprintf('\n[R3-5] sweep done: %d/%d combinations succeeded (%d reused from cache, %d newly computed) -> %s\n', ...
        nnz(okRun), nCombo, nnz([SW.fromCache]), nnz(okRun & ~[SW.fromCache]), sweepFile);
fprintf('[R3-5] run section (6) to inspect the results.\n');

%% ################################################################
%  (6) BLEACH-DR SWEEP RESULTS  (viewer for what section (5) produced)
%
%  Section (5) runs the sweep and saves everything to
%  <save_figto>/R3_bleachdr_sweep_f<f>.mat. This section only READS that and
%  can be run on its own after a restart -- it reloads the .mat if the sweep
%  variables are not in the workspace.
%
%   6a. summary table + metric scatters + blood spectra
%   6b. DID THE BRANCH-LOCAL FAST SIGNAL SURVIVE?  (the question that matters)
%   6c. UNSMOOTHED corrected matrix vs the SMOOTHED original -- heatmaps
%   6d. best combination -> correlation matrices vs uncorrected, corr vs distance
%   6e. stacked movie of the first bdr.nPreviewFrames frames, chosen conditions
%   6f. interactive viewer
%   6g. high-pass variance image of every swept condition
% ################################################################

% ---------------------------------------------------------------------------
% 6-PARAMS. ALL VIEWER PARAMETERS IN ONE PLACE
%     Nothing here is recomputed from the movies -- section (5) already did the
%     scoring with scr.*, which is reloaded from the sweep .mat. These only
%     control selection thresholds and how the results are drawn.
% ---------------------------------------------------------------------------
f = 20;

vw = struct();
% -- which setting counts as acceptable -------------------------------------
vw.bAPfloor      = 0.9;          % worst-branch bAP kept must exceed this. THE GATE:
                                 % a setting below it lost branch-local fast signal.
vw.slopeRange    = [0.9 1.1];    % acceptable somatic spike_slope window
vw.sortBy        = 'bAP_min';    % summary-table sort key ('bAP_min' | 'blood_removed' | ...)
% -- plots -------------------------------------------------------------------
vw.corrClim      = [-0.25 1];    % colour limits of the correlation-matrix panels
vw.ratioYlim     = [0 1.3];      % y-limits of the bAP-kept panels
vw.distSmooth    = 15;           % moving-average window (pairs) for correlation-vs-distance
vw.psdSmooth     = 21;           % moving-average window (bins) for the blood spectra
vw.diffClim      = [-0.3 0.3];   % colour limits of the correlation DIFFERENCE maps (6c)
% -- stacked movie (6e) ------------------------------------------------------
vw.mvCombos      = [];           % WHICH CONDITIONS to show, as combination indices, e.g. [1 3].
                                 % [] = every combination that ran. Panels are stacked
                                 % vertically: Original on top, then corrected for each.
vw.mvShowBlood   = false;        % append each condition's blood movie under its corrected panel
vw.mvClimPrc     = [0.5 99.5];   % percentiles of the INPUT movie setting the fixed display scale
vw.mvSepPx       = 2;            % thickness (px) of the separator drawn between stacked panels
vw.mvWriteMp4    = false;        % also write an .mp4 next to the sweep .mat
% -- high-pass variance images per condition (6g) ---------------------------
vw.hpSource      = 'cache';      % where the movies come from:
                                 %  'cache'   : re-read vw.hpFrames frames from each
                                 %              condition's saved .h5 (better statistics)
                                 %  'preview' : the bdr.nPreviewFrames already in the .mat
vw.hpFrames      = 2000;         % frames read per condition when hpSource = 'cache'
vw.hpClim        = [];           % [] = shared range from the UNCORRECTED image; or set [lo hi]
vw.hpPerImage    = false;        % true = every panel scaled to itself (structure visible,
                                 % brightness no longer comparable between panels)
vw.hpContour     = false;        % draw the bvMask outline on them
vw.hpRows        = 7;            % rows in the 6g grid. Columns follow from the panel count,
                                 % so 7 rows with <=7 panels means one full-width image per
                                 % row -- which suits the 1:3 aspect of these frames.
vw.mvFrameRate   = 20;           % playback frame rate of that .mp4

sweepFile = fullfile(save_figto, sprintf('R3_bleachdr_sweep_f%d.mat', f));
if ~exist('SW','var') || ~exist('Cmat_sw','var')
    fprintf('[R3-6] loading %s ...\n', sweepFile);
    load(sweepFile);
end
nCombo = numel(SW);
runIdx = find(okRun);                    % combinations that completed
if ~exist('cmapCorr','var'), cmapCorr = gen_colormap([0 0.5 1; 1 1 1; 1 0 0], 256); end

% ---------------------------------------------------------------------------
% 6a. SUMMARY
%     Target: HIGH blood_removed / HF_removed with spike_slope ~ 1 (spikes
%     untouched) and a HIGH blood_in/out ratio (what was removed sits on the
%     vessels). corr_vs_in tells how much the correlation structure moved.
% ---------------------------------------------------------------------------
if ~all(okRun), fprintf('[R3-6] %d/%d combination(s) failed and are dropped\n', nnz(~okRun), nCombo); end
sweepTable = struct2table(rmfield(SW(okRun), {'bAP_ratio','fast_keep','cacheFile'}));
sweepTable = sortrows(sweepTable, vw.sortBy, 'descend');   % safest branch first
disp(sweepTable(:, {'tag','bAP_median','bAP_min','bAP_worstclass','bAP_worstclass_val', ...
                    'blood_removed','hf_removed','bloodratio','spike_slope','corr_vs_in','seconds'}));

figure('Name','R3-BVsweep-metrics','NumberTitle','off'); clf
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; hold on
scatter([SW.blood_removed], [SW.spike_slope], 60, 'filled');
text([SW.blood_removed], [SW.spike_slope], strrep({SW.tag},'_','\_'), 'FontSize',7, 'VerticalAlignment','bottom');
yline(1,'k--'); xlabel('blood removed (1 - var ratio)'); ylabel('spike slope'); grid on
title('blood removal vs spike preservation');
nexttile; hold on
scatter([SW.hf_removed], [SW.bloodratio], 60, 'filled');
yline(hfratio_in,'k--','uncorrected HF in/out'); xlabel('HF power removed'); ylabel('removed power in/out bvMask'); grid on
title('is the removed power vascular?');
nexttile; hold on
scatter([SW.rem_smooth], [SW.corr_vs_in], 60, 'filled');
xline(mean(rem_in,'omitnan'),'k--','uncorrected'); yline(1,'k:');
xlabel(sprintf('%% var removed by %d-ms smoothing', scr.sub_win)); ylabel('r(Cmat_{corr}, Cmat_{uncorr})'); grid on
title('effect on the subthreshold statistics');
sgtitle(sprintf('f=%d  bleach-dr parameter sweep (frames %d:%d)', f, sweepFrames(1), sweepFrames(end)));

% blood band check: where does each swept band sit in the soma spectrum?
figure('Name','R3-BVsweep-spectra','NumberTitle','off'); clf; hold on
Lp   = 2^nextpow2(nTS);  fax = (0:Lp/2-1) * (bdr.frameRate/Lp);
psd_ = @(x) movmean(abs(fft(x(:) - mean(x(:)), Lp)).^2, vw.psdSmooth);
Pin  = psd_(soma_in);
plot(fax, Pin(1:Lp/2), 'k', 'LineWidth',1.5, 'DisplayName','uncorrected');
for k = 1:numel(blood_tr)
    if isempty(blood_tr{k}), continue; end
    Pb = psd_(blood_tr{k});
    plot(fax, Pb(1:Lp/2), '-', 'DisplayName', strrep(SW(k).tag,'_','\_'));
end
set(gca,'XScale','log','YScale','log'); xlim([1 bdr.frameRate/2]); grid on
xlabel('frequency (Hz)'); ylabel('power (a.u.)'); legend('Location','sw','Box','off');
title(sprintf('f=%d  soma-projected blood component spectrum per parameter set', f));

% ---------------------------------------------------------------------------
% 6b. DID THE BRANCH-LOCAL FAST SIGNAL SURVIVE?
%     bAP_ratio(roi) = bAP amplitude after / before. The bAP is a real, fast,
%     branch-local voltage event, so it is the positive control: a setting that
%     removes blood without touching the neuron keeps this at 1 on EVERY branch,
%     including the distal ones where the signal is smallest.
%     fast_keep(roi) = shot-noise-free fast (raw - 20 ms) power after / before.
%     This one SHOULD drop -- that drop is the blood being removed. The
%     diagnostic pattern is:
%        bAP_ratio ~ 1 everywhere  +  fast_keep < 1        -> blood only. Good.
%        bAP_ratio falling on distal/thin branches         -> signal being eaten.
%        bAP_ratio and fast_keep falling together per ROI  -> signal being eaten.
% ---------------------------------------------------------------------------
figure('Name','R3-BVsweep-branchSurvival','NumberTitle','off'); clf
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
cmapK = lines(nCombo);

ax1 = nexttile; hold(ax1,'on');           % bAP kept vs dendritic position
for k = runIdx
    plot(ax1, dendaxisS, SW(k).bAP_ratio, '.-', 'Color',cmapK(k,:), 'MarkerSize',10, ...
         'DisplayName', strrep(SW(k).tag,'_','\_'));
end
yline(ax1,1,'k--'); yline(ax1,vw.bAPfloor,'r:');
xlabel(ax1,'dendritic position (\mum)'); ylabel(ax1,'bAP amplitude kept'); grid(ax1,'on');
title(ax1,'bAP survival along the dendrite (want flat at 1)'); ylim(ax1,vw.ratioYlim);

ax2 = nexttile; hold(ax2,'on');           % bAP kept per dendrite class
cls = unique(labelvecS(~isnan(labelvecS)))';
M = nan(numel(runIdx), numel(cls));
for i = 1:numel(runIdx)
    for c = 1:numel(cls)
        M(i,c) = median(SW(runIdx(i)).bAP_ratio(labelvecS==cls(c)), 'omitnan');
    end
end
bar(ax2, M'); yline(ax2,1,'k--'); yline(ax2,vw.bAPfloor,'r:');
set(ax2,'XTick',1:numel(cls),'XTickLabel',classNm(cls)); ylabel(ax2,'median bAP kept');
grid(ax2,'on'); title(ax2,'bAP survival per dendrite class'); ylim(ax2,vw.ratioYlim);
legend(ax2, strrep({SW(runIdx).tag},'_','\_'), 'Location','southoutside','Box','off','FontSize',7);

ax3 = nexttile; hold(ax3,'on');           % the diagnostic scatter, per ROI
for k = runIdx
    scatter(ax3, SW(k).fast_keep, SW(k).bAP_ratio, 22, cmapK(k,:), 'filled', ...
            'MarkerFaceAlpha',0.6, 'DisplayName', strrep(SW(k).tag,'_','\_'));
end
yline(ax3,1,'k--'); xlabel(ax3,'fast power kept (blood should lower this)');
ylabel(ax3,'bAP kept (must stay 1)'); grid(ax3,'on'); ylim(ax3,vw.ratioYlim);
title(ax3,'per-ROI: blood removed vs signal lost');

ax4 = nexttile; hold(ax4,'on');           % the summary trade-off, per combo
for k = runIdx
    scatter(ax4, SW(k).blood_removed, SW(k).bAP_min, 70, cmapK(k,:), 'filled');
    text(ax4, SW(k).blood_removed, SW(k).bAP_min, ['  ' strrep(SW(k).tag,'_','\_')], 'FontSize',7);
end
yline(ax4,1,'k--'); yline(ax4,vw.bAPfloor,'r:','acceptable floor');
xlabel(ax4,'blood removed'); ylabel(ax4,'WORST-branch bAP kept'); grid(ax4,'on');
title(ax4,'the trade-off that decides the setting');
sgtitle(sprintf('f=%d  did the branch-local fast signal survive?', f));

fprintf('\n[R3-6] branch-local fast signal survival:\n');
for k = runIdx
    fprintf('   %-34s bAP kept med=%.3f min=%.3f (worst: %s %.3f) | fast power kept med=%.3f\n', ...
        SW(k).tag, SW(k).bAP_median, SW(k).bAP_min, SW(k).bAP_worstclass, ...
        SW(k).bAP_worstclass_val, median(SW(k).fast_keep,'omitnan'));
end

% ---------------------------------------------------------------------------
% 6c. DOES BLOOD-CORRECTION REPLACE THE 20-ms SMOOTHING?
%     The paper smooths at 20 ms to kill the blood. If bleach-dr removes the
%     blood properly, then the correlation matrix of the CORRECTED traces with
%     NO temporal smoothing at all should already equal the 20-ms-smoothed
%     matrix of the ORIGINAL traces -- i.e. the smoothing was only ever doing
%     the blood's job, and the paper's structure is not a smoothing artifact.
%
%       row 1: the matrices themselves (original 20 ms | original raw | each
%              condition's corrected raw), common leaf order
%       row 2: difference from the original 20-ms matrix
%       then : a compact heatmap of pairwise off-diagonal r between all of them
% ---------------------------------------------------------------------------
% leaf order from the original 20-ms matrix, reused everywhere
DdS = 1 - Cmat_in;  DdS = (DdS+DdS')/2;  DdS(1:nROIs+1:end) = 0;  DdS(isnan(DdS)) = 1;
dvecS = squareform(DdS);  ZcS = linkage(dvecS,'average');
try
    ordS = optimalleaforder(ZcS, dvecS);            % Stats toolbox
catch
    fT = figure('Visible','off'); [~,~,ordS] = dendrogram(ZcS, 0); close(fT);
end
tri1S = triu(true(nROIs),1);

% the set of matrices being compared: 2 references + one per condition
matsC  = [{Cmat_in}, {Cmat_raw_in}, Cmat_raw_sw(runIdx)];
namesC = [{sprintf('ORIGINAL %d ms', scr.sub_win)}, {'ORIGINAL raw'}, ...
          arrayfun(@(k) sprintf('%s\nraw', SW(k).tag), runIdx, 'UniformOutput', false)];
nM = numel(matsC);

figure('Name','R3-BVsweep-rawVsSmoothed','NumberTitle','off'); clf
tlR = tiledlayout(2, nM, 'TileSpacing','compact','Padding','compact');
for m = 1:nM
    nexttile(m);
    imagesc(matsC{m}(ordS,ordS), vw.corrClim); axis image; colormap(gca,cmapCorr);
    if m > 2
        title(sprintf('%s\nr=%.3f', strrep(namesC{m},'_','\_'), SW(runIdx(m-2)).corrRaw_vs_in20), 'FontSize',8);
    else
        title(strrep(namesC{m},'_','\_'), 'FontSize',8, 'FontWeight','bold');
    end
    if m == 1, ylabel('ROI (leaf order)'); end
    set(gca,'XTick',[],'YTick',[]);
end
for m = 1:nM
    nexttile(nM + m);
    D = matsC{m} - Cmat_in;                       % difference from the 20-ms original
    imagesc(D(ordS,ordS), vw.diffClim); axis image; colormap(gca,cmapCorr);
    title(sprintf('\\Delta vs %d ms   mean|\\Delta|=%.3f', scr.sub_win, mean(abs(D(tri1S)),'omitnan')), 'FontSize',7);
    set(gca,'XTick',[],'YTick',[]);
end
title(tlR, sprintf(['f=%d  unsmoothed corrected vs %d-ms-smoothed original ' ...
                    '(top: matrices, bottom: difference)'], f, scr.sub_win));

% compact similarity heatmap: pairwise off-diagonal r between every matrix above
Rsim = nan(nM);
for i = 1:nM
    for j = 1:nM
        a = matsC{i}(tri1S);  b = matsC{j}(tri1S);  ok = ~isnan(a) & ~isnan(b);
        Rsim(i,j) = corr(a(ok), b(ok));
    end
end
shortNm = cellfun(@(x) strrep(strrep(x, newline, ' '), '_','\_'), namesC, 'UniformOutput', false);
figure('Name','R3-BVsweep-corrSimilarity','NumberTitle','off'); clf
imagesc(Rsim, [min(Rsim(:)) 1]); axis image; colormap(gca, flipud(gray(256))); colorbar
set(gca,'XTick',1:nM,'XTickLabel',shortNm,'XTickLabelRotation',35, ...
        'YTick',1:nM,'YTickLabel',shortNm,'FontSize',8);
for i = 1:nM
    for j = 1:nM
        text(j, i, sprintf('%.3f', Rsim(i,j)), 'HorizontalAlignment','center', ...
             'FontSize',8, 'Color', [1 1 1]*(Rsim(i,j) > (min(Rsim(:))+1)/2));
    end
end
title(sprintf('f=%d  off-diagonal Pearson r between correlation matrices', f));

fprintf('\n[R3-6] unsmoothed corrected matrix vs the %d-ms original:\n', scr.sub_win);
fprintf('   %-34s %8s %8s %8s\n', 'condition', 'r vs 20ms', 'r vs raw', 'mean|d|');
fprintf('   %-34s %8.3f %8s %8.3f\n', 'ORIGINAL raw (no correction)', Rsim(2,1), '-', ...
        mean(abs(Cmat_raw_in(tri1S) - Cmat_in(tri1S)),'omitnan'));
for k = runIdx
    fprintf('   %-34s %8.3f %8.3f %8.3f\n', SW(k).tag, ...
            SW(k).corrRaw_vs_in20, SW(k).corrRaw_vs_inRaw, SW(k).corrRaw_absdiff);
end
fprintf(['   Read it as: a condition whose "r vs 20ms" beats the ORIGINAL raw row has been\n' ...
         '   pushed TOWARDS the smoothed original by the correction -- blood-removal is doing\n' ...
         '   what the smoothing did. Cross-check bAP_min in 6b before believing it.\n']);

% ---------------------------------------------------------------------------
% 6d. BEST COMBINATION: corrected correlation matrix vs the 20-ms matrix
%     "Best" = removes the most blood while keeping the WORST branch's bAP
%     above bAPfloor. Signal preservation gates the choice, not blood removal.
% ---------------------------------------------------------------------------
okSel = okRun & [SW.bAP_min] > vw.bAPfloor ...
              & [SW.spike_slope] > vw.slopeRange(1) & [SW.spike_slope] < vw.slopeRange(2);
if ~any(okSel)
    okSel = okRun;
    warning('[R3-6] NO setting kept every branch above bAP %.2f -- all of them eat branch signal', vw.bAPfloor);
end
cand  = find(okSel);
[~,b] = max([SW(cand).blood_removed]);  kBest = cand(b);
fprintf('\n[R3-6] best combination: %s (blood_rm=%.3f, worst-branch bAP=%.3f, spk_slope=%.3f)\n', ...
        SW(kBest).tag, SW(kBest).blood_removed, SW(kBest).bAP_min, SW(kBest).spike_slope);
% (ordS / tri1S come from 6c)

figure('Name','R3-BVsweep-corrMats','NumberTitle','off'); clf
nPan = 1 + numel(Cmat_sw(~cellfun(@isempty,Cmat_sw)));
tlS  = tiledlayout(1, nPan, 'TileSpacing','compact','Padding','compact');
nexttile; imagesc(Cmat_in(ordS,ordS), vw.corrClim); axis image; colormap(gca,cmapCorr);
title(sprintf('uncorrected (%d ms)', scr.sub_win)); ylabel('ROI (leaf order)');
for k = 1:nCombo
    if isempty(Cmat_sw{k}), continue; end
    nexttile; imagesc(Cmat_sw{k}(ordS,ordS), vw.corrClim); axis image; colormap(gca,cmapCorr);
    title(sprintf('%s\nr=%.3f', strrep(SW(k).tag,'_','\_'), SW(k).corr_vs_in), 'FontSize',8);
end
cb = colorbar; cb.Layout.Tile = 'east';
title(tlS, sprintf('f=%d  %d-ms subthreshold correlation: uncorrected vs bleach-dr settings', f, scr.sub_win));

figure('Name','R3-BVsweep-corrVsDist','NumberTitle','off'); clf; hold on
xdS   = Dpair(tri1S);
ycS   = Cmat_in(tri1S);  okS = ~isnan(xdS) & ~isnan(ycS);
[xs,si] = sort(xdS(okS));  ys = ycS(okS);  ys = ys(si);
plot(xs, movmean(ys,vw.distSmooth), 'k-', 'LineWidth',2, 'DisplayName','uncorrected');
for k = 1:nCombo
    if isempty(Cmat_sw{k}), continue; end
    ycS = Cmat_sw{k}(tri1S);  okS = ~isnan(xdS) & ~isnan(ycS);
    [xs,si] = sort(xdS(okS));  ys = ycS(okS);  ys = ys(si);
    plot(xs, movmean(ys,vw.distSmooth), '-', 'LineWidth',1.2, 'DisplayName', strrep(SW(k).tag,'_','\_'));
end
xlabel('dendritic distance (\mum)'); ylabel('pairwise correlation'); grid on; legend('Box','off');
title(sprintf('f=%d  correlation vs distance, blood-corrected vs uncorrected', f));

fprintf('\n[R3-6] Interpretation:\n');
fprintf('   - bAP_min ~ 1 => no branch lost its fast signal. THIS GATES EVERYTHING ELSE.\n');
fprintf('   - fast_keep < 1 while bAP_ratio = 1 => only blood was taken out of the fast band.\n');
fprintf('   - blood_in/out >> 1 and r(HF, blood) high => what was removed is spatially vascular.\n');
fprintf('   - corr_vs_in ~ 1 => the subthreshold correlation structure (the paper''s result) is\n');
fprintf('     unchanged by the blood correction.\n');

%% ---------------------------------------------------------------------------
% 6e. STACKED MOVIE, first bdr.nPreviewFrames frames
%     Panels are stacked VERTICALLY: Original on top, then the corrected movie
%     of each condition you picked in vw.mvCombos (and its blood movie too, if
%     vw.mvShowBlood). Every panel shares one fixed intensity scale taken from
%     the ORIGINAL, so "less flicker further down" reads as real removal rather
%     than as rescaling. A bright separator line marks the panel boundaries.
% ---------------------------------------------------------------------------
assert(exist('prev_in','var')==1, ...
    '[R3-6] no preview movies -- rerun section (5) with bdr.savePreview = true');

mvCombos = vw.mvCombos;
if isempty(mvCombos), mvCombos = runIdx; end
mvCombos = mvCombos(ismember(mvCombos, runIdx));      % drop failed / out-of-range picks
assert(~isempty(mvCombos), '[R3-6] vw.mvCombos selects no combination that ran');

mvClim = prctile(prev_in(:), vw.mvClimPrc);           % fixed scale from the ORIGINAL
sepRow = [];
if vw.mvSepPx > 0
    sepRow = repmat(single(mvClim(2)), vw.mvSepPx, size(prev_in,2), size(prev_in,3));
end

mont      = prev_in;
mvPanels  = {'ORIGINAL'};
for k = mvCombos(:)'
    mont = cat(1, mont, sepRow, prev_corr{k});
    mvPanels{end+1} = sprintf('corrected: %s', SW(k).tag);            %#ok<SAGROW>
    if vw.mvShowBlood
        mont = cat(1, mont, sepRow, prev_blood{k});
        mvPanels{end+1} = sprintf('blood:     %s', SW(k).tag);        %#ok<SAGROW>
    end
end

fprintf('\n[R3-6] stacked movie: %d panels x %d frames (%d x %d px)\n', ...
        numel(mvPanels), size(mont,3), size(mont,1), size(mont,2));
for i = 1:numel(mvPanels), fprintf('        panel %d (top->bottom): %s\n', i, mvPanels{i}); end

figure('Name','R3-BVsweep-movie','NumberTitle','off'); clf
moviefixsc(mont, mvClim);
title(sprintf('f=%d  top->bottom: %s   [scale %.1f %.1f]', f, ...
      strrep(strjoin(mvPanels, ' / '), '_', '\_'), mvClim(1), mvClim(2)), 'FontSize', 8);

if vw.mvWriteMp4
    mp4fn = fullfile(save_figto, sprintf('R3_BVsweep_f%d_stack%s.mp4', f, ...
                     sprintf('_%d', mvCombos)));
    writeMov(mp4fn, mont, [], vw.mvFrameRate, 1, mvClim);
    fprintf('[R3-6] wrote %s\n', mp4fn);
end

%% ---------------------------------------------------------------------------
% 6f. INTERACTIVE VIEWER  -- everything above, in one window
%     Pick conditions in the list (ctrl/shift for several) and the movie AND
%     every plot update together, so more blood removed can be weighed against
%     the branch-signal cost in one glance.
%
%       list      : the conditions that ran; multi-select stacks them
%       Show      : Corrected / Blood / Corrected-Original / Original only
%       checkboxes: add a blood panel, keep ORIGINAL on top, tile horizontally
%       scale     : shrink/expand the shared display range (all panels share it)
%       Play/fps  : animate; MP4 writes the current stack out
%       right     : primary condition's unsmoothed correlation matrix and its
%                   difference from the ORIGINAL sub_win matrix, plus bAP
%                   survival and the per-ROI trade-off for ALL selections
% ---------------------------------------------------------------------------
BVsweep_viewer(sweepFile);          % or: BVsweep_viewer(struct(...)) from the workspace

%% ---------------------------------------------------------------------------
% 6g. HIGH-PASS VARIANCE IMAGE FOR EVERY SWEPT CONDITION
%     var( mov - movmean(mov, sub_win, 3) ) per pixel, for the uncorrected movie
%     and for each parameter combination. This is the spatial version of the
%     sweep: it shows WHERE each setting removed the fast power, so a setting
%     that removes it off the vessels as well as on them is visible directly
%     rather than only through the summary numbers.
% ---------------------------------------------------------------------------
runIdxG = find(okRun);
assert(~isempty(runIdxG), '[R3-6g] no combination completed');

movG = cell(1, 1+numel(runIdxG));    % {uncorrected, cond 1, cond 2, ...}
nmG  = [{'UNCORRECTED'}, {SW(runIdxG).tag}];
switch vw.hpSource
    case 'cache'
        okCache = arrayfun(@(k) ~isempty(SW(k).cacheFile) && isfile(SW(k).cacheFile), runIdxG);
        if ~all(okCache)
            warning('[R3-6g] %d cached movie(s) missing -- falling back to the saved previews', nnz(~okCache));
            vw.hpSource = 'preview';
        end
    case 'preview'
    otherwise, error('vw.hpSource must be ''cache'' or ''preview''');
end

if strcmp(vw.hpSource,'preview')
    assert(exist('prev_in','var')==1, '[R3-6g] no previews in the sweep file');
    movG{1} = prev_in;
    for i = 1:numel(runIdxG), movG{i+1} = prev_corr{runIdxG(i)}; end
    srcTagG = sprintf('%d preview frames', size(prev_in,3));
else
    % the uncorrected movie is not cached, so it still comes from the preview
    assert(exist('prev_in','var')==1, '[R3-6g] no uncorrected preview in the sweep file');
    movG{1} = prev_in;
    [nrG, ncG] = size(prev_in, [1 2]);
    for i = 1:numel(runIdxG)
        k  = runIdxG(i);
        nF = min(vw.hpFrames, h5info(SW(k).cacheFile, '/corrected').Dataspace.Size(3));
        movG{i+1} = permute(h5read(SW(k).cacheFile, '/corrected', [1 1 1], [ncG nrG nF]), [2 1 3]);
    end
    srcTagG = sprintf('%d frames from the cached .h5 (uncorrected: %d preview frames)', ...
                      min(vw.hpFrames, size(movG{2},3)), size(prev_in,3));
end

hpG = cell(1, numel(movG));  ratG = nan(1, numel(movG));
for i = 1:numel(movG)
    hp    = movG{i} - movmean(movG{i}, scr.sub_win, 3);
    hpG{i} = var(hp, 0, 3);
    ratG(i) = mean(hpG{i}(bvpixS)) / mean(hpG{i}(~bvpixS));
    movG{i} = [];                                   % free as we go
end
clear movG hp

if isempty(vw.hpClim), climG = prctile(hpG{1}(:), [1 99.5]); else, climG = vw.hpClim; end
nPanG = numel(hpG);
nRowG = max(1, vw.hpRows);
nColG = ceil(nPanG / nRowG);
get_figure('R3-6g-hpVariance');
tlG = tiledlayout(nRowG, nColG, 'TileSpacing','compact','Padding','compact');
if nPanG < nRowG*nColG
    fprintf('[R3-6g] %d panels in a %d x %d grid (%d tiles left empty)\n', ...
            nPanG, nRowG, nColG, nRowG*nColG - nPanG);
end
for i = 1:nPanG
    ax = nexttile(i);
    if vw.hpPerImage, cl = prctile(hpG{i}(:), [1 99.5]); else, cl = climG; end
    if cl(2) <= cl(1), cl(2) = cl(1) + eps; end
    imagesc(ax, hpG{i}, cl); axis(ax,'image','off'); colormap(ax, turbo); hold(ax,'on');
    if vw.hpContour, contour(ax, bvBWS, [0.5 0.5], 'w', 'LineWidth', 0.8); end
    if i == 1, ttl = nmG{i}; else, ttl = strrep(nmG{i},'_','\_'); end
    title(ax, sprintf('%s\nin/out bvMask = %.2f', ttl, ratG(i)), 'FontSize', 8);
    if vw.hpPerImage || i == nPanG, colorbar(ax); end
end
if vw.hpPerImage, scaleTagG = 'each panel scaled to itself';
else,             scaleTagG = sprintf('shared scale [%.4g %.4g] from the uncorrected image', climG(1), climG(2)); end
title(tlG, sprintf('f=%d  high-pass variance (raw - %d ms) per condition -- %s; %s', ...
      f, scr.sub_win, srcTagG, scaleTagG), 'Interpreter','none');

fprintf('\n[R3-6g] high-pass variance in/out bvMask, by condition (%s):\n', srcTagG);
for i = 1:nPanG
    fprintf('   %-34s %.2f', nmG{i}, ratG(i));
    if i > 1, fprintf('   (%.0f%% of the uncorrected mean HF power remains)', ...
                      100*mean(hpG{i}(:))/mean(hpG{1}(:))); end
    fprintf('\n');
end
fprintf(['   A GOOD setting lowers the mean HF power while RAISING nothing off the\n' ...
         '   vessels: the in/out ratio should fall because the numerator fell, not\n' ...
         '   because the denominator rose. Cross-check bAP_min in 6b.\n']);

%% ################################################################
%  (7) FULL 10,000-FRAME BLOOD CORRECTION  (f = 20, frames 78000:88000)
%
%  Section (5)/(6) picked the parameters on a short window. This runs the SAME
%  settings on the whole 10k segment and produces the deliverables:
%
%     * the corrected movie          -> <cacheDir>/f20_<tag>_seg<digest>.h5
%     * corrected traces             -> BVcorr10k_traces.mat
%       (combined + odd/even checkerboard channels, so downstream correlation
%        matrices keep a shot-noise-free diagonal, exactly as elsewhere)
%     * the verification metrics recomputed on all 10k frames
%
%  WHY THIS MATTERS BEYOND JUST BEING LONGER:
%  the 1100-frame window held only a handful of somatic spikes, so bAP_ratio
%  there was noise-dominated (that is why part of it came out above 1). With
%  ~10x the spikes the bAP survival number becomes trustworthy -- so this run
%  is also the real test of whether these settings preserve the branch-local
%  fast signal. 7e prints the verdict.
% ################################################################
f  = 20;
Pr = PCresult{f};
if ~exist('Result','var') || ~isfield(Result,'spike')
    load(fullfile(fpath{f},'PC_Result.mat'),'Result');
end

% ---------------------------------------------------------------------------
% 7a. PARAMETERS  (the setting chosen in section 6, now fixed)
% ---------------------------------------------------------------------------
b7 = struct();
b7.venvpy   = '/Users/bhlee1117/Documents/GitHub/bleach-dr/.venv/bin/python';
b7.script   = '/Users/bhlee1117/Documents/GitHub/bleach-dr/run_correction.py';

% -- the chosen correction parameters ---------------------------------------
b7.bandpass       = [200 220];    % --bandpass-low/high
b7.strategyPreset = 2;            % --strategy-preset   -> blood blocks {20,40,60}
b7.blockSize      = 20;           % --block-size
b7.maxComponents  = 30;           % --max-components
b7.rank1Method    = 'skip';       % --rank1-method
b7.frameRate      = 1000;         % --frame-rate
b7.device         = 'cpu';

% -- the segment -------------------------------------------------------------
b7.frames    = 78000:88000;       % 10,001 frames
b7.segSource = 'mov_res';         % 'mov_res' : reuse section (4)'s movie -- IDENTICAL
                                  %             preprocessing to what section (5) validated
                                  % 'chunk'   : reload via loadVoltChunk (adds the photobleach
                                  %             normalisation that mov_res skips; matches the
                                  %             whole-movie production run instead)
b7.mc_prefix = 'mc_ShutterReg';
b7.seg_len   = 15000;
b7.rank1Source = 'STAmap';        % same rank1 image as the sweep

% -- memory / IO -------------------------------------------------------------
% A single pass over 10k frames makes python hold the movie (~4 GB float32) plus
% the blood basis (v and a are rank x nFrames; preset 2 at block 20 reaches
% ~37k rows -> several GB). Chunking bounds that. Each processed window is
% chunkSize + 2*overlap frames and must stay >= 1024; the overlap is trimmed
% after correction so chunk boundaries do not show.
b7.chunkSize    = 5000;           % 0 = single pass (faithful but memory-hungry)
b7.chunkOverlap = 500;
b7.cacheDir  = fullfile(fpath{f}, 'BVsweep_corrected');   % same store as the sweep
b7.useCache  = true;
% fpath{f} is on the cohen_lab NFS mount, so the ~3.7 GB input h5 goes to a LOCAL
% disk instead (it is deleted as soon as python finishes). Only the corrected
% result is written to the network store.
b7.scratch   = '/Volumes/BHL18TB_D1';
if ~isfolder(b7.scratch), b7.scratch = tempdir; end
b7.traceBlk  = 2000;              % frames per block when streaming traces back out
b7.qcFrames  = 1:3000;            % sub-window used for the movie-LEVEL variance maps
b7.nPreview  = 200;               % frames kept for a before/after movie

% scoring uses the SAME scr.* as section (5) so the numbers are comparable
if ~exist('scr','var')
    scr = struct('sub_win',20,'sub_dilate',7,'minPairs',100,'bAP_pre',5, ...
                 'bAP_post',15,'bAP_baseprc',30,'bvThresh',0.2,'blueDilate',201);
end

tag7 = sprintf('bp%g-%g_p%d_%s_b%d_c%d', b7.bandpass(1), b7.bandpass(2), ...
               b7.strategyPreset, b7.rank1Method, b7.blockSize, b7.maxComponents);
assert((b7.strategyPreset+1)*b7.blockSize <= min(size(Pr.ref_im)), ...
       '[R3-7] blocks up to %d px exceed the %d px FOV', ...
       (b7.strategyPreset+1)*b7.blockSize, min(size(Pr.ref_im)));

% ---------------------------------------------------------------------------
% 7b. INPUT MOVIE + rank1 IMAGE
% ---------------------------------------------------------------------------
switch b7.segSource
    case 'mov_res'
        assert(exist('mov_res','var')==1, '[R3-7] run section (4) first (mov_res / segFrames)');
        [tf, loc] = ismember(b7.frames, segFrames);
        assert(all(tf), '[R3-7] b7.frames must lie inside segFrames (%d:%d)', ...
               segFrames(1), segFrames(end));
        mov7 = single(-mov_res(:,:,loc));                     % spike-positive
    case 'chunk'
        load(fullfile(fpath{f},'output_data.mat'),'Device_Data');
        sz7  = double(Device_Data{1,3}.ROI([2 4]));
        j7   = ceil(b7.frames(1)/b7.seg_len);
        assert(ceil(b7.frames(end)/b7.seg_len) == j7, ...
               '[R3-7] the window must stay inside one %d-frame mc chunk', b7.seg_len);
        loc7 = b7.frames - (j7-1)*b7.seg_len;
        fprintf('[R3-7] loading mc chunk %d ...\n', j7);
        mc7  = loadVoltChunk(fpath{f}, j7, sz7, b7.seg_len, b7.mc_prefix);
        mov7 = single(mc7(:,:,loc7));
        clear mc7
    otherwise
        error('b7.segSource must be ''mov_res'' or ''chunk''');
end
[nR7, nC7, nT7] = size(mov7);
fprintf('[R3-7] segment %d x %d x %d  (frames %d:%d, %s)\n', ...
        nR7, nC7, nT7, b7.frames(1), b7.frames(end), b7.segSource);

% rank1 image: SS spike-triggered-average movie, peak over time, spike-positive
STAmovie7 = importdata(fullfile(fpath{f},'SS_STAmovie.mat'));
STAmap7   = max(-STAmovie7, [], 3);
r07 = floor((nR7 - size(STAmap7,1))/2);  c07 = floor((nC7 - size(STAmap7,2))/2);
spikeImg7 = zeros(nR7, nC7, 'single');
spikeImg7(r07 + (1:size(STAmap7,1)), c07 + (1:size(STAmap7,2))) = STAmap7;

% footprints (sorted ROI order) + odd/even checkerboard channels
oddPix7  = mod((1:nR7)' + (1:nC7), 2) == 1;
ft7      = Pr.Ftprnts(:,:,Pr.sortDist);
ftv7     = tovec(ft7);   ftO7 = tovec(ft7 .* oddPix7);   ftE7 = tovec(ft7 .* ~oddPix7);
ftS7     = tovec(Pr.Ftprnts(:,:,1));            % soma
nROI7    = size(ftv7, 2);

spk7     = Result.spike(1, b7.frames) > 0;
good7    = ~imdilate(Pr.BlueStim(b7.frames) > 0, ones(1,scr.blueDilate));
bvMap7   = max(Pr.bvMask, [], 3);
bvpix7   = bvMap7(:) > 0;
fprintf('[R3-7] %d good frames, %d somatic spikes  (sweep window had far fewer)\n', ...
        nnz(good7), nnz(spk7));

% ---------------------------------------------------------------------------
% 7c. RUN (or reuse a cached result with matching provenance)
% ---------------------------------------------------------------------------
if ~isfolder(b7.cacheDir), mkdir(b7.cacheDir); end
step7 = max(1, floor(numel(mov7)/1e5));
md7   = java.security.MessageDigest.getInstance('MD5');
md7.update(typecast(single([reshape(mov7(1:step7:end),[],1); spikeImg7(:)]), 'int8'));
sig7  = sprintf('%02x', mod(double(md7.digest()), 256));  sig7 = sig7(1:10);
in7   = fullfile(b7.scratch, 'bv_full_in.h5');
out7  = fullfile(b7.cacheDir, sprintf('f%d_%s_seg%s.h5', f, tag7, sig7));
fprintf('[R3-7] output: %s\n', out7);

% ---- disk pre-flight: input (deleted later) + corrected + blood ------------
needIn7  = nT7*nR7*nC7*4  / 2^30;
needOut7 = 2*needIn7;
for chk = {b7.scratch needIn7; b7.cacheDir needOut7}'
    dd = chk{1};  while ~isfolder(dd) && ~isempty(dd), dd = fileparts(dd); end
    fo = java.io.File(dd);
    freeGB = double(fo.getUsableSpace())/2^30;
    if freeGB < chk{2}*1.2
        warning('[R3-7] only %.1f GB free on %s but %.1f GB is needed', freeGB, dd, chk{2});
    else
        fprintf('[R3-7] %-40s %.0f GB free (need %.1f GB)\n', dd, freeGB, chk{2});
    end
end

use7 = b7.useCache && isfile(out7);
if use7
    try
        a7 = @(n) h5readatt(out7,'/',n);
        use7 = a7('bandpass_low')==b7.bandpass(1) && a7('bandpass_high')==b7.bandpass(2) && ...
               a7('block_size')==b7.blockSize && a7('max_components')==b7.maxComponents && ...
               a7('strategy_preset')==b7.strategyPreset && a7('frame_rate')==b7.frameRate && ...
               strcmp(strtrim(char(a7('rank1_method'))), b7.rank1Method) && ...
               any(strcmp({h5info(out7).Datasets.Name},'blood'));
        if ~use7, fprintf('   cached file parameters differ -> recomputing\n'); end
    catch
        use7 = false;
    end
end

if use7
    d7 = dir(out7);
    fprintf('[R3-7] REUSING cached result (%.1f GB, written %s)\n', d7.bytes/2^30, ...
            string(datetime(d7.datenum,'ConvertFrom','datenum','Format','yyyy-MM-dd HH:mm')));
else
    if isfile(in7),  delete(in7);  end
    if isfile(out7), delete(out7); end
    fprintf('[R3-7] writing input h5 (%.1f GB) ...\n', nT7*nR7*nC7*4/2^30);
    mo7 = permute(mov7, [2 1 3]);   si7 = permute(spikeImg7, [2 1]);
    h5create(in7,'/mov',        size(mo7),'Datatype','single'); h5write(in7,'/mov',        mo7);
    h5create(in7,'/rank1_image',size(si7),'Datatype','single'); h5write(in7,'/rank1_image',si7);
    clear mo7 si7

    cmd7 = sprintf(['"%s" "%s" --input "%s" --output "%s" ' ...
        '--frame-rate %g --bandpass-low %g --bandpass-high %g ' ...
        '--block-size %d --max-components %d --strategy-preset %d ' ...
        '--rank1-method %s --device %s --chunk-size %d --chunk-overlap %d'], ...
        b7.venvpy, b7.script, in7, out7, b7.frameRate, b7.bandpass(1), b7.bandpass(2), ...
        b7.blockSize, b7.maxComponents, b7.strategyPreset, b7.rank1Method, b7.device, ...
        b7.chunkSize, b7.chunkOverlap);
    fprintf('[R3-7] %s\n', cmd7);
    t7 = tic;
    [st7, o7] = system(cmd7, '-echo');
    if st7 ~= 0
        if isfile(out7), delete(out7); end          % never keep a partial result
        error('[R3-7] python failed:\n%s', o7);
    end
    fprintf('[R3-7] correction done in %.1f min\n', toc(t7)/60);
    if isfile(in7), delete(in7); end
end

% ---------------------------------------------------------------------------
% 7d. TRACES  -- streamed block by block so the 10k movie is never fully in RAM
% ---------------------------------------------------------------------------
% Corrected, blood AND uncorrected are projected in the SAME blocked pass, so
% no full 10k movie (4 GB each) is ever materialised alongside its transpose.
tr7C = zeros(nROI7, nT7);  tr7_1 = zeros(nROI7, nT7);  tr7_2 = zeros(nROI7, nT7);
bl7C = zeros(nROI7, nT7);  in7C  = zeros(nROI7, nT7);
in7_1 = zeros(nROI7, nT7); in7_2 = zeros(nROI7, nT7);
tr7_soma = zeros(nT7,1);   bl7_soma = zeros(nT7,1);   in7_soma = zeros(nT7,1);
fprintf('[R3-7] extracting traces in blocks of %d frames ...\n', b7.traceBlk);
for t0 = 1:b7.traceBlk:nT7
    t1 = min(t0 + b7.traceBlk - 1, nT7);
    n  = t1 - t0 + 1;
    cc = permute(h5read(out7,'/corrected',[1 1 t0],[nC7 nR7 n]), [2 1 3]);
    bb = permute(h5read(out7,'/blood',    [1 1 t0],[nC7 nR7 n]), [2 1 3]);
    vc = tovec(cc)';   vb = tovec(bb)';   v0 = tovec(mov7(:,:,t0:t1))';
    tr7C(:,t0:t1)  = (vc * ftv7)';    tr7_1(:,t0:t1) = (vc * ftO7)';
    tr7_2(:,t0:t1) = (vc * ftE7)';    tr7_soma(t0:t1) = vc * ftS7;
    bl7C(:,t0:t1)  = (vb * ftv7)';    bl7_soma(t0:t1) = vb * ftS7;
    in7C(:,t0:t1)  = (v0 * ftv7)';    in7_1(:,t0:t1)  = (v0 * ftO7)';
    in7_2(:,t0:t1) = (v0 * ftE7)';    in7_soma(t0:t1) = v0 * ftS7;
    clear cc bb vc vb v0
    fprintf('   %d/%d frames\n', t1, nT7);
end

traceFile7 = fullfile(fpath{f}, 'BVcorr10k_traces.mat');
save(traceFile7, 'tr7C','tr7_1','tr7_2','bl7C','in7C','in7_1','in7_2', ...
     'tr7_soma','in7_soma','bl7_soma','b7','tag7','sig7','scr', ...
     'spk7','good7','nROI7','nT7','out7','-v7.3');
fprintf('[R3-7] traces -> %s\n', traceFile7);

% ---------------------------------------------------------------------------
% 7e. VERIFICATION on all 10k frames -- the properly-powered bAP check
% ---------------------------------------------------------------------------
% bAP survival (this is the number the 1100-frame window could not resolve)
bAPi7 = get_STA(in7C, spk7, scr.bAP_pre, scr.bAP_post);
bAPi7 = bAPi7 - prctile(bAPi7, scr.bAP_baseprc, 2);   ampI7 = max(bAPi7, [], 2);
bAPc7 = get_STA(tr7C, spk7, scr.bAP_pre, scr.bAP_post);
bAPc7 = bAPc7 - prctile(bAPc7, scr.bAP_baseprc, 2);   ampC7 = max(bAPc7, [], 2);
bAPkeep7 = ampC7 ./ ampI7;

% shot-noise-free fast power (odd x even covariance of the raw - sub_win part)
% NOTE the order: the moving average is taken over the FULL trace and only then
% restricted to good frames -- identical to section (5), so the numbers compare.
fi1 = in7_1 - movmean(in7_1, scr.sub_win, 2);   fi2 = in7_2 - movmean(in7_2, scr.sub_win, 2);
fc1 = tr7_1 - movmean(tr7_1, scr.sub_win, 2);   fc2 = tr7_2 - movmean(tr7_2, scr.sub_win, 2);
fastI7 = mean(fi1(:,good7).*fi2(:,good7), 2, 'omitnan');
fastC7 = mean(fc1(:,good7).*fc2(:,good7), 2, 'omitnan');
fastKeep7 = fastC7 ./ fastI7;

% correlation matrices: unsmoothed and sub_win-smoothed, before and after
sub = @(x,w) get_subthreshold(x, spk7, scr.sub_dilate, w);
Ci_raw = get_corrTensor(sub(in7_1,1), sub(in7_2,1), sub(in7C,1), good7, 0, scr.minPairs);
Ci_sm  = get_corrTensor(sub(in7_1,scr.sub_win), sub(in7_2,scr.sub_win), sub(in7C,scr.sub_win), good7, 0, scr.minPairs);
Cc_raw = get_corrTensor(sub(tr7_1,1), sub(tr7_2,1), sub(tr7C,1), good7, 0, scr.minPairs);
Cc_sm  = get_corrTensor(sub(tr7_1,scr.sub_win), sub(tr7_2,scr.sub_win), sub(tr7C,scr.sub_win), good7, 0, scr.minPairs);
C7 = struct('in_raw',Ci_raw(:,:,1),'in_sm',Ci_sm(:,:,1),'c_raw',Cc_raw(:,:,1),'c_sm',Cc_sm(:,:,1));

% ---- ROI-ORDER CHECK -------------------------------------------------------
% All six trace arrays are filled in one loop from the same ftv7/ftO7/ftE7, so
% row i is the same ROI everywhere and the four C7 matrices share one ordering.
% That is asserted rather than assumed: with matching rows, corr(in7C(i,:),
% tr7C(i,:)) on the diagonal must stand far above the off-diagonal pairs, and
% the off-diagonal structure of the two correlation matrices must agree. A row
% permutation collapses both (measured on synthetic data: diagonal separation
% 5.0 SD -> 0.0 SD, and r between the matrices 0.89 -> -0.06).
Rid  = corr(in7C(:,good7)', tr7C(:,good7)');
dId  = diag(Rid);  oId = Rid(~eye(nROI7));
sepId = (median(dId) - median(oId)) / std(oId);
tri7  = triu(true(nROI7),1);
rStruct = corr(C7.in_sm(tri7), C7.c_sm(tri7), 'rows','complete');
fprintf(['[R3-7] ROI-order check: same-ROI r = %.3f vs other-ROI %.3f (%.1f SD apart); ' ...
         'r(C7.in_sm, C7.c_sm) off-diagonal = %.3f\n'], ...
        median(dId), median(oId), sepId, rStruct);
if sepId < 2 || rStruct < 0.2
    warning(['[R3-7] the uncorrected and corrected traces do NOT look row-aligned. ' ...
             'Every C7 matrix would then be in a different ROI order -- stop and check ' ...
             'ftv7 / Pr.sortDist before using these matrices.']);
end

tri7 = triu(true(nROI7),1);
rOf  = @(A,B) corr(A(tri7 & ~isnan(A) & ~isnan(B)), B(tri7 & ~isnan(A) & ~isnan(B)));
r_rawVsSm   = rOf(C7.c_raw, C7.in_sm);      % the key one (see 6c)
r_inRawVsSm = rOf(C7.in_raw, C7.in_sm);     % baseline: no correction
r_smVsSm    = rOf(C7.c_sm,  C7.in_sm);

% movie-level maps on the QC sub-window only (bounded memory)
q7 = b7.qcFrames(b7.qcFrames <= nT7);
cq = permute(h5read(out7,'/corrected',[1 1 q7(1)],[nC7 nR7 numel(q7)]), [2 1 3]);
bq = permute(h5read(out7,'/blood',    [1 1 q7(1)],[nC7 nR7 numel(q7)]), [2 1 3]);
mq = mov7(:,:,q7);
% Spikes are fast and sit ON THE NEURON, so leaving them in makes the "high-pass"
% image a spike map rather than a blood map. Mask them out (and the blue-on
% frames) exactly as the trace statistics do: moving average over the CONTIGUOUS
% movie first, analysed frames selected afterwards.
gq7   = good7(q7);
hpT   = mq - movmean(mq, scr.sub_win, 3);   hf_i7 = var(hpT(:,:,gq7), 0, 3);
hpT   = cq - movmean(cq, scr.sub_win, 3);   hf_c7 = var(hpT(:,:,gq7), 0, 3);
clear hpT
bvar7 = var(bq(:,:,gq7) - mean(bq(:,:,gq7),3), 0, 3);
fprintf('[R3-7] high-pass images from %d/%d frames (spike + blue-on excluded)\n', nnz(gq7), numel(q7));
ratio_i7 = mean(hf_i7(bvpix7))/mean(hf_i7(~bvpix7));
ratio_c7 = mean(hf_c7(bvpix7))/mean(hf_c7(~bvpix7));
ratio_b7 = mean(bvar7(bvpix7))/mean(bvar7(~bvpix7));
prev7_in = mq(:,:,1:min(b7.nPreview,numel(q7)));
prev7_c  = cq(:,:,1:min(b7.nPreview,numel(q7)));
prev7_b  = bq(:,:,1:min(b7.nPreview,numel(q7)));
clear cq bq mq

save(traceFile7, 'bAPkeep7','ampI7','ampC7','fastKeep7','C7','hf_i7','hf_c7','bvar7', ...
     'prev7_in','prev7_c','prev7_b','r_rawVsSm','r_inRawVsSm','r_smVsSm','-append');

% ---- verdict ----
labs7 = Pr.labelvec(:);   clsNm7 = {'Basal','Soma','Trunk','Oblique','Distal'};
fprintf('\n[R3-7] ===== VERIFICATION ON %d FRAMES (%d spikes) =====\n', nT7, nnz(spk7));
fprintf('   bAP kept : median %.3f   min %.3f   (n=%d ROIs)\n', ...
        median(bAPkeep7,'omitnan'), min(bAPkeep7), nROI7);
for c = unique(labs7(~isnan(labs7)))'
    fprintf('      %-8s median %.3f   min %.3f   (%d ROIs)\n', clsNm7{c}, ...
            median(bAPkeep7(labs7==c),'omitnan'), min(bAPkeep7(labs7==c)), nnz(labs7==c));
end
fprintf('   fast power kept : median %.3f\n', median(fastKeep7,'omitnan'));
fprintf('   correlation( fast power kept , bAP kept ) = %.3f\n', ...
        corr(fastKeep7(~isnan(fastKeep7)&~isnan(bAPkeep7)), bAPkeep7(~isnan(fastKeep7)&~isnan(bAPkeep7))));
fprintf('      ~0 => only blood was removed.  strongly positive => branch signal was eaten.\n');
fprintf('   HF var in/out bvMask : %.2f -> %.2f   removed power in/out : %.2f\n', ...
        ratio_i7, ratio_c7, ratio_b7);
fprintf('   corr matrices: r(corrected raw, ORIGINAL %dms) = %.3f   [no-correction baseline %.3f]\n', ...
        scr.sub_win, r_rawVsSm, r_inRawVsSm);
fprintf('                  r(corrected %dms, ORIGINAL %dms) = %.3f\n', scr.sub_win, scr.sub_win, r_smVsSm);
if min(bAPkeep7) > 0.9
    fprintf('   VERDICT: every branch kept >90%% of its bAP -- settings preserve branch-local signal.\n');
else
    fprintf(['   VERDICT: %d/%d ROIs fell below 0.9 bAP kept (worst %.3f). With %d spikes this is\n' ...
             '            now a well-powered estimate, so it reflects real signal loss.\n'], ...
             nnz(bAPkeep7 < 0.9), nROI7, min(bAPkeep7), nnz(spk7));
end

% ---- figures ----
if ~exist('cmapCorr','var'), cmapCorr = gen_colormap([0 0.5 1; 1 1 1; 1 0 0], 256); end
figure('Name','R3-7-verify','NumberTitle','off'); clf
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nexttile; hold on
[xs,si] = sort(Pr.dendaxis(:));  plot(xs, bAPkeep7(si), '.-','MarkerSize',12);
yline(1,'k--'); yline(0.9,'r:'); ylim([0 1.4]); grid on
xlabel('dendritic position (\mum)'); ylabel('bAP kept'); title('bAP survival (10k frames)');
nexttile; hold on
scatter(fastKeep7, bAPkeep7, 26, 'filled','MarkerFaceAlpha',.6);
yline(1,'k--'); ylim([0 1.4]); grid on
xlabel('fast power kept'); ylabel('bAP kept'); title('per ROI: blood removed vs signal lost');
nexttile; imagesc(hf_i7); axis image off; colorbar; hold on
contour(bvMap7 > scr.bvThresh*max(bvMap7(:)), [.5 .5], 'g','LineWidth',1);
title(sprintf('HF var BEFORE (in/out=%.2f)', ratio_i7));
nexttile; imagesc(hf_c7); axis image off; colorbar; hold on
contour(bvMap7 > scr.bvThresh*max(bvMap7(:)), [.5 .5], 'g','LineWidth',1);
title(sprintf('HF var AFTER (in/out=%.2f)', ratio_c7));
nexttile; imagesc(C7.in_sm, [-0.25 1]); axis image; colormap(gca,cmapCorr);
set(gca,'XTick',[],'YTick',[]); title(sprintf('ORIGINAL %d ms', scr.sub_win));
nexttile; imagesc(C7.c_raw, [-0.25 1]); axis image; colormap(gca,cmapCorr);
set(gca,'XTick',[],'YTick',[]);
title(sprintf('CORRECTED unsmoothed\nr vs %dms = %.3f', scr.sub_win, r_rawVsSm));
sgtitle(sprintf('f=%d  %s  frames %d:%d', f, tag7, b7.frames(1), b7.frames(end)));

figure('Name','R3-7-movie','NumberTitle','off'); clf
sep7 = repmat(single(prctile(prev7_in(:),99.5)), 2, nC7, size(prev7_in,3));
moviefixsc(cat(1, prev7_in, sep7, prev7_c, sep7, prev7_b), ...
           prctile(prev7_in(:), [0.5 99.5]));
title(sprintf('f=%d  top->bottom: ORIGINAL / corrected / blood   (%s)', f, tag7), ...
      'Interpreter','none');

%% ################################################################
%  (8) FINAL FIGURE  (f = 20, frames b7.frames)
%
%  The reviewer's argument in one figure. Sub-ROIs (~30 px, a few um) are tiled
%  along the SWC skeleton, so each one sees a short stretch of dendrite -- some
%  sit on a blood vessel, most do not.
%
%   8a  two example traces, ON- and OFF-vessel: raw + 20 ms average, and
%       (raw - 20 ms) underneath                          [panels 1-2]
%   8b  violin: power of the raw trace vs the high-pass (raw - 20 ms), with the
%       on-vessel sub-ROIs marked                         [panel 3]
%   8c  high-pass variance image + bvMask                 [panel 4]
%   8d  the same after bleach-dr correction + bvMask      [panel 5]
%   8e  correlation matrix of the corrected movie         [panel 6]
%   8f  PCA explained variance inside the SWC dendrite mask, for the raw /
%       20-ms-smoothed / bleach-dr-corrected movie        [panel 7]
%
%  Requires section (4) (mov_res, segFrames) and section (7) (out7, b7).
% ################################################################
f  = 20;
Pr = PCresult{f};
if ~exist('Result','var') || ~isfield(Result,'spike')
    load(fullfile(fpath{f},'PC_Result.mat'),'Result');
end
assert(exist('mov_res','var')==1 && exist('segFrames','var')==1, '[R3-8] run section (4) first');
assert(exist('out7','var')==1 && isfile(out7), '[R3-8] run section (7) first');
% section 8 reads these from sections (5)/(7); name them explicitly so a missing
% one fails here with a clear message rather than deep inside a plotting call
for req8 = {'b7','C7','in7C','bl7C','nT7','r_rawVsSm','r_smVsSm','r_inRawVsSm'}
    assert(exist(req8{1},'var')==1, '[R3-8] %s is missing -- run section (7) first', req8{1});
end
for req8 = {'frames','traceBlk'}
    assert(isfield(b7, req8{1}), '[R3-8] b7.%s is missing -- rerun section (7)', req8{1});
end
if ~exist('scr','var') || ~isfield(scr,'sub_win')      % section (5) normally sets this
    scr = struct('sub_win',20,'sub_dilate',7,'minPairs',100,'bAP_pre',5, ...
                 'bAP_post',15,'bAP_baseprc',30,'bvThresh',0.2,'blueDilate',201);
    warning('[R3-8] scr was missing -- using the section (5) defaults (sub_win = %d ms)', scr.sub_win);
end

% ---------------------------------------------------------------------------
% 8-PARAMS
% ---------------------------------------------------------------------------
p8 = struct();
p8.Wsm        = 20;      % ms. "smoothed" window; high-pass := raw - movmean(raw,Wsm)
p8.subPx      = 30;      % pixels per sub-ROI (divide_ROIs target size)
p8.swcRadius  = 3;       % px radius around each SWC node -> dendrite mask
p8.Wsp        = 25;      % frames dilated around each spike -> excluded
p8.bvFracOn   = 0.30;    % an ROI counts as ON-vessel if this fraction of its pixels are
p8.bvFracOff  = 0.00;    % ... and OFF-vessel if at most this fraction are
p8.subCorrected = false; % add a third box "high-pass after bleach-dr" to the variance panel.
                         % FALSE also SKIPS projecting the corrected 10k movie onto the
                         % sub-ROIs -- that is a ~4 GB read from the h5, so leaving this
                         % off makes section (8) noticeably faster.
p8.bloodBand  = 'total'; % what the blood axis of the blood-vs-high-pass panel shows:
                         %  'total'    : variance of the whole blood component
                         %  'highpass' : variance of its (raw - Wsm) part only, i.e. the
                         %               band the y axis lives in -- the stricter test
p8.bvFracMark = 0;       % scatter/map marking threshold: an ROI is drawn as "vessel" when
                         % its overlap EXCEEDS this. 0 = mark ANY overlap at all. Whole
                         % Result.ftprnt ROIs are big, so they rarely reach bvFracOn=0.30 --
                         % that threshold is meant for the ~30 px sub-ROIs.
p8.exROI      = [21 6];      % [onIdx offIdx] to override the automatic example pick
p8.baseWin    = 501;    % ms, display-only running-median baseline
p8.showMs     = 2000;    % length of the displayed stretch
p8.showWindow = 'spiking';  % which stretch the example traces show:
                         %  'spiking' : longest blue-off run, then the sub-window with the
                         %              MOST spikes -- so the reader sees the cell firing
                         %  'silent'  : longest run with no spikes at all (the old default;
                         %              it made the traces look like they had no action
                         %              potentials, which invites the wrong question)
p8.showClip   = true;    % set the y-limits from the SUBTHRESHOLD trace, so bAPs run off the
                         % top rather than compressing the subthreshold comparison that the
                         % panel exists to make. false = fit everything, spikes included.
p8.pcaCov      = 'lag1';          % covariance the movie PCA diagonalises:
                                  %  'lag1' : X(t)X(t+1), symmetrised. White shot noise has
                                  %           zero lag-1 covariance, so the spectrum is
                                  %           noise-free -- no flat floor inflating the tail.
                                  %           Eigenvalues may go NEGATIVE once the signal is
                                  %           exhausted; that is the unbiased estimator
                                  %           fluctuating about zero, and the code normalises
                                  %           and plots only the positive ones.
                                  %  'same' : ordinary X(t)X(t) covariance (biased by noise)
p8.pcShow      = {1:5, 100:105};  % individual PC spatial maps to draw (figure R3-8-PCimages)
p8.pcShowMovie = [2 3];           % WHICH MOVIES those maps come from (vector):
                                  % 1 = raw, 2 = 20 ms smoothed, 3 = bleach-dr corrected.
                                  % One row per (movie, PC group) pair.
p8.pcBands    = {1:5, 6:10, 11:20, 21:100};   % PC bands for the per-band variance images
                          % in the R3-8-PCA figure. Contiguous and NON-overlapping; each
                          % band is capped at the available rank. Change freely, e.g.
                          % {1:3, 4:10, 11:50} -- the figure sizes itself to numel(p8.pcBands).
p8.imgFrames  = 1:3000;  % frames used for the high-pass variance IMAGES.
                          % Set to b7.qcFrames to reproduce R3-7-verify exactly.
p8.imgUseGood = true;    % exclude spike + blue-on frames from those images. Keep TRUE:
                          % spikes are fast and live on the neuron, so including them
                          % turns the "high-pass" map into a spike map (measured: it
                          % inflates neuron pixels ~16x and flips the vessel/neuron
                          % contrast from 8.2 to 0.5). Section (7) now does the same.
p8.nPcaFrames = 5000;    % frames used for the PCA. The component count is
                         % min(#dendrite pixels, #analysed frames) -- with ~21k pixels the
                         % FRAMES set the rank, so this is what fixes how many PCs exist.
                         % Measured: svd(20959 x 3800) = 8.8 s, ~26 s for the three movies.
                         % Cost grows with frames^2, so raising this gets expensive fast.
p8.exColr     = [0.85 0.33 0.10; 0.00 0.45 0.74];   % on-vessel / off-vessel example colours
p8.vesCol     = [0.90 0.00 0.00];                   % blood vessel (map contour + histogram)
p8.traceROIset = 'ftprnt';
p8.traceNorm  = 'robust'; % how the ROI traces are turned into dF/F:
                          %  'robust' : Pr.normTraces ./ Pr.scalingfactor -- the pipeline's
                          %             own robust dF/F fit (RobustdFFfit.mat). Verified:
                          %             a bAP is ~0.10 in these units, i.e. FRACTIONAL dF/F,
                          %             so x100 gives %%. Only defined for Result.ftprnt ROIs.
                          %  'F0'     : divide the mov_res trace by the mean-image F0 seen
                          %             through the same footprint (works for any ROI set)
p8.f0Source   = 'auto';   % baseline F0 for the 'F0' normalisation:
                          %  'auto'  : mean(mov_mc,3) if section (4) left it in memory,
                          %            otherwise Pr.ref_im
                          %  'movmc' : force mean(mov_mc,3)   'refim' : force Pr.ref_im
p8.varClim     = [2000 55000];      % colour limits of the RAW high-pass variance image.
                          % [] = auto (1st/99.5th percentile of that image).
p8.varClimCorr = [0 1000];                % colour limits of the CORRECTED one, set separately.
                          % [] = reuse p8.varClim (shared scale, so the drop is visible);
                          % give it its own range, e.g. [0 900], to see structure inside
                          % the corrected movie that the raw scale would crush.
p8.varContour  = false;             % draw the bvMask outline on the high-pass variance images
p8.pcaContour  = false;             % ... on the per-band PC variance images (R3-8-PCA)
p8.pcaClimMode = 'perimage';        % colour scaling of the per-band PC variance images:
                                    %  'perimage' : every panel gets its own range, so the
                                    %               structure inside each is visible. Brightness
                                    %               is then NOT comparable between panels.
                                    %  'perband'  : one range per band, set by the raw movie,
                                    %               shared by all three -- brightness IS
                                    %               comparable, so a band losing power shows.
p8.pcaClimPrc  = 99.5;              % upper percentile defining those ranges
p8.pcImgContour = true;             % ... on the individual PC maps (R3-8-PCimages), where it
                                    % is the point: it shows whether a PC is vessel-shaped
p8.corrShow   = 'c_raw';  % which matrix the correlation panel draws, a field of C7:
                          %  'c_raw' : corrected, NO temporal smoothing  <- the paper's claim
                          %  'c_sm'  : corrected, sub_win smoothed  (this is what 6d shows)
                          %  'in_sm' : the uncorrected sub_win matrix (the reference itself)
p8.leafOrder  = 'fig2d';    % ROW/COLUMN ORDER of the correlation matrix:
                          % 'fig2d'    : reproduce Figure 2D of Analysis_20260522_Figure2all
                          %              exactly -- cluster the SILENT-period 20 ms
                          %              cross-channel matrix of the WHOLE recording
                          %              (Pr.Subthreshold_ch via get_corrMat), flip so
                          %              basal sits on top, then take dendrogram's outperm
                          % 'dendrite' : NO reordering -- keep Pr.sortDist, which is
                          %              Result.BrancDist_order, i.e. sorted by distance
                          %              along the dendrite. Near-diagonal then means
                          %              "close together on the neuron". Interpretable.
                          % 'own'      : hierarchical-clustering leaf order of THIS window's
                          %              uncorrected 20 ms matrix (blocks, but anatomy is
                          %              scrambled -- this is what looked "weird")
                          % 'sweep'    : the clustering leaf order from section (5)/(6),
                          %              for row-for-row comparison with R3-BVsweep-corrMats

% ---------------------------------------------------------------------------
% 8-PREP: sub-ROIs on the SWC skeleton, traces from raw and corrected movies
% ---------------------------------------------------------------------------
[nrow8, ncol8] = size(Pr.ref_im);
dMask8 = point2img(Result.SWC, p8.swcRadius, [nrow8 ncol8]) > 0;   % SWC dendrite mask
if ~exist('ftprintsmallest','var') || size(ftprintsmallest,1) ~= nrow8
    Rs8 = importdata(fullfile(fpath{f}, 'PC_Result_smallerROI.mat'));
    ftprintsmallest = divide_ROIs(Rs8.ftprint_smaller .* dMask8, p8.subPx, dMask8);
end
M8 = size(ftprintsmallest, 3);

% unit-sum footprints so traces are dF in camera counts
F8   = double(ftprintsmallest);
F8   = F8 ./ max(sum(sum(F8,1),2), eps);
Fv8  = tovec(F8);                                   % pixels x M8

% frames: the section (7) window, so raw and corrected are the same frames
fr8   = b7.frames;
[~, loc8] = ismember(fr8, segFrames);
X8    = -tovec(mov_res(:,:,loc8));                  % pixels x T, spike-positive
T_raw = (X8' * Fv8)';                               % M8 x T   raw
nT8   = size(T_raw, 2);
clear X8

% corrected traces on the SAME sub-ROIs, streamed from the section (7) output.
% Only read when the panel actually shows them -- this is a ~4 GB h5 read.
T_cor = [];
if p8.subCorrected
    T_cor = zeros(M8, nT8);
    for t0 = 1:b7.traceBlk:nT8
        t1 = min(t0 + b7.traceBlk - 1, nT8);
        cc = permute(h5read(out7,'/corrected',[1 1 t0],[ncol8 nrow8 t1-t0+1]), [2 1 3]);
        T_cor(:, t0:t1) = (tovec(cc)' * Fv8)';
        clear cc
    end
end

spk8  = full(Result.spike(1, fr8)) > 0;
bad8  = imdilate(spk8, ones(1, 2*p8.Wsp+1)) > 0 | imdilate(Pr.BlueStim(fr8) > 0, ones(1,201));
gd8   = ~bad8;                                       % frames used for every statistic

% vessel overlap per sub-ROI
bvMap8 = max(Pr.bvMask, [], 3);
bvBW8  = bvMap8 > 0.2*max(bvMap8(:));
ftBin8 = ftprintsmallest > 0;
npix8  = squeeze(sum(sum(ftBin8,1),2));
nbv8   = squeeze(sum(sum(ftBin8 & bvBW8,1),2));
bvFrac8 = nbv8 ./ max(npix8,1);
onVes  = bvFrac8 >= p8.bvFracOn;
offVes = bvFrac8 <= p8.bvFracOff;
fprintf('[R3-8] %d sub-ROIs: %d on-vessel (>=%.0f%% overlap), %d off-vessel\n', ...
        M8, nnz(onVes), 100*p8.bvFracOn, nnz(offVes));

% ---- baseline F0 per sub-ROI, for the dF/F0 normalisation -----------------
% F0 is the resting fluorescence seen through the SAME unit-sum footprint that
% produced the trace, so dF/F0 is the ordinary fractional change and traces from
% different sub-ROIs (and different brightness) become directly comparable.
switch p8.f0Source
    case 'movmc', assert(exist('mov_mc','var')==1, '[R3-8] mov_mc not in memory'); F0img = mean(mov_mc,3);
    case 'refim', F0img = double(Pr.ref_im);
    otherwise
        if exist('mov_mc','var'), F0img = mean(mov_mc,3); f0src='mean(mov_mc,3)';
        else,                     F0img = double(Pr.ref_im); f0src='Pr.ref_im'; end
end
if ~exist('f0src','var'), f0src = p8.f0Source; end
F0sub = (tovec(F0img)' * Fv8)';                      % M8 x 1, counts
assert(all(F0sub > 0), '[R3-8] non-positive F0 on %d sub-ROI(s)', nnz(F0sub<=0));
fprintf('[R3-8] F0 from %s: median %.0f counts (range %.0f-%.0f)\n', ...
        f0src, median(F0sub), min(F0sub), max(F0sub));

% the three trace versions, as dF/F0 in PERCENT
% (a slow running median is removed first; that is a display choice, not a filter
%  used by any statistic -- every number below is a variance of a fast component)
dF   = T_raw - movmedian(T_raw, p8.baseWin, 2, 'omitnan');
dF   = 100 * dF ./ F0sub;                            % -> %dF/F0
dF20 = movmean(dF, p8.Wsm, 2);
dFhp = dF - dF20;                                    % what the 20 ms average discards

% the same three versions for the bleach-dr CORRECTED movie, on the SAME
% sub-ROIs / F0 / frames, so the third column of the variance panel is directly
% comparable with the first two
dFhpc = [];
if p8.subCorrected
    dFc   = T_cor - movmedian(T_cor, p8.baseWin, 2, 'omitnan');
    dFc   = 100 * dFc ./ F0sub;
    dFhpc = dFc - movmean(dFc, p8.Wsm, 2);
end

% ---- the SAME three versions for Result.ftprnt, the paper's own ROIs ------
% in7C (section 7) is already Result.ftprnt projected on exactly these frames,
% in Pr.sortDist order, so nothing has to be re-extracted. F0 uses the SAME
% (un-normalised) footprint, so the arbitrary footprint scaling cancels in dF/F0.
assert(nT7 == nT8, '[R3-8] section (7) and (8) frame counts differ (%d vs %d)', nT7, nT8);
ftFull  = Pr.Ftprnts(:,:,Pr.sortDist);
ftFullV = tovec(ftFull);
F0full = (tovec(F0img)' * ftFullV)';                  % nROI7 x 1, counts. Always needed:
                                                     % the blood-vs-high-pass panel converts
                                                     % counts^2 to (%dF/F0)^2 with it.
assert(all(F0full > 0), '[R3-8] non-positive F0 on %d Result.ftprnt ROI(s)', nnz(F0full<=0));
switch p8.traceNorm
    case 'robust'
        % the pipeline's own robust dF/F: normTraces divided by the per-ROI
        % ScaleFactor that RobustdFFfit.mat solved for. This is the quantity the
        % paper analyses, so it is the honest thing to show. It is NOT the same
        % movie as in7C (-mov_res), which is why the ScaleFactor must be applied
        % to normTraces and never to in7C: the two are in different units.
        assert(isfield(Pr,'normTraces') && isfield(Pr,'scalingfactor'), ...
               '[R3-8] Pr.normTraces / Pr.scalingfactor missing -- use p8.traceNorm=''F0''');
        sfF  = Pr.scalingfactor(Pr.sortDist);
        dFf  = 100 * Pr.normTraces(Pr.sortDist, fr8) ./ sfF;      % -> %dF/F
        normNm = 'robust \DeltaF/F';
    case 'F0'
        dFf    = 100 * in7C ./ F0full;
        normNm = '\DeltaF/F_0';
    otherwise, error('p8.traceNorm must be ''robust'' or ''F0''');
end
dFf     = dFf - movmedian(dFf, p8.baseWin, 2, 'omitnan');         % display baseline only
dF20f   = movmean(dFf, p8.Wsm, 2);
dFhpf   = dFf - dF20f;
binFull = ftFull > 0;
bvFracF = squeeze(sum(sum(binFull & bvBW8,1),2)) ./ max(squeeze(sum(sum(binFull,1),2)),1);

% ---- which set the two example traces are drawn from ---------------------
switch p8.traceROIset
    case 'ftprnt'
        TR = struct('dF',dFf,'dF20',dF20f,'dFhp',dFhpf,'bvFrac',bvFracF, ...
                    'ft',ftFull,'name','Result.ftprnt');
    case 'sub'
        if strcmp(p8.traceNorm,'robust')
            warning('[R3-8] robust dF/F exists only for Result.ftprnt ROIs -- sub-ROIs use F0');
        end
        TR = struct('dF',dF,'dF20',dF20,'dFhp',dFhp,'bvFrac',bvFrac8, ...
                    'ft',ftprintsmallest,'name','sub-ROI');
        normNm = '\DeltaF/F_0';
    otherwise, error('p8.traceROIset must be ''ftprnt'' or ''sub''');
end
TR.norm   = normNm;
TR.onVes  = TR.bvFrac >= p8.bvFracOn;
TR.offVes = TR.bvFrac <= p8.bvFracOff;
fprintf('[R3-8] example traces from %s (%d ROIs: %d on-vessel, %d off-vessel)\n', ...
        TR.name, size(TR.dF,1), nnz(TR.onVes), nnz(TR.offVes));

% ---------------------------------------------------------------------------
% 8a. TWO EXAMPLE TRACES: one ON a vessel, one OFF it
% ---------------------------------------------------------------------------
if ~isempty(p8.exROI) && (max(p8.exROI) > size(TR.dF,1) || min(p8.exROI) < 1)
    warning('[R3-8] p8.exROI = [%s] is out of range for %s (%d ROIs) -- picking automatically', ...
            num2str(p8.exROI), TR.name, size(TR.dF,1));
    p8.exROI = [];
end
if isempty(p8.exROI)
    hpVar = var(TR.dFhp(:,gd8), 0, 2, 'omitnan');
    cand_on = find(TR.onVes);
    if isempty(cand_on), [~,cand_on] = max(TR.bvFrac); end          % most-overlapping if none pass
    [~, i1] = max(hpVar(cand_on));                                  % worst-contaminated
    cand_of = find(TR.offVes);
    if isempty(cand_of), [~,cand_of] = min(TR.bvFrac); end
    [~, i2] = min(hpVar(cand_of));                                  % cleanest
    p8.exROI = [cand_on(i1) cand_of(i2)];
end
exLab = {sprintf('%s #%d ON vessel (%.0f%% overlap)',  TR.name, p8.exROI(1), 100*TR.bvFrac(p8.exROI(1))), ...
         sprintf('%s #%d OFF vessel (%.0f%% overlap)', TR.name, p8.exROI(2), 100*TR.bvFrac(p8.exROI(2)))};
fprintf('[R3-8] examples: %s | %s\n', exLab{1}, exLab{2});

% Display window. Statistics always use gd8 (peri-spike frames removed); this is
% only about what the example traces show.
if strcmp(p8.showWindow, 'silent')
    dsp8 = gd8;                                          % no spikes anywhere in the window
else
    dsp8 = ~imdilate(Pr.BlueStim(fr8) > 0, ones(1,201)); % blue-off only, spikes allowed
end
bw8 = bwlabel(dsp8);
rl8 = accumarray(bw8(bw8>0)', 1);
[~, best8] = max(rl8);
runIdx8 = find(bw8 == best8);
nShow8  = min(p8.showMs, numel(runIdx8));
if strcmp(p8.showWindow,'silent') || nShow8 >= numel(runIdx8)
    tShow8 = runIdx8(1:nShow8);
else
    % slide the window inside that run and take the position with the most spikes
    sc = movsum(double(spk8(runIdx8)), nShow8);
    [~, ctr] = max(sc);
    st = min(max(1, ctr - floor(nShow8/2)), numel(runIdx8) - nShow8 + 1);
    tShow8 = runIdx8(st:st+nShow8-1);
end
tms8   = tShow8 - tShow8(1);
fprintf('[R3-8] example traces: %s window, %d ms, %d spikes shown\n', ...
        p8.showWindow, numel(tShow8), nnz(spk8(tShow8)));
pale8  = 1 - 0.45*(1 - p8.exColr);

% ---------------------------------------------------------------------------
% FIGURE LAYOUT  (4 x 4 tiles)
%    1, 2, 5, 6 : neuron image + the two example ROIs
%    3, 4, 7, 8 : the traces -- ON raw+subth, OFF raw+subth, both diffs overlaid
%    9,10       : high-pass variance image, RAW            (stacked above ...)
%   13,14       : high-pass variance image, bleach-dr      (... its corrected twin)
%   11,15       : variance vs vessel overlap, raw | high-pass columns
%   12,16       : correlation matrix of the corrected movie
%  (the PCA moved to its own figure, R3-8-PCA)
% ---------------------------------------------------------------------------
get_figure('R3-8-final');
tl8 = tiledlayout(4, 4, 'TileSpacing','compact','Padding','compact');

% ---- panel 1: the neuron, the vessel, and the two example sub-ROIs --------
% Only what the figure argues about is drawn: the vessel (red) and the two ROIs
% whose traces are next to it. The other sub-ROIs are left out as clutter.
ax = nexttile(1, [2 2]);
imagesc(ax, Pr.ref_im, prctile(Pr.ref_im(:), [20 99.8]));
colormap(ax, gray); axis(ax,'image','off'); hold(ax,'on');
contour(ax, bvBW8, [0.5 0.5], 'Color', p8.vesCol, 'LineWidth', 1.2);   % blood vessel
coord8 = get_coord(TR.ft);
for k = 1:2                                           % the two example ROIs
    B = bwboundaries(TR.ft(:,:,p8.exROI(k)) > 0);
    for b = 1:numel(B), plot(ax, B{b}(:,2), B{b}(:,1), '-', 'Color', p8.exColr(k,:), 'LineWidth', 2.5); end
    text(ax, coord8(p8.exROI(k),1)+5, coord8(p8.exROI(k),2)-6, sprintf('#%d', p8.exROI(k)), ...
         'Color', p8.exColr(k,:), 'FontSize', 11, 'FontWeight','bold');
end
sbpx8 = 100 / Pixelsize(f);
plot(ax, ncol8-12-[sbpx8 0], (nrow8-10)*[1 1], 'w-', 'LineWidth', 3);
text(ax, ncol8-12-sbpx8/2, nrow8-100, '100 µm', 'Color','w','HorizontalAlignment','center','FontSize',9);
title(ax, sprintf('f=%d   red = blood vessel;  traces from %s', f, TR.name), 'FontSize',10);

% ---- panels 2-3: one axes per example ROI ---------------------------------
% raw and subthreshold OVERLAID at the same baseline (so the average is seen
% sitting inside the raw trace), the difference offset below.
% All three are divided by the RAW trace's own SD over the analysed frames, so
% the raw trace is a true z-score AND the other two stay in the same units --
% z-scoring each separately would force all three to unit variance and destroy
% exactly the amplitude comparison this panel exists to make.
% Three rows in ONE axes so every trace is on the same scale:
%   row 1  ON  vessel : raw (pale) + subthreshold (bold)
%   row 2  OFF vessel : raw (pale) + subthreshold (bold)
%   row 3  the two (raw - 20 ms) difference traces OVERLAID, so the amplitude
%          gap between an on- and an off-vessel ROI is read directly
ax = nexttile(3, [2 2]); hold(ax,'on');
Yb = cell(1,2);
for k = 1:2
    Yb{k} = [TR.dF(p8.exROI(k), tShow8); TR.dF20(p8.exROI(k), tShow8); TR.dFhp(p8.exROI(k), tShow8)];
    Yb{k}(1:2,:) = Yb{k}(1:2,:) - mean(Yb{k}(1,:));
    Yb{k}(3,:)   = Yb{k}(3,:)   - mean(Yb{k}(3,:));
end
if p8.showClip, kAmp = 2; else, kAmp = 1; end                  % 2 = subthreshold, 1 = raw
rowAmp = max(cellfun(@(Y) max(Y(kAmp,:)) - min(Y(kAmp,:)), Yb));   % common row height
rowAmp = max(rowAmp, max(cellfun(@(Y) max(Y(3,:)) - min(Y(3,:)), Yb)));
dRow   = -1.25 * rowAmp;
rowOff = [0, dRow, 2*dRow];
vesTag = {'ON vessel', 'OFF vessel'};
for k = 1:2                                                     % rows 1 and 2
    plot(ax, tms8, Yb{k}(1,:) + rowOff(k), '-', 'Color', [.65 .65 .65], 'LineWidth', 0.4);
    plot(ax, tms8, Yb{k}(2,:) + rowOff(k), '-', 'Color', p8.exColr(k,:), 'LineWidth', 1.6);
    text(ax, 6, rowOff(k) + 0.55*rowAmp, ...
         sprintf('%s  #%d (%.0f%% overlap)   raw + subthreshold (%d ms)', ...
                 vesTag{k}, p8.exROI(k), 100*TR.bvFrac(p8.exROI(k)), p8.Wsm), ...
         'Color', p8.exColr(k,:), 'FontSize', 8);
end
for k = 1:2                                                     % row 3: both diffs overlaid
    plot(ax, tms8, Yb{k}(3,:) + rowOff(3), '-', 'Color', p8.exColr(k,:), 'LineWidth', 0.6);
end
text(ax, 6, rowOff(3) + 0.55*rowAmp, sprintf('diff (raw - %d ms), both ROIs overlaid', p8.Wsm), ...
     'Color', [.2 .2 .2], 'FontSize', 8);
for r = 2:3, yline(ax, rowOff(r) + 0.62*rowAmp, ':', 'Color',[.85 .85 .85], 'HandleVisibility','off'); end

sb  = round(rowAmp/3, 1, 'significant');                        % scale bar, no y axis
xsb = tms8(end) - 0.02*tms8(end);
ysb = rowOff(3) - 0.55*rowAmp;
plot(ax, xsb*[1 1], ysb + [0 sb], 'k-', 'LineWidth', 2.5);
text(ax, xsb - 0.012*tms8(end), ysb + sb/2, sprintf('%g%% %s', sb, TR.norm), ...
     'HorizontalAlignment','right', 'VerticalAlignment','middle', 'FontSize', 8);
set(ax,'YTick',[],'YColor','none'); box(ax,'off');
xlim(ax,[0 tms8(end)]); ylim(ax, [ysb - 0.25*rowAmp, 0.75*rowAmp]);
xlabel(ax,'time (ms)');
title(ax, sprintf('%s traces, common scale', TR.name), 'FontSize', 9);

fprintf('[R3-8] variance decomposition of the two examples:\n');
for k = 1:2
    vR = var(TR.dF  (p8.exROI(k), gd8), 0, 2, 'omitnan');
    vS = var(TR.dF20(p8.exROI(k), gd8), 0, 2, 'omitnan');
    vH = var(TR.dFhp(p8.exROI(k), gd8), 0, 2, 'omitnan');
    fprintf('   %-38s raw=%.3g  %dms=%.3g (%.0f%%)  high-pass=%.3g (%.0f%%)   [(%%dF/F0)^2]\n', ...
            exLab{k}, vR, p8.Wsm, vS, 100*vS/vR, vH, 100*vH/vR);
end

% ---- panels 4-5: high-pass variance images, RAW above bleach-dr -----------
q8  = p8.imgFrames(p8.imgFrames <= nT8);
mq8 = single(-mov_res(:,:,loc8(q8)));
cq8 = permute(h5read(out7,'/corrected',[1 1 q8(1)],[ncol8 nrow8 numel(q8)]), [2 1 3]);
% high-pass := raw - movmean(raw, Wsm, 3), pixel-wise along time. The moving
% average runs over the CONTIGUOUS movie and only then are the analysed frames
% selected -- the same order as the traces, and it keeps spikes and blue-on
% frames from inflating the variance on the neuron.
gq8 = gd8(q8);
if ~p8.imgUseGood, gq8 = true(size(gq8)); end
hpTmp      = mq8 - movmean(mq8, p8.Wsm, 3);
hp_raw_img = var(hpTmp(:,:,gq8), 0, 3);
hpTmp      = cq8 - movmean(cq8, p8.Wsm, 3);
hp_cor_img = var(hpTmp(:,:,gq8), 0, 3);
clear hpTmp
bvpix8 = bvMap8(:) > 0;

% ONE colour scale for both images so "less high-pass power after correction"
% reads off the picture. p8.varClim overrides the automatic range.
if isempty(p8.varClim), climHP8 = prctile(hp_raw_img(:), [1 99.5]); else, climHP8 = p8.varClim; end
if isempty(p8.varClimCorr), climHP8c = climHP8; else, climHP8c = p8.varClimCorr; end
fprintf('[R3-8] high-pass variance clim: raw [%.4g %.4g], corrected [%.4g %.4g]%s\n', ...
        climHP8(1), climHP8(2), climHP8c(1), climHP8c(2), ...
        repmat('  (shared)', 1, isempty(p8.varClimCorr)));
for pnl = 1:2
    ax = nexttile(9 + 4*(pnl-1), [1 2]);
    if pnl == 1, IM = hp_raw_img; ttl = 'high-pass variance: RAW';
    else,        IM = hp_cor_img; ttl = 'high-pass variance: bleach-dr corrected'; end
    if pnl == 1, clNow = climHP8; else, clNow = climHP8c; end
    imagesc(ax, IM, clNow); axis(ax,'image','off');
    colormap(ax, turbo); hold(ax,'on');
    if p8.varContour, contour(ax, bvBW8, [0.5 0.5], 'w', 'LineWidth', 0.8); end
    title(ax, ttl, 'FontSize', 9);
    colorbar(ax);
end
if p8.imgUseGood, frTag8 = 'spike + blue-on excluded'; else, frTag8 = 'ALL frames'; end
fprintf('[R3-8] high-pass images from %d/%d frames (%s)\n', nnz(gq8), numel(q8), frTag8);
fprintf('[R3-8] high-pass variance in/out bvMask: raw %.2f -> corrected %.2f\n', ...
        mean(hp_raw_img(bvpix8))/mean(hp_raw_img(~bvpix8)), ...
        mean(hp_cor_img(bvpix8))/mean(hp_cor_img(~bvpix8)));

% ---- panel 6: variance in two columns, raw vs high-pass ------------------
% One dot per Result.ftprnt ROI in each column, faint line joining the same ROI
% across columns. Red = the ROI overlaps a blood vessel. The vessel ROIs sit at
% the top of BOTH columns, but the separation is far larger in the high-pass one:
% that gap is the contamination the 20 ms average is removing.
% This panel deliberately uses the ~30 px SUB-ROIs (dF/F0-normalised), not the
% whole Result.ftprnt ROIs: a whole ROI averages over far more dendrite than a
% vessel covers, which dilutes the contamination. The small footprints resolve it.
% Boxes rather than one dot per sub-ROI: with a few hundred sub-ROIs the cloud
% is unreadable, and the comparison that matters is between the vessel and
% non-vessel DISTRIBUTIONS. Each column is split into two boxes by that.
ax = nexttile(11, [2 1]); hold(ax,'on');
vRawF = var(dF   (:,gd8), 0, 2, 'omitnan');     % raw, dF/F0, sub-ROIs
vHpF  = var(dFhp (:,gd8), 0, 2, 'omitnan');     % what 20 ms smoothing discards
setNm = sprintf('sub-ROI (%d px)', p8.subPx);
setNorm = '\\DeltaF/F_0';
onF   = bvFrac8(:) > p8.bvFracMark;        % ANY vessel overlap among the sub-ROIs
fprintf('[R3-8] vessel overlap over %s: max %.1f%%, median %.1f%%, %d of %d with any overlap\n', ...
        setNm, 100*max(bvFrac8), 100*median(bvFrac8), nnz(onF), numel(bvFrac8));
if ~any(onF)
    warning('[R3-8] no sub-ROI overlaps bvMask at all -- the two boxes will be identical');
end

colV  = {vRawF, vHpF};
colNm = {'raw', sprintf('high-pass (raw-%dms)', p8.Wsm)};
if p8.subCorrected && ~isempty(dFhpc)
    vHpC  = var(dFhpc(:,gd8), 0, 2, 'omitnan');
    colV{end+1}  = vHpC;
    colNm{end+1} = sprintf('high-pass (bleach-dr-%dms)', p8.Wsm);
end
nCol8 = numel(colV);

xg8 = categorical(repelem((1:nCol8)', numel(vRawF), 1), 1:nCol8, colNm);
cg8 = categorical(repmat(onF, nCol8, 1), [false true], ...
                  {sprintf('no vessel overlap (n=%d)', nnz(~onF)), ...
                   sprintf('overlaps vessel (n=%d)',   nnz(onF))});
yy8 = vertcat(colV{:});
bh8 = boxchart(ax, xg8, yy8, 'GroupByColor', cg8, 'MarkerStyle','none');   % no outlier markers
if numel(bh8) == 2
    bh8(1).BoxFaceColor = [.55 .55 .55];  bh8(1).MarkerColor = [.55 .55 .55];
    bh8(2).BoxFaceColor = p8.vesCol;      bh8(2).MarkerColor = p8.vesCol;
else
    bh8(1).BoxFaceColor = [.55 .55 .55];  bh8(1).MarkerColor = [.55 .55 .55];
end
% Linear axis. Hiding the outlier markers does NOT stop them driving the axis
% limits (measured: ylim top 1.4x the whisker extent), which flattens the boxes
% on a linear scale -- so the limits are set to the whisker range explicitly.
wLo8 = inf;  wHi8 = -inf;
for c = 1:nCol8
    for g = [false true]
        v = colV{c}(onF == g);  v = v(~isnan(v));
        if numel(v) < 2, continue; end
        q = quantile(v, [.25 .75]);  iq = diff(q);
        wHi8 = max(wHi8, max(v(v <= q(2) + 1.5*iq)));
        wLo8 = min(wLo8, min(v(v >= q(1) - 1.5*iq)));
    end
end
set(ax,'XTickLabelRotation',12); grid(ax,'on');
if isfinite(wLo8) && isfinite(wHi8) && wHi8 > wLo8
    pad8 = 0.05*(wHi8 - wLo8);
    ylim(ax, [max(0, wLo8 - pad8), wHi8 + pad8]);    % variance is non-negative
end
ylabel(ax, sprintf('variance  ((%%%s)^2)', setNorm));
legend(ax,'Location','best','Box','off','FontSize',7);
title(ax, setNm, 'FontSize', 9);

if any(onF) && any(~onF)
    fprintf('[R3-8] on/off-vessel median variance ratio:');
    for c = 1:nCol8
        fprintf('  %s x%.2f', colNm{c}, median(colV{c}(onF))/median(colV{c}(~onF)));
    end
    fprintf('\n       ranksum p vs off-vessel:');
    for c = 1:nCol8
        fprintf('  %s %.2g', colNm{c}, ranksum(colV{c}(onF), colV{c}(~onF)));
    end
    fprintf('\n');
end
if p8.subCorrected && nCol8 == 3
    fprintf('[R3-8] high-pass variance kept after bleach-dr (median ratio): on-vessel %.3f, off-vessel %.3f\n', ...
            median(colV{3}(onF)./colV{2}(onF)), median(colV{3}(~onF)./colV{2}(~onF)));
end
xOv = 100*bvFrac8(:);  okOv = ~isnan(xOv) & ~isnan(vHpF);
if nnz(okOv) > 3
    fprintf('[R3-8] Spearman vs vessel overlap: high-pass rho=%.3f, raw rho=%.3f\n', ...
            corr(xOv(okOv), vHpF(okOv), 'type','Spearman'), ...
            corr(xOv(okOv), vRawF(okOv),'type','Spearman'));
end

% ---- panel 7: correlation matrix of the corrected movie ------------------
% Leaf order comes from the ORIGINAL (uncorrected, sub_win-smoothed) matrix, so
% the blocks sit where the original data put them and any change is visible as a
% change WITHIN that ordering rather than as a re-shuffle.
ax = nexttile(12, [2 1]);
if ~exist('cmapCorr','var'), cmapCorr = gen_colormap([0 0.5 1; 1 1 1; 1 0 0], 256); end
% leaf order. 'own' clusters THIS window's uncorrected sub_win matrix; 'sweep'
% reproduces section (6)'s order. The sweep order is rebuilt from Cmat_in in the
% sweep .mat rather than from a stale `ordS` in the workspace, so it works even
% if section (6) was never run in this session -- and if it cannot be rebuilt it
% says so loudly instead of quietly falling back.
if strcmp(p8.leafOrder,'fig2d')
    % --- verbatim reproduction of the Figure 2D ordering ---------------------
    % Source: Analysis_20260522_Figure2all.m, cell "Figure 2D". It clusters a
    % DIFFERENT matrix from anything else here: the silent-period 20 ms
    % cross-channel correlation of the WHOLE recording, not this 10k window. The
    % basal-on-top flip and the final dendrogram('Reorder') step are part of the
    % order, so both are reproduced rather than approximated.
    spkA = Pr.SpikeMat(1,:) > 0;
    goodA = ~imdilate(Pr.BlueStim > 0, ones(1,201));
    perispkA = unique([tovec(find(spkA)' + (-3:10)); find(Pr.CStrace)']);
    perispkA(perispkA < 1 | perispkA > numel(spkA)) = [];
    silentA = true(1, numel(spkA));  silentA(perispkA) = false;
    silentMaskA = goodA & silentA;
    sub1A = Pr.Subthreshold_ch{1};  sub1A(isnan(Pr.NormalizedTrace_ch{1})) = NaN;
    sub2A = Pr.Subthreshold_ch{2};  sub2A(isnan(Pr.NormalizedTrace_ch{2})) = NaN;
    CA = get_corrMat(sub1A, sub2A, find(silentMaskA));
    nA = size(CA,1);
    assert(nA == size(C7.in_sm,1), ...
           '[R3-8] Figure 2D matrix has %d ROIs but C7 has %d -- cannot reuse that order', ...
           nA, size(C7.in_sm,1));
    DA = 1 - CA;  DA = (DA+DA')/2;  DA(1:nA+1:end) = 0;  DA(isnan(DA)) = 1;
    dvA = squareform(DA);  ZA = linkage(dvA, 'average');
    try
        ordA = optimalleaforder(ZA, dvA);
    catch
        ordA = 1:nA;
    end
    labA = Pr.labelvec;
    if mean(find(labA(ordA)==1)) > (nA+1)/2, ordA = fliplr(ordA); end   % basal on top
    fTmp = figure('Visible','off');
    [~,~,ord8] = dendrogram(ZA, 0, 'Orientation','left', 'Reorder', ordA);
    close(fTmp);
    ordSrc = 'Figure 2D leaf order (silent 20 ms, whole recording)';
    fprintf('[R3-8] correlation matrix in FIGURE 2D order (%d ROIs, basal-on-top flip applied)\n', nA);
elseif strcmp(p8.leafOrder,'dendrite')
    % No reordering at all. The rows of C7 are already in Pr.sortDist order,
    % which is Result.BrancDist_order -- distance along the dendrite. Keeping it
    % means position i of the matrix IS the i-th ROI along the neuron.
    ord8 = 1:size(C7.in_sm,1);
    ordSrc = 'dendritic distance (Pr.sortDist, no clustering)';
    fprintf('[R3-8] correlation matrix in DENDRITIC DISTANCE order (no clustering)\n');
end
leafSrcMat = C7.in_sm;  ordSrc0 = 'own window, uncorrected 20 ms';
if ~any(strcmp(p8.leafOrder, {'dendrite','fig2d'}))
if strcmp(p8.leafOrder,'sweep')
    Csweep = [];
    if exist('Cmat_in','var') && size(Cmat_in,1) == size(C7.in_sm,1)
        Csweep = Cmat_in;                                  % already in the workspace
    elseif exist('sweepFile','var') && isfile(sweepFile)
        Ssw = load(sweepFile, 'Cmat_in');                  % rebuild it from the sweep file
        if isfield(Ssw,'Cmat_in') && size(Ssw.Cmat_in,1) == size(C7.in_sm,1)
            Csweep = Ssw.Cmat_in;
        end
    end
    if isempty(Csweep)
        warning(['[R3-8] p8.leafOrder=''sweep'' requested but the sweep Cmat_in is not ' ...
                 'available (run section (5)/(6), or check the ROI count) -- using this ' ...
                 'window''s own order instead']);
    else
        leafSrcMat = Csweep;  ordSrc0 = 'section-(6) sweep leaf order';
    end
end
Dd8 = 1 - leafSrcMat;  Dd8 = (Dd8+Dd8')/2;  Dd8(1:size(Dd8,1)+1:end) = 0;  Dd8(isnan(Dd8)) = 1;
dv8 = squareform(Dd8);  Z8 = linkage(dv8,'average');
try
    ord8 = optimalleaforder(Z8, dv8);
catch
    fT = figure('Visible','off'); [~,~,ord8] = dendrogram(Z8,0); close(fT);
end
ordSrc = ordSrc0;
fprintf('[R3-8] leaf order = %s; it keeps %d/%d ROIs in dendritic-distance position\n', ...
        ordSrc, sum(ord8(:) == (1:numel(ord8))'), numel(ord8));
end
fprintf('[R3-8] matrix row order (first 12 ROI indices): %s ...\n', ...
        num2str(ord8(1:min(12,numel(ord8)))));
Cshow8 = C7.(p8.corrShow);
switch p8.corrShow
    case 'c_raw', shTtl = 'corrected, unsmoothed';   shR = r_rawVsSm;
    case 'c_sm',  shTtl = sprintf('corrected, %d ms', scr.sub_win);  shR = r_smVsSm;
    case 'in_sm', shTtl = sprintf('UNcorrected, %d ms', scr.sub_win); shR = 1;
end
imagesc(ax, Cshow8(ord8,ord8), [-0.25 1]); axis(ax,'image'); colormap(ax, cmapCorr);
set(ax,'XTick',[],'YTick',[]); colorbar(ax);
title(ax, sprintf('%s\n(%s)  r = %.3f', shTtl, ordSrc, shR), 'FontSize', 9);

% Why this can look unlike 6d's R3-BVsweep-corrMats: that figure draws the
% sub_win-SMOOTHED matrix, on the short sweep window, in the sweep window's own
% leaf order. All three differ here. p8.corrShow='c_sm' and p8.leafOrder='sweep'
% make them comparable (the frame count still differs).
fprintf(['[R3-8] correlation panel: %s, %s, %d frames.\n' ...
         '       r(corrected unsmoothed, ORIGINAL %dms) = %.3f | r(corrected %dms, ORIGINAL %dms) = %.3f\n' ...
         '       baseline r(UNcorrected raw, ORIGINAL %dms) = %.3f\n'], ...
        p8.corrShow, ordSrc, nT8, scr.sub_win, r_rawVsSm, scr.sub_win, scr.sub_win, r_smVsSm, ...
        scr.sub_win, r_inRawVsSm);

%% ---------------------------------------------------------------------------
% 8f. PCA -- ITS OWN FIGURE (R3-8-PCA)
%   rows: raw / 20 ms smoothed / bleach-dr corrected
%   cols: the spatial variance carried by each PC band, plus the full spectrum
%
% The variance image of a band needs no reconstruction: with A = U S V' and the
% columns of V orthonormal, the temporal variance of pixel i restricted to the
% components in band B is
%       var_i(B) = 1/(T-1) * sum_{k in B} ( U(i,k) * S_k )^2
% (cross terms vanish because the V columns are orthogonal). Each band uses ONE
% colour scale across the three movies, so losing power in a band is visible.
% ---------------------------------------------------------------------------
dpix8 = find(dMask8);
kPCA  = 1:min(p8.nPcaFrames, numel(q8));
Praw  = reshape(mq8, [], numel(q8));  Praw = Praw(dpix8, kPCA);
Psm   = movmean(Praw, p8.Wsm, 2);
Pcor  = reshape(cq8,  [], numel(q8)); Pcor = Pcor(dpix8, kPCA);
clear mq8 cq8
gPCA  = gd8(q8(kPCA));  nUse8 = nnz(gPCA);
nMax  = max([cellfun(@max, p8.pcBands), cellfun(@max, p8.pcShow)]);
fprintf('[R3-8] PCA: %d dendrite-mask pixels x %d frames -> %d components\n', ...
        numel(dpix8), nUse8, min(numel(dpix8), nUse8));

pcaNm  = {'raw', sprintf('%d ms smoothed',p8.Wsm), 'bleach-dr corrected'};
evList = cell(1,3);  vbList = cell(3, numel(p8.pcBands));  UkList = cell(1,3);  skList = cell(1,3);
for i = 1:3
    switch i
        case 1, A = Praw;
        case 2, A = Psm;
        case 3, A = Pcor;
    end
    A = double(A(:, gPCA));  A = A - mean(A, 2);
    [U8, S8, V8] = svd(A, 'econ');
    sv = diag(S8);
    if strcmp(p8.pcaCov, 'lag1')
        % Lag-1 covariance C = A*G*A'. Forming C is nPix^2 (3.5 GB at ~21k px),
        % but with A = U*S*V' it collapses to a T x T eigenproblem:
        %     C = U * [ S*(V1'V2 + V2'V1)/2*S ] * U'
        % which is exact -- verified to 4e-15 against the direct eigendecomposition.
        V1 = V8(1:end-1,:);  V2 = V8(2:end,:);
        G8 = (V1'*V2 + V2'*V1) / 2 / max(size(V8,1)-2, 1);
        K8 = (sv .* G8) .* sv';   K8 = (K8 + K8')/2;
        [W8, lam] = eig(K8, 'vector');
        [lam, si8] = sort(lam, 'descend');   W8 = W8(:, si8);
        Uall = U8 * W8;                                   % orthonormal, eigvecs of C
        clear V1 V2 G8 K8 W8
    else
        lam  = sv.^2 / max(nUse8-1, 1);                   % ordinary covariance eigenvalues
        Uall = U8;
    end
    clear U8 S8 V8 A
    posSum = sum(lam(lam > 0));                           % lam can be negative for 'lag1'
    evList{i} = lam / posSum;
    kKeep = 1:min(nMax, numel(lam));
    Uk = Uall(:, kKeep);  sk = lam(kKeep);                % sk now holds EIGENVALUES
    clear Uall
    for b = 1:numel(p8.pcBands)                       % per-band variance image
        B = p8.pcBands{b};  B = B(B <= numel(sk));
        B = B(sk(B) > 0);                             % negative dims carry no variance
        vimg = zeros(nrow8, ncol8);
        if ~isempty(B)
            vimg(dpix8) = sum((Uk(:,B).^2) .* sk(B)', 2);
        end
        vbList{i,b} = vimg;
    end
    UkList{i} = Uk;  skList{i} = sk;          % kept for the individual-PC maps below
    clear Uk sk
    fprintf('   %-22s %d components, %d with positive eigenvalue\n', ...
            pcaNm{i}, numel(lam), nnz(lam > 0));
end
clear Praw Psm Pcor

get_figure('R3-8-PCA');
tlP = tiledlayout(4, numel(p8.pcBands), 'TileSpacing','compact','Padding','compact');
bandLab = cellfun(@(B) sprintf('PC %d-%d', min(B), max(B)), p8.pcBands, 'UniformOutput', false);
fprintf('[R3-8] PC band images scaled %s (upper percentile %.3g)\n', p8.pcaClimMode, p8.pcaClimPrc);
for b = 1:numel(p8.pcBands)
    climB = [0 prctile(vbList{1,b}(dpix8), p8.pcaClimPrc)];   % 'perband': set by the RAW movie
    if climB(2) <= 0, climB(2) = eps; end
    for i = 1:3
        ax = nexttile((i-1)*numel(p8.pcBands) + b);
        if strcmp(p8.pcaClimMode, 'perimage')
            climB = [0 prctile(vbList{i,b}(dpix8), p8.pcaClimPrc)];
            if climB(2) <= 0, climB(2) = eps; end
        end
        imagesc(ax, vbList{i,b}, climB); axis(ax,'image','off'); colormap(ax, turbo);
        hold(ax,'on');
        if p8.pcaContour, contour(ax, bvBW8, [0.5 0.5], 'w', 'LineWidth', 0.5); end
        if i == 1, title(ax, bandLab{b}, 'FontSize', 9); end
        if b == 1
            text(ax, -0.02, 0.5, pcaNm{i}, 'Units','normalized', 'Rotation',90, ...
                 'HorizontalAlignment','center', 'FontSize', 9, 'FontWeight','bold');
        end
        if strcmp(p8.pcaClimMode,'perimage') || b == numel(p8.pcBands), colorbar(ax); end
    end
end
ax = nexttile(3*numel(p8.pcBands) + 1, [1 numel(p8.pcBands)]); hold(ax,'on');
cmap8 = [0 0 0; 0.85 0.33 0.10; 0.00 0.45 0.74];
hP8 = gobjects(1,3);
for i = [3 2 1]                                        % raw drawn LAST so it stays visible
    ev = evList{i};
    hP8(i) = plot(ax, 1:numel(ev), 100*ev, '-', 'Color', cmap8(i,:), ...
                  'LineWidth', 1.5 + 0.6*(i==1), 'DisplayName', pcaNm{i});
end
for b = 1:numel(p8.pcBands)
    xline(ax, min(p8.pcBands{b}), ':', 'Color',[.7 .7 .7], 'HandleVisibility','off');
end
set(ax,'XScale','log','YScale','log'); grid(ax,'on');
xlim(ax, [1 max(cellfun(@numel, evList))]);
xlabel(ax,'PC #'); ylabel(ax,'variance explained (%)');
legend(ax, hP8, pcaNm, 'Location','southwest','Box','off');
if strcmp(p8.pcaClimMode,'perimage'), climTag = 'each panel scaled to itself';
else,                                  climTag = 'one colour scale per band'; end
title(tlP, sprintf(['f=%d  PCA inside the SWC dendrite mask (%d px, %d frames -> %d PCs)   ' ...
                    'per-band spatial variance; %s'], ...
                    f, numel(dpix8), nUse8, min(numel(dpix8), nUse8), climTag), 'Interpreter','none');

fprintf('[R3-8] variance explained per band (%% of total):\n');
fprintf('   %-22s %s\n', 'movie', strjoin(bandLab, '  '));
for i = 1:3
    ev = evList{i};
    vals = cellfun(@(B) 100*sum(ev(B(B<=numel(ev)))), p8.pcBands);
    fprintf('   %-22s %s\n', pcaNm{i}, sprintf('%8.2f', vals));
end
fprintf('[R3-8] band variance inside bvMask / outside (raw -> corrected):\n');
bvIn = bvMap8 > 0;
for b = 1:numel(p8.pcBands)
    r1 = mean(vbList{1,b}(bvIn & dMask8)) / mean(vbList{1,b}(~bvIn & dMask8));
    r3 = mean(vbList{3,b}(bvIn & dMask8)) / mean(vbList{3,b}(~bvIn & dMask8));
    fprintf('   %-10s  %.2f -> %.2f\n', bandLab{b}, r1, r3);
end

% ---------------------------------------------------------------------------
% 8g. dF/F VARIANCE OF EACH ROI, PAINTED ON THE DENDRITE MAP  (R3-8-ROImap)
%     Every Result.ftprnt ROI is filled with its own variance. Overlapping
%     footprints are resolved by giving each pixel to the ROI with the largest
%     weight there, so the map is a clean partition of the dendrite.
% ---------------------------------------------------------------------------
[wMax8, aMax8] = max(TR.ft, [], 3);      % which ROI owns each pixel
supp8 = wMax8 > 0;
vRawM = var(TR.dF  (:,gd8), 0, 2, 'omitnan');
vHpM  = var(TR.dFhp(:,gd8), 0, 2, 'omitnan');
mapsM = cell(1,3);
for m = 1:3
    switch m
        case 1, vv = vRawM;
        case 2, vv = vHpM;
        case 3, vv = 100*vHpM./vRawM;
    end
    IMm = nan(nrow8, ncol8);  IMm(supp8) = vv(aMax8(supp8));
    mapsM{m} = IMm;
end
namesM = {'raw variance', sprintf('high-pass variance (raw - %d ms)', p8.Wsm), ...
          'high-pass / raw  (%)'};

get_figure('R3-8-ROImap');
tlM = tiledlayout(3, 1, 'TileSpacing','compact','Padding','compact');
for m = 1:3
    ax = nexttile(m);
    IM = mapsM{m};
    if m < 3, cl = [0 prctile(IM(supp8), 99)]; else, cl = [0 100]; end
    if cl(2) <= cl(1), cl(2) = cl(1) + eps; end
    imagesc(ax, IM, 'AlphaData', supp8); set(ax,'Color',[.12 .12 .12]);
    clim(ax, cl); axis(ax,'image','off'); colormap(ax, turbo); hold(ax,'on');
    contour(ax, bvBW8, [0.5 0.5], 'Color', p8.vesCol, 'LineWidth', 1.2);
    colorbar(ax);
    title(ax, sprintf('%s   (%s, %s)', namesM{m}, TR.name, TR.norm), 'FontSize', 9);
end
title(tlM, sprintf('f=%d  per-ROI %s variance on the dendrite (red = blood vessel)', f, TR.norm));

%% ---------------------------------------------------------------------------
% 8h. INDIVIDUAL PC SPATIAL MAPS  (R3-8-PCimages)
%     The columns of U reshaped back onto the dendrite mask. Sign is arbitrary
%     in an SVD, so each map is flipped to make its strongest pixel positive.
%     Low PCs should look like the neuron; if the blood is a distinct source it
%     shows up as vessel-shaped structure in some of them.
% ---------------------------------------------------------------------------
iShow  = unique(min(max(1, p8.pcShowMovie(:)'), 3), 'stable');
nColPC = max(cellfun(@numel, p8.pcShow));
nRowPC = numel(iShow) * numel(p8.pcShow);
get_figure('R3-8-PCimages');
tlPC = tiledlayout(nRowPC, nColPC, 'TileSpacing','compact','Padding','compact');
rowPC = 0;
for mi = 1:numel(iShow)
    i   = iShow(mi);
    UkS = UkList{i};
    for g = 1:numel(p8.pcShow)
        rowPC = rowPC + 1;
        pcs = p8.pcShow{g};  pcs = pcs(pcs <= size(UkS,2));
        for c = 1:nColPC
            ax = nexttile((rowPC-1)*nColPC + c);
            if c > numel(pcs), axis(ax,'off'); continue; end
            k  = pcs(c);
            u  = UkS(:,k);
            [~, imx] = max(abs(u));  u = u * sign(u(imx));       % fix the arbitrary sign
            im = nan(nrow8, ncol8);  im(dpix8) = u;
            cl = prctile(abs(u), 99.5);  if cl <= 0, cl = eps; end
            imagesc(ax, im, 'AlphaData', dMask8); set(ax,'Color',[.12 .12 .12]);
            clim(ax, [-cl cl]); axis(ax,'image','off'); colormap(ax, cmapCorr); hold(ax,'on');
            if p8.pcImgContour
                contour(ax, bvBW8, [0.5 0.5], 'Color', p8.vesCol, 'LineWidth', 0.7);
            end
            title(ax, sprintf('PC %d  (%.2f%%)', k, 100*evList{i}(k)), 'FontSize', 8);
            if c == 1
                text(ax, -0.03, 0.5, pcaNm{i}, 'Units','normalized', 'Rotation',90, ...
                     'HorizontalAlignment','center', 'FontSize',8, 'FontWeight','bold');
            end
        end
    end
end
title(tlPC, sprintf(['f=%d  individual PC maps (red = blood vessel, SVD sign fixed so the ' ...
                     'strongest pixel is positive)'], f), 'Interpreter','none');

figure(findobj('Type','figure','Name','R3-8-final'));
title(tl8, sprintf('f=%d  frames %d:%d -- the fast subthreshold content is vascular (%s)', ...
      f, fr8(1), fr8(end), tag7), 'Interpreter','none');

%% ################################################################
%  (9) cvPCA -- NOISE-CORRECTED EIGENSPECTRUM  (Stringer et al. 2019, Fig. 2a)
%
%  Section (8) plots ordinary PCA of a single movie. That estimator is BIASED:
%  independent shot noise adds an approximately FLAT floor to the eigenvalue
%  spectrum, so the long shallow tail there is largely noise, and the raw-vs-
%  corrected comparison is confounded because the corrected movie has been
%  through PMD (a low-rank denoiser) while the raw one has not.
%
%  cvPCA fixes both. It needs two measurements of the SAME signal with
%  INDEPENDENT noise; here that is the odd/even checkerboard pixel halves of
%  each footprint, which section (7) already projected:
%        in7_1 / in7_2   uncorrected, odd / even
%        tr7_1 / tr7_2   bleach-dr corrected, odd / even
%  The two halves are simultaneous, so every real optical signal -- INCLUDING
%  the blood -- is shared and survives; only shot noise cancels. That is exactly
%  what is wanted: the spectrum of real variance, with the noise floor removed,
%  so raw and corrected can be compared on equal terms.
%
%  READ THE CAVEAT AT THE END BEFORE QUOTING AN EXPONENT.
% ################################################################
f  = 20;
Pr = PCresult{f};
for req9 = {'in7_1','in7_2','tr7_1','tr7_2','b7','scr'}
    assert(exist(req9{1},'var')==1, '[R3-9] %s is missing -- run section (7) first', req9{1});
end

% ---------------------------------------------------------------------------
% 9a. PARAMETERS
% ---------------------------------------------------------------------------
p9 = struct();
p9.roiSet     = 'ftprnt';   % 'ftprnt' : Result.ftprnt ROIs (what was asked for).
                            %            NOTE this caps the spectrum at ~nROI points.
                            % 'sub'    : the ~30 px sub-ROIs -- hundreds of dimensions,
                            %            a far longer spectrum, at the cost of noisier rows.
p9.subPx      = 30;         % sub-ROI size when roiSet = 'sub'
p9.swcRadius  = 3;
p9.useSub     = true;       % run on the SUBTHRESHOLD traces (spikes interpolated out),
                            % matching the correlation matrices. false = full traces.
p9.smoothMs   = [1 20];     % temporal windows to compare (1 = none)
p9.nSplit     = 10;         % random sample half-splits inside cvPCA
p9.fitRange   = [2 Inf];    % dimensions used for the power-law fit (Inf = last positive one)
p9.subsetFrac = [0.25 0.5 0.75 1];   % ROI fractions for the Stringer Fig. 2g-j check
p9.nSubsetRep = 5;
p9.frames     = [];         % [] = all of b7.frames

% ---------------------------------------------------------------------------
% 9b. TRACES: odd / even channels, dF/F0, for uncorrected and corrected
% ---------------------------------------------------------------------------
[nrow9, ncol9] = size(Pr.ref_im);
if exist('mov_mc','var'), F0img9 = mean(mov_mc,3); else, F0img9 = double(Pr.ref_im); end
oddPix9 = mod((1:nrow9)' + (1:ncol9), 2) == 1;

switch p9.roiSet
    case 'ftprnt'
        ftA = Pr.Ftprnts(:,:,Pr.sortDist);
        A1 = in7_1;  A2 = in7_2;  B1 = tr7_1;  B2 = tr7_2;      % already projected in (7)
    case 'sub'
        dM9 = point2img(Result.SWC, p9.swcRadius, [nrow9 ncol9]) > 0;
        if ~exist('ftprintsmallest','var') || size(ftprintsmallest,1) ~= nrow9
            Rs9 = importdata(fullfile(fpath{f},'PC_Result_smallerROI.mat'));
            ftprintsmallest = divide_ROIs(Rs9.ftprint_smaller .* dM9, p9.subPx, dM9);
        end
        ftA = double(ftprintsmallest);
        error(['[R3-9] roiSet=''sub'' needs the odd/even sub-ROI traces, which section (7) ' ...
               'does not build. Project mov_res and the corrected movie onto ' ...
               'tovec(ftA.*oddPix9) / tovec(ftA.*~oddPix9) first.']);
    otherwise
        error('p9.roiSet must be ''ftprnt'' or ''sub''');
end
ftO9 = tovec(ftA .*  oddPix9);   ftE9 = tovec(ftA .* ~oddPix9);
F0o9 = (tovec(F0img9)' * ftO9)';  F0e9 = (tovec(F0img9)' * ftE9)';
assert(all(F0o9 > 0) && all(F0e9 > 0), '[R3-9] non-positive F0 on some channel/ROI');

fr9 = 1:size(A1,2);  if ~isempty(p9.frames), fr9 = p9.frames; end
spk9  = full(Result.spike(1, b7.frames(fr9))) > 0;
good9 = ~imdilate(Pr.BlueStim(b7.frames(fr9)) > 0, ones(1,scr.blueDilate));

% dF/F0 per channel (each channel divided by ITS OWN F0, so the two are comparable)
mk9 = @(X, F0) 100 * X(:, fr9) ./ F0;
condNm9 = {}; XX = {};
for w = p9.smoothMs
    for c = 1:2
        if c == 1, Y1 = mk9(A1, F0o9);  Y2 = mk9(A2, F0e9);  nm = 'uncorrected';
        else,      Y1 = mk9(B1, F0o9);  Y2 = mk9(B2, F0e9);  nm = 'bleach-dr';
        end
        if p9.useSub
            Y1 = get_subthreshold(Y1, spk9, scr.sub_dilate, w);
            Y2 = get_subthreshold(Y2, spk9, scr.sub_dilate, w);
        elseif w > 1
            Y1 = movmean(Y1, w, 2);  Y2 = movmean(Y2, w, 2);
        end
        Y1(:, ~good9) = NaN;  Y2(:, ~good9) = NaN;            % analysed frames only
        XX{end+1} = {Y1, Y2};                                  %#ok<SAGROW>
        if w == 1, condNm9{end+1} = sprintf('%s, no smoothing', nm);        %#ok<SAGROW>
        else,      condNm9{end+1} = sprintf('%s, %d ms', nm, w);            %#ok<SAGROW>
        end
    end
end
nC9 = numel(XX);
fprintf('[R3-9] %s: %d ROIs x %d frames (%d analysed), %d conditions\n', ...
        p9.roiSet, size(A1,1), numel(fr9), nnz(good9), nC9);

% ---------------------------------------------------------------------------
% 9c. cvPCA, and ordinary PCA on the same data for comparison
% ---------------------------------------------------------------------------
ssCV = cell(1,nC9);  evOrd = cell(1,nC9);  alphaCV = nan(1,nC9);  alphaOrd = nan(1,nC9);
nPosCV = nan(1,nC9);
rng(0);
for c = 1:nC9
    [ssCV{c}, ~, ~] = cvPCA(XX{c}{1}, XX{c}{2}, p9.nSplit);
    % ordinary PCA on the COMBINED (odd+even) trace, i.e. what section (8) does
    Z = (XX{c}{1} + XX{c}{2}) / 2;
    Z = Z(:, all(isfinite(Z),1));
    sv = svd(Z - mean(Z,2), 'econ');
    evOrd{c} = sv.^2 / sum(sv.^2);
    % power-law fit over the dimensions where cvPCA is still positive
    lastPos = find(ssCV{c} > 0, 1, 'last');
    k1 = max(1, p9.fitRange(1));
    k2 = min([numel(ssCV{c}), p9.fitRange(2), lastPos]);
    nPosCV(c) = lastPos;
    if k2 > k1 + 2
        kk = (k1:k2)';  kk = kk(ssCV{c}(kk) > 0);
        pf = [ones(numel(kk),1) log10(kk)] \ log10(ssCV{c}(kk));
        alphaCV(c) = -pf(2);
        kko = kk(kk <= numel(evOrd{c}));
        pfo = [ones(numel(kko),1) log10(kko)] \ log10(evOrd{c}(kko));
        alphaOrd(c) = -pfo(2);
    end
end

% ---------------------------------------------------------------------------
% 9d. FIGURE
% ---------------------------------------------------------------------------
get_figure('R3-9-cvPCA');
tl9 = tiledlayout(1, 2, 'TileSpacing','compact','Padding','compact');
cmap9 = lines(nC9);

ax = nexttile(1); hold(ax,'on');
h9 = gobjects(1,nC9);
for c = 1:nC9
    y = ssCV{c};  kpos = find(y > 0);
    h9(c) = plot(ax, kpos, y(kpos), '-o', 'Color', cmap9(c,:), 'MarkerSize',3, ...
                 'MarkerFaceColor', cmap9(c,:), 'LineWidth',1.4, 'DisplayName', condNm9{c});
    plot(ax, 1:numel(evOrd{c}), evOrd{c}, ':', 'Color', cmap9(c,:), 'LineWidth',1, ...
         'HandleVisibility','off');                          % ordinary PCA, dotted
end
kref = 1:max(cellfun(@numel, ssCV));
plot(ax, kref, kref.^-1 / sum(kref.^-1), 'k--', 'LineWidth',1, 'DisplayName','1/n');
set(ax,'XScale','log','YScale','log'); grid(ax,'on');
xlabel(ax,'PC dimension'); ylabel(ax,'fraction of reliable variance');
legend(ax,'Location','southwest','Box','off','FontSize',7);
title(ax, sprintf('cvPCA (solid) vs ordinary PCA (dotted), %s', p9.roiSet), 'FontSize',9);

% Stringer Fig. 2g-j analogue: does the exponent depend on how many ROIs went in?
ax = nexttile(2); hold(ax,'on');
nR9 = size(XX{1}{1},1);
aSub = nan(nC9, numel(p9.subsetFrac), p9.nSubsetRep);
for c = 1:nC9
    for q = 1:numel(p9.subsetFrac)
        m = max(4, round(p9.subsetFrac(q)*nR9));
        for r = 1:p9.nSubsetRep
            idx = randperm(nR9, m);
            s2 = cvPCA(XX{c}{1}(idx,:), XX{c}{2}(idx,:), max(3,round(p9.nSplit/2)));
            lp = find(s2 > 0, 1, 'last');
            kk = (2:min(lp, numel(s2)))';  kk = kk(s2(kk) > 0);
            if numel(kk) > 3
                pf = [ones(numel(kk),1) log10(kk)] \ log10(s2(kk));
                aSub(c,q,r) = -pf(2);
            end
        end
    end
    mA = mean(aSub(c,:,:), 3, 'omitnan');  sA = std(aSub(c,:,:), 0, 3, 'omitnan');
    errorbar(ax, p9.subsetFrac*nR9, mA, sA, '-o', 'Color', cmap9(c,:), 'LineWidth',1.4, ...
             'MarkerFaceColor', cmap9(c,:), 'DisplayName', condNm9{c});
end
yline(ax, 1, 'k--', 'HandleVisibility','off');
set(ax,'XScale','log'); grid(ax,'on');
xlabel(ax,'number of ROIs analysed'); ylabel(ax,'fitted power-law exponent \alpha');
legend(ax,'Location','best','Box','off','FontSize',7);
title(ax, 'exponent vs sample size (Stringer Fig. 2g-j check)', 'FontSize',9);
title(tl9, sprintf('f=%d  cvPCA on the odd/even channels, frames %d:%d', ...
      f, b7.frames(fr9(1)), b7.frames(fr9(end))));

% ---------------------------------------------------------------------------
% 9e. NUMBERS AND THE CAVEAT
% ---------------------------------------------------------------------------
fprintf('\n[R3-9] ===== cvPCA =====\n');
fprintf('   %-28s %8s %8s %10s %10s\n', 'condition','alpha_cv','alpha_pca','#pos dims','PC1 frac');
for c = 1:nC9
    fprintf('   %-28s %8.2f %8.2f %10d %10.3f\n', condNm9{c}, alphaCV(c), alphaOrd(c), ...
            nPosCV(c), ssCV{c}(1));
end
fprintf(['\n   alpha_cv  : power-law exponent from the NOISE-CORRECTED spectrum\n' ...
         '   alpha_pca : the same fit on ordinary PCA of the combined trace -- the biased\n' ...
         '               estimator section (8) uses. The gap between the two columns is\n' ...
         '               how much of that spectrum was shot noise.\n' ...
         '   #pos dims : dimensions where cvPCA is still positive, i.e. where reliable\n' ...
         '               signal remains. Beyond it the estimator fluctuates about zero.\n']);
fprintf(['\n   CAVEAT -- do not quote an exponent from %d ROIs as a power law. Stringer\n' ...
         '   et al. Fig. 2g-j show the usable range and the fitted exponent only converge\n' ...
         '   as unit count grows; they used ~10,000 neurons over ~2 decades. With %d ROIs\n' ...
         '   the spectrum has at most %d points. The right-hand panel is that check on\n' ...
         '   this data: if alpha is still moving at the largest subset, it has not\n' ...
         '   converged and the number is a lower bound on the true exponent, nothing more.\n' ...
         '   Use p9.roiSet = ''sub'' for a longer spectrum.\n'], nR9, nR9, nR9);

%% ################################################################
%  (10) THREE COVARIANCE ESTIMATORS -- UNCORRECTED TRACE
%
%  All three come from the SAME ROIs, frames and normalisation, so any
%  difference is the estimator, not the data:
%
%   (i)  X(t)X(t)      same-time, one channel. Signal + shot noise. The noise
%                      sits ENTIRELY on the diagonal, inflating every variance
%                      and therefore deflating every correlation.
%   (ii) X_odd X_even  same time, two independent pixel halves. Shot noise
%                      cancels in expectation, diagonal included. Keeps ALL real
%                      optical signal at every timescale -- blood included.
%   (iii) X(t)X(t+1)   one-frame lag. White noise has zero lag-1 covariance, so
%                      this is unbiased too, BUT it also attenuates real signal
%                      that decorrelates within one frame, by about exp(-1/tau).
%
%  (ii) and (iii) are both noise-free; (ii) minus (iii) is the real signal at
%  sub-frame timescales. Section (10b) repeats all of it on the bleach-dr trace.
%  Both call covEstimators3 / plotCovEstimators3 so they cannot drift apart.
% ################################################################
f  = 20;
Pr = PCresult{f};
for req10 = {'in7_1','in7_2','tr7_1','tr7_2','b7','scr'}
    assert(exist(req10{1},'var')==1, '[R3-10] %s is missing -- run section (7) first', req10{1});
end

% ---------------------------------------------------------------------------
% 10a. PARAMETERS  (shared by 10 and 10b)
% ---------------------------------------------------------------------------
p10 = struct();
p10.useSub   = true;           % subthreshold traces (spikes interpolated out)
p10.smoothMs = 1;              % 1 = no temporal smoothing
p10.nPC      = 5;              % PC maps to draw per estimator
p10.corrClim = [-0.25 1];
p10.paintPrc = 99;             % percentile setting the PC-map colour range

% ---- traces, shared by both sections --------------------------------------
[nrow10, ncol10] = size(Pr.ref_im);
if exist('mov_mc','var'), F0img10 = mean(mov_mc,3); else, F0img10 = double(Pr.ref_im); end
oddPix10 = mod((1:nrow10)' + (1:ncol10), 2) == 1;
ftT      = Pr.Ftprnts(:,:,Pr.sortDist);
ftO10    = tovec(ftT .*  oddPix10);   ftE10 = tovec(ftT .* ~oddPix10);
F0o10    = (tovec(F0img10)' * ftO10)';  F0e10 = (tovec(F0img10)' * ftE10)';
assert(all(F0o10 > 0) && all(F0e10 > 0), '[R3-10] non-positive F0 on some channel/ROI');

spk10  = full(Result.spike(1, b7.frames)) > 0;
good10 = ~imdilate(Pr.BlueStim(b7.frames) > 0, ones(1,scr.blueDilate));

% one place that turns a raw/corrected pair into the analysed odd/even traces
srcRaw10  = {in7_1, in7_2};
srcCorr10 = {tr7_1, tr7_2};
prep10 = cell(1,2);
for q = 1:2
    if q == 1, R = srcRaw10; else, R = srcCorr10; end
    Y1 = 100 * R{1} ./ F0o10;      Y2 = 100 * R{2} ./ F0e10;    % %dF/F0 per channel
    if p10.useSub
        Y1 = get_subthreshold(Y1, spk10, scr.sub_dilate, p10.smoothMs);
        Y2 = get_subthreshold(Y2, spk10, scr.sub_dilate, p10.smoothMs);
    elseif p10.smoothMs > 1
        Y1 = movmean(Y1, p10.smoothMs, 2);  Y2 = movmean(Y2, p10.smoothMs, 2);
    end
    prep10{q} = {Y1, Y2};
end
srcNm10 = {'uncorrected', 'bleach-dr corrected'};

% painting info for the PC maps, and the ROI order used elsewhere
[wMax10, aMax10] = max(ftT, [], 3);
PT10 = struct('aMax', aMax10, 'supp', wMax10 > 0, 'siz', [nrow10 ncol10]);
if exist('ord8','var') && numel(ord8) == size(ftT,3), ordC10 = ord8; else, ordC10 = 1:size(ftT,3); end
if p10.useSub, subTag10 = 'subthreshold'; else, subTag10 = 'full trace'; end
if p10.smoothMs > 1, smTag10 = sprintf('%d ms', p10.smoothMs); else, smTag10 = 'no smoothing'; end

% ---------------------------------------------------------------------------
% 10b'. RUN -- uncorrected
% ---------------------------------------------------------------------------
E10 = cell(1,2);
E10{1} = covEstimators3(prep10{1}{1}, prep10{1}{2}, good10);
plotCovEstimators3(E10{1}, PT10, struct( ...
    'figName','R3-10-estimators-uncorrected', 'nPC',p10.nPC, 'ord',ordC10, ...
    'corrClim',p10.corrClim, 'paintPrc',p10.paintPrc, ...
    'titleStr', sprintf('f=%d  UNCORRECTED, %s, %s -- three covariance estimators', ...
                        f, subTag10, smTag10)));
report_cov3(E10{1}, srcNm10{1});

%% ################################################################
%  (10b) THE SAME THREE ESTIMATORS ON THE BLEACH-DR CORRECTED TRACE
%
%  Identical ROIs, frames, normalisation and code path as section (10) -- only
%  the movie differs -- so the two figures are directly comparable panel for
%  panel. At the end the two runs are compared to each other:
%    * how much variance the correction removed, per estimator
%    * whether the correlation STRUCTURE survived it
%  Estimator (ii) is the one to trust for that comparison: it is noise-free, so
%  it is not confounded by PMD having denoised the corrected movie.
% ################################################################
assert(exist('E10','var')==1 && ~isempty(E10{1}), '[R3-10b] run section (10) first');

E10{2} = covEstimators3(prep10{2}{1}, prep10{2}{2}, good10);
plotCovEstimators3(E10{2}, PT10, struct( ...
    'figName','R3-10b-estimators-bleachdr', 'nPC',p10.nPC, 'ord',ordC10, ...
    'corrClim',p10.corrClim, 'paintPrc',p10.paintPrc, ...
    'titleStr', sprintf('f=%d  BLEACH-DR CORRECTED, %s, %s -- three covariance estimators', ...
                        f, subTag10, smTag10)));
report_cov3(E10{2}, srcNm10{2});

% ---------------------------------------------------------------------------
% 10b''. UNCORRECTED vs CORRECTED, estimator by estimator
% ---------------------------------------------------------------------------
nR10  = E10{1}.nROI;
tri10 = triu(true(nR10),1);
fprintf('\n[R3-10b] ===== uncorrected vs bleach-dr, per estimator =====\n');
fprintf('   %-34s %10s %10s %12s\n', 'estimator', 'var kept', 'r(corr)', 'mean|dCorr|');
for e = 1:3
    vk = mean(diag(E10{2}.Cov{e})) / mean(diag(E10{1}.Cov{e}));
    u = E10{1}.Corr{e}(tri10);  v = E10{2}.Corr{e}(tri10);
    ok = isfinite(u) & isfinite(v);
    fprintf('   %-34s %10.3f %10.3f %12.3f\n', E10{1}.names{e}, vk, ...
            corr(u(ok), v(ok)), mean(abs(u(ok)-v(ok))));
end
fprintf(['\n   var kept   : mean diagonal after / before. Read it from estimator (ii):\n' ...
         '                (i) is confounded because PMD removed shot noise as well.\n' ...
         '   r(corr)    : off-diagonal correlation between the two matrices. Near 1 means\n' ...
         '                the correlation STRUCTURE survived; well below 1 means the\n' ...
         '                correction reshaped it, not just scaled it down.\n']);

% how far apart are the two noise-free estimators, within each movie?
fprintf('\n[R3-10b] noise-free estimators (ii) vs (iii), within each movie:\n');
for q = 1:2
    u = E10{q}.Corr{2}(tri10);  v = E10{q}.Corr{3}(tri10);
    ok = isfinite(u) & isfinite(v);
    fprintf('   %-22s r = %.3f | diag ratio (iii)/(ii) = %.3f  -> %.1f%% of the real\n', ...
            srcNm10{q}, corr(u(ok), v(ok)), ...
            mean(diag(E10{q}.Cov{3}))/mean(diag(E10{q}.Cov{2})), ...
            100*(1 - mean(diag(E10{q}.Cov{3}))/mean(diag(E10{q}.Cov{2}))));
    fprintf('   %-22s variance decorrelates within one frame\n', '');
end

% ===================== local functions (must stay last) =====================
function report_cov3(E, srcName)
% Console summary for one covEstimators3 result.
fprintf('\n[R3-10] ===== %s =====\n', srcName);
fprintf('   %d ROIs, %d analysed frames, %d consecutive pairs\n', E.nROI, E.nGood, E.nPair);
fprintf('   %-34s %10s %10s %10s %10s\n', 'estimator', 'mean diag', '#pos eig', 'PC1 %', 'PC1-5 %');
for e = 1:3
    D = E.Eval{e};  pos = D > 0;  tot = sum(D(pos));
    fprintf('   %-34s %10.4g %10d %10.1f %10.1f\n', E.names{e}, mean(diag(E.Cov{e})), ...
            E.nPosEig(e), 100*D(1)/tot, 100*sum(D(1:min(5,end)))/tot);
    if E.nBadDiag(e) > 0
        fprintf('        (%d ROI(s) had non-positive estimated variance -> NaN in the matrix)\n', E.nBadDiag(e));
    end
end
tri = triu(true(E.nROI),1);
fprintf('   off-diagonal agreement:');
lbl = {'i','ii','iii'};
for a = 1:3
    for b = a+1:3
        u = E.Corr{a}(tri);  v = E.Corr{b}(tri);  ok = isfinite(u) & isfinite(v);
        fprintf('  r(%s,%s)=%.3f', lbl{a}, lbl{b}, corr(u(ok), v(ok)));
    end
end
fprintf('\n   diag ratio (ii)/(i) = %.3f  -> %.1f%% of the ordinary variance was shot noise\n', ...
        mean(diag(E.Cov{2}))/mean(diag(E.Cov{1})), ...
        100*(1 - mean(diag(E.Cov{2}))/mean(diag(E.Cov{1}))));
end

%% ################################################################
%  (10h) THE TWO REVIEWER FIGURES
%
%  Everything here is ASSEMBLED from what sections (7)/(8) already computed --
%  nothing is recalculated -- so these publication figures and the diagnostic
%  ones cannot disagree.
%
%  FIGURE 1  "the fast component is vascular"
%     i    neuron image with the two example ROIs marked
%     ii   their traces: raw + subthreshold, and both difference traces overlaid
%     iii  high-pass variance image with the bvMask contour
%     iv   box plots of raw vs high-pass variance, split by vessel overlap
%
%  FIGURE 2  "the high PCs collapse after smoothing / correction"
%     i    one example trace before and after bleach-dr
%     ii   PC variance, log-log        (raw / 20 ms smoothed / corrected)
%     iii  cumulative variance         (same three)
%     iv   eigenvector maps of PC1-3 for raw and corrected, titled with the
%          variance each explains
% ################################################################
for reqH = {'TR','dF','dFhp','hp_raw_img','hp_cor_img','bvBW8','bvFrac8','p8', ...
            'evList','UkList','dpix8','dMask8','tShow8','tms8','in7C','tr7C','F0full'}
    assert(exist(reqH{1},'var')==1, '[R3-10h] %s is missing -- run section (8) first', reqH{1});
end
for reqH = {'in7_1','in7_2','tr7_1','tr7_2','F0o10','F0e10','spk10','good10','PT10','ordC10','p10'}
    assert(exist(reqH{1},'var')==1, '[R3-10h] %s is missing -- run section (10) first', reqH{1});
end

pH = struct();
pH.exROI      = [21];        % ROI for the before/after trace in figure 2; [] = p8.exROI(1)
pH.pcaOn      = 'ftprnt';  % what figure 2's PCA runs on:
                           %  'ftprnt' : the Result.ftprnt ROI traces, exactly as section (10)
                           %             builds them, via covEstimators3. Interpretable and
                           %             noise-free, but the spectrum is only nROI long.
                           %  'pixels' : section (8)'s dendrite-mask pixel PCA (~21k px, so a
                           %             far longer spectrum, but one point per pixel mode).
pH.pcaEst     = 3;         % which covEstimators3 estimator when pcaOn = 'ftprnt':
                           %  1 = X(t)X(t) same-time (biased by shot noise)
                           %  2 = X_odd X_even       (noise-free IF the checkerboard halves
                           %                          really are independent -- see the
                           %                          sub-pixel interpolation caveat)
                           %  3 = X(t)X(t+1) lag-1   (noise-free, immune to that caveat)
pH.varNorm    = 'self';     % how panels F and G scale the eigenvalues:
                           %  'raw'  : all three divided by the RAW total. Each curve then
                           %           ENDS at the fraction of raw variance it retains, and
                           %           its curvature is still the concentration -- both
                           %           readable at once.
                           %  'self' : each divided by its OWN total (every curve ends at
                           %           100%%). Shows concentration only, and HIDES the fact
                           %           that the corrected movie has far less total variance --
                           %           which is the result being claimed.
                           %  'abs'  : no normalisation, in (%%dF/F)^2.
pH.nPCmap     = 3;         % PCs mapped in figure 2; also the number regressed out for
                           % the residual image in the 4th column
pH.mapConds   = [1 2 3];   % which of the three conditions get a PC-map row in figure 2:
                           %  1 = raw, 2 = Wsm smoothed, 3 = bleach-dr corrected.
                           % One row per entry, in this order. [1 3] gives the
                           % two-row version (raw vs corrected only).
pH.residROI    = 'ftprnt';    % what the residual is measured on, then painted on the dendrite:
                           %  'ftprnt' : the same Result.ftprnt ROIs the PCA used
                           %  'sub'    : the ~30 px sub-ROIs -- same idea, finer spatial map
pH.residBlk    = 2000;     % frames per block when a movie has to be streamed for 'sub'
pH.residClim   = [];       % colour limits SHARED by every residual panel, so the drop is
                           % readable off the picture. [] = 0.1/99.9 percentile over ALL
                           % of the mapped conditions together; set [lo hi] to fix by hand.
pH.residContour = true;    % outline bvMask on the residual panels, so any surviving
                           % variance can be read as on- or off-vessel
pH.residEst    = 'oddeven'; % how the residual is measured:
                           %  'oddeven' : covariance of the odd x even checkerboard channels.
                           %          Their shot noise is independent whatever temporal filter
                           %          has been applied, so this is the only one of the three
                           %          that stays unbiased across the smoothed conditions.
                           %  'lag1'  : X(t)X(t+1). Unbiased ONLY on unsmoothed traces. After a
                           %          w-point moving average white noise acquires a lag-1
                           %          autocovariance of sigma^2*(w-1)/w^2 against a variance of
                           %          sigma^2/w, i.e. a fraction (w-1)/w of the noise SURVIVES
                           %          -- 0% at w=1 but 95% at w=20. That bias is differential
                           %          across exactly the conditions this figure compares, which
                           %          is why it is no longer the default. Diagnostic only.
                           %  'var'   : ordinary variance. Dominated by shot noise once the top
                           %          components are gone (~78x inflation), so both movies
                           %          then look the same -- diagnostic only.
pH.hfClim     = [0 35000];        % [] = 1st/99.5th percentile of the raw HF image
pH.corrTrace  = 'robust';  % units of the before/after trace in figure 2:
                           %  'robust' : the pipeline's robust dF/F. The ScaleFactor is fit on
                           %             Result.normTraces, which exists only for the
                           %             UNCORRECTED movie, so it cannot be applied to the
                           %             corrected trace directly. Instead the movie trace is
                           %             CALIBRATED onto the robust dF/F on the uncorrected
                           %             data, and the SAME per-ROI scalar is applied to the
                           %             corrected one -- so both sides are in robust dF/F and
                           %             stay comparable. The calibration r is printed; if it
                           %             is poor the code falls back to F0 and says so.
                           %  'F0'     : divide both by the mean-image F0 through the same
                           %             footprint (always available, no calibration step)
pH.showMs     = min(p8.showMs, numel(tShow8));

% ===========================================================================
% FIGURE 1
% ===========================================================================
get_figure('R3-10h-fig1');
tlH1 = tiledlayout(5, 1, 'TileSpacing','compact','Padding','compact');

% ---- i. neuron + example ROIs --------------------------------------------
ax = nexttile(1);
imagesc(ax, Pr.ref_im, prctile(Pr.ref_im(:), [20 99.8]));
colormap(ax, gray); axis(ax,'image','off'); hold(ax,'on');
contour(ax, bvBW8, [0.5 0.5], 'Color', p8.vesCol, 'LineWidth', 1.2);
coordH = get_coord(TR.ft);
for k = 1:2
    B = bwboundaries(TR.ft(:,:,p8.exROI(k)) > 0);
    for b = 1:numel(B), plot(ax, B{b}(:,2), B{b}(:,1), '-', 'Color', p8.exColr(k,:), 'LineWidth', 2.5); end
    text(ax, coordH(p8.exROI(k),1)+5, coordH(p8.exROI(k),2)-6, sprintf('#%d', p8.exROI(k)), ...
         'Color', p8.exColr(k,:), 'FontSize', 11, 'FontWeight','bold');
end
sbH = 100 / Pixelsize(f);
plot(ax, size(Pr.ref_im,2)-12-[sbH 0], (size(Pr.ref_im,1)-10)*[1 1], 'w-', 'LineWidth', 3);
title(ax, sprintf('f=%d   red = blood vessel,  20 \\mum bar', f), 'FontSize', 10);

% ---- ii. traces: ON raw+subth / OFF raw+subth / both diffs overlaid -------
ax = nexttile(2, [2 1]); hold(ax,'on');
tSh = tShow8(1:pH.showMs);  tMs = tms8(1:pH.showMs);
YbH = cell(1,2);
for k = 1:2
    YbH{k} = [TR.dF(p8.exROI(k), tSh); TR.dF20(p8.exROI(k), tSh); TR.dFhp(p8.exROI(k), tSh)];
    YbH{k}(1:2,:) = YbH{k}(1:2,:) - mean(YbH{k}(1,:));
    YbH{k}(3,:)   = YbH{k}(3,:)   - mean(YbH{k}(3,:));
end
if p8.showClip, kAmpH = 2; else, kAmpH = 1; end            % 2 = subthreshold, 1 = raw
ampH = max([cellfun(@(Y) max(Y(kAmpH,:))-min(Y(kAmpH,:)), YbH), ...
            cellfun(@(Y) max(Y(3,:))-min(Y(3,:)), YbH)]);
dRowH = -1.25*ampH;  offH = [0 dRowH 2*dRowH];
vesTagH = {'ON vessel','OFF vessel'};
for k = 1:2
    plot(ax, tMs, YbH{k}(1,:) + offH(k), '-', 'Color', [.65 .65 .65], 'LineWidth', 0.4);
    plot(ax, tMs, YbH{k}(2,:) + offH(k), '-', 'Color', p8.exColr(k,:), 'LineWidth', 1.6);
    text(ax, 6, offH(k) + 0.55*ampH, sprintf('%s  #%d (%.0f%% overlap)   raw + subthreshold (%d ms)', ...
         vesTagH{k}, p8.exROI(k), 100*TR.bvFrac(p8.exROI(k)), p8.Wsm), ...
         'Color', p8.exColr(k,:), 'FontSize', 8);
end
for k = 1:2
    plot(ax, tMs, YbH{k}(3,:) + offH(3), '-', 'Color', p8.exColr(k,:), 'LineWidth', 0.6);
end
text(ax, 6, offH(3) + 0.55*ampH, sprintf('diff (raw - %d ms), both ROIs overlaid', p8.Wsm), ...
     'Color', [.2 .2 .2], 'FontSize', 8);
sbV = round(ampH/3, 1, 'significant');
xsb = tMs(end) - 0.02*tMs(end);  ysb = offH(3) - 0.55*ampH;
plot(ax, xsb*[1 1], ysb + [0 sbV], 'k-', 'LineWidth', 2.5);
text(ax, xsb - 0.012*tMs(end), ysb + sbV/2, sprintf('%g%% %s', sbV, TR.norm), ...
     'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',8);
set(ax,'YTick',[],'YColor','none'); box(ax,'off');
xlim(ax,[0 tMs(end)]); ylim(ax,[ysb - 0.25*ampH, 0.75*ampH]); xlabel(ax,'time (ms)');

% ---- iii. high-pass variance image with the vessel contour ---------------
ax = nexttile(4);
if isempty(pH.hfClim), clH = prctile(hp_raw_img(:), [1 99.5]); else, clH = pH.hfClim; end
imagesc(ax, hp_raw_img, clH); axis(ax,'image','off'); colormap(ax, turbo); hold(ax,'on');
contour(ax, bvBW8, [0.5 0.5], 'w', 'LineWidth', 1);
colorbar(ax);
title(ax, sprintf('high-pass variance  var(raw - %d ms),  white = bvMask', p8.Wsm), 'FontSize', 9);

% ---- iv. box plots -------------------------------------------------------
ax = nexttile(5); hold(ax,'on');
vR = var(dF  (:,gd8), 0, 2, 'omitnan');
vH = var(dFhp(:,gd8), 0, 2, 'omitnan');
onH = bvFrac8(:) > p8.bvFracMark;
cNmH = {'raw', sprintf('high-pass (raw-%dms)', p8.Wsm)};
xgH = categorical(repelem((1:2)', numel(vR), 1), 1:2, cNmH);
cgH = categorical(repmat(onH, 2, 1), [false true], ...
                  {sprintf('no vessel overlap (n=%d)', nnz(~onH)), ...
                   sprintf('overlaps vessel (n=%d)',   nnz(onH))});
bhH = boxchart(ax, xgH, [vR; vH], 'GroupByColor', cgH, 'MarkerStyle','none');
bhH(1).BoxFaceColor = [.55 .55 .55];
if numel(bhH) == 2, bhH(2).BoxFaceColor = p8.vesCol; end
wLoH = inf; wHiH = -inf;                                  % axis to the whiskers, not outliers
for cc = 1:2
    v0 = vR; if cc == 2, v0 = vH; end
    for g = [false true]
        v = v0(onH == g);  v = v(~isnan(v));
        if numel(v) < 2, continue; end
        q = quantile(v,[.25 .75]);  iq = diff(q);
        wHiH = max(wHiH, max(v(v <= q(2)+1.5*iq)));
        wLoH = min(wLoH, min(v(v >= q(1)-1.5*iq)));
    end
end
if isfinite(wLoH) && wHiH > wLoH
    padH = 0.05*(wHiH-wLoH);  ylim(ax, [max(0,wLoH-padH), wHiH+padH]);
end
grid(ax,'on'); ylabel(ax, sprintf('variance  ((%%\\DeltaF/F_0)^2)'));
legend(ax,'Location','best','Box','off','FontSize',7);

% ---- statistics for this panel -------------------------------------------
% Wilcoxon rank-sum (Mann-Whitney U), two-sided, on-vessel vs off-vessel. The
% ROI variances are strongly right-skewed, so a rank test is the honest choice.
% The THIRD test is the one that carries the argument: it asks whether vessel
% ROIs put a larger FRACTION of their variance in the high-pass band, which is
% not answered by either column on its own -- an ROI could simply be brighter.
statH = struct('name',{},'p',{},'ratio',{},'medOn',{},'medOff',{});
for cc = 1:2
    v0 = vR;  if cc == 2, v0 = vH; end
    a = v0(onH);  b = v0(~onH);
    a = a(isfinite(a));  b = b(isfinite(b));
    statH(cc).name  = cNmH{cc};
    statH(cc).p     = ranksum(a, b);
    statH(cc).medOn = median(a);  statH(cc).medOff = median(b);
    statH(cc).ratio = median(a)/median(b);
end
fr = vH ./ vR;                                    % fraction of variance in the high-pass band
fa = fr(onH & isfinite(fr));  fb = fr(~onH & isfinite(fr));
statH(3).name  = 'high-pass FRACTION (hp/raw)';
statH(3).p     = ranksum(fa, fb);
statH(3).medOn = median(fa);  statH(3).medOff = median(fb);
statH(3).ratio = median(fa)/median(fb);

fprintf('\n[R3-10h] figure 1 box panel -- %s, %d ROIs (%d on vessel, %d off)\n', ...
        sprintf('sub-ROI %d px', p8.subPx), numel(vR), nnz(onH), nnz(~onH));
fprintf('   %-30s %10s %10s %8s %12s\n', 'quantity', 'med ON', 'med OFF', 'ratio', 'ranksum p');
for q = 1:numel(statH)
    if statH(q).p < 1e-4, pStr = sprintf('%.1e', statH(q).p);
    else,                 pStr = sprintf('%.4f',  statH(q).p);  end
    fprintf('   %-30s %10.4g %10.4g %8.2f %12s\n', statH(q).name, statH(q).medOn, ...
            statH(q).medOff, statH(q).ratio, pStr);
end
fprintf(['   (two-sided Wilcoxon rank-sum; the third row is the one that separates\n' ...
         '    "vessel ROIs are brighter" from "vessel ROIs are relatively FASTER")\n']);

if statH(1).p < 1e-4, p1 = sprintf('p=%.0e', statH(1).p); else, p1 = sprintf('p=%.3f', statH(1).p); end
if statH(2).p < 1e-4, p2 = sprintf('p=%.0e', statH(2).p); else, p2 = sprintf('p=%.3f', statH(2).p); end
title(ax, sprintf('per sub-ROI (%d px), %d ROIs   |   raw x%.2f, %s   high-pass x%.2f, %s', ...
      p8.subPx, numel(vR), statH(1).ratio, p1, statH(2).ratio, p2), 'FontSize', 8);
title(tlH1, sprintf('f=%d  the fast subthreshold component is vascular', f));

% ===========================================================================
% FIGURE 2
% ===========================================================================
if isempty(pH.exROI), exH = p8.exROI(1); else, exH = pH.exROI; end
% Both sides must share ONE scaling or the comparison is meaningless. Whatever
% that scaling is, it cancels out of the before/after difference -- which is why
% calibrating onto the robust dF/F is legitimate even though the ScaleFactor
% itself was never fit on the corrected movie.
normTagH = '\DeltaF/F_0';
kCal = 100 / F0full(exH);                                  % the F0 route, as a fallback
if strcmp(pH.corrTrace, 'robust')
    robTr = 100 * Pr.normTraces(Pr.sortDist(exH), b7.frames) ./ Pr.scalingfactor(Pr.sortDist(exH));
    gCal  = good10 & isfinite(robTr) & isfinite(in7C(exH,:));
    aC = in7C(exH,gCal) - mean(in7C(exH,gCal));
    rC = robTr(gCal)    - mean(robTr(gCal));
    rCal = corr(aC(:), rC(:));
    if abs(rCal) > 0.5
        kCal = (aC*rC') / (aC*aC');                        % counts -> robust dF/F, this ROI
        normTagH = 'robust \DeltaF/F';
        fprintf('[R3-10h] trace calibrated to robust dF/F on ROI #%d: r = %.3f, scale = %.4g\n', ...
                exH, rCal, kCal);
    else
        warning(['[R3-10h] the movie trace and the pipeline robust dF/F agree only at r = %.2f ' ...
                 'on ROI #%d -- calibration is not trustworthy, falling back to dF/F0'], rCal, exH);
    end
end
dFbefore = kCal * in7C(exH,:);
dFafter  = kCal * tr7C(exH,:);
dFbefore = dFbefore - movmedian(dFbefore, p8.baseWin, 2, 'omitnan');
dFafter  = dFafter  - movmedian(dFafter,  p8.baseWin, 2, 'omitnan');

% ---- bAP PRESERVATION: the number the panel-E claim rests on -------------
% Recomputed here from in7C / tr7C rather than reused from section (7), so it
% cannot drift from what this figure shows. The bAP is fast (~2 ms) and
% branch-local, so if it survives the correction, branch-local fast signal did.
% The ratio is scale-free, which is why the raw counts traces can be used
% directly with no normalisation.
spkH  = full(Result.spike(1, b7.frames)) > 0;
bapI  = get_STA(in7C, spkH, scr.bAP_pre, scr.bAP_post);
bapI  = bapI - prctile(bapI, scr.bAP_baseprc, 2);   ampI = max(bapI, [], 2);
bapC  = get_STA(tr7C, spkH, scr.bAP_pre, scr.bAP_post);
bapC  = bapC - prctile(bapC, scr.bAP_baseprc, 2);   ampC = max(bapC, [], 2);
% Two measures, because they fail differently:
%   amplitude ratio  max(bAP_corr)/max(bAP_raw) -- "did the bAP SHRINK?"
%   waveform r       corr of the two STA waveforms -- "did its SHAPE change?"
% Correlation alone would be misleading: it is SCALE-INVARIANT and reads 1.000
% even for a bAP attenuated to 40%, i.e. it is blind to the failure this claim
% is about. Amplitude alone misses distortion and latency shifts. Report both.
bapKeep = 100 * ampC ./ ampI;                              % amplitude, %
bapCorr = nan(size(bapKeep));                              % waveform shape
for n = 1:size(bapI,1)
    bapCorr(n) = corr(bapI(n,:)', bapC(n,:)');
end
labH  = Pr.labelvec(:);  clsNmH = {'Basal','Soma','Trunk','Oblique','Distal'};

fprintf('\n[R3-10h] bAP PRESERVATION after blood correction (%d spikes, %d ROIs)\n', ...
        nnz(spkH), numel(bapKeep));
fprintf('   AMPLITUDE kept  max(bAP_corr)/max(bAP_raw):\n');
fprintf('      median %.1f%%   mean %.1f%%   worst branch %.1f%% (ROI #%d)   IQR %.1f-%.1f%%\n', ...
        median(bapKeep,'omitnan'), mean(bapKeep,'omitnan'), min(bapKeep), ...
        find(bapKeep == min(bapKeep), 1), prctile(bapKeep,25), prctile(bapKeep,75));
fprintf('      (in VARIANCE units that is %.1f%% median, %.1f%% worst -- amplitude squared)\n', ...
        median(bapKeep,'omitnan')^2/100, min(bapKeep)^2/100);
fprintf('   WAVEFORM r  corr(bAP_raw, bAP_corr):\n');
fprintf('      median %.3f   worst %.3f\n', median(bapCorr,'omitnan'), min(bapCorr));
fprintf('   per dendrite class:\n');
fprintf('      %-8s %9s %9s %9s\n', 'class', 'amp med', 'amp worst', 'r med');
for c = unique(labH(~isnan(labH)))'
    v = bapKeep(labH == c);  q = bapCorr(labH == c);
    fprintf('      %-8s %8.1f%% %8.1f%% %9.3f   (%d ROIs)\n', ...
            clsNmH{c}, median(v,'omitnan'), min(v), median(q,'omitnan'), numel(v));
end
fprintf('   example ROI #%d (panel E): amplitude %.1f%%, waveform r %.3f\n', ...
        exH, bapKeep(exH), bapCorr(exH));
fprintf(['   -> for the legend: "the correction preserves %.0f%% of the bAP AMPLITUDE on\n' ...
         '      average (median; no less than %.0f%% on any branch), with the bAP waveform\n' ...
         '      unchanged (median r = %.2f)"\n'], ...
        median(bapKeep,'omitnan'), min(bapKeep), median(bapCorr,'omitnan'));
if nnz(spkH) < 20
    warning('[R3-10h] only %d spikes -- the bAP ratio is noisy; widen b7.frames before quoting it', nnz(spkH));
end

get_figure('R3-10h-fig2');
% rows 1-2 are the trace and the two spectrum panels; every mapped condition
% then gets a row of PC images + its residual.
assert(~isempty(pH.mapConds) && all(ismember(pH.mapConds, 1:3)), ...
       '[R3-10h] pH.mapConds must be a non-empty subset of [1 2 3]');
nRowH2 = 2 + numel(pH.mapConds);
tlH2 = tiledlayout(nRowH2, 8, 'TileSpacing','compact','Padding','compact');

% ---- i. example trace before / after --------------------------------------
% Before and after OVERLAID on one baseline, so what the correction removed is
% the visible gap between the two rather than something the reader has to
% compare across an offset.
ax = nexttile(1, [1 8]); hold(ax,'on');
yB = dFbefore(tSh) - mean(dFbefore(tSh));
yA = dFafter(tSh)  - mean(dFafter(tSh));
ampB = max(yB) - min(yB);
hB = plot(ax, tMs, yB, '-', 'Color', [.62 .62 .62], 'LineWidth', 0.5, ...
          'DisplayName', sprintf('before  (ROI #%d)', exH));
hA = plot(ax, tMs, yA, '-', 'Color', [0 .45 .74], 'LineWidth', 0.9, ...
          'DisplayName', 'after bleach-dr');
sbB = round(ampB/4, 1, 'significant');
yLo = min([yB yA]);
plot(ax, (tMs(end)-25)*[1 1], yLo + [0 sbB], 'k-', 'LineWidth', 2.5);
text(ax, tMs(end)-35, yLo + sbB/2, sprintf('%g%% %s', sbB, normTagH), ...
     'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',8);
legend(ax, [hB hA], 'Location','northwest', 'Box','off', 'FontSize',8);
set(ax,'YTick',[],'YColor','none'); box(ax,'off');
xlim(ax,[0 tMs(end)]); ylim(ax, [yLo - 0.08*ampB, max([yB yA]) + 0.22*ampB]);
xlabel(ax,'time (ms)');
title(ax, sprintf(['example trace before and after blood correction  ' ...
                   '(variance %.3g -> %.3g, %.0f%% removed;  bAP amplitude kept: median ' ...
                   '%.0f%%, worst branch %.0f%%;  waveform r = %.2f)'], ...
                  var(yB), var(yA), 100*(1-var(yA)/var(yB)), ...
                  median(bapKeep,'omitnan'), min(bapKeep), median(bapCorr,'omitnan')), 'FontSize', 9);

% ---- the three spectra --------------------------------------------------
pcaNmH = {'raw', sprintf('%d ms smoothed', p8.Wsm), 'bleach-dr corrected'};
cmH    = [0 0 0; 0.85 0.33 0.10; 0.00 0.45 0.74];
if strcmp(pH.pcaOn, 'ftprnt')
    % Same ROI traces, normalisation and frames as section (10); the three
    % conditions differ only in movie and smoothing window.
    srcH2 = { {in7_1,in7_2}, 1        ; ...      % raw, no smoothing
              {in7_1,in7_2}, p8.Wsm   ; ...      % raw, Wsm smoothed
              {tr7_1,tr7_2}, 1        };         % bleach-dr corrected
    evH = cell(1,3);  UH = cell(1,3);  EH = cell(1,3);  lamAll = cell(1,3);
    for i = 1:3
        Y1 = 100 * srcH2{i,1}{1} ./ F0o10;
        Y2 = 100 * srcH2{i,1}{2} ./ F0e10;
        Y1 = get_subthreshold(Y1, spk10, scr.sub_dilate, srcH2{i,2});
        Y2 = get_subthreshold(Y2, spk10, scr.sub_dilate, srcH2{i,2});
        EH{i} = covEstimators3(Y1, Y2, good10);
        lamAll{i} = EH{i}.Eval{pH.pcaEst};
        UH{i}     = EH{i}.Evec{pH.pcaEst};
    end
    % scale AFTER all three exist, so 'raw' can use the raw total as the reference
    totH = cellfun(@(L) sum(L(L>0)), lamAll);
    switch pH.varNorm
        case 'self', den = totH;                 vUnit = '% of own total';
        case 'raw',  den = repmat(totH(1),1,3);  vUnit = '% of RAW total';
        case 'abs',  den = ones(1,3);            vUnit = '(%\DeltaF/F)^2';
        otherwise,   error('pH.varNorm must be ''raw'', ''self'' or ''abs''');
    end
    for i = 1:3, evH{i} = lamAll{i} / den(i); end
    fprintf('[R3-10h] total variance: raw %.4g | %d ms %.4g (%.0f%% of raw) | corrected %.4g (%.0f%% of raw)\n', ...
            totH(1), p8.Wsm, totH(2), 100*totH(2)/totH(1), totH(3), 100*totH(3)/totH(1));
    fprintf('         panels F/G scaled as ''%s'' (%s)\n', pH.varNorm, vUnit);
    estNmH = EH{1}.names{pH.pcaEst};
    pcaTagH = sprintf('%s ROI traces, %s', 'Result.ftprnt', estNmH);
    fprintf('[R3-10h] figure 2 PCA on %d Result.ftprnt ROIs, %s\n', EH{1}.nROI, estNmH);
    if EH{1}.nROI < 100
        fprintf(['         NOTE the spectrum is only %d points long (one per ROI), against\n' ...
                 '         ~%d for the pixel version. The "high PCs" argument is much weaker\n' ...
                 '         over %d dimensions -- set pH.pcaOn = ''pixels'' for the long one.\n'], ...
                 EH{1}.nROI, numel(dpix8), EH{1}.nROI);
    end
else
    evH = evList;  UH = UkList;  vUnit = '% of own total';
    pcaTagH = sprintf('SWC dendrite mask, %d px, %s covariance', numel(dpix8), p8.pcaCov);
end
ax = nexttile(9, [1 4]); hold(ax,'on');
hH = gobjects(1,3);
for i = [3 2 1]
    ev = evH{i};  kp = find(ev > 0);      % noise-free eigenvalues go negative past the signal
    scF = 100;  if strcmp(pH.varNorm,'abs'), scF = 1; end
    hH(i) = plot(ax, kp, scF*ev(kp), '-', 'Color', cmH(i,:), ...
                 'LineWidth', 1.5 + 0.6*(i==1), 'DisplayName', pcaNmH{i});
end
set(ax,'XScale','log','YScale','log'); grid(ax,'on');
xlim(ax,[1 max(cellfun(@numel, evH))]);
xlabel(ax,'PC #'); ylabel(ax, sprintf('variance  %s', vUnit));
legend(ax, hH, pcaNmH, 'Location','southwest','Box','off','FontSize',7);
title(ax,'PC variance spectrum','FontSize',9);

ax = nexttile(13, [1 4]); hold(ax,'on');
sc = 100;  if strcmp(pH.varNorm,'abs'), sc = 1; end
for i = 1:3
    ev = evH{i};  kp = ev > 0;               % cumulate the positive spectrum only
    plot(ax, find(kp), sc*cumsum(ev(kp)), '-', 'Color', cmH(i,:), 'LineWidth', 1.5, ...
         'DisplayName', pcaNmH{i});
end
set(ax,'XScale','log'); grid(ax,'on');
xlim(ax,[1 max(cellfun(@numel, evH))]);
if strcmp(pH.varNorm,'self'), ylim(ax,[0 100]); end   % only 'self' is bounded at 100
xlabel(ax,'PC #'); ylabel(ax, sprintf('cumulative variance  %s', vUnit));
legend(ax,'Location','southeast','Box','off','FontSize',7);
title(ax,'cumulative variance','FontSize',9);

% ---- iv. PC1-3 maps, then the residual after regressing them out ---------
% The three components become TIME COURSES (from the ftprnt traces the PCA ran
% on) and those are regressed out of every ROI of the chosen set. Only the time
% courses matter, so the regression target need not live in the same space as
% the eigenvectors -- which is what lets 'sub' give a finer map.
srcH = pH.mapConds(:)';                        % one row per mapped condition,
rowH = 17 + 8*(0:numel(srcH)-1);               % starting on row 3 of the layout
nMapH  = numel(srcH);
residH = cell(1,nMapH);  residFrac = nan(1,nMapH);

switch pH.residROI
    case 'ftprnt'
        ftR = Pr.Ftprnts(:,:,Pr.sortDist);
        setTagH = 'Result.ftprnt';
    case 'sub'
        assert(exist('ftprintsmallest','var')==1, '[R3-10h] ftprintsmallest missing -- run section (8)');
        ftR = double(ftprintsmallest);
        ftR = ftR ./ max(sum(sum(ftR,1),2), eps);               % unit-sum, as in section (8)
        setTagH = sprintf('sub-ROI (%d px)', p8.subPx);
    otherwise
        error('pH.residROI must be ''ftprnt'' or ''sub''');
end
[wR, aR] = max(ftR, [], 3);  suppR = wR > 0;                    % pixel -> owning ROI
% Split each footprint by checkerboard parity. The odd and even halves see the
% same voltage but INDEPENDENT shot noise, and that independence survives any
% temporal filtering -- unlike the lag-1 trick, which a moving average breaks.
% The price is that a very small sub-ROI can contain no pixels of one parity, so
% those ROIs are dropped rather than fatal.
ftRv = tovec(ftR);
ftRo = tovec(ftR .*  oddPix10);   ftRe = tovec(ftR .* ~oddPix10);
F0R  = (tovec(F0img10)' * ftRv)';
F0Ro = (tovec(F0img10)' * ftRo)'; F0Re = (tovec(F0img10)' * ftRe)';
okPar = F0Ro > 0 & F0Re > 0;                                    % usable on both channels
if ~any(okPar)
    error('[R3-10h] no %s ROI has pixels of both checkerboard parities', setTagH);
end
if ~all(okPar)
    warning(['[R3-10h] %d of %d %s ROI(s) lack one checkerboard parity (or have F0 <= 0 ' ...
             'on it) and are left blank in the residual map'], nnz(~okPar), numel(okPar), setTagH);
end
F0Ro(~okPar) = NaN;  F0Re(~okPar) = NaN;
fprintf('[R3-10h] residual on %s (%d of %d ROIs usable), regressing out PC1-%d\n', ...
        setTagH, nnz(okPar), size(ftR,3), pH.nPCmap);

for r = 1:nMapH
    i  = srcH(r);
    U3 = UH{i}(:, 1:pH.nPCmap);

    % --- component time courses, from the traces the PCA used --------------
    Y1p = get_subthreshold(100*srcH2{i,1}{1}./F0o10, spk10, scr.sub_dilate, srcH2{i,2});
    Y2p = get_subthreshold(100*srcH2{i,1}{2}./F0e10, spk10, scr.sub_dilate, srcH2{i,2});
    Ccomp = U3' * ((Y1p + Y2p)/2);                              % nPC x nT

    % --- the traces the components get regressed OUT of, one per parity ----
    if strcmp(pH.residROI,'ftprnt')
        Zo = 100 * srcH2{i,1}{1} ./ F0Ro;                       % odd  channel
        Ze = 100 * srcH2{i,1}{2} ./ F0Re;                       % even channel
    else
        Zo = zeros(size(ftR,3), nT8);  Ze = zeros(size(ftR,3), nT8);
        for t0 = 1:pH.residBlk:nT8                              % sub-ROI traces are not
            t1 = min(t0 + pH.residBlk - 1, nT8);                % precomputed for either movie
            if i == 3
                Mv = permute(h5read(out7,'/corrected',[1 1 t0], ...
                                    [ncol8 nrow8 t1-t0+1]), [2 1 3]);
            else
                Mv = single(-mov_res(:,:, loc8(t0:t1)));
            end
            Mvv = tovec(Mv);
            Zo(:, t0:t1) = (Mvv' * ftRo)';
            Ze(:, t0:t1) = (Mvv' * ftRe)';  clear Mv Mvv
        end
        Zo = 100 * Zo ./ F0Ro;  Ze = 100 * Ze ./ F0Re;
    end
    Zo = get_subthreshold(Zo, spk10, scr.sub_dilate, srcH2{i,2});
    Ze = get_subthreshold(Ze, spk10, scr.sub_dilate, srcH2{i,2});

    % --- regress the component time courses out, separately per channel ---
    % (finiteness is checked on the USABLE ROIs only: one unusable ROI would
    %  otherwise mark every frame bad, since all() runs down the ROI dimension)
    gH = good10 & all(isfinite(Zo(okPar,:)),1) & all(isfinite(Ze(okPar,:)),1) ...
                & all(isfinite(Ccomp),1);
    assert(nnz(gH) > pH.nPCmap + 2, '[R3-10h] too few analysed frames for the regression');
    nGh = nnz(gH);
    Cg  = Ccomp(:,gH) - mean(Ccomp(:,gH),2);
    prj = @(X) X - ((X*Cg')/(Cg*Cg'))*Cg;

    Zco = Zo(okPar,gH) - mean(Zo(okPar,gH),2);   Rro = prj(Zco);   % centred / residual
    Zce = Ze(okPar,gH) - mean(Ze(okPar,gH),2);   Rre = prj(Zce);

    vRes = nan(size(ftR,3),1);  vTot = nan(size(ftR,3),1);
    switch pH.residEst
        case 'oddeven'      % independent noise cancels at ANY smoothing window
            vRes(okPar) = sum(Rro .* Rre, 2) / (nGh-1);
            vTot(okPar) = sum(Zco .* Zce, 2) / (nGh-1);
            residTagH = 'odd \times even var';
        case 'lag1'         % valid only unsmoothed -- see pH.residEst
            if srcH2{i,2} > 1
                warning(['[R3-10h] ''lag1'' on the %d ms smoothed condition keeps ~%.0f%% of the ' ...
                         'shot noise ((w-1)/w) -- the residual is not comparable to the raw row'], ...
                         srcH2{i,2}, 100*(srcH2{i,2}-1)/srcH2{i,2});
            end
            Rrc = (Rro + Rre)/2;  Zcc = (Zco + Zce)/2;
            vRes(okPar) = mean(Rrc(:,1:end-1) .* Rrc(:,2:end), 2);
            vTot(okPar) = mean(Zcc(:,1:end-1) .* Zcc(:,2:end), 2);
            residTagH = 'lag-1 var';
        case 'var'          % includes shot noise
            Rrc = (Rro + Rre)/2;  Zcc = (Zco + Zce)/2;
            vRes(okPar) = var(Rrc, 0, 2);
            vTot(okPar) = var(Zcc, 0, 2);
            residTagH = 'ordinary var';
        otherwise
            error('pH.residEst must be ''oddeven'', ''lag1'' or ''var''');
    end
    residFrac(r) = mean(vRes,'omitnan') / mean(vTot,'omitnan');  % same measure both sides
    residH{r} = nan(size(wR));  residH{r}(suppR) = vRes(aR(suppR));
    rIm = nan(size(wR));  rIm(suppR) = vRes(aR(suppR));
    onBv = isfinite(rIm) & bvBW8;   offBv = isfinite(rIm) & ~bvBW8;
    fprintf('   %-22s beyond PC1-%d: %.4g = %.1f%% of total  (on bvMask %.4g / off %.4g = %.2f)\n', ...
            pcaNmH{i}, pH.nPCmap, mean(vRes,'omitnan'), 100*residFrac(r), ...
            mean(rIm(onBv)), mean(rIm(offBv)), mean(rIm(onBv))/mean(rIm(offBv)));
end

% ONE scale for every residual panel: with separate ranges each panel is
% re-normalised to itself and "less is left after correction" becomes invisible.
allRes = cell2mat(cellfun(@(R) R(isfinite(R)), residH(:), 'UniformOutput', false));
if isempty(pH.residClim), climRes = prctile(allRes, [0.1 99.9]); else, climRes = pH.residClim; end
if ~(isfinite(climRes(1)) && climRes(2) > climRes(1)), climRes = [min(allRes) max(allRes)+eps]; end
% Units of the residual maps. The traces are 100*X./F0, F0 = mean(mov_mc,3) over
% the footprint, so these are percent dF/F0 SQUARED -- NOT the robust dF/F that
% panel (i)'s trace is calibrated to (that calibration, kCal, is fit on one ROI
% and is never applied here).
residUnit = '(% \DeltaF/F_0)^2';
fprintf('   shared residual colour scale [%.4g %.4g]  in (%%dF/F0)^2\n', climRes(1), climRes(2));

for r = 1:nMapH
    i = srcH(r);
    for k = 1:pH.nPCmap
        ax = nexttile(rowH(r) + 2*(k-1), [1 2]);
        v = UH{i}(:,k);
        [~, im] = max(abs(v));  v = v * sign(v(im));       % eigenvector sign is arbitrary
        if strcmp(pH.pcaOn,'ftprnt')
            IM = nan(PT10.siz);  IM(PT10.supp) = v(PT10.aMax(PT10.supp));  alphaH = PT10.supp;
        else
            IM = nan(size(dMask8));  IM(dpix8) = v;        alphaH = dMask8;
        end
        cl = prctile(abs(v), 99.5);  if cl <= 0, cl = eps; end
        imagesc(ax, IM, 'AlphaData', alphaH); set(ax,'Color',[.12 .12 .12]);
        clim(ax, [-cl cl]); axis(ax,'image','off'); colormap(ax, cmapCorr);
        title(ax, sprintf('%s  PC%d  (%.2f%%)', pcaNmH{i}, k, 100*evH{i}(k)), 'FontSize', 8);
    end
    ax = nexttile(rowH(r) + 2*pH.nPCmap, [1 2]);          % the 4th column
    imagesc(ax, residH{r}, 'AlphaData', isfinite(residH{r})); set(ax,'Color',[.12 .12 .12]);
    clim(ax, climRes);                                    % SHARED by both rows
    axis(ax,'image','off'); colormap(ax, turbo); hold(ax,'on');
    if pH.residContour
        [~, hBv] = contour(ax, bvBW8, [0.5 0.5], 'Color', p8.vesCol, 'LineWidth', 1);
        if r == 1
            legend(ax, hBv, 'bvMask', 'Location','northeast', 'Box','off', ...
                   'TextColor','w', 'FontSize', 7);
        end
    end
    cbR = colorbar(ax);
    cbR.Label.String = residUnit;  cbR.Label.FontSize = 7;
    title(ax, sprintf('%s\nvariance beyond PC1-%d  (%.1f%% of total)\n%s, %s', ...
          pcaNmH{i}, pH.nPCmap, 100*residFrac(r), setTagH, residTagH), 'FontSize', 8);
end

fprintf('\n[R3-10h] figure 2 numbers:\n');
for i = 1:3
    ev = evH{i};
    evp = ev(ev > 0);
    fprintf('   %-22s PC1-3 = %5.1f%%   PC1-10 = %5.1f%%   PC>100 = %5.2f%%   (%d positive dims)\n', ...
            pcaNmH{i}, 100*sum(evp(1:min(3,end))), 100*sum(evp(1:min(10,end))), ...
            100*sum(evp(min(101,end+1):end)), numel(evp));
end
