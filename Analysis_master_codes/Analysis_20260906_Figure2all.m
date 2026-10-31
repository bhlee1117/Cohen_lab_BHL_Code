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
%foi=20;
%% Load Result files
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

%% SWC image (Figure 2a)
figure(100); clf;
f=20;
cmap_ExTr=gen_colormap(Plasma,10);
SWCpoints=PCresult{f}.swc;
SWCpoints(:,3)=SWCpoints(:,3)+5;
SWCpoints(1,3)=50;
SomaROI=[18 16 17]; ApicalROI=[32 33 34]; BasalROI=[7 8 9];
nROI=size(PCresult{f}.NormalizedTrace_dirt,1);
ftprint=PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist);
ROIvec=[1:nROI]; %ROIvec(ShowROI([35 37 39]))=2; ROIvec(ShowROI([3]))=1;
ftprint_s=max(ftprint(:,:,SomaROI),[],3);
ftprint_ap=max(ftprint(:,:,ApicalROI),[],3);
ftprint_bs=max(ftprint(:,:,BasalROI),[],3);
boundary_s=bwboundaries(ftprint_s); boundary_ap=cell2mat(bwboundaries(ftprint_ap)); boundary_bs=bwboundaries(ftprint_bs);
showScaleScatter(ROIvec, SWCpoints, ftprint,repmat([0 0 0],256,1)); hold all
plot(boundary_s{1}(:,1),boundary_s{1}(:,2),'color',cmap_ExTr(5,:),'linewidth',2);
plot(boundary_bs{1}(:,1),boundary_bs{1}(:,2),'color',cmap_ExTr(1,:),'linewidth',2);
plot(boundary_ap(:,1),boundary_ap(:,2),'color',cmap_ExTr(end,:),'linewidth',2);
drawScaleBar(100/Pixelsize(f),'vertical')

%% Voltage trace example, (Figure 2B–C)
f=20; cmap=gen_colormap(Plasma,3);

basalind=[7 8 9];
apicalind=[32 33 34];
somaind=[16 17 18];
panelBt2show=[45000:160000];

NormalizedTrace_noNaN=interpolateNaN(PCresult{f}.NormalizedTrace_dirt);
filteredNormTr2 = pcafilterTrace(NormalizedTrace_noNaN, setdiff([1:5],[4 5 6 7 10]));
ftprnt=PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist);
[nROI nTime]=size(filteredNormTr2);

Show_tr=[mean(filteredNormTr2(basalind,:),1,'omitnan'); mean(filteredNormTr2(somaind,:),1,'omitnan'); mean(filteredNormTr2(apicalind,:),1,'omitnan')];
Show_tr=interpolateNaN(Show_tr);
Show_subtr=get_subthreshold(Show_tr,PCresult{f}.SpikeMat(1,:)>0,7,20);
Show_ftprnt=cat(3,max(ftprnt(:,:,basalind)>0,[],3),max(ftprnt(:,:,somaind)>0,[],3),max(ftprnt(:,:,apicalind)>0,[],3));
roi_show=setdiff([1:nROI],[9 5 13 16 19 22]); %10 12 13 14 15 21
figure(5); clf; tiledlayout(4,1);
tax=[1:nTime]/1000; lscale=[0.2];

ax3=nexttile([3 1]);
l=plot(tax(panelBt2show),Show_tr(:,panelBt2show)'-[1:3]*lscale);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap,2));
hold all;
l2=plot(tax(panelBt2show),Show_subtr(:,panelBt2show)'-[1:3]*lscale);
arrayfun(@(l,c) set(l,'Color',c{:}),l2,num2cell(cmap/2,2));
%legend(l,{'Basal','Soma','Distal'})
axis off
drawScaleBar(lscale/2,'vertical','position',[162 -0.25],'color','k');
ylim([-0.67 0]);

ax2=nexttile([1 1]);
plot(tax(panelBt2show),PCresult{f}.VR(5,panelBt2show),'color',[0 0 0],'linewidth',2);
axis off
linkaxes([ax2 ax3],'x')
drawScaleBar(10,'horizontal','position',[162 -15],'color','k','Linewidth',3);
drawScaleBar(115,'vertical','position',[162 -15],'color','k','Linewidth',3);
xlim([45 162.5]); ylim([-17 115]);
set_figsize(200,170)

figure(6); clf; tiledlayout(5,1); %  zoom in version of figure(5);
tax=[1:size(Result.normTraces,2)]/1000;
zoomintime=[76900:79100];
NormalizedTrace_noNaN_zoomin=NormalizedTrace_noNaN(:,zoomintime);
filteredNormTr2_zoom=pcafilterTrace(NormalizedTrace_noNaN_zoomin,setdiff([1:37],[3 4 9 11 13 15 16]));

ax3=nexttile([3 1]);
Show_tr2show=Show_tr(:,zoomintime);
Show_subtr2show=Show_subtr(:,zoomintime);
l=plot(zoomintime-zoomintime(1),Show_tr2show'-[1:3]*0.2,'linewidth',1);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap,2));
hold all;
l=plot(zoomintime-zoomintime(1),Show_subtr2show'-[1:3]*lscale,'linewidth',1);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap/2,2));
axis off; ylim([-0.65 -0.1]);
drawScaleBar(100,'horizontal','position',[2200 -0.5],'color','k','Linewidth',3);
drawScaleBar(0.1,'vertical','position',[2200 -0.5],'color','k','Linewidth',3);

ax1=nexttile([2 1]);
kymo2show=filteredNormTr2_zoom(roi_show,:);
[~, dsort]=sort(PCresult{f}.dendaxis(roi_show));
kymo2show=kymo2show(dsort,:);
imagesc(zoomintime-zoomintime(1),[1:length(roi_show)],kymo2show,[-0.03 0.15]);
%set(gca,'YTick',[1:length(roi_show)],'YTickLabel',num2str(dsort'))
set_kymoYtick(PCresult{f}.dendaxis(roi_show))
colormap("turbo");
cb=colorbar;
cb.Label.String = '∆F/F';
ylabel('Distance from soma (µm)');
xlabel('Time (ms)');
set_font('Arial')
linkaxes([ax1 ax3],'x')
axis tight;
set_fontsize(16);
set_figsize(270,170);

%% ===== Correlation between ROIs — calculation (spatial & temporal) =====
% For each neuron and each frame condition (condName: Good, Silent, Run, PF) a
% single spatiotemporal correlation tensor Cmat = N x N x (2*maxLag+1) is built
% with get_corrTensor (raw traces, nanXCorrFFT on a continuous time axis).
%   * spatial correlation matrix = LAG-1 correlation (get_lag1cov, see below)
%   * temporal correlation       = the lag profile  Cmat(i,j,:)
% The spatial matrix is now the LAG-1 correlation of the combined trace rather
% than the zero-lag cross-channel (odd/even) matrix. Shot noise is temporally
% white, so it contributes nothing to X(t)X(t+1)' -- including on the diagonal,
% which is the only place the odd/even split was needed. One trace per ROI
% therefore suffices, and the same estimator can be used on movies that have no
% two independent-noise halves. The correlogram TENSOR is unchanged and still
% uses tr1 x tr2 on the diagonal. Spatial decay (corr vs distance)
% and temporal decay (corr vs lag) are both extracted from this tensor. Auto
% class pairs (dClass2temporal{dc,1}==dClass2temporal{dc,2}) use diagonal pairs
% only; dClass2temporal row 5 ({1..5} auto) is the all-ROI pooled autocorr.
% Separate cache from the cross-channel version: Correlation_Fig2_data.mat holds
% the zero-lag odd/even corrMat and is still read by Analysis_2026????_Response2reviewer3
% and Analysis_20260522_Figure2all. Do not overwrite it with the lag-1 matrices.
corrDataFile = fullfile(save_figto, 'Correlation_Fig2_Lag1data.mat');   % tensors + fit inputs
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
        % spatial matrix from the LAG-1 covariance of the combined trace: white
        % shot noise drops out of X(t)X(t+1)' on the diagonal as well, so this
        % needs no channel split. Raw (unsmoothed) traces, as required.
        [~, corrMat{f,c}] = get_lag1cov(trC, condMask{c});
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

%% ===== Exponential fits (edit fit function / range / bounds here, then re-run this cell) =====
% Cheap and re-runnable: the expensive tensors above are cached in corrDataFile,
% so tweak the settings below and re-run only this cell. Model has NO offset
% (c0 = 0): y = sum_i c_i*exp(-x/tau_i), decaying to 0. Single vs double
% exponential is set by the length of fitGuess (scalar -> single, [g1 g2] -> double).
fitGuess = struct('spatial', 50, 'temporal', [10 100]);   % 1/e guess per FitType (scalar/2-vec)
fitXmax  = struct('spatial', Inf, 'temporal', Inf);       % include only x <= this in the fit
fitTauLB = struct('spatial', 0,   'temporal', 0);         % lower bound on tau
fitTauUB = struct('spatial', Inf, 'temporal', Inf);       % upper bound on tau (set <=0 -> growing exp; then also set LB<0)
expfit   = @(x, y, xs, g, lb, ub) expfit_bnd(x, y, xs, g, lb, ub);   % swap in a different fit here

nRow = height(CorrelationData);
Tau = nan(nRow,1);  R2 = nan(nRow,1);  Yfit = cell(nRow,1);  TauAll = cell(nRow,1);
for r = 1:nRow
    ft    = char(CorrelationData.FitType(r));
    xfull = CorrelationData.Xdata{r};   yfull = CorrelationData.Ydata{r};
    inR   = xfull <= fitXmax.(ft);
    x = xfull(inR);  y = yfull(inR);
    if numel(x) <= 5, continue; end
    [yhatFull, tcon] = expfit(x, y, xfull, fitGuess.(ft), fitTauLB.(ft), fitTauUB.(ft));   % fit on range, eval at full x
    Tau(r)    = max(tcon);              % slowest 1/e component (for multi-exp guesses)
    TauAll{r} = tcon;
    R2(r)     = 1 - sum((y - yhatFull(inR)).^2) / sum((y - mean(y)).^2);
    Yfit{r}   = yhatFull;
end
CorrelationFitResult = CorrelationData;
CorrelationFitResult.Tau    = Tau;
CorrelationFitResult.TauAll = TauAll;
CorrelationFitResult.R2     = R2;
CorrelationFitResult.Yfit   = Yfit;

%% ===== Interactive review / refinement of each exponential fit =====
% Step through the fits in CorrelationFitResult and refine each one with the
% interactive GUI (lasso-select the points to include, live a*exp(-x/b) fit).
% Refined a / b(=tau) / R^2 are written into new table columns; the automatic
% Tau/R2 columns are left untouched. Set reviewRows to a subset to review only
% some fits, e.g. find(CorrelationFitResult.FitType=="temporal").

%reviewRows = 1:height(CorrelationFitResult);
reviewRows = find(CorrelationFitResult.ClassPair=="All, All (auto)");

nRow = height(CorrelationFitResult);
A_manual = nan(nRow,1);  Tau_manual = nan(nRow,1);  R2_manual = nan(nRow,1);
FitManual = cell(nRow,1);
for r = reviewRows(:)'
    X = CorrelationFitResult.Xdata{r};
    Y = CorrelationFitResult.Ydata{r};
    fprintf('Review fit %d/%d | f=%d | %s | %s | %s\n', ...
        find(reviewRows==r,1), numel(reviewRows), CorrelationFitResult.FileInd(r), ...
        CorrelationFitResult.FitType(r), CorrelationFitResult.ClassPair(r), ...
        CorrelationFitResult.Condition(r));

    ft2 = char(CorrelationFitResult.FitType(r));
    opt = struct();
    opt.x_fit   = 0:ceil(max(X));
    opt.t_guess = fitGuess.(ft2);   % single/double, same guess/bounds as the batch fit
    opt.tau_lb  = fitTauLB.(ft2);
    opt.tau_ub  = fitTauUB.(ft2);
    res = interactive_expfit(X, Y, opt);
    if isfield(res,'aborted') && res.aborted, break; end   % window aborted -> keep auto fit

    FitManual{r}  = res;
    A_manual(r)   = res.a;
    Tau_manual(r) = res.b;
    R2_manual(r)  = res.Rsq;
end
CorrelationFitResult.A_manual   = A_manual;
CorrelationFitResult.Tau_manual = Tau_manual;
CorrelationFitResult.R2_manual  = R2_manual;
CorrelationFitResult.FitManual  = FitManual;

%% ===== Save fitted results =====
save(fullfile(save_figto, 'Correlation_Fig2_Lag1results.mat'), ...
    'corrMat', 'condName', 'heatCond', ...
    'temporalCorr', 'temporalCorrMat', 'lagsCorr', ...
    'CorrelationData', 'CorrelationFitResult', 'dClass2temporal', 'title_str', 'classNames', ...
    'nPair', 'nCond', 'maxLag', 'foi', '-v7.3');

% Decay-constant summary (1/e) from the interactively refined fits (Tau_manual), silent condition
T = CorrelationFitResult;
validSp = rmmissing(T.Tau_manual(T.FitType=="spatial" & T.ClassPair=="All, All (auto)" & T.Condition=="Silent" & T.R2_manual>0.3));
fprintf('Spatial decay length (all pairs) = %.0f +/- %.0f um (n=%d)\n', ...
    mean(validSp), std(validSp), numel(validSp));
for dc = 1:nPair
    ss = rmmissing(T.Tau_manual(T.FitType=="spatial" & T.ClassPair==string(title_str{dc}) & T.Condition=="Silent" & T.R2_manual>0.3));
    fprintf('Spatial length [%s] = %.0f +/- %.0f um (n=%d)\n', ...
        title_str{dc}, mean(ss), std(ss), numel(ss));
end
for dc = 1:nPair
    tt = rmmissing(T.Tau_manual(T.FitType=="temporal" & T.ClassPair==string(title_str{dc}) & T.Condition=="Silent" & T.R2_manual>0.3));
    fprintf('Temporal decay [%s] = %.0f +/- %.0f frames (n=%d)\n', ...
        title_str{dc}, mean(tt), std(tt), numel(tt));
end
%% Show example Spatio-temporal correlation
showf=23; showind=[3 39];
figure(117); clf;
for i=showind
    dax=PCresult{showf}.dendaxis;
    [dax dsort]=sort(dax);
    tax=size(temporalCorrMat{showf},3)-round(size(temporalCorrMat{showf},3)/2);
    nexttile([1 1]);
    imagesc(tax,[1:length(dax)],squeeze(temporalCorrMat{showf}(i,dsort,:)))
    set_kymoYtick(dax);
end


%% ===== Figure 2D — subthreshold correlation matrix (lag-1), leaf-ordered + dendrogram =====
% LAG-1 correlation matrix for one neuron over the SILENT period (blue-off AND
% not peri-spike). Replaces the zero-lag cross-channel (ch1 x ch2) matrix: white
% shot noise contributes nothing to X(t)X(t+1)', on the diagonal as well, so no
% channel split is needed.
%
% NOTE the input is the SPIKE-REMOVED but UNSMOOTHED trace (fig2dWin = 1), not
% the 20-ms subthreshold trace. After a k-frame moving average shot noise has
% lag-1 autocorrelation ~ (k-1)/k, so smoothing would put most of it straight
% back onto the diagonal and defeat the estimator. Set fig2dWin = 17 to
% reproduce the old smoothed input, but the diagonal is then no longer clean.
%
% Reordered by hierarchical-clustering leaf order, with the dendrogram on the
% left and dendrite-class colour dots on the top and left.
f = 20;
Pr  = PCresult{f};
trC = Pr.NormalizedTrace_dirt;          % lag-1 needs one trace per ROI, not the two channels
spk = Pr.SpikeMat(1,:) > 0;
goodvec = ~imdilate(Pr.BlueStim>0, ones(1,201));
% silent period = blue-off AND not peri-spike ( +/- [-3:10] around spikes, + CS )
perispk = unique([tovec(find(spk)' + (-3:10)); find(Pr.CStrace)']);
perispk(perispk < 1 | perispk > numel(spk)) = [];
silentvec = true(1, numel(spk));  silentvec(perispk) = false;
silentMask = goodvec & silentvec;
fig2dWin = 1;                 % get_subthreshold averaging window: 1 = none (see note above)
subU = get_subthreshold(trC, spk, 7, fig2dWin);   % spikes removed, NOT smoothed
subU(isnan(trC)) = NaN;
[~, C, nPr2D] = get_lag1cov(subU, silentMask);
fprintf('Fig 2D: lag-1 correlation, window %d, %d consecutive silent frame pairs\n', ...
        fig2dWin, nPr2D);
n = size(C,1);
lab = PCresult{f}.labelvec;                        % dendrite class per ROI (1..5, NaN=unclassified)
cmap_label = hsv(5);  cmap_label = cmap_label([4 2 3 5 1],:);

% hierarchical clustering + optimal leaf order (distance = 1 - correlation)
D = 1 - C;  D = (D+D')/2;  D(1:n+1:end) = 0;  D(isnan(D)) = 1;
dvec = squareform(D);  Z = linkage(dvec, 'average');
try, ord = optimalleaforder(Z, dvec); catch, ord = 1:n; end

% simple 3-colour blue-white-red map + caxis (as before; white at the caxis midpoint)
cmapCorr = gen_colormap([0 0.5 1; 1 1 1; 1 0 0], 256);
climc = [-0.25 1];

% put basal (class 1) at the top of the ordering
if mean(find(lab(ord)==1)) > (n+1)/2, ord = fliplr(ord); end

figure(118); clf;
axD = axes('Position',[0.06  0.13 0.13  0.72]);    % dendrogram (left)
axL = axes('Position',[0.200 0.13 0.02  0.72]);    % left class dots (rows)
axM = axes('Position',[0.225 0.13 0.55  0.72]);    % correlation matrix
axT = axes('Position',[0.225 0.86 0.55  0.02]);    % top class dots (columns)

% --- dendrogram: flipped L-R, basal on top; leaf k at data-y=k ---
axes(axD);
[Hd,~,outperm] = dendrogram(Z, 0, 'Orientation','left', 'Reorder', ord, 'ColorThreshold','default');
set(Hd,'LineWidth',1);
set(axD,'YLim',[0.5 n+0.5], 'YDir','reverse');  axis(axD,'off');
ordUse = outperm;                                  % leaf order along the axis (data-y 1..n)

% --- correlation matrix: fills its Position box (no axis image), row 1 at top ---
axes(axM);
imagesc(C(ordUse,ordUse), climc);
set(axM,'XLim',[0.5 n+0.5], 'YLim',[0.5 n+0.5], 'YDir','reverse');  axis(axM,'off');
colormap(axM, cmapCorr);
cb = colorbar(axM);  cb.Label.String = 'Correlation coefficient';
cb.Position = [0.79 0.13 0.025 0.72];  axM.Position = [0.225 0.13 0.55 0.72];   % restore after colorbar

% --- class colour dots (share axM row/column limits -> aligned; NaN class -> grey) ---
labc = lab(ordUse);  dotcol = repmat([0.7 0.7 0.7], n, 1);
vv = ~isnan(labc);   dotcol(vv,:) = cmap_label(labc(vv),:);
axes(axT);  scatter(1:n, ones(1,n), 26, dotcol, 'filled');
set(axT,'XLim',[0.5 n+0.5], 'YLim',[0.5 1.5]);  axis(axT,'off');
axes(axL);  scatter(ones(1,n), 1:n, 26, dotcol, 'filled');
set(axL,'XLim',[0.5 1.5], 'YLim',[0.5 n+0.5], 'YDir','reverse');  axis(axL,'off');

set_font('Arial'); set_fontsize(14);  set_figsize(133,100);

%% ===== Figure 119 — Fig 2D matrix: three ways =====
% Same neuron, same silent frames, same ROI order. Three matrices, arranged so
% that smoothing and estimator can be told apart:
%   (a) 20 ms smoothed, ODD/EVEN cross-channel -- the original Figure 2D.
%   (b) 20 ms smoothed, LAG-1                  -- same data as (a), other estimator.
%   (c) unsmoothed,     LAG-1                  -- the new Figure 2D.
% (a) vs (b) isolates the ESTIMATOR at fixed data; (b) vs (c) isolates the
% SMOOTHING at fixed estimator; (a) vs (c) compares the two defensible choices.
%
% CAVEAT on (b): lag-1 cancels shot noise only on unsmoothed traces. After a
% k-frame moving average the noise keeps a lag-1 autocorrelation of ~(k-1)/k, so
% most of it survives into the lag-1 product and inflates the diagonal, which
% deflates every correlation. (b) is here to show the estimator's effect, not as
% a defensible matrix in its own right -- (a) and (c) are the clean ones.
%
% All three are symmetrised for display (the ch1 x ch2 matrix is not symmetric);
% the leaf order comes from (c) and is shared, so differences show up as colour
% rather than as reshuffling.
f = 20;
Pr  = PCresult{f};
smWin = 20;                        % moving-average window in frames (~20 ms at 1 kHz)
trC = Pr.NormalizedTrace_dirt;
ch1 = Pr.NormalizedTrace_ch{1};  ch2 = Pr.NormalizedTrace_ch{2};
spk = Pr.SpikeMat(1,:) > 0;
goodvec   = ~imdilate(Pr.BlueStim>0, ones(1,201));
perispk   = unique([tovec(find(spk)' + (-3:10)); find(Pr.CStrace)']);
perispk(perispk < 1 | perispk > numel(spk)) = [];
silentvec = true(1, numel(spk));  silentvec(perispk) = false;
silentMask = goodvec(:)' & silentvec;

% (a) smoothed, odd x even
s1 = get_subthreshold(ch1, spk, 7, smWin);  s1(isnan(trC)) = NaN;
s2 = get_subthreshold(ch2, spk, 7, smWin);  s2(isnan(trC)) = NaN;
Csm  = get_corrMat(s1, s2, find(silentMask));   Csm = (Csm + Csm')/2;

% (b) smoothed, lag-1  (NOT shot-noise free -- see caveat above)
sC = get_subthreshold(trC, spk, 7, smWin);  sC(isnan(trC)) = NaN;
[~, Csl] = get_lag1cov(sC, silentMask);

% (c) unsmoothed, lag-1
subU = get_subthreshold(trC, spk, 7, 1);  subU(isnan(trC)) = NaN;
[~, Clag] = get_lag1cov(subU, silentMask);

nC = size(Clag,1);
Dc = 1 - Clag;  Dc = (Dc+Dc')/2;  Dc(1:nC+1:end) = 0;  Dc(isnan(Dc)) = 1;
dvecC = squareform(Dc);  Zc = linkage(dvecC, 'average');
try, ordC = optimalleaforder(Zc, dvecC); catch, ordC = 1:nC; end
labC = Pr.labelvec;
if mean(find(labC(ordC)==1)) > (nC+1)/2, ordC = fliplr(ordC); end   % basal on top

cmapCorr2 = gen_colormap([0 0.5 1; 1 1 1; 1 0 0], 256);
climc2    = [-0.25 1];

matsC  = {Csm, Csl, Clag};
namesC = {sprintf('(a) %d ms smoothed, odd x even', smWin), ...
          sprintf('(b) %d ms smoothed, lag-1', smWin), ...
          '(c) unsmoothed, lag-1'};
mUp = triu(true(nC),1);                       % off-diagonal pairs; the diagonal is 1 by construction

figure(119); clf; tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for q = 1:3
    ax = nexttile(q);
    imagesc(ax, matsC{q}(ordC,ordC), climc2);  axis(ax,'image');
    set(ax,'XTick',[],'YTick',[]);  colormap(ax, cmapCorr2);
    title(ax, namesC{q}, 'FontSize', 9);
    if q == 3, cb = colorbar(ax); cb.Label.String = 'Correlation coefficient'; end
end

% pairwise scatters: estimator effect, smoothing effect, and the two clean ones
cmpC = [1 2; 2 3; 1 3];
cmpL = {'estimator (data fixed)', 'smoothing (estimator fixed)', 'the two clean matrices'};
for q = 1:3
    a = cmpC(q,1);  b = cmpC(q,2);
    xs = matsC{a}(mUp);  ys = matsC{b}(mUp);
    okC = isfinite(xs) & isfinite(ys);
    rC  = corr(xs(okC), ys(okC));
    ax = nexttile(3+q);  hold(ax,'on');
    plot(ax, xs(okC), ys(okC), '.', 'Color', [.5 .5 .5], 'MarkerSize', 8);
    plot(ax, climc2, climc2, 'k--', 'LineWidth', 1);
    xlabel(ax, namesC{a}); ylabel(ax, namesC{b});
    xlim(ax, climc2); ylim(ax, climc2); axis(ax,'square'); box(ax,'off'); grid(ax,'on');
    title(ax, sprintf('%s:  r = %.3f', cmpL{q}, rC), 'FontSize', 9);
end

fprintf('\nFig 2D correlation matrix, f = %d, %d ROIs, %d silent frames\n', f, nC, nnz(silentMask));
for q = 1:3
    vq = matsC{q}(mUp);  vq = vq(isfinite(vq));
    fprintf('  %-32s mean off-diag %.3f   (diag mean %.3f)\n', ...
            namesC{q}, mean(vq), mean(diag(matsC{q}), 'omitnan'));
end
for q = 1:3
    a = cmpC(q,1);  b = cmpC(q,2);
    xs = matsC{a}(mUp);  ys = matsC{b}(mUp);
    okC = isfinite(xs) & isfinite(ys);
    fprintf('  (%d) vs (%d)  r = %.3f  slope %.3f   [%s]\n', ...
            a, b, corr(xs(okC), ys(okC)), xs(okC)\ys(okC), cmpL{q});
end
set_font('Arial'); set_fontsize(11);

%% ===== Figure 217 — cumulative variance vs PC number (lag-1 or odd x even covariance) =====
% How many dimensions the subthreshold dynamics actually occupy, across all
% neurons in foi. The eigenspectrum comes from the LAG-1 covariance rather than
% the ordinary same-time covariance: shot noise inflates every same-time
% variance and lands entirely on the diagonal, which flattens the tail of the
% spectrum and makes the dynamics look higher-dimensional than they are. Lag-1
% is unbiased for signal in every dimension.
%
% As in Fig 2D the input is the SPIKE-REMOVED but UNSMOOTHED trace (cumWin = 1):
% a k-frame moving average would give the noise a lag-1 autocorrelation of
% ~(k-1)/k and put it straight back in. This also means the curve is NOT
% contaminated by the 20-ms smoothing.
%
% Because the estimator is unbiased it produces NEGATIVE eigenvalues once the
% signal is exhausted -- that is the estimator fluctuating about zero, not a
% failure. The cumulative fraction is taken over the POSITIVE eigenvalues only,
% and nPosD records how many each neuron has.
%
% cumEst = 'oddeven' uses the ODD x EVEN cross-channel covariance instead
% (NormalizedTrace_ch{1} x {2}, symmetrised). Shot noise is independent between
% the two channels, so it cancels on the diagonal too -- and unlike lag-1 it
% stays cancelled after smoothing, so cumWin > 1 (e.g. 20) is legitimate there.
% Same positive-eigenvalue rule applies.
cumEst = 'oddeven';                  % 'lag1' or 'oddeven'
cumWin = 20;                       % get_subthreshold averaging window; 1 = none
if strcmp(cumEst, 'lag1') && cumWin > 1
    warning('Fig217:lag1Smoothed', 'lag-1 on a %d-frame smoothed trace is NOT shot-noise free', cumWin);
end
cumsumD = cell(1, max(foi));  SubV = cell(1, max(foi));  SubD = cell(1, max(foi));
nPosD   = nan(1, max(foi));

figure(217); clf; hold on;
for f = foi
    Prc = PCresult{f};
    trCc = Prc.NormalizedTrace_dirt;
    spkc = Prc.SpikeMat(1,:) > 0;
    goodc = ~imdilate(Prc.BlueStim > 0, ones(1,201));   % blue off

    switch cumEst
        case 'lag1'
            subc = get_subthreshold(trCc, spkc, 7, cumWin);
            subc(isnan(trCc)) = NaN;
            goodc = goodc(:)' & all(isfinite(subc), 1);      % finite in every ROI
            Ccov = get_lag1cov(subc, goodc);
        case 'oddeven'
            sub1 = get_subthreshold(Prc.NormalizedTrace_ch{1}, spkc, 7, cumWin);  sub1(isnan(trCc)) = NaN;
            sub2 = get_subthreshold(Prc.NormalizedTrace_ch{2}, spkc, 7, cumWin);  sub2(isnan(trCc)) = NaN;
            goodc = goodc(:)' & all(isfinite(sub1), 1) & all(isfinite(sub2), 1);
            A1 = sub1(:, goodc) - mean(sub1(:, goodc), 2);
            A2 = sub2(:, goodc) - mean(sub2(:, goodc), 2);
            Ccov = (A1*A2' + A2*A1') / 2 / (nnz(goodc) - 1);
        otherwise
            error('Fig217:cumEst', 'cumEst must be ''lag1'' or ''oddeven''');
    end
    [V, D] = eig((Ccov + Ccov')/2, 'vector');
    [D, si] = sort(D, 'descend');   V = V(:, si);
    SubV{f} = V;  SubD{f} = D;

    pos = D > 0;  nPosD(f) = nnz(pos);
    cumsumD{f} = [(1:nPosD(f))', cumsum(D(pos))/sum(D(pos))];
    fprintf('  f=%d: %d ROIs, %d positive %s eigenvalues, PC1-3 = %.1f%%\n', ...
            f, numel(D), nPosD(f), cumEst, 100*cumsumD{f}(min(3,end),2));

    plot(cumsumD{f}(:,1), cumsumD{f}(:,2), 'color', [0.7 0.7 0.7],'LineWidth',1);
    drawnow;
end

VarianceResult = binning_data(cumsumD(foi), 0.5:1:50.5);
errorbar(VarianceResult.centers, VarianceResult.mean, VarianceResult.std, 'color', 'k','LineWidth',1.5);
set(gca, 'xscale', 'linear');%  xlim([1 50]);  ylim([0 1]);  box off;
xlabel('PC component #');  ylabel('Fraction of variance');
%title(sprintf('%s covariance, window %d, n = %d neurons', cumEst, cumWin, numel(foi)));
xlim([0.5 50]); ylim([0 1.1]);
set(gca,'XTick',[1 10 20 30 40 50],'LineWidth',1);
set_font('Arial'); set_fontsize(16); set_figsize(90,115);

%% ===== Figure 218 — moving average vs Savitzky-Golay: CorrMat (f=20) and cumulative variance =====
% Same pipeline as get_subthreshold (spikes +/-3 frames removed, linear
% interpolation), only the final smoothing step differs:
%   MA       : movmean, maWin frames                 (= get_subthreshold(..., maWin))
%   SG same  : sgolayfilt, order sgOrd, same span     (maWin+1, framelen must be odd)
%   SG match : sgolayfilt, order sgOrd, span chosen so its -3 dB cutoff matches MA
% At equal span SG passes far more high frequency than MA (it preserves peak
% shape by design), so "same span" and "same cutoff" answer different questions.
% Cutoffs: MA ~0.443/k * fs; SG ~(N+1)/(3.2M-4.6) * fs (Schafer 2011).
%
% ESTIMATOR: smoothed traces cannot go through lag-1 (the smoothing gives the
% shot noise a lag-1 autocorrelation, see Fig 119 caveat), so both the matrix
% and the eigenspectrum use the ODD x EVEN cross-channel product, whose noise
% cancels regardless of smoothing. Matrix = Fig 119 (a); spectrum = eig of the
% symmetrised cross-channel covariance, positive eigenvalues only (as Fig 217).
fsSub = 1000;                     % frame rate (Hz)
maWin = 20;                       % moving-average window (frames)
sgOrd = 2;                        % Savitzky-Golay polynomial order
fcMA  = 0.443/maWin;                                     % fraction of fs
fcSGf = @(N,M) (N+1)/(3.2*M - 4.6);
sgSame  = maWin + ~mod(maWin,2);
sgMatch = round(((sgOrd+1)/fcMA + 4.6)/3.2);  sgMatch = sgMatch + ~mod(sgMatch,2);

smSpec = struct('type', {'ma', 'sg', 'sg'}, 'len', {maWin, sgSame, sgMatch}, 'ord', {NaN, sgOrd, sgOrd});
nSm = numel(smSpec);
smName = cell(1, nSm);  smFc = nan(1, nSm);
for m = 1:nSm
    if strcmp(smSpec(m).type, 'ma')
        smFc(m)   = fcMA*fsSub;
        smName{m} = sprintf('MA %d', smSpec(m).len);
    else
        smFc(m)   = fcSGf(smSpec(m).ord, smSpec(m).len)*fsSub;
        smName{m} = sprintf('SG %d/%d', smSpec(m).ord, smSpec(m).len);
    end
    fprintf('  %-10s -3 dB cutoff ~ %.1f Hz\n', smName{m}, smFc(m));
end
smCol = [0 0 0; 0.85 0.33 0.10; 0 0.45 0.74];

% ---- (1) representative correlation matrix, f = 20 ----
f = 20;
Pr  = PCresult{f};
trC = Pr.NormalizedTrace_dirt;
spk = Pr.SpikeMat(1,:) > 0;
goodvec   = ~imdilate(Pr.BlueStim>0, ones(1,201));
perispk   = unique([tovec(find(spk)' + (-3:10)); find(Pr.CStrace)']);
perispk(perispk < 1 | perispk > numel(spk)) = [];
silentvec = true(1, numel(spk));  silentvec(perispk) = false;
silentMask = goodvec(:)' & silentvec;

Csm = cell(1, nSm);
for m = 1:nSm
    s1 = smooth_subthreshold(Pr.NormalizedTrace_ch{1}, spk, smSpec(m));  s1(isnan(trC)) = NaN;
    s2 = smooth_subthreshold(Pr.NormalizedTrace_ch{2}, spk, smSpec(m));  s2(isnan(trC)) = NaN;
    Cm = get_corrMat(s1, s2, find(silentMask));
    Csm{m} = (Cm + Cm')/2;
end

% shared leaf order from the MA matrix, basal on top
nC = size(Csm{1}, 1);
Dc = 1 - Csm{1};  Dc = (Dc+Dc')/2;  Dc(1:nC+1:end) = 0;  Dc(isnan(Dc)) = 1;
dvecC = squareform(Dc);  Zc = linkage(dvecC, 'average');
try, ordC = optimalleaforder(Zc, dvecC); catch, ordC = 1:nC; end
if mean(find(Pr.labelvec(ordC)==1)) > (nC+1)/2, ordC = fliplr(ordC); end
cmapCorr3 = gen_colormap([0 0.5 1; 1 1 1; 1 0 0], 256);
climc3    = [-0.25 1];
mUp = triu(true(nC), 1);

% ---- (2) cumulative variance, all neurons ----
cumsumSm = cell(nSm, max(foi));
nPC80    = nan(nSm, max(foi));                % PCs needed to reach 80 % of the positive variance
for f = foi
    Prc  = PCresult{f};
    trCc = Prc.NormalizedTrace_dirt;
    spkc = Prc.SpikeMat(1,:) > 0;
    goodc0 = ~imdilate(Prc.BlueStim > 0, ones(1,201));
    for m = 1:nSm
        s1 = smooth_subthreshold(Prc.NormalizedTrace_ch{1}, spkc, smSpec(m));  s1(isnan(trCc)) = NaN;
        s2 = smooth_subthreshold(Prc.NormalizedTrace_ch{2}, spkc, smSpec(m));  s2(isnan(trCc)) = NaN;
        g  = goodc0(:)' & all(isfinite(s1), 1) & all(isfinite(s2), 1);
        A1 = s1(:, g) - mean(s1(:, g), 2);
        A2 = s2(:, g) - mean(s2(:, g), 2);
        Cx = (A1*A2' + A2*A1') / 2 / (nnz(g) - 1);       % cross-channel covariance, noise-free in expectation
        D  = sort(eig((Cx + Cx')/2), 'descend');
        D  = D(D > 0);
        cumsumSm{m, f} = [(1:numel(D))', cumsum(D)/sum(D)];
        nPC80(m, f) = find(cumsumSm{m, f}(:,2) >= 0.8, 1);
    end
    fprintf('  f=%d: PCs to 80%% = %s\n', f, sprintf('%d  ', nPC80(:, f)));
end

% ---- plot ----
figure(218); clf; tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for m = 1:nSm
    ax = nexttile(m);
    imagesc(ax, Csm{m}(ordC, ordC), climc3);  axis(ax, 'image');
    set(ax, 'XTick', [], 'YTick', []);  colormap(ax, cmapCorr3);
    vq = Csm{m}(mUp);
    title(ax, sprintf('%s (f=20)\nmean off-diag %.3f', smName{m}, mean(vq, 'omitnan')), 'FontSize', 9);
    if m == nSm, cb = colorbar(ax); cb.Label.String = 'Correlation coefficient'; end
end

% filter frequency responses (what each smoother actually passes)
ax = nexttile(4);  hold(ax, 'on');
for m = 1:nSm
    if strcmp(smSpec(m).type, 'ma')
        h = ones(1, smSpec(m).len) / smSpec(m).len;
    else
        B = sgolay(smSpec(m).ord, smSpec(m).len);  h = B((smSpec(m).len+1)/2, :);
    end
    [H, w] = freqz(h, 1, 2048, fsSub);
    plot(ax, w, 20*log10(abs(H)), 'Color', smCol(m,:), 'LineWidth', 1.5, 'DisplayName', smName{m});
end
yline(ax, -3, 'k:');  xlim(ax, [0 200]);  ylim(ax, [-40 3]);  box(ax, 'off');
xlabel(ax, 'Frequency (Hz)');  ylabel(ax, 'Gain (dB)');  legend(ax, 'Location', 'northeast', 'Box', 'off');
title(ax, 'Filter response', 'FontSize', 9);

% MA vs each SG, off-diagonal pairs
for m = 2:nSm
    ax = nexttile(3 + m);  hold(ax, 'on');
    xs = Csm{1}(mUp);  ys = Csm{m}(mUp);  ok = isfinite(xs) & isfinite(ys);
    plot(ax, xs(ok), ys(ok), '.', 'Color', smCol(m,:), 'MarkerSize', 6);
    plot(ax, climc3, climc3, 'k--');
    xlim(ax, climc3);  ylim(ax, climc3);  axis(ax, 'square');  box(ax, 'off');
    xlabel(ax, smName{1});  ylabel(ax, smName{m});
    title(ax, sprintf('r = %.3f, slope %.3f', corr(xs(ok), ys(ok)), xs(ok)\ys(ok)), 'FontSize', 9);
end

% cumulative variance, mean +/- std across neurons
ax = nexttile(7);  hold(ax, 'on');
for m = 1:nSm
    VR = binning_data(cumsumSm(m, foi), 0.5:1:50.5);
    errorbar(ax, VR.centers, VR.mean, VR.std, 'Color', smCol(m,:), 'LineWidth', 1, 'DisplayName', smName{m});
end
set(ax, 'XScale', 'log');  xlim(ax, [1 50]);  ylim(ax, [0 1]);  box(ax, 'off');
xlabel(ax, 'PC component #');  ylabel(ax, 'Fraction of variance');
legend(ax, 'Location', 'southeast', 'Box', 'off');
title(ax, sprintf('odd x even covariance, n = %d neurons', numel(foi)), 'FontSize', 9);

% PCs to 80 %, paired per neuron
ax = nexttile(8);  hold(ax, 'on');
plot(ax, 1:nSm, nPC80(:, foi), '-', 'Color', [0.75 0.75 0.75]);
for m = 1:nSm
    plot(ax, m, nPC80(m, foi), 'o', 'Color', smCol(m,:), 'MarkerFaceColor', smCol(m,:), 'MarkerSize', 4);
end
set(ax, 'XTick', 1:nSm, 'XTickLabel', smName);  xlim(ax, [0.5 nSm+0.5]);  box(ax, 'off');
ylabel(ax, '# PCs to 80% variance');
for m = 2:nSm
    pS = signrank(nPC80(1, foi), nPC80(m, foi));
    fprintf('  %s vs %s: median PCs to 80%% %g vs %g, signrank p = %.3g\n', ...
            smName{1}, smName{m}, median(nPC80(1, foi)), median(nPC80(m, foi)), pS);
end
title(ax, 'Dimensionality per neuron', 'FontSize', 9);
set_font('Arial'); set_fontsize(11);

%% ===== Figure 219 — INTERACTIVE: pick global / see-saw PCs per neuron =====
% Uses the eigenvectors/eigenvalues from Figure 217 (SubV, SubD), so run that
% cell first; the estimator it used (cumEst, cumWin) is stored with the result.
% For each neuron the top modeNShow PCs are drawn as loadings on the morphology
% (sign flipped so the loadings sum >= 0; blue-white-red, symmetric range).
% modeView = 'dendrite' colours the SWC trace inside each ROI footprint with that
% ROI's loading (grey = traced dendrite outside any ROI); 'footprint' paints the
% smoothed ROI footprints instead.
%   left click  : toggle GLOBAL  (red frame)
%   right click : toggle SEE-SAW (blue frame)     (ctrl-click on a Mac trackpad)
%   Enter / n / ->  : accept, next neuron      b / <-  : back one neuron
%   c : clear this neuron                      q / Esc / close window : save & stop
% More than one PC can be given the same mode; their variances are summed.
% Variance explained = eigenvalue / sum of ALL positive eigenvalues (as in the
% Methods). u = |sum v| / (sqrt(n)*||v||) is printed per PC as a guide: 1 = all
% ROIs equal and same sign (global-like), ~0 = loadings cancel (see-saw-like).
% The result is saved after every neuron, so stopping halfway loses nothing;
% re-running resumes at the first neuron with no selection (modeRedo = true
% starts over, keeping the saved selections as the starting point).
assert(exist('SubV', 'var') == 1 && exist('SubD', 'var') == 1, 'Run Figure 217 first (SubV / SubD)');
modeNShow = 15;                   % PCs shown per neuron
modeFOI   = foi;                  % neurons to review
modeRedo  = false;                % true: review every neuron in modeFOI again
modeView  = 'dendrite';           % 'dendrite' (SWC trace) or 'footprint'
modeFile  = fullfile(save_figto, 'PCmode_selection_Fig2.mat');

modeTmpl = struct('globalPC', [], 'seesawPC', [], 'varGlobal', NaN, 'varSeesaw', NaN, ...
                  'varFrac', [], 'uIdx', [], 'globalVec', [], 'seesawVec', [], 'done', false);
if exist(modeFile, 'file')
    S = load(modeFile, 'ModeSel', 'modeMeta');  ModeSel = S.ModeSel;
    if ~strcmp(S.modeMeta.cumEst, cumEst) || S.modeMeta.cumWin ~= cumWin
        warning('Fig219:est', 'saved selection used %s / window %d; SubV now comes from %s / window %d', ...
                S.modeMeta.cumEst, S.modeMeta.cumWin, cumEst, cumWin);
    end
    fprintf('Loaded %s\n', modeFile);
else
    ModeSel = repmat(modeTmpl, 1, max(foi));
end
if numel(ModeSel) < max(modeFOI), ModeSel(max(modeFOI)) = modeTmpl; end
cmapPC = gen_colormap([0 0.5 1; 1 1 1; 1 0 0], 256);

ii = 1;
if ~modeRedo
    ii = find(~[ModeSel(modeFOI).done], 1);
    if isempty(ii), ii = numel(modeFOI) + 1; disp('All neurons already selected (set modeRedo = true to review).'); end
end
while ii <= numel(modeFOI)
    f = modeFOI(ii);
    V = SubV{f};  D = SubD{f};
    nShow   = min(modeNShow, nnz(D > 0));
    varFrac = D / sum(D(D > 0));
    sgn = sign(sum(V(:, 1:nShow), 1));  sgn(sgn == 0) = 1;
    Vs  = V(:, 1:nShow) .* sgn;                          % display sign: loadings sum >= 0
    uIdx = abs(sum(Vs, 1)) ./ (sqrt(size(Vs,1)) * vecnorm(Vs));

    % starting selection: saved one, else the most uniform PC as global if u > 0.7
    sel = zeros(1, nShow);
    if ModeSel(f).done
        sel(ModeSel(f).globalPC(ModeSel(f).globalPC <= nShow)) = 1;
        sel(ModeSel(f).seesawPC(ModeSel(f).seesawPC <= nShow)) = 2;
    else
        [uMax, kMax] = max(uIdx);
        if uMax > 0.7, sel(kMax) = 1; end
    end

    Ft   = PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist) > 0;
    bnd  = max(Ft, [], 3);
    swcP = zeros(size(PCresult{f}.swc,1), 4);           % [row col ~ markersize], as in Fig 2F
    swcP(:,1:3) = PCresult{f}.swc;
    swcP(:,4) = 6;  swcP(1,4) = 30;                      % first point = soma
    fig  = figure(219); clf(fig);
    fig.UserData = struct('sel', sel, 'action', '');
    fig.KeyPressFcn     = @pcmode_key;
    fig.CloseRequestFcn = @(h,~) pcmode_close(h);
    tl = tiledlayout(fig, 3, ceil(modeNShow/3), 'TileSpacing', 'tight', 'Padding', 'compact');
    for k = 1:nShow
        ax = nexttile(tl, k);
        m  = max(abs(Vs(:,k)));
        if strcmp(modeView, 'dendrite')
            fig.CurrentAxes = ax;
            showScaleScatter(Vs(:,k), swcP, Ft, cmapPC, [-m m]);
            set(ax.Children, 'HitTest', 'off');
        else
            RGB = showScaleImageSmooth(Ft, Vs(:,k), cmapPC, [-m m]);
            image(ax, RGB .* bnd, 'HitTest', 'off');  axis(ax, 'image');  axis(ax, 'off');
        end
        hold(ax, 'on');
        xl = ax.XLim + [-1 1]*0.04*diff(ax.XLim);  yl = ax.YLim + [-1 1]*0.04*diff(ax.YLim);
        ax.XLim = xl;  ax.YLim = yl;                     % margin so the frame clears the trace
        ax.UserData = rectangle(ax, 'Position', [xl(1) yl(1) diff(xl) diff(yl)], ...
                                'LineWidth', 4, 'EdgeColor', 'none', 'HitTest', 'off');
        title(ax, sprintf('PC%d  %.1f%%  u=%.2f', k, 100*varFrac(k), uIdx(k)), 'FontSize', 9);
        ax.PickableParts = 'all';
        ax.ButtonDownFcn = @(a,~) pcmode_click(a, k);
        ax.Tag = 'pcmode';
    end
    title(tl, sprintf(['f=%d (%d/%d)   left = GLOBAL (red), right = SEE-SAW (blue)   ' ...
                       'Enter: next   b: back   c: clear   q: save & stop'], f, ii, numel(modeFOI)));
    pcmode_redraw(fig);
    uiwait(fig);
    act = fig.UserData.action;  sel = fig.UserData.sel;

    if ~strcmp(act, 'quit') || any(sel)          % record unless quitting with nothing picked
        gPC = find(sel == 1);  sPC = find(sel == 2);
        ModeSel(f).globalPC  = gPC;               ModeSel(f).seesawPC  = sPC;
        ModeSel(f).varGlobal = sum(varFrac(gPC)); ModeSel(f).varSeesaw = sum(varFrac(sPC));
        ModeSel(f).varFrac   = varFrac(1:nShow);  ModeSel(f).uIdx      = uIdx;
        ModeSel(f).globalVec = Vs(:, gPC);        ModeSel(f).seesawVec = Vs(:, sPC);
        ModeSel(f).done      = true;
        fprintf('  f=%-3d global PC [%s] = %5.1f%%   see-saw PC [%s] = %5.1f%%   (top %d PCs = %.1f%%)\n', ...
                f, num2str(gPC), 100*ModeSel(f).varGlobal, num2str(sPC), 100*ModeSel(f).varSeesaw, ...
                nShow, 100*sum(varFrac(1:nShow)));
    end
    modeMeta = struct('cumEst', cumEst, 'cumWin', cumWin, 'nShow', modeNShow, 'date', datestr(now));
    save(modeFile, 'ModeSel', 'modeMeta');

    switch act
        case 'quit', break;
        case 'back', ii = max(1, ii - 1);
        otherwise,   ii = ii + 1;
    end
end
if exist('fig', 'var') && ishandle(fig), fig.CloseRequestFcn = 'closereq'; fig.KeyPressFcn = ''; end

% ---- summary over every neuron with a selection ----
fDone = modeFOI([ModeSel(modeFOI).done]);
vG = [ModeSel(fDone).varGlobal];  vS = [ModeSel(fDone).varSeesaw];
fprintf('\nGlobal / see-saw modes (%s covariance, window %d), saved to %s\n', cumEst, cumWin, modeFile);
fprintf('  %-5s %-10s %9s %-10s %9s\n', 'f', 'global PC', 'var (%)', 'seesaw PC', 'var (%)');
for f = fDone
    fprintf('  %-5d %-10s %9.1f %-10s %9.1f\n', f, num2str(ModeSel(f).globalPC), 100*ModeSel(f).varGlobal, ...
            num2str(ModeSel(f).seesawPC), 100*ModeSel(f).varSeesaw);
end
fprintf('  Global : %.1f +/- %.1f %% (mean +/- std, n = %d neurons)\n', 100*mean(vG), 100*std(vG), numel(vG));
hasS = ~cellfun(@isempty, {ModeSel(fDone).seesawPC});
fprintf('  See-saw: %.1f +/- %.1f %% (n = %d neurons with a see-saw PC; %d without)\n', ...
        100*mean(vS(hasS)), 100*std(vS(hasS)), nnz(hasS), nnz(~hasS));
fprintf('  Global + see-saw: %.1f +/- %.1f %%\n', 100*mean(vG + vS), 100*std(vG + vS));

%% ===== Figure 125 — silent correlation of example ROIs on morphology =====
f = 20;
figure(125); clf; tiledlayout('flow');
cmap_corr = gen_colormap([0 0.5 1; 1 1 1; 1 0 0]);
Ft = PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist);
for j = [2 7 20 33]   % example seed ROIs
    nexttile([1 1]);
    showScaleImageSmooth(Ft > 0, corrMat{f,2}(j,:), cmap_corr, [-0.25 1]); hold all
    bd = cell2mat(bwboundaries(imdilate(Ft(:,:,j) > 0, strel('diamond',5))));
    plot(bd(:,2), bd(:,1), 'color', [1 1 1]);
end
cb = colorbar;  colormap(cmap_corr);  cb.Ticks = [0 1];  cb.TickLabels = [-0.4 1];

%% ===== Figure 126 — spatial decay of correlation with distance =====
figure(126); clf; tiledlayout(1,2);
nexttile([1 1]); hold on;   % pooled correlation vs distance + per-session fits
dat2show_cat=[]; gg=1; Lc2show=[]; amp2show=[];
for f = foi
    Dpair   = PCresult{f}.interDendDist(PCresult{f}.sortDist, PCresult{f}.sortDist);
    Dpair   = max(cat(3,Dpair,Dpair'),[],3);
    triuIdx = triu(true(size(Dpair)));   % upper triangle incl. diagonal (self)
    %xd = Dpair(triuIdx);  yc = corrMat{f,2}(triuIdx);
    xd = T.Xdata{row}(T.FitManual{row}.idxSel);  yc = T.Ydata{row}(T.FitManual{row}.idxSel);
    ok = ~isnan(xd) & ~isnan(yc);
    plot(xd(ok), yc(ok), '.', 'color', [0.85 0.85 0.85], 'markersize', 10);
    row = find(T.FitType=="spatial" & T.ClassPair=="All, All (auto)" & T.Condition=="Silent" & T.R2_manual>0.4 & T.FileInd==f, 1);
    Lc2show=[Lc2show T.Tau_manual(row)];
    amp2show=[amp2show T.A_manual(row)];
    dat2show_cat{gg}=[xd yc];
    gg=gg+1;
end
binResult=binning_data(dat2show_cat,[[0:30:200]-15 400 500]);
errorbar_shade(binResult.centers,binResult.mean,binResult.sem,[1 0 0]);
xlabel('Distance (µm)'); ylabel('Correlation');
title(sprintf('\\lambda: %2.0f ± %2.0f µm',mean(Lc2show),std(Lc2show)))
meanfit=@(x) mean(amp2show)*exp(-x/mean(Lc2show));
plot([1:500],meanfit(1:500),'k--','linewidth',2);
xlim([0 450]); ylim([0 1]); box off;

setLbl   = title_str;                 % dc==5 ('All, All (auto)') = pooled autocorrelation
cmap_temp = gen_colormap(Plasma, numel(setLbl));
T = CorrelationFitResult;

firstTile = true;
for k = 5   % all-classes autocorrelation; use 1:numel(setLbl) to show every pair
    cc = squeeze(temporalCorr(foi, heatCond, k));   % silent condition
    valid = find(~cellfun(@isempty, cc));
    if isempty(valid), continue; end

    lags      = lagsCorr{foi(valid(1))};
    Corr2show = cell2mat(cc(valid));                   % nSess x nLags
    pos       = lags >= 0;                             % positive side only

    nexttile; hold on;
    plot(lags(pos), Corr2show(:,pos), 'color', [0.75 0.75 0.75]);
    mCorr = mean(Corr2show(:,pos), 1, 'omitnan');
    sCorr = std(Corr2show(:,pos), 0, 1, 'omitnan');
    errorbar_shade(lags(pos), mCorr, sCorr./sqrt(size(Corr2show,1)), [0.5 0 1]);

    % overlay the interactively refined mean fit and annotate tau
    m  = T.FitType=="temporal" & T.ClassPair==string(setLbl{k}) & T.Condition=="Silent";
    tt = rmmissing(cell2mat(cellfun(@(x) x.tau(1), T.FitManual(m),'UniformOutput',0))); 
    tt2 = rmmissing(cell2mat(cellfun(@(x) x.tau(2), T.FitManual(m),'UniformOutput',0)));
    aa = rmmissing(cell2mat(cellfun(@(x) x.amp(1), T.FitManual(m),'UniformOutput',0))); 
    aa2 = rmmissing(cell2mat(cellfun(@(x) x.amp(2), T.FitManual(m),'UniformOutput',0))); 
    title(sprintf('\\tau_{fast} = %.0f ± %.0f ms,\n \\tau_{slow} = %.0f ± %.0f ms', mean(tt), std(tt), mean(tt2), std(tt2)));
    hold all;
    meanfit2=@(x) mean(aa)*exp(-x/mean(tt))+mean(aa2)*exp(-x/mean(tt2));
    plot([1:500],meanfit2(1:500),'k--','linewidth',2);
    xlabel('Lag (ms)'); xlim([0 300]); box off;
    ylim([0 1]);
    if firstTile, ylabel('Correlation'); firstTile = false; end
end
set_font('Arial'); set_fontsize(12);
set_figsize(185,90);

%% ===== Figure 127 — spatiotemporal correlation heatmap (distance x lag) =====
% Bin every ROI pair by dendritic distance (dendbin) with binning_data, then
% average its correlogram within each distance bin -> mean correlation as a
% function of (distance, lag), pooled over all sessions.
dendbin = 0:20:500;                       % dendritic-distance bin edges (µm)
lags    = lagsCorr{foi(find(~cellfun(@isempty, lagsCorr(foi)), 1))};
nLagH   = numel(lags);

distCells = cell(numel(foi),1);           % {session}: [nPairs x 1] pair distances
corrCells = cell(numel(foi),1);           % {session}: [nPairs x nLag] correlograms
for fi = 1:numel(foi)
    f = foi(fi);
    if isempty(temporalCorrMat{f}), continue; end
    C  = double(temporalCorrMat{f});
    nR = size(C,1);
    D  = PCresult{f}.interDendDist(PCresult{f}.sortDist, PCresult{f}.sortDist);
    D  = max(cat(3,D,D'),[],3);
    up = triu(true(nR));                  % unique pairs incl. diagonal (distance 0)
    Cr = reshape(C, nR*nR, nLagH);
    distCells{fi} = D(up);
    corrCells{fi} = Cr(up, :);
end

% spatial binning via binning_data (membership), then average correlograms per bin
binRes = binning_data(cellfun(@(d) [d, zeros(size(d))], distCells, 'UniformOutput', false), dendbin);
heat = nan(numel(dendbin)-1, nLagH);
for b = 1:numel(dendbin)-1
    pooled = [];
    for fi = 1:numel(foi)
        mb = binRes.membership{b, fi};
        if isempty(mb), continue; end
        pooled = [pooled; corrCells{fi}(mb, :)];   %#ok<AGROW>
    end
    if ~isempty(pooled), heat(b,:) = mean(pooled, 1, 'omitnan'); end
end

% row-normalized version: each distance bin divided by its own value at lag 0
[~, i0]  = min(abs(lags));                % index of zero lag
heatNorm = heat ./ heat(:, i0);           % row-wise normalization by t = 0
heatNorm(heat(:, i0) <= 0 | isnan(heat(:, i0)), :) = NaN;   % drop unusable rows

% diffusive spread reference: d^2 = 2*Dcoef*t  ->  d = sqrt(2*Dcoef*|t|)
Dcoef  = 807;                             % µm^2/ms
tRef   = linspace(min(lags), max(lags), 501);
dRef   = sqrt(2*Dcoef*abs(tRef));

figure(127); clf; tl = tiledlayout(1,3);
cmap_st = gen_colormap([0 0.4 1; 1 1 1; 1 0 0], 256);

ax1 = nexttile; %posL = lags >= 0;         % positive lag side
imagesc(lags(:), binRes.centers, heat(:, :)); axis xy; hold on;
plot(tRef, dRef, 'k:', 'linewidth', 1.5);
colormap(ax1, cmap_st); caxis(ax1, [-0.05 0.25]);
cb = colorbar; %cb.Label.String = 'Correlation';
xlim(ax1, [min(lags) max(lags)]); ylim(ax1, [min(binRes.centers) max(binRes.centers)]);
xlabel('Lag (ms)'); ylabel('Distance (µm)');
title('Correlation');

ax2 = nexttile;
imagesc(lags(:), binRes.centers, heatNorm(:, :)); axis xy; hold on;
plot(tRef, dRef, 'k:', 'linewidth', 1.5);
colormap(ax2, cmap_st); caxis(ax2, [-0.2 1]);
cb2 = colorbar; %cb2.Label.String = 'Correlation (norm. to lag 0)';
xlim(ax2, [min(lags) max(lags)]); ylim(ax2, [min(binRes.centers) max(binRes.centers)]);
xlabel('Lag (ms)');
title(sprintf('Row-normalized (t=0), d^2 = 2\\cdot%d\\mum^2/ms\\cdott', Dcoef));

ax3 = nexttile; hold on;                   % iso-correlation lines only, no background image
lvlsC   = [0.3 0.6 0.8];                   % contour levels (fraction of the lag-0 value)
cmapLvl = gen_colormap([0 0.4 1; 1 0 1; 1 0 0], numel(lvlsC));   % low -> high level
hLvl    = gobjects(1, numel(lvlsC));
for ii = 1:numel(lvlsC)                    % one call per level so each line gets its own color
    [~, hLvl(ii)] = contour(lags(:), binRes.centers, heatNorm, lvlsC(ii)*[1 1], ...
                            'color', cmapLvl(ii,:), 'linewidth', 1.5);
end
hRef = plot(tRef, dRef, 'k:', 'linewidth', 1.5);
set(ax3, 'color', 'none', 'box', 'off');
legend([hLvl hRef], [compose('%.1f', lvlsC(:)'), {'d^2 = 2Dt'}], ...
       'box', 'off', 'location', 'northeast');
xlim(ax3, [min(lags) max(lags)]); ylim(ax3, [min(binRes.centers) max(binRes.centers)]);
xlabel('Lag (ms)');
title('Iso-correlation contours');
linkaxes([ax1 ax2 ax3],'xy');
ylim([0 350]); xlim([-400 400]);

set_font('Arial'); set_fontsize(15);
set_figsize(345,90);
%% Peak triggered average (Figure 2F–G)

if isfile(fullfile(save_figto, 'PeakTriggerAverage_Fig2_results.mat'))
    PeakSTA=[]; peakvec=[]; TroughSTA=[]; troughvec=[]; PeakTrough_FR=[];
    PeakSTA_ch=[];
    nTauPeak=[1000 1000];
    PeakTroughMat=[]; PeakSTASpClass=[]; TroughSTASpClass=[];
    for f=[foi]
        PeakSTA{f}=[];
        for ch=1:2
            PeakSTA_ch{f,ch}=[];
        end
        TroughSTA{f}=[];
        StopFreq=[12 200]; %filter high frequency
        [nROI nTime]=size(PCresult{f}.Subthreshold);
        Subthreshold_int=interpolateNaN(PCresult{f}.Subthreshold/PCresult{f}.SpikeHeight);
        FilterTr=[];
        for n=1:nROI
            %[PhaseTrace BasalSubFilt BasalthetaPower] = get_phase(Subthreshold_int(n,:), 1000, FilterFreq);
            perispikefrm=unique([tovec(find(PCresult{f}.SpikeMat(1,:))'+[-10:40]); find(PCresult{f}.CStrace)']);
            perispikefrm(perispikefrm<1 | perispikefrm>nTime)=[];
            FilterTr(n,:)=get_bandstop(Subthreshold_int(n,:),1000,StopFreq);
            FilterTr(n,:)=FilterTr(n,:)-movmedian(FilterTr(n,:),300,2);
            [pks, locs] = findpeaks(FilterTr(n,:),'MinPeakHeight', 0.3,'MinPeakDistance',50, ...
                'MinPeakProminence', 0.1);  % Prominence can also be tuned
            [trough, loc_trgh] = findpeaks(-FilterTr(n,:),'MinPeakHeight', 0.3,'MinPeakDistance',50, ...
                'MinPeakProminence', 0.1);  % Prominence can also be tuned
            % Turn to 0 during blue Stim and peri-spike frame
            peakvec{f}(n,:)=ind2vec(nTime,locs,1);
            peakvec{f}(n,ind2vec(nTime,perispikefrm,1)>0)=0;
            troughvec{f}(n,:)=ind2vec(nTime,loc_trgh,1);
            troughvec{f}(n,ind2vec(nTime,perispikefrm,1)>0)=0;

            if ~isempty(find(peakvec{f}(n,:)>0))
                [PeakSTA{f}(:,:,n), peakMat, peakloc]=get_STA(PCresult{f}.Subthreshold,peakvec{f}(n,:),nTauPeak(1),nTauPeak(2));
                [~, PeakSTASpClass{f,n}]=get_STA(PCresult{f}.SpikeClassMat,peakvec{f}(n,:),nTauPeak(1),nTauPeak(2));
                for ch=1:2
                    [PeakSTA_ch{f,ch}(:,:,n), ~, ~]=get_STA(PCresult{f}.Subthreshold_ch{ch},peakvec{f}(n,:),nTauPeak(1),nTauPeak(2));
                end
            else
                PeakSTA{f}(:,:,n)=NaN(nROI,sum(nTauPeak)+1);
                for ch=1:2
                    PeakSTA_ch{f,ch}(:,:,n)=NaN(nROI,sum(nTauPeak)+1);
                end
                peakMat=NaN(nROI,sum(nTauPeak)+1);
                peakloc=[];
                PeakSTASpClass{f,n}=NaN(nROI,sum(nTauPeak)+1);
            end

            if ~isempty(find(troughvec{f}(n,:)>0))
                [TroughSTA{f}(:,:,n), troughMat, troughloc]=get_STA(PCresult{f}.Subthreshold,troughvec{f}(n,:),nTauPeak(1),nTauPeak(2));
                [~, TroughSTASpClass{f,n}]=get_STA(PCresult{f}.SpikeClassMat,troughvec{f}(n,:),nTauPeak(1),nTauPeak(2));
            else
                TroughSTA{f}(:,:,n)=NaN(nROI,sum(nTauPeak)+1);
                troughMat=NaN(nROI,sum(nTauPeak)+1);
                troughloc=[];
                TroughSTASpClass{f,n}=NaN(nROI,sum(nTauPeak)+1);
            end
        end
        disp(['Peak triggered average finished, file:' num2str(f)])
    end

    save(fullfile(save_figto, 'PeakTriggerAverage_Fig2_results.mat'), ...
        'PeakSTASpClass', 'troughvec', 'peakvec', 'nTauPeak', 'PeakSTA_ch', ...
        'TroughSTA', 'PeakSTA', 'TroughSTASpClass', '-v7.3');
else
    load(fullfile(save_figto, 'PeakTriggerAverage_Fig2_results.mat'))
end
%% Peak triggering average kymo (Figure 2F)
figure(132); clf; tiledlayout(2,1);
showf=[23]; cax=[-0.01 0.03];
[nROI, nTime]=size(PCresult{showf}.Subthreshold);
showInd=[3 8 39];
showSTA=PeakSTA{showf}-median(PeakSTA{showf}(:,1:1000,:),2);
for i=[1 3]
    nexttile([1 1]);
    [dendaxissorted, dsort]=sort(PCresult{showf}.dendaxis);
    imagesc([-nTauPeak(1):nTauPeak(2)],[1:nROI],showSTA(dsort,:,showInd(i)),cax);
    set_kymoYtick(dendaxissorted);
    xlabel('Time (ms)'); %ylabel('Distance from soma (µm)');
    xlim([-500 500]);
    cb=colorbar;
    cb.Label.String='∆F/F';
end
colormap(turbo);
set_fontsize(18); set_font('Arial');
set_figsize(160,220);

SWCpoints=zeros(size(PCresult{showf}.swc,1),4);
SWCpoints(:,1:3)=PCresult{showf}.swc;
SWCpoints(:,4)=SWCpoints(:,4)+6; SWCpoints(1,4)=30;
ftprints=PCresult{showf}.Ftprnts(:,:,PCresult{showf}.sortDist);
[H, W, T] = size(ftprints);
t2show=[-120:60:120]; shift_x=220;

figure(133); clf;
tiledlayout(2,1,'TileSpacing','tight');
for i=[1 3]
    nexttile([1 1]);
    [dendaxissorted, dsort]=sort(PCresult{showf}.dendaxis);
    tax=[-nTauPeak(1):nTauPeak(2)];
    g=1;
    for t=1:length(t2show)
        peakpattern=showSTA(dsort,find(tax==t2show(t)),showInd(i));
        SWCpoints2show=SWCpoints+[0 g 0 0]*shift_x;
        newH = H + shift_x*g;
        ftprints_shifted = zeros(newH, W, T);
        ftprints_shifted(shift_x*g+1 : shift_x*g+H,:, :) = ftprints;
        showScaleScatter(peakpattern,SWCpoints2show,ftprints_shifted(:,:,dsort),'turbo',cax);
        g=g+1;
    end
    if t==length(t2show) & i==3
        drawScaleBar(100*Pixelsize(showf),'vertical','color','k','position',[shift_x*(g-0.1) 50],'Linewidth',3);
    end
end
set_figsize(250,220);

%% Peak triggering average by ROIs (Figure 2G)
PeakTraceROI=[];
figure(134); clf; tiledlayout(2,2,'TileSpacing','compact'); cmap_temp=gen_colormap(Plasma,4); 
title_str={'Basal peak, Basal voltage','Apical peak, Apical voltage','Apical peak, Basal voltage','Basal peak, Apical voltage'};
dClass2temporal = [[{1},{1}];[{[3 5]} {[3 5]}];[{[1]} {[5]}];[{[5]} {[1]}]]; % Reading and triggering
ax1=[]; g=1; cmap_seq=[1 4 1 4];
for dc=[1 4 3 2]
    for f=foi
        ROI2plot=ismember(PCresult{f}.labelvec,dClass2temporal{dc,2})'*ismember(PCresult{f}.labelvec,dClass2temporal{dc,1});
        % uppertri_ind=1-triu(ones(length(labelvec{f})));
        % ROI2plot=ROI2plot & uppertri_ind;
        STAsub=permute(PeakSTA{f}./PCresult{f}.SpikeHeight,[3 1 2]); % triggering ROI, measuring ROI, peri-peak time
        STAsub=STAsub-median(STAsub,3);
        PeakTraceROI{f}(:,dc)=tovec(STAsub)'*tovec(ROI2plot)/sum(ROI2plot(:));
    end

    ax1=[ax1 nexttile([1 1])];
    PeakTrace2show=cell2mat(cellfun(@(x) x(:,dc),PeakTraceROI(foi),'UniformOutput',false));
    plot([-nTauPeak(1):nTauPeak(2)],PeakTrace2show,'color',[0.7 0.7 0.7]); hold all;
    errorbar_shade([-nTauPeak(1):nTauPeak(2)],mean(PeakTrace2show,2,'omitnan'),std(PeakTrace2show,0,2,'omitnan'),cmap_temp(cmap_seq(g),:));
    xlabel('Time (ms)'); ylabel('Subthreshold (spike height)');
    title(title_str{dc})
    box off;
    g=g+1;
end
linkaxes(ax1,'x'); xlim([-500 500]);
set_font('Arial'); set_fontsize(20);
set_figsize(300,270)

%% ===== ALL PEAKS — one unified scoring pass, NO THRESHOLDS =====
% Replaces the four cells that used to live here (L50 spatial extent / oscillation
% score / suprathreshold length & duration / scatter viewer). Those pipelines each
% picked their own amplitude window, applied their own cutoffs and threw events away
% before anything could be checked; the result was hard to follow and the cutoffs were
% doing the discriminating rather than the data.
%
% Here EVERY peak in peakvec is scored with EVERY parameter and NOTHING is rejected on
% the basis of its VALUES: no amplitude band, no length/duration/oscillation cutoff, no
% spike/blue exclusion, no hand-drawn polygons. Peri-spike proximity is recorded as a
% COLUMN (nearSpk) so it can be tested as a discriminator like anything else, instead of
% being silently removed.
%
% The ONE thing that is collapsed is multiple counting of the SAME event. A subthreshold
% depolarization spreads over many ROIs, so peakvec fires once per ROI for what is a
% single event and every parameter would enter the table several times over. Within a
% dedupWin_PK ms window only the LARGEST peak (the ROI with the highest window-mean
% subthreshold voltage — the very amplitude the parameters are measured from) is kept
% and the co-active ROIs are dropped. This is a de-duplication of one event, not a
% quality filter: nothing is judged, and peaks farther apart than the window are
% untouched, so genuinely separate events are never merged.
%
% Length and duration are measured the established way: the extent over which the peak
% stays above Vthr_sh = 0.15 spike height. A peak that never reaches 0.15 gets 0, not
% NaN, because the region above 0.15 is then genuinely empty — that keeps the measure
% defined for every peak instead of silently dropping the small ones. NOTE this measure
% is amplitude-coupled by construction (a bigger peak clears 0.15 for longer and
% further), so lenSet / durSet will correlate with amp; keep that in mind when reading
% the classifier's coefficients. The old relative version measured at 50% of each
% event's OWN amplitude, so they are defined for every peak however small. The old
% absolute-threshold versions (lenSet / durSet at Vthr_sh) are kept alongside them for
% comparison — they are NaN for any peak below Vthr_sh, which is most of them, and that
% is the point.
%
% Workflow:
%   this cell        -> PkTab (nPeak x nParam) + PkMeta (f, ROI, frame)
%   next cell        -> label a random sample of peaks by eye
%   cell after that  -> ask the labels which parameter actually discriminates
%
% Kept for the downstream cells that still expect them: exampleList / nEx / exAMP /
% exL50 / exCOMP / exR0 / exProf / exAMP_dff, and an EventTab shim at the end.
% Helpers at the end of the file: event_extent, event_len_dur, event_peak_stats,
% event_local_psd, peak_classifier, param_auc.

% --- TEST RUN ----------------------------------------------------------------
% Scoring every neuron is the slow step, so do a test pass on one file first: set
% foi_PK to a single index (e.g. 20), check the printout and the figures, then set it
% back to foi for the real run. maxPkPerFile_PK additionally caps how many peaks are
% scored per neuron, which turns the pass into a few-second smoke test.
% Everything downstream keys off runTag_PK, so a partial run writes its labels to
% PeakLabels_f20.mat and CANNOT overwrite the labels from a full run.
foi_PK = foi;                    % <-- set to foi for the full run
maxPkPerFile_PK = Inf;          % Inf = every peak; e.g. 2000 for a quick smoke test
isTestRun_PK = ~isequal(sort(foi_PK(:)'), sort(foi(:)')) || isfinite(maxPkPerFile_PK);
if isTestRun_PK
    if numel(foi_PK) > 4
        runTag_PK = sprintf('n%dfiles', numel(foi_PK));
    else
        runTag_PK = ['f' strjoin(arrayfun(@num2str, foi_PK(:)', 'UniformOutput',false), '_')];
    end
    if isfinite(maxPkPerFile_PK)
        runTag_PK = [runTag_PK sprintf('_max%d', maxPkPerFile_PK)];
    end
    fprintf(['*** TEST RUN: %d of %d neurons%s. Label/ranking/classifier files are ' ...
        'tagged "%s". The EventTab shim and every cell below the classifier will be ' ...
        'PARTIAL — do not read them as results. ***\n'], numel(foi_PK), numel(foi), ...
        repmat_str(isfinite(maxPkPerFile_PK), sprintf(', <=%d peaks each', maxPkPerFile_PK), ''), ...
        runTag_PK);
else
    runTag_PK = 'all';
end
% Labels are keyed by peak IDENTITY (file, ROI, frame), never by row number, so a label
% file written by ANY run can be read back by any other. That is what makes a test run
% useful: if you already labelled peaks on the full set, the ones belonging to the
% neurons in foi_PK are found and carried over automatically. It works in reverse too —
% point importFrom_PK at a test run's file when you later score everything.
labFile_PK    = fullfile(save_figto, ['PeakLabels_' runTag_PK '.mat']);
importFrom_PK = fullfile(save_figto, 'PeakLabels_all.mat');   % {} or '' to disable

% --- geometry / amplitude ---
avgwin_ext = 10;      % +/- ms averaged at the peak for every window amplitude
Dmax_ext   = 300;     % um: ROIs farther than this from the trigger are ignored
Vthr_sh    = 0.15;    % spike-height units: LENGTH and DURATION are the extent over
                      % which the peak stays above this level — the measure already in
                      % use. Absolute, not a fraction of each event's own amplitude, so
                      % lengths are comparable across events in real voltage terms.
baseWin_PK = 1000;    % +/- ms local median baseline (the event window excluded)
durMax_PK  = 100;     % +/- ms cap on the duration search
smoothDur_PK = 5;     % movmean (ms) on the trace before the duration search. When
                      % bpLenDur_PK is on this is FORCED TO 1: a 60 Hz low-pass and a
                      % 5 ms moving average are two low-passes in series, and stacking
                      % them broadens every duration. The effective value is printed.
% --- the trace LENGTH and DURATION are measured on --------------------------------
bpLenDur_PK = [0.1 60];  % Hz, zero-phase bandpass applied before the threshold search.
                      % [] = off, i.e. the unfiltered trace the old code used.
                      % The high corner removes the 40-250 Hz blood noise that makes a
                      % raw trace cross 0.15 spike height for a millisecond at a time;
                      % the low corner removes the slow drift that shifts the whole
                      % trace relative to an ABSOLUTE threshold. Both directly bias a
                      % threshold-crossing measurement, which is why it is worth doing
                      % here even though amp / bandRatio / oscRatio deliberately keep
                      % using the unfiltered trace (psdP MUST see 40-250 Hz).
                      % Widening this changes lenSet and durSet, so re-score and say so.
bpOrd_PK   = 2;       % Butterworth order per band edge, applied as second-order
                      % sections with filtfilt (zero phase, so no timing shift: a
                      % duration measured on the filtered trace is still centred on the
                      % peak). Order 2 is deliberate — a steeper filter rings, and a
                      % ringing filter invents threshold crossings.
spkWin_PK  = 30;      % +/- ms around somatic spikes / CS frames -> nearSpk (a COLUMN)
dedupWin_PK = 50;     % ms: FULL width of the one-event-one-peak window (so +/- 25 ms).
                      % Peaks in different ROIs closer in time than this belong to the
                      % same spreading depolarization, so only the biggest one is scored.
                      % Set to 0 to score every ROI's peak separately (old behaviour).
% --- blood noise (fast) vs event shape (slow): TWO timescales, kept separate ------
% This is the one thing not to mix up. Blood noise is fast (40-250 Hz) and the events
% we care about are slow (tens of ms), so the two families of measurement need OPPOSITE
% filtering and must not share settings:
%   psdP  works on the UNSMOOTHED raw trace and only ever looks inside fVess_PK.
%         -> bandRatio, bandFrac, fLoc
%   oscP  works on a trace smoothed to the EVENT timescale, which deliberately
%         suppresses the blood band so that oscRatio compares the peak against slow
%         fluctuation rather than against high-frequency noise.
%         -> amp, oscRatio, cmCorr
% Deriving oscP.smooth from fVess_PK (which is what I tried first) breaks this: with a
% 1 ms smooth, oscAmp is dominated by the 40-250 Hz noise itself, and a real 0.35
% event sitting in a vessel-contaminated ROI scored oscRatio 1.24 against 0.93 for pure
% vessel noise — i.e. no separation at all. With the 5 ms smooth below the same event
% scores well clear of it.
fVess_PK = [40 250];  % Hz: blood / red-blood-cell flow noise. THE only frequency knob.
                      % Set [2 12.5] instead to test the cardiac-pulsation hypothesis.
oscP = struct('smooth',5, 'Wosc',300, 'Wac',700, 'Wplot',400);   % ms, EVENT timescale
psdP = struct();
psdP.Wnar  = 100;     % +/- ms narrow window = the peak itself
psdP.Wwide = 2000;    % +/- ms wide window   = the background it sits in
psdP.nfft  = 2^12;    % common zero-padded grid, 0.24 Hz bins
psdP.fVess = fVess_PK;
psdP.fBand = [max(1,fVess_PK(1)/4) min(450,2*fVess_PK(2))];  % fLoc search range, kept
                      % wider than fVess_PK so "fLoc lands in band" stays a real check
fprintf(['Blood band %g-%g Hz (psdP, raw trace) | event timescale %g ms smooth (oscP) ' ...
    '| fLoc searched %g-%g Hz\n'], fVess_PK(1), fVess_PK(2), oscP.smooth, ...
    psdP.fBand(1), psdP.fBand(2));
bvDil_PK = 3;         % px dilation of bvMask when measuring ROI overlap
% --- the ROI-SET length and duration (lenSet / durSet) ----------------------------
% Every ROI above Vthr_sh at the peak instant is collected, the above-threshold ROIs
% are grouped into sets connected along the real dendrite, and lenSet/durSet describe
% the set the detecting ROI belongs to. See event_roi_sets for the full rule; the
% dedicated ROI-SET cell further down shows individual events and their maps.
ciZ_PK     = 1.96;    % the CI a single separating ROI is judged against. The gap
                      % between two sets is believed only when Vthr_sh is above that
                      % ROI's A + ciZ*sigma, i.e. even its plausible maximum stays
                      % below threshold. 1.96 = 95%, 2.58 = 99%.
repTrace_PK = 'peak'; % the set's representative trace, which BOTH height and duration
                      % are read from: 'peak' (largest member) or 'center' (1-centre)
adjMax_PK  = 60;      % µm cap on branch-point links in the ROI graph ([] = no cap).
                      % Intra-branch links are never capped — the dendrite is
                      % continuous between consecutive ROIs however far apart they are.
alsoRaw_PK = true;    % also measure lenSet/durSet on the UNFILTERED trace, as a
                      % control that the band-pass is not creating the extent. Costs
                      % one extra set computation per peak; set false to skip it and
                      % leave the lenSetRaw/durSetRaw/nROIsetRaw columns NaN.
dtDirtMax_PK = 1000;  % ms: cap on dtDirt. Beyond this the exact distance stops
                      % meaning anything, and an uncapped Inf would break any
                      % correlation or model that uses the column.
useFastBase_PK = true;% SPEED. true = the local baselines are precomputed for every
                      % frame with movmedian, once per neuron, instead of a median per
                      % peak. ~3x faster overall and the only change that is not
                      % bit-exact: movmedian cannot cut each event out of its own
                      % baseline window (measured effect <= 4e-4 spike height on a
                      % 2001-frame window). Set false to restore the exact per-peak
                      % medians — everything else in the speed-up is exact either way.
% --- area under the peak (aucSum / aucArea) ---
aucWin_PK  = [-50 50];  % ms relative to the peak: the integration window. Wide enough
                        % to contain a slow subthreshold event, narrow enough that the
                        % integral is not dominated by baseline wander.
ftpThr_PK  = 0.2;       % fraction of each footprint's own maximum weight that counts as
                        % "inside the ROI" when its area is measured. The footprints are
                        % weighted, not binary, so an area needs SOME cut; 0.2 is where
                        % the weights have fallen to the background. aucArea scales
                        % roughly as 1/area, so changing this rescales the parameter —
                        % it does not reorder the peaks.
% --- hand-picked local candidates (kept: they are the closest thing to ground truth) ---
ampTol_ext = 0.25;    % used by the downstream amplitude-matched control cell
nCtrl_ext  = 100;     % same
exampleList = [26 493523; 26 539246; 26 562118; 26 576267; 18 422263];  % [fileInd frame]
candROIman  = [17; 17; 17; 17; 12];   % trigger ROI per candidate, sortDist numbering
nEx = size(exampleList,1);

% --- THE PARAMETER SET ------------------------------------------------------
% Twelve parameters, one per distinct physical question. This used to be 37; knowing
% the blood-noise band is 40-250 Hz collapsed most of them, because the whole
% machinery that ESTIMATED an unknown oscillation frequency (fLoc gating, pkProm,
% statScore, fMatch, fACF, acMin/acLag/acReb) existed only to find the band. With the
% band given, in-band power is measured directly and all of that goes away. Dropped
% for redundancy at the same time: ampRaw / pkHeight (same thing as amp), L50 / nROI50
% (same thing as lenAbs), lenHalf / durHalf (the 50%-relative twins), asymHalf, dipSym / dipAsym / bipol / nZC (the
% slow down-up-down-up family — at 40-250 Hz that shape is a 4-25 ms ripple, which
% bandRatio measures far more directly, and on synthetic data dipSym alone threw away
% 60% of real events), oscAmp, cmFrac, hbPhase, compart, nNbrPk, nearBlue, censSpace,
% censTime. event_len_dur and event_extent still compute the dropped extras if you
% want any of them back.
pkNames = { ...
  'amp'       % 1  peak amplitude, spike-height units (subthreshold, +/-avgwin mean)
  'lenSet'    % 2  LENGTH of the ROI SET the detecting ROI belongs to: the full contour
              %    distance between the set's two farthest members (= 2 x lenHalf; the
              %    full length, not the half). The set is built by event_roi_sets: ROIs
              %    above Vthr_sh at this one instant, connected along the ACTUAL
              %    dendrite, with a single separating ROI able to split the set only
              %    when it is confidently below threshold (Vthr above the top of its
              %    CI). 0 = the set is one ROI, or the detecting ROI is not above
              %    threshold at all, so there is no extent to measure.
              %    This REPLACES lenAbs, which walked outward in sorted distance and
              %    was truncated by the first ROI below threshold — one noisy ROI could
              %    take a 50 µm event to 0.
  'durSet'    % 3  DURATION of that same set, measured on its representative trace
              %    (the largest member by default): ms above Vthr_sh. Height and
              %    duration therefore describe the same signal.
  'bandRatio' % 4  40-250 Hz power AT the peak (+/-100 ms) / same band in the background
              %    (+/-2 s). ~1 = the vessel noise is always there in this ROI;
              %    >>1 = a burst of flow noise landed here. THE blood-noise measure.
  'bandFrac'  % 5  share of the peak window's total power that is in 40-250 Hz:
              %    "is this peak MADE of blood noise?"
  'fLoc'      % 6  diagnostic: argmax of PSD(+/-2 s) - PSD(+/-100 ms) in Hz. Not a
              %    discriminator any more, kept so you can confirm the in-band power
              %    really sits where fVess_PK says it does.
  'oscRatio'  % 7  peak height / amplitude of the local fluctuation in the flanks.
              %    Band-agnostic local signal-to-fluctuation ratio: ~1 = the "peak" is
              %    no bigger than what fluctuates around it, >>1 = it stands out.
  'cmCorr'    % 8  correlation with the mean-over-ROI trace: motion and whole-field
              %    artifacts are common-mode, a local dendritic event is not
  'bvFrac'    % 9  trigger ROI footprint lying on the blood-vessel mask
  'roiPkRate' % 10 peaks per second detected in this ROI over the whole recording
              %    (a vessel-contaminated ROI fires constantly)
  'distSoma'  % 11 |contour distance| of the trigger ROI from the soma (um)
  'nearSpk'   % 12 within +/-spkWin_PK ms of a somatic spike or CS frame — a COVARIATE,
              %    recorded so it can be tested, never used to reject
  'aucSum'    % 13 SUM over ROIs of the area under the peak: for every ROI within
              %    Dmax_ext of the trigger, integrate the baseline-subtracted
              %    subthreshold trace over aucWin_PK and add the integrals up.
              %    Units: spike-height x ms. The total depolarization the event
              %    delivered to the imaged dendrite, amplitude x duration x spread in
              %    one number.
  'aucArea'   % 14 MEAN AUC per peak = aucSum / (total footprint area of those same
              %    ROIs), spike-height x ms / um^2. Dividing by area rather than by ROI
              %    count makes it an area-weighted mean, so a neuron imaged with a few
              %    big ROIs and one with many small ones are comparable. Signed and NOT
              %    rectified: an ROI the event never reached contributes ~0, so this is
              %    deliberately diluted by how much quiet dendrite is inside Dmax_ext —
              %    that dilution IS the "mean over the dendrite" the measure asks for.
              %    Report Dmax_ext and aucWin_PK with any value taken from it.
  'nROIset'   % 15 ROIs in the detecting ROI's set. 0 = it was not above threshold.
              %    1 = above threshold but alone, which is why lenSet is 0 there: a
              %    genuine measurement limited by the ROI spacing, NOT a failure.
  'nSetEv'    % 16 how many separate sets existed anywhere on the dendrite at this
              %    instant. >1 means other domains were simultaneously above
              %    threshold; they are NOT necessarily independent events, so this is
              %    a flag to look at the event, not a classification.
  'dtDirt'    % 17 ms from the peak to the nearest DIRT-masked frame IN ITS OWN ROI.
              %    Dirt was blanked to NaN in NormalizedTrace_dirt and then filled by
              %    interpolateNaN, so by the time anything is measured the trace looks
              %    continuous and the removal is invisible — a peak sitting next to, or
              %    on the shoulder of, an interpolated stretch is exactly where a
              %    spurious one is likely. This records how close that is.
              %    0 = the peak frame itself was dirt. Capped at dtDirtMax_PK, so the
              %    cap value means "no dirt within that window", including neurons
              %    with no dirt removal at all (ifdirtRemov false).
              %    Dirt is per ROI x frame; whole-frame motion rejection is a separate
              %    thing and is excluded from this measure.
  'lenSetRaw' % 18 CONTROL: lenSet measured on the UNFILTERED dF/F trace instead of
              %    the 0.1-60 Hz one. Same threshold, same adjacency, same gap rule —
              %    the only difference is the trace. If the two agree, the extent is
              %    not something the band-pass created; where they disagree, the raw
              %    trace is the one carrying slow drift and 40-250 Hz blood noise
              %    across an absolute threshold. NaN when alsoRaw_PK is false.
  'durSetRaw' % 19 the same control for durSet
  'nROIsetRaw'% 20 ROIs in the raw-trace set, for comparison with nROIset
  }';
nP = numel(pkNames);
iP = @(s) find(strcmp(pkNames,s), 1);      % column index by name

PkMeta = [];  PkTab = [];
nRawAll_PK = [];  nDedupAll_PK = [];   % peaks in peakvec before / after the dedup, per file
% length/duration are measured on the bandpassed trace, so the extra moving average
% would be a second low-pass in series — see the note next to smoothDur_PK
smoothDurEff_PK = smoothDur_PK;
if ~isempty(bpLenDur_PK), smoothDurEff_PK = 1; end
if isfinite(maxPkPerFile_PK), rng(0); end       % reproducible peak subsample
tic
for f = foi_PK(:)'                              % row vector: one iteration per neuron
    RawSH = interpolateNaN(PCresult{f}.NormalizedTrace_dirt / PCresult{f}.SpikeHeight);
    SubSH = interpolateNaN(PCresult{f}.Subthreshold      / PCresult{f}.SpikeHeight);
    Dmat  = PCresult{f}.interDendDist(PCresult{f}.sortDist, PCresult{f}.sortDist);
    Dmat  = max(cat(3,Dmat,Dmat'),[],3);
    lab   = PCresult{f}.labelvec(:);
    dsom  = abs(PCresult{f}.dendaxis(:));            % contour distance from the soma
    [nR, nT] = size(RawSH);
    Gl = mean(SubSH, 1, 'omitnan');                  % global / common-mode trace
    % the trace lenSet/durSet are read off. RawSH itself is left alone: amp,
    % bandRatio, bandFrac, fLoc, oscRatio and cmCorr all still come from the unfiltered
    % trace, and the PSD measures MUST (their band is above this filter's corner).
    [BpSH, bpTxt_PK] = bandpass_rows(RawSH, bpLenDur_PK, 1000, bpOrd_PK);
    % everything the ROI-set measurement needs, once per neuron
    MwPK  = movmean(BpSH, 2*avgwin_ext+1, 2, 'omitnan');   % the window-mean statistic
    sigPK = roi_noise(MwPK, [], baseWin_PK, avgwin_ext, 'global');
    % LOCAL BASELINES FOR EVERY FRAME, computed once instead of per peak. A median over
    % +/-baseWin for all ROIs costs ~0.8 ms per peak; movmedian over the whole
    % recording costs ~0.1 s for the entire neuron and answers every peak by lookup.
    % The one difference: movmedian cannot exclude each event's own +/-3*avgwin. Over a
    % 2001-frame window that changed the baseline by at most 4e-4 spike height in
    % testing, against a 0.15 threshold — set useFastBase_PK = false to go back to the
    % exact per-peak medians and confirm that on your own data.
    if useFastBase_PK
        baseMwPK  = movmedian(MwPK,  2*baseWin_PK+1, 2, 'omitnan');
        baseSubPK = movmedian(SubSH, 2*baseWin_PK+1, 2, 'omitnan');
        baseRawPK = movmedian(RawSH, 2*baseWin_PK+1, 2, 'omitnan');
        baseGlPK  = movmedian(Gl,    2*baseWin_PK+1, 2, 'omitnan');
    else
        baseMwPK = [];  baseSubPK = [];  baseRawPK = [];  baseGlPK = [];
    end
    blPK  = PCresult{f}.BranchLabel(:);  blPK = blPK(PCresult{f}.sortDist);
    [AdjPK, adjInfo] = dend_adjacency(blPK, abs(PCresult{f}.dendaxis(:)), Dmat, adjMax_PK);
    Pset  = struct('Vthr',Vthr_sh, 'ciZ',ciZ_PK, 'durMax',durMax_PK, ...
        'avgwin',avgwin_ext, 'repMode',repTrace_PK);
    % the same machinery on the UNFILTERED trace, as a control on the band-pass
    MwRaw = [];  sigRaw = [];  baseMwRaw = [];
    if alsoRaw_PK
        MwRaw  = movmean(RawSH, 2*avgwin_ext+1, 2, 'omitnan');
        sigRaw = roi_noise(MwRaw, [], baseWin_PK, avgwin_ext, 'global');
        if useFastBase_PK
            baseMwRaw = movmedian(MwRaw, 2*baseWin_PK+1, 2, 'omitnan');
        end
    end

    % footprint overlap with the blood-vessel mask, and each ROI's footprint AREA
    ftp = abs(PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist));
    bvF = nan(nR,1);
    if isfield(PCresult{f},'bvMask') && ~isempty(PCresult{f}.bvMask)
        bvm = imdilate(max(PCresult{f}.bvMask,[],3) > 0, strel('disk',bvDil_PK));
        bvF = squeeze(sum(sum(ftp .* double(bvm),1),2)) ./ ...
              max(squeeze(sum(sum(ftp,1),2)), eps);
    end
    aROI = roi_areas(ftp, PCresult{f}.pixelsize, ftpThr_PK);   % um^2 per ROI, for aucArea
    % distance to the nearest dirt-masked frame, per ROI, once per neuron
    motRej = [];
    if isfield(PCresult{f},'motionReject'), motRej = PCresult{f}.motionReject; end
    [dtDirtPK, dirtInfo] = dirt_distance(PCresult{f}.NormalizedTrace_dirt, motRej, ...
        dtDirtMax_PK);

    % somatic-spike / CS proximity — recorded, NOT used to reject
    badT = false(1,nT);
    badT(PCresult{f}.SpikeMat(1,:) > 0) = true;
    badT(PCresult{f}.CStrace > 0)       = true;
    badT = imdilate(badT, ones(1, 2*spkWin_PK+1));

    pv = peakvec{f} > 0;
    roiRate = sum(pv,2) / (nT/1000);                 % peaks per second, per ROI
    [r0all, t0all] = find(pv);
    nPkRaw = numel(r0all);

    % --- ONE EVENT, ONE PEAK: non-maximum suppression in time --------------------
    % peakvec fires in every ROI a depolarization reaches, so a single spreading event
    % appears as a cluster of peaks a few ms apart in neighbouring ROIs and would be
    % scored once per ROI. Keep only the largest peak of each cluster, where "largest"
    % is the window-mean subthreshold voltage at the peak — exactly the amplitude
    % event_extent reports, so the winner is the ROI the parameters would call the
    % trigger anyway. Greedy: the biggest peak in the file claims its +/- dedupWin_PK/2
    % ms window, the biggest peak outside every claimed window claims the next, and so
    % on, so a run of overlapping windows cannot chain into one long exclusion zone.
    % roiPkRate above is computed from ALL peaks and is unaffected by this.
    if dedupWin_PK > 0 && nPkRaw > 0
        SubWin = movmean(SubSH, 2*avgwin_ext+1, 2, 'omitnan');   % = mean over t0 +/- avgwin
        ampAll = SubWin(sub2ind([nR nT], r0all, t0all));
        [~, ordAmp] = sort(ampAll, 'descend');       % NaN amplitudes sort last
        halfDedup = round(dedupWin_PK/2);
        blocked = false(1, nT + 2*halfDedup);        % padded so edge peaks need no clipping
        keepPk  = false(nPkRaw, 1);
        for k = ordAmp(:)'
            tb = t0all(k) + halfDedup;               % index into the padded mask
            if blocked(tb), continue; end            % a bigger peak already owns this time
            keepPk(k) = true;
            blocked(tb-halfDedup : tb+halfDedup) = true;
        end
        r0all = r0all(keepPk);  t0all = t0all(keepPk);
    end
    nDedupAll_PK(end+1) = numel(r0all);                        %#ok<AGROW>

    if numel(r0all) > maxPkPerFile_PK          % test-run cap; roiPkRate still uses ALL peaks
        kp = sort(randperm(numel(r0all), maxPkPerFile_PK));
        r0all = r0all(kp);  t0all = t0all(kp);
    end
    Tb = nan(numel(r0all), nP);
    for e = 1:numel(r0all)
        r0 = r0all(e);  t0 = t0all(e);

        if useFastBase_PK
            bStat = [baseRawPK(r0,t0) baseGlPK(t0)];
            bMw   = baseMwPK(:,t0);
            bSub  = baseSubPK(:,t0);
        else
            bStat = [];  bMw = [];  bSub = [];
        end
        S = event_peak_stats(RawSH, Gl, r0, t0, avgwin_ext, baseWin_PK, oscP, bStat);
        Q = event_local_psd(RawSH, r0, S.tc, psdP);
        [lenH, durH, nRs, nSe] = peak_set_measure(MwPK, BpSH, sigPK, Dmat, AdjPK, ...
            r0, t0, baseWin_PK, avgwin_ext, Pset, bMw);
        lenR = NaN;  durR = NaN;  nRsR = NaN;
        if alsoRaw_PK
            bRw = [];
            if ~isempty(baseMwRaw), bRw = baseMwRaw(:,t0); end
            [lenR, durR, nRsR] = peak_set_measure(MwRaw, RawSH, sigRaw, Dmat, AdjPK, ...
                r0, t0, baseWin_PK, avgwin_ext, Pset, bRw);
        end
        aSub = event_extent(SubSH, Dmat, r0, t0, avgwin_ext, Dmax_ext);
        [aucS, aucA] = event_auc(SubSH, Dmat, aROI, r0, t0, aucWin_PK, Dmax_ext, ...
            baseWin_PK, avgwin_ext, bSub);

        Tb(e,:) = [aSub lenH durH Q.bandRatio Q.bandFrac Q.fLoc ...
                   S.oscRatio S.cmCorr bvF(r0) roiRate(r0) dsom(r0) badT(t0) ...
                   aucS aucA nRs nSe double(dtDirtPK(r0,t0)) lenR durR nRsR];
    end
    PkMeta = [PkMeta; repmat(f,numel(r0all),1) r0all t0all];   %#ok<AGROW>
    PkTab  = [PkTab;  Tb];                                     %#ok<AGROW>
    nRawAll_PK(end+1) = nPkRaw;                                %#ok<AGROW>
    nDone = find(foi_PK(:)' == f, 1);
    fprintf('  f%-3d %6d peaks scored of %6d in peakvec (%.0f%% collapsed as co-active ROIs) (%7d total, %.0f s elapsed, ~%.0f s left)\n', ...
        f, numel(r0all), nPkRaw, 100*(1 - numel(r0all)/max(nPkRaw,1)), size(PkTab,1), ...
        toc, toc/max(nDone,1)*(numel(foi_PK)-max(nDone,1)));
end
nPk = size(PkTab,1);
fprintf('ALL PEAKS scored: %d peaks x %d parameters from %d neurons in %.0f s\n', ...
    nPk, nP, numel(foi_PK), toc);
if dedupWin_PK > 0
    fprintf(['  one event = one peak: %d of %d peakvec peaks (%.1f%%) were co-active ROIs ' ...
        'within %g ms of a bigger peak and are represented by it\n'], ...
        sum(nRawAll_PK) - sum(nDedupAll_PK), sum(nRawAll_PK), ...
        100*(1 - sum(nDedupAll_PK)/max(sum(nRawAll_PK),1)), dedupWin_PK);
end
fprintf('  no peak was rejected on its values; nearSpk %.1f%% is a COLUMN, not a filter\n', ...
    100*mean(PkTab(:,iP('nearSpk'))>0));
fprintf(['  lenSet/durSet: ROI sets on %s, threshold %.2f, gap ROI judged at %.0f%% CI,\n' ...
    '     dendritic adjacency (%d branches, %d intra-branch links, %d branch points),\n' ...
    '     representative trace = %s ROI. Every other parameter uses the unfiltered trace.\n'], ...
    bpTxt_PK, Vthr_sh, 100*erf(ciZ_PK/sqrt(2)), adjInfo.nBranch, adjInfo.nIntra, ...
    adjInfo.nJunction, repTrace_PK);
% A fragmented ROI graph splits every event into one set per piece regardless of the
% amplitudes, so it has to be visible. adjMax_PK is what usually causes it.
if adjInfo.nCompRaw > 1
    fprintf(['     NOTE: the branch labels + adjMax_PK %s left the graph in %d pieces; ' ...
        '%d link(s) were\n     added to connect it (longest %.0f um). Raise adjMax_PK ' ...
        'to make the junction\n     search find these itself.\n'], mat2str(adjMax_PK), ...
        adjInfo.nCompRaw, adjInfo.nRepair, max(adjInfo.repairs(:,3)));
end
fprintf('  blood noise (%g-%g Hz) elevated at the peak: bandRatio median %.2f, %.1f%% above 1.5\n', ...
    fVess_PK(1), fVess_PK(2), median(PkTab(:,iP('bandRatio')),'omitnan'), ...
    100*mean(PkTab(:,iP('bandRatio')) > 1.5));
fprintf('  amplitude range %.3f - %.3f spike-height (median %.3f)\n', ...
    min(PkTab(:,1)), max(PkTab(:,1)), median(PkTab(:,1),'omitnan'));
fprintf(['  area under the peak over %g..%g ms, ROIs within %g um: aucSum median %.1f ' ...
    '(range %.1f - %.1f) spike-height*ms\n     mean AUC per unit ROI area (aucArea): ' ...
    'median %.3g (range %.3g - %.3g) spike-height*ms/um^2\n'], ...
    aucWin_PK(1), aucWin_PK(2), Dmax_ext, ...
    median(PkTab(:,iP('aucSum')),'omitnan'), min(PkTab(:,iP('aucSum'))), ...
    max(PkTab(:,iP('aucSum'))), median(PkTab(:,iP('aucArea')),'omitnan'), ...
    min(PkTab(:,iP('aucArea'))), max(PkTab(:,iP('aucArea'))));
% Three distinct outcomes at length 0, counted separately because they mean different
% things. This is the diagnostic that tells you whether lenSet is usable as a feature.
noSet  = PkTab(:,iP('nROIset')) == 0;            % the detecting ROI never reached Vthr
alone  = PkTab(:,iP('nROIset')) == 1;            % above threshold, but on its own
spread = PkTab(:,iP('nROIset')) >= 2;
fprintf(['  ROI sets: %d/%d peaks (%.1f%%) had the detecting ROI BELOW %.2f -> lenSet 0\n' ...
    '            %d (%.1f%%) above threshold but ALONE -> lenSet 0, limited by ROI spacing\n' ...
    '            %d (%.1f%%) spread over >=2 ROIs: length %.0f-%.0f um (median %.0f),\n' ...
    '               duration %.0f-%.0f ms (median %.0f), %.0f ROIs median\n'], ...
    sum(noSet), nPk, 100*mean(noSet), Vthr_sh, sum(alone), 100*mean(alone), ...
    sum(spread), 100*mean(spread), ...
    min(PkTab(spread,iP('lenSet'))), max(PkTab(spread,iP('lenSet'))), ...
    median(PkTab(spread,iP('lenSet')),'omitnan'), ...
    min(PkTab(spread,iP('durSet'))), max(PkTab(spread,iP('durSet'))), ...
    median(PkTab(spread,iP('durSet')),'omitnan'), ...
    median(PkTab(spread,iP('nROIset')),'omitnan'));
fprintf('  more than one set on the dendrite at the same instant: %.1f%% of peaks\n', ...
    100*mean(PkTab(:,iP('nSetEv')) > 1));
% The band-pass control. If the two agree the extent is not an artifact of filtering;
% where they disagree the RAW trace is the one carrying drift and 40-250 Hz noise
% across an absolute threshold, so a raw value LARGER than the filtered one is the
% expected direction of the bias, not evidence of a longer event.
if alsoRaw_PK
    lf = PkTab(:,iP('lenSet'));  lr = PkTab(:,iP('lenSetRaw'));
    gg = ~isnan(lf) & ~isnan(lr);
    both = gg & (lf > 0 | lr > 0);
    fprintf(['  CONTROL, same measurement on the unfiltered trace:\n' ...
        '     lenSet  filtered median %.0f um vs raw median %.0f um  (over the %d peaks\n' ...
        '        where either is non-zero); Spearman %.3f, identical on %.1f%% of them\n' ...
        '     raw LONGER on %.1f%%, raw SHORTER on %.1f%%\n' ...
        '     durSet  filtered median %.0f ms vs raw median %.0f ms\n'], ...
        median(lf(both)), median(lr(both)), sum(both), ...
        corr(lf(gg), lr(gg), 'Type','Spearman'), 100*mean(lf(both)==lr(both)), ...
        100*mean(lr(both)>lf(both)), 100*mean(lr(both)<lf(both)), ...
        median(PkTab(both,iP('durSet'))), median(PkTab(both,iP('durSetRaw'))));
end
dtd = PkTab(:,iP('dtDirt'));
fprintf(['  dirt: %.2f%% of ROI-frames masked in the last neuron (%d of %d ROIs affected).\n' ...
    '     Peaks by distance to the nearest dirt frame in their own ROI:\n' ...
    '        on dirt (0 ms) %d (%.1f%%) | <10 ms %d (%.1f%%) | <100 ms %d (%.1f%%) | ' ...
    'none within %g ms %d (%.1f%%)\n'], 100*dirtInfo.fracDirt, dirtInfo.nROIwithDirt, ...
    dirtInfo.nROI, sum(dtd==0), 100*mean(dtd==0), sum(dtd<10), 100*mean(dtd<10), ...
    sum(dtd<100), 100*mean(dtd<100), dtDirtMax_PK, sum(dtd>=dtDirtMax_PK), ...
    100*mean(dtd>=dtDirtMax_PK));
if mean(noSet | alone) > 0.5
    fprintf(['     NOTE: over half the peaks give lenSet = 0. They are not lost, but the ' ...
        'length x duration\n     plane cannot separate them from each other — only ' ...
        'the other parameters can.\n']);
end
fprintf('  parameters defined for every peak: %s\n', ...
    strjoin(pkNames(sum(isnan(PkTab),1) < 0.02*nPk), ', '));
fprintf('  parameters often NaN: %s\n', ...
    strjoin(pkNames(sum(isnan(PkTab),1) > 0.2*nPk), ', '));

% --- the hand-picked local candidates, scored identically ---
% These are the only peaks with an outside opinion attached, so they are marked in
% every histogram of the classifier and in the discrimination analysis.
exAMP = nan(nEx,1); exL50 = nan(nEx,1); exCOMP = nan(nEx,1);
exR0 = nan(nEx,1);  exProf = cell(nEx,1);  exPk = nan(nEx,nP);
for i = 1:nEx
    fe = exampleList(i,1);  te = exampleList(i,2);
    RawSH = interpolateNaN(PCresult{fe}.NormalizedTrace_dirt / PCresult{fe}.SpikeHeight);
    SubSH = interpolateNaN(PCresult{fe}.Subthreshold      / PCresult{fe}.SpikeHeight);
    Dmat = PCresult{fe}.interDendDist(PCresult{fe}.sortDist, PCresult{fe}.sortDist);
    Dmat = max(cat(3,Dmat,Dmat'),[],3);
    Gl   = mean(SubSH, 1, 'omitnan');
    dsom = abs(PCresult{fe}.dendaxis(:));
    r0 = candROIman(i);                        % hand-identified trigger ROI
    if isnan(r0), r0 = trigger_roi(SubSH, te, avgwin_ext); end
    exR0(i) = r0;
    [exAMP(i), exL50(i), exProf{i}] = event_extent(SubSH, Dmat, r0, te, avgwin_ext, Dmax_ext);
    l = PCresult{fe}.labelvec(r0);
    if ~isnan(l), exCOMP(i) = 1 + (l>=3); end
    S = event_peak_stats(RawSH, Gl, r0, te, avgwin_ext, baseWin_PK, oscP);
    Q = event_local_psd(RawSH, r0, S.tc, psdP);
    BpSH = bandpass_rows(RawSH, bpLenDur_PK, 1000, bpOrd_PK);
    MwE  = movmean(BpSH, 2*avgwin_ext+1, 2, 'omitnan');
    sigE = roi_noise(MwE, [], baseWin_PK, avgwin_ext, 'global');
    blE  = PCresult{fe}.BranchLabel(:);  blE = blE(PCresult{fe}.sortDist);
    AdjE = dend_adjacency(blE, dsom, Dmat, adjMax_PK);
    PsetE = struct('Vthr',Vthr_sh, 'ciZ',ciZ_PK, 'durMax',durMax_PK, ...
        'avgwin',avgwin_ext, 'repMode',repTrace_PK);
    [lenH, durH, nRsE, nSeE] = peak_set_measure(MwE, BpSH, sigE, Dmat, AdjE, r0, te, ...
        baseWin_PK, avgwin_ext, PsetE);
    lenRE = NaN;  durRE = NaN;  nRsRE = NaN;
    if alsoRaw_PK
        MwRawE  = movmean(RawSH, 2*avgwin_ext+1, 2, 'omitnan');
        sigRawE = roi_noise(MwRawE, [], baseWin_PK, avgwin_ext, 'global');
        [lenRE, durRE, nRsRE] = peak_set_measure(MwRawE, RawSH, sigRawE, Dmat, AdjE, ...
            r0, te, baseWin_PK, avgwin_ext, PsetE);
    end
    % the ROI-level and timing parameters too, so a multi-parameter classifier can
    % actually score the candidates (they used to be left NaN here, which silently
    % excluded all five from any model using bvFrac / roiPkRate / nearSpk)
    nTe = size(RawSH,2);
    ftpE = abs(PCresult{fe}.Ftprnts(:,:,PCresult{fe}.sortDist));
    bvFe = NaN;
    if isfield(PCresult{fe},'bvMask') && ~isempty(PCresult{fe}.bvMask)
        bvmE = imdilate(max(PCresult{fe}.bvMask,[],3) > 0, strel('disk',bvDil_PK));
        bvFe = sum(sum(ftpE(:,:,r0) .* double(bvmE))) / max(sum(sum(ftpE(:,:,r0))), eps);
    end
    rateE = sum(peakvec{fe}(r0,:) > 0) / (nTe/1000);
    badE = false(1,nTe);
    badE(PCresult{fe}.SpikeMat(1,:) > 0) = true;
    badE(PCresult{fe}.CStrace > 0)       = true;
    badE = imdilate(badE, ones(1, 2*spkWin_PK+1));
    aROIe = roi_areas(ftpE, PCresult{fe}.pixelsize, ftpThr_PK);
    [aucSe, aucAe] = event_auc(SubSH, Dmat, aROIe, r0, te, aucWin_PK, Dmax_ext, ...
        baseWin_PK, avgwin_ext);
    motRejE = [];
    if isfield(PCresult{fe},'motionReject'), motRejE = PCresult{fe}.motionReject; end
    dtDirtE = dirt_distance(PCresult{fe}.NormalizedTrace_dirt, motRejE, dtDirtMax_PK);
    exPk(i,:) = [exAMP(i) lenH durH Q.bandRatio Q.bandFrac Q.fLoc ...
                 S.oscRatio S.cmCorr bvFe rateE dsom(r0) badE(te) aucSe aucAe ...
                 nRsE nSeE double(dtDirtE(r0,te)) lenRE durRE nRsRE];
    fprintf(['Candidate %d (f%d ROI%d frame%d): amp %.3f  lenSet %.0f um (%d ROIs, %d set) ' ...
        'durSet %.0f ms  bandRatio %.2f  bandFrac %.2f  oscRatio %.2f  fLoc %.0f Hz  ' ...
        'aucSum %.1f  aucArea %.3g\n'], i, fe, r0, te, ...
        exAMP(i), lenH, nRsE, nSeE, durH, Q.bandRatio, Q.bandFrac, S.oscRatio, Q.fLoc, ...
        aucSe, aucAe);
end
exAMP_dff = arrayfun(@(i) exAMP(i) * PCresult{exampleList(i,1)}.SpikeHeight, 1:nEx);

% --- percentile of each candidate within the full peak population -------------
% If a parameter really separates local events from artifacts, the five candidates
% should sit together and away from the middle of the distribution.
fprintf('\nWhere the hand-picked candidates sit in the all-peak distribution (percentile):\n');
exPct = nan(nEx,nP);
for p = 1:nP
    col = PkTab(~isnan(PkTab(:,p)), p);
    if isempty(col), continue; end
    for i = 1:nEx
        if ~isnan(exPk(i,p)), exPct(i,p) = 100*mean(col <= exPk(i,p)); end
    end
end
[~, ordSep] = sort(abs(mean(exPct,1,'omitnan') - 50), 'descend');
for p = ordSep(1:min(12,nP))
    fprintf('  %-10s  %s  (median pct %.0f)\n', pkNames{p}, ...
        num2str(exPct(:,p)','%5.0f'), median(exPct(:,p),'omitnan'));
end
fprintf('  ^ suggestive only — 5 candidates is not a labelled set. That is the next cell.\n');

% --- backward-compatibility shim ---------------------------------------------
% The cells further down (amplitude-matched controls, figure 142, figure 145) still
% expect the old EventTab = [f r0 t0 amp L50 comp]. It is now just a VIEW of PkTab,
% built with the old 0.75 x mean-candidate-amplitude cut so those cells behave as
% before. Nothing above this line is thresholded, and PkTab is untouched.
ampThresh = 0.75 * mean(exAMP);
shim = find(PkTab(:,iP('amp')) >= ampThresh);
shimL50 = nan(numel(shim),1);  shimComp = nan(numel(shim),1);
for f = unique(PkMeta(shim,1))'
    SubSH = interpolateNaN(PCresult{f}.Subthreshold / PCresult{f}.SpikeHeight);
    Dmat = PCresult{f}.interDendDist(PCresult{f}.sortDist, PCresult{f}.sortDist);
    Dmat = max(cat(3,Dmat,Dmat'),[],3);
    lab = PCresult{f}.labelvec(:);
    for k = find(PkMeta(shim,1)==f)'
        [~, shimL50(k)] = event_extent(SubSH, Dmat, PkMeta(shim(k),2), PkMeta(shim(k),3), ...
            avgwin_ext, Dmax_ext);
        l = lab(PkMeta(shim(k),2));
        if ~isnan(l), shimComp(k) = 1 + (l>=3); end
    end
end
EventTab = [PkMeta(shim,1:3) PkTab(shim,iP('amp')) shimL50 shimComp];
fprintf(['\nEventTab shim for the downstream cells: %d/%d peaks above the old ' ...
    'ampThresh %.3f (PkTab itself is not filtered)\n'], numel(shim), nPk, ampThresh);
if isTestRun_PK
    warning(['TEST RUN: EventTab holds only %d of the %d neurons, so every cell below ' ...
        'the classifier is incomplete.'], numel(foi_PK), numel(foi));
end


%% ===== CLASSIFY a random sample of peaks by eye =====
% One peak at a time: its kymograph (all ROIs, sorted by distance from the soma, with
% the ROI that peakvec fired on marked), the trigger ROI's raw trace, and a histogram
% of each parameter over ALL peaks with this peak's value drawn on it. Keys:
%     r = real     a = artifact     u = unsure/skip
%     <- / ->      previous / next            (classifying auto-advances)
%     s = save     q = save and close
% Labels are written to labFile_PK after every keypress, so the session is resumable
% and can be stopped at any time. Re-running this cell reloads the same sample (the
% seed is stored in the file) and jumps to the first unlabelled peak.
%
% ALREADY-DONE WORK IS NEVER REPEATED. The labels live in labFile_PK and are keyed by
% peak identity (file, ROI, frame), so re-running this cell reads them back rather than
% asking again. If nothing is left to label, the cell loads the file, prints what is in
% it and RETURNS WITHOUT OPENING THE VIEWER — the kymographs are not even cut, which is
% the slow part. The viewer opens only when there are unlabelled peaks (because you
% raised nSample_PK, or scored more neurons) or when you ask for it with
% forceRelabel_PK. Nothing downstream needs the viewer: the ranking and classifier
% cells read labFile_PK themselves.
nSample_PK   = 300;        % how many peaks you want labelled in total
rngSeed_PK   = 1;          % change for a fresh independent sample
sampleMode_PK = 'random';  % 'random' = uniform over all peaks (what you asked for)
                           % 'stratified' = equal numbers per amp decile, which labels
                           % far more large peaks per hour of your time
topUpSample_PK = true;     % after importing old labels, draw fresh peaks up to nSample_PK
forceRelabel_PK = false;   % true = open the viewer even when every peak is labelled,
                           % e.g. to revisit and correct calls you already made
kymoHalf_PK  = 500;        % +/- ms shown in the kymograph and the trace
% every parameter gets a histogram; nothing is summarised away
histParams_PK = pkNames;      % all of them fit in a 4 x 4 grid — nothing hidden

% --- labels, matched by peak identity (see the note next to labFile_PK) ------------
[sampleIdx, pkLabel, srcLab] = load_peak_labels({labFile_PK, importFrom_PK}, PkMeta);
if isempty(sampleIdx)
    fprintf('No previous labels apply to these %d peaks — starting fresh.\n', nPk);
else
    fprintf('Carried over %d labels (%d real, %d artifact, %d unsure) from %s\n', ...
        numel(sampleIdx), sum(pkLabel==1), sum(pkLabel==0), sum(pkLabel==-1), srcLab);
end

% --- top the sample up to nSample_PK with fresh, never-seen peaks ---
nNeed = nSample_PK - numel(sampleIdx);
if topUpSample_PK && nNeed > 0
    rng(rngSeed_PK);
    pool = setdiff((1:nPk)', sampleIdx);
    switch lower(sampleMode_PK)
        case 'stratified'
            ed = prctile(PkTab(pool,iP('amp')), 0:10:100);  ed(end) = inf;
            per = ceil(nNeed/10);  add = [];
            for d = 1:10
                inD = pool(PkTab(pool,iP('amp')) >= ed(d) & PkTab(pool,iP('amp')) < ed(d+1));
                add = [add; inD(randperm(numel(inD), min(per,numel(inD))))]; %#ok<AGROW>
            end
            add = add(randperm(numel(add)));
        otherwise
            add = pool(randperm(numel(pool), min(nNeed, numel(pool))));
    end
    sampleIdx = [sampleIdx; add];              % already-labelled first, new ones after
    pkLabel   = [pkLabel;   nan(numel(add),1)];
    fprintf('Topped up with %d fresh %s peaks -> %d to label (%d still blank)\n', ...
        numel(add), sampleMode_PK, numel(sampleIdx), sum(isnan(pkLabel)));
end
labMeta = PkMeta(sampleIdx,:);                 % the identity key that makes this portable
nPkAll = nPk;
save(labFile_PK, 'sampleIdx','pkLabel','labMeta','rngSeed_PK','sampleMode_PK','nPkAll');

% --- is there anything left to do? --------------------------------------------------
% Everything above is cheap bookkeeping. The kymograph precompute and the viewer below
% are the expensive parts, so they are skipped entirely when every sampled peak already
% carries a label. Set forceRelabel_PK = true to open the viewer anyway.
nBlank_PK = sum(isnan(pkLabel));
openViewer_PK = nBlank_PK > 0 || forceRelabel_PK;
if ~openViewer_PK
    fprintf(['Labels are COMPLETE — nothing to classify, viewer not opened.\n' ...
        '  %d labelled peaks (%d real, %d artifact, %d unsure) loaded from %s\n' ...
        '  Revise them with forceRelabel_PK = true; raise nSample_PK to add fresh ' ...
        'peaks to label.\n'], numel(pkLabel), sum(pkLabel==1), sum(pkLabel==0), ...
        sum(pkLabel==-1), labFile_PK);
else
    fprintf('%d of %d sampled peaks unlabelled — cutting kymographs and opening the viewer.\n', ...
        nBlank_PK, numel(pkLabel));

    % --- precompute one kymograph slice per sampled peak (fast, robust viewer) ---
    % Rows are ALL ROIs of that neuron ordered by distance from the soma (basal first),
    % so the trigger ROI's position in the row order is itself informative.
    tAxis_PK = -kymoHalf_PK:kymoHalf_PK;
    nRmax = max(arrayfun(@(f) size(PCresult{f}.NormalizedTrace_dirt,1), unique(PkMeta(sampleIdx,1))));
    KymoPK = nan(nRmax, numel(tAxis_PK), numel(sampleIdx), 'single');
    TrigPK = nan(numel(sampleIdx), numel(tAxis_PK));      % trigger ROI raw trace
    rowTrigPK = nan(numel(sampleIdx),1);                  % its row in the kymograph
    rowSomaPK = nan(numel(sampleIdx),1);                  % where the soma sits in that order
    dendPK    = nan(nRmax, numel(sampleIdx));             % signed distance per row
    for f = unique(PkMeta(sampleIdx,1))'
        RawSH = interpolateNaN(PCresult{f}.NormalizedTrace_dirt / PCresult{f}.SpikeHeight);
        [dsorted, dord] = sort(PCresult{f}.dendaxis);     % basal -> soma -> apical
        [nR, nT] = size(RawSH);
        for j = find(PkMeta(sampleIdx,1)==f)'
            r0 = PkMeta(sampleIdx(j),2);  t0 = PkMeta(sampleIdx(j),3);
            W = t0 + tAxis_PK;  inb = W >= 1 & W <= nT;
            % the local baseline the measurements use, not a whole-recording median,
            % so a colour here is on the same scale as the numbers beside it
            base = median(RawSH(:, baseline_idx(t0, nT, baseWin_PK, avgwin_ext)), ...
                2, 'omitnan');
            KymoPK(1:nR, inb, j) = RawSH(dord, W(inb)) - base(dord);
            TrigPK(j, inb)  = RawSH(r0, W(inb)) - base(r0);
            rowTrigPK(j)    = find(dord == r0, 1);
            [~, rowSomaPK(j)] = min(abs(dsorted));
            dendPK(1:nR, j) = dsorted;
        end
    end
    fprintf('Kymographs ready: %d peaks, %d ROI rows, +/-%d ms (%.0f MB)\n', ...
        numel(sampleIdx), nRmax, kymoHalf_PK, numel(KymoPK)*4/2^20);

    cfgPK = struct('tAxis',tAxis_PK, 'histParams',{histParams_PK}, 'labFile',labFile_PK, ...
        'pkNames',{pkNames}, 'nPkAll',nPk, 'labMeta',labMeta);
    peak_classifier(PkTab, PkMeta, sampleIdx, pkLabel, KymoPK, TrigPK, rowTrigPK, ...
        rowSomaPK, dendPK, exPk, cfgPK);
end


%% ===== WHICH parameter separates the peaks I called real? =====
% Reads the labels back and ranks every parameter by how well it separates them, with
% no cutoff chosen in advance: AUC is the probability that a randomly chosen "real"
% peak scores above a randomly chosen "artifact" one, so 0.5 is useless and 0/1 is
% perfect separation (AUC < 0.5 just means the parameter runs the other way).
% Only then is a cutoff read off, at the point that maximises Youden's J.
minLab_PK = 20;        % refuse to rank parameters on fewer labels than this per class
nTopShow_PK = 8;       % parameters detailed in the figure
[sampleIdx, pkLabel, srcLab] = load_peak_labels({labFile_PK, importFrom_PK}, PkMeta);
if isempty(sampleIdx), error('No labels apply to these peaks — run the classifier cell first.'); end
fprintf('Labels from %s\n', srcLab);
isReal = pkLabel == 1;  isArt = pkLabel == 0;
fprintf('Labels: %d real, %d artifact, %d unsure, %d unlabelled (of %d sampled)\n', ...
    sum(isReal), sum(isArt), sum(pkLabel==-1), sum(isnan(pkLabel)), numel(pkLabel));
if sum(isReal) < minLab_PK || sum(isArt) < minLab_PK
    error(['Need at least %d of each class (have %d real, %d artifact) — ' ...
        'label more peaks.'], minLab_PK, sum(isReal), sum(isArt));
end

X = PkTab(sampleIdx,:);
aucP = nan(nP,1); nUse = nan(nP,1); dirP = nan(nP,1);
cutP = nan(nP,1); accP = nan(nP,1);
for p = 1:nP
    [aucP(p), nUse(p), cutP(p), accP(p)] = param_auc(X(:,p), isReal, isArt);
    dirP(p) = sign(aucP(p) - 0.5);
end
sepP = abs(aucP - 0.5) + 0.5;                       % direction-free separation
[~, ord] = sort(sepP, 'descend', 'MissingPlacement','last');

fprintf('\n%-11s %6s %6s %8s %9s %7s   %s\n', 'parameter','AUC','sep','cutoff','balanced','n','direction');
for p = ord(1:min(20,nP))
    if isnan(aucP(p)), continue; end
    fprintf('%-11s %6.3f %6.3f %8.3g %8.1f%% %7d   real is %s\n', pkNames{p}, aucP(p), ...
        sepP(p), cutP(p), 100*accP(p), nUse(p), repmat_str(dirP(p)>0,'HIGHER','LOWER'));
end
fprintf('\nBest single parameter: %s (AUC %.3f, balanced accuracy %.1f%% at %s %s %.3g)\n', ...
    pkNames{ord(1)}, aucP(ord(1)), 100*accP(ord(1)), pkNames{ord(1)}, ...
    repmat_str(dirP(ord(1))>0,'>=','<='), cutP(ord(1)));

% --- does a COMBINATION beat the best single parameter? ---
% Logistic regression on the top parameters, scored by 5-fold cross-validation so the
% comparison against the single-parameter AUC is honest. Parameters are used only
% where they are defined for nearly every labelled peak.
cvAUC = NaN;  useCombo = [];
if exist('fitglm','file')
    cand = ord(~isnan(aucP(ord)))';
    cand = cand(arrayfun(@(p) mean(isnan(X(isReal|isArt,p))) < 0.05, cand));
    useCombo = cand(1:min(5,numel(cand)));
    yy = isReal(isReal|isArt);  XX = X(isReal|isArt, useCombo);
    ok = all(~isnan(XX),2);  XX = XX(ok,:);  yy = yy(ok);
    rng(7);  cvp = mod(randperm(numel(yy)), 5) + 1;  sc = nan(numel(yy),1);
    for k = 1:5
        tr = cvp ~= k;  te = cvp == k;
        if numel(unique(yy(tr))) < 2, continue; end
        mdl = fitglm(XX(tr,:), yy(tr), 'Distribution','binomial');
        sc(te) = predict(mdl, XX(te,:));
    end
    cvAUC = param_auc(sc, yy==1, yy==0);
    fprintf(['5-fold CV logistic regression on {%s}: AUC %.3f  ' ...
        '(best single %.3f -> %s)\n'], strjoin(pkNames(useCombo),', '), cvAUC, ...
        sepP(ord(1)), repmat_str(cvAUC > sepP(ord(1))+0.02, ...
        'the combination helps', 'no better than the single parameter'));
end

% ===== figure 151: the ranking and what the top parameters look like =====
figure(151); clf; tiledlayout(3, max(4,nTopShow_PK/2), 'TileSpacing','compact');

nexttile([1 max(4,nTopShow_PK/2)]); hold all;      % A: the ranking
top = ord(~isnan(aucP(ord)));
barh(sepP(top(end:-1:1)) - 0.5, 'FaceColor',[0.4 0.4 0.4], 'EdgeColor','none');
set(gca, 'YTick',1:numel(top), 'YTickLabel',pkNames(top(end:-1:1)), 'FontSize',8);
xline(0.1,'k:'); xlabel('|AUC - 0.5|   (0 = useless, 0.5 = perfect)');
if ~isnan(cvAUC), xline(cvAUC-0.5, 'r--', 'Label','CV combination'); end
title('A  every parameter, ranked by how well it separates YOUR labels'); box off;

for j = 1:nTopShow_PK                              % B..: the top parameters
    p = ord(j);  if isnan(aucP(p)), continue; end
    nexttile([1 1]); hold all;
    vr = X(isReal,p);  va = X(isArt,p);
    ed = linspace(min(X(isReal|isArt,p)), max(X(isReal|isArt,p)), 24);
    if ed(1) == ed(end), ed = ed(1) + [-1 1]; end
    histogram(va, ed, 'Normalization','probability', 'FaceColor',[0.85 0.33 0.1], ...
        'EdgeColor','none', 'FaceAlpha',0.6);
    histogram(vr, ed, 'Normalization','probability', 'FaceColor',[0.2 0.2 0.2], ...
        'EdgeColor','none', 'FaceAlpha',0.5);
    xline(cutP(p), 'b--', 'linewidth',1.5);
    yl = ylim;
    for i = 1:nEx, xline(exPk(i,p), 'g-', 'linewidth',1); end
    ylim(yl);
    xlabel(pkNames{p}); ylabel('Fraction'); box off;
    title(sprintf('%s  AUC %.3f', char('A'+j), aucP(p)), 'FontWeight','normal');
end
set_font('Arial'); set_fontsize(10); set_figsize(560,300);

% ===== figure 152: ROC of the top parameters, and the best pair =====
figure(152); clf; tiledlayout(1,2,'TileSpacing','compact');
nexttile([1 1]); hold all;
cmapR = gen_colormap(Plasma, nTopShow_PK);
for j = 1:min(nTopShow_PK,numel(top))
    p = top(j);  s = dirP(p)*X(:,p);                % oriented so "real" is always high
    g = (isReal|isArt) & ~isnan(s);
    thr = [inf; sort(unique(s(g)),'descend'); -inf];
    tpr = nan(numel(thr),1);  fpr = tpr;
    for q = 1:numel(thr)
        pred = s >= thr(q);
        tpr(q) = mean(pred(isReal & g));  fpr(q) = mean(pred(isArt & g));
    end
    plot(fpr, tpr, '-', 'color',cmapR(j,:), 'linewidth',1.5);
end
plot([0 1],[0 1],'k:');
xlabel('False positive rate (artifacts called real)'); ylabel('True positive rate');
legend(pkNames(top(1:min(nTopShow_PK,numel(top)))), 'box','off', 'location','southeast', 'FontSize',8);
title('ROC, top parameters'); box off; axis square;

nexttile([1 1]); hold all;                          % the best pair
p1 = ord(1);  p2 = ord(2);
scatter(PkTab(:,p1), PkTab(:,p2), 6, [0.85 0.85 0.85], 'filled', 'MarkerFaceAlpha',0.15);
scatter(X(isArt,p1),  X(isArt,p2),  26, [0.85 0.33 0.1], 'filled', 'MarkerFaceAlpha',0.7);
scatter(X(isReal,p1), X(isReal,p2), 26, [0.2 0.2 0.2],  'filled', 'MarkerFaceAlpha',0.7);
plot(exPk(:,p1), exPk(:,p2), 'go', 'markersize',10, 'linewidth',2);
xline(cutP(p1),'b--'); yline(cutP(p2),'b--');
xlabel(pkNames{p1}); ylabel(pkNames{p2}); box off;
legend({'all peaks','artifact','real','hand-picked'}, 'box','off', 'location','best', 'FontSize',8);
title(sprintf('best pair: %s x %s', pkNames{p1}, pkNames{p2}));
set_font('Arial'); set_fontsize(11); set_figsize(400,180);

% --- apply the winning rule to every peak (reported, still not enforced) ---
score1 = dirP(ord(1)) * PkTab(:,ord(1));
passes = score1 >= dirP(ord(1))*cutP(ord(1));
fprintf(['\nIf %s %s %.3g were applied to all %d peaks it would keep %d (%.1f%%). ' ...
    'Nothing is filtered by this cell — set the cutoff where you want it.\n'], ...
    pkNames{ord(1)}, repmat_str(dirP(ord(1))>0,'>=','<='), cutP(ord(1)), nPk, ...
    sum(passes), 100*mean(passes));
save(fullfile(save_figto,['PeakParamRanking_' runTag_PK '.mat']), 'pkNames','aucP','sepP','cutP', ...
    'accP','dirP','ord','cvAUC','useCombo','sampleIdx','pkLabel');

%% ===== Multi-parameter classifier -> the length x duration plane =====
% Trains one logistic model on the hand labels using SEVERAL parameters at once, then
% applies it to every peak and draws the result on the length x duration plane — the
% same plane the retired L&D section used, except the colour now comes from a
% classifier fitted to your labels instead of from cutoffs chosen by hand.
%
% ONE DELIBERATE CHOICE: by default lenSet and durSet are NOT given to the
% classifier (excludeAxes_CL). If they were, the decision drawn on the length x
% duration plane would just be the model's own boundary in that plane and the figure
% could not tell you anything — it would be circular by construction. Trained without
% them, the plane becomes an INDEPENDENT check: do peaks the other parameters call
% real actually occupy a distinct region of length x duration? Set excludeAxes_CL =
% false to get the best possible classifier instead, and both AUCs are printed so the
% cost of leaving them out is visible.
%
% cmCorr is included by default but flagged: it was the top single parameter, yet its
% sign was the OPPOSITE of the physical prediction (real peaks correlate MORE with the
% global mean, i.e. "real" = spatially widespread) and its own cutoff rejects one of the
% five hand-picked LOCAL candidates. Since this paper is about local events, a model
% leaning on cmCorr may be learning "widespread", not "genuine". Drop it via
% dropParams_CL and compare — the candidate probabilities printed below are the test.
% aucSum / aucArea are dropped for the same reason as lenSet/durSet: they are
% amplitude x duration x spread rolled into one number, so a model trained on them
% would be judging events by their size, and any later "how big is the AUC of the real
% peaks" readout would be circular. They are measurements to be READ OUT of the
% classified peaks, not predictors. Take them out of dropParams_CL to test them.
excludeAxes_CL = true;                 % keep lenSet/durSet out of the predictors
dropParams_CL = {'nearSpk','lenSet','durSet','aucSum','aucArea', ...
                 'lenSetRaw','durSetRaw','nROIsetRaw'};   % e.g. add 'cmCorr' to test without it
nFold_CL   = 5;                        % cross-validation folds
rngSeed_CL = 11;
nShowScatter_CL = 20000;               % peaks drawn in the scatter panel
nBin_CL    = 36;                       % bins per axis for the binned-probability map
minBin_CL  = 15;                       % peaks needed in a bin before it is coloured
[sampleIdx, pkLabel, srcLab] = load_peak_labels({labFile_PK, importFrom_PK}, PkMeta);
if isempty(sampleIdx), error('No labels apply to these peaks — run the classifier cell first.'); end
fprintf('Labels from %s\n', srcLab);
labSel = pkLabel == 1 | pkLabel == 0;
y = pkLabel(labSel) == 1;
fprintf('Classifier on %d labelled peaks (%d real, %d artifact)\n', sum(labSel), sum(y), sum(~y));
if sum(y) < 15 || sum(~y) < 15
    error('Need >=15 of each class (have %d real, %d artifact).', sum(y), sum(~y));
end

% --- predictors ---
drop = dropParams_CL;
if excludeAxes_CL, drop = [drop, {'lenSet','durSet','lenSetRaw','durSetRaw'}]; end
pidx = find(~ismember(pkNames, drop));
Xlab = PkTab(sampleIdx(labSel), pidx);
badP = mean(isnan(Xlab),1) >= 0.05 | std(Xlab,0,1,'omitnan') == 0;   % unusable predictors
if any(badP)
    fprintf('  dropped as unusable (>5%% NaN or zero variance among the labels): %s\n', ...
        strjoin(pkNames(pidx(badP)), ', '));
end
pidx = pidx(~badP);
Xlab = PkTab(sampleIdx(labSel), pidx);
fprintf('  predictors (%d): %s\n', numel(pidx), strjoin(pkNames(pidx), ', '));

% z-score on the labelled set, then apply those SAME numbers everywhere, so the
% coefficients are comparable and the all-peak scores are on the training scale
mu = mean(Xlab,1,'omitnan');  sd = std(Xlab,0,1,'omitnan');  sd(sd==0) = 1;
Z = (Xlab - mu) ./ sd;
ok = all(~isnan(Z),2);
Zok = Z(ok,:);  yok = y(ok);
fprintf('  %d/%d labelled peaks have all predictors defined\n', sum(ok), numel(y));

% --- cross-validated performance (out-of-fold, so it is not self-graded) ---
rng(rngSeed_CL);
cvg = mod(randperm(numel(yok)), nFold_CL) + 1;
oof = nan(numel(yok),1);
ws = warning('off','all');
for k = 1:nFold_CL
    tr = cvg ~= k;  te = cvg == k;
    if numel(unique(yok(tr))) < 2, continue; end
    mk = fitglm(Zok(tr,:), yok(tr), 'Distribution','binomial');
    oof(te) = predict(mk, Zok(te,:));
end
mdlCL = fitglm(Zok, yok, 'Distribution','binomial');     % final model, all labels
warning(ws);
[aucCL, ~, thrCL, accCL] = param_auc(oof, yok==1, yok==0);
aucIn = param_auc(predict(mdlCL, Zok), yok==1, yok==0);
% best single parameter, for reference (same labelled peaks)
aucSingle = nan(nP,1);
for p = 1:nP, aucSingle(p) = param_auc(PkTab(sampleIdx(labSel),p), y==1, y==0); end
sepSingle = abs(aucSingle-0.5)+0.5;
[bestSep, bestP] = max(sepSingle);
seCL = sqrt(aucCL*(1-aucCL)/min(sum(yok),sum(~yok)));    % rough SE, for a sanity check
fprintf(['  cross-validated AUC %.3f (+/-~%.3f)   in-sample %.3f   best single (%s) %.3f\n' ...
    '  -> %s\n'], aucCL, seCL, aucIn, pkNames{bestP}, bestSep, ...
    repmat_str(aucCL > bestSep + 2*seCL, 'the combination genuinely beats any one parameter', ...
    'NOT clearly better than the best single parameter — the extra parameters may be redundant'));
fprintf('  operating point: p >= %.3f (Youden on the out-of-fold scores), balanced accuracy %.1f%%\n', ...
    thrCL, 100*accCL);

% --- standardized coefficients: which parameters the model actually leans on ---
cf = mdlCL.Coefficients;
bStd = cf.Estimate(2:end);  pVal = cf.pValue(2:end);
[~, ordC] = sort(abs(bStd), 'descend');
fprintf('\n  %-11s %9s %9s   %s\n', 'predictor','coef(z)','p','pushes towards');
for j = ordC'
    fprintf('  %-11s %+9.2f %9.1g   %s\n', pkNames{pidx(j)}, bStd(j), pVal(j), ...
        repmat_str(bStd(j)>0, 'REAL when high', 'REAL when low'));
end

% --- apply to every peak, and to the hand-picked candidates ---
Zall = (PkTab(:,pidx) - mu) ./ sd;
pReal = nan(nPk,1);
okAll = all(~isnan(Zall),2);
ws = warning('off','all');
pReal(okAll) = predict(mdlCL, Zall(okAll,:));
Zex = (exPk(:,pidx) - mu) ./ sd;
pEx = nan(nEx,1);  okEx = all(~isnan(Zex),2);
if any(okEx), pEx(okEx) = predict(mdlCL, Zex(okEx,:)); end
warning(ws);
fprintf('\n  scored %d/%d peaks (%d had a missing predictor); %.1f%% above p=%.3f\n', ...
    sum(okAll), nPk, sum(~okAll), 100*mean(pReal(okAll)>=thrCL), thrCL);
fprintf('  hand-picked LOCAL candidates: p = %s\n', num2str(pEx','%.3f  '));
fprintf('  -> %d/%d candidates classified REAL. %s\n', sum(pEx>=thrCL), nEx, ...
    repmat_str(sum(pEx>=thrCL)==nEx, 'all of them survive, which is what we want.', ...
    'SOME ARE REJECTED — the model is learning something these events do not have.'));

% ===== figure 153: the length x duration plane, coloured by the classifier =====
iL = iP('lenSet');  iD = iP('durSet');
L = PkTab(:,iL);  D = PkTab(:,iD);
edL = linspace(0, prctile(L,99.5), nBin_CL+1);
edD = linspace(0, prctile(D,99.5), nBin_CL+1);
figure(153); clf; tiledlayout(2,2,'TileSpacing','compact');

% A: mean classifier probability in each length x duration bin — the whole plane at once
nexttile([1 1]); hold all;
gL = discretize(L, edL);  gD = discretize(D, edD);
g = ~isnan(gL) & ~isnan(gD) & ~isnan(pReal);
Mp = accumarray([gD(g) gL(g)], pReal(g), [nBin_CL nBin_CL], @mean, NaN);
Nb = accumarray([gD(g) gL(g)], 1,        [nBin_CL nBin_CL], @sum,  0);
Mp(Nb < minBin_CL) = NaN;
imagesc(edL(1:end-1), edD(1:end-1), Mp, 'AlphaData', ~isnan(Mp));
set(gca,'YDir','normal','Color',[0.94 0.94 0.94]);  clim([0 1]);  colormap(gca, turbo);
contour(edL(1:end-1), edD(1:end-1), Mp, thrCL*[1 1], 'k-', 'linewidth',1.5);
plot(exPk(:,iL), exPk(:,iD), 'wo', 'markersize',10, 'linewidth',2);
plot(exPk(:,iL), exPk(:,iD), 'ko', 'markersize',10, 'linewidth',1);
cb = colorbar;  cb.Label.String = 'p(real)';
xlabel(sprintf('Length above %.2f spike height (µm)', Vthr_sh));
ylabel(sprintf('Duration above %.2f spike height (ms)', Vthr_sh));
title(sprintf('A  mean p(real) per bin (>=%d peaks), black = p %.2f', minBin_CL, thrCL));
axis tight; box off;

% B: the peaks themselves, coloured by p(real) — the scatter you asked for.
% This panel is static. The next cell draws the same scatter in its own figure with the
% dots CLICKABLE: click one and its kymograph appears beside it.
nexttile([1 1]); hold all;
sel = find(okAll & ~isnan(L) & ~isnan(D));
if numel(sel) > nShowScatter_CL
    rng(rngSeed_CL);  sel = sel(randperm(numel(sel), nShowScatter_CL));
end
[~, so] = sort(pReal(sel));  sel = sel(so);          % draw likely-real on top
scatter(L(sel), D(sel), 7, pReal(sel), 'filled', 'MarkerFaceAlpha',0.45);
clim([0 1]);  colormap(gca, turbo);
plot(exPk(:,iL), exPk(:,iD), 'ko', 'markersize',11, 'linewidth',2, 'MarkerFaceColor','w');
cb = colorbar;  cb.Label.String = 'p(real)';
xlabel(sprintf('Length above %.2f spike height (µm)', Vthr_sh));
ylabel(sprintf('Duration above %.2f spike height (ms)', Vthr_sh));
xlim([0 edL(end)]);  ylim([0 edD(end)]);
title(sprintf('B  %d of %d peaks (black rings = hand-picked local)', numel(sel), sum(okAll)));
box off;

% C: the labelled peaks in the same plane — does the classifier agree with your eye?
nexttile([1 1]); hold all;
scatter(L(sel), D(sel), 5, [0.85 0.85 0.85], 'filled', 'MarkerFaceAlpha',0.25);
si = sampleIdx(labSel);
scatter(L(si(~y)), D(si(~y)), 34, [0.85 0.33 0.1], 'filled', 'MarkerFaceAlpha',0.75);
scatter(L(si( y)), D(si( y)), 34, [0.15 0.15 0.15], 'filled', 'MarkerFaceAlpha',0.75);
contour(edL(1:end-1), edD(1:end-1), Mp, thrCL*[1 1], 'b-', 'linewidth',1.5);
plot(exPk(:,iL), exPk(:,iD), 'go', 'markersize',11, 'linewidth',2);
xlabel(sprintf('Length above %.2f spike height (µm)', Vthr_sh));
ylabel(sprintf('Duration above %.2f spike height (ms)', Vthr_sh));
xlim([0 edL(end)]);  ylim([0 edD(end)]);
legend({'all peaks','you: artifact','you: real','classifier boundary','hand-picked local'}, ...
    'box','off','location','northeast','FontSize',8);
title('C  your labels on the same plane'); box off;

% D: what the model leans on, and how well it does
nexttile([1 1]); hold all;
barh(bStd(ordC(end:-1:1)), 'FaceColor',[0.4 0.4 0.4], 'EdgeColor','none');
set(gca,'YTick',1:numel(ordC), 'YTickLabel',pkNames(pidx(ordC(end:-1:1))), 'FontSize',9);
plot([0 0], ylim, 'k-');
xlabel('standardized logistic coefficient   (+ = real when high)');
title(sprintf('D  CV AUC %.3f vs best single %.3f (%s)', aucCL, bestSep, pkNames{bestP}));
box off;
set_font('Arial'); set_fontsize(11); set_figsize(420,340);

% --- where the two disagree, so it can be looked at rather than trusted ---
pLab = pReal(si);
disagree = (pLab >= thrCL) ~= y;
fprintf('\n  classifier disagrees with your label on %d/%d labelled peaks (%.0f%%)\n', ...
    sum(disagree & ~isnan(pLab)), sum(~isnan(pLab)), 100*mean(disagree(~isnan(pLab))));
dr = find(disagree & ~isnan(pLab));
[~, dord] = sort(abs(pLab(dr)-thrCL), 'descend');     % most confident mistakes first
fprintf('  most confident disagreements (f, ROI, frame, p, your label):\n');
for j = dr(dord(1:min(8,numel(dord))))'
    fprintf('     f%-3d ROI%-3d frame%-8d p=%.3f  you said %s\n', PkMeta(si(j),1), ...
        PkMeta(si(j),2), PkMeta(si(j),3), pLab(j), repmat_str(y(j),'REAL','ARTIFACT'));
end
fprintf('  (sampleIdx positions %s — browse them in the classifier cell with the arrow keys)\n', ...
    num2str(find(ismember(sampleIdx, si(dr(dord(1:min(8,numel(dord))))))')));
save(fullfile(save_figto,['PeakClassifier_' runTag_PK '.mat']), 'mdlCL','pidx','pkNames','mu','sd', ...
    'thrCL','aucCL','pReal','pEx','excludeAxes_CL','dropParams_CL');

%% ===== ROI SETS: split an event into sets, then measure each one ====================
% The rule, in order:
%   1. a peak is detected in some ROI (peakvec, as before)
%   2. every ROI above Vthr_sh AT THAT ONE INSTANT is collected. Amplitude is the
%      +/-avgwin window mean of the 0.1-60 Hz BANDPASSED trace, baseline-subtracted,
%      read at the SAME t0 for every ROI. One common instant is essential: if each ROI
%      contributed its own maximum over a search window, signals happening at different
%      times would be swept into one set.
%   3. those ROIs form SETS, connected along the ACTUAL DENDRITE (dend_adjacency:
%      consecutive within a BranchLabel, plus one link per branch point), not by
%      distance ordering. When two sets are separated by EXACTLY ONE ROI, that ROI
%      decides: the separation is believed only when the ROI is confidently below
%      threshold, i.e. Vthr is above the TOP of its 95% CI (Vthr >= A + ciZ*sigma). A
%      noisy ROI whose CI reaches past 0.15 is no evidence of a gap, so the sets merge.
%   4. separated by two or more ROIs: always different sets, no test
%   5. each set's LENGTH is half the contour distance between its two farthest members.
%      HEIGHT and DURATION both come from ONE representative trace (repTrace_ST), so
%      they describe the same signal.
%
% TWO THINGS THAT ARE NOT ASSUMED HERE
%   * One physical event fires peakvec in several ROIs, which would enter this cell as
%     several detections of the same thing. Detections in the same neuron that overlap
%     IN TIME (within dedupWin_ST) and SHARE AT LEAST ONE ROI are therefore collapsed
%     into one EVENT, and only its representative detection is measured. eventID and
%     nDet in SetTab record the collapse, so nothing is counted twice.
%   * Sets that are spatially separate at one instant are NOT assumed to be independent
%     biological events. They are recorded as separate DOMAINS of the same event
%     (same eventID, different setNo). Whether two domains are one event or two is a
%     question this cell reports on rather than answers — nSet per event is in the
%     printout, and figure 164 shows examples so the cases can be judged by eye.
setVthr_ST  = Vthr_sh;  % 0.15 spike height, the established threshold
ciZ_ST      = 1.96;     % 95% CI. 1.645 = one-sided 95%, 2.58 = 99%
repTrace_ST = 'peak';   % the set's representative trace: 'peak' (largest member) or
                        % 'center' (the 1-centre). Height AND duration come from it.
dedupWin_ST = 50;       % ms: detections this close in time that share an ROI are the
                        % same event. 0 = treat every peakvec detection separately
nTest_ST    = 3000;     % peaks to test (random sample of PkTab; Inf = all)
adjMax_ST   = 60;       % µm cap on branch-point links only (see dend_adjacency). [] = none
noiseMode_ST = 'global';% sigma for the CI: 'global' (whole recording, better estimated)
                        % or 'local' (each event's own baseline window)
durMax_ST   = 100;      % +/- ms cap on the duration walk
rngSeed_ST  = 5;
if ~exist('PkTab','var'), error('Run the all-peaks scoring cell first.'); end

rng(rngSeed_ST);
selST = (1:nPk)';
if numel(selST) > nTest_ST, selST = sort(selST(randperm(nPk, nTest_ST))); end
stNames = {'eventID','nDet','f','r0','t0','setNo','nSet','hasTrigger','nROI', ...
           'lenHalf','lenSpan','height','duration','centerROI','repROI', ...
           'nGapTested','nGapMerged'};
iS = @(s) find(strcmp(stNames,s), 1);
SetTab = [];  evCount = 0;  adjTxt = '';
Pst = struct('Vthr',setVthr_ST, 'ciZ',ciZ_ST, 'durMax',durMax_ST, ...
    'avgwin',avgwin_ext, 'repMode',repTrace_ST, 'wantPath',true);   % the maps draw it
tic
for f = unique(PkMeta(selST,1))'
    RawSH = interpolateNaN(PCresult{f}.NormalizedTrace_dirt / PCresult{f}.SpikeHeight);
    BpSH  = bandpass_rows(RawSH, bpLenDur_PK, 1000, bpOrd_PK);   % 0.1-60 Hz
    Mw    = movmean(BpSH, 2*avgwin_ext+1, 2, 'omitnan');
    Dmat  = PCresult{f}.interDendDist(PCresult{f}.sortDist, PCresult{f}.sortDist);
    Dmat  = max(cat(3,Dmat,Dmat'),[],3);
    bl    = PCresult{f}.BranchLabel(:);  bl = bl(PCresult{f}.sortDist);
    dsom  = abs(PCresult{f}.dendaxis(:));
    [Adj, ainfo] = dend_adjacency(bl, dsom, Dmat, adjMax_ST);
    sigG  = roi_noise(Mw, [], baseWin_PK, avgwin_ext, 'global');
    rows  = find(PkMeta(selST,1) == f);

    % --- pass 1: sets for every detection, nothing collapsed yet ---
    nD = numel(rows);
    detS = cell(nD,1);  detInfo = cell(nD,1);
    detMem = false(nD, size(Mw,1));  detT = nan(nD,1);  detH = nan(nD,1);
    for q = 1:nD
        k = selST(rows(q));  r0 = PkMeta(k,2);  t0 = PkMeta(k,3);
        if strcmpi(noiseMode_ST,'global'), sg = sigG;
        else, sg = roi_noise(Mw, t0, baseWin_PK, avgwin_ext, 'local');
        end
        B = t0 + (-baseWin_PK:baseWin_PK);  B = B(B >= 1 & B <= size(Mw,2));
        B = setdiff(B, t0 + (-3*avgwin_ext:3*avgwin_ext));
        if isempty(B), B = 1:size(Mw,2); end
        base = median(Mw(:,B), 2, 'omitnan');
        A    = Mw(:,t0) - base;                 % SAME instant for every ROI
        [S, info] = event_roi_sets(A, sg, BpSH, base, Dmat, Adj, r0, t0, Pst);
        detS{q} = S;  detInfo{q} = info;  detT(q) = t0;
        if ~isempty(S)
            detMem(q, [S.roi]) = true;
            detH(q) = max([S.height]);
        end
    end

    % --- pass 2: collapse detections of the SAME event ---
    % Same neuron, within dedupWin_ST ms, and sharing at least one above-threshold ROI.
    % Transitive, so a chain of detections across a spreading event becomes one event.
    grpD = 1:nD;
    if dedupWin_ST > 0
        [~, ordT] = sort(detT);
        for ii = 1:nD
            a = ordT(ii);
            if ~any(detMem(a,:)), continue; end
            for jj = ii+1:nD
                b = ordT(jj);
                if detT(b) - detT(a) > dedupWin_ST, break; end
                if any(detMem(a,:) & detMem(b,:))
                    grpD(grpD == grpD(b)) = grpD(a);
                end
            end
        end
    end

    % --- pass 3: one representative detection per event ---
    for g = unique(grpD)
        mem = find(grpD == g);
        mem = mem(arrayfun(@(q) ~isempty(detS{q}), mem));
        if isempty(mem), continue; end
        [~, best] = max(detH(mem));             % the detection with the tallest set
        q = mem(best);
        S = detS{q};  info = detInfo{q};
        k = selST(rows(q));
        evCount = evCount + 1;
        for c = 1:numel(S)
            SetTab(end+1,:) = [evCount numel(mem) f PkMeta(k,2) PkMeta(k,3) ...
                c numel(S) S(c).hasTrigger S(c).nROI S(c).lenHalf S(c).lenSpan ...
                S(c).height S(c).duration S(c).center S(c).repROI ...
                size(info.gaps,1) sum(info.gaps(:,4))]; %#ok<SAGROW>
        end
    end
    adjTxt = sprintf('%d branches, %d intra-branch links, %d branch points', ...
        ainfo.nBranch, ainfo.nIntra, ainfo.nJunction);
    fprintf('  f%-3d %4d detections -> %4d events -> %4d sets   (%s)  [%.0f s]\n', ...
        f, nD, numel(unique(grpD(arrayfun(@(q) ~isempty(detS{q}), 1:nD)))), ...
        sum(SetTab(:,iS('f'))==f), adjTxt, toc);
end

% --- what came out ---
nDet = numel(selST);
if isempty(SetTab), error('No ROI reached %.2f in any sampled peak.', setVthr_ST); end
ev1  = SetTab(SetTab(:,iS('setNo'))==1, :);      % one row per event
nEvS = size(ev1,1);
fprintf('\nROI SETS from %d peakvec detections (Vthr %.2f on the %g-%g Hz trace, %.0f%% CI, %s sigma)\n', ...
    nDet, setVthr_ST, bpLenDur_PK(1), bpLenDur_PK(2), 100*erf(ciZ_ST/sqrt(2)), noiseMode_ST);
fprintf('  detections collapsed into %d EVENTS (same neuron, <=%g ms apart, sharing an ROI)\n', ...
    nEvS, dedupWin_ST);
fprintf('    %.1f detections per event on average — that many would have been double counted\n', ...
    mean(ev1(:,iS('nDet'))));
for q = 1:max(ev1(:,iS('nSet')))
    n = sum(ev1(:,iS('nSet')) == q);
    if n == 0, continue; end
    fprintf('  %d spatial set(s) per event: %5d events (%.1f%%)\n', q, n, 100*n/max(nEvS,1));
end
fprintf('  ^ separate sets at one instant are DOMAINS of one detection, not proven to be\n');
fprintf('    independent events — figure 164 shows examples to judge that by eye\n');
fprintf('  single-ROI gaps tested: %d, of which MERGED (gap not trustworthy): %d (%.1f%%)\n', ...
    sum(ev1(:,iS('nGapTested'))), sum(ev1(:,iS('nGapMerged'))), ...
    100*sum(ev1(:,iS('nGapMerged'))) / max(sum(ev1(:,iS('nGapTested'))),1));
trg = SetTab(:,iS('hasTrigger')) == 1;
fprintf('  sets containing the detecting ROI: %d of %d (%.1f%%)\n', ...
    sum(trg), size(SetTab,1), 100*mean(trg));
fprintf('  trigger sets : lenHalf median %.0f µm (span %.0f), height %.2f, dur %.0f ms, %.0f ROIs\n', ...
    median(SetTab(trg,iS('lenHalf'))), median(SetTab(trg,iS('lenSpan'))), ...
    median(SetTab(trg,iS('height'))), median(SetTab(trg,iS('duration'))), ...
    median(SetTab(trg,iS('nROI'))));
if any(~trg)
    fprintf('  other sets   : lenHalf median %.0f µm (span %.0f), height %.2f, dur %.0f ms, %.0f ROIs\n', ...
        median(SetTab(~trg,iS('lenHalf'))), median(SetTab(~trg,iS('lenSpan'))), ...
        median(SetTab(~trg,iS('height'))), median(SetTab(~trg,iS('duration'))), ...
        median(SetTab(~trg,iS('nROI'))));
end

% ===== figure 163: what the sets look like =====
figure(163); clf; tiledlayout(2,2,'TileSpacing','compact');
nexttile([1 1]); hold all;
histogram(SetTab(trg,iS('lenHalf')), 20, 'FaceColor',[0.2 0.2 0.2], 'EdgeColor','none');
if any(~trg)
    histogram(SetTab(~trg,iS('lenHalf')), 20, 'FaceColor',[0.85 0.33 0.1], ...
        'EdgeColor','none', 'FaceAlpha',0.6);
    legend({'set with the detecting ROI','other sets'}, 'box','off');
end
xlabel('lenHalf = span/2 (µm)'); ylabel('sets'); box off;
title('A  length of every set');

nexttile([1 1]); hold all;
scatter(SetTab(:,iS('lenHalf')), SetTab(:,iS('height')), 8, SetTab(:,iS('nROI')), 'filled', ...
    'MarkerFaceAlpha',0.4);
colormap(gca, turbo); cb = colorbar; cb.Label.String = 'ROIs in set';
xlabel('lenHalf (µm)'); ylabel('height (spike height)'); box off;
title(sprintf('B  length vs height  (both from the %s ROI)', repTrace_ST));

nexttile([1 1]); hold all;
scatter(SetTab(:,iS('lenHalf')), SetTab(:,iS('duration')), 8, [0.3 0.3 0.3], 'filled', ...
    'MarkerFaceAlpha',0.3);
xlabel('lenHalf (µm)'); ylabel('duration (ms)'); box off;
title('C  length vs duration');

nexttile([1 1]); hold all;
histogram(ev1(:,iS('nSet')), 0.5:1:(max(ev1(:,iS('nSet')))+0.5), 'FaceColor',[0.3 0.3 0.3], ...
    'EdgeColor','none');
xlabel('spatial sets per event'); ylabel('events'); box off;
title(sprintf('D  %.1f%% of single-ROI gaps merged; %.1f det/event', ...
    100*sum(ev1(:,iS('nGapMerged')))/max(sum(ev1(:,iS('nGapTested'))),1), ...
    mean(ev1(:,iS('nDet')))));
sgtitle(sprintf('ROI sets: Vthr %.2f, one instant, dendritic adjacency, gap ROI at %.0f%% CI', ...
    setVthr_ST, 100*erf(ciZ_ST/sqrt(2))), 'FontSize',11);
set_font('Arial'); set_fontsize(10); set_figsize(460,340);

% ===== figure 164: individual examples, one row each =====
% Chosen to cover the cases: one set; two sets kept apart by a trustworthy gap ROI;
% two sets MERGED because the gap ROI was too noisy to believe; three or more sets.
pickST = [];
i1 = find(ev1(:,iS('nSet'))==1 & ev1(:,iS('nROI'))>2, 2);
i2 = find(ev1(:,iS('nSet'))==2 & ev1(:,iS('nGapMerged'))==0, 2);
i3 = find(ev1(:,iS('nGapMerged'))>0, 2);
i4 = find(ev1(:,iS('nSet'))>=3, 2);
pickST = unique([i1; i2; i3; i4], 'stable');
% Figure 165 shows the SAME events on the morphology, three panels each:
%   i)   dendrite amplitude map — every ROI painted with its amplitude at the event
%        instant; ROIs below Vthr are grey, so "who is above threshold" is visible
%   ii)  the ROI sets — the same footprints coloured by which set they ended up in,
%        with the tested gap ROIs outlined
%   iii) the connected sets — the dendritic adjacency graph itself, nodes at the ROI
%        centroids and an edge for every neighbour link, so you can see WHY two blobs
%        are one set or two: the edges are the connectivity the rule walked on
figure(164); clf;
figure(165); clf;
if isempty(pickST)
    figure(164); text(0.5,0.5,'no events to show','HorizontalAlignment','center');
else
    figure(164); tl164 = tiledlayout(numel(pickST), 1, 'TileSpacing','compact');
    figure(165); tl165 = tiledlayout(numel(pickST), 3, 'TileSpacing','compact', ...
        'Padding','compact');
    for rr = pickST'
        f = ev1(rr,iS('f'));  r0 = ev1(rr,iS('r0'));  t0 = ev1(rr,iS('t0'));
        RawSH = interpolateNaN(PCresult{f}.NormalizedTrace_dirt / PCresult{f}.SpikeHeight);
        BpSH  = bandpass_rows(RawSH, bpLenDur_PK, 1000, bpOrd_PK);
        Mw    = movmean(BpSH, 2*avgwin_ext+1, 2, 'omitnan');
        Dmat  = PCresult{f}.interDendDist(PCresult{f}.sortDist, PCresult{f}.sortDist);
        Dmat  = max(cat(3,Dmat,Dmat'),[],3);
        bl    = PCresult{f}.BranchLabel(:);  bl = bl(PCresult{f}.sortDist);
        dsom  = abs(PCresult{f}.dendaxis(:));
        Adj   = dend_adjacency(bl, dsom, Dmat, adjMax_ST);
        if strcmpi(noiseMode_ST,'global')
            sg = roi_noise(Mw, [], baseWin_PK, avgwin_ext, 'global');
        else
            sg = roi_noise(Mw, t0, baseWin_PK, avgwin_ext, 'local');
        end
        B = t0 + (-baseWin_PK:baseWin_PK);  B = B(B>=1 & B<=size(Mw,2));
        B = setdiff(B, t0 + (-3*avgwin_ext:3*avgwin_ext));
        base = median(Mw(:,B), 2, 'omitnan');
        A = Mw(:,t0) - base;
        [S, info] = event_roi_sets(A, sg, BpSH, base, Dmat, Adj, r0, t0, Pst);
        ttl = sprintf('f%d ROI%d fr%d (%d det)', f, r0, t0, ev1(rr,iS('nDet')));
        ax = nexttile(tl164);
        draw_roi_sets(ax, A, sg, S, info, PCresult{f}.dendaxis(:), r0, setVthr_ST, ...
            ciZ_ST, ttl);
        % --- the same event on the morphology ---
        ftM  = abs(PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist)) > 0;
        cooR = get_coord(PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist));
        draw_roi_set_maps(nexttile(tl165), nexttile(tl165), nexttile(tl165), ...
            ftM, cooR, A, S, info, Adj, r0, setVthr_ST, ttl);
    end
    figure(164); xlabel('position along dendrite (µm from soma)');
end
figure(164); set_font('Arial'); set_fontsize(9);
set_figsize(430, 95*max(numel(pickST),1)+50);
figure(165); set_font('Arial'); set_fontsize(9);
set_figsize(620, 190*max(numel(pickST),1)+30);
fprintf('\nSetTab: %d sets from %d events (%d detections). Columns:\n  %s\n', ...
    size(SetTab,1), nEvS, nDet, strjoin(stNames, ', '));

%% ===== Length distribution of the peaks the classifier calls real =====
% THE HISTOGRAM. lenSet — the FULL length of the ROI set the detecting ROI belongs to
% — for every peak the classifier calls real, inside an amplitude band YOU choose.
%
% The length axis is a FAIR readout here only because excludeAxes_CL kept lenSet out of
% the classifier's predictors: the model decided "real" from the blood-noise, shape and
% ROI parameters, and the length distribution below is something it was never shown.
% Set excludeAxes_CL = false and this figure becomes circular — a warning is printed.
%
% AMPLITUDE — the knob for choosing which events are in the histogram is the 'amp' row
% of extraFilt_CL. Length is measured at an ABSOLUTE threshold (0.15 spike height), so
% it is coupled to amplitude by construction: the same event at twice the amplitude
% clears 0.15 over more ROIs. A distribution pooled over all amplitudes therefore
% partly reports the amplitude distribution. Restricting to a band compares like with
% like. The coupling does not disappear on its own: within any band the lengths are
% still spread, and two bands are two different measurements, so quote the band
% alongside any length you report. Figure 154 panel 2 draws the band so you can see
% what the choice selected.
%
% WHAT SITS AT LENGTH 0 — two different things, so the counts are printed separately:
%   (a) the detecting ROI never reaches 0.15 (nROIset = 0). There is no set to measure,
%       so this is not "length zero", it is "not measurable by this definition".
%   (b) the ROI clears 0.15 but its set is that ROI alone (nROIset = 1, lenSet = 0) — a
%       real measurement, limited by the ROI spacing.
% ONE PLACE FOR EVERY CUT. minLen_CL, strictP_CL and ampRange_CL used to be separate
% variables, so the selection was spread over four places and no single line told you
% what had been removed. They are all rows of extraFilt_CL now:
%     old  minLen_CL   = 35      ->  {'lenSet', [35 inf]}
%     old  ampRange_CL = [0.7 1] ->  {'amp',    [0.7 1]}
%     old  thrCL / strictP_CL    ->  {'pReal',  [0.5 inf]}
% Bounds are INCLUSIVE and a peak whose value is NaN for a constrained variable cannot
% be checked, so it is dropped and counted separately. The full list of variables you
% can filter on is PRINTED when the cell runs, with each one's range over the scored
% peaks, so you never have to guess a sensible bound. The cut-by-cut waterfall is
% printed too: how many peaks each row removes on its own, and how many survive after
% it is applied in order.
%
% Note that constraining a variable the classifier itself used (bandFrac, oscRatio,
% cmCorr...) truncates that variable's own scatter panel by construction — those panels
% are marked "(cut)" and their rho should not be read.
refLen_CL   = 37;       % µm — reference length drawn as a red line (display only)
extRep_CL   = [50 100 150 200 250];  % µm — how many of the filtered events extend
                        % beyond each of these is printed. One number if you only care
                        % about one, [] to skip the report.
capLen_CL   = 250;      % µm — lengths at or above this are AGGREGATED into one
                        % overflow bin at the right edge, labelled ">=capLen_CL". Display
                        % only: nothing is removed and every statistic printed below
                        % (median, IQR, mean) still uses the true unclipped lengths.
                        % Inf or [] = no aggregation, draw the full range.
binLen_CL   = 15;       % µm per histogram bin (display only)
extraFilt_CL = { ...    % EVERY cut, one row per variable: {name, [lo hi]}
    'pReal',   [0.8 inf]      % classifier probability; thrCL = the operating point
    'amp',     [0.2 inf]        % amplitude band (spike height)
    'lenSet',  [5 inf]         % lower bound on the length itself
    'aucArea', [0 inf]          % e.g. drop negative-area peaks
    'dtDirt', [1000 inf]          % e.g. drop negative-area peaks
    };
% --- the ARTIFACT BOUNDARY on the length x duration plane -------------------------
% extraFilt_CL can only cut one variable at a time — a rectangle in the plane. Some
% artifact families are a REGION: long but brief, or brief but broad. Draw that region
% once on the density panel and it becomes a cut like any other, applied to every run
% afterwards and reported in the same waterfall.
%   'load'   use the saved boundary if there is one (default, and what you want once
%            it is drawn). Silently does nothing when no file exists yet.
%   'draw'   draw a NEW boundary on the density panel at the end of this cell and save
%            it. The drawing is interactive: click the vertices, double-click to
%            close, drag them to adjust, then press Enter. The CURRENT figure is not
%            re-filtered — re-run the cell to apply it.
%   'off'    ignore any saved boundary
% Peaks INSIDE the polygon are treated as ARTIFACTS and removed. The file records
% which two parameters the axes were, so a boundary drawn on one plane can never be
% applied to a different one by mistake.
artifactPoly_CL = 'load';
polyFile_CL = fullfile(save_figto, ['ArtifactPoly_' runTag_PK '.mat']);
if ~exist('pReal','var'), error('Run the multi-parameter classifier cell first.'); end
if ~excludeAxes_CL
    warning(['excludeAxes_CL is false, so lenSet was one of the classifier''s own ' ...
        'predictors — this length distribution is circular and should not be read.']);
end
iLn = iP('lenSet');  iAm = iP('amp');  iDu = iP('durSet');

% --- what you can filter on, and what the values look like ---
% pReal is not a column of PkTab, so the two are concatenated into one table here and
% everything downstream indexes that. Add any other derived quantity the same way.
filtNames_CL = [pkNames(:); {'pReal'}];
FiltVal_CL   = [PkTab, pReal(:)];
fprintf('\n=== variables you can filter on (extraFilt_CL), over all %d scored peaks ===\n', nPk);
fprintf('  %-11s %9s %9s %9s %9s %9s %7s\n', 'name','min','p25','median','p75','max','NaN%');
for q = 1:numel(filtNames_CL)
    v = FiltVal_CL(:,q);  vv = v(~isnan(v));
    if isempty(vv)
        fprintf('  %-11s %9s %9s %9s %9s %9s %6.1f%%\n', filtNames_CL{q}, '-','-','-','-','-', 100);
    else
        fprintf('  %-11s %9.3g %9.3g %9.3g %9.3g %9.3g %6.1f%%\n', filtNames_CL{q}, ...
            min(vv), prctile(vv,25), median(vv), prctile(vv,75), max(vv), ...
            100*mean(isnan(v)));
    end
end

% --- apply every cut, and say what each one did ---
inKeep = true(nPk,1);  filtTxt = '';  filtIdx = [];  filtRng = [];
loLen_CL = 0;                                  % histogram x-axis start, from the cuts
fprintf('\n=== filter waterfall ===\n');
fprintf('  %-11s %-18s %9s %9s %9s\n', 'variable','range','alone','NaN','remaining');
fprintf('  %-11s %-18s %9s %9s %9d\n', '(start)','', '', '', nPk);
for kf = 1:size(extraFilt_CL,1)
    nm = extraFilt_CL{kf,1};
    qf = find(strcmp(filtNames_CL, nm), 1);
    if isempty(qf)
        warning('extraFilt_CL: "%s" is not a filterable variable — ignored.', nm);
        continue
    end
    rg = extraFilt_CL{kf,2};
    v  = FiltVal_CL(:,qf);
    isNaNp = isnan(v);
    keep = v >= rg(1) & v <= rg(2);            % NaN -> false
    inKeep = inKeep & keep;
    fprintf('  %-11s %-18s %9d %9d %9d\n', nm, sprintf('[%g %g]', rg(1), rg(2)), ...
        sum(~keep & ~isNaNp), sum(isNaNp), sum(inKeep));
    filtTxt = [filtTxt sprintf(', %s in [%g %g]', nm, rg(1), rg(2))]; %#ok<AGROW>
    pf = find(strcmp(pkNames, nm), 1);         % only PkTab columns get a scatter bound
    if ~isempty(pf), filtIdx(end+1) = pf;  filtRng(end+1,:) = rg(:)'; end %#ok<SAGROW>
    if strcmp(nm,'lenSet'), loLen_CL = max(0, rg(1)); end
end
% --- the artifact boundary, applied as one more cut ---
artPoly = [];  artTxt = '';
if ~strcmpi(artifactPoly_CL,'off') && isfile(polyFile_CL)
    Lp = load(polyFile_CL);
    if isfield(Lp,'artPoly') && isfield(Lp,'xName') && isfield(Lp,'yName')
        if strcmp(Lp.xName, pkNames{iLn}) && strcmp(Lp.yName, pkNames{iDu})
            artPoly = Lp.artPoly;
        else
            warning(['%s was drawn on %s vs %s, not %s vs %s — ignored. Delete it or ' ...
                'redraw with artifactPoly_CL = ''draw''.'], polyFile_CL, Lp.xName, ...
                Lp.yName, pkNames{iLn}, pkNames{iDu});
        end
    end
end
if ~isempty(artPoly)
    inPoly = inpolygon(FiltVal_CL(:,iLn), FiltVal_CL(:,iDu), artPoly(:,1), artPoly(:,2));
    keepA = ~inPoly;
    inKeep = inKeep & keepA;
    fprintf('  %-11s %-18s %9d %9d %9d\n', 'artifact', 'polygon (len x dur)', ...
        sum(~keepA), 0, sum(inKeep));
    artTxt = sprintf(', outside the artifact boundary (%d vertices)', size(artPoly,1));
end
fprintf('  %-11s %-18s %9s %9s %9d  (%.1f%% of all peaks)\n', '(final)','','','', ...
    sum(inKeep), 100*mean(inKeep));
% 'alone' is how many that one cut would remove from the FULL set, so the column does
% not add up to the total: cuts overlap, and the 'remaining' column is the one that
% tracks the selection actually used.

iNR = iP('nROIset');
unmeasured = PkTab(:,iNR) == 0;                        % (a) never reached threshold
confined   = PkTab(:,iNR) == 1;                        % (b) set is one ROI: real, but
                                                       %     below the ROI spacing
% The selection: every cut in extraFilt_CL, applied together.
selShown = inKeep & ~isnan(PkTab(:,iLn));
Lr = PkTab(selShown, iLn);
if isempty(Lr)
    error(['No peak survives extraFilt_CL. The waterfall above shows which row emptied ' ...
        'the set — widen that one.']);
end
% The like-for-like reference: same cuts EXCEPT the classifier probability, split on it.
qP = find(strcmp(filtNames_CL,'pReal'), 1);
kP = find(strcmp(extraFilt_CL(:,1),'pReal'), 1);
inNoP = true(nPk,1);
for kf = 1:size(extraFilt_CL,1)
    if kf == kP, continue; end
    qf = find(strcmp(filtNames_CL, extraFilt_CL{kf,1}), 1);
    if isempty(qf), continue; end
    rg = extraFilt_CL{kf,2};
    inNoP = inNoP & FiltVal_CL(:,qf) >= rg(1) & FiltVal_CL(:,qf) <= rg(2);
end
inNoP = inNoP & ~isnan(PkTab(:,iLn));
if isempty(kP)
    isArtCL = false(nPk,1);
else
    pLo = extraFilt_CL{kP,2}(1);
    isArtCL = inNoP & pReal < pLo;              % failed only the probability cut
end
La = PkTab(isArtCL, iLn);

nInBand = sum(selShown);
fprintf('\n%d peaks are in the histogram. Of those, at length 0:\n', nInBand);
fprintf('   %6d (%.1f%%) UNMEASURABLE — the detecting ROI never reaches %.2f spike height\n', ...
    sum(selShown & unmeasured), 100*sum(selShown & unmeasured)/max(nInBand,1), Vthr_sh);
fprintf('   %6d (%.1f%%) set of ONE ROI — a real measurement, below the ROI spacing\n', ...
    sum(selShown & confined), 100*sum(selShown & confined)/max(nInBand,1));
fprintf('   %6d (%.1f%%) had another set elsewhere on the dendrite at the same instant\n', ...
    sum(selShown & PkTab(:,iP('nSetEv'))>1), ...
    100*sum(selShown & PkTab(:,iP('nSetEv'))>1)/max(nInBand,1));
fprintf('  median %.0f µm (IQR %.0f-%.0f), mean %.0f +/- %.0f\n', ...
    median(Lr), prctile(Lr,25), prctile(Lr,75), mean(Lr), std(Lr));
fprintf('  below %g µm: %d (%.1f%%)   |   at or above: %d (%.1f%%)\n', ...
    refLen_CL, sum(Lr<refLen_CL), 100*mean(Lr<refLen_CL), ...
    sum(Lr>=refLen_CL), 100*mean(Lr>=refLen_CL));
% How far the events that survived every cut actually reach. Strictly greater than, so
% "> 200" excludes an event measured at exactly 200 µm.
if ~isempty(extRep_CL)
    fprintf('  of the %d filtered events, how many extend beyond:\n', numel(Lr));
    for xt = extRep_CL(:)'
        nx = sum(Lr > xt);
        fprintf('     > %4g µm : %5d  (%.1f%%)%s\n', xt, nx, 100*nx/numel(Lr), ...
            repmat_str(nx == 0, '   <- none', ''));
    end
end
if ~isempty(capLen_CL) && isfinite(capLen_CL)
    nOv = sum(Lr >= capLen_CL);
    if nOv > 0
        fprintf(['  histogram overflow bin: %d peaks (%.1f%%) are >= %g µm and are drawn ' ...
            'in one bar.\n     Their real lengths span %.0f-%.0f µm (median %.0f) and are ' ...
            'INCLUDED in every statistic above.\n'], nOv, 100*nOv/numel(Lr), capLen_CL, ...
            min(Lr(Lr>=capLen_CL)), max(Lr), median(Lr(Lr>=capLen_CL)));
    else
        fprintf('  histogram overflow bin: no peak reaches %g µm, so it is empty\n', capLen_CL);
    end
end
if ~isempty(La)
    aucLen = param_auc(PkTab(:,iLn), selShown, isArtCL);
    fprintf(['  (peaks failing ONLY the pReal cut are not drawn; for reference n = %d, ' ...
        'median %.0f µm, and length alone separates the two groups with AUC %.3f)\n'], ...
        numel(La), median(La), aucLen);
end
% How sensitive is the answer to where the probability cut sits? Reported for a few
% levels rather than through a separate strictP_CL variable.
if ~isempty(kP)
    fprintf('  sensitivity to the pReal cut (all other cuts unchanged):\n');
    for pc = unique([extraFilt_CL{kP,2}(1) 0.5 0.75 0.9])
        Lp = PkTab(inNoP & pReal >= pc, iLn);
        if isempty(Lp), continue; end
        fprintf('     p >= %.2f : n = %6d, median %3.0f µm, %.1f%% below %g\n', ...
            pc, numel(Lp), median(Lp), 100*mean(Lp<refLen_CL), refLen_CL);
    end
end
fprintf('  hand-picked candidates: length %s µm at amplitude %s\n', ...
    num2str(exPk(:,iLn)','%.0f  '), num2str(exPk(:,iAm)','%.2f  '));

% ===== figure 154: the length distribution, plus the scatters YOU ask for =====
% Every point in this figure is one peak that survived extraFilt_CL — exactly the
% peaks counted in the waterfall above. selShown is that one selection, used
% throughout, so no panel can be showing a different set from another.
%
% Panel 1 is the length histogram and panel 2 takes ONE peak apart to show how its
% lenSet was measured. After them comes one DENSITY HEATMAP per row of scatterXY_CL
% and nothing else.
%   colour     = peaks per 2-D bin (log10 by default); empty bins are left white, so
%                the support of the distribution is visible instead of being painted
%                as "low"
%   blue line  = running median of y in bins of x (the trend through the density)
%   black dash = the passive-cable prediction, on a length x duration panel only
%   red line   = refLen_CL, drawn only on an axis that IS length
%   magenta    = an extraFilt_CL bound on that axis — the cloud is cut there BY
%                CONSTRUCTION, so the rho of a cut parameter should not be read
%   rho        = Spearman correlation between the two plotted parameters
% A length axis reuses the histogram's bin edges, overflow bin included, so the two
% panels line up; any other axis is binned over its 0.5-99.5 percentile range.
%
% scatterXY_CL: one row per scatter, {xParameter, yParameter}. Any two names from
% pkNames, in any order — 'amp' 'lenSet' 'durSet' 'bandRatio' 'bandFrac' 'fLoc'
% 'oscRatio' 'cmCorr' 'bvFrac' 'roiPkRate' 'distSoma' 'nearSpk' 'aucSum' 'aucArea'
% 'nROIset' 'nSetEv'.
% Give it {} to draw the histogram alone.
scatterXY_CL = {'lenSet','durSet'};       % e.g. {'amp','lenSet'; 'lenSet','durSet'}
densMk_CL    = 8;       % marker size in the scatter_density panels. Every peak is
                        % drawn at its own coordinates and coloured by local kernel
                        % density — nothing is binned, so nothing saturates and the
                        % long tail stays visible. capLen_CL applies to the HISTOGRAM
                        % only; these axes follow the data.
exampleSel_CL = [];  % which peak the worked example takes apart: 'median' (the
                        % one closest to the median length, with >=3 ROIs), 'longest',
                        % or a row index into PkTab to show one specific peak
kymoHalfEx_CL = 300;    % +/- ms shown in the example's kymograph and trace
setName_CL  = 'depolarized ROIs';   % what the coloured ROIs are CALLED in the map
                        % legend, and the term to use in the figure caption. They are
                        % the ROIs above Vthr_sh that are CONTIGUOUS with the peak
                        % along the dendrite, so "above threshold" on its own would
                        % misdescribe them: a far ROI above threshold but not connected
                        % is deliberately excluded. Alternatives that keep that meaning:
                        % 'event ROIs', 'involved ROIs', 'depolarized domain'.
labelROI_CL = false;    % print each ROI's index on its footprint. Useful only when
                        % cross-checking one ROI against the profile panel or PkMeta.
% --- the passive-cable comparison drawn on a length x duration panel ---
cableD_CL   = 807;      % µm^2/ms. For a passive cable D = lambda^2 / tau, so e.g.
                        % lambda 200 µm with tau 50 ms gives ~800.  [] = no line.
cableFun_CL = @(t,D) 2*sqrt(2*D*t);   % FULL extent reached after diffusing for t ms.
                        % The factor 2 is because lenSet is the full span and an event
                        % spreads both ways from its trigger, while sqrt(2*D*t) is the
                        % one-sided RMS displacement. Change this line if you want a
                        % different convention, e.g. @(t,D) sqrt(4*D*t) — the panel
                        % legend and the printed table follow whatever is written here.
cReal = [0.2 0.2 0.2];  cTrend = [0 0.45 0.74];
xL = PkTab(selShown, iLn);
topL = max([prctile(Lr,99.5), 2*refLen_CL, loLen_CL+binLen_CL]);
% no subsampling: a density heatmap uses every point by construction, which is one of
% the reasons to prefer it over a scatter here.

% Spearman rho of every parameter against length, for the table below. This used to be
% a side effect of drawing the panels; it is computed here so the table stays complete
% however few panels are drawn.
rhoP = nan(nP,1);
for pj = 1:nP
    yP = PkTab(selShown, pj);  gd = ~isnan(yP) & ~isnan(xL);
    if sum(gd) > 2, rhoP(pj) = corr(xL(gd), yP(gd), 'Type','Spearman'); end
end

if isempty(scatterXY_CL), xyPairs = cell(0,2); else, xyPairs = scatterXY_CL; end
if size(xyPairs,2) ~= 2
    error('scatterXY_CL must be an N x 2 cell array of {xParameter, yParameter} rows.');
end
nSc = size(xyPairs,1);

% Layout: the histogram and the density panels are wide (2 columns each) and sit on
% the top row(s); the worked example's four panels fill the row underneath.
figure(154); clf;
ncP = 4;
nWide = 1 + nSc;                                 % histogram + one per density panel
nrTop = ceil(nWide/2);
nrP = nrTop + 1;
tiledlayout(nrP, ncP, 'TileSpacing','compact', 'Padding','compact');

% --- panel 1: THE HISTOGRAM — length of every peak that survived the filters ---
% Everything at or above capLen_CL goes into ONE overflow bin drawn at the right edge.
% The long tail is usually a handful of peaks spread over hundreds of µm, which
% stretches the axis and squashes the part you actually want to read. Only the drawing
% is clipped: Lr itself is untouched, so the median, IQR and mean in the title and the
% printout are the true values, and the overflow count is stated explicitly.
nexttile([1 2]); hold all;
doCap_CL = ~isempty(capLen_CL) && isfinite(capLen_CL) && capLen_CL > loLen_CL;
if doCap_CL
    edL = loLen_CL:binLen_CL:capLen_CL;
    if edL(end) < capLen_CL
        % the cap does not sit on the bin grid, so the bin before it is narrower than
        % the rest. With 'probability' normalization a narrow bin holds proportionally
        % fewer peaks and therefore draws SHORTER, which reads as a dip that is not in
        % the data. Say so rather than letting it pass as a feature of the histogram.
        fprintf(['  NOTE: capLen_CL %g is not on the bin grid (start %g, width %g), so ' ...
            'the bin\n     %g-%g µm is only %g µm wide and draws short. Use binLen_CL ' ...
            'that divides %g, e.g. %s.\n'], capLen_CL, loLen_CL, binLen_CL, edL(end), ...
            capLen_CL, capLen_CL-edL(end), capLen_CL-loLen_CL, ...
            strjoin(compose('%g', divisors_near(capLen_CL-loLen_CL, binLen_CL)), ' or '));
        edL(end+1) = capLen_CL;
    end
    edL(end+1) = capLen_CL + binLen_CL;                    % the overflow bin
    Lplot = min(Lr, capLen_CL);                            % display only
    nOver = sum(Lr >= capLen_CL);
else
    edL = loLen_CL:binLen_CL:ceil(topL/binLen_CL)*binLen_CL;
    Lplot = Lr;  nOver = 0;
end
histogram(Lplot, edL, 'Normalization','probability', 'FaceColor',cReal, ...
    'EdgeColor','none', 'FaceAlpha',0.85);
yl = ylim;
% for i = 1:nEx
%     if ~isnan(exPk(i,iLn))
%         plot(min(exPk(i,iLn), edL(end-1))*[1 1], yl, 'g-', 'linewidth',1);
%     end
% end
plot(refLen_CL*[1 1], yl, 'r-', 'linewidth',2);
text(refLen_CL, yl(2), sprintf(' %g µm', refLen_CL), 'Color','r', ...
    'VerticalAlignment','top', 'FontWeight','bold');
if doCap_CL
    % mark the overflow bin and label the tick as ">=cap" so the bar cannot be read
    % as "a lot of peaks at exactly 200 µm"
    plot(capLen_CL*[1 1], yl, '-', 'Color',[0.5 0.5 0.5], 'linewidth',1);
    xt = get(gca,'XTick');  xt = xt(xt < capLen_CL);
    xt = [xt, capLen_CL + binLen_CL/2];
    xtl = [compose('%g', xt(1:end-1)), {sprintf('\\geq%g', capLen_CL)}];
    set(gca, 'XTick',xt, 'XTickLabel',xtl);
end
ylim(yl);  xlim([loLen_CL edL(end)]);
xlabel('lenSet: full length of the ROI set (µm)');  ylabel('fraction');
ttlExtra = '';
if doCap_CL && nOver > 0
    ttlExtra = sprintf(', %d (%.1f%%) \\geq%g', nOver, 100*nOver/numel(Lr), capLen_CL);
end
title(sprintf('LENGTH  n=%d  med %.0f µm, %.0f%% < %g%s', numel(Lr), median(Lr), ...
    100*mean(Lr<refLen_CL), refLen_CL, ttlExtra), 'FontWeight','bold');
box off;

% --- the density panels listed in scatterXY_CL ---
% scatter_density: every peak is drawn at its OWN coordinates and coloured by the local
% kernel density, so nothing is binned, nothing saturates, and the long tail is visible
% instead of being folded into an overflow column. capLen_CL applies to the histogram
% only; here the axis follows the data, which is also what lets the cable line show.
for j = 1:nSc
    px = find(strcmp(pkNames, xyPairs{j,1}), 1);
    py = find(strcmp(pkNames, xyPairs{j,2}), 1);
    if isempty(px) || isempty(py)
        warning('scatterXY_CL row %d: "%s" vs "%s" — not a parameter name, skipped.', ...
            j, xyPairs{j,1}, xyPairs{j,2});
        continue
    end
    axD = nexttile([1 2]); hold(axD,'on');
    xP = PkTab(selShown, px);  yP = PkTab(selShown, py);
    gd = ~isnan(xP) & ~isnan(yP);
    if sum(gd) < 3
        axis(axD,'off');
        title(axD, sprintf('%s vs %s — no data', pkNames{py}, pkNames{px}));
        continue
    end
    rhoXY = corr(xP(gd), yP(gd), 'Type','Spearman');
    xv = xP(gd);  yv = yP(gd);
    axes(axD);
    scatter_density(xv, yv, densMk_CL, turbo(256));
    cbD = findobj(gcf,'Type','ColorBar');
    if ~isempty(cbD), cbD(1).Label.String = 'local density'; end
    % axis limits follow the data, clipped at the 99.5th percentile so a couple of
    % extreme peaks cannot squash everything else; the count outside is reported
    xhi = max(prctile(xv,99.5), min(xv)+eps);  yhi = max(prctile(yv,99.5), min(yv)+eps);
    xlo = min(0, min(xv));  ylo = min(0, min(yv));
    nOut = sum(xv > xhi | yv > yhi);
    % --- the passive cable / diffusion prediction ---
    % D = lambda^2/tau for a passive cable; cableFun_CL turns a duration into the
    % extent such a spread reaches. The convention (factor of 2 or not) is written out
    % in cableFun_CL — change it there, and the legend and printout follow.
    % isLD = ~isempty(cableD_CL) && cableD_CL > 0 && ...
    %        ((px == iLn && py == iDu) || (py == iLn && px == iDu));
    % if isLD
    %     if px == iLn
    %         tt = linspace(max(ylo,eps), yhi, 400);
    %         cx = cableFun_CL(tt, cableD_CL);  cy = tt;
    %     else
    %         tt = linspace(max(xlo,eps), xhi, 400);
    %         cx = tt;  cy = cableFun_CL(tt, cableD_CL);
    %     end
    %     inBox = cx >= xlo & cx <= xhi & cy >= ylo & cy <= yhi;
    %     plot(axD, cx, cy, 'w-', 'linewidth',3.5, 'HandleVisibility','off');
    %     hC = plot(axD, cx, cy, 'k--', 'linewidth',2);
    %     if any(inBox)
    %         legend(hC, sprintf('cable, D = %g µm^2/ms', cableD_CL), 'box','off', ...
    %             'Location','northwest', 'FontSize',8, 'AutoUpdate','off');
    %         fprintf('     cable line visible over %s = %.0f-%.0f\n', pkNames{px}, ...
    %             min(cx(inBox)), max(cx(inBox)));
    %     else
    %         delete(hC);
    %         text(axD, xhi, yhi, sprintf(' cable D=%g is off-scale here ', cableD_CL), ...
    %             'HorizontalAlignment','right', 'VerticalAlignment','top', 'FontSize',7);
    %     end
    % end
    % if px == iLn, plot(axD, refLen_CL*[1 1], [ylo yhi], 'r-', 'linewidth',1.5, 'HandleVisibility','off'); end
    % if py == iLn, plot(axD, [xlo xhi], refLen_CL*[1 1], 'r-', 'linewidth',1.5, 'HandleVisibility','off'); end
    isCut = false;
    kf = find(filtIdx == px, 1);
    if ~isempty(kf)
        isCut = true;
        for e = filtRng(kf,:)
            if isfinite(e) && e > xlo && e < xhi
                plot(axD, e*[1 1], [ylo yhi], '--', 'color',[0.6 0 0.6], ...
                    'linewidth',1.2, 'HandleVisibility','off');
            end
        end
    end
    kf = find(filtIdx == py, 1);
    if ~isempty(kf)
        isCut = true;
        for e = filtRng(kf,:)
            if isfinite(e) && e > ylo && e < yhi
                plot(axD, [xlo xhi], e*[1 1], '--', 'color',[0.6 0 0.6], ...
                    'linewidth',1.2, 'HandleVisibility','off');
            end
        end
    end
    % the artifact boundary: outline what was removed, on the plane it was drawn on
    isLDplane = (px == iLn && py == iDu);
    if isLDplane && ~isempty(artPoly)
        plot(axD, artPoly([1:end 1],1), artPoly([1:end 1],2), '-', ...
            'Color',[0.85 0.33 0.1], 'linewidth',2, 'HandleVisibility','off');
        patch(axD, 'XData',artPoly(:,1), 'YData',artPoly(:,2), ...
            'FaceColor',[0.85 0.33 0.1], 'FaceAlpha',0.10, 'EdgeColor','none', ...
            'HandleVisibility','off');
        text(axD, mean(artPoly(:,1)), mean(artPoly(:,2)), 'artifact', ...
            'Color',[0.6 0.2 0.05], 'FontSize',8, 'FontWeight','bold', ...
            'HorizontalAlignment','center');
    end
    if isLDplane, axArt_CL = axD; end     % remembered for the draw step below
    xlim(axD, [xlo xhi]);  ylim(axD, [ylo yhi]);
    xlabel(axD, pkNames{px});  ylabel(axD, pkNames{py});
    ttlCut = repmat_str(isCut, ' (cut)', '');
    title(axD, sprintf('%s vs %s%s   \\rho = %+.2f   n = %d', pkNames{py}, pkNames{px}, ...
        ttlCut, rhoXY, sum(gd)), 'FontWeight','normal', ...
        'Color', repmat_str(isCut, [0.6 0 0.6], [0 0 0]));
    box(axD,'off');
    fprintf('  density %s vs %s: n = %d, Spearman rho = %+.3f, %d peaks beyond the 99.5%% axis limit\n', ...
        pkNames{py}, pkNames{px}, sum(gd), rhoXY, nOut);
    % recomputed here rather than reused from the drawing block above, so the cable
    % numbers are still reported even when that block is commented out
    isLD = ~isempty(cableD_CL) && cableD_CL > 0 && ...
           ((px == iLn && py == iDu) || (py == iLn && px == iDu));
    if isLD
        dSam = [10 25 50 100 200];
        fprintf('     cable D = %g µm^2/ms predicts length %s µm at %s ms\n', ...
            cableD_CL, num2str(cableFun_CL(dSam, cableD_CL), '%.0f  '), ...
            num2str(dSam, '%d  '));
        if px == iLn, lOK = xv;  dOK = yv;  else, lOK = yv;  dOK = xv;  end
        fprintf('     %.1f%% of peaks are SHORTER than that line at their own duration\n', ...
            100*mean(lOK < cableFun_CL(max(dOK,eps), cableD_CL)));
    end
end
if mod(nWide,2) == 1                             % keep the example on its own row
    axPad = nexttile([1 2]);  axis(axPad,'off');
end

% --- the worked example: FOUR panels showing how this peak's lenSet was measured ---
% One peak, taken apart. Re-measured through the same peak_set_measure the scoring
% pass called, so these are the numbers that went into panel 1.
kEx = pick_example_CL(exampleSel_CL, find(selShown), PkTab(:,iLn), ...
    PkTab(:,iP('nROIset')), PkMeta);
axE = gobjects(1,4);
for q = 1:4, axE(q) = nexttile([1 1]); end
if isempty(kEx)
    for q = 1:4, axis(axE(q),'off'); end
    title(axE(1),'no example available');
else
    cfgEx = struct('Vthr',Vthr_sh, 'ciZ',ciZ_PK, 'durMax',durMax_PK, ...
        'avgwin',avgwin_ext, 'repMode',repTrace_PK, 'baseWin',baseWin_PK, ...
        'band',bpLenDur_PK, 'bpOrd',bpOrd_PK, 'adjMax',adjMax_PK, ...
        'fastBase',useFastBase_PK, 'kymoHalf',kymoHalfEx_CL, ...
        'setName',setName_CL, 'labelROIs',labelROI_CL);
    [lEx, dEx, nrEx] = draw_example_event(axE(1), axE(2), axE(3), axE(4), ...
        PkMeta(kEx,1), PkMeta(kEx,2), PkMeta(kEx,3), PCresult, cfgEx);
    fprintf(['\nExample in figure 154 (4 panels): f%d ROI%d frame%d -> lenSet %.0f µm ' ...
        'over %d ROIs, durSet %.0f ms\n'], PkMeta(kEx,1), PkMeta(kEx,2), PkMeta(kEx,3), ...
        lEx, nrEx, dEx);
end

sgtitle(sprintf('%d peaks after%s%s', numel(Lr), filtTxt, artTxt), 'FontSize',11);
set_font('Arial'); set_fontsize(9); %set_figsize(240*ncP, 170*nrP);

% --- draw and save a new artifact boundary ---
% Interactive, so it only happens when you ask for it. The peaks on the panel are the
% ones that survived the OTHER cuts, which is the right background to judge against.
% The current figure is deliberately not re-filtered: re-run the cell so the boundary
% goes through the waterfall like every other cut and the counts stay honest.
if strcmpi(artifactPoly_CL,'draw')
    if ~exist('axArt_CL','var') || ~isgraphics(axArt_CL)
        warning(['artifactPoly_CL = ''draw'' needs a %s vs %s panel to draw on — put ' ...
            'that pair in scatterXY_CL.'], pkNames{iLn}, pkNames{iDu});
    elseif ~exist('drawpolygon','file')
        warning(['drawpolygon is not available (Image Processing Toolbox). Define ' ...
            'artPoly yourself and save it to %s with xName/yName.'], polyFile_CL);
    else
        figure(154);  drawnow;  axes(axArt_CL);
        fprintf(['\n>>> DRAW THE ARTIFACT BOUNDARY on the %s vs %s panel of figure 154.\n' ...
            '    Click the vertices; double-click (or click the first vertex) to close.\n' ...
            '    Drag any vertex to adjust, then press Enter in this window to save.\n'], ...
            pkNames{iLn}, pkNames{iDu});
        roi = drawpolygon(axArt_CL, 'Color',[0.85 0.33 0.1], 'FaceAlpha',0.12, ...
            'LineWidth',2);
        input('    press Enter when the boundary is where you want it ', 's');
        if isvalid(roi) && size(roi.Position,1) >= 3
            artPoly = roi.Position;
            xName = pkNames{iLn};  yName = pkNames{iDu};   %#ok<NASGU>
            drawnWhen = datetime('now');                   %#ok<NASGU>
            drawnFilt = filtTxt;                           %#ok<NASGU>
            save(polyFile_CL, 'artPoly','xName','yName','drawnWhen','drawnFilt');
            nIn = sum(inpolygon(PkTab(:,iLn), PkTab(:,iDu), artPoly(:,1), artPoly(:,2)));
            fprintf(['    saved %d vertices to %s\n    it encloses %d of %d scored ' ...
                'peaks (%.1f%%).\n    Set artifactPoly_CL = ''load'' and re-run this ' ...
                'cell to apply it.\n'], size(artPoly,1), polyFile_CL, nIn, nPk, ...
                100*nIn/nPk);
        else
            fprintf('    no polygon with at least 3 vertices — nothing saved.\n');
        end
    end
elseif ~isempty(artPoly)
    fprintf('\nArtifact boundary: %d vertices from %s, removing %d peaks.\n', ...
        size(artPoly,1), polyFile_CL, sum(inPoly));
elseif ~strcmpi(artifactPoly_CL,'off')
    fprintf(['\nNo artifact boundary saved yet. Set artifactPoly_CL = ''draw'' and ' ...
        're-run this cell to draw one.\n']);
end

% --- the same peaks, as a table, with how strongly each tracks length ---
% Every parameter appears here whether or not it was plotted: dropping the scatter grid
% removed the pictures, not the numbers, so nothing stops being checked.
ordP = [iLn, setdiff(1:nP, iLn, 'stable')];        % length first, then the rest
fprintf('\n  parameters of those %d peaks (rho = Spearman correlation WITH length):\n', numel(Lr));
fprintf('  %-11s %9s %9s %9s %8s %9s\n', 'parameter','median','p25','p75','rho','n');
for j = 1:numel(ordP)
    pj = ordP(j);
    v = PkTab(selShown, pj);  v = v(~isnan(v));
    if isempty(v)
        fprintf('  %-11s %9s %9s %9s %8s %9d\n', pkNames{pj}, '-','-','-','-', 0);
    else
        fprintf('  %-11s %9.3g %9.3g %9.3g %8s %9d\n', pkNames{pj}, median(v), ...
            prctile(v,25), prctile(v,75), ...
            repmat_str(isnan(rhoP(pj)), '-', sprintf('%+.2f',rhoP(pj))), numel(v));
    end
end
[~, ordR] = sort(abs(rhoP), 'descend', 'MissingPlacement','last');
fprintf('  strongest length dependence: %s\n', strjoin(arrayfun(@(k) ...
    sprintf('%s (%+.2f)', pkNames{k}, rhoP(k)), ordR(1:min(3,sum(~isnan(rhoP)))), ...
    'UniformOutput',false), ', '));
ufN  = unique(PkMeta(selShown,1))';
medN = arrayfun(@(f) median(PkTab(selShown & PkMeta(:,1)==f, iLn),'omitnan'), ufN);
fprintf(['  per-neuron check: median of the %d per-neuron median lengths %.0f µm ' ...
    '(spread %.0f-%.0f) vs pooled %.0f µm\n'], numel(ufN), median(medN), min(medN), ...
    max(medN), median(Lr));

%% ===== INTERACTIVE length x duration: click a dot, see its kymograph =====
% FIGURE 154 MADE CLICKABLE. This runs straight after the length cell and browses
% EXACTLY the peaks that survived extraFilt_CL — the same selShown that panel 1's
% histogram and the density panel are built from — so the plane you click on is the
% plane you just read, not a different population. The window opens beside figure 154.
%
% Left panel: lenSet vs durSet for the filtered peaks. Click anywhere and the nearest
% drawn peak is selected, all parameters appear in the status bar, and one pasteable
% line per selection goes to the command window.
%
% The four panels on the right SHOW THE MEASUREMENT. Each one is recomputed for the
% clicked peak by calling the same peak_set_measure the scoring pass called, with the
% settings passed in cfgLD below, so these are the numbers in PkTab and not a
% lookalike. (If they ever disagree — because a scoring parameter was edited without
% re-scoring — the status bar says so instead of quietly showing two answers.)
%   kymograph   every ROI by distance from the soma; white dots on the left edge mark
%               the members of the detecting ROI's SET, i.e. the ROIs lenSet spans
%   trace       the SET'S representative ROI: unfiltered in grey, the 0.1-60 Hz trace
%               the crossings were found on in black, the threshold as a red dotted
%               line, and the shaded span between the crossings, which IS durSet
%   set map     the morphology painted by set membership — grey below threshold, one
%               colour per set, and a magenta outline on any single-ROI gap that was
%               tested (filled = not trusted, so the sets either side were merged)
%   profile     amplitude vs position along the dendrite with every ROI's CI bar, the
%               sets coloured and their spans drawn: LENGTH is read off this panel
%   click            select the nearest peak to the pointer
%   <- / ->          step out from / back towards the point you clicked. This is how you
%                    get through the peaks stacked on the same dot: they come out
%                    nearest-first, so repeated -> walks the pile.
%   up / down        jump 10 at a time
%   p                re-print the current peak
%   c                copy "[f frame]" to the clipboard, ready for exampleList
%   x                PICK THIS PEAK for the worked example: sets exampleSel_CL in the
%                    base workspace, then re-run the length cell and figure 154's four
%                    example panels take this peak apart
% Only the drawn subsample is selectable — a click cannot reach a peak that is not on
% screen. Raise nShowBrowse_LD to draw more (slower to render, no other cost).
% Set useFilter_LD = false to browse every scored peak instead of the filtered set,
% and fShow_LD to restrict it to particular neurons.
% Kymographs are cut ON DEMAND, so the number of peaks browsed is unlimited and nothing
% is precomputed; the three most recently visited neurons stay cached, so clicking
% around within one neuron is instant and moving to another costs one normalisation.
% This cell only READS PkTab/PkMeta — it writes nothing and changes no labels.
% NOTE: MATLAB's zoom and pan tools swallow mouse clicks. Zoom in, then click the tool
% off (or press the arrow cursor in the toolbar) to go back to selecting.
useFilter_LD = true;       % true = browse only the peaks in figure 154 (selShown);
                           % false = every scored peak, as before
fShow_LD = [20 23];             % WHICH NEURONS to browse, e.g. [20 23]. [] = all of them.
                           % Named _LD rather than f2show because f2show is a live
                           % variable in several other cells of this file and reusing
                           % it here would silently change what they plot.
nShowBrowse_LD = 30000;    % peaks drawn, and therefore clickable
rngSeed_LD  = 1;           % the subsample and the jitter are both reproducible
jitter_LD   = [4 1];       % [um ms] uniform +/- jitter added to the DOT POSITIONS ONLY.
                           % lenSet is quantised by the ROI spacing and durSet by the
                           % 1 kHz frame rate, so thousands of peaks land on exactly the
                           % same coordinates — most of them at (0,0) — and without
                           % jitter only one of them could ever be picked. Every number
                           % printed or shown comes from PkTab and is never jittered.
kymoHalf_LD = 500;         % +/- ms shown in the kymograph and the trace
colorBy_LD  = 'auto';      % 'auto' = p(real) if the classifier cell has run, else
                           % 'amp'. Any parameter name also works, e.g. 'aucArea'.
if ~exist('PkTab','var'), error('Run the all-peaks scoring cell first.'); end
if useFilter_LD && ~exist('selShown','var')
    error(['useFilter_LD is true but selShown does not exist — run the length ' ...
        'distribution cell (figure 154) first, or set useFilter_LD = false.']);
end
iL = iP('lenSet');  iD = iP('durSet');

% --- which peaks to draw: the SAME selection figure 154 is built from ---
rng(rngSeed_LD);
if useFilter_LD
    selLD = find(selShown & ~isnan(PkTab(:,iL)) & ~isnan(PkTab(:,iD)));
    filtTxt_LD = filtTxt;
else
    selLD = find(~isnan(PkTab(:,iL)) & ~isnan(PkTab(:,iD)));
    filtTxt_LD = ' (no filter — every scored peak)';
end
if ~isempty(fShow_LD)
    avail = unique(PkMeta(selLD,1))';
    keepN = ismember(PkMeta(selLD,1), fShow_LD);
    if ~any(keepN)
        error(['fShow_LD = %s selects no peak. Neurons present in this selection: ' ...
            '%s.'], mat2str(fShow_LD), mat2str(avail));
    end
    missing = setdiff(fShow_LD, avail);
    if ~isempty(missing)
        warning('fShow_LD: no peak from neuron(s) %s in this selection — skipped.', ...
            mat2str(missing));
    end
    selLD = selLD(keepN);
    filtTxt_LD = [filtTxt_LD sprintf(', neurons %s', mat2str(intersect(fShow_LD, avail)))];
end
if isempty(selLD)
    error('No peak to browse. Widen extraFilt_CL or set useFilter_LD = false.');
end
if numel(selLD) > nShowBrowse_LD
    selLD = sort(selLD(randperm(numel(selLD), nShowBrowse_LD)));
end
uf = unique(PkMeta(selLD,1))';
fprintf('Browsing %d peaks from %d neuron(s): %s\n', numel(selLD), numel(uf), ...
    strjoin(compose('f%d (%d)', uf(:), arrayfun(@(f) sum(PkMeta(selLD,1)==f), uf(:))), ', '));

% --- what to colour them by ---
cLimLD = [];
switch lower(colorBy_LD)
    case {'auto','preal','p(real)'}
        if exist('pReal','var') && numel(pReal) == nPk
            cLD = pReal(selLD);  cLabLD = 'p(real)';  cLimLD = [0 1];
        else
            cLD = PkTab(selLD,iP('amp'));  cLabLD = 'amp (spike height)';
        end
    otherwise
        pc = iP(colorBy_LD);
        if isempty(pc), error('colorBy_LD "%s" is not a parameter name.', colorBy_LD); end
        cLD = PkTab(selLD,pc);  cLabLD = colorBy_LD;
end

% --- hand labels and classifier scores, so a selected peak can report them ---
labAllLD = nan(nPk,1);
[liLD, lvLD] = load_peak_labels({labFile_PK, importFrom_PK}, PkMeta);
if ~isempty(liLD), labAllLD(liLD) = lvLD; end
pAllLD = nan(nPk,1);
if exist('pReal','var') && numel(pReal) == nPk, pAllLD = pReal; end

% The browser re-measures each clicked peak with event_len_dur, so it needs the same
% settings the scoring pass used. Passing them explicitly (rather than letting the
% viewer pick its own) is what makes the panels the measurement and not a picture of
% something similar; if you edit any of these without re-scoring, the viewer notices
% and says so on the peak it redraws.
if isempty(bpLenDur_PK)
    bpShort_LD = 'unfiltered';
else
    bpShort_LD = sprintf('%g-%g Hz', bpLenDur_PK(1), bpLenDur_PK(2));
end
cfgLD = struct('tAxis', -kymoHalf_LD:kymoHalf_LD, 'pkNames',{pkNames}, 'iL',iL, 'iD',iD, ...
    'Vthr',Vthr_sh, 'jitter',jitter_LD, 'jitterSeed',rngSeed_LD, 'cLab',cLabLD, ...
    'cLim',cLimLD, 'aucWin',aucWin_PK, 'exPk',exPk, 'labAll',labAllLD, 'pAll',pAllLD, ...
    'avgwin',avgwin_ext, 'Dmax',Dmax_ext, 'baseWin',baseWin_PK, 'durMax',durMax_PK, ...
    'smoothDur',smoothDurEff_PK, 'band',bpLenDur_PK, 'bpOrd',bpOrd_PK, ...
    'bpTxtShort',bpShort_LD, 'ciZ',ciZ_PK, 'repMode',repTrace_PK, 'adjMax',adjMax_PK, ...
    'fastBase',useFastBase_PK, 'nextTo',154, 'filtTxt',filtTxt_LD, ...
    'setName',setName_CL, 'labelROIs',labelROI_CL);
browse_len_dur(PkTab, PkMeta, PCresult, selLD, cLD, cfgLD);
fprintf(['Interactive browser open (figure 160, beside figure 154): %d peaks drawn' ...
    '%s,\n  coloured by %s. Click a dot; arrows walk through its neighbours; ' ...
    'c copies the identity for exampleList.\n'], numel(selLD), filtTxt_LD, cLabLD);

%% ===== Amplitude-matched control selection (SAME neuron) =====
% For each local candidate, pick control events from the SAME neuron with peak
% amplitude within +/- ampTol_ext and the same broad compartment (basal vs
% apical). If too few, relax to amplitude-only. Keep up to nCtrl_ext controls
% (closest in amplitude). These are amplitude-matched "broad" comparison events.
ctrlRows = cell(nEx,1);          % EventTab row indices of matched controls per candidate
for i = 1:nEx
    fe = exampleList(i,1);  te = exampleList(i,2);  a = exAMP(i);
    same  = EventTab(:,1)==fe & abs(EventTab(:,3)-te)>200;        % same neuron, not the candidate
    ampok = EventTab(:,4)>=a*(1-ampTol_ext) & EventTab(:,4)<=a*(1+ampTol_ext);
    m = find(same & ampok & EventTab(:,6)==exCOMP(i));           % + same compartment
    if numel(m) < 5, m = find(same & ampok); end                 % relax if too few
    [~, ord] = sort(abs(EventTab(m,4)-a));                        % closest in amplitude first
    ctrlRows{i} = m(ord(1:min(nCtrl_ext,numel(m))));
    fprintf('Candidate %d (f%d): %d amplitude-matched controls (same neuron)\n', i, fe, numel(ctrlRows{i}));
end

% profiles: candidates (have exProf) + pooled matched controls (recompute)
allCtrlRows = unique(cell2mat(ctrlRows));
ctrlProfAll = get_profiles(allCtrlRows, EventTab, PCresult, avgwin_ext, Dmax_ext);

% one example broad control = matched control with the largest L50
[~, bi] = max(EventTab(allCtrlRows,5));      % col 5 = L50
broadRow = allCtrlRows(bi);
fprintf('Example broad control: f%d frame %d ROI %d, L50=%.0f um\n', ...
    EventTab(broadRow,1), EventTab(broadRow,3), EventTab(broadRow,2), EventTab(broadRow,5));

%% ===== Figure 142 — aggregate (A/B top, C bottom-left, D bottom-right) =====
% A: example local candidate; B: example amplitude-matched broad event (each =
% morphology colored by normalized amplitude + amplitude-vs-distance profile).
% C: mean +/- sem normalized amplitude vs contour distance, local candidates
% (red) vs same-neuron amplitude-matched controls (gray). D: histogram of L50
% over the matched controls, with each candidate marked. E: histogram of L50 over
% ALL events (all neurons) with the candidates' mean L50 marked.
globalThresh = 50;    % L50 (um) above which an event is called global/broad
L50min       = 20;     % L50 (um) below this is likely blood-vessel noise -> excluded

% pooled L50 over all neurons (exclude L50 < L50min) + total recording time
allL50 = EventTab(:,5);  allL50 = allL50(~isnan(allL50) & allL50 >= L50min);
totalSec = 0;
for f = foi, totalSec = totalSec + size(PCresult{f}.Subthreshold,2)/1000; end
nGlobal = sum(allL50 > globalThresh);
rate10k = nGlobal / totalSec * 1e4;

% candidate mean amplitude in dF/F ( = spike-height units x per-neuron SpikeHeight )
exAMP_dff = arrayfun(@(i) exAMP(i) * PCresult{exampleList(i,1)}.SpikeHeight, 1:nEx);

% candidate L50 measured from the voltage movie (ExampleL50_movie.mat); the ROIs
% are too sparse to resolve L50, so use the polyline-kymograph value where present,
% falling back to the ROI-based exL50 for any not-yet-measured candidate.
exL50_mov = exL50;
L50movFile = fullfile(save_figto, 'ExampleL50_movie.mat');
if isfile(L50movFile)
    Sv = load(L50movFile);
    for i = 1:min(nEx, numel(Sv.MovieL50))
        if ~isempty(Sv.MovieL50(i).rec) && isfield(Sv.MovieL50(i).rec,'L50_half')
            exL50_mov(i) = Sv.MovieL50(i).rec.L50_half;
        end
    end
end
fprintf('Candidate L50 (movie): %s um  (mean %.0f)\n', ...
    num2str(round(exL50_mov(:)')), mean(exL50_mov,'omitnan'));

figure(142); clf;
tiledlayout(3, 2, 'TileSpacing','compact');

% --- A (top-left): example local candidate ---
axM=nexttile([1 1]);  axP=nexttile([1 1]);
draw_event_footprint(axM, axP, exampleList(1,1), exR0(1), exampleList(1,2), ...
    PCresult, Pixelsize, avgwin_ext, Dmax_ext, exL50_mov(1), ...
    sprintf('Local  (L_{50} = %.0f µm)', exL50_mov(1)));
title('');

% --- B (top-middle): matched broad control ---
axM=nexttile([1 1]);  axP=nexttile([1 1]);
draw_event_footprint(axM, axP, EventTab(broadRow,1), EventTab(broadRow,2), EventTab(broadRow,3), ...
    PCresult, Pixelsize, avgwin_ext, Dmax_ext, EventTab(broadRow,5), ...
    sprintf('Broad  (L_{50} > %.0f µm)', EventTab(broadRow,5)));
title('');

% --- C (bottom-left): normalized amplitude vs contour distance ---
nexttile([1 1]); hold all;
edgesC  = 0:20:Dmax_ext;
binCtrl = binning_data(ctrlProfAll, edgesC);
binCand = binning_data(exProf,      edgesC);
hC = errorbar_shade(binCtrl.centers, binCtrl.mean, binCtrl.sem, [0.5 0.5 0.5]);
hL = errorbar_shade(binCand.centers, binCand.mean, binCand.sem, [1 0 0]);
xlabel('Contour distance from trigger (µm)'); ylabel('Normalized amplitude');
legend([hL hC], {'Local','All'}, 'box','off','location','northeast');
xlim([0 150]); ylim([-0.2 1.05]); box off;

% --- E (top-right): L50 histogram over ALL events (all neurons) ---
nexttile([1 1]); hold all;
histogram(allL50, 0:10:Dmax_ext, 'FaceColor',[0.4 0.4 0.4], 'EdgeColor','none');
xline(mean(exL50_mov,'omitnan'), 'r-', 'linewidth',2, ...
    'Label',sprintf('cand mean %.0f µm', mean(exL50_mov,'omitnan')));
xlabel('L_{50} (µm)'); ylabel('Event count'); box off;
title(sprintf('E  all events (n=%d, %.0f s)', numel(allL50), totalSec));


%title('C  normalized amplitude vs distance');
% 
% % --- D (bottom-right): matched-control L50 histogram; candidates marked ---
% nexttile([1 1]); hold all;
% ctrlL50 = EventTab(allCtrlRows,5);  ctrlL50 = ctrlL50(ctrlL50 >= L50min);   % drop bv noise
% histogram(ctrlL50, 0:20:2*Dmax_ext, 'Normalization','probability', ...
%     'FaceColor',[0.6 0.6 0.6], 'EdgeColor','none');
% yl = ylim;
% xline(median(ctrlL50,'omitnan'), 'k--', 'linewidth',1.5, ...
%     'Label',sprintf('ctrl median %.0f µm', median(ctrlL50,'omitnan')));
% xline(mean(exL50,'omitnan'), 'r-', 'linewidth',2, ...
%     'Label',sprintf('cand mean %.0f µm', mean(exL50,'omitnan')));
% xlabel('L_{50} (µm)'); ylabel('Probability'); box off; ylim(yl);
% title(sprintf('D  L_{50}: %d matched controls vs %d candidates', numel(ctrlL50), nEx));

set_font('Arial'); set_fontsize(12);
set_figsize(190,202);

fprintf('Candidate mean amplitude: %.3f spike-height = %.4f dF/F\n', mean(exAMP), mean(exAMP_dff));
fprintf('Amplitude-decay plot (C) averages %d local-candidate and %d matched-control peaks\n', ...
    sum(~cellfun(@isempty,exProf)), sum(~cellfun(@isempty,ctrlProfAll)));
fprintf('Recording time: %.0f s over %d neurons (L50>=%d um kept)\n', totalSec, numel(foi), L50min);
fprintf('Global events (L50>%d um): %d  ->  %.1f per 10^4 s\n', globalThresh, nGlobal, rate10k);

%% ===== Figure 145 — event footprints sorted by L50 (one neuron) =====
% Montage of event morphologies (colored by normalized window amplitude, trigger
% ROI outlined) from one neuron, ordered small L50 (local) -> large L50 (broad).
% nShow_mon events are sampled evenly across the sorted L50 range; each panel is
% titled with its L50 (red title = one of the nominated local candidates).
fMon      = exampleList(1,1);      % neuron to montage (default: candidate neuron)
nShow_mon = 20;                    % number of footprints to show

rowsMon = find(EventTab(:,1)==fMon & ~isnan(EventTab(:,5)));
[~, ord] = sort(EventTab(rowsMon,5));                 % ascending L50: local -> broad
rowsMon = rowsMon(ord);
pick = unique(round(linspace(1, numel(rowsMon), min(nShow_mon,numel(rowsMon)))));
rowsMon = rowsMon(pick);

SubNorm = interpolateNaN(PCresult{fMon}.Subthreshold / PCresult{fMon}.SpikeHeight);
ftprnts = PCresult{fMon}.Ftprnts(:,:,PCresult{fMon}.sortDist);
candFrames = exampleList(exampleList(:,1)==fMon, 2);  % candidate frames in this neuron

nc = ceil(sqrt(numel(rowsMon)));  nr = ceil(numel(rowsMon)/nc);
figure(145); clf; tiledlayout(nr, nc, 'TileSpacing','none', 'Padding','compact');
for k = 1:numel(rowsMon)
    rr = rowsMon(k);  r0 = EventTab(rr,2);  t0 = EventTab(rr,3);  L50 = EventTab(rr,5);
    W = t0 + (-avgwin_ext:avgwin_ext);  W(W<1 | W>size(SubNorm,2)) = [];
    V = mean(SubNorm(:,W), 2, 'omitnan');  Anorm = V / V(r0);
    nexttile;
    showScaleImageSmooth(ftprnts > 0, Anorm', turbo(256), [0 1]); hold all;
    bd = cell2mat(bwboundaries(imdilate(ftprnts(:,:,r0) > 0, strel('diamond',3))));
    plot(bd(:,2), bd(:,1), 'w', 'linewidth',1);
    isCand = any(abs(candFrames - t0) <= 200);
    tcol = [0 0 0];  if isCand, tcol = [1 0 0]; end
    title(sprintf('%.0f µm', L50), 'FontSize',9, 'Color',tcol);
    axis equal tight off;
end
sgtitle(sprintf('f%d — footprints sorted by L_{50} (local \\rightarrow broad)', fMon));
set_figsize(140*nc, 140*nr);

%% ===== Peri-spike time fraction (+/- 3 ms) — somatic & dendritic spikes =====
% Fraction of the total recording within +/- spikeWin ms of a somatic spike
% (SpikeMat row 1) or a dendritic spike (SpikeClassMat row 4 = dSP), pooled over
% all neurons in foi (1 kHz -> #frames == ms).
spikeWin = 3;
periSoma = 0;  periDsp = 0;  totFrames = 0;
for f = foi
    nT   = size(PCresult{f}.SpikeMat,2);
    soma = find(PCresult{f}.SpikeMat(1,:) > 0);                 % somatic spikes
    if size(PCresult{f}.SpikeClassMat,1) >= 4
        dsp = find(PCresult{f}.SpikeClassMat(4,:) > 0);        % dendritic spikes (dSP)
    else
        dsp = [];
    end
    iS = unique(soma(:) + (-spikeWin:spikeWin));  iS = iS(iS>=1 & iS<=nT);
    iD = unique(dsp(:)  + (-spikeWin:spikeWin));  iD = iD(iD>=1 & iD<=nT);
    periSoma  = periSoma  + numel(iS);
    periDsp   = periDsp   + numel(iD);
    totFrames = totFrames + nT;
end
fprintf('Total recording: %.1f s\n', totFrames/1000);
fprintf('Peri-somatic-spike   (+/-%d ms): %.1f s = %.2f%%\n', spikeWin, periSoma/1000, 100*periSoma/totFrames);
fprintf('Peri-dendritic-spike (+/-%d ms): %.1f s = %.2f%%\n', spikeWin, periDsp/1000,  100*periDsp/totFrames);

%% ===== Measure L50 from the voltage movie along a hand-drawn polyline =====
% ROIs are too sparse to resolve the ~30 um decay, so L50 for the example events
% is measured directly from the movie. Set iMov to the example (row of
% exampleList), run this cell, then draw a polyline along the active dendrite in
% the depolarization image (double-click to finish). The amplitude is sampled
% along the polyline (a polyline kymograph), and L50 = half-decay distance is
% measured from the profile at the peak and appended to ExampleL50_movie.mat.
% NOTE: loads a full ~15000-frame movie segment (like the LocalActivity script) —
% memory-heavy; process one example at a time. cross_half() is at end of file.
iMov        = 5;                            % which example (row of exampleList)
avgwin_mov  = avgwin_ext;                   % +/- ms window to average the event image
baseWin_mov = [1 1000];                     % pre-event baseline window (ms before event)
kymoHalf    = 100;                          % +/- ms shown in the polyline kymograph
dL_mov      = 12;                            % polyLineKymo3 box length along the line (pixels)
dP_mov      = 12;                            % polyLineKymo3 box width across the line (pixels)
saveL50file = fullfile(save_figto, 'ExampleL50_movie.mat');
time_segment = 15000;  overlap = 200;

f = exampleList(iMov,1);  fr = exampleList(iMov,2);

% --- load the movie segment containing the event; motion-residual (as LocalActivity) ---
load(fullfile(fpath{f}, 'PC_Result.mat'), 'Result');
Sdev = load(fullfile(fpath{f}, 'output_data.mat'));   Device_Data = Sdev.Device_Data;
sz = double(Device_Data{1,3}.ROI([2 4]));
frm_end = EndFrame(f);
f_seg = [1:time_segment:frm_end frm_end+1];  f_seg(2:end) = f_seg(2:end)-1;
j   = find(f_seg < fr, 1, 'last');           % segment index containing fr
loc = fr - f_seg(j);                         % local frame index in that segment

mov_mc  = double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'], sz(2), sz(1)));
mov_mc  = mov_mc(:,:,1:end-overlap);
mcseg   = Result.mcTrace(f_seg(j)+1:f_seg(j+1), :);
mov_res = mov_mc - mean(mov_mc,3);
mov_res = SeeResiduals(mov_res, mcseg);
mov_res = SeeResiduals(mov_res, mcseg.^2);
mov_res = SeeResiduals(mov_res, mcseg(:,1).*mcseg(:,end));
nT = size(mov_res,3);

% --- depolarization image: event window (+/-avgwin) minus pre-event baseline ---
Wt = loc + (-avgwin_mov:avgwin_mov);           Wt = Wt(Wt>=1 & Wt<=nT);
Bt = loc - (baseWin_mov(1):baseWin_mov(2));    Bt = Bt(Bt>=1 & Bt<=nT);
depolImg = imgaussfilt(mean(-mov_res(:,:,Wt),3) - mean(-mov_res(:,:,Bt),3), 1);

% --- extract the polyline kymograph with polyLineKymo3 ---
% (left-click to start the line, click vertices, right-click to end)
tK = loc + (-kymoHalf:kymoHalf);   tK = tK(tK>=1 & tK<=nT);   tRel = tK - loc;
dispImg = depolImg;  dispImg(dispImg<0) = 0;                  % display image for the drawing
[kymoTx, kymoROI] = polyLineKymo3(-mov_res(:,:,tK), dL_mov, dP_mov, dispImg);
K = kymoTx';                                                  % [nBox x nTime] = space x time

% arc length from the box centers returned in kymoROI
pos  = cell2mat(cellfun(@(r) mean(r(1:4,:),1), kymoROI(:), 'UniformOutput',false));  % [X Y] per box
s_um = [0; cumsum(sqrt(sum(diff(pos).^2,2)))] * Pixelsize(f);

baseCol = tRel <= -baseWin_mov(1) & tRel >= -baseWin_mov(2);
if ~any(baseCol), baseCol = tRel < 0; end                     % fallback if window out of range
K = K - mean(K(:,baseCol), 2, 'omitnan');                     % per-position baseline

% --- spatial profile at the event & L50 (half-decay distance) ---
evCol = abs(tRel) <= avgwin_mov;
Aprof = movmean(mean(K(:,evCol), 2, 'omitnan'), 3);           % light spatial smoothing
An = Aprof / max(Aprof);
[~, pk] = max(An);
sL = cross_half(s_um, An, pk, -1);   sR = cross_half(s_um, An, pk, +1);
halfL = s_um(pk) - sL;   halfR = sR - s_um(pk);
FWHM  = halfL + halfR;   L50_half = FWHM/2;

% --- record the polyline on the depolarization image ---
figure(150); clf;
imagesc(depolImg); axis image off; colormap(turbo); hold all;
caxis([0 prctile(depolImg(:),99.5)]);  cb = colorbar;  cb.Label.String = '\DeltaF (a.u.)';
plot(pos(:,1), pos(:,2), 'w.-', 'markersize',6);
title(sprintf('ex%d (f%d frame %d): polyline', iMov, f, fr));

% --- show kymograph + measured profile ---
figure(151); clf; tiledlayout(1,2,'TileSpacing','compact');
nexttile; imagesc(tRel, s_um, K); axis xy; colormap(turbo); hold all;
plot([0 0],[s_um(1) s_um(end)],'w--');
xlabel('Time from peak (ms)'); ylabel('Distance along polyline (µm)');
title('Polyline kymograph'); cb=colorbar; cb.Label.String='\DeltaF';
nexttile; plot(An, s_um, 'k-'); hold all; axis ij;
plot([0.5 0.5],[sL sR],'r-'); plot(0.5,sL,'ro'); plot(0.5,sR,'ro'); plot(An(pk),s_um(pk),'b*');
xlabel('Norm. amplitude'); ylabel('Distance along polyline (µm)'); xlim([-0.2 1.1]);
title(sprintf('L_{50}(half)=%.0f µm, FWHM=%.0f µm', L50_half, FWHM));
set_font('Arial'); set_fontsize(13); set_figsize(500,260);

% --- append/update and save ---
rec = struct('f',f, 'frame',fr, 'iMov',iMov, 'polyline',pos, 'arclen_um',s_um, ...
    'profile',Aprof, 'profileNorm',An, 'peak_um',s_um(pk), ...
    'L50_half',L50_half, 'FWHM',FWHM, 'halfLeft',halfL, 'halfRight',halfR, ...
    'kymo',K, 'tRel',tRel, 'dL',dL_mov, 'dP',dP_mov);
rec.kymoROI = kymoROI;                        % box polygons (cell) stored separately
if isfile(saveL50file), Sv = load(saveL50file);  MovieL50 = Sv.MovieL50;  else, MovieL50 = struct([]); end
MovieL50(iMov).rec = rec;
save(saveL50file, 'MovieL50');
fprintf('ex%d (f%d frame %d): L50(half)=%.1f um, FWHM=%.1f um  (halfL=%.1f, halfR=%.1f)  saved -> %s\n', ...
    iMov, f, fr, L50_half, FWHM, halfL, halfR, saveL50file);

clear mov_mc mov_res                          % free movie memory

%% Peak decay over distance (figure 2)
dClass2plot=[[{1},{1}];[{[3 5]} {[3 5]}];[{[4]} {[4]}];[{[3 5]} {[1 2 3 4 5]}]; [{[1 2 3 4 5]},{[1 2 3 4 5]}]];
%dClass2plot=[[{[1 2 3 4 5]},{[1 2 3 4 5]}]]; % Trigger, Read
pair_str=[];
cmap_corrscatter=[0.5 0.8 1; 1 0.5 0.5];
cmap_neuron=hsv(max(foi));
M=[]; S=[]; M1d=[]; S1d=[]; dendbin=[0 15:45:450 500 600];

saveto='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism';
load([saveto '/PeakTriggeredAverageResult_20260620.mat']);

% PeakDecayData_Pair=[]; PeakDecayData_Cat=[];
% BadROIs={[26 3 19 20 24 25 26 36],[24 23 24 22 26 25 29]};
% for f=foi
%     DistMat=PCresult{f}.interDendDist(PCresult{f}.sortDist,PCresult{f}.sortDist);
%     DistMat=max(cat(3,DistMat,DistMat'),[],3);
%     DistMat1d = abs(PCresult{f}.dendaxis1d - PCresult{f}.dendaxis1d');
%     PTA_tmp=PeakSTA{f};
%     PTA_tmp=PTA_tmp-prctile(PTA_tmp(:,nTauPeak(1)+[-800:800],:),30,2);
%     PeakZeroMat=squeeze(PTA_tmp(:,nTauPeak(1)+1,:)); %row: Read, column: Trigger
%     PeakZeroMat_ch=cellfun(@(x) squeeze(x(:,nTauPeak(1)+1,:)),PeakSTA_ch(f,:),'UniformOutput',false); %row: Read, column: Trigger
%     ft=PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist);
%     ftcoord=get_coord(ft);
%     BadFs=(cellfun(@(x) x(1),BadROIs));
%     if any(ismember(BadFs,f))
%     PCresult{f}.labelvec(BadROIs{find(ismember(BadFs,f))}(2:end))=NaN;
%     end
% 
%     for dc=1:size(dClass2plot,1) %region class to plot
% 
%         TrigROIs=find(ismember(PCresult{f}.labelvec,dClass2plot{dc,1}));
%         ReadROIs=find(ismember(PCresult{f}.labelvec,dClass2plot{dc,2}));
% 
%         for tr=TrigROIs
%             dax=DistMat(ReadROIs,tr);
%             dax1d=DistMat1d(ReadROIs,tr);
%             Volt=PeakZeroMat(ReadROIs,tr);
%             Voltodd=PeakZeroMat_ch{1}(ReadROIs,tr);
%             Volteven=PeakZeroMat_ch{2}(ReadROIs,tr);
%             [~, trInd_inReadROI]=min(dax);            
%             PeakDecayData_Cat=[PeakDecayData_Cat; {f, dClass2plot{dc,1}, dClass2plot{dc,2}, tr, ReadROIs, dax, dax1d, Volt, ...
%                                                    ftcoord(tr,:)}, Voltodd, Volteven];
%         end
% 
%         %uppertri_ind=triu(ones(size(DistMat)))';
%         label_ind=ismember(PCresult{f}.labelvec,dClass2plot{dc,2})'*ismember(PCresult{f}.labelvec,dClass2plot{dc,1})>0;
%         daxpair=tovec(DistMat(label_ind));
%         dax1dpair=tovec(DistMat1d(label_ind));
%         peakdat2fit=tovec(PeakZeroMat(label_ind));
%         peakVoltage=mean(peakdat2fit(daxpair<20),'omitnan');
%         N_points=sum(~isnan(peakdat2fit));
%         prc_triggeringpoint=sum(peakdat2fit<peakVoltage)./N_points;
% 
%         if N_points>5
%             [fitted_exp params R2]=expfit_wBd(daxpair,peakdat2fit,[0:500],[peakVoltage -500],[peakVoltage*0.9 -485],...
%                                                                            [peakVoltage*1.2 -50],1./(daxpair+1).^.2);
%             [fitted_exp1d params1d R21d]=expfit_wBd(dax1dpair,peakdat2fit,[0:500],[peakVoltage -500],[peakVoltage*0.7 -max(dax1dpair)*2],[peakVoltage*1.5 -50]);
%             PeakDecayData_Pair=[PeakDecayData_Pair; {f, dClass2plot{dc,1}, dClass2plot{dc,2}, TrigROIs, ReadROIs, -params(2), R2, -params1d(2), R21d, N_points, ...
%                 prc_triggeringpoint,daxpair, dax1dpair,peakdat2fit, fitted_exp, fitted_exp1d}];
%         end
%     end
% end
% 
% fieldName={'NeuronID','TriggerClass','ReadClass','TriggerROI','ReadROIs','dax','dax1d','TriggerVoltage','Trigger_coord','TriggerVoltageodd','TriggerVoltageeven'};
% PeakDecayData_Cat=array2table(PeakDecayData_Cat,'VariableNames',fieldName);
% 
% fieldName={'NeuronID','TriggerClass','ReadClass','TriggerROI','ReadROIs','LengthConstant','Rsquare','LengthConstant1d','Rsquare1d','NFitpoints','Percentile','dax','dax1d','TriggerVoltage','Expfit','Expfit1d'};
% PeakDecayData_Pair=array2table(PeakDecayData_Pair,'VariableNames',fieldName);
% 
% CorrDistMat=[];
% ClassInd=cellfun(@(x,y) isequal(x,dClass2plot{5,1}) & isequal(y,dClass2plot{5,2}),PeakDecayData_Cat.TriggerClass,PeakDecayData_Cat.ReadClass,'UniformOutput',0);
% for f=foi
%     CorrDistMat{f}=[];
%     showInd=find(cell2mat(PeakDecayData_Cat.NeuronID)==f & cell2mat(ClassInd))';
%     DistMat=PCresult{f}.interDendDist(PCresult{f}.sortDist,PCresult{f}.sortDist);
%     DistMat=max(cat(3,DistMat,DistMat'),[],3);
%     DistMat1d = abs(PCresult{f}.dendaxis1d - PCresult{f}.dendaxis1d');
% 
%     for j=1:length(showInd)
%         for jj=j:length(showInd)
%             if isequal(PeakDecayData_Cat.TriggerClass{showInd(j)},PeakDecayData_Cat.TriggerClass{showInd(jj)}) &...
%                isequal(PeakDecayData_Cat.ReadClass{showInd(j)},PeakDecayData_Cat.ReadClass{showInd(jj)}) % If read and trigger class is same for both paris,
%                 if j==jj
%                     CorrDistMat{f}=[CorrDistMat{f}; [DistMat(PeakDecayData_Cat.TriggerROI{showInd(j)},PeakDecayData_Cat.TriggerROI{showInd(jj)}) ...
%                         DistMat1d(PeakDecayData_Cat.TriggerROI{showInd(j)},PeakDecayData_Cat.TriggerROI{showInd(jj)}) ...
%                         corr(PeakDecayData_Cat.TriggerVoltageeven{showInd(j)},PeakDecayData_Cat.TriggerVoltageodd{showInd(jj)})...
%                         PeakDecayData_Cat.TriggerROI{showInd(j)} PeakDecayData_Cat.TriggerROI{showInd(jj)}]];
%                 else
%                     CorrDistMat{f}=[CorrDistMat{f}; [DistMat(PeakDecayData_Cat.TriggerROI{showInd(j)},PeakDecayData_Cat.TriggerROI{showInd(jj)}) ...
%                         DistMat1d(PeakDecayData_Cat.TriggerROI{showInd(j)},PeakDecayData_Cat.TriggerROI{showInd(jj)}) ...
%                         corr(PeakDecayData_Cat.TriggerVoltage{showInd(j)},PeakDecayData_Cat.TriggerVoltage{showInd(jj)})...
%                         PeakDecayData_Cat.TriggerROI{showInd(j)} PeakDecayData_Cat.TriggerROI{showInd(jj)}]];
%                 end
%             end
%         end
%     end
%     if ~isempty(CorrDistMat{f})
%         fieldName={'ContourDist','LaminarDist','Correlation','Pair_i','Pair_j'};
%         CorrDistMat{f}=array2table(CorrDistMat{f},'VariableNames',fieldName);
%     end
% end

% Do interactive exp fitting
for f=[11]
    dc=4; %trigger: Apical, Read: all
    fprintf('Interactive fitting f # %d\n', f);
    ClassInd=cellfun(@(x,y) isequal(x,dClass2plot{dc,1}) & isequal(y,dClass2plot{dc,2}),PeakDecayData_Pair.TriggerClass,PeakDecayData_Pair.ReadClass,'UniformOutput',0);
    PTA2analyze=find(cell2mat(PeakDecayData_Pair.NeuronID)==f & cell2mat(ClassInd));
    opt = struct();
    opt.x_fit   = 0:ceil(max(PeakDecayData_Pair.dax{PTA2analyze})); opt.t_guess = 100;
    opt.tau_lb  = 0;  opt.tau_ub  = 1000;
    fitResult=interactive_expfit(PeakDecayData_Pair.dax{PTA2analyze},PeakDecayData_Pair.TriggerVoltage{PTA2analyze},opt);
    if fitResult.aborted
        break;
    end
    PeakDecayData_Pair.fitresult{PTA2analyze}=fitResult;
    PeakDecayData_Pair.Expfit{PTA2analyze}=fitResult.yFit;
    PeakDecayData_Pair.Rsquare{PTA2analyze}=fitResult.Rsq;
    PeakDecayData_Pair.LengthConstant{PTA2analyze}=fitResult.p(2);
end
saveto='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism';
save([saveto '/PeakTriggeredAverageResult_20260620.mat'],'PeakDecayData_Cat','PeakDecayData_Pair','dClass2plot','-v7.3');
%% Show example of contour decay and PTA correlation 
f2show=23;
ClassInd=cellfun(@(x,y) isequal(x,dClass2plot{5,1}) & isequal(y,dClass2plot{5,2}),PeakDecayData_Cat.TriggerClass,PeakDecayData_Cat.ReadClass,'UniformOutput',0);
showInd=find(cell2mat(PeakDecayData_Cat.NeuronID)==f2show & cell2mat(ClassInd))';
Ftcoord=get_coord(PCresult{f2show}.Ftprnts(:,:,PCresult{f2show}.sortDist));
skelDend=point2img(PCresult{f2show}.swc,2,size(PCresult{f2show}.ref_im));
idx=SelectFromScatter(CorrDistMat{f2show}.LaminarDist,CorrDistMat{f2show}.ContourDist);
pth=[]; g=1;
for ss=idx'
    [~, pp]= geodesic_distance(skelDend,Ftcoord(CorrDistMat{f2show}.Pair_i(ss),:),Ftcoord(CorrDistMat{f2show}.Pair_j(ss),:));
    bw_skel = bwskel(pp);       % or bwmorph(bw, 'skel', Inf)
    [r, c]  = find(bw_skel);
    pth{g}=[r c];
    pth{g}(1,:)=Ftcoord(CorrDistMat{f2show}.Pair_i(ss),[2 1]);
    pth{g}(end+1,:)=Ftcoord(CorrDistMat{f2show}.Pair_j(ss),[2 1]);
g=g+1;
end
figure(131); clf; tiledlayout(2,12);
ax2=nexttile(1,[1 4]);
dendaxis1d=PCresult{f2show}.dendaxis1d;
PTA2show=PeakDecayData_Cat.TriggerVoltage(showInd);
ReadROI2show=PeakDecayData_Cat.ReadROIs(showInd);
distcmap=vec2cmap(dendaxis1d(cell2mat(PeakDecayData_Cat.TriggerROI(showInd))),'turbo');
drange1d=round([min(dendaxis1d(cell2mat(PeakDecayData_Cat.TriggerROI(showInd)))) max(dendaxis1d(cell2mat(PeakDecayData_Cat.TriggerROI(showInd))))]);
[~, dsort]=cellfun(@(x) sort(dendaxis1d(x)),ReadROI2show,'UniformOutput',false);
hold all;
l=cellfun(@(x,y,z) plot(dendaxis1d(y(z)),x(z),'linewidth',1),PTA2show,ReadROI2show,dsort,'UniformOutput',false);
cellfun(@(x,y) set(x,'Color',y),l,num2cell(distcmap,2));
xlabel('Somatodendritic coordinate (µm)'); ylabel('∆F/F');
title('Peak-triggered average (PTA) voltage');
cb=colorbar; colormap(ax2,turbo); 
cb.Label.String=sprintf('Somatodendritic location \n of triggering ROI (µm)');
cb.Ticks=cb.Limits; cb.TickLabels=[drange1d];
box off; axis tight;

ax1=nexttile([1 4]);
scatter(CorrDistMat{f2show}.LaminarDist,CorrDistMat{f2show}.ContourDist,10,vec2cmap(CorrDistMat{f2show}.Correlation,'turbo',[-1 1]),'filled'); hold all
scatter(CorrDistMat{f2show}.LaminarDist(idx),CorrDistMat{f2show}.ContourDist(idx),60,hsv(length(idx)),'linewidth',2); hold all
xlabel('Somatodendritic distance (µm)'); ylabel('Contour distance (µm)'); 
title('PTA correlation');
cb=colorbar; colormap(ax1,turbo); 
cb.Label.String='Correlation';
cb.Ticks=cb.Limits; cb.TickLabels=[-1 1];

nexttile([1 4]);
cmap_line=hsv(length(pth));
%imshow2(PCresult{f2show}.ref_im,[600 3500]); hold all
scatter(PCresult{f2show}.swc(:,1),PCresult{f2show}.swc(:,2),7,[0.6 0.6 0.6],'filled'); hold all;
for g=1:length(pth)
scatter(pth{g}(:,2),pth{g}(:,1),5,cmap_line(g,:),'filled');
scatter(pth{g}([1 end],2),pth{g}([1 end],1),70,cmap_line(g,:),'linewidth',2);
end
quiver(Ftcoord(dendaxis1d==0,1), Ftcoord(dendaxis1d==0,2), PCresult{f2show}.principal_axis(1)*300, 300*PCresult{f2show}.principal_axis(2), 0, 'Color', 'k', 'LineWidth', 2, 'MaxHeadSize', 0.5);
drawScaleBar(100/Pixelsize(f2show),'Horizontal','color','k','position',[100 160]);
axis equal tight off;
title('Example contour paths between ROI pairs')

nexttile([1 3]);
for j=showInd
    d2show=PeakDecayData_Cat.dax{j};
    v2show=PeakDecayData_Cat.TriggerVoltage{j};
    [d2show, dsort]=sort(d2show);
    v2show=v2show(dsort);
    %v2show=v2show./v2show(1);
    plot(d2show,v2show,'color',[0.5 0.5 0.5]); hold all;
end
v2showCat=cellfun(@(x,y) [x y],PeakDecayData_Cat.dax(showInd),PeakDecayData_Cat.TriggerVoltage(showInd),'UniformOutput',false)';
v2showCat_binned=binning_data(v2showCat,[0:30:450]);
errorbar_shade(v2showCat_binned.centers,v2showCat_binned.mean,v2showCat_binned.std,[1 0 0]);
xlabel('Contour distance from triggering ROI (µm)'); ylabel('Peak-triggered ∆F/F');
xlim([0 400]); box off;

nexttile([1 3]);
for j=showInd
    d2show=PeakDecayData_Cat.dax1d{j};
    v2show=PeakDecayData_Cat.TriggerVoltage{j};
    [d2show, dsort]=sort(d2show);
    v2show=v2show(dsort);
    %v2show=v2show./v2show(1);
    plot(d2show,v2show,'color',[0.5 0.5 0.5]); hold all;
end
v2showCat=cellfun(@(x,y) [x y],PeakDecayData_Cat.dax1d(showInd),PeakDecayData_Cat.TriggerVoltage(showInd),'UniformOutput',false)';
v2showCat_binned=binning_data(v2showCat,[0:30:450]);
errorbar_shade(v2showCat_binned.centers,v2showCat_binned.mean,v2showCat_binned.std,[1 0 0]);
xlabel('Somatodendritic distance from triggering ROI (µm)'); ylabel('Peak-triggered ∆F/F');
xlim([0 400]); box off;

nexttile([1 3]);
scatter_density(CorrDistMat{f2show}.ContourDist,CorrDistMat{f2show}.Correlation,40);
xlabel('Contour distance (µm)'); ylabel('PTA Correlation');
xlim([0 600]); box off; delete(colorbar);

nexttile([1 3]);
scatter_density(CorrDistMat{f2show}.LaminarDist,CorrDistMat{f2show}.Correlation,40);
xlabel('Somatodendritic distance (µm)'); ylabel('PTA Correlation');
xlim([0 500]); box off; delete(colorbar);
set_font('Arial'); set_fontsize(16);
set_figsize(450,200);

figure(132); clf;
shortbd_L=[1 100]; longbd_L=[100 200]; shortbd_C=[1 100]; longbd_C=[100 200];
Ratiothreshold=1.7;
% ShortCLongL_ind=cellfun(@(x) x.ContourDist > shortbd_C(1) & x.ContourDist < shortbd_C(2) & x.LaminarDist > longbd_L(1) & x.LaminarDist < longbd_L(2),CorrDistMat(foi),'UniformOutput',false);
% LongCShortL_ind=cellfun(@(x) x.ContourDist > longbd_C(1) & x.ContourDist < longbd_C(2) & x.LaminarDist > shortbd_L(1) & x.LaminarDist < shortbd_L(2),CorrDistMat(foi),'UniformOutput',false);

ShortCLongL_ind=cellfun(@(x) x.ContourDist./x.LaminarDist<Ratiothreshold & x.LaminarDist<200 & x.LaminarDist>0,CorrDistMat(foi),'UniformOutput',false);
LongCShortL_ind=cellfun(@(x) x.ContourDist./x.LaminarDist>Ratiothreshold & x.LaminarDist<200 & x.LaminarDist>0,CorrDistMat(foi),'UniformOutput',false);

Corr_SCLL=cellfun(@(x,y) mean(x.Correlation(y),'omitnan'),CorrDistMat(foi),ShortCLongL_ind,'UniformOutput',false);
Corr_LCSL=cellfun(@(x,y) mean(x.Correlation(y),'omitnan'),CorrDistMat(foi),LongCShortL_ind,'UniformOutput',false);

p=Boxplot_wPoints2(cell2mat([Corr_SCLL' Corr_LCSL']),hsv(2));
set(gca,'xtick',[1 2],'xticklabel',{['Long contour ,short laminar'], 'Short contour, long laminar'});
drawPValueLines(p,0.05,'TextYOffset',0.05);
ylabel('Mean correlation'); box off;
set_font('Arial'); set_fontsize(16); 
set_figsize(130,230)

%% Plot local stimulation data, figure 2J, K
load(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism/LocalStimulationResult_20260620.mat']);
t2show=[-5:5:60];
cmap_errorbar=[0.3 0.3 0.3; 0.5 0.5 0.5; 0.7 0.7 0.7; 0 0.6 1];
tInd=cell2mat(Lstaall.ReadBlueTime)>=20 & cell2mat(Lstaall.ReadBlueTime)<=40;
path2show=find(cell2mat(Lstaall.StimROI_percentile)>0.9 & cell2mat(Lstaall.Fit_R2)>0.2 & tInd & cell2mat(Lstaall.Fit_length)>20 & cell2mat(Lstaall.Fit_length)<300); % Percentile, R^2
Stim2show=unique(cell2mat(Lstaall.StimInd(path2show)));
ClassInd=cellfun(@(x,y) isequal(x,dClass2plot{4,1}) & isequal(y,dClass2plot{4,2}),PeakDecayData_Pair.TriggerClass,PeakDecayData_Pair.ReadClass,'UniformOutput',0);
Lc2Spont=cell2mat(PeakDecayData_Pair.LengthConstant(cell2mat(ClassInd)));

DecayLengthCat=[]; DecayCurveCat=[];
for i=1:length(Stim2show)
    Dat2showInd=find(ismember(cell2mat(Lstaall.StimInd),Stim2show(i)));
    [dax, dd]=cellfun(@(x) sort(x(:,1)), Lstaall.dendaxis(Dat2showInd),'UniformOutput',0);
    SortedVoltageSTA=cellfun(@(x,y) x(y), Lstaall.VoltageSTA(Dat2showInd),dd,'UniformOutput',0);
    t2read=cell2mat(cellfun(@(x) sort(x(:,1)), Lstaall.ReadBlueTime(Dat2showInd),'UniformOutput',0));
    [t2read, tsort]=sort(t2read);
    KymoSta=cell2mat(SortedVoltageSTA(tsort)');
    Lc2print=[]; g=1;
    ind2print=ismember(Dat2showInd,find(ismember(cell2mat(Lstaall.ReadBlueTime),[15:35]))); %StimInd + Time
    Dat2showInd_sub=Dat2showInd(ind2print);
    [maxR2, mi]=max(cell2mat(Lstaall.Fit_R2(Dat2showInd_sub)));
    DecayLengthCat= [DecayLengthCat; [Stim2show(i) Lstaall.Fit_length{Dat2showInd_sub(mi)} maxR2]];
    DecayCurveCat = [DecayCurveCat {[Lstaall.dendaxis{Dat2showInd_sub(mi)}, Lstaall.VoltageSTA{Dat2showInd_sub(mi)}]}];
end

Lc2Stim=DecayLengthCat(:,2);
[~, dsort]=cellfun(@(x) sort(x(:,1)),DecayCurveCat,'UniformOutput',0);
data2show_ave2=cellfun(@(x,y) [x(y,1) x(y,2)/mean(x(x(:,1)<10,2))],DecayCurveCat,dsort,'UniformOutput',false); %normalize
LocalstimVoltage_bin=binning_data(data2show_ave2,dendbin);

figure(158); clf; %check PTAs
ClassInd=cellfun(@(x,y) isequal(x,dClass2plot{4,1}) & isequal(y,dClass2plot{4,2}),PeakDecayData_Pair.TriggerClass,PeakDecayData_Pair.ReadClass,'UniformOutput',0);
ClassInd=find(cell2mat(ClassInd));
for c=ClassInd'
    nexttile([1 1]);
    scatter(PeakDecayData_Pair.dax{c},PeakDecayData_Pair.TriggerVoltage{c},10,'filled'); hold all;
    plot(PeakDecayData_Pair.Expfit{c},'r');
    title(PeakDecayData_Pair.LengthConstant{c});
end

StimN=unique(cell2mat(Lstaall.StimInd(Stim2show))); 
MNlist=[Mouse(cell2mat(Lstaall.fileInd(Stim2show))) NeuronInd(cell2mat(Lstaall.fileInd(Stim2show)))];
Mlist=unique(Mouse(cell2mat(Lstaall.fileInd(Stim2show))));
[UniqueMNlist]=unique(MNlist);
fprintf('Stim-triggered average, # of stimulation: %2.0f, # of Neuron: %2.0f, # of mouse: %2.0f \n',length(StimN),size(UniqueMNlist,1),size(Mlist,1));
%%
figure(156); clf; tiledlayout(2,2);
nexttile([1 1]);
PTA2show=[]; PTA2show_all=[];
for d=4%size(dClass2plot,1)
    ClassInd=cellfun(@(x,y) isequal(x,dClass2plot{d,1}) & isequal(y,dClass2plot{d,2}),PeakDecayData_Pair.TriggerClass,PeakDecayData_Pair.ReadClass,'UniformOutput',0);  
    for f=foi
        idx=cell2mat(PeakDecayData_Pair.NeuronID)==f & cell2mat(ClassInd);
        voltdat_tmp=cellfun(@(x,y,z) [x(z.idxSel) y(z.idxSel)],PeakDecayData_Pair.dax(idx), PeakDecayData_Pair.TriggerVoltage(idx),PeakDecayData_Pair.fitresult(idx),'UniformOutput',0)';
        %voltdat_tmp=cellfun(@(x,y,z) [x y],PeakDecayData_Pair.dax(idx), PeakDecayData_Pair.TriggerVoltage(idx),'UniformOutput',0)';
        voltdat_tmp{1}(:,2)=voltdat_tmp{1}(:,2)./mean(voltdat_tmp{1}(voltdat_tmp{1}(:,1)<5,2));
        voltdat_bin_tmp=binning_data(voltdat_tmp,dendbin);
        PTA2show{f}=[voltdat_bin_tmp.centers' voltdat_bin_tmp.mean'];
        %PTA2show{f}=voltdat_tmp{1};
    end
    PTA2show_all{d}=binning_data(PTA2show(foi),dendbin);
    errorbar_shade(PTA2show_all{d}.centers,PTA2show_all{d}.mean,PTA2show_all{d}.sem,cmap_errorbar(2,:)); hold all;
end
l=errorbar_shade(LocalstimVoltage_bin.centers,LocalstimVoltage_bin.mean,LocalstimVoltage_bin.sem,cmap_errorbar(4,:));
hold all;
meanfitPTA=@(x) exp(-x/mean(Lc2Spont));
meanfitSTA=@(x) exp(-x/mean(Lc2Stim));
% plot([0:300],meanfitSTA([0:300]),'color',cmap_errorbar(4,:),'linestyle','--','linewidth',1.5);
% plot([0:300],meanfitPTA([0:300]),'color',[0 0 0],'linestyle','--','linewidth',1.5);
xlim([0 250]); box off;
ylabel('Normalized voltage'); xlabel('Distance from triggering ROI (µm)');

nexttile([1 1])
cmap=[0 0.6 1;0 0 0; 0.7 0 0.9];
LengthC2show={Lc2Stim, Lc2Spont};
p=Violin_wPoints(LengthC2show,cmap([1 2],:)); axis tight; hold all
line([1 2],[270 270],'LineWidth',2,'color','k');
text([1.5],[280],'N.S.','HorizontalAlignment','center','FontSize',18);
set(gca,'xtick',[1 2],'xticklabel',{'Opto. Stim.','Spontaneous'})
ylabel('Decay length (µm)'); ylim([0 300]); xlim([0.5 2.5]); box off;

nexttile([1 1])
l=errorbar_shade(LocalstimVoltage_bin.centers,LocalstimVoltage_bin.mean,LocalstimVoltage_bin.sem,cmap_errorbar(4,:));
hold all;
meanfitSTA=@(x) exp(-x/mean(Lc2Stim));
plot([0:300],meanfitSTA([0:300]),'color',cmap_errorbar(4,:),'linestyle','--','linewidth',1.5);
ylabel('Normalized voltage'); xlabel('Distance from stimulation ROI (µm)');
title(sprintf('\\lambda : %2.0f ± %2.0f µm',mean(Lc2Stim),std(Lc2Stim)));
xlim([0 200]); box off;

set_font('Arial'); set_fontsize(20);
set_figsize(240,270);
%% Figure 2E, Trigger, vertex, endpoint
R_th=1.5;
Laminar_Maxdist=50;
PairDistVertex=[];
for f=foi
    DistMat=PCresult{f}.interDendDist(PCresult{f}.sortDist,PCresult{f}.sortDist);
    DistMat=max(cat(3,DistMat,DistMat'),[],3);
    DistMat1d = abs(PCresult{f}.dendaxis1d - PCresult{f}.dendaxis1d');
    uppertri_ind=triu(ones(size(DistMat)));
    RMat=DistMat./DistMat1d;
    [pairr pairc]=find(DistMat1d<Laminar_Maxdist & RMat>R_th & uppertri_ind>0);
    swcmask=point2img(PCresult{f}.swc,2,size(PCresult{f}.ref_im));
    ftprnts=PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist);
    ftprnts_vec=tovec(ftprnts>0);
    ftcoord=get_coord(ftprnts);
    labelvec=PCresult{f}.labelvec;
    for p=1:length(pairr)
        [~, pth]=geodesic_distance(swcmask,ftcoord(pairr(p),:),ftcoord(pairc(p),:));
        RoiOnpth=find(sum(ftprnts_vec(tovec(pth>0),:),1)>0);
        if PCresult{f}.dendaxis1d(pairr(p))<=0 % basal
            [~, maxind]=max(PCresult{f}.dendaxis1d(RoiOnpth));
        else %apical
            [~, maxind]=min(PCresult{f}.dendaxis1d(RoiOnpth));            
        end
        PairDistVertex=[PairDistVertex; [f, pairr(p), pairc(p), RoiOnpth(maxind), DistMat1d(pairr(p),RoiOnpth(maxind)), DistMat1d(pairc(p),RoiOnpth(maxind)), ...
            DistMat(pairr(p),RoiOnpth(maxind)),DistMat(pairc(p),RoiOnpth(maxind)), RMat(pairr(p),pairc(p)), length(RoiOnpth), labelvec(pairr(p)), ...
            labelvec(pairc(p)), labelvec(RoiOnpth(maxind)), DistMat1d(pairr(p),pairc(p))]];
    end
end

fieldName={'FileInd','StartROI','EndROI','VertexROI','Dist1D_SV','Dist1D_VE','Dist_SV','Dist_VE','Cont2LamRatio','NROIonPth','LabelS','LabelE','LabelV','Dist1D_SE'};
PairDistVertex=array2table(PairDistVertex,'VariableNames',fieldName);

Dist1D_thresSVE=50;
Dist1D_thresSE=30;
Vertex2show_cat=[];
CorrPath2show=[]; PeakPath2show=[]; DistPath2show=[]; Dist1DPath2show=[];
for f=foi
    CorrPath2show{f}=[]; PeakPath2show{f}=[]; DistPath2show{f}=[]; Dist1DPath2show{f}=[];
    vert2show=find(PairDistVertex.FileInd==f & PairDistVertex.NROIonPth>2 & ...
        PairDistVertex.Dist1D_SV>Dist1D_thresSVE & PairDistVertex.Dist1D_VE>Dist1D_thresSVE & PairDistVertex.Dist1D_SE<Dist1D_thresSE);
    PTA_tmp=PeakSTA{f};
    PTA_tmp=PTA_tmp-prctile(PTA_tmp(:,nTauPeak(1)+[-800:800],:),30,2);
    PeakZeroMat=squeeze(PTA_tmp(:,nTauPeak(1)+1,:))./PCresult{f}.SpikeHeight; %row: Read, column: Trigger
    Corrtmp=corrMat{f,1};

    swcmask=point2img(PCresult{f}.swc,2,size(PCresult{f}.ref_im));
    ftprnts=PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist);
    ftcoord=get_coord(ftprnts);
    ROIlist=[PairDistVertex.StartROI(vert2show) PairDistVertex.EndROI(vert2show) PairDistVertex.VertexROI(vert2show)];
    usedVert=zeros(length(vert2show),1);
    for v=1:length(vert2show)
        PeakPath=[PeakZeroMat(ROIlist(v,1),ROIlist(v,1)) PeakZeroMat(ROIlist(v,1),ROIlist(v,3)) PeakZeroMat(ROIlist(v,1),ROIlist(v,2))];
        if max(PeakPath)<0.35
        CorrPath=[Corrtmp(ROIlist(v,1),ROIlist(v,1)) Corrtmp(ROIlist(v,1),ROIlist(v,3)) Corrtmp(ROIlist(v,1),ROIlist(v,2))];
        DistPath=[PairDistVertex.Dist_SV(vert2show(v)) PairDistVertex.Dist_VE(vert2show(v))];
        Dist1dPath=[PairDistVertex.Dist1D_SV(vert2show(v)) PairDistVertex.Dist1D_VE(vert2show(v))];
        Vertex2show_cat=[Vertex2show_cat; {DistPath, Dist1dPath, CorrPath, PeakPath, f, vert2show(v)}];
        usedVert(v)=1;
        end
    end

    DistPath2show{f}=[-PairDistVertex.Dist_SV(vert2show(usedVert>0)) zeros(sum(usedVert),1) PairDistVertex.Dist_VE(vert2show(usedVert>0))];
    Dist1DPath2show{f}=[-PairDistVertex.Dist1D_SV(vert2show(usedVert>0)) zeros(sum(usedVert),1) PairDistVertex.Dist1D_VE(vert2show(usedVert>0))];
    PeakPath2show{f}=cell2mat(Vertex2show_cat(cell2mat(Vertex2show_cat(:,5))==f,4));
    CorrPath2show{f}=cell2mat(Vertex2show_cat(cell2mat(Vertex2show_cat(:,5))==f,3));
end
fieldName={'DistPath','Dist1DPath','CorrPath','PeakPath','FileInd','VertexInd'};
Vertex2show_cat=array2table(Vertex2show_cat,'VariableNames',fieldName);
%% Figure 2E, Trigger, vertex, endpoint
figure(120); clf; dendbin=[-300 -10 10 300];
tiledlayout(1,2,'TileSpacing','tight');
nexttile([1 1]); hold all;
ValidF=cell2mat(cellfun(@(x) size(x,1),DistPath2show,'UniformOutput',0))>0;
dat2bin=cellfun(@(x,y) [tovec(x) tovec(y)],DistPath2show(ValidF),PeakPath2show(ValidF),'UniformOutput',0);
binResult=binning_data(dat2bin,dendbin);
cellfun(@(x,y) plot(x',y','color',[0.6 0.6 0.6]),DistPath2show(ValidF),PeakPath2show(ValidF),'UniformOutput',0);
errorbar(binResult.centers,binResult.mean,binResult.sem,'color',[1 0 0],'linewidth',2)
p=get_pValue(cell2mat(PeakPath2show(ValidF)'),1,'Method','ttest');
ylim([-0.05 0.35]);
drawPValueLines(p,0,'XCoord',binResult.centers,'TextYOffset',0.015,'StepHeight',0.03);
ylim([-0.05 0.45]);
xlabel(['Contour distance' newline 'from vertex (µm)']);
ylabel(['Peak-triggered voltage' newline '(spike height)']);

nexttile([1 1]); hold all;
ValidF=cell2mat(cellfun(@(x) size(x,1),DistPath2show,'UniformOutput',0))>0;
dat2bin=cellfun(@(x,y) [tovec(x) tovec(y)],DistPath2show(ValidF),CorrPath2show(ValidF),'UniformOutput',0);
dat2bin2D=cellfun(@(x,y,z) [tovec(x) abs(tovec(y)) tovec(z)],DistPath2show(ValidF),Dist1DPath2show(ValidF),CorrPath2show(ValidF),'UniformOutput',0);
dat2bin2D=cell2mat(dat2bin2D');
binResult=binning_data(dat2bin,dendbin);
bin2DResult=binning_data2D(dat2bin2D(:,1),dat2bin2D(:,2),dat2bin2D(:,3),dendbin,dendbin);

cellfun(@(x,y) plot(x',y','color',[0.6 0.6 0.6]),DistPath2show(ValidF),CorrPath2show(ValidF),'UniformOutput',0);
errorbar(binResult.centers,binResult.mean,binResult.sem,'color',[1 0 0],'linewidth',2)
p=get_pValue(cell2mat(CorrPath2show(ValidF)'),1,'Method','ttest');
ylim([0 1]);
drawPValueLines(p,0,'XCoord',binResult.centers,'TextYOffset',0.035,'StepHeight',0.07);
xlabel(['Contour distance' newline 'from vertex (µm)']);
ylabel('Correlation'); ylim([0 1.2]);
set_font('Arial'); set_fontsize(18);
set_figsize(200,140)

[fileTriggerList]=unique([PairDistVertex.FileInd(cell2mat(Vertex2show_cat.VertexInd)) PairDistVertex.StartROI(cell2mat(Vertex2show_cat.VertexInd))],'row');
fprintf('n = %2.0f trigger ROI, %2.0f contours, %2.0f neurons, %2.0f mice \n',...
        size(fileTriggerList,1),size(Vertex2show_cat,1),...
        length(unique(cell2mat(Vertex2show_cat.FileInd))),length(unique(mouse(unique(cell2mat(Vertex2show_cat.FileInd))))))

%% Figure S?, show contour/laminar example
f2show=27;
vert2show=find(PairDistVertex.FileInd==f2show & PairDistVertex.NROIonPth>2 & ...
    PairDistVertex.Dist1D_SV>Dist1D_thres & PairDistVertex.Dist1D_VE>Dist1D_thres & PairDistVertex.Dist1D_SE<Dist1D_thresSE);
swcmask=point2img(PCresult{f2show}.swc,2,size(PCresult{f2show}.ref_im));
ftprnts=PCresult{f2show}.Ftprnts(:,:,PCresult{f2show}.sortDist);
ftcoord=get_coord(ftprnts);
ROIlist=[PairDistVertex.StartROI(vert2show) PairDistVertex.EndROI(vert2show) PairDistVertex.VertexROI(vert2show)];
figure(11); clf; cmap_SVE=nebula(3);
imshow2(PCresult{f2show}.ref_im,[250 1000]); hold all;
showV=[5 12 13]; cmap_V=[1 1 0.2];
for v=showV
    [~, pth]=geodesic_distance(swcmask,ftcoord(ROIlist(v,1),:),ftcoord(ROIlist(v,2),:));
    bw_skel = bwskel(pth);
    [r, c]  = find(bw_skel);
    scatter(c,r,5,cmap_V,'filled');
    scatter(ftcoord(ROIlist(v,:),1),ftcoord(ROIlist(v,:),2),100,cmap_SVE,'filled')
end
drawScaleBar(100/Pixelsize(f2show),'horizontal','color',[1 1 1])
set_figsize(200,100)

%%
f2show=20;
vert2show=find(cell2mat(PairDistVertex.FileInd)==f2show);
swcmask=point2img(PCresult{f2show}.swc,2,size(PCresult{f2show}.ref_im));
ftprnts=PCresult{f2show}.Ftprnts(:,:,PCresult{f2show}.sortDist);
ftcoord=get_coord(ftprnts);
ROIlist=[cell2mat(PairDistVertex.StartROI(vert2show)) cell2mat(PairDistVertex.EndROI(vert2show)) cell2mat(PairDistVertex.VertexROI(vert2show))];
figure(11); clf;
imshow2(PCresult{f2show}.ref_im,[]); hold all;
for v=1:length(vert2show)
    [~, pth]=geodesic_distance(swcmask,ftcoord(ROIlist(v,1),:),ftcoord(ROIlist(v,2),:));
    bw_skel = bwskel(pth);
    [r, c]  = find(bw_skel);
    scatter(c,r,10,'filled');
    scatter(ftcoord(ROIlist(v,:),1),ftcoord(ROIlist(v,:),2),'ro')
end


%% Plot decay length over time (Figure S10 H, I)
path2show=find(cell2mat(Lstaall.StimROI_percentile)>0.9 & cell2mat(Lstaall.Fit_R2)>0.2 & tInd); % Percentile, R^2
Stim2show=unique(cell2mat(Lstaall.StimInd(path2show)));
Lsta2show=ismember(cell2mat(Lstaall.StimInd),Stim2show) & cell2mat(Lstaall.Fit_length)<300 & cell2mat(Lstaall.Fit_length)>20;
fprintf('# of sample: %2.0f\n',length(Stim2show));

figure(155); clf;
for i=1:length(Stim2show)
    Dat2showInd=find(ismember(cell2mat(Lstaall.StimInd),Stim2show(i)));
    [dax, dd]=cellfun(@(x) sort(x(:,1)), Lstaall.dendaxis(Dat2showInd),'UniformOutput',0);
    SortedVoltageSTA=cellfun(@(x,y) x(y), Lstaall.VoltageSTA(Dat2showInd),dd,'UniformOutput',0);
    t2read=cell2mat(cellfun(@(x) sort(x(:,1)), Lstaall.ReadBlueTime(Dat2showInd),'UniformOutput',0));
    [t2read, tsort]=sort(t2read);
    KymoSta=cell2mat(SortedVoltageSTA(tsort)');
    Lc2print=[]; g=1;
    for t2print=15:10:35
        ind2print=ismember(Dat2showInd,find(cell2mat(Lstaall.ReadBlueTime)==t2print));
        Lc2print=[Lc2print Lstaall.Fit_length{Dat2showInd(ind2print)}];
    end
    nexttile([1 1]);
    imagesc(t2read,[1:size(KymoSta,1)],KymoSta);
    set_kymoYtick(dax{1});
    title(sprintf('StimInd: %2.0f, Lc: %3.0f %3.0f %3.0f',Stim2show(i),Lc2print(1),Lc2print(2),Lc2print(3)));
end

data2showall=[];
for i=1:length(Stim2show)
    Dat2showInd=find(ismember(cell2mat(Lstaall.StimInd),Stim2show(i)));
    [dax, dd]=cellfun(@(x) sort(x(:,1)), Lstaall.dendaxis(Dat2showInd),'UniformOutput',0);
    SortedVoltageSTA=cellfun(@(x,y) x(y), Lstaall.VoltageSTA(Dat2showInd),dd,'UniformOutput',0);
    t2read=cell2mat(cellfun(@(x) sort(x(:,1)), Lstaall.ReadBlueTime(Dat2showInd),'UniformOutput',0));
    [t2read, tsort]=sort(t2read);
    KymoSta=cell2mat(SortedVoltageSTA(tsort)');
    Lc2print=[]; g=1;
    for t=1:length(t2read)
        ind2print=ismember(Dat2showInd,find(cell2mat(Lstaall.ReadBlueTime)==t2read(t)));
        Lc2print=[Lc2print Lstaall.Fit_length{Dat2showInd(ind2print)}];
        data2showall{i,t}=[dax{1} KymoSta(:,t)];
    end
end

timeevolutionMat=[];
for t=1:size(data2showall,2)
    localstim_binResult=binning_data(data2showall(:,t),dendbin);
    timeevolutionMat=[timeevolutionMat localstim_binResult.mean'];
end

figure(154); clf;
tiledlayout(2,1,'TileSpacing','loose')
cmap_t2=gen_colormap([0 0.6 1; 1 0 0; 0 0 0],length(t2read));
cb=colorbar; cb.Label.String='Time after blue onset (ms)';
cb.Ticks=cb.Limits; cb.TickLabels=[t2read(1) t2read(end)];

ax2=nexttile([1 1]);
pcolor(t2show,localstim_binResult.centers,timeevolutionMat);
ylim([0 500]); xlim([-5 45]); axis tight;
cb=colorbar; cb.Label.String='∆F/F';
xlabel('Time after blue onset (ms)')
ylabel('Distance from stimulation ROI (µm)')
set(gca,'ydir','reverse'); colormap(turbo);

d2show=[cell2mat(Lstaall.ReadBlueTime(Lsta2show)) cell2mat(Lstaall.Fit_length(Lsta2show))];
LcData_binning=binning_data({d2show},[-7.5:5:60]);
d2showCell=[];
for m=1:length(LcData_binning.membership)
    d2showCell{m}=[d2show(LcData_binning.membership{m},2)];
end
ax3=nexttile([1 1]);
showtInd=[4:8]; cmap_t=hsv(length(showtInd));
p=Violin_wPoints(d2showCell(showtInd),cmap_t);
line([1 5],[305 305],'LineWidth',2,'color','k');
text([3],[325],'N.S.','HorizontalAlignment','center','FontSize',20);
set(gca,'xtick',[1:length(showtInd)],'XTickLabel',LcData_binning.centers(showtInd))
ylabel('Decay length (µm)'); xlabel('Time after blue onset (ms)'); ylim([0 325]);
set_font('Arial'); set_fontsize(22);
set_figsize(170,320)

%% ===== Reviewer: 20-ms subthreshold correlation over time (thirds, f=20) =====
% Bleaching / drift control: split the recording into three equal time-thirds and
% recompute the 20-ms cross-channel subthreshold correlation matrix in each (silent
% frames, same leaf order). Stable matrices / correlation-vs-distance across thirds
% => the dynamics are not a photobleaching or drift artifact.
f = 20;
Pr  = PCresult{f};
ch1 = Pr.NormalizedTrace_ch{1};  ch2 = Pr.NormalizedTrace_ch{2};  trC = Pr.NormalizedTrace_dirt;
[nROI, nTime] = size(trC);
spk = Pr.SpikeMat(1,:) > 0;
Dpair = Pr.interDendDist(Pr.sortDist, Pr.sortDist);  Dpair = max(cat(3,Dpair,Dpair'),[],3);
minPairs = 100;  win20 = 20;
goodvec = ~imdilate(Pr.BlueStim>0, ones(1,201));
% silent period = blue-off AND not peri-spike ( +/- [-3:10] around spikes, + CS )
perispk = unique([tovec(find(spk)' + (-3:10)); find(Pr.CStrace)']);
perispk(perispk < 1 | perispk > nTime) = [];
silentvec = true(1, nTime);  silentvec(perispk) = false;
silentMask = goodvec & silentvec;

% 20-ms subthreshold (per channel + combined); keep motion/dirt as NaN
sub1 = get_subthreshold(ch1, spk, 7, win20);  sub1(isnan(ch1)) = NaN;
sub2 = get_subthreshold(ch2, spk, 7, win20);  sub2(isnan(ch2)) = NaN;
subC = get_subthreshold(trC, spk, 7, win20);  subC(isnan(trC)) = NaN;

% full-recording matrix -> leaf order (applied to every third so blocks line up)
cmFull = get_corrTensor(sub1, sub2, subC, silentMask, 0, minPairs);
C20 = cmFull(:,:,1);
Dd  = 1 - C20;  Dd = (Dd+Dd')/2;  Dd(1:nROI+1:end) = 0;  Dd(isnan(Dd)) = 1;
Zc  = linkage(squareform(Dd),'average');
try, ord = optimalleaforder(Zc, squareform(Dd)); catch, fT=figure('Visible','off'); [~,~,ord]=dendrogram(Zc,0); close(fT); end
cmapCorr = gen_colormap([0 0.5 1; 1 1 1; 1 0 0], 256);

% per-third correlation matrices
edgesT = round(linspace(0, nTime, 4));  Cmat_third = cell(1,3);  thirdMask = false(3,nTime);
for t = 1:3
    thirdMask(t, edgesT(t)+1:edgesT(t+1)) = true;
    cm = get_corrTensor(sub1, sub2, subC, silentMask & thirdMask(t,:), 0, minPairs);
    Cmat_third{t} = cm(:,:,1);
end
tlab3 = {'1st third','2nd third','3rd third'};

figure(210); clf;
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for t = 1:3
    nexttile(t);
    imagesc(Cmat_third{t}(ordUse,ordUse), [-0.25 1]); axis image; colormap(gca,cmapCorr);
    title(sprintf('%s (%d frm)', tlab3{t}, nnz(silentMask & thirdMask(t,:))));
    if t==1, ylabel('ROI (leaf order)'); end
end
cb = colorbar; cb.Layout.Tile = 'east';
title(tl, sprintf('f=%d  20-ms subthreshold correlation over time (thirds)', f));
set_figsize(420,160);

triu1 = triu(true(nROI),1);
fprintf('\n[Rev] 20-ms correlation-matrix similarity across time-thirds (f=%d, off-diag Pearson r):\n', f);
prs = [1 2; 2 3; 1 3];
for p = 1:size(prs,1)
    a = Cmat_third{prs(p,1)}(triu1);  b = Cmat_third{prs(p,2)}(triu1);  ok = ~isnan(a)&~isnan(b);
    fprintf('   third %d vs %d : r = %.4f , mean|Δ| = %.4f\n', prs(p,1),prs(p,2), corr(a(ok),b(ok)), mean(abs(a(ok)-b(ok))));
end

figure(211); clf; hold on;
xd = Dpair(triu1);
for t = 1:3
    yc = Cmat_third{t}(triu1);  ok = ~isnan(xd)&~isnan(yc);
    [xs,si] = sort(xd(ok));  ys = yc(ok);  ys = ys(si);
    plot(xs, movmean(ys,15), '-','LineWidth',1.5,'DisplayName',tlab3{t});
end
xlabel('dendritic distance (µm)'); ylabel('pairwise correlation'); legend; box off;
title(sprintf('f=%d  20-ms correlation vs distance over time', f));
set_font('Arial'); set_fontsize(14); set_figsize(220,180);

%% ===== Reviewer: peak-triggered average KYMOGRAPH over time (thirds, f=23) =====
% Bleaching / drift control: recompute the subthreshold peak-triggered-average
% kymograph (measuring ROIs sorted by distance x time, triggered on one ROI's
% peaks) in each time-third. Peaks are detected as in the Peak-triggered-average
% cell (band-stop-filtered per-ROI subthreshold; peri-spike and blue-on excluded).
% A stable kymograph across thirds => no bleaching/drift change in the spatial
% peak dynamics.
fp   = 23;
Prp  = PCresult{fp};
[nRp, nTp] = size(Prp.Subthreshold);
SubN = interpolateNaN(Prp.Subthreshold ./ Prp.SpikeHeight);       % spike-height units
StopFreq = [12 200];

perispk = unique([tovec(find(Prp.SpikeMat(1,:))' + (-10:40)); find(Prp.CStrace)']);
perispk(perispk < 1 | perispk > nTp) = [];
badp = ind2vec(nTp, perispk, 1) > 0 | imdilate(Prp.BlueStim > 0, ones(1,201));

peakvec = zeros(nRp, nTp);
for n = 1:nRp
    ft = get_bandstop(SubN(n,:), 1000, StopFreq);  ft = ft - movmedian(ft, 300, 2);
    [~, locs] = findpeaks(ft, 'MinPeakHeight',0.2, 'MinPeakDistance',50, 'MinPeakProminence',0.1);
    pv = ind2vec(nTp, locs, 1);  pv(badp) = 0;  peakvec(n,:) = pv;
end

% STA measured on the RAW subthreshold (as PeakSTA in Fig 2F); detection above used SubN
SubRaw = interpolateNaN(Prp.Subthreshold);
[dendaxissorted, dsort] = sort(PCresult{fp}.dendaxis);

trigROIs = [3 39];                        % same triggering ROIs as Fig 2F (line 837, showInd([1 3]))
edgesTp  = round(linspace(0, nTp, 4));
preP = 1000;  postP = 1000;  tauP = (-preP:postP);   % same +/-1 s window as Fig 2F
cax  = [-0.01 0.03];                       % same color range as Fig 2F
labs3 = {'1st third','2nd third','3rd third'};

figure(213); clf;
tl = tiledlayout(numel(trigROIs), 3, 'TileSpacing','compact','Padding','compact');
nPk = zeros(numel(trigROIs),3);
for ri = 1:numel(trigROIs)
    r0 = trigROIs(ri);
    for t = 1:3
        inThird = false(1,nTp);  inThird(edgesTp(t)+1:edgesTp(t+1)) = true;
        trigVec = double(peakvec(r0,:) > 0 & inThird);
        nPk(ri,t) = nnz(trigVec);
        sta = get_STA(SubRaw, trigVec, preP, postP);            % [nRp x L], all ROIs
        sta = sta - median(sta(:, 1:preP), 2);                  % pre-peak baseline per ROI
        nexttile((ri-1)*3 + t);
        imagesc(tauP, 1:nRp, sta(dsort,:), cax);  set_kymoYtick(dendaxissorted);
        xlim([-500 500]);
        if ri==1, title(sprintf('%s (%d pk)', labs3{t}, nPk(ri,t))); end
        if t==1,  ylabel(sprintf('ROI %d\ndist. from soma (µm)', r0)); end
        if ri==numel(trigROIs), xlabel('time from peak (ms)'); end
    end
end
colormap(turbo);  cb = colorbar;  cb.Layout.Tile = 'east';  cb.Label.String = '∆F/F';
title(tl, sprintf('f=%d  peak-triggered average kymograph over time (Fig 2F ROIs)', fp));
set_font('Arial'); set_fontsize(12); set_figsize(480, 160*numel(trigROIs));

fprintf('\n[Rev] f=%d peak-triggered kymograph peak counts per time-third:\n', fp);
for ri = 1:numel(trigROIs)
    fprintf('   trigger ROI %d : %s\n', trigROIs(ri), mat2str(nPk(ri,:)));
end

%% ===== Local functions =====
function show_peak_footprint(axM, axT, axP, f, r0, t0, ratVal, PCresult, Pixelsize, ...
    avgwin_ext, dref_ext, dband_ext)
% Draw the 3-panel peak footprint into the provided axes (morphology axM, traces
% axT, voltage-vs-distance axP) for one peak: neuron f, trigger ROI r0, time t0,
% with far/trigger voltage ratio ratVal.
ftprnts = PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist);
nROI    = size(ftprnts,3);
Dmat = PCresult{f}.interDendDist(PCresult{f}.sortDist, PCresult{f}.sortDist);
Dmat = max(cat(3, Dmat, Dmat'), [], 3);
SubNorm = interpolateNaN(PCresult{f}.Subthreshold / PCresult{f}.SpikeHeight); % raw subthreshold

W    = t0 + (-avgwin_ext:avgwin_ext);  W(W<1 | W>size(SubNorm,2)) = [];
Vall = mean(SubNorm(:,W), 2, 'omitnan');                 % window-mean voltage per ROI
dall = Dmat(r0,:)';                                      % contour distance from trigger
band = abs(dall - dref_ext) <= dband_ext;               % far reference ROIs

% --- (1) morphology colored by window-mean voltage; trigger ROI outlined ---
axes(axM); cla(axM);
showScaleImageSmooth(ftprnts > 0, Vall', turbo(256), [0 prctile(Vall,99)]); hold all;
bd = cell2mat(bwboundaries(imdilate(ftprnts(:,:,r0) > 0, strel('diamond',3))));
plot(bd(:,2), bd(:,1), 'color',[1 1 1], 'linewidth',1.5);
title(sprintf('f%d, frame %d, ROI %d', f, t0, r0));
axis equal tight off;
drawScaleBar(50/Pixelsize(f), 'horizontal', 'color',[1 1 1]);

% --- (2) all traces sorted by distance from trigger; averaging window shaded ---
axes(axT); cla(axT);
win = 150;  off = 0.4;
tt  = t0 + (-win:win);  tt(tt<1 | tt>size(SubNorm,2)) = [];
[~, dorder] = sort(dall);
yb  = -(nROI-1)*off - 0.4;
fill([-avgwin_ext avgwin_ext avgwin_ext -avgwin_ext], [yb yb 0.6 0.6], ...
    [1 0.9 0.9], 'EdgeColor','none'); hold all;
for k = 1:nROI
    r   = dorder(k);
    col = [0.6 0.6 0.6];
    if r==r0, col = [1 0 0]; elseif band(r), col = [0 0.4 1]; end   % trigger red, far band blue
    plot(tt-t0, SubNorm(r,tt) - (k-1)*off, 'color',col);
end
title('Traces (near \rightarrow far)');
xlabel('Time from peak (ms)');
axis tight; box off; set(gca,'ytick',[]);
drawScaleBar(0.5, 'vertical', 'color','k', 'position',[win-90 yb]);
drawScaleBar(50,  'horizontal','color','k', 'position',[win-90 yb]);

% --- (3) window-mean voltage vs distance; trigger, far band, and ratio ---
axes(axP); cla(axP);
yspan = [min(Vall) max(Vall)];
fill([dref_ext-dband_ext dref_ext+dband_ext dref_ext+dband_ext dref_ext-dband_ext], ...
    [yspan(1) yspan(1) yspan(2) yspan(2)], [0.85 0.92 1], 'EdgeColor','none'); hold all;
plot(dall, Vall, 'k.', 'markersize',12);
Vfar = mean(Vall(band), 'omitnan');
yline(Vall(r0), 'r-', 'linewidth',1.5, 'Label','V_{trig}');
yline(Vfar,     'b-', 'linewidth',1.5, 'Label',sprintf('V_{~%dµm}',dref_ext));
yline(0, 'k:');
xlabel('Distance from trigger (µm)'); ylabel('Window-mean voltage');
title(sprintf('ratio = %.2f', ratVal));
box off; xlim([0 max(dall)*1.05]);
set_font('Arial'); set_fontsize(12);
drawnow;
end

function r0 = trigger_roi(SubNorm, t0, avgwin)
% Trigger ROI = ROI with the largest window-mean subthreshold voltage at frame t0.
nTime = size(SubNorm,2);
W = t0 + (-avgwin:avgwin);  W(W<1 | W>nTime) = [];
V = mean(SubNorm(:,W), 2, 'omitnan');
[~, r0] = max(V);
end

function [amp, L50, prof] = event_extent(SubNorm, Dmat, r0, t0, avgwin, Dmax)
% amp  = trigger window-mean voltage (spike-height units)
% prof = [contour-distance-from-r0 , normalized-amplitude] for ROIs within Dmax
% L50  = contour radius at which the normalized profile first crosses 0.5, i.e.
%        the half-decay distance from the trigger. Bounded to [0, Dmax]. Small
%        L50 -> local; large L50 -> broad.
nTime = size(SubNorm,2);
W = t0 + (-avgwin:avgwin);  W(W<1 | W>nTime) = [];
V = mean(SubNorm(:,W), 2, 'omitnan');
amp = V(r0);
prof = [];  L50 = NaN;
if ~(amp > 0), return; end
Anorm = V / amp;
d = Dmat(r0,:)';
sel = d <= Dmax & ~isnan(Anorm);
prof = [d(sel) Anorm(sel)];
% L50 = contour radius where the profile first crosses 0.5
[ds, si] = sort(d(sel));  asel = Anorm(sel);  as = asel(si);
below = find(as <= 0.5 & ds > 0, 1);
if isempty(below)
    L50 = Dmax;                             % never drops below 0.5 within Dmax
elseif below == 1
    L50 = 0;
else
    i1 = below-1;  i2 = below;              % i1 closer (>0.5), i2 first <=0.5
    dn = as(i2) - as(i1);
    if dn == 0, L50 = ds(i1);
    else,       L50 = ds(i1) + (0.5-as(i1))*(ds(i2)-ds(i1))/dn;
    end
end
end

function [Xf, txt] = bandpass_rows(X, band, fs, ord)
% Zero-phase bandpass of every ROW of X (nROI x nTime). Returns X unchanged when band
% is empty, so the caller can switch the filter off with [].
% Butterworth as SECOND-ORDER SECTIONS, not b/a: a 0.1 Hz corner at 1 kHz is a
% normalized frequency of 2e-4, where direct-form coefficients lose their precision and
% the filter can return garbage or NaN. sos + filtfilt is stable there, and filtfilt
% makes it zero-phase, so nothing in the trace moves in time.
% If the Signal Processing Toolbox is missing, falls back to a moving-average
% approximation (long mean subtracted = high pass, short mean = low pass) and says so
% in txt, which the caller prints — the measurement is never silently a different one.
Xf = X;  txt = 'no bandpass (raw trace)';
if isempty(band), return; end
band = sort(band(:))';
if exist('butter','file') && exist('filtfilt','file') && exist('zp2sos','file')
    Wn = band / (fs/2);
    Wn = min(max(Wn, 1e-6), 1-1e-6);
    [z, p, k] = butter(ord, Wn, 'bandpass');
    [sos, g] = zp2sos(z, p, k);
    Xf = filtfilt(sos, g, double(X)')';
    txt = sprintf('%g-%g Hz Butterworth (order %d/edge, sos, zero phase)', ...
        band(1), band(2), ord);
else
    nLo = max(3, round(fs/band(1)));            % long window ~ the high-pass period
    nHi = max(1, round(fs/band(2)));            % short window ~ the low-pass period
    Xf = movmean(double(X) - movmean(double(X), nLo, 2, 'omitnan'), nHi, 2, 'omitnan');
    txt = sprintf('%g-%g Hz moving-average approximation (no Signal Processing Toolbox)', ...
        band(1), band(2));
end
end

function draw_roi_set_maps(axA, axS, axC, ftMask, coords, A, S, info, Adj, r0, Vthr, ttl)
% One event drawn on the morphology, three panels.
%   axA  AMPLITUDE MAP: each footprint painted with its amplitude at the event instant.
%        The colormap is grey below Vthr and turbo above it, so the panel answers "who
%        crossed threshold" and "how big was it" at the same time.
%   axS  ROI SETS: the same footprints coloured by which set they ended up in. Grey =
%        below threshold. A tested single-ROI gap is outlined in magenta, filled when
%        it was NOT trusted (so the sets either side were merged). The dark line is
%        THE PATH lenSet measures — the dendritic route between the set's two
%        farthest members — with rings on the two ends and the length beside it.
%   axC  CONNECTED SETS: the dendritic adjacency graph. Every neighbour link is drawn,
%        coloured when it joins two members of one set and left thin grey otherwise.
%        This is the panel that explains the other two — two blobs are one set exactly
%        when a coloured path runs between them.
% The images are painted as RGB directly rather than through a scaled-colour call, so
% the discrete set colours stay exact and ROIs cannot be silently blended.
nR = size(ftMask,3);
[H, W, ~] = size(ftMask);
setOf = zeros(nR,1);
for c = 1:numel(S), setOf(S(c).roi) = c; end
cmapSet = lines(max(numel(S),1));
cmx = max([max(A), Vthr*1.5, eps]);
grey = [0.72 0.72 0.72];

% --- i) amplitude map ---
cmapA = turbo(256);
nLow = max(1, min(256, round(Vthr/cmx*256)));
cmapA(1:nLow,:) = repmat(grey, nLow, 1);          % below threshold reads as grey
rgbA = zeros(nR,3);
for r = 1:nR
    v = A(r);
    if isnan(v), rgbA(r,:) = [1 1 1]*0.95;  continue; end
    ix = max(1, min(256, round(v/cmx*255)+1));
    rgbA(r,:) = cmapA(ix,:);
end
image(axA, paint_rois(ftMask, rgbA));
colormap(axA, cmapA);  clim(axA, [0 cmx]);
cb = colorbar(axA, 'Location','southoutside');   % horizontal: does not squeeze the map
cb.Label.String = 'amp';  cb.Ticks = [0 Vthr cmx];
cb.TickLabels = compose('%.2f', [0 Vthr cmx]);
axis(axA, 'image', 'off');
title(axA, sprintf('%s\ni) amplitude (grey < %.2f)', ttl, Vthr), 'FontSize',8, ...
    'FontWeight','normal');

% --- ii) the sets ---
% no ROI numbers here: figure 165 puts three panels per row and they collide.
% Change the last argument to true if you want them.
map_sets_panel(axS, ftMask, coords, S, info, r0, false);
title(axS, sprintf('ii) %d set(s), %d gap merged', numel(S), sum(info.gaps(:,4))), ...
    'FontSize',8, 'FontWeight','normal');

% --- iii) the connectivity that produced them ---
bgAny = any(ftMask, 3);
image(axC, repmat(1 - 0.20*double(bgAny), 1, 1, 3));  hold(axC,'on');
[ii, jj] = find(triu(Adj, 1));
for e = 1:numel(ii)
    a = ii(e);  b = jj(e);
    sameSet = setOf(a) > 0 && setOf(a) == setOf(b);
    if sameSet
        plot(axC, coords([a b],1), coords([a b],2), '-', 'LineWidth',3, ...
            'Color',cmapSet(setOf(a),:));
    else
        plot(axC, coords([a b],1), coords([a b],2), '-', 'LineWidth',1.2, ...
            'Color',[0.55 0.55 0.55]);
    end
end
for r = 1:nR
    if setOf(r) > 0
        plot(axC, coords(r,1), coords(r,2), 'o', 'markersize',5, 'LineWidth',0.8, ...
            'MarkerFaceColor',cmapSet(setOf(r),:), 'MarkerEdgeColor','k');
    else
        plot(axC, coords(r,1), coords(r,2), 'o', 'markersize',3.5, ...
            'MarkerFaceColor',grey, 'MarkerEdgeColor','none');
    end
end
for g = 1:size(info.gaps,1)
    w = info.gaps(g,1);
    mk = 'o';  fc = 'none';
    if info.gaps(g,4), fc = [0.8 0 0.8]; end
    plot(axC, coords(w,1), coords(w,2), mk, 'markersize',11, 'LineWidth',1.8, ...
        'Color',[0.8 0 0.8], 'MarkerFaceColor',fc);
end
plot(axC, coords(r0,1), coords(r0,2), 'k^', 'MarkerFaceColor','k', 'markersize',6);
axis(axC, 'image', 'off');  set(axC, 'YDir','reverse');
title(axC, sprintf('iii) connectivity (%d links)', numel(ii)), 'FontSize',8, ...
    'FontWeight','normal');
end

function map_sets_panel(ax, ftMask, coords, S, info, r0, labelROIs, setName, showLegend)
% The morphology painted by SET membership. Grey = below threshold. A tested
% single-ROI gap is outlined magenta, filled when it was not trusted (sets merged).
% Black triangle = the ROI the peak was detected in.
%   labelROIs : print each set member's ROI number on its footprint. Default FALSE —
%               it clutters the panel and is only worth it when you are cross-checking
%               a specific ROI against the profile panel or PkMeta. Pass true then.
if nargin < 7 || isempty(labelROIs), labelROIs = false; end
if nargin < 8 || isempty(setName),  setName = 'depolarized ROIs'; end
if nargin < 9 || isempty(showLegend), showLegend = false; end
%   setName    : what to CALL the coloured ROIs in the legend. They are the ROIs above
%                Vthr that are contiguous with the peak along the dendrite, so a label
%                of "above threshold" alone would be wrong — a far ROI above threshold
%                but not connected is deliberately NOT one of them.
%   showLegend : draw the legend (default false; the small multi-panel figures do not
%                have room for it)
nR = size(ftMask,3);
setOf = zeros(nR,1);
for c = 1:numel(S), setOf(S(c).roi) = c; end
cmapSet = lines(max(numel(S),1));
rgbS = repmat([0.72 0.72 0.72], nR, 1);
for r = 1:nR
    if setOf(r) > 0, rgbS(r,:) = cmapSet(setOf(r),:); end
end
image(ax, paint_rois(ftMask, rgbS));  hold(ax,'on');
for g = 1:size(info.gaps,1)
    w = info.gaps(g,1);
    b = bwboundaries(ftMask(:,:,w));
    if isempty(b), continue; end
    b = b{1};
    if info.gaps(g,4)
        fill(ax, b(:,2), b(:,1), [0.8 0 0.8], 'FaceAlpha',0.45, ...
            'EdgeColor',[0.8 0 0.8], 'LineWidth',1.5);
    else
        plot(ax, b(:,2), b(:,1), '-', 'Color',[0.8 0 0.8], 'LineWidth',1.5);
    end
end
% THE PATH lenSet MEASURES: the dendritic route between the set's two farthest
% members, drawn through the ROI centroids with its ends marked. Without it the
% length is a number beside a picture; with it you can see the route and check that
% it goes where you think it does.
if ~isempty(coords)
    for c = 1:numel(S)
        if ~isfield(S(c),'path') || numel(S(c).path) < 2, continue; end
        xy = coords(S(c).path, :);
        plot(ax, xy(:,1), xy(:,2), '-', 'Color','w', 'LineWidth',4);
        plot(ax, xy(:,1), xy(:,2), '-', 'Color',cmapSet(c,:)*0.55, 'LineWidth',1.8);
        plot(ax, xy([1 end],1), xy([1 end],2), 'o', 'MarkerSize',7, 'LineWidth',1.8, ...
            'MarkerFaceColor','w', 'MarkerEdgeColor',cmapSet(c,:)*0.55);
        % label at the middle OF THE PATH, not between its endpoints: on a route that
        % bends round a branch point the straight-line midpoint is off the dendrite
        mid = xy(max(1,round(size(xy,1)/2)), :);
        % nudged off the path so it does not sit on top of an ROI number
        off = 0.035 * size(ftMask,1);
        text(ax, mid(1), mid(2) - off, sprintf(' %.0f µm ', S(c).lenSpan), ...
            'Color',cmapSet(c,:)*0.55, 'FontSize',8, 'FontWeight','bold', ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'BackgroundColor','w', 'Margin',0.5);
    end
    % ROI numbers on the members, drawn BEFORE the trigger marker so that stays on
    % top. These are the sortDist indices used everywhere else — the status bar, the
    % printouts, PkMeta — so a blob here can be matched to a point in the profile
    % panel and to a peak identity.
    if labelROIs
        for c = 1:numel(S)
            % white on a dark set colour, black on a light one — 'yellow' and 'green'
            % from lines() are light enough that white text disappears on them
            lum = 0.299*cmapSet(c,1) + 0.587*cmapSet(c,2) + 0.114*cmapSet(c,3);
            txtCol = [1 1 1];
            if lum > 0.55, txtCol = [0 0 0]; end
            for rr = S(c).roi
                text(ax, coords(rr,1), coords(rr,2), sprintf('%d', rr), ...
                    'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                    'FontSize',7, 'FontWeight','bold', 'Color',txtCol, 'Clipping','on');
            end
        end
    end
    hTrig = plot(ax, coords(r0,1), coords(r0,2), 'k^', 'MarkerFaceColor','k', ...
        'markersize',6);
end
axis(ax, 'image', 'off');
% A legend naming what the colours are. setName is the term that goes in the figure
% caption too, so it is written once here rather than being retyped per panel.
if showLegend && ~isempty(S)
    hSet  = patch(ax, 'XData',nan, 'YData',nan, 'FaceColor',cmapSet(1,:), 'EdgeColor','none');
    hOth  = patch(ax, 'XData',nan, 'YData',nan, 'FaceColor',[0.72 0.72 0.72], 'EdgeColor','none');
    hPath = plot(ax, nan, nan, '-', 'Color',cmapSet(1,:)*0.55, 'LineWidth',1.8);
    hs = [hSet hOth hPath];
    ls = {setName, 'other ROIs', sprintf('spatial extent (%.0f µm)', S(1).lenSpan)};
    if exist('hTrig','var'), hs(end+1) = hTrig;  ls{end+1} = 'detecting ROI'; end
    if ~isempty(info.gaps)
        hGap = plot(ax, nan, nan, 's', 'Color',[0.8 0 0.8], 'LineWidth',1.5, ...
            'MarkerFaceColor','none', 'LineStyle','none');
        hs(end+1) = hGap;  ls{end+1} = 'tested gap ROI';
    end
    legend(ax, hs, ls, 'Box','off', 'Location','southoutside', 'FontSize',7, ...
        'AutoUpdate','off');
end
end

function img = paint_rois(ftMask, rgbPerROI)
% Paint each ROI's footprint mask with its own RGB colour on a white background.
% Later ROIs win where footprints overlap, which is only a display choice.
[H, W, nR] = size(ftMask);
img = ones(H, W, 3);
for r = 1:nR
    m = ftMask(:,:,r);
    if ~any(m(:)), continue; end
    for ch = 1:3
        tmp = img(:,:,ch);  tmp(m) = rgbPerROI(r,ch);  img(:,:,ch) = tmp;
    end
end
end

function draw_roi_sets(ax, A, sig, S, info, dendaxis, r0, Vthr, ciZ, ttl)
% One event as the set-splitting rule sees it. x = position along the dendrite so the
% sets are where they physically are; y = each ROI's amplitude at the event instant.
%   grey dots + bars : every ROI with its +/-ciZ*sigma CI
%   coloured         : the members of each set, one colour per set, with a bar under
%                      them spanning the set and labelled with its FULL length
%                      (lenSpan), the same number that goes into lenSet
%   black rings      : the two members the span is measured BETWEEN — the ends of the
%                      path drawn on the dendrite map
%   red dotted       : Vthr
%   magenta square   : a single-ROI gap that was TESTED. Filled = the gap was not
%                      trustworthy and the two sets were merged; open = the gap was
%                      believed and the sets stayed apart
%   black triangle   : the ROI peakvec fired in
hold(ax,'on');
x = dendaxis(:);
errorbar(ax, x, A, ciZ*sig, '.', 'Color',[0.65 0.65 0.65], 'markersize',8, 'CapSize',0);
plot(ax, xlim(ax), Vthr*[1 1], 'r:', 'linewidth',1);
cmap = lines(max(numel(S),1));
ylo = min([A - ciZ*sig; -0.05]);  yhi = max([A + ciZ*sig; Vthr*1.5]);
for c = 1:numel(S)
    m = S(c).roi;
    plot(ax, x(m), A(m), 'o', 'markersize',6, 'linewidth',1.2, 'Color',cmap(c,:), ...
        'MarkerFaceColor',cmap(c,:));
    yb = ylo + 0.06*(yhi-ylo)*c;
    plot(ax, [min(x(m)) max(x(m))], yb*[1 1], '-', 'Color',cmap(c,:), 'linewidth',2.5);
    text(ax, mean([min(x(m)) max(x(m))]), yb, sprintf(' %.0f µm ', S(c).lenSpan), ...
        'Color',cmap(c,:), 'FontSize',8, 'VerticalAlignment','bottom', ...
        'HorizontalAlignment','center');
    plot(ax, x(S(c).repROI), A(S(c).repROI), 'k+', 'markersize',9, 'linewidth',1.2);
    if isfield(S(c),'ends') && numel(S(c).ends) == 2   % the two ends of the span
        plot(ax, x(S(c).ends), A(S(c).ends), 'o', 'markersize',10, 'linewidth',1.5, ...
            'Color',cmap(c,:)*0.55);
    end
end
for g = 1:size(info.gaps,1)
    w = info.gaps(g,1);
    if info.gaps(g,4)                       % merged: the gap was not believed
        plot(ax, x(w), A(w), 's', 'markersize',10, 'linewidth',1.5, ...
            'Color',[0.8 0 0.8], 'MarkerFaceColor',[0.8 0 0.8]);
    else
        plot(ax, x(w), A(w), 's', 'markersize',10, 'linewidth',1.5, 'Color',[0.8 0 0.8]);
    end
end
plot(ax, x(r0), A(r0), 'k^', 'markersize',7, 'MarkerFaceColor','k');
ylim(ax, [ylo yhi]);  box(ax,'off');
ylabel(ax, 'amp');
title(ax, sprintf('%s — %d set(s), %d gap(s) tested, %d merged', ttl, numel(S), ...
    size(info.gaps,1), sum(info.gaps(:,4))), 'FontSize',9, 'FontWeight','normal');
end

function [Adj, info] = dend_adjacency(branchLab, dsoma, Dmat, maxEdge)
% The ROI neighbour graph built from the ACTUAL dendritic structure, not from distance
% ordering. Two ROIs are linked only when they are consecutive on the same branch, or
% when they are the two ends meeting at a branch point.
%   1. within each branch (BranchLabel), sort the ROIs by distance from the soma and
%      link consecutive ones. These are unambiguous neighbours: nothing lies between.
%   2. each branch's soma-proximal end is linked to the nearest ROI outside the branch,
%      which is the junction it grows from. This is the only place branches connect.
% Sorting by distance INSIDE a branch is safe — a branch is a line — whereas sorting
% globally is what interleaves separate branches and lets a walk step across the tree.
%   maxEdge : optional µm PREFERENCE for the cross-branch (junction) search. An
%             intra-branch link is real however far apart the two ROIs are — the
%             dendrite is continuous between them, just sparsely sampled — but a long
%             junction link usually means the two branches were never sampled to their
%             meeting point. [] = no cap.
% The returned graph is ALWAYS CONNECTED. maxEdge can reject junction links, and two
% sibling branches can pick each other, either of which would leave the graph in
% pieces — and since a set is a connected component, a graph in N pieces splits every
% event into N sets whatever the amplitudes are. Any remaining pieces are therefore
% joined by their closest ROI pair and reported in info.repairs, so a long repair is
% visible rather than silent.
%   info.nCompRaw : components before repair (1 = the labels alone were sufficient)
%   info.repairs  : [roiA roiB distance] for each link added to connect the graph
% Falls back to a nearest-neighbour graph if branch labels are unavailable.
nR = numel(dsoma);
Adj = false(nR);
info = struct('nIntra',0, 'nJunction',0, 'nBranch',0, 'fallback',false);
if isempty(branchLab) || all(isnan(branchLab)) || numel(branchLab) ~= nR
    D = Dmat;  D(1:nR+1:end) = inf;
    [~, nn] = min(D, [], 2);
    Adj(sub2ind([nR nR], (1:nR)', nn)) = true;
    Adj = Adj | Adj';
    info.fallback = true;
    return
end
ub = unique(branchLab(~isnan(branchLab)));
info.nBranch = numel(ub);
for b = ub(:)'
    idx = find(branchLab == b);
    if numel(idx) < 2, continue; end
    [~, o] = sort(dsoma(idx));                 % along the branch, soma outward
    idx = idx(o);
    for q = 1:numel(idx)-1
        Adj(idx(q), idx(q+1)) = true;
        info.nIntra = info.nIntra + 1;
    end
end
for b = ub(:)'                                  % branch point: proximal end -> outside
    idx = find(branchLab == b);
    if isempty(idx), continue; end
    [~, im] = min(dsoma(idx));
    pp = idx(im);
    % A branch attaches to something CLOSER TO THE SOMA. Restricting the search that
    % way is what stops two sibling branches from choosing each other and forming a
    % pair that never reaches the trunk.
    out = find(branchLab ~= b & dsoma < dsoma(pp));
    if isempty(out), out = find(branchLab ~= b); end
    if isempty(out), continue; end
    [dm, jm] = min(Dmat(pp, out));
    if ~isempty(maxEdge) && isfinite(maxEdge) && dm > maxEdge, continue; end
    Adj(pp, out(jm)) = true;
    info.nJunction = info.nJunction + 1;
end
Adj = Adj | Adj';
Adj(1:nR+1:end) = false;

% --- the graph MUST end up connected -------------------------------------------
% Sets are connected components of this graph, so a fragmented graph splits every
% event into pieces no matter what the amplitudes are: with all ROIs above threshold
% the number of sets is exactly the number of components. That is a property of the
% morphology bookkeeping, not of the data, so it is repaired here rather than being
% allowed to masquerade as structure. Repeatedly join the two closest components by
% their nearest ROI pair. maxEdge is a PREFERENCE for the junction search above, never
% a licence to leave the graph in pieces — every repair is recorded so a long one is
% visible instead of silent.
comp = graph_components(Adj);
info.nCompRaw = max(comp);
info.repairs = zeros(0,3);                      % [roiA roiB distance]
while max(comp) > 1
    best = inf;  bi = 0;  bj = 0;
    for ci = 1:max(comp)-1
        ii = find(comp == ci);
        jj = find(comp > ci);
        [v, k1] = min(Dmat(ii,jj), [], 'all', 'linear');
        if v < best
            [r1, c1] = ind2sub([numel(ii) numel(jj)], k1);
            best = v;  bi = ii(r1);  bj = jj(c1);
        end
    end
    if ~isfinite(best) || bi == 0, break; end   % nothing can be joined
    Adj(bi,bj) = true;  Adj(bj,bi) = true;
    info.repairs(end+1,:) = [bi bj best];       %#ok<AGROW>
    comp = graph_components(Adj);
end
info.nComp = max(comp);
info.nRepair = size(info.repairs,1);
end

function pth = tree_path(Adj, Dmat, a, b)
% Shortest route from ROI a to ROI b along the dendrite graph, as a list of ROIs.
% On a tree this is the only route, and it is the path the set's lenSpan measures.
% Plain Dijkstra: the graphs here are tens of nodes, so this is instant and keeps the
% function free of any dependency for something this small. Note the route may pass
% through ROIs that are NOT set members — a bridged-over gap ROI, for instance — which
% is correct: the dendrite between the two ends includes them.
n = size(Adj,1);
dist = inf(n,1);  prev = zeros(n,1);  seen = false(n,1);
dist(a) = 0;
for it = 1:n
    u = 0;  best = inf;
    for k = 1:n
        if ~seen(k) && dist(k) < best, best = dist(k);  u = k; end
    end
    if u == 0 || u == b, break; end
    seen(u) = true;
    for w = find(Adj(u,:))
        alt = dist(u) + Dmat(u,w);
        if alt < dist(w), dist(w) = alt;  prev(w) = u; end
    end
end
if ~isfinite(dist(b)), pth = [];  return; end
pth = b;
while pth(1) ~= a && prev(pth(1)) ~= 0, pth = [prev(pth(1)) pth]; end %#ok<AGROW>
if pth(1) ~= a, pth = []; end
end

function comp = graph_components(Adj)
% Connected-component labels of a small symmetric adjacency matrix, by breadth-first
% search. Written out rather than calling conncomp so this file does not depend on a
% toolbox for something this simple.
n = size(Adj,1);
comp = zeros(n,1);  c = 0;
for s0 = 1:n
    if comp(s0), continue; end
    c = c + 1;  comp(s0) = c;  q = s0;
    while ~isempty(q)
        u = q(1);  q(1) = [];
        w = find(Adj(u,:) & comp(:)' == 0);
        comp(w) = c;  q = [q w]; %#ok<AGROW>
    end
end
end

function [lenFull, dur, nROIset, nSetEv, S, info, A] = peak_set_measure(Mw, Bp, sig, Dmat, Adj, r0, t0, baseWin, avgwin, P, baseIn)
% lenSet / durSet for ONE peak: build the ROI sets at t0 and report the set the
% detecting ROI r0 belongs to. This is the wrapper the scoring pass calls, so the
% per-peak columns and the dedicated ROI-SET cell cannot drift apart — both go through
% event_roi_sets with the same arguments.
%   lenFull = the set's FULL length (both ends), i.e. 2 x lenHalf = lenSpan
%   dur     = that set's duration, from its representative trace
%   nROIset = ROIs in the set (0 if r0 is not above threshold)
%   nSetEv  = how many sets existed anywhere at this instant
% A peak whose own ROI is below threshold gets 0, not NaN: there is genuinely no
% extent to measure, which is the convention the rest of this file uses.
%   baseIn : optional nROI x 1 local baseline for this peak, precomputed by the caller
%            (movmedian once per neuron instead of a median per peak — the single
%            biggest cost in the scoring pass). [] or omitted = compute it here.
if nargin < 11 || isempty(baseIn)
    B = baseline_idx(t0, size(Mw,2), baseWin, avgwin);
    base = median(Mw(:,B), 2, 'omitnan');
else
    base = baseIn(:);
end
A = Mw(:,t0) - base;                        % every ROI at the SAME instant
[S, info] = event_roi_sets(A, sig, Bp, base, Dmat, Adj, r0, t0, P);
nSetEv = numel(S);
lenFull = 0;  dur = 0;  nROIset = 0;
if nSetEv == 0, return; end
k = find([S.hasTrigger], 1);
if isempty(k), return; end                  % r0 itself never reached threshold
lenFull = S(k).lenSpan;                     % FULL length, not the half
dur     = S(k).duration;
nROIset = S(k).nROI;
end

function [S, info] = event_roi_sets(A, sig, Bp, base, Dmat, Adj, r0, t0, P)
% Split one event into SETS of above-threshold ROIs, deciding each gap on whether the
% ROI making the gap can be trusted, and measure every set.
%
%   1. ROIs with A >= Vthr form the candidate sets: connected components on the ROI
%      neighbour graph, so "connected" means adjacent along the dendrite, not merely
%      close in sorted distance.
%   2. Two sets separated by EXACTLY ONE ROI: that ROI is below Vthr, but is it
%      RELIABLY below? Trust the gap only when even the top of its 95% CI fails to
%      reach the threshold, i.e. Vthr >= A(gap) + ciZ*sigma(gap). If the ROI is noisy
%      enough that Vthr sits inside its CI, the gap is not evidence of anything and the
%      two sets are merged.
%   3. Separated by two or more ROIs: always separate sets, no test. One unreliable ROI
%      can fake a gap; two in a row is a real interruption.
%
% Inputs
%   A    : nROI x 1 window amplitude of EVERY ROI AT THE SAME INSTANT — the +/-avgwin
%          mean of the bandpassed trace centred on t0, baseline-subtracted. It has to
%          be one common instant (or one common short window): if each ROI contributed
%          its own maximum over some search window, signals occurring at DIFFERENT
%          times would be collected into one set.
%   sig  : nROI x 1 noise of that same statistic (robust sigma)
%   Bp   : nROI x nTime bandpassed traces and base : nROI x 1 baselines, kept SEPARATE
%          so only the representative ROI's own window is baseline-subtracted for the
%          duration walk. Passing Bp - base in allocated a full copy per peak.
%   Adj  : ROI neighbour graph along the dendrite (dend_adjacency)
%   P    : Vthr, ciZ, durMax, avgwin, repMode
%
% Returns one struct per set:
%   roi        member ROIs (above threshold; a merged-over gap ROI is NOT a member,
%              it only fails to separate)
%   center     the set's 1-centre: the member whose farthest other member is nearest
%   lenSpan    contour distance between the two farthest members
%   lenHalf    lenSpan/2. THE length
%   ends       the two member ROIs that distance is measured between
%   path       the ROIs along the dendritic route from one end to the other — the
%              path lenSpan is the length of, for drawing on the morphology. Empty
%              unless P.wantPath is set, because the scoring pass never draws.
%   repROI     the set's representative trace (P.repMode), which height and duration
%              are BOTH measured on
%   height     that trace's amplitude;  duration: ms it stays above Vthr
%   maxAmp     the set's largest amplitude (equals height when repMode = 'peak')
%   hasTrigger whether this set contains the ROI the peak was detected in
% info.gaps records every single-ROI gap that was tested: [gapROI A hiCI merged].
nR = numel(A);
S = struct('roi',{}, 'center',{}, 'lenSpan',{}, 'lenHalf',{}, 'ends',{}, 'path',{}, ...
           'repROI',{}, 'height',{}, 'duration',{}, 'maxAmp',{}, 'hasTrigger',{}, ...
           'nROI',{});
info = struct('above',false(nR,1), 'gaps',zeros(0,4), 'nRaw',0);
above = A >= P.Vthr & ~isnan(A);
info.above = above;
if ~any(above), return; end

% --- 1. connected components of the above-threshold ROIs ---
lab = zeros(nR,1);  nc = 0;
for seed = find(above)'
    if lab(seed), continue; end
    nc = nc + 1;  lab(seed) = nc;  q = seed;
    while ~isempty(q)
        u = q(1);  q(1) = [];
        for w = find(Adj(u,:))
            if above(w) && ~lab(w), lab(w) = nc;  q(end+1) = w; end %#ok<AGROW>
        end
    end
end
info.nRaw = nc;

% --- 2. single-ROI gaps: merge unless the gap ROI is confidently below threshold ---
grp = 1:nc;
for w = find(~above)'
    nb = find(Adj(w,:));
    cs = unique(lab(nb(above(nb))));
    cs = cs(cs > 0);
    if numel(cs) < 2, continue; end          % this ROI does not separate anything
    hiCI = A(w) + P.ciZ * sig(w);            % top of the gap ROI's CI
    trustGap = ~isnan(hiCI) && P.Vthr >= hiCI;
    info.gaps(end+1,:) = [w A(w) hiCI double(~trustGap)]; %#ok<AGROW>
    if ~trustGap                             % the gap could be an artifact -> merge
        for j = 2:numel(cs)
            grp(grp == grp(cs(j))) = grp(cs(1));
        end
    end
end
lab(above) = grp(lab(above));

% --- 3. measure every surviving set ---
% LENGTH = half the contour distance between the set's two farthest members. A single
% ROI gives 0. HEIGHT and DURATION both come from ONE representative trace, chosen by
% P.repMode, so the two numbers describe the same signal rather than being assembled
% from different ROIs:
%   'peak'   the largest member (default) — height is then the set's maximum
%   'center' the set's 1-centre, the member whose farthest other member is nearest
for c = unique(lab(above))'
    mem = find(lab == c & above);
    D = Dmat(mem,mem);
    [~, ci] = min(max(D, [], 2));            % 1-centre of the set
    [~, pi_] = max(A(mem));
    [mxD, li] = max(D(:));
    [ia, ib] = ind2sub(size(D), li);
    sr.roi        = mem(:)';
    sr.center     = mem(ci);
    sr.lenSpan    = mxD;                     % the two farthest members
    sr.lenHalf    = sr.lenSpan / 2;          % THE length
    sr.ends       = [mem(ia) mem(ib)];       % WHICH two members that is
    % The dendritic route between them — the path lenSpan is the length OF. Only
    % computed when a viewer asks (P.wantPath): the scoring pass measures hundreds of
    % thousands of peaks and never draws any of them.
    sr.path = [];
    if isfield(P,'wantPath') && P.wantPath
        sr.path = tree_path(Adj, Dmat, sr.ends(1), sr.ends(2));
    end
    switch lower(P.repMode)
        case 'center', sr.repROI = mem(ci);
        otherwise,     sr.repROI = mem(pi_);
    end
    sr.height     = A(sr.repROI);
    sr.duration   = trace_span_above(Bp, base, sr.repROI, t0, P.Vthr, P.durMax, P.avgwin);
    sr.maxAmp     = A(mem(pi_));             % kept so 'center' can be sanity-checked
    sr.hasTrigger = ismember(r0, mem);
    sr.nROI       = numel(mem);
    S(end+1) = sr; %#ok<AGROW>
end
[~, ord] = sort([S.height], 'descend');      % biggest set first
S = S(ord);
end

function B = baseline_idx(t0, nTime, baseWin, avgwin)
% The local-baseline sample indices: t0 +/- baseWin with the event's own +/-3*avgwin
% cut out. Identical to the setdiff form it replaces (the range is already sorted and
% unique, so setdiff was only paying for a sort) and about 20x cheaper, which matters
% because this is built three times for every peak.
lo = max(1, t0 - baseWin);  hi = min(nTime, t0 + baseWin);
if hi < lo, B = 1:nTime;  return; end
keep = true(1, hi - lo + 1);
e1 = max(lo, t0 - 3*avgwin);  e2 = min(hi, t0 + 3*avgwin);
if e2 >= e1, keep(e1-lo+1 : e2-lo+1) = false; end
B = lo:hi;  B = B(keep);
if isempty(B), B = 1:nTime; end
end

function dur = trace_span_above(Bp, base, rep, t0, Vthr, durMax, avgwin)
% Contiguous ms that ROI rep's baseline-subtracted trace stays >= Vthr around t0. The
% centre is the trace's own maximum within +/-avgwin of t0, because the peak frame came
% from a different (subthreshold) trace and can sit a frame or two off. 0 when it never
% reaches Vthr, matching the convention everywhere else in this file.
% Only the +/-(avgwin+durMax) samples the walk can possibly reach are extracted, and
% the subtraction happens on that slice: passing a whole baseline-subtracted nROI x
% nTime matrix in (which is what this used to receive) allocated ~46 MB per peak.
nT = size(Bp,2);
lo0 = max(1, t0 - avgwin - durMax - 1);
hi0 = min(nT, t0 + avgwin + durMax + 1);
v = Bp(rep, lo0:hi0) - base(rep);                 % local coordinates from here on
nv = numel(v);
W = (t0 - avgwin : t0 + avgwin) - lo0 + 1;
W = W(W >= 1 & W <= nv);
if isempty(W), dur = NaN; return; end
[~, im] = max(v(W));  tc = W(im);
if ~(v(tc) >= Vthr), dur = 0; return; end
lo = max(1, tc-durMax);  hi = min(nv, tc+durMax);
i = tc;  while i > lo && v(i-1) >= Vthr, i = i-1; end
j = tc;  while j < hi && v(j+1) >= Vthr, j = j+1; end
dur = j - i;
end

function sig = roi_noise(Mw, t0, baseWin, avgwin, mode)
% Robust noise of the WINDOW-MEAN statistic, per ROI, in the same units as the
% amplitudes it will be compared against.
% This is the point of the whole exercise: the quantity being thresholded is a
% +/-avgwin mean, not a single sample, so the relevant noise is how much THAT
% fluctuates when nothing is happening — which is what the MAD of Mw measures. Using a
% sample-wise noise estimate here would overestimate it by roughly sqrt(2*avgwin+1) and
% using one number for every ROI would ignore exactly the vessel-contaminated ROIs this
% is meant to catch.
% 1.4826 x MAD rather than std, so a few real events inside the baseline window do not
% inflate the noise they are supposed to be measured against.
%   mode 'global' : the whole recording (t0 ignored)
%   mode 'local'  : t0 +/- baseWin with the event's own window cut out
nT = size(Mw,2);
if strcmpi(mode,'global')
    X = Mw;
else
    B = t0 + (-baseWin:baseWin);  B = B(B >= 1 & B <= nT);
    B = setdiff(B, t0 + (-3*avgwin:3*avgwin));
    if isempty(B), B = 1:nT; end
    X = Mw(:,B);
end
med = median(X, 2, 'omitnan');
sig = 1.4826 * median(abs(X - med), 2, 'omitnan');
sig(~(sig > 0)) = NaN;                          % a flat ROI has no scale: not testable
end

function [dt, info] = dirt_distance(NormTrDirt, motionRej, capMs)
% Milliseconds from every frame to the nearest DIRT-masked frame, per ROI.
% Dirt was blanked to NaN in NormalizedTrace_dirt and later filled by interpolateNaN,
% so downstream the trace looks continuous and nothing marks where the removal was.
% The mask is recovered here as "NaN that whole-frame motion rejection does not
% explain": motionReject blanks entire columns, dirtTrace blanks individual
% ROI x frame entries, and only the second is what this measures.
%   motionRej : the 1 x nTime motion-rejection flag, or [] when the neuron has none
%   dt   : nROI x nTime single, 0 on a dirt frame, capped at capMs
%   info : fracDirt, nROIwithDirt — printed so a neuron with no dirt is obvious
[nR, nT] = size(NormTrDirt);
mot = false(1, nT);
if ~isempty(motionRej)
    m = motionRej(:)' > 0;
    if numel(m) == nT, mot = m; end
end
dirtMask = isnan(NormTrDirt) & ~repmat(mot, nR, 1);
dt = capMs * ones(nR, nT, 'single');
for r = 1:nR
    if any(dirtMask(r,:))
        dt(r,:) = min(single(bwdist(dirtMask(r,:))), capMs);
    end
end
info = struct('fracDirt', mean(dirtMask(:)), 'nROIwithDirt', sum(any(dirtMask,2)), ...
    'nROI', nR);
end

function A = roi_areas(ftp, pixelsize, thr)
% Footprint area of every ROI in um^2. The footprints are WEIGHTED (they come out of
% the demixing, not from a hand-drawn mask), so an area needs a cut: each footprint is
% normalized to its own maximum and the pixels at or above thr are counted. Returns a
% column vector in the ROI order of ftp.
%   ftp       : H x W x nROI, already abs() and in the ROI order you want
%   pixelsize : um per pixel (PCresult{f}.pixelsize)
mx = max(max(ftp,[],1),[],2);                    % 1 x 1 x nROI
mx(mx <= 0) = eps;
A = squeeze(sum(sum(ftp >= thr.*mx, 1), 2)) * pixelsize^2;
A = A(:);
A(A <= 0) = NaN;                                 % an empty footprint has no area
end

function [aucSum, aucArea, nUse] = event_auc(Sub, Dmat, areaROI, r0, t0, aucWin, Dmax, baseWin, avgwin, baseIn)
% Area under one peak, summed over ROIs, and the same sum divided by the total ROI
% area — the "mean AUC per peak".
%   aucSum  = sum over ROIs within Dmax of the trigger of  integral( Sub - base ) dt
%             over t0 + aucWin, in spike-height x ms.
%   aucArea = aucSum / sum(areaROI of those same ROIs), spike-height x ms / um^2.
%             Dividing by AREA rather than by ROI COUNT is what makes this an
%             area-weighted mean: two neurons whose ROIs differ in size still give
%             comparable numbers.
%   nUse    = how many ROIs entered the sum.
% The integral is SIGNED and not rectified, so an ROI the event never reached
% contributes approximately zero rather than adding noise with a positive bias. That
% also means aucArea is diluted by however much quiet dendrite lies within Dmax, which
% is the intended reading of "mean over the dendrite" — but it does make the value
% depend on Dmax and on the field of view, so quote both.
% base is the per-ROI local median over t0 +/- baseWin with the event's own window
% excluded — the same baseline convention as event_len_dur and event_peak_stats, so
% these amplitudes are comparable with amp / lenSet / durSet.
aucSum = NaN;  aucArea = NaN;  nUse = 0;
nTime = size(Sub,2);
sel = find(Dmat(r0,:)' <= Dmax & ~isnan(Dmat(r0,:)'));
if isempty(sel), return; end
W = t0 + (round(aucWin(1)):round(aucWin(2)));   W(W<1 | W>nTime) = [];
if isempty(W), return; end
if nargin < 10 || isempty(baseIn)
    B = baseline_idx(t0, nTime, baseWin, avgwin);
    base = median(Sub(sel,B), 2, 'omitnan');
else
    base = baseIn(sel);  base = base(:);        % precomputed by the caller
end
X = Sub(sel,W) - base;                          % baseline-subtracted event window
% An ROI counts only if BOTH its integral and its area are measurable, so the same set
% of ROIs is behind the numerator and the denominator.
ok = ~any(isnan(X),2) & ~isnan(base) & ~isnan(areaROI(sel));
if ~any(ok), return; end
% 1 kHz sampling, so the frame index IS the time in ms; trapz integrates over it
auc = trapz(double(W(:)), double(X(ok,:))', 1)';       % one integral per usable ROI
aucSum  = sum(auc);
aucArea = aucSum / sum(areaROI(sel(ok)));
nUse    = numel(auc);
end

function profs = get_profiles(rows, EventTab, PCresult, avgwin, Dmax)
% Normalized amplitude-vs-distance profiles for a set of EventTab rows; SubNorm/
% Dmat are computed once per involved neuron for speed.
profs = cell(numel(rows),1);
for f = unique(EventTab(rows,1))'
    SubNorm = interpolateNaN(PCresult{f}.Subthreshold / PCresult{f}.SpikeHeight);
    Dmat = PCresult{f}.interDendDist(PCresult{f}.sortDist, PCresult{f}.sortDist);
    Dmat = max(cat(3,Dmat,Dmat'),[],3);
    for k = find(EventTab(rows,1)==f)'
        rr = rows(k);
        [~,~,profs{k}] = event_extent(SubNorm, Dmat, EventTab(rr,2), EventTab(rr,3), avgwin, Dmax);
    end
end
end

function draw_event_footprint(axMorph, axProf, f, r0, t0, PCresult, Pixelsize, avgwin, Dmax, L50val, ttl)
% Two-panel event footprint: (left) morphology colored by normalized window
% amplitude with the trigger ROI outlined; (right) normalized amplitude vs contour
% distance with the 0.5 line and L50/2 radius marked.
ftprnts = PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist);
Dmat = PCresult{f}.interDendDist(PCresult{f}.sortDist, PCresult{f}.sortDist);
Dmat = max(cat(3,Dmat,Dmat'),[],3);
SubNorm = interpolateNaN(PCresult{f}.Subthreshold / PCresult{f}.SpikeHeight);
nTime = size(SubNorm,2);
W = t0 + (-avgwin:avgwin);  W(W<1 | W>nTime) = [];
V = mean(SubNorm(:,W), 2, 'omitnan');
Anorm = V / V(r0);
d = Dmat(r0,:)';

axes(axMorph); cla(axMorph);
showScaleImageSmooth(ftprnts > 0, Anorm', turbo(256), [0 1]); hold all;
bd = cell2mat(bwboundaries(imdilate(ftprnts(:,:,r0) > 0, strel('diamond',3))));
plot(bd(:,2), bd(:,1), 'color',[1 1 1], 'linewidth',1.5);
colormap(axMorph, turbo(256));  caxis(axMorph, [0 1]);
cb = colorbar(axMorph);  cb.Label.String = 'Norm. amplitude';  cb.Ticks = [0 0.5 1];
title(ttl);  axis equal tight off;
drawScaleBar(50/Pixelsize(f), 'horizontal', 'color',[1 1 1]);

axes(axProf); cla(axProf);
plot(d, Anorm, 'k.', 'markersize',12); hold all;
yline(0.5, 'r--', '0.5');
xline(L50val, 'b--', 'Label','L_{50}');
yline(0, 'k:');
xlabel('Contour distance (µm)'); ylabel('Normalized amplitude');
title(sprintf('L_{50} = %.0f µm', L50val));
box off; xlim([0 Dmax]); ylim([-0.3 1.1]);
end

function sX = cross_half(s, A, pk, dir)
% Interpolated position along arc-length s where the normalized profile A crosses
% 0.5, searching from the peak index pk in direction dir (+1 right, -1 left). If
% it never drops below 0.5, returns the corresponding polyline endpoint.
n = numel(A);  i = pk;
if dir < 0
    while i > 1 && A(i) >= 0.5, i = i-1; end
    if A(i) < 0.5 && A(i+1) ~= A(i)
        sX = interp1([A(i) A(i+1)], [s(i) s(i+1)], 0.5);
    else
        sX = s(1);
    end
else
    while i < n && A(i) >= 0.5, i = i+1; end
    if A(i) < 0.5 && A(i) ~= A(i-1)
        sX = interp1([A(i-1) A(i)], [s(i-1) s(i)], 0.5);
    else
        sX = s(end);
    end
end
end

function [len, dur, ampTrig, tPre, tPost, cens, prof, Vall, tc, vTrig] = event_len_dur(Sub, Raw, Dmat, r0, t0, avgwin, Dmax, Vthr, baseWin, durMax, smoothDur)
% Threshold-crossing length and duration of one event, measured outward from the
% trigger point (peak ROI r0, peak frame t0).
% The peak was located on the SUBTHRESHOLD trace (Sub, used here only to seed the
% time search); LENGTH and DURATION are both measured on Raw, which is not broadened
% by the subthreshold filter and therefore gives shorter durations. Callers pass the
% BANDPASSED trace as Raw (see bpLenDur_PK) so the crossings are not set by blood
% noise or by slow drift against an absolute threshold.
% The last three outputs exist so a viewer can DRAW the measurement rather than just
% report it:
%   Vall  : baseline-subtracted window amplitude of EVERY ROI (the spatial profile the
%           length was read off; prof is the subset within Dmax, with distances)
%   tc    : the frame the duration was measured around, so tc+tPre and tc+tPost are the
%           crossings in absolute frames
%   vTrig : the trigger ROI's baseline-subtracted, smoothed trace — the exact signal
%           the crossings were found on
%   Sub       : nROI x nTime subthreshold traces, spike-height units (sortDist order)
%   Raw       : nROI x nTime RAW traces, same units and ROI order (NormalizedTrace_dirt)
%   Dmat      : nROI x nROI contour distance matrix (µm, same ROI order)
%   r0, t0    : trigger ROI and peak frame (from the subthreshold detection)
%   avgwin    : +/- frames of Raw averaged for the per-ROI window amplitude (profile)
%   Dmax      : ROIs farther than this from the trigger are ignored (µm)
%   Vthr      : absolute voltage threshold, same units as the traces (spike height)
%   baseWin   : +/- frames of local per-ROI baseline (event window excluded)
%   durMax    : +/- frames searched for the temporal crossings
%   smoothDur : movmean window (frames) on Raw before the duration search (1 = none)
% dur     = tPost - tPre, the outermost frames on either side of the peak at which
%           the trigger ROI's RAW trace is still above Vthr (ms; 1 kHz, no interp)
% len     = contour distance of the farthest ROI still above Vthr before the RAW
%           peak-window profile drops below it (µm, a sampled ROI distance)
% ampTrig = trigger window-mean amplitude of Raw (baseline subtracted)
% cens    = [spaceFlag timeFlag]. spaceFlag: 1 = no crossing within the sampled
%           ROIs (len = outermost ROI / Dmax), 2 =
%           only the trigger ROI is above threshold, so the extent is shorter than
%           the ROI spacing and therefore unresolved (len = 0).
%           timeFlag: 1 = no crossing within the durMax cap.
% prof    = [distance amplitude] of the peak-window profile (ROIs within Dmax)
if nargin < 11 || isempty(smoothDur), smoothDur = 1; end
nTime = size(Raw,2);
len = NaN;  dur = NaN;  tPre = NaN;  tPost = NaN;  cens = [0 0];  prof = [];
W = t0 + (-avgwin:avgwin);   W(W<1 | W>nTime) = [];
B = t0 + (-baseWin:baseWin); B(B<1 | B>nTime) = [];
B = setdiff(B, t0 + (-3*avgwin:3*avgwin));     % keep the event out of its own baseline
if isempty(B), B = 1:nTime; end
base = median(Raw(:,B), 2, 'omitnan');         % baseline from the measurement trace
V = mean(Raw(:,W), 2, 'omitnan') - base;       % +/-avgwin mean of Raw at the peak
ampTrig = V(r0);
% the drawing outputs, assigned before any early return so a viewer can always show
% what was (or was not) measurable
Vall  = V;
vTrig = movmean(Raw(r0,:) - base(r0), smoothDur, 'omitnan');
tc    = t0;
% A peak that never reaches Vthr has an EMPTY region above it, so its extent above
% Vthr is genuinely zero — not unknown. Returning 0 rather than NaN keeps the
% measurement defined for every peak in the all-peak pass (an absolute threshold that
% returned NaN would silently drop every small peak, i.e. pre-select large events),
% while staying faithful to the definition. cens = [3 3] marks these so they can be
% told apart from a peak that IS above threshold but confined to one ROI (len = 0,
% cens(1) = 2, but dur > 0).
if isnan(ampTrig)
    return;                                     % no usable trace here: stays NaN,
end                                             % which is NOT the same as "below Vthr"
if ampTrig < Vthr
    len = 0;  dur = 0;  tPre = 0;  tPost = 0;  cens = [3 3];
    return;
end

% --- LENGTH: walk out along contour distance until the profile crosses Vthr ---
d = Dmat(r0,:)';
sel = d <= Dmax & ~isnan(V) & ~isnan(d);
prof = [d(sel) V(sel)];
[ds, si] = sort(d(sel));  Vsel = V(sel);  Vs = Vsel(si);
below = find(Vs < Vthr & ds > 0, 1);           % first ROI below threshold
if isempty(below)
    len = min(max(ds), Dmax);                  % never drops within the sampled ROIs
    cens(1) = 1;
else
    i1 = find(Vs >= Vthr & ds < ds(below), 1, 'last');  % farthest ROI still above
    if isempty(i1)
        len = 0;            % only the trigger ROI is above -> shorter than ROI spacing
        cens(1) = 2;        % unresolved: below the sampling resolution
    else
        len = ds(i1);       % measured ROI distance, no interpolation
    end
end

% --- DURATION: crossings of Vthr in the trigger ROI's RAW trace, either side of the peak ---
% The peak frame is taken from the SUBTHRESHOLD trace (its maximum inside the peak
% window); only the crossings are read off the raw trace. If the raw trace happens
% to sit below threshold at that frame (noise), the raw maximum inside the same
% window is used instead.
v  = vTrig;                                    % already smoothed and baseline-subtracted
vs = Sub(r0,:) - median(Sub(r0,B), 2, 'omitnan');
[~, im] = max(vs(W));  tc = W(im);             % subthreshold peak frame
if v(tc) < Vthr
    [~, im] = max(v(W));  tc = W(im);          % fall back to the raw maximum
end
if v(tc) < Vthr, cens(2) = 1; tPre = 0; tPost = 0; dur = 0; return; end
lo = max(1, tc-durMax);  hi = min(nTime, tc+durMax);
i = tc;  while i > lo && v(i-1) >= Vthr, i = i-1; end
tPre = i - tc;                                 % last frame above threshold (ms, <=0)
if i == lo, cens(2) = 1; end                   % hit the search cap
i = tc;  while i < hi && v(i+1) >= Vthr, i = i+1; end
tPost = i - tc;                                % last frame above threshold (ms, >=0)
if i == hi, cens(2) = 1; end
dur = tPost - tPre;                            % sampled at 1 kHz, no interpolation
end

function S = event_peak_stats(Raw, Gl, r0, t0, avgwin, baseWin, P, baseIn)
% Size, local signal-to-fluctuation ratio and common-mode correlation of one peak, all
% on the RAW trace with a local median baseline (the event's own window excluded, the
% same baseline convention event_len_dur uses, so the amplitudes are comparable).
%   amp      +/-avgwin mean amplitude at the trigger ROI, baseline subtracted
%   pkAmp    largest value inside that window;  tc = the frame where it occurs
%            (tc is reused as the centre of the local-PSD windows)
%   oscRatio pkAmp / oscAmp, where oscAmp is half the 5-95 percentile range of the trace
%            in the FLANKS only (+/-P.Wac, the event itself cut out at its baseline
%            crossings). A local signal-to-fluctuation ratio: ~1 means the "peak" is no
%            bigger than whatever fluctuates around it — noise, or an ongoing ripple —
%            and >>1 means it genuinely stands out. Deliberately band-agnostic, so it
%            does not care what frequency the artifact happens to be at.
%   cmCorr   correlation of the peri-peak waveform with Gl, the mean-over-ROI trace.
%            Motion and whole-field artifacts are common-mode; a local dendritic event
%            contributes almost nothing to a mean over 20-40 ROIs.
%   wave     peri-peak waveform on +/-P.Wplot, normalized to the peak (for the viewer)
S = struct('amp',NaN, 'pkAmp',NaN, 'tc',NaN, 'oscAmp',NaN, 'oscRatio',NaN, ...
    'cmCorr',NaN, 'wave',nan(1,2*P.Wplot+1));
nTime = size(Raw,2);
% Only a window around t0 is ever read (the furthest reach is avgwin + the largest of
% Wplot / Wac / Wosc), so the smoothing is done on that slice instead of on the whole
% recording. It used to movmean the full trace TWICE per peak — 1.5 ms each on a 200 s
% recording — and then use about 1% of the result. The slice is padded by the smoothing
% width so every sample actually used still sees a full window, which makes this
% numerically identical, not an approximation.
reach = avgwin + max([P.Wplot, P.Wac, P.Wosc]) + P.smooth + 2;
lo0 = max(1, t0 - reach);  hi0 = min(nTime, t0 + reach);
off = lo0 - 1;                                  % local index = global index - off
if nargin < 8 || isempty(baseIn)
    B = baseline_idx(t0, nTime, baseWin, avgwin);
    bRaw = median(Raw(r0,B), 'omitnan');
    bGl  = NaN;
    if ~isempty(Gl), bGl = median(Gl(B), 'omitnan'); end
else
    bRaw = baseIn(1);  bGl = baseIn(min(2,end));
end
v = movmean(Raw(r0,lo0:hi0) - bRaw, P.smooth, 'omitnan');
nv = numel(v);

W = (t0-avgwin : t0+avgwin) - off;  W = W(W >= 1 & W <= nv);
if isempty(W), return; end
S.amp = mean(v(W),'omitnan');
[pk, im] = max(v(W));  tcL = W(im);             % local
S.pkAmp = pk;  S.tc = tcL + off;                % tc is reported in GLOBAL frames
if ~(pk > 0), return; end

Wp = tcL + (-P.Wplot:P.Wplot);
inW = Wp >= 1 & Wp <= nv & (Wp+off) >= 1 & (Wp+off) <= nTime;
S.wave(inW) = v(Wp(inW)) / pk;

% event edges = nearest baseline crossings, so the flanks exclude the event itself
lo = max(1, tcL-P.Wosc);  hi = min(nv, tcL+P.Wosc);
iL = tcL;  while iL > lo && v(iL-1) > 0, iL = iL-1; end
iR = tcL;  while iR < hi && v(iR+1) > 0, iR = iR+1; end
Fl = [max(1,tcL-P.Wac):max(1,iL-1), min(nv,iR+1):min(nv,tcL+P.Wac)];
Fl = Fl(Fl >= 1 & Fl <= nv);
vf = v(Fl);  vf = vf(~isnan(vf));
if numel(vf) > 50
    S.oscAmp   = 0.5*(prctile(vf,95) - prctile(vf,5));
    S.oscRatio = pk / max(S.oscAmp, eps);
end

if ~isempty(Gl)
    g = movmean(Gl(lo0:hi0) - bGl, P.smooth, 'omitnan');
    Wb = tcL + (-P.Wosc:P.Wosc);  Wb = Wb(Wb >= 1 & Wb <= nv);
    a = v(Wb);  b = g(Wb);  ok = ~isnan(a) & ~isnan(b);
    if sum(ok) > 10 && std(a(ok)) > 0 && std(b(ok)) > 0
        S.cmCorr = corr(a(ok)', b(ok)');
    end
end
end

function S = event_local_psd(Raw, r0, tc, P)
% BLOOD-NOISE content of one peak, from two power spectra centred on the same frame:
%   Pn = PSD of the NARROW window (+/-P.Wnar ms)  — the peak itself
%   Pw = PSD of the WIDE   window (+/-P.Wwide ms) — the background it sits in
% Both are power DENSITIES on one zero-padded grid, so for a stationary process they
% estimate the same thing whatever the window length — which is what makes their ratio
% meaningful. Red-blood-cell flow noise occupies P.fVess (40-250 Hz).
%   bandRatio in-band power of the narrow window / in-band power of the wide window.
%             ~1 = the vessel noise is simply always present in this ROI; >>1 = a burst
%             of flow noise landed on this peak, which is what creates a false peak.
%   bandFrac  in-band share of the narrow window's total power (1 Hz - Nyquist):
%             "is this peak MADE of blood noise?" Diluted by broadband shot noise, so
%             it separates groups better than it judges single peaks.
%   fLoc      argmax of Pw-Pn over P.fBand (Hz). Diagnostic only: a stationary rhythm's
%             periodogram peak density grows with window length, so the difference peaks
%             at the frequency of whatever is oscillating around this peak. Its
%             resolution comes from the WIDE window (0.25 Hz), not the narrow one.
% Only the MEAN is removed from each segment; a linear detrend over 201 samples would
% act as a strong high-pass and bias the low end of the spectrum.
S = struct('bandRatio',NaN, 'bandFrac',NaN, 'fLoc',NaN);
nTime = size(Raw,2);
Wn = tc + (-P.Wnar:P.Wnar);    Wn(Wn<1 | Wn>nTime) = [];
Ww = tc + (-P.Wwide:P.Wwide);  Ww(Ww<1 | Ww>nTime) = [];
xn = Raw(r0,Wn);  xn = xn(~isnan(xn));
xw = Raw(r0,Ww);  xw = xw(~isnan(xw));
if numel(xn) < P.Wnar || numel(xw) < 4*P.Wnar, return; end     % window too truncated

% Hann-windowed one-sided periodogram, computed directly. This is what
% periodogram() returns (checked to 5e-16 relative) but ~10x faster: periodogram's
% argument handling dominated its runtime here, and hann() was being rebuilt from
% scratch on every call. The windows and the frequency axis depend only on the window
% LENGTH, which is the same for every peak away from the recording edges, so they are
% cached and reused.
[Pn, fx] = hann_psd(xn, P.nfft, 1000);
Pw       = hann_psd(xw, P.nfft, 1000);

bl = fx >= P.fVess(1) & fx <= P.fVess(2);       % inside the blood-noise band
ba = fx >= 1;                                   % 1 Hz - Nyquist
inN = trapz(fx(bl), Pn(bl));
S.bandRatio = inN / max(trapz(fx(bl), Pw(bl)), eps);
S.bandFrac  = inN / max(trapz(fx(ba), Pn(ba)), eps);

b = fx >= P.fBand(1) & fx <= P.fBand(2);
ib = find(b);  [~, i1] = max(Pw(b) - Pn(b));  S.fLoc = fx(ib(i1));
end

function [Pxx, fx] = hann_psd(x, nfft, fs)
% One-sided Hann-windowed periodogram of x, identical to
% periodogram(x(:)-mean(x), hann(numel(x)), nfft, fs) but without the argument
% parsing, and with the window / frequency axis cached per window length.
persistent wCache fCache
if isempty(wCache), wCache = containers.Map('KeyType','double','ValueType','any'); end
if isempty(fCache), fCache = containers.Map('KeyType','double','ValueType','any'); end
n = numel(x);
if isKey(wCache, n)
    w = wCache(n);
else
    w = hann(n);  wCache(n) = w;
end
key = nfft*1e6 + fs;
if isKey(fCache, key)
    fx = fCache(key);
else
    fx = (0:floor(nfft/2))' * (fs/nfft);  fCache(key) = fx;
end
y = (x(:) - mean(x)) .* w;
Y = fft(y, nfft);
Pxx = abs(Y(1:numel(fx))).^2 / (fs * sum(w.^2));
Pxx(2:end-1) = 2*Pxx(2:end-1);              % fold the negative frequencies
end

function peak_classifier(PkTab, PkMeta, sampleIdx, pkLabel, Kymo, Trig, rowTrig, rowSoma, dend, exPk, cfg)
% Hand-classify one peak at a time. Left: the peri-peak kymograph of every ROI of that
% neuron, ordered by distance from the soma (basal at the bottom, apical at the top),
% with the ROI peakvec fired on drawn as a white line and the soma as a dotted one;
% below it, that ROI's raw trace. Right: a histogram of each parameter over ALL peaks
% with this peak's value in red and the five hand-picked candidates in green; far
% All twelve parameters are shown, so nothing is hidden behind a summary.
% Keys:  r real   a artifact   u unsure   <- -> browse   s save   q save and close.
% Labels are written to cfg.labFile after every keypress, so this is resumable.
st.PkTab = PkTab;  st.PkMeta = PkMeta;  st.sampleIdx = sampleIdx;  st.pkLabel = pkLabel;
st.Kymo = Kymo;  st.Trig = Trig;  st.rowTrig = rowTrig;  st.rowSoma = rowSoma;
st.dend = dend;  st.exPk = exPk;  st.cfg = cfg;
st.hp = cellfun(@(s) find(strcmp(cfg.pkNames,s),1), cfg.histParams);
st.i = find(isnan(pkLabel), 1);
if isempty(st.i), st.i = 1; end

fig = figure('Name','Peak classifier — r real / a artifact / u unsure / arrows browse / q quit', ...
    'NumberTitle','off', 'Color','w', 'KeyPressFcn',@pc_key);
st.axK = axes(fig, 'Position',[0.045 0.46 0.30 0.47]);
st.axT = axes(fig, 'Position',[0.045 0.11 0.30 0.26]);
nH = numel(st.hp);  nc = 4;  nr = ceil(nH/nc);
st.axH = gobjects(nH,1);
for k = 1:nH
    c = mod(k-1,nc);  r = floor((k-1)/nc);
    w = 0.58/nc;  h = 0.88/nr;
    st.axH(k) = axes(fig, 'Position',[0.385+c*w, 0.95-(r+1)*h, w-0.04, h-0.055]);
end
st.txt = uicontrol(fig, 'Style','text', 'Units','normalized', 'BackgroundColor','w', ...
    'Position',[0.045 0.005 0.80 0.045], 'HorizontalAlignment','left', 'FontSize',10);
uicontrol(fig,'Style','pushbutton','String','REAL (r)','Units','normalized', ...
    'Position',[0.86 0.005 0.045 0.04], 'Callback',@(s,e) pc_set(fig,1));
uicontrol(fig,'Style','pushbutton','String','ARTIFACT (a)','Units','normalized', ...
    'Position',[0.907 0.005 0.05 0.04], 'Callback',@(s,e) pc_set(fig,0));
uicontrol(fig,'Style','pushbutton','String','? (u)','Units','normalized', ...
    'Position',[0.959 0.005 0.035 0.04], 'Callback',@(s,e) pc_set(fig,-1));
guidata(fig, st);
pc_draw(fig);
end

function pc_draw(fig)
st = guidata(fig);
i = st.i;  k = st.sampleIdx(i);
f = st.PkMeta(k,1);  r0 = st.PkMeta(k,2);  t0 = st.PkMeta(k,3);
K = st.Kymo(:,:,i);
nR = find(any(~isnan(K),2), 1, 'last');
if isempty(nR), nR = size(K,1); end
d = st.dend(1:nR, i);
lblStr = {'ARTIFACT','','REAL'};

% --- kymograph: every ROI, ordered by distance from the soma ---
cla(st.axK);
sub = K(1:nR,:);  sub = sub(~isnan(sub));
if isempty(sub), cax = [-1 1]; else, cax = double(prctile(sub,[2 99.5])); end
if ~(cax(2) > cax(1)), cax = cax(1) + [-1 1]; end
imagesc(st.axK, st.cfg.tAxis, 1:nR, K(1:nR,:), cax);
set(st.axK, 'YDir','normal');  hold(st.axK,'on');
plot(st.axK, st.cfg.tAxis([1 end]), st.rowTrig(i)*[1 1], 'w-', 'linewidth',1.5);
plot(st.axK, st.cfg.tAxis(1)+0.03*range(st.cfg.tAxis), st.rowTrig(i), 'w>', ...
    'MarkerFaceColor','w', 'markersize',9);
if ~isnan(st.rowSoma(i))
    plot(st.axK, st.cfg.tAxis([1 end]), st.rowSoma(i)*[1 1], 'w:', 'linewidth',1);
end
plot(st.axK, [0 0], [0.5 nR+0.5], 'w--');
colormap(st.axK, turbo);  axis(st.axK,'tight');
yt = unique(round(linspace(1,nR,7)));
set(st.axK, 'YTick',yt, 'YTickLabel',compose('%.0f', d(yt)));
ylabel(st.axK, 'distance from soma (µm)');
title(st.axK, sprintf('f%d   ROI %d (%.0f µm, row %d/%d)   frame %d', ...
    f, r0, d(st.rowTrig(i)), st.rowTrig(i), nR, t0));

% --- the trigger ROI's own raw trace ---
cla(st.axT);
plot(st.axT, st.cfg.tAxis, st.Trig(i,:), 'k-', 'linewidth',1);  hold(st.axT,'on');
plot(st.axT, st.cfg.tAxis([1 end]), [0 0], 'k:');
plot(st.axT, [0 0], ylim(st.axT), 'r--');
xlim(st.axT, st.cfg.tAxis([1 end]));  box(st.axT,'off');
xlabel(st.axT, 'Time from peak (ms)');  ylabel(st.axT, 'Raw (spike height)');

% --- one histogram per parameter, this peak in red, candidates in green ---
v = st.PkTab(k,:);
for h = 1:numel(st.hp)
    p = st.hp(h);  ax = st.axH(h);  cla(ax);
    col = st.PkTab(~isnan(st.PkTab(:,p)), p);
    if isempty(col), axis(ax,'off'); continue; end
    lo = prctile(col,0.5);  hi = prctile(col,99.5);
    if ~(hi > lo), lo = min(col);  hi = max(col); end
    if ~(hi > lo), hi = lo + 1; end
    histogram(ax, col, linspace(lo,hi,30), 'FaceColor',[0.55 0.55 0.55], 'EdgeColor','none');
    hold(ax,'on');
    for e = 1:size(st.exPk,1)
        if ~isnan(st.exPk(e,p)), plot(ax, st.exPk(e,p)*[1 1], ylim(ax), 'g-', 'linewidth',1); end
    end
    if ~isnan(v(p))
        plot(ax, v(p)*[1 1], ylim(ax), 'r-', 'linewidth',2);
        xlim(ax, [min(lo,v(p)) max(hi,v(p))]);
    else
        xlim(ax, [lo hi]);
    end
    set(ax,'YTick',[]);  box(ax,'off');
    title(ax, sprintf('%s = %.3g', st.cfg.pkNames{p}, v(p)), 'FontSize',9, 'FontWeight','normal');
end

% --- status line ---
cur = st.pkLabel(i);
if isnan(cur), curStr = 'unlabelled'; elseif cur == -1, curStr = 'UNSURE';
else, curStr = lblStr{cur+2}; end
set(st.txt, 'String', sprintf(['peak %d/%d   this one: %s        labelled so far: ' ...
    '%d real, %d artifact, %d unsure, %d left        [r]eal  [a]rtifact  [u]nsure  ' ...
    '<- ->  [s]ave  [q]uit'], i, numel(st.sampleIdx), curStr, sum(st.pkLabel==1), ...
    sum(st.pkLabel==0), sum(st.pkLabel==-1), sum(isnan(st.pkLabel))));
end

function pc_key(fig, ev)
switch lower(ev.Key)
    case 'r', pc_set(fig, 1);   return
    case 'a', pc_set(fig, 0);   return
    case 'u', pc_set(fig, -1);  return
    case 's', pc_save(fig);     return
    case 'q', pc_save(fig);  close(fig);  return
    case 'rightarrow', pc_move(fig, +1);
    case 'leftarrow',  pc_move(fig, -1);
end
end

function pc_move(fig, step)
st = guidata(fig);
st.i = min(max(st.i + step, 1), numel(st.sampleIdx));
guidata(fig, st);  pc_draw(fig);
end

function pc_set(fig, val)
st = guidata(fig);
st.pkLabel(st.i) = val;
guidata(fig, st);
pc_save(fig);
st = guidata(fig);
nxt = find(isnan(st.pkLabel(st.i+1:end)), 1);      % jump to the next unlabelled one
if isempty(nxt), st.i = min(st.i+1, numel(st.sampleIdx));
else,            st.i = st.i + nxt;
end
guidata(fig, st);  pc_draw(fig);
end

function pc_save(fig)
st = guidata(fig);
sampleIdx = st.sampleIdx;  pkLabel = st.pkLabel;
labMeta = st.cfg.labMeta;  nPkAll = st.cfg.nPkAll;      % labMeta = the portable key
rngSeed_PK = NaN;  sampleMode_PK = '';
if isfile(st.cfg.labFile)
    old = load(st.cfg.labFile);
    if isfield(old,'rngSeed_PK'),    rngSeed_PK = old.rngSeed_PK; end
    if isfield(old,'sampleMode_PK'), sampleMode_PK = old.sampleMode_PK; end
end
save(st.cfg.labFile, 'sampleIdx','pkLabel','labMeta','nPkAll','rngSeed_PK','sampleMode_PK');
end

function v = jitter_few(v, gd)
% Spread a few-valued parameter out just enough that a scatter shows density instead of
% a line of dots on top of each other. Anything with more than 5 distinct values is
% returned untouched, so a continuous parameter is never moved.
if sum(gd) < 2, return; end
u = unique(v(gd));
if numel(u) > 5, return; end
sp = range(v(gd));  if sp == 0, sp = 1; end
v = v + 0.12*sp*(rand(size(v))-0.5);
end

function [lenFull, dur, nROIset] = draw_example_event(axK, axM, axT, axP, f, r0, t0, PCresult, P)
% ONE event, taken apart in four panels, re-measured through the same
% peak_set_measure the scoring pass uses so the numbers drawn are the ones in PkTab.
%   axK  KYMOGRAPH, basal at the bottom through the soma to apical at the top, over
%        +/-kymoHalf ms. The members of the detecting ROI's set are ticked on the left
%        edge, so the vertical extent of the ticks IS lenSet.
%   axM  DENDRITE MAP, footprints painted by set membership
%   axT  DURATION, the set's representative trace with the span above Vthr shaded
%   axP  AMPLITUDE vs POSITION, every ROI with its CI, the sets coloured and their
%        spans drawn — the panel lenSet is read off
RawE = interpolateNaN(PCresult{f}.NormalizedTrace_dirt / PCresult{f}.SpikeHeight);
BpE  = bandpass_rows(RawE, P.band, 1000, P.bpOrd);
MwE  = movmean(BpE, 2*P.avgwin+1, 2, 'omitnan');
sigE = roi_noise(MwE, [], P.baseWin, P.avgwin, 'global');
DmE  = PCresult{f}.interDendDist(PCresult{f}.sortDist, PCresult{f}.sortDist);
DmE  = max(cat(3,DmE,DmE'),[],3);
dax  = PCresult{f}.dendaxis(:);
blE  = PCresult{f}.BranchLabel(:);  blE = blE(PCresult{f}.sortDist);
AdjE = dend_adjacency(blE, abs(dax), DmE, P.adjMax);
Pset = struct('Vthr',P.Vthr, 'ciZ',P.ciZ, 'durMax',P.durMax, 'avgwin',P.avgwin, ...
    'repMode',P.repMode, 'wantPath',true);     % a viewer: draw the path lenSet measures
bE = [];
if P.fastBase
    bm = movmedian(MwE, 2*P.baseWin+1, 2, 'omitnan');  bE = bm(:,t0);
end
[lenFull, dur, nROIset, ~, S, info, A] = peak_set_measure(MwE, BpE, sigE, DmE, AdjE, ...
    r0, t0, P.baseWin, P.avgwin, Pset, bE);
kT = find([S.hasTrigger], 1);
[nR, nT] = size(RawE);
tAx = -P.kymoHalf:P.kymoHalf;
W = t0 + tAx;  inb = W >= 1 & W <= nT;

% --- kymograph, basal -> soma -> apical ---
[dsorted, dord] = sort(dax);
% the SAME local baseline the measurements use — median over t0 +/- baseWin with the
% event's own window cut out — so a colour in the kymograph and a measured amplitude
% are on one scale. A whole-recording median (what this used to be) drifts away from
% the local level and makes the two disagree.
baseK = median(RawE(:, baseline_idx(t0, nT, P.baseWin, P.avgwin)), 2, 'omitnan');
K = nan(nR, numel(tAx));
K(:, inb) = RawE(dord, W(inb)) - baseK(dord);
cla(axK);
sub = K(~isnan(K));
if isempty(sub), cax = [-1 1]; else, cax = double(prctile(sub,[2 99.5])); end
if ~(cax(2) > cax(1)), cax = cax(1) + [-1 1]; end
imagesc(axK, tAx, 1:nR, K, cax);
set(axK,'YDir','normal');  hold(axK,'on');
rowTrig = find(dord == r0, 1);
plot(axK, tAx([1 end]), rowTrig*[1 1], 'w-', 'linewidth',1.2);
[~, rowSoma] = min(abs(dsorted));
plot(axK, tAx([1 end]), rowSoma*[1 1], 'w:', 'linewidth',1);
plot(axK, [0 0], [0.5 nR+0.5], 'w--');
if ~isempty(kT)
    rowsSet = arrayfun(@(rr) find(dord == rr, 1), S(kT).roi);
    plot(axK, repmat(tAx(1)+0.06*range(tAx), size(rowsSet)), rowsSet, 'w.', 'markersize',10);
end
colormap(axK, turbo);  axis(axK,'tight');
yt = unique(round(linspace(1,nR,7)));
set(axK, 'YTick',yt, 'YTickLabel',compose('%.0f', dsorted(yt)));
xlabel(axK, 'time from peak (ms)');  ylabel(axK, 'basal <- soma -> apical (µm)');
title(axK, sprintf('f%d ROI%d fr%d — kymograph', f, r0, t0), 'FontSize',8, ...
    'FontWeight','normal');

% --- dendrite map ---
axes(axM);  cla(axM);
ftM = abs(PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist)) > 0;
coo = get_coord(PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist));
map_sets_panel(axM, ftM, coo, S, info, r0, P.labelROIs, P.setName, true);
title(axM, sprintf('dendrite map — %d set(s)', numel(S)), 'FontSize',8, 'FontWeight','normal');

% --- duration ---
cla(axT);
repR = r0;  if ~isempty(kT), repR = S(kT).repROI; end
Bb = baseline_idx(t0, nT, P.baseWin, P.avgwin);
b0 = median(MwE(repR,Bb), 'omitnan');
trR = nan(1,numel(tAx));  trR(inb) = RawE(repR, W(inb)) - baseK(repR);
trB = nan(1,numel(tAx));  trB(inb) = BpE(repR, W(inb)) - b0;
plot(axT, tAx, trR, '-', 'Color',[0.75 0.75 0.75], 'linewidth',0.5);  hold(axT,'on');
plot(axT, tAx, trB, 'k-', 'linewidth',1.2);
plot(axT, tAx([1 end]), [0 0], 'k:');
plot(axT, tAx([1 end]), P.Vthr*[1 1], 'r:', 'linewidth',1);
yl = ylim(axT);
if dur > 0
    WW = t0 + (-P.avgwin:P.avgwin);  WW = WW(WW>=1 & WW<=nT);
    vv = BpE(repR,:) - b0;
    [~, im] = max(vv(WW));  tc = WW(im);
    i1 = tc;  while i1 > 1 && vv(i1-1) >= P.Vthr, i1 = i1-1; end
    j1 = tc;  while j1 < nT && vv(j1+1) >= P.Vthr, j1 = j1+1; end
    xa = i1 - t0;  xb = j1 - t0;
    patch(axT, [xa xb xb xa], [P.Vthr P.Vthr yl(2) yl(2)], [1 0.85 0.85], ...
        'EdgeColor','none', 'FaceAlpha',0.45);
    plot(axT, [xa xb], P.Vthr*[1 1], '-', 'Color',[0.85 0.2 0.2], 'linewidth',2.5);
    text(axT, (xa+xb)/2, P.Vthr, sprintf(' %.0f ms ', dur), 'Color',[0.6 0 0], ...
        'FontSize',8, 'FontWeight','bold', 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom');
    uistack(findobj(axT,'Type','patch'), 'bottom');
end
ylim(axT, yl);  xlim(axT, tAx([1 end]));  box(axT,'off');
xlabel(axT, 'time from peak (ms)');  ylabel(axT, 'spike height');
title(axT, sprintf('durSet %.0f ms (rep ROI %d)', dur, repR), 'FontSize',8, ...
    'FontWeight','normal');

% --- amplitude vs position: where lenSet comes from ---
cla(axP);
draw_roi_sets(axP, A, sigE, S, info, dax, r0, P.Vthr, P.ciZ, '');
xlabel(axP, 'position along dendrite (µm from soma)');
title(axP, sprintf('lenSet %.0f µm over %d ROIs', lenFull, nROIset), 'FontSize',8, ...
    'FontWeight','normal');
end

function k = pick_example_CL(sel, rows, lenAll, nroiAll, PkMeta)
% Choose the peak the worked example takes apart. sel may be:
%   'median'    the peak closest to the median length, preferring one spread over
%               >=3 ROIs so the panels actually show a set being cut (default)
%   'longest' / 'shortest'   the extremes of length within the selection
%   'widest'    the one involving the most ROIs
%   [f t0]      BY IDENTITY: neuron and frame, exactly the form the browser prints as
%               "[exampleList: 7 13967]" and the status bar shows. Copy one straight
%               across.
%   [f r0 t0]   the same, pinned to a particular detecting ROI, for the rare frame
%               where two ROIs both fired
%   scalar      a row index into PkTab, if you already have one
% An identity that is not in the current selection is still shown, with a warning:
% you asked for that peak, so you get it, but the panels will then be displaying
% something the histogram beside them excluded.
k = [];
if isempty(rows), return; end
L = lenAll(rows);  nr = nroiAll(rows);

if isnumeric(sel) && numel(sel) >= 2
    f = sel(1);
    if numel(sel) >= 3, r0 = sel(2);  t0 = sel(3);  else, r0 = [];  t0 = sel(2);  end
    hit = PkMeta(:,1) == f & PkMeta(:,3) == t0;
    if ~isempty(r0), hit = hit & PkMeta(:,2) == r0; end
    k = find(hit, 1);
    if isempty(k)
        % the frame may have moved by a few ms: the peak detector's own dedup keeps
        % the largest peak of a cluster, which is not always the frame you noted
        cand = find(PkMeta(:,1) == f);
        if ~isempty(r0), cand = cand(PkMeta(cand,2) == r0); end
        if ~isempty(cand)
            [dt, i1] = min(abs(PkMeta(cand,3) - t0));
            if dt <= 50
                k = cand(i1);
                warning(['exampleSel_CL: no peak at f%d frame %d; using the nearest ' ...
                    'one in that neuron, frame %d (%d ms away).'], f, t0, PkMeta(k,3), dt);
            end
        end
    end
    if isempty(k)
        warning('exampleSel_CL: no peak matches f%d frame %d — falling back to the median.', f, t0);
    elseif ~ismember(k, rows)
        warning(['exampleSel_CL: f%d ROI%d frame%d is NOT in the filtered set the ' ...
            'histogram shows.\n         It is drawn anyway, but it is not one of the ' ...
            'peaks counted beside it.'], PkMeta(k,1), PkMeta(k,2), PkMeta(k,3));
    end
    if ~isempty(k), return; end
elseif isnumeric(sel) && isscalar(sel)
    k = sel;  return
end

switch lower(char(string(sel)))
    case 'longest',  [~, i1] = max(L);  k = rows(i1);
    case 'shortest', [~, i1] = min(L);  k = rows(i1);
    case 'widest',   [~, i1] = max(nr); k = rows(i1);
    otherwise
        cand = rows(nr >= 3);  Lc = L(nr >= 3);
        if isempty(cand), cand = rows;  Lc = L; end
        [~, i1] = min(abs(Lc - median(L)));
        k = cand(i1);
end
end

function d = divisors_near(n, want)
% The few whole divisors of n closest to `want`, for suggesting a bin width that makes
% a histogram's last bin come out full-width. Returns at most three.
n = round(n);
if n < 1, d = want; return; end
dv = find(mod(n, 1:n) == 0);
dv = dv(dv >= 2 & dv <= n);
if isempty(dv), d = want; return; end
[~, o] = sort(abs(dv - want));
d = sort(dv(o(1:min(3,numel(o)))));
end

function [lo, hi] = panel_range(v, isLengthAxis, minLen, topLen)
% Axis limits for one scatter axis: a LENGTH axis keeps the bounds the rest of the cell
% uses (so panels are comparable and the refLen line lands where it does in the
% histogram), any other parameter is clipped to its own 0.5-99.5 percentile so a few
% outliers cannot squash the cloud.
if isLengthAxis
    lo = minLen;  hi = topLen;
else
    lo = prctile(v,0.5);  hi = prctile(v,99.5);
    if ~(hi > lo), lo = min(v);  hi = max(v); end
end
if ~(hi > lo), hi = lo + 1; end
end

function browse_len_dur(PkTab, PkMeta, PCresult, sel, cdata, cfg)
% Clickable length x duration scatter that shows HOW the clicked peak was measured.
% Left: one dot per drawn peak at (lenSet, durSet), jittered by cfg.jitter so that
% peaks sharing a coordinate can still be picked apart. Clicking sorts every drawn peak
% by distance to the pointer, so the arrow keys then walk outward through the
% neighbours — that is what makes the pile at (0,0) browsable.
% Right, four panels, all recomputed for the clicked peak by calling the SAME
% event_len_dur the scoring pass used, so what is drawn is the measurement itself and
% not a redrawing of it:
%   kymograph  every ROI by distance from the soma; the ROIs whose amplitude cleared
%              Vthr — the ones LENGTH counted — are ticked on the left edge
%   trace      the set's representative ROI, filtered over raw, with the threshold and
%              the span between the crossings that IS durSet
%   set map    the morphology painted by set membership, gaps marked
%   profile    amplitude against position withevery ROI's CI and the set spans drawn —
%              length is read off this panel, so it shows the measurement directly
% Read-only: no labels are written and PkTab is not touched.
st.PkTab = PkTab;  st.PkMeta = PkMeta;  st.PCresult = PCresult;
st.sel = sel(:);  st.cfg = cfg;
rs = RandStream('twister','Seed',cfg.jitterSeed);
st.X = PkTab(st.sel,cfg.iL) + cfg.jitter(1)*2*(rand(rs,numel(st.sel),1)-0.5);
st.Y = PkTab(st.sel,cfg.iD) + cfg.jitter(2)*2*(rand(rs,numel(st.sel),1)-0.5);
st.ord = [];  st.j = 1;

fig = figure(160);  clf(fig);
set(fig, 'Name','length x duration browser — click a dot for its kymograph', ...
    'NumberTitle','off', 'Color','w', 'KeyPressFcn',@bld_key);
place_beside(fig, cfg);   % open next to the figure this browser belongs to
st.axS = axes(fig, 'Position',[0.045 0.13 0.235 0.78]);
hs = scatter(st.axS, st.X, st.Y, 7, cdata, 'filled', 'MarkerFaceAlpha',0.45);
hold(st.axS,'on');
if ~isempty(cfg.cLim), clim(st.axS, cfg.cLim); end
colormap(st.axS, turbo);
cb = colorbar(st.axS);  cb.Label.String = cfg.cLab;
% hand-labelled peaks and the hand-picked candidates, marked but not selectable
inSel = ismember(st.sel, find(cfg.labAll == 1));
hR = plot(st.axS, st.X(inSel), st.Y(inSel), 'ko', 'markersize',7, 'linewidth',1);
inSel = ismember(st.sel, find(cfg.labAll == 0));
hA = plot(st.axS, st.X(inSel), st.Y(inSel), 'o', 'markersize',7, 'linewidth',1, ...
    'Color',[0.85 0.33 0.1]);
hE = plot(st.axS, cfg.exPk(:,cfg.iL), cfg.exPk(:,cfg.iD), 'go', 'markersize',11, ...
    'linewidth',2);
st.hSel = plot(st.axS, nan, nan, 'mo', 'markersize',14, 'linewidth',2);
% every child must be click-through, or clicking a marker never reaches the axes
set([hs hR hA hE st.hSel], 'HitTest','off', 'PickableParts','none');
set(st.axS, 'ButtonDownFcn',@(s,e) bld_click(fig));
xlabel(st.axS, sprintf('Length above %.2f spike height (µm)', cfg.Vthr));
ylabel(st.axS, sprintf('Duration above %.2f spike height (ms)', cfg.Vthr));
ttlS = sprintf('%d peaks — click one  (jitter +/-%g µm, +/-%g ms)', ...
    numel(st.sel), cfg.jitter(1), cfg.jitter(2));
if isfield(cfg,'filtTxt') && ~isempty(cfg.filtTxt)
    ttlS = sprintf('%s\n%s', ttlS, strtrim(regexprep(cfg.filtTxt, '^,\s*', '')));
end
title(st.axS, ttlS, 'FontSize',9);
box(st.axS,'off');

st.axK = axes(fig, 'Position',[0.365 0.57 0.235 0.34]);   % kymograph
st.axT = axes(fig, 'Position',[0.365 0.13 0.235 0.33]);   % trace + DURATION
st.axM = axes(fig, 'Position',[0.680 0.57 0.180 0.34]);   % voltage footprint
st.axP = axes(fig, 'Position',[0.680 0.13 0.245 0.33]);   % profile + LENGTH
st.txt = uicontrol(fig, 'Style','text', 'Units','normalized', 'BackgroundColor','w', ...
    'Position',[0.02 0.005 0.96 0.055], 'HorizontalAlignment','left', 'FontSize',9, ...
    'String','Click a dot in the scatter.');
guidata(fig, st);
end

function place_beside(fig, cfg)
% Put this figure immediately to the right of the figure it belongs to, so the plane
% and the browser of that plane are on screen together. Falls back to leaving the
% position alone if that figure is not open or there is no room beside it.
if ~isfield(cfg,'nextTo') || isempty(cfg.nextTo), return; end
ref = findobj(groot, 'Type','figure', 'Number', cfg.nextTo);
if isempty(ref), return; end
try
    pr = get(ref(1), 'Position');
    pf = get(fig, 'Position');
    scr = get(groot, 'ScreenSize');
    x = pr(1) + pr(3) + 12;
    if x + pf(3) > scr(3)                 % no room on the right: go underneath
        x = pr(1);
        y = max(40, pr(2) - pf(4) - 70);
    else
        y = pr(2);
    end
    set(fig, 'Position', [x y pf(3) pf(4)]);
catch
    % window management is a convenience, never a reason to fail
end
end

function bld_click(fig)
st = guidata(fig);
bld_select(fig, st.axS.CurrentPoint(1,1:2));
end

function bld_select(fig, cp)
% Rank every drawn peak by distance to the point cp = [length duration] and select the
% nearest. Distances are normalized by the axis ranges, so a click weighs the two axes
% the way they look on screen rather than in µm vs ms.
st = guidata(fig);
rx = max(diff(xlim(st.axS)), eps);  ry = max(diff(ylim(st.axS)), eps);
d2 = ((st.X - cp(1))/rx).^2 + ((st.Y - cp(2))/ry).^2;
[~, st.ord] = sort(d2);
st.j = 1;
guidata(fig, st);
bld_draw(fig, true);
end

function bld_key(fig, ev)
st = guidata(fig);
if isempty(st.ord), return; end
switch lower(ev.Key)
    case 'rightarrow', st.j = min(st.j+1, numel(st.ord));
    case 'leftarrow',  st.j = max(st.j-1, 1);
    case 'uparrow',    st.j = min(st.j+10, numel(st.ord));
    case 'downarrow',  st.j = max(st.j-10, 1);
    case 'p',          guidata(fig, st);  bld_draw(fig, true);  return
    case 'c'
        k = st.sel(st.ord(st.j));
        s = sprintf('%d %d', st.PkMeta(k,1), st.PkMeta(k,3));
        try, clipboard('copy', s); fprintf('  copied "[%s]" for exampleList\n', s);
        catch, fprintf('  (clipboard unavailable) [%s]\n', s); end
        return
    case 'x'
        % hand this peak to the worked example: set exampleSel_CL in the base
        % workspace so re-running the length cell takes THIS peak apart. Nothing is
        % redrawn here — the figure-154 panels are built by that cell, not this one.
        k = st.sel(st.ord(st.j));
        id = st.PkMeta(k,[1 2 3]);
        assignin('base', 'exampleSel_CL', id);
        fprintf(['  exampleSel_CL = [%d %d %d]  (f%d ROI%d frame%d).\n' ...
            '  Re-run the length-distribution cell to show this peak in figure 154.\n'], ...
            id(1), id(2), id(3), id(1), id(2), id(3));
        return
    otherwise, return
end
guidata(fig, st);
bld_draw(fig, true);
end

function bld_draw(fig, doPrint)
st = guidata(fig);
if isempty(st.ord), return; end
q = st.ord(st.j);  k = st.sel(q);
f = st.PkMeta(k,1);  r0 = st.PkMeta(k,2);  t0 = st.PkMeta(k,3);
v = st.PkTab(k,:);  P = st.cfg;
set(st.hSel, 'XData', st.X(q), 'YData', st.Y(q));

C = peak_cache(st.PCresult, f, P);
[nR, nT] = size(C.Raw);

% --- REDO the measurement with the SAME wrapper the scoring pass used, so the panels
% --- below are the measurement itself rather than an illustration of it
Pset = struct('Vthr',P.Vthr, 'ciZ',P.ciZ, 'durMax',P.durMax, 'avgwin',P.avgwin, ...
    'repMode',P.repMode, 'wantPath',true);     % a viewer: draw the path lenSet measures
bIn = [];
if ~isempty(C.baseMw), bIn = C.baseMw(:,t0); end
[lenN, durN, nRs, nSe, S, info, A] = peak_set_measure(C.Mw, C.Bp, C.sig, C.Dmat, ...
    C.Adj, r0, t0, P.baseWin, P.avgwin, Pset, bIn);
kTrig = find([S.hasTrigger], 1);

% --- kymograph: every ROI, ordered by distance from the soma ---
W = t0 + P.tAxis;  inb = W >= 1 & W <= nT;
K = nan(nR, numel(P.tAxis));
baseLoc = median(C.Raw(:, baseline_idx(t0, nT, P.baseWin, P.avgwin)), 2, 'omitnan');
K(:, inb) = C.Raw(C.dord, W(inb)) - baseLoc(C.dord);   % local baseline, as measured
rowTrig = find(C.dord == r0, 1);
[~, rowSoma] = min(abs(C.dsorted));
d = C.dsorted(:);
cla(st.axK);
sub = K(~isnan(K));
if isempty(sub), cax = [-1 1]; else, cax = double(prctile(sub,[2 99.5])); end
if ~(cax(2) > cax(1)), cax = cax(1) + [-1 1]; end
imagesc(st.axK, P.tAxis, 1:nR, K, cax);
set(st.axK, 'YDir','normal');  hold(st.axK,'on');
plot(st.axK, P.tAxis([1 end]), rowTrig*[1 1], 'w-', 'linewidth',1.5);
plot(st.axK, P.tAxis(1)+0.03*range(P.tAxis), rowTrig, 'w>', ...
    'MarkerFaceColor','w', 'markersize',9);
if ~isnan(rowSoma)
    plot(st.axK, P.tAxis([1 end]), rowSoma*[1 1], 'w:', 'linewidth',1);
end
plot(st.axK, [0 0], [0.5 nR+0.5], 'w--');
% the members of the detecting ROI's SET — the ROIs lenSet was measured across
if ~isempty(kTrig)
    rowsSet = arrayfun(@(rr) find(C.dord == rr, 1), S(kTrig).roi);
    plot(st.axK, repmat(P.tAxis(1)+0.09*range(P.tAxis), size(rowsSet)), rowsSet, ...
        'w.', 'markersize',9);
end
if ~isempty(P.aucWin)                          % what aucSum integrated over
    plot(st.axK, P.aucWin(1)*[1 1], [0.5 nR+0.5], 'w-', 'linewidth',0.5);
    plot(st.axK, P.aucWin(2)*[1 1], [0.5 nR+0.5], 'w-', 'linewidth',0.5);
end
colormap(st.axK, turbo);  axis(st.axK,'tight');
yt = unique(round(linspace(1,nR,7)));
set(st.axK, 'YTick',yt, 'YTickLabel',compose('%.0f', d(yt)));
ylabel(st.axK, 'distance from soma (\mum)');
title(st.axK, sprintf('f%d  ROI %d (%.0f \\mum, row %d/%d)  frame %d   [%d/%d]', ...
    f, r0, d(rowTrig), rowTrig, nR, t0, st.j, numel(st.ord)), 'FontSize',9);

% --- trace: the SET'S representative ROI, and the span that is durSet ---
cla(st.axT);
repR = r0;
if ~isempty(kTrig), repR = S(kTrig).repROI; end
base0 = mean(C.Mw(repR, max(1,t0-P.baseWin):min(nT,t0+P.baseWin)), 'omitnan');
trRaw = nan(1, numel(P.tAxis));  trRaw(inb) = C.Raw(repR, W(inb)) - baseLoc(repR);
trBp  = nan(1, numel(P.tAxis));  trBp(inb)  = C.Bp(repR, W(inb)) - base0;
plot(st.axT, P.tAxis, trRaw, '-', 'Color',[0.72 0.72 0.72], 'linewidth',0.5);
hold(st.axT,'on');
plot(st.axT, P.tAxis, trBp, 'k-', 'linewidth',1.2);
plot(st.axT, P.tAxis([1 end]), [0 0], 'k:');
plot(st.axT, P.tAxis([1 end]), P.Vthr*[1 1], 'r:', 'linewidth',1);
yl = ylim(st.axT);
plot(st.axT, [0 0], yl, 'r--');
if durN > 0
    % the walk starts at the trace's own maximum within +/-avgwin of t0
    WW = t0 + (-P.avgwin:P.avgwin);  WW = WW(WW>=1 & WW<=nT);
    [~, im] = max(C.Bp(repR,WW) - base0);  tc = WW(im);
    xa = (tc - t0) - 0;  vv = C.Bp(repR,:) - base0;
    i = tc;  while i > 1 && vv(i-1) >= P.Vthr, i = i-1; end
    j2 = tc; while j2 < nT && vv(j2+1) >= P.Vthr, j2 = j2+1; end
    xa = i - t0;  xb = j2 - t0;
    patch(st.axT, [xa xb xb xa], [P.Vthr P.Vthr yl(2) yl(2)], [1 0.85 0.85], ...
        'EdgeColor','none', 'FaceAlpha',0.45);
    plot(st.axT, [xa xa], yl, '-', 'Color',[0.85 0.2 0.2], 'linewidth',1);
    plot(st.axT, [xb xb], yl, '-', 'Color',[0.85 0.2 0.2], 'linewidth',1);
    plot(st.axT, [xa xb], P.Vthr*[1 1], '-', 'Color',[0.85 0.2 0.2], 'linewidth',2.5);
    text(st.axT, (xa+xb)/2, P.Vthr, sprintf(' %.0f ms ', durN), 'Color',[0.6 0 0], ...
        'FontSize',9, 'FontWeight','bold', 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom');
    uistack(findobj(st.axT,'Type','patch'), 'bottom');
end
ylim(st.axT, yl);  xlim(st.axT, P.tAxis([1 end]));  box(st.axT,'off');
xlabel(st.axT, 'Time from peak (ms)');
ylabel(st.axT, sprintf('set rep ROI %d (spike height)', repR));
title(st.axT, sprintf('DURATION %.0f ms   (grey = unfiltered, black = %s)', durN, ...
    P.bpTxtShort), 'FontSize',9, 'FontWeight','normal');

% --- the sets on the morphology ---
axes(st.axM);  cla(st.axM);
map_sets_panel(st.axM, C.ftMask, C.coords, S, info, r0, P.labelROIs, P.setName, true);
title(st.axM, sprintf('%d set(s), %d gap merged', numel(S), sum(info.gaps(:,4))), ...
    'FontSize',9, 'FontWeight','normal');

% --- the profile the sets were cut from: LENGTH is read off here ---
cla(st.axP);
draw_roi_sets(st.axP, A, C.sig, S, info, C.dendRaw, r0, P.Vthr, P.ciZ, '');
xlabel(st.axP, 'position along dendrite (\mum)');
title(st.axP, sprintf('LENGTH %.0f \\mum (full), %d ROIs in set', lenN, nRs), ...
    'FontSize',9, 'FontWeight','normal');

% --- every parameter, plus the classifier score and your label if they exist ---
txt = strjoin(compose('%s=%.3g', string(P.pkNames(:)), v(:))', '   ');
extra = '';
if ~isnan(P.pAll(k)), extra = sprintf('   p(real)=%.3f', P.pAll(k)); end
if ~isnan(P.labAll(k))
    lbl = {'ARTIFACT','','REAL'};
    if P.labAll(k) == -1, extra = [extra '   your label: UNSURE'];
    else, extra = [extra '   your label: ' lbl{P.labAll(k)+2}];
    end
end
% If the recomputation disagrees with the stored value, the scoring parameters have
% been edited since PkTab was built. Say so rather than showing two different numbers
% in two panels and letting them be read as one result.
reTxt = '';
if abs(lenN - v(P.iL)) > 1e-6 || abs(durN - v(P.iD)) > 1e-6
    reTxt = sprintf(['   ** redrawn len %.0f / dur %.0f DIFFER from PkTab — the ' ...
        'scoring parameters changed, re-run the scoring cell **'], lenN, durN);
end
set(st.txt, 'String', sprintf('f%d ROI%d frame%d%s%s\n%s', f, r0, t0, extra, reTxt, txt));
if doPrint
    fprintf('  f%-3d ROI%-3d frame%-8d %s%s%s   [exampleList: %d %d]\n', f, r0, t0, ...
        txt, extra, reTxt, f, t0);
end
guidata(fig, st);
end

function C = peak_cache(PCresult, f, P)
% Everything the browser needs about ONE neuron, built on demand and kept for the last
% few neurons visited: normalising a whole recording and bandpassing it is the slow
% part, and browsing tends to stay within one neuron.
%   Raw  : normalized traces, spike-height units, sortDist order (what the kymograph
%          and the grey trace show)
%   Bp   : the SAME traces bandpassed exactly as the scoring pass did — this is what
%          length and duration are measured on
%   Sub  : subthreshold traces (event_len_dur uses them to find the peak frame)
%   base : per-ROI median of Raw over the whole recording, for display only
%   Dmat : symmetric contour-distance matrix
%   ftMask : footprint mask per ROI for the voltage-footprint panel. Only the LOGICAL
%          mask is kept — it is what showScaleImageSmooth wants and it is ~8x smaller
%          than the weighted stack, which matters with three neurons cached. The
%          trigger ROI's outline is traced per click instead of all of them per neuron.
persistent ck cv
if isempty(ck), ck = [];  cv = {}; end
h = find(ck == f, 1);
if ~isempty(h)
    C = cv{h};                                            % move to the front (LRU)
    ck = [f ck([1:h-1 h+1:end])];  cv = [{C} cv([1:h-1 h+1:end])];
    return
end
C.Raw = interpolateNaN(PCresult{f}.NormalizedTrace_dirt / PCresult{f}.SpikeHeight);
C.Sub = interpolateNaN(PCresult{f}.Subthreshold      / PCresult{f}.SpikeHeight);
C.Bp  = bandpass_rows(C.Raw, P.band, 1000, P.bpOrd);
[C.dsorted, C.dord] = sort(PCresult{f}.dendaxis);         % basal -> soma -> apical
C.base = median(C.Raw, 2, 'omitnan');   % whole-recording level, kept only as a fallback;
                                        % the panels use a LOCAL baseline per peak
C.Dmat = PCresult{f}.interDendDist(PCresult{f}.sortDist, PCresult{f}.sortDist);
C.Dmat = max(cat(3,C.Dmat,C.Dmat'),[],3);
C.pixelsize = PCresult{f}.pixelsize;
C.ftMask = abs(PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist)) > 0;
C.coords = get_coord(PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist));
C.dendRaw = PCresult{f}.dendaxis(:);                  % signed position, sortDist order
% what the ROI-set measurement needs, built once per neuron
C.Mw  = movmean(C.Bp, 2*P.avgwin+1, 2, 'omitnan');
C.sig = roi_noise(C.Mw, [], P.baseWin, P.avgwin, 'global');
% the viewer must use the SAME baseline convention the scoring pass used, or its
% recomputation disagrees by a frame or two and the "differs from PkTab" warning fires
% on every peak for no reason
C.baseMw = [];
if isfield(P,'fastBase') && P.fastBase
    C.baseMw = movmedian(C.Mw, 2*P.baseWin+1, 2, 'omitnan');
end
bl    = PCresult{f}.BranchLabel(:);  bl = bl(PCresult{f}.sortDist);
C.Adj = dend_adjacency(bl, abs(C.dendRaw), C.Dmat, P.adjMax);
ck = [f ck];  cv = [{C} cv];
if numel(ck) > 3, ck(4:end) = [];  cv(4:end) = []; end
end

function [sampleIdx, pkLabel, srcTxt] = load_peak_labels(files, PkMeta)
% Read hand labels back by PEAK IDENTITY (file, ROI, frame) rather than by row number.
% Row numbers into PkTab are meaningless across runs — scoring one neuron produces a
% completely different PkTab from scoring twenty — so the label files store labMeta,
% the [f r0 t0] of every labelled peak, and this function matches on that instead.
% Consequences: a test run on one neuron automatically inherits whatever you already
% labelled for that neuron in a full run; a full run can inherit a test run's labels;
% and re-scoring with changed parameters keeps the labels, since peak identity is
% unchanged.
%   files  : cell array of label files, highest priority FIRST. The first file to
%            supply a label for a given peak wins; later files only fill gaps.
%   PkMeta : nPk x 3 [f r0 t0] of the run being labelled now.
% Returns the rows of PkMeta that carry a label, and those labels (1/0/-1).
if ischar(files) || isstring(files), files = {char(files)}; end
files = files(~cellfun(@(x) isempty(x) || ~(ischar(x) || isstring(x)), files));
lblByRow = nan(size(PkMeta,1), 1);
srcTxt = '';  parts = {};
for k = 1:numel(files)
    if ~isfile(files{k}), continue; end
    Lp = load(files{k});
    if ~isfield(Lp,'labMeta') || ~isfield(Lp,'pkLabel')
        warning(['%s predates identity-keyed labels (no labMeta) and cannot be ' ...
            'matched to this run — skipped.'], files{k});
        continue
    end
    [tf, loc] = ismember(PkMeta, Lp.labMeta, 'rows');
    v = nan(size(PkMeta,1),1);
    v(tf) = Lp.pkLabel(loc(tf));
    take = ~isnan(v) & isnan(lblByRow);        % first file wins, later ones fill gaps
    lblByRow(take) = v(take);
    if any(take)
        [~, nm, ex] = fileparts(files{k});
        parts{end+1} = sprintf('%s%s (%d of its %d labels apply here)', nm, ex, ...
            sum(take), sum(~isnan(Lp.pkLabel))); %#ok<AGROW>
    end
end
srcTxt = strjoin(parts, ' + ');
sampleIdx = find(~isnan(lblByRow));
pkLabel   = lblByRow(sampleIdx);
end

function [auc, n, cut, acc] = param_auc(v, isReal, isArt)
% Mann-Whitney AUC of one parameter against the hand labels, plus the cutoff that
% maximises Youden's J and the balanced accuracy there. AUC > 0.5 means "real" peaks
% score HIGHER; < 0.5 means they score LOWER; 0.5 means the parameter is useless.
% Nothing here assumes a direction or a threshold in advance.
v = v(:);  isReal = isReal(:);  isArt = isArt(:);
auc = NaN;  n = 0;  cut = NaN;  acc = NaN;
g1 = isReal & ~isnan(v);  g0 = isArt & ~isnan(v);
n1 = sum(g1);  n0 = sum(g0);  n = n1 + n0;
if n1 < 3 || n0 < 3, return; end
rk = tiedrank([v(g1); v(g0)]);
auc = (sum(rk(1:n1)) - n1*(n1+1)/2) / (n1*n0);
s = sign(auc - 0.5);  if s == 0, s = 1; end        % orient so "real" is high
thr = unique(v(g1|g0));
J = nan(numel(thr),1);
for q = 1:numel(thr)
    pred = s*v >= s*thr(q);
    J(q) = mean(pred(g1)) - mean(pred(g0));
end
[~, jb] = max(J);  cut = thr(jb);
pred = s*v >= s*cut;
acc = 0.5*(mean(pred(g1)) + 1 - mean(pred(g0)));   % balanced accuracy
end

function [len, FW] = movie_len(rec, ampTrig, Vthr)
% Threshold-crossing length of an example event measured from the MOVIE, using the
% polyline kymograph record saved in ExampleL50_movie.mat (rec = MovieL50(i).rec).
% rec.profileNorm is the spatial profile normalized to its peak (∆F, a.u.), so it
% is scaled by ampTrig (the event's ROI-measured peak amplitude, spike-height) to make the
% absolute threshold Vthr meaningful. Crossings are found on both sides of the peak
% along the polyline; FW = full width, len = FW/2 (same convention as L50_half, so
% it is comparable to a crossing radius measured from the trigger).
len = NaN;  FW = NaN;
if isempty(rec) || ~isfield(rec,'profileNorm') || ~isfield(rec,'arclen_um'), return; end
s = rec.arclen_um(:);  V = rec.profileNorm(:) * ampTrig;
if numel(s) < 3 || ~(max(V) >= Vthr), return; end
[~, pk] = max(V);
sL = cross_level(s, V, pk, -1, Vthr);
sR = cross_level(s, V, pk, +1, Vthr);
FW = sR - sL;  len = FW/2;
end

function sX = cross_level(s, A, pk, dir, lev)
% Position along arc-length s of the LAST sample still at or above the level lev,
% searching from the peak index pk in direction dir (+1 right, -1 left). Sampled
% positions only — no interpolation between polyline boxes. If A never drops below
% lev, the corresponding polyline endpoint is returned (censored).
% cross_half() is the interpolated lev = 0.5 version used by the L50 analysis.
n = numel(A);  i = pk;
if dir < 0
    while i > 1 && A(i-1) >= lev, i = i-1; end
else
    while i < n && A(i+1) >= lev, i = i+1; end
end
sX = s(i);
end

function lbl = viewLbl_one(ViewTab, k)
% One-line identity of event k for the kymograph titles and the interactive viewer.
lbl = sprintf('f%d  ROI %d  frame %d  |  len %.0f µm  dur %.0f ms%s', ...
    ViewTab(k,1), ViewTab(k,2), ViewTab(k,3), ViewTab(k,5), ViewTab(k,6), ...
    repmat('  [CANDIDATE]', 1, ViewTab(k,9)==1));
end

function r0 = movie_trigger_roi(rec, ftprnts)
% Trigger ROI defined by the MOVIE: the peak of the polyline profile in
% MovieL50(i).rec is a pixel on the drawn line, and the returned ROI is the one
% whose footprint carries the most weight there (nearest footprint centre of mass
% if no footprint covers the pixel). NaN if the record has no polyline.
r0 = NaN;
if isempty(rec) || ~isfield(rec,'polyline') || ~isfield(rec,'profileNorm'), return; end
pos = rec.polyline;                       % [X Y] pixel centre of each kymograph box
if size(pos,1) < 2, return; end
[~, pk] = max(rec.profileNorm(:));
pk = min(max(pk,1), size(pos,1));
xy = round(pos(pk,:));                    % pixel of the profile peak
[H, W, nR] = size(ftprnts);
xy(1) = min(max(xy(1),1), W);  xy(2) = min(max(xy(2),1), H);
wgt = squeeze(ftprnts(xy(2), xy(1), :));  % footprint weight of every ROI at that pixel
if any(wgt > 0)
    [~, r0] = max(wgt);  return
end
cen = nan(nR,2);                          % no footprint there -> nearest centre of mass
for r = 1:nR
    f = ftprnts(:,:,r);  tot = sum(f(:));
    if tot <= 0, continue; end
    [yy, xx] = ndgrid(1:H, 1:W);
    cen(r,:) = [sum(xx(:).*f(:)) sum(yy(:).*f(:))] / tot;
end
[~, r0] = min(sum((cen - xy).^2, 2));
end

function s = repmat_str(tf, sTrue, sFalse)
% Pick one of two strings (keeps the fprintf calls on one line).
if tf, s = sTrue; else, s = sFalse; end
end

function Xs = smooth_subthreshold(X, spk, spec)
% Spike-removed, interpolated subthreshold trace (get_subthreshold, +/-3 frames),
% then smoothed by spec.type: 'ma' = movmean(spec.len), 'sg' = sgolayfilt(spec.ord, spec.len).
% sgolayfilt cannot skip NaN, so gaps are linearly filled first; callers re-mask NaN.
[~, Xu] = get_subthreshold(X, spk, 7, 1);
switch spec.type
    case 'ma'
        Xs = movmean(Xu, spec.len, 2, 'omitnan');
    case 'sg'
        Xf = fillmissing(Xu, 'linear', 2, 'EndValues', 'nearest');
        Xs = sgolayfilt(Xf, spec.ord, spec.len, [], 2);
    otherwise
        error('smooth_subthreshold:type', 'unknown smoothing type ''%s''', spec.type);
end
end

function pcmode_click(ax, k)
% Figure 219: left click toggles GLOBAL, right click toggles SEE-SAW for PC k.
fig = ancestor(ax, 'figure');
st  = fig.UserData;
switch fig.SelectionType
    case 'normal', want = 1;
    case 'alt',    want = 2;
    otherwise,     return;                  % double click etc.
end
if st.sel(k) == want, st.sel(k) = 0; else, st.sel(k) = want; end
fig.UserData = st;
pcmode_redraw(fig);
end

function pcmode_redraw(fig)
% Figure 219: frame colour per PC tile from fig.UserData.sel (0 none, 1 global, 2 see-saw).
col = {'none', [1 0 0], [0 0.3 1]};
axs = findobj(fig, 'Type', 'axes', 'Tag', 'pcmode');
sel = fig.UserData.sel;
tcol = {get(groot, 'defaultTextColor'), [1 0 0], [0 0.3 1]};
for a = axs(:)'
    k = sscanf(a.Title.String, 'PC%d');
    a.UserData.EdgeColor = col{sel(k) + 1};
    a.Title.Color = tcol{sel(k) + 1};
end
end

function pcmode_key(fig, ev)
% Figure 219 keys: Enter/n/-> next, b/<- back, c clear, q/Esc save & stop.
switch ev.Key
    case {'return', 'n', 'rightarrow'}, act = 'next';
    case {'b', 'leftarrow'},            act = 'back';
    case {'q', 'escape'},               act = 'quit';
    case 'c'
        fig.UserData.sel(:) = 0;  pcmode_redraw(fig);  return;
    otherwise, return;
end
fig.UserData.action = act;
uiresume(fig);
end

function pcmode_close(fig)
% Figure 219: closing the window = save & stop (the figure itself stays until the loop ends).
fig.UserData.action = 'quit';
uiresume(fig);
end
