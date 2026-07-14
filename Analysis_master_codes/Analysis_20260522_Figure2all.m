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
save(fullfile(save_figto, 'Correlation_Fig2_results.mat'), ...
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


%% ===== Figure 2D — spatial correlation matrix per session =====
figure(118); clf; tiledlayout('flow');
cmap_label = hsv(5);  cmap_label = cmap_label([4 2 3 5 1],:);
for f = foi
    nexttile([1 1]);
    imshow_label(corrMat{f,1}, PCresult{f}.labelvec, cmap_label, classNames);
    axis equal tight off
    caxis([-0.3 1]);  title(num2str(f));
    if f ~= foi(end), legend off; end
    drawnow;
end
colormap(gen_colormap([0 0.4 1; 1 1 1; 1 0 0], 256));

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

figure(127); clf;
%posL = lags >= 0;                         % positive lag side
imagesc(lags(:), binRes.centers, heat(:, :)); axis xy;
colormap(gen_colormap([0 0.4 1; 1 1 1; 1 0 0], 256)); caxis([-0.05 0.25]);
cb = colorbar; %cb.Label.String = 'Correlation';
xlabel('Lag (ms)'); ylabel('Distance from triggering ROI (µm)');
set_font('Arial'); set_fontsize(15);
set_figsize(115,90);
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
            [pks, locs] = findpeaks(FilterTr(n,:),'MinPeakHeight', 0.2,'MinPeakDistance',50, ...
                'MinPeakProminence', 0.1);  % Prominence can also be tuned
            [trough, loc_trgh] = findpeaks(-FilterTr(n,:),'MinPeakHeight', 0.2,'MinPeakDistance',50, ...
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
% %%
% figure(135); clf; tiledlayout(1,2); cmap_temp=gen_colormap(Plasma,4);
% ddd={[1 4],[3 2]};
% title_str={'Basal peak-TA','Apical peak-TA'};
% for dc2=1:length(ddd)
%     ldc=[]; g=1;
%     for dc=ddd{dc2}
%         nexttile(dc2,[1 1]);
%         PeakTrace2show=cell2mat(cellfun(@(x) x(:,dc),PeakTraceROI(foi),'UniformOutput',false));
%         %plot([-nTauPeak(1):nTauPeak(2)],PeakTrace2show,'color',[0.7 0.7 0.7]); hold all;
%         ldc(g)=errorbar_shade([-nTauPeak(1):nTauPeak(2)],mean(PeakTrace2show,2,'omitnan'),std(PeakTrace2show,0,2,'omitnan'),cmap_errorbar(dc,:));
%         xlabel('Time (ms)'); ylabel('Subthreshold (spike height)');
%         title(title_str{dc2})
%         g=g+1;
%     end
%     legend(ldc,{'Basal','Apical'})
% end


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