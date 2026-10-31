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


%% ===== Figure 2D — subthreshold correlation matrix (ch1 x ch2), leaf-ordered + dendrogram =====
% Cross-channel (ch1 x ch2) subthreshold correlation matrix for one neuron,
% computed here from the 20-ms subthreshold traces (get_subthreshold) over the
% SILENT period (blue-off AND not peri-spike); the ch1 x ch2 diagonal is
% shot-noise-free. Same calculation as the reviewer thirds cell, but on silent
% frames. Reordered by hierarchical-clustering leaf order, with the dendrogram on
% the left and dendrite-class colour dots on the top and left.
f = 20;
Pr  = PCresult{f};
ch1 = Pr.NormalizedTrace_ch{1};  ch2 = Pr.NormalizedTrace_ch{2};  trC = Pr.NormalizedTrace_dirt;
spk = Pr.SpikeMat(1,:) > 0;
win20 = 20;  minPairs = 100;
goodvec = ~imdilate(Pr.BlueStim>0, ones(1,201));
% silent period = blue-off AND not peri-spike ( +/- [-3:10] around spikes, + CS )
perispk = unique([tovec(find(spk)' + (-3:10)); find(Pr.CStrace)']);
perispk(perispk < 1 | perispk > numel(spk)) = [];
silentvec = true(1, numel(spk));  silentvec(perispk) = false;
silentMask = goodvec & silentvec;
sub1 = Pr.Subthreshold_ch{1};  sub1(isnan(ch1)) = NaN;
sub2 = Pr.Subthreshold_ch{2};  sub2(isnan(ch2)) = NaN;
subC = Pr.Subthreshold;  subC(isnan(trC)) = NaN;
%cm = get_corrTensor(sub1, sub2, subC, silentMask, 0, minPairs);   % zero-lag cross-channel matrix (silent)
C = get_corrMat(sub1, sub2, find(silentMask));
%C = cm(:,:,1);
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

%% ===== Local vs amplitude-matched events — L50 spatial extent (Supp. Fig) =====
% (1) Among the nominated local candidates (exampleList), 0.75 x their mean peak
%     amplitude sets an amplitude threshold. Amplitudes are in SPIKE-HEIGHT units
%     (subthreshold trace divided by the somatic spike height).
% (2) Every detected subthreshold peak (peakvec) with peak amplitude >= threshold,
%     across the whole recording, is kept as an event.
% (3) Each event's spatial extent is L50 = the contour radius at which the
%     normalized amplitude vs distance profile first crosses 0.5 (half-decay
%     distance from the trigger). Small L50 -> local; large L50 -> broad/global
%     (bounded at Dmax_ext).
% Voltage comes from the raw subthreshold trace (spike-height units); the trigger
% ROI is the max-window-voltage ROI at the peak. Helper functions (trigger_roi /
% event_extent / get_profiles / draw_event_footprint) are at the end of this file.
avgwin_ext = 10;      % +/- averaging window around each peak (ms == frames at 1 kHz)
Dmax_ext   = 300;     % contour distance range included / profiles (um); L50 <= Dmax
ampTol_ext = 0.25;    % amplitude-match tolerance for controls (+/- fraction)
nCtrl_ext  = 100;      % max amplitude-matched controls per candidate
exampleList = [26 493523; 26 539246; 26 562118; 26 576267; 18 422263];   % [fileInd  frame] local candidates
% --- (1) candidate extents + amplitude threshold ---
nEx = size(exampleList,1);
exR0=nan(nEx,1); exAMP=nan(nEx,1); exL50=nan(nEx,1);
exCOMP=nan(nEx,1); exProf=cell(nEx,1);
for i = 1:nEx
    fe = exampleList(i,1);  te = exampleList(i,2);
    SubNorm = interpolateNaN(PCresult{fe}.Subthreshold / PCresult{fe}.SpikeHeight);
    Dmat = PCresult{fe}.interDendDist(PCresult{fe}.sortDist, PCresult{fe}.sortDist);
    Dmat = max(cat(3,Dmat,Dmat'),[],3);
    r0e = trigger_roi(SubNorm, te, avgwin_ext);
    [exAMP(i), exL50(i), exProf{i}] = event_extent(SubNorm, Dmat, r0e, te, avgwin_ext, Dmax_ext);
    exR0(i) = r0e;
    lab = PCresult{fe}.labelvec(r0e);
    if ~isnan(lab), exCOMP(i) = 1 + (lab>=3); end     % 1 basal, 2 apical-ish (labels 3/4/5)
    fprintf('Candidate %d: f=%d frame=%d ROI=%d amp=%.3f spike-height = %.4f dF/F  L50=%.0f um\n', ...
        i, fe, te, r0e, exAMP(i), exAMP(i)*PCresult{fe}.SpikeHeight, exL50(i));
end
ampThresh = 0.75 * mean(exAMP);      % 0.75 x mean candidate amplitude (spike-height units)
% same threshold in dF/F for the candidate neuron(s) (spike-height x SpikeHeight)
ampThresh_dff = 0.75 * mean(arrayfun(@(i) exAMP(i)*PCresult{exampleList(i,1)}.SpikeHeight, 1:nEx));
fprintf('Amplitude threshold (0.75 x mean of %d candidates) = %.3f spike-height = %.4f dF/F (cand neuron)\n', ...
    nEx, ampThresh, ampThresh_dff);

% --- (2) detect all events above threshold; compute extent ---
% EventTab columns: [f  r0  t0  amp  L50  comp]
EventTab = [];
for f = foi
    SubNorm = interpolateNaN(PCresult{f}.Subthreshold / PCresult{f}.SpikeHeight);
    Dmat = PCresult{f}.interDendDist(PCresult{f}.sortDist, PCresult{f}.sortDist);
    Dmat = max(cat(3,Dmat,Dmat'),[],3);
    lab  = PCresult{f}.labelvec(:);
    [r0all, t0all] = find(peakvec{f} > 0);
    for e = 1:numel(r0all)
        [a, L50] = event_extent(SubNorm, Dmat, r0all(e), t0all(e), avgwin_ext, Dmax_ext);
        if a >= ampThresh
            l = lab(r0all(e));  comp = NaN;  if ~isnan(l), comp = 1 + (l>=3); end
            EventTab = [EventTab; f r0all(e) t0all(e) a L50 comp]; %#ok<AGROW>
        end
    end
    fprintf('  file %d done (%d events so far)\n', f, size(EventTab,1));
end
fprintf('Total events above threshold: %d   |   local candidates: %d\n', ...
    size(EventTab,1), nEx);

%% ===== Suprathreshold LENGTH & DURATION of amplitude-matched events =====
% Alternative to L50: the (peak ROI, peak frame) is the trigger point, and the
% event is measured by where the voltage crosses BELOW a fixed absolute threshold
% Vthr_sh (in SPIKE-HEIGHT units, i.e. subthreshold trace / somatic spike height —
% the same units as EventTab(:,4); 0.1 = 10% of a spike), walking away from that
% trigger point along each axis:
%   duration = along TIME, from the trigger ROI's RAW trace: walk back to the last
%              frame still above threshold (tPre, e.g. -100 ms) and forward to the
%              last one (tPost, e.g. +200 ms); duration = tPost - tPre (e.g. 300 ms)
%   length   = along SPACE, from the RAW trace averaged over +/- avgwin_ext ms at the
%              peak: the contour distance of the FARTHEST ROI still above Vthr_sh
%              before the profile drops below it. Bounded to Dmax_ext.
% NO INTERPOLATION: both values are reported at the sampled grid (ROI distances,
% 1 kHz frames). A length of 0 means only the trigger ROI was above threshold — the
% ROI spacing is the resolution limit, so anything shorter is simply unresolved.
% WHICH TRACE: the peak (trigger ROI r0, frame t0) is located on the SUBTHRESHOLD
% trace as before (peakvec / EventTab), but both measurements are made on the RAW
% trace (NormalizedTrace_dirt / SpikeHeight). The subthreshold filter broadens an
% event in time, so raw-trace crossings give the shorter, less filtered duration.
% Inside the +/- avgwin_ext window around t0 the raw maximum is used as the seed of
% the crossing walk (the event is still anchored to the subthreshold peak).
% NOTE: the raw trace contains spikes — a spike inside the event window keeps the
% trace above threshold and lengthens the measured duration.
% If no crossing is found the value is censored (time -> durMax_LD cap, space ->
% outermost sampled ROI / Dmax_ext) and
% counted in the printout. Amplitudes are baseline-subtracted per ROI with a local
% median (the event window itself is excluded from the baseline). Events are
% amplitude-matched to the local candidates: only EventTab rows within
% ±ampTol_ext of the candidates' mean peak amplitude are kept.
% Blood-vessel / motion artifacts show up as events that are both spatially and
% temporally tiny (one ROI, a few ms), so events below lenMin_LD or durMin_LD are
% rejected as noise (same idea as L50min in the L50 analysis) — the rejected counts
% are printed and the cutoffs are drawn in panel A. The cutoffs apply ONLY to the
% automatically detected control pool; the hand-picked candidates are never
% rejected by them (they are known-good events, and dropping one would bias the
% comparison towards whatever the cutoff selects for).
% Helper event_len_dur() is at the end of this file.
Vthr_sh    = 0.15;     % voltage threshold in SPIKE-HEIGHT units (0.10 = 10% of a spike)
baseWin_LD = 1000;      % +/- ms local baseline window per ROI (event window excluded)
durMax_LD  = 500;      % +/- ms cap on the duration search
smoothDur_LD = 3;      % movmean (ms) applied to the RAW trace before the duration search (1 = none)
lenMin_LD  = 25;        % NOISE CUTOFF (controls only): length below this (µm) -> blood vessel
durMin_LD  = 0;       % NOISE CUTOFF (controls only): duration below this (ms) -> noise
nMaxPerN_LD   = 500;    % BALANCE: max matched events per neuron (a few neurons dominate)
minEvtPerN_LD = 0;     % BALANCE: neurons with fewer kept events are dropped from the stats
drawReject_LD  = true;
% --- PEAK AMPLITUDE RANGE of the events to analyse (spike-height units) ---
% [] -> derived from the candidates: mean(exAMP) * (1 +/- ampTol_ext).
% Otherwise set an explicit [lo hi] range here, e.g. ampRange_LD = [0.15 0.40].
% EventTab only holds events above ampThresh = 0.75*mean(exAMP), so a lower bound
% below that cannot add events back (a warning is printed if you ask for one).
ampRange_LD = [];
if isempty(ampRange_LD)
    ampBand_LD = [1-ampTol_ext 1+ampTol_ext] * mean(exAMP);
else
    ampBand_LD = ampRange_LD(:)';
end
if ampBand_LD(1) < ampThresh
    warning(['Requested lower amplitude bound %.3f is below the detection threshold ' ...
        'ampThresh=%.3f — no events exist below it.'], ampBand_LD(1), ampThresh);
end
if isempty(ampRange_LD), srcStr = 'from candidates'; else, srcStr = 'explicit ampRange_LD'; end
fprintf('Peak amplitude range: %.3f – %.3f spike-height (%s)\n', ampBand_LD(1), ampBand_LD(2), srcStr);
rowsLD = find(EventTab(:,4) >= ampBand_LD(1) & EventTab(:,4) <= ampBand_LD(2));
if numel(rowsLD) < 20                       % relax if the band is too tight
    rowsLD = (1:size(EventTab,1))';
    fprintf('Amplitude band too tight (<20 events) -> using all %d events above threshold\n', numel(rowsLD));
end

% --- candidate trigger ROI ---------------------------------------------------
% exR0 (from the L50 section) is trigger_roi(): the ROI with the largest window-mean
% SUBTHRESHOLD voltage at the candidate frame — an estimate from the traces, which
% can land on the wrong ROI when the active dendrite is poorly covered. Since the
% polyline was drawn on the movie through the active dendrite, the movie itself
% defines the location: the peak of the polyline profile is mapped to the ROI whose
% footprint covers that pixel (nearest footprint centre if it covers none).
%   'movie'  : take it from the polyline peak (default; falls back to the trace)
%   'trace'  : keep the trigger_roi() estimate
%   'manual' : use candROIman, one ROI index per candidate
% candROIman_numbering says how candROIman is indexed: 'sorted' = this section's
% indexing (traces/kymographs/LDTab are all in PCresult{f}.sortDist order, which is
% what the viewer and figure 146 titles print), 'original' = the ROI id before
% sorting. Both indices are printed so the choice can be checked at a glance.
candROIsrc = 'manual';
candROIman = [17; 17; 17; 17; 12];   % trigger ROI per candidate (NaN = fall back)
candROIman_numbering = 'sorted';

L50movFile_LD = fullfile(save_figto, 'ExampleL50_movie.mat');
movRec = cell(nEx,1);
if isfile(L50movFile_LD)
    Sv_LD = load(L50movFile_LD);
    for i = 1:min(nEx, numel(Sv_LD.MovieL50))
        if ~isempty(Sv_LD.MovieL50(i).rec), movRec{i} = Sv_LD.MovieL50(i).rec; end
    end
else
    warning('%s not found — candidate length and trigger ROI fall back to the traces', L50movFile_LD);
end

exR0_use = exR0;  r0src = repmat({'trace'}, nEx, 1);
for i = 1:nEx
    fe = exampleList(i,1);
    ftp = PCresult{fe}.Ftprnts(:,:,PCresult{fe}.sortDist);
    switch candROIsrc
        case 'manual'
            if ~isnan(candROIman(i))
                if strcmpi(candROIman_numbering, 'original')
                    j = find(PCresult{fe}.sortDist == candROIman(i), 1);
                else
                    j = candROIman(i);
                end
                if isempty(j) || j < 1 || j > size(ftp,3)
                    warning('Candidate %d: ROI %d (%s numbering) not found in f%d — keeping the trace estimate', ...
                        i, candROIman(i), candROIman_numbering, fe);
                else
                    exR0_use(i) = j;  r0src{i} = 'manual';
                end
            end
        case 'movie'
            rm = movie_trigger_roi(movRec{i}, ftp);
            if ~isnan(rm), exR0_use(i) = rm; r0src{i} = 'movie'; end
    end
    origID = PCresult{fe}.sortDist(exR0_use(i));      % ROI id before sorting
    if exR0_use(i) ~= exR0(i)
        Dm = PCresult{fe}.interDendDist(PCresult{fe}.sortDist, PCresult{fe}.sortDist);
        Dm = max(cat(3,Dm,Dm'),[],3);
        fprintf(['Candidate %d: trigger ROI %d (%s; original id %d) instead of %d (trace) ' ...
            '— %.0f µm apart\n'], i, exR0_use(i), r0src{i}, origID, exR0(i), Dm(exR0(i), exR0_use(i)));
    else
        fprintf('Candidate %d: trigger ROI %d (%s; original id %d)\n', ...
            i, exR0_use(i), r0src{i}, origID);
    end
end

% --- candidate measurements (at exR0_use) ---
exLEN=nan(nEx,1); exDUR=nan(nEx,1); exAmpSH_LD=nan(nEx,1);
exTPRE=nan(nEx,1); exTPOST=nan(nEx,1);
for i = 1:nEx
    fe = exampleList(i,1);  te = exampleList(i,2);
    SubSH = interpolateNaN(PCresult{fe}.Subthreshold / PCresult{fe}.SpikeHeight);  % peak location
    RawSH = interpolateNaN(PCresult{fe}.NormalizedTrace_dirt / PCresult{fe}.SpikeHeight);  % measurement
    Dmat = PCresult{fe}.interDendDist(PCresult{fe}.sortDist, PCresult{fe}.sortDist);
    Dmat = max(cat(3,Dmat,Dmat'),[],3);
    [exLEN(i), exDUR(i), exAmpSH_LD(i), exTPRE(i), exTPOST(i)] = event_len_dur(SubSH, RawSH, Dmat, ...
        exR0_use(i), te, avgwin_ext, Dmax_ext, Vthr_sh, baseWin_LD, durMax_LD, smoothDur_LD);
    fprintf(['Candidate %d (ROI-based): amp=%.3f spike-height (%.4f ∆F/F)  length=%.0f µm  ' ...
        'duration=%.0f ms (%.0f -> %+.0f ms)  [thr %.3f sh = %.4f ∆F/F]\n'], ...
        i, exAmpSH_LD(i), exAmpSH_LD(i)*PCresult{fe}.SpikeHeight, exLEN(i), exDUR(i), ...
        exTPRE(i), exTPOST(i), Vthr_sh, Vthr_sh*PCresult{fe}.SpikeHeight);
    if isnan(exLEN(i)) || isnan(exDUR(i))       % explain the NaN instead of hiding it
        if isnan(exAmpSH_LD(i))
            fprintf('   -> NaN: ROI %d of f%d has no usable raw trace at frame %d (all NaN)\n', ...
                exR0_use(i), fe, te);
        else
            fprintf(['   -> NaN: the raw +/-%d ms amplitude at the trigger (%.3f spike-height = ' ...
                '%.4f ∆F/F) is below Vthr_sh = %.3f. Lower Vthr_sh to measure this event.\n'], ...
                avgwin_ext, exAmpSH_LD(i), exAmpSH_LD(i)*PCresult{fe}.SpikeHeight, Vthr_sh);
        end
    end
end

% --- candidate LENGTH: measured exactly like every other event -----------------
% The candidates are treated as ordinary events: length comes from the ROI profile
% (event_len_dur above), the same way the control events are measured, so the two
% groups are compared on identical footing. Setting useMovieLen_LD = true instead
% takes the candidate length from the polyline movie profile in ExampleL50_movie.mat
% (scaled by the candidate's amplitude, FW/2 like L50_half) — finer than the ROI
% grid, but then the candidates and the controls are no longer measured alike.
useMovieLen_LD = true;
exLEN_mov = nan(nEx,1);  exFW_mov = nan(nEx,1);
if useMovieLen_LD
    for i = 1:nEx
        if isempty(movRec{i}), continue; end
        [exLEN_mov(i), exFW_mov(i)] = movie_len(movRec{i}, exAmpSH_LD(i), Vthr_sh);
        fprintf('Candidate %d (movie)    : length=%.0f µm (full width %.0f µm)\n', ...
            i, exLEN_mov(i), exFW_mov(i));
        if isnan(exLEN_mov(i))
            fprintf(['   -> movie NaN: the polyline profile scaled by amp=%.3f never reaches ' ...
                'Vthr_sh=%.3f (peak of the scaled profile = %.3f)\n'], ...
                exAmpSH_LD(i), Vthr_sh, exAmpSH_LD(i));
        end
    end
end
useMov  = ~isnan(exLEN_mov);
exLEN_use = exLEN;  exLEN_use(useMov) = exLEN_mov(useMov);
fprintf('Candidate length source: ROI profile for %d/%d, movie for %d\n', ...
    sum(~useMov), nEx, sum(useMov));

% --- all amplitude-matched events ---
% LDTab columns: [f r0 t0 amp_sh length_um duration_ms tPre tPost isCand censSpace censTime]
LDTab = nan(numel(rowsLD), 11);
for f = unique(EventTab(rowsLD,1))'
    SubSH = interpolateNaN(PCresult{f}.Subthreshold / PCresult{f}.SpikeHeight);          % peak location
    RawSH = interpolateNaN(PCresult{f}.NormalizedTrace_dirt / PCresult{f}.SpikeHeight);  % measurement
    Dmat = PCresult{f}.interDendDist(PCresult{f}.sortDist, PCresult{f}.sortDist);
    Dmat = max(cat(3,Dmat,Dmat'),[],3);
    candFr = exampleList(exampleList(:,1)==f, 2);
    for k = find(EventTab(rowsLD,1)==f)'
        rr = rowsLD(k);  r0 = EventTab(rr,2);  t0 = EventTab(rr,3);
        [len, dur, aSH, tPre, tPost, cens] = event_len_dur(SubSH, RawSH, Dmat, r0, t0, ...
            avgwin_ext, Dmax_ext, Vthr_sh, baseWin_LD, durMax_LD, smoothDur_LD);
        isC = any(abs(candFr - t0) <= 200);
        LDTab(k,:) = [f r0 t0 aSH len dur tPre tPost isC cens(1) cens(2)];
    end
    fprintf('  file %d done (%d/%d events)\n', f, sum(~isnan(LDTab(:,5))), size(LDTab,1));
end
okLD = ~isnan(LDTab(:,5)) & ~isnan(LDTab(:,6));       % trigger above threshold

% --- reject peaks near spikes / CS bursts / blue stimulation ---------------
% peakvec already dropped peaks at -10..+40 ms around somatic spikes and on the
% CS-trace frames themselves, but NOT a margin around the CS bursts and NOT the
% blue-stimulation periods. Here every event whose peak frame falls within
% spkExcl_LD ms of a somatic spike or of a CS-trace frame is removed, and (with
% blueOff_LD) so is everything within blueDil_LD ms of blue light — the same
% dilation the rest of this script uses (imdilate(BlueStim>0, ones(1,201))).
% Candidates are not removed; a warning is printed if one falls inside a window.
spkExcl_LD = 30;       % +/- ms around somatic spikes and CS-trace frames
blueOff_LD = true;     % keep blue-OFF periods only
blueDil_LD = 100;      % +/- ms dilation around blue-on frames
nearSpk = false(size(okLD));  nearBlue = false(size(okLD));
for f = unique(LDTab(~isnan(LDTab(:,1)),1))'
    nT  = size(PCresult{f}.Subthreshold,2);
    bad = false(1,nT);
    bad(PCresult{f}.SpikeMat(1,:) > 0) = true;         % somatic spikes
    bad(find(PCresult{f}.CStrace > 0)) = true;         %#ok<FNDSB> % CS bursts
    bad = imdilate(bad, ones(1, 2*spkExcl_LD+1));      % +/- spkExcl_LD ms around both
    if blueOff_LD
        blu = imdilate(PCresult{f}.BlueStim > 0, ones(1, 2*blueDil_LD+1));
    else
        blu = false(1,nT);
    end
    idx = find(LDTab(:,1)==f & ~isnan(LDTab(:,3)));
    nearSpk(idx)  = bad(LDTab(idx,3));
    nearBlue(idx) = blu(LDTab(idx,3));
end
isCandRow = LDTab(:,9)==1;
if any((nearSpk | nearBlue) & isCandRow)
    warning(['%d candidate event(s) fall within the spike/CS/blue exclusion window ' ...
        '— they are KEPT, but check them.'], sum((nearSpk|nearBlue) & isCandRow));
end
dropEvt = (nearSpk | nearBlue) & ~isCandRow;           % controls only
fprintf('Excluded near spikes/CS (+/-%d ms): %d;  blue-on (+/-%d ms): %d;  total removed %d/%d\n', ...
    spkExcl_LD, sum(nearSpk & okLD & ~isCandRow), blueDil_LD, ...
    sum(nearBlue & okLD & ~isCandRow), sum(dropEvt & okLD), sum(okLD));
okLD = okLD & ~dropEvt;

% --- MANUAL artifact rejection: hand-drawn regions on the length-duration plane ---
% Floating bubbles and other artifacts that survive the automatic cutoffs cluster in
% the length-duration plane, so they can be lassoed by hand. Set drawReject_LD =
% true once, draw one polygon per artifact cluster (double-click to close each one,
% then answer the dialog), and the regions are saved to ManualReject_LD.mat and
% re-applied automatically on later runs — set the flag back to false afterwards.
% Delete that file (or set rejPoly = {}) to start over. Like the noise cutoffs, the
% manual regions reject CONTROL events only; the hand-picked candidates are kept.
manualRejFile  = fullfile(save_figto, 'ManualReject_LD.mat');
rejPoly = {};
if isfile(manualRejFile), Sr = load(manualRejFile); rejPoly = Sr.rejPoly; end
if drawReject_LD
    figure(144); clf; hold all;
    scatter(LDTab(okLD,5), LDTab(okLD,6), 16, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha',0.4);
    plot(exLEN_use, exDUR, 'ro', 'markersize',9, 'linewidth',1.5, 'MarkerFaceColor','r');
    xlabel('Length: crossing distance (µm)'); ylabel('Duration: t_{post} - t_{pre} (ms)');
    title('Enclose ARTIFACT clusters — double-click to close each polygon');
    for p = 1:numel(rejPoly)                        % show what is already saved
        plot([rejPoly{p}(:,1); rejPoly{p}(1,1)], [rejPoly{p}(:,2); rejPoly{p}(1,2)], 'r-');
    end
    more2draw = true;
    while more2draw
        if exist('drawpolygon','file')              % Image Processing Toolbox
            h = drawpolygon('Color','r');
            pos = [];  if isvalid(h), pos = h.Position; end
        else                                        % fallback: click, Enter to finish
            title('Click the vertices, press Enter to close this polygon');
            [px, py] = ginput;  pos = [px py];
        end
        if size(pos,1) >= 3
            rejPoly{end+1} = pos;
            plot([pos(:,1); pos(1,1)], [pos(:,2); pos(1,2)], 'r-', 'linewidth',1.5);
            fprintf('Region %d: %d events inside\n', numel(rejPoly), ...
                sum(inpolygon(LDTab(:,5), LDTab(:,6), pos(:,1), pos(:,2)) & okLD));
        end
        more2draw = strcmp(questdlg('Draw another artifact region?','Manual rejection', ...
            'Yes','No','No'), 'Yes');
    end
    save(manualRejFile, 'rejPoly');
    fprintf('Saved %d rejection region(s) -> %s  (set drawReject_LD = false now)\n', ...
        numel(rejPoly), manualRejFile);
end
manRej = false(size(okLD));
for p = 1:numel(rejPoly)
    manRej = manRej | inpolygon(LDTab(:,5), LDTab(:,6), rejPoly{p}(:,1), rejPoly{p}(:,2));
end
manRej = manRej & okLD & LDTab(:,9)==0;               % controls only, candidates untouched

% blood-vessel / motion noise rejection
tooShortL = okLD & LDTab(:,5) <  lenMin_LD;
tooShortD = okLD & LDTab(:,6) <  durMin_LD;
keepLD    = okLD & LDTab(:,5) >= lenMin_LD & LDTab(:,6) >= durMin_LD & ~manRej;
rawCtrl   = keepLD & LDTab(:,9)==0;                   % matched controls (candidates excluded)
fprintf('Manual rejection: %d region(s), %d events removed\n', numel(rejPoly), sum(manRej));
fprintf('Amplitude-matched events: %d (band %.3f–%.3f spike-height), %d above %.0f%% spike height\n', ...
    numel(rowsLD), ampBand_LD(1), ampBand_LD(2), sum(okLD), 100*Vthr_sh);
fprintf('Noise cutoff (length<%g µm or duration<%g ms): rejected %d (%d short-length, %d short-duration) -> %d kept\n', ...
    lenMin_LD, durMin_LD, sum(okLD)-sum(keepLD), sum(tooShortL), sum(tooShortD), sum(keepLD));
fprintf('Censored (no crossing found): space %d/%d, time %d/%d\n', ...
    sum(LDTab(keepLD,10)==1), sum(keepLD), sum(LDTab(keepLD,11)==1), sum(keepLD));
% events confined to the trigger ROI: extent below the ROI sampling resolution
fprintf('Unresolved length (only the trigger ROI above threshold, len=0): %d/%d events\n', ...
    sum(LDTab(keepLD,10)==2), sum(keepLD));
for f = unique(LDTab(okLD,1))'
    Dm = PCresult{f}.interDendDist(PCresult{f}.sortDist, PCresult{f}.sortDist);
    Dm = max(cat(3,Dm,Dm'),[],3);  Dm(Dm<=0) = NaN;
    fprintf('   f%-3d median nearest-ROI spacing %.0f µm -> lengths below that are unresolved\n', ...
        f, median(min(Dm,[],2,'omitnan'),'omitnan'));
end

% --- balance the control pool across neurons -------------------------------
% Event counts are very uneven (a single neuron can supply half of all matched
% events), so a pooled event average mostly reports that one neuron. Each neuron
% is therefore capped at nMaxPerN_LD events - the ones closest in peak amplitude
% to the candidates' mean, i.e. the best-matched ones - and every summary below is
% computed NEURON-LEVEL: mean per neuron first, then mean +/- sem across neurons
% (neurons with < minEvtPerN_LD kept events are dropped).
fneu  = LDTab(:,1);
ampSH = EventTab(rowsLD,4);                            % peak amplitude, spike-height units
aTarget = mean(exAMP);
ufRaw = unique(fneu(rawCtrl))';
fprintf('Per-neuron matched events (before balancing, %d total):\n', sum(rawCtrl));
for f = ufRaw
    n = sum(rawCtrl & fneu==f);
    fprintf('   f%-3d %5d  (%4.1f%%)\n', f, n, 100*n/sum(rawCtrl));
end
ctrlLD = false(size(rawCtrl));
for f = ufRaw
    idx = find(rawCtrl & fneu==f);
    [~, ord] = sort(abs(ampSH(idx) - aTarget));        % best amplitude match first
    ctrlLD(idx(ord(1:min(nMaxPerN_LD,numel(idx))))) = true;
end
fprintf('After capping at %d events/neuron: %d events from %d neurons (max share %.1f%%)\n', ...
    nMaxPerN_LD, sum(ctrlLD), numel(unique(fneu(ctrlLD))), ...
    100*max(arrayfun(@(f) sum(ctrlLD & fneu==f), unique(fneu(ctrlLD))))/sum(ctrlLD));

% neuron-level summaries (the unit of replication is the neuron, not the event)
ufK  = unique(fneu(ctrlLD))';
nCnt = arrayfun(@(f) sum(ctrlLD & fneu==f), ufK);
neuL = arrayfun(@(f) mean(LDTab(ctrlLD & fneu==f,5),'omitnan'), ufK);
neuD = arrayfun(@(f) mean(LDTab(ctrlLD & fneu==f,6),'omitnan'), ufK);
selN = nCnt >= minEvtPerN_LD;
neuL = neuL(selN);  neuD = neuD(selN);  ufK_sel = ufK(selN);
fprintf('Neuron-level pool: %d/%d neurons with >=%d events\n', ...
    numel(ufK_sel), numel(ufK), minEvtPerN_LD);

% --- group means (individual candidates are not plotted separately) ---
% The noise cutoff is NOT applied to the candidates: they are hand-picked local
% events, so only candidates that could not be measured at all (NaN) are dropped.
exKeep = ~isnan(exLEN_use) & ~isnan(exDUR);
if any(~exKeep)
    warning('%d/%d candidate(s) could not be measured (NaN), excluded: %s', ...
        sum(~exKeep), nEx, num2str(find(~exKeep)'));
end
belowCut = exKeep & (exLEN_use < lenMin_LD | exDUR < durMin_LD);
if any(belowCut)
    fprintf('Note: candidate(s) %s fall below the control noise cutoff but are KEPT\n', ...
        num2str(find(belowCut)'));
end
exLENk = exLEN_use(exKeep);  exDURk = exDUR(exKeep);
ctrlL = LDTab(ctrlLD,5);  ctrlD = LDTab(ctrlLD,6);     % pooled (balanced) events
mL_c = mean(exLENk,'omitnan');  eL_c = std(exLENk,'omitnan')/sqrt(numel(exLENk));
mD_c = mean(exDURk,'omitnan');  eD_c = std(exDURk,'omitnan')/sqrt(numel(exDURk));
mL_k = mean(neuL,'omitnan');    eL_k = std(neuL,'omitnan')/sqrt(numel(neuL));   % across neurons
mD_k = mean(neuD,'omitnan');    eD_k = std(neuD,'omitnan')/sqrt(numel(neuD));

% --- figure ---
figure(143); clf; tiledlayout(1,3,'TileSpacing','compact');

% A: length vs duration — control events (gray) with the two group means +/- sem
nexttile([1 1]); hold all;
% rejected blood-vessel-like events, shown faintly behind the cutoff lines
scatter(LDTab(okLD & ~keepLD,5), LDTab(okLD & ~keepLD,6), 10, [0.85 0.85 0.85], 'filled', ...
    'MarkerFaceAlpha',0.5, 'MarkerEdgeColor','none');
scatter(ctrlL, ctrlD, 14, [0.6 0.6 0.6], 'filled', ...
    'MarkerFaceAlpha',0.3, 'MarkerEdgeColor','none');
xline(lenMin_LD, ':', 'color',[0.3 0.3 0.3], 'Label','cutoff');
yline(durMin_LD, ':', 'color',[0.3 0.3 0.3]);
for p = 1:numel(rejPoly)                              % hand-drawn artifact regions
    plot([rejPoly{p}(:,1); rejPoly{p}(1,1)], [rejPoly{p}(:,2); rejPoly{p}(1,2)], ...
        '-', 'color',[0.9 0.5 0.5], 'linewidth',1);
end
hK = errorbar(mL_k, mD_k, eD_k, eD_k, eL_k, eL_k, 'ko', ...
    'markersize',9, 'linewidth',1.5, 'MarkerFaceColor','k', 'CapSize',0);
hC = errorbar(mL_c, mD_c, eD_c, eD_c, eL_c, eL_c, 'ro', ...
    'markersize',9, 'linewidth',1.5, 'MarkerFaceColor','r', 'CapSize',0);
xlabel('Length: crossing distance (µm)'); ylabel('Duration: t_{post} - t_{pre} (ms)');
legend([hC hK], {sprintf('Local (n=%d)',numel(exLENk)), ...
    sprintf('Matched mean \pm sem (%d neurons)',numel(neuL))}, ...
    'box','off','location','northeast');
title(sprintf('A  raw-trace crossings at %.0f%% spike height', 100*Vthr_sh)); box off;

% B: length distribution + candidate mean +/- sem
nexttile([1 1]); hold all;
histogram(ctrlL, 0:20:Dmax_ext, 'FaceColor',[0.4 0.4 0.4], 'EdgeColor','none');
yl = ylim;
xline(mL_c, 'r-', 'linewidth',2, 'Label',sprintf('local %.0f µm', mL_c));
xline(mL_k, 'k--', 'linewidth',1.5);
ylim(yl);
xlabel('Length: crossing distance (µm)'); ylabel('Event count'); box off;
title(sprintf('B  length: local %.0f±%.0f vs matched %.0f±%.0f µm (per neuron)', mL_c, eL_c, mL_k, eL_k));
fprintf('%2.0f %% of events exceed 50 µm\n',sum(ctrlL>50)./length(ctrlL)*100);

% C: duration distribution + candidate mean +/- sem
nexttile([1 1]); hold all;
histogram(ctrlD, 0:5:2*durMax_LD, 'FaceColor',[0.4 0.4 0.4], 'EdgeColor','none');
yl = ylim;
xline(mD_c, 'r-', 'linewidth',2, 'Label',sprintf('local %.0f ms', mD_c));
xline(mD_k, 'k--', 'linewidth',1.5);
ylim(yl);
xlabel('Duration: t_{post} - t_{pre} (ms)'); ylabel('Event count'); box off;
title(sprintf('C  duration: local %.0f±%.0f vs matched %.0f±%.0f ms (per neuron)', mD_c, eD_c, mD_k, eD_k));

set_font('Arial'); set_fontsize(12);
set_figsize(420,140);

% --- Figure 147: valid events only (clean version of the same three panels) ---
figure(147); clf; tiledlayout(1,3,'TileSpacing','compact');

nexttile([1 1]); hold all;
scatter(ctrlL, ctrlD, 14, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha',0.35, 'MarkerEdgeColor','none');
hK7 = errorbar(mL_k, mD_k, eD_k, eD_k, eL_k, eL_k, 'ko', ...
    'markersize',9, 'linewidth',1.5, 'MarkerFaceColor','k', 'CapSize',0);
hC7 = errorbar(mL_c, mD_c, eD_c, eD_c, eL_c, eL_c, 'ro', ...
    'markersize',9, 'linewidth',1.5, 'MarkerFaceColor','r', 'CapSize',0);
xlabel('Length: crossing distance (µm)'); ylabel('Duration: t_{post} - t_{pre} (ms)');
legend([hC7 hK7], {sprintf('Local (n=%d)',numel(exLENk)), ...
    sprintf('Amplitude-matched (n=%d)',numel(ctrlL))}, 'box','off','location','northeast');
title(sprintf('Crossings at %.0f%% spike height', 100*Vthr_sh)); box off;

nexttile([1 1]); hold all;
histogram(ctrlL, 0:10:Dmax_ext, 'FaceColor',[0.4 0.4 0.4], 'EdgeColor','none');
yl = ylim;
xline(mL_c, 'r-', 'linewidth',2, 'Label',sprintf('local %.0f µm', mL_c));
xline(mL_k, 'k--', 'linewidth',1.5);
ylim(yl);
xlabel('Length: crossing distance (µm)'); ylabel('Event count'); box off;
title(sprintf('Length (median %.0f µm)', median(ctrlL,'omitnan')));

nexttile([1 1]); hold all;
histogram(ctrlD, 0:5:2*durMax_LD, 'FaceColor',[0.4 0.4 0.4], 'EdgeColor','none');
yl = ylim;
xline(mD_c, 'r-', 'linewidth',2, 'Label',sprintf('local %.0f ms', mD_c));
xline(mD_k, 'k--', 'linewidth',1.5);
ylim(yl);
xlabel('Duration: t_{post} - t_{pre} (ms)'); ylabel('Event count'); box off;
title(sprintf('Duration (median %.0f ms)', median(ctrlD,'omitnan')));

set_font('Arial'); set_fontsize(12);
set_figsize(420,140);

% --- how many events extend beyond lenPass_LD ---
lenPass_LD = 50;       % µm
nPass = sum(ctrlL > lenPass_LD);
fprintf('Events longer than %g µm: %d / %d valid events (%.1f%%)\n', ...
    lenPass_LD, nPass, numel(ctrlL), 100*nPass/numel(ctrlL));
for f = unique(LDTab(ctrlLD,1))'
    cl = LDTab(ctrlLD & LDTab(:,1)==f, 5);
    fprintf('   f%-3d %d / %d (%.1f%%)\n', f, sum(cl > lenPass_LD), numel(cl), ...
        100*sum(cl > lenPass_LD)/numel(cl));
end
fprintf('Local candidates longer than %g µm: %d / %d\n', ...
    lenPass_LD, sum(exLENk > lenPass_LD), numel(exLENk));

% --- group comparison (means, not individual candidates) ---
fprintf('Local candidates    : length %.0f ± %.0f µm (%s), duration %.0f ± %.0f ms (n=%d events)\n', ...
    mL_c, eL_c, repmat_str(any(useMov),'movie','ROI profile'), mD_c, eD_c, numel(exLENk));
fprintf('Amp-matched (neuron): length %.0f ± %.0f µm, duration %.0f ± %.0f ms (n=%d neurons, %d events)\n', ...
    mL_k, eL_k, mD_k, eD_k, numel(neuL), sum(ctrlLD));
% same-neuron control: the cleanest comparison, since all candidates come from one neuron
for fe = unique(exampleList(:,1))'
    s = ctrlLD & fneu==fe;
    if any(s)
        fprintf('Same neuron f%-3d    : length %.0f ± %.0f µm, duration %.0f ± %.0f ms (n=%d events)\n', ...
            fe, mean(LDTab(s,5),'omitnan'), std(LDTab(s,5),'omitnan')/sqrt(sum(s)), ...
            mean(LDTab(s,6),'omitnan'), std(LDTab(s,6),'omitnan')/sqrt(sum(s)), sum(s));
    end
end
fprintf('Candidate mean vs neuron means: %.0fth pct of length, %.0fth pct of duration\n', ...
    100*mean(neuL <= mL_c), 100*mean(neuD <= mD_c));
if exist('ranksum','file')
    fprintf('Wilcoxon rank sum (candidate events vs neuron means): length p=%.3g, duration p=%.3g\n', ...
        ranksum(exLENk, neuL), ranksum(exDURk, neuD));
end

%% ===== Interactive: click a scatter point -> that event's kymograph =====
% Same events as the panel-A scatter. interactive_scatter_matrix_viewer(Mat1,A,B)
% takes A/B as the scatter axes (here length and duration) and Mat1{i} as the
% matrix shown when point i is clicked — here the event's peri-event kymograph:
% RAW trace (spike-height units, local baseline subtracted), ROIs sorted by contour
% distance from the trigger ROI (row 1 = trigger), over +/- kymoHalf_LD ms.
% ViewTab holds the event identity of every point: [f r0 t0 amp len dur ... isCand].
% Slices are passed as an N x T x nEvent array (all slices must share one size).
kymoHalf_LD = 1200;     % +/- ms shown per event
memMaxGB_LD = 2;       % memory ceiling for the kymograph stack; only then is it thinned

% every point drawn in panel A (all kept controls) + the candidates
viewRows = find(ctrlLD);
nRmax = max(arrayfun(@(f) size(PCresult{f}.NormalizedTrace_dirt,1), unique(LDTab(ctrlLD,1))));
bytesPer = nRmax * (2*kymoHalf_LD+1) * 4;                    % single precision
nFit = floor(memMaxGB_LD*2^30 / bytesPer);
if numel(viewRows) > nFit                                    % thin evenly, keep the spread
    warning(['Viewer would need %.1f GB for all %d events — showing %d of them ' ...
        '(evenly spaced). Raise memMaxGB_LD or lower kymoHalf_LD to show all.'], ...
        numel(viewRows)*bytesPer/2^30, numel(viewRows), nFit);
    viewRows = viewRows(round(linspace(1, numel(viewRows), nFit)));
end
% candidates are appended from the candidate arrays (their trigger ROI is exR0_use,
% which is not necessarily the ROI that peakvec fired on for the same event)
% ALL candidates are kept here, including those whose length/duration came out NaN
% — their kymograph is exactly what is needed to see why the measurement failed.
% (A NaN candidate has no scatter position, so it appears in figure 146 but cannot
% be clicked in the viewer.)
candTab = [exampleList(:,1) exR0_use exampleList(:,2) exAmpSH_LD exLEN_use exDUR ...
    exTPRE exTPOST ones(nEx,1) zeros(nEx,2)];
ViewTab = [LDTab(viewRows,:); candTab];
viewRows = [viewRows; nan(size(candTab,1),1)];
% One fixed-size slice per event (the viewer needs every slice the same size):
% rows = max ROI count over the neurons involved, columns = the full +/-kymoHalf_LD
% window with the peak always in the centre column. Missing rows/frames stay NaN.
nRmax = max(arrayfun(@(f) size(PCresult{f}.NormalizedTrace_dirt,1), unique(ViewTab(:,1))));
tAxis = -kymoHalf_LD:kymoHalf_LD;
KymoArr = nan(nRmax, numel(tAxis), numel(viewRows), 'single');
viewRowLen = nan(numel(viewRows),1);
for f = unique(ViewTab(:,1))'
    RawSH = interpolateNaN(PCresult{f}.NormalizedTrace_dirt / PCresult{f}.SpikeHeight);
    Dmat = PCresult{f}.interDendDist(PCresult{f}.sortDist, PCresult{f}.sortDist);
    Dmat = max(cat(3,Dmat,Dmat'),[],3);
    [nR, nT] = size(RawSH);
    for k = find(ViewTab(:,1)==f)'
        r0 = ViewTab(k,2);  t0 = ViewTab(k,3);
        W = t0 + tAxis;  inb = W >= 1 & W <= nT;    % keep the peak centred
        B = t0 + (-baseWin_LD:baseWin_LD);    B(B<1 | B>nT) = [];
        B = setdiff(B, t0 + (-3*avgwin_ext:3*avgwin_ext));
        if isempty(B), B = 1:nT; end
        base = median(RawSH, 2, 'omitnan');
        [dsrt, dord] = sort(Dmat(r0,:));            % ROIs by distance from the trigger
        KymoArr(1:nR, inb, k) = RawSH(dord,W(inb)) - base(dord);
        viewRowLen(k) = sum(dsrt <= ViewTab(k,5));  % last ROI row inside the measured length
    end
end
% --- candidates in the same kymograph format, as a static figure ---------------
% Same slices the viewer shows (row 1 = trigger ROI, rows ordered by contour
% distance from it; time centred on the peak). Dashed white lines mark the measured
% duration (tPre..tPost) and the last ROI row inside the measured length.
candIdx = find(ViewTab(:,9)==1);
if ~isempty(candIdx)
    sub = KymoArr(1:37:end);                         % subsample for a fast color axis
    caxK = [prctile(sub,5) prctile(sub,99)];
    figure(146); clf; tiledlayout(1, numel(candIdx), 'TileSpacing','compact');
    for k = 1:numel(candIdx)
        kk = candIdx(k);  K = KymoArr(:,:,kk);
        nRk = find(any(~isnan(K),2), 1, 'last');
        nexttile; hold all;
        imagesc(tAxis, 1:nRk, K(1:nRk,:), caxK);
        plot([ViewTab(kk,7) ViewTab(kk,7)], [0.5 nRk+0.5], 'w--', 'linewidth',1);
        plot([ViewTab(kk,8) ViewTab(kk,8)], [0.5 nRk+0.5], 'w--', 'linewidth',1);
        if ~isnan(viewRowLen(kk)) && viewRowLen(kk) >= 1
            plot([tAxis(1) tAxis(end)], (viewRowLen(kk)+0.5)*[1 1], 'w:', 'linewidth',1);
        end
        axis xy tight;
        xlabel('Time from peak (ms)');
        if k == 1, ylabel('ROI (sorted by distance from trigger)'); end
        title(viewLbl_one(ViewTab, kk), 'Interpreter','none', 'FontSize',9);
    end
    colormap(turbo);
    cb = colorbar;  cb.Label.String = 'Raw (spike height)';
    set_font('Arial'); set_fontsize(11); set_figsize(190*numel(candIdx), 210);
end

% one identifying label per point: neuron, trigger ROI, peak frame, measurements
viewLbl = arrayfun(@(k) viewLbl_one(ViewTab,k), 1:size(ViewTab,1), 'UniformOutput', false);
interactive_scatter_matrix_viewer(KymoArr, ViewTab(:,5), ViewTab(:,6), viewLbl);
% label the scatter axes of the viewer (the function labels them 'A'/'B')
for a = findall(gcf,'Type','axes')'
    if ~isempty(findobj(a,'Type','scatter'))
        xlabel(a,'Length (µm)'); ylabel(a,'Duration (ms)');
        title(a,'Click a point — kymograph on the right');
    end
end
fprintf('Viewer: %d events (%d controls + %d candidates), kymo +/-%d ms\n', ...
    size(ViewTab,1), sum(ViewTab(:,9)==0), sum(ViewTab(:,9)==1), kymoHalf_LD);
if any(isnan(ViewTab(ViewTab(:,9)==1,5)) | isnan(ViewTab(ViewTab(:,9)==1,6)))
    fprintf('  %d candidate(s) have NaN length/duration: drawn in figure 146, not clickable in the viewer\n', ...
        sum(isnan(ViewTab(:,5)) | isnan(ViewTab(:,6))));
end
fprintf('  clicking a point prints f, ROI and frame (also shown in the panel title)\n');
fprintf('  ViewTab(i,:) = [f r0 t0 amp len dur tPre tPost isCand censSpace censTime]\n');
fprintf('  panel-A control events: %d, in viewer: %d%s\n', sum(ctrlLD), sum(ViewTab(:,9)==0), ...
    repmat(' (thinned — see warning)', 1, sum(ViewTab(:,9)==0) < sum(ctrlLD)));
if any(ViewTab(:,9)==1)
    fprintf('  candidate points: %s\n', num2str(find(ViewTab(:,9)==1)'));
end

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

function [len, dur, ampTrig, tPre, tPost, cens, prof] = event_len_dur(Sub, Raw, Dmat, r0, t0, avgwin, Dmax, Vthr, baseWin, durMax, smoothDur)
% Threshold-crossing length and duration of one event, measured outward from the
% trigger point (peak ROI r0, peak frame t0).
% The peak was located on the SUBTHRESHOLD trace (Sub, used here only to seed the
% time search); LENGTH and DURATION are both measured on the RAW trace (Raw), which
% is not broadened by the subthreshold filter and therefore gives shorter durations.
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
base = median(Raw(:,B), 2, 'omitnan');         % baseline from the RAW trace
V = mean(Raw(:,W), 2, 'omitnan') - base;       % +/-avgwin mean of Raw at the peak
ampTrig = V(r0);
if ~(ampTrig >= Vthr), return; end              % trigger point itself below threshold

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
v  = movmean(Raw(r0,:) - base(r0), smoothDur, 'omitnan');
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
