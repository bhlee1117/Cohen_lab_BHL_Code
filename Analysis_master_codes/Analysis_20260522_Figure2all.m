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
    end

    PCresult{f}.BlueStim = Result.Blue;
    PCresult{f}.CStrace = CS_trace;
    PCresult{f}.BStrace = BS_trace;
    PCresult{f}.avgImg = Result.ref_im;
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

NormalizedTrace_noNaN=interpolateNaN(PCresult{f}.NormalizedTrace_dirt);
filteredNormTr2 = pcafilterTrace(NormalizedTrace_noNaN, setdiff([1:5],[4 5 6 7 10]));
ftprnt=PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist);
[nROI nTime]=size(filteredNormTr2);

Show_tr=[mean(filteredNormTr2(basalind,:),1,'omitnan'); mean(filteredNormTr2(somaind,:),1,'omitnan'); mean(filteredNormTr2(apicalind,:),1,'omitnan')];
Show_tr=interpolateNaN(Show_tr);
Show_subtr=get_subthreshold(Show_tr,PCresult{f}.SpikeMat(1,:)>0,7,20);
Show_ftprnt=cat(3,max(ftprnt(:,:,basalind)>0,[],3),max(ftprnt(:,:,somaind)>0,[],3),max(ftprnt(:,:,apicalind)>0,[],3));
roi_show=setdiff([1:nROI],[10 12 13 14 15 21]); %12 13 21 12 21

figure(5); clf; tiledlayout(4,1);
tax=[1:nTime]/1000; lscale=[0.2];

ax3=nexttile([3 1]);
l=plot(tax(45000:160000),Show_tr(:,45000:160000)'-[1:3]*lscale);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap,2));
hold all;
l2=plot(tax(45000:160000),Show_subtr(:,45000:160000)'-[1:3]*lscale);
arrayfun(@(l,c) set(l,'Color',c{:}),l2,num2cell(cmap/2,2));
%legend(l,{'Basal','Soma','Distal'})
axis off
drawScaleBar(lscale/2,'vertical','position',[162 -0.25],'color','k');
ylim([-0.67 0]);

ax2=nexttile([1 1]);
plot(tax(45000:160000),PCresult{f}.VR(5,45000:160000),'color',[0 0 0],'linewidth',2);
axis off
linkaxes([ax2 ax3],'x')
drawScaleBar(10,'horizontal','position',[162 -15],'color','k','Linewidth',3);
drawScaleBar(115,'vertical','position',[162 -15],'color','k','Linewidth',3);
xlim([45 162.5]); ylim([-17 115]);
set_figsize(300,160)

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
%set(gca,'YTick',[1:34],'YTickLabel',num2str(dsort'))
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
set_figsize(270,160);

%% Figure 2D, Correlation matrix between ROIs
figure(118); clf
perispike_time=[-3:10];
Corrcoeff=[]; corrMat=[];
for f=20
    [nROI nTime]=size(PCresult{f}.NormalizedTrace_dirt);
    sub_ch=[];
    for ch=1:2
        sub_ch{ch}=get_subthreshold(PCresult{f}.NormalizedTrace_ch{ch},PCresult{f}.SpikeMat(1,:)>0,7,20);
    end
    perispike_frame=unique([tovec(find(double(PCresult{f}.SpikeMat(1,:)>0))'+perispike_time); find(PCresult{f}.CStrace)']);
    perispike_frame(perispike_frame<=0 | perispike_frame>nTime)=[];
    silenttime_vec=ind2vec(nTime,perispike_frame,0,1);
    nonvalid_frame=find(sum(isnan(PCresult{f}.Subthreshold),1)>0);

    Blue_on_frame=find(imdilate(PCresult{f}.BlueStim>0, [ones(1, 1), 1, ones(1, 200)]));
    Badframe=unique([Blue_on_frame nonvalid_frame]);
    Goodframe=setdiff([1:nTime],Badframe);
    roisD_order_ind=cellfun(@find,PCresult{f}.roisD_order,'UniformOutput',false);

    nexttile([1 1]);
    PFinframe=find(ind2vec(nTime, Goodframe,1)>0 & PCresult{f}.PFvec); PFoutframe=find(ind2vec(nTime, Goodframe,1)>0 & (1-PCresult{f}.PFvec));
    Runframe=find(ind2vec(nTime, Goodframe,1)>0 & PCresult{f}.runVec); Restframe=find(ind2vec(nTime, Goodframe,1)>0 & (1-PCresult{f}.runVec));
    silentGoodframe=find(ind2vec(nTime, Goodframe,1)>0 & silenttime_vec);
    corrMat{f,1}=get_corrMat(sub_ch{1},sub_ch{2},Goodframe);
    corrMat{f,2}=get_corrMat(sub_ch{1},sub_ch{2},silentGoodframe); %silent
    corrMatPF{f,1}=get_corrMat(sub_ch{1},sub_ch{2},PFinframe);
    corrMatPF{f,2}=get_corrMat(sub_ch{1},sub_ch{2},PFoutframe);
    corrMatRun{f,1}=get_corrMat(sub_ch{1},sub_ch{2},Runframe);
    corrMatRun{f,2}=get_corrMat(sub_ch{1},sub_ch{2},Restframe);

    cmap_label=hsv(5);
    cmap_label=cmap_label([4 2 3 5 1],:);
    distMatrix = sqrt(2*(1 - corrMat{f,1}));
    distMatrix(logical(eye(size(distMatrix)))) = 0;
    Z_ref = linkage(squareform(distMatrix), 'average');
    leafOrder = optimalleaforder(Z_ref,squareform(distMatrix));
    Cluster_ref= switchlabel(cluster(Z_ref, 'maxclust', 5));
    [~, order_cluster]=sort(Cluster_ref,'ascend');

    l=imshow_label(corrMat{f,1}(:,:),PCresult{f}.labelvec,cmap_label,{'Basal','Soma','Trunk','Oblique','Apical (distal)'});
    axis equal tight off
    if f==foi(end)
    else
        legend off
    end
    caxis([-0.3 1])
    title(num2str(f))
    drawnow;
end
colormap(gen_colormap([0 0.4 1; 1 1 1; 1 0 0],256));

f=20;
figure(125);
clf;
for j=[2 7 20 33];%1:size(Subthreshold{f},1)
    nexttile([1 1]);
    Ft=PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist);
    bd=cell2mat(bwboundaries(imdilate(Ft(:,:,j)>0,strel('diamond',5))));
    showScaleImageSmooth(permute(Ft>0,[1 2 3]), corrMat{f,2}(j,:),colormap(gen_colormap([0 0.5 1; 1 1 1; 1 0 0])),[-0.25 1]); hold all
    plot(bd(:,2),bd(:,1),'color',[1 1 1])
end
cb=colorbar;
colormap(gen_colormap([0 0.5 1; 1 1 1; 1 0 0]));
cb.Ticks=[0 1];
cb.TickLabels=[-0.4 1];

%% Auto/Cross correlogram (Figure S)
dClass2temporal = [[{1},{1}];[{[3 5]} {[3 5]}];[{[1]} {[5]}];[{[5]} {[5]}]]; % example

maxLag = 500;   % max lag in frames (change as you like)

temporalCorr = cell(length(foi), size(dClass2temporal,1));  % {f, pair}
lagsCorr      = cell(length(foi), 1);

for fi = 1:length(foi)
    f = foi(fi);

    % Use the same subthreshold + silent frames as above
    [nROI nTime] = size(PCresult{f}.NormalizedTrace_dirt);
    sub_ch = cell(1,2);
    for ch = 1:2
        sub_ch{ch} = get_subthreshold(PCresult{f}.NormalizedTrace_ch{ch}, max(PCresult{f}.SpikeMat(1,:),[],1)>0, 7, 17);
    end

    % Rebuild silentGoodframe for this loop (or reuse if still in scope)
    perispike_time  = -3:10;
    perispike_frame = unique([tovec(find(double(PCresult{f}.SpikeMat(1,:)==1))' + perispike_time); find(PCresult{f}.CStrace)']);
    perispike_frame(perispike_frame <= 0 | perispike_frame > nTime) = [];
    silenttime_vec  = ind2vec(nTime, perispike_frame, 0, 1);
    nonvalid_frame  = find(sum(isnan(PCresult{f}.Subthreshold),1) > 0);
    Blue_on_frame   = find(imdilate(PCresult{f}.BlueStim > 0, [ones(1,1), 1, ones(1,200)]));
    Badframe        = unique([Blue_on_frame nonvalid_frame]);
    Goodframe       = setdiff(1:nTime, Badframe);
    silentGoodframe = find(ind2vec(nTime, Goodframe, 1) > 0 & silenttime_vec);

    % Loop over class pairs defined in dClass2temporal
    for dc = 1:size(dClass2temporal,1)
        classA = dClass2temporal{dc,1};   % e.g. [1] or [3 5]
        classB = dClass2temporal{dc,2};   % e.g. [1] or [5]

        roiA = ismember(PCresult{f}.labelvec, classA);
        roiB = ismember(PCresult{f}.labelvec, classB);

        % Average subthreshold across ROIs in each class and across channels
        % (N_rois x T → 1 x T)
        traceA_ch1 = mean(sub_ch{1}(roiA, silentGoodframe), 1, 'omitnan');
        traceA_ch2 = mean(sub_ch{2}(roiA, silentGoodframe), 1, 'omitnan');
        traceB_ch1 = mean(sub_ch{1}(roiB, silentGoodframe), 1, 'omitnan');
        traceB_ch2 = mean(sub_ch{2}(roiB, silentGoodframe), 1, 'omitnan');

        traceA = (traceA_ch1) / 2;
        traceB = (traceB_ch2) / 2;

        % Remove NaNs and mean-center
        valid  = ~isnan(traceA) & ~isnan(traceB);
        tA     = traceA(valid) - mean(traceA(valid));
        tB     = traceB(valid) - mean(traceB(valid));

        % Normalized auto/cross-correlation
        [xc, lags] = xcorr(tA, tB, maxLag, 'coeff');

        temporalCorr{fi, dc} = xc;   % correlation vs lag
        lagsCorr{fi}         = lags; % same for all dc in this f
    end
end

param=[];
for ff=1:size(temporalCorr,1)
    for dd=1:size(temporalCorr,2)
        [~, param(ff,dd)]=expfitDM_2(abs(lagsCorr{ff})',temporalCorr{ff,dd}',[0:500]',100);
    end
end

M_b=mean(param(param(:,1)~=100,1)); S_b=std(param(param(:,1)~=100,1));
M_a=mean(param(param(:,2)~=100,2)); S_a=std(param(param(:,2)~=100,2));
fprintf(['1/e const. of autocorrelogram, basal = %2.0f ± %2.0f, apical = %2.0f ± %2.0f \n'],M_b,S_b,M_a,S_a);

figure(124); clf; tiledlayout(2,2);
cmap_temp   = gen_colormap(Plasma,4);

title_str = { ...
    'Basal, Basal (auto)', ...
    'Apical, Apical (auto)', ...
    'Basal, Distal (cross)',...
    'Distal, Distal (auto)'};

panel_order = [1 2 4 3];
tileID = 1;

for dc = panel_order
    % find sessions that have BOTH classes in this pair
    validSessions = [];
    for fi = 1:length(foi)
        f = foi(fi);
        hasA = any(ismember(PCresult{f}.labelvec, dClass2temporal{dc,1}));
        hasB = any(ismember(PCresult{f}.labelvec, dClass2temporal{dc,2}));
        if hasA && hasB
            validSessions(end+1) = fi;
        end
    end

    % skip if no valid sessions
    if isempty(validSessions), cofntinue; end

    % collect correlograms for this pair across valid sessions
    lags       = lagsCorr{1};
    Corr2show  = cell2mat(temporalCorr(validSessions, dc));  % [nLags x nSess]
    t_lag      = lags;

    nexttile(tileID); tileID = tileID + 1;

    % individual sessions in gray
    plot(t_lag, Corr2show, 'color', [0.7 0.7 0.7]); hold on;

    % mean ± std across valid sessions
    mCorr = mean(Corr2show, 1, 'omitnan');
    sCorr = std(Corr2show, 0, 1, 'omitnan');
    errorbar_shade(t_lag, mCorr, sCorr, ...
        cmap_temp(dc, :));

    xlabel('Lag (ms)');
    title(title_str{dc});
    box off;
    set(gca,'YLabel',[])
    if dc~=3
        ylim([-0.2 1])
    end
end
nexttile(1,[1 1]);
ylabel('Correlation');

%% Peak triggered average (Figure 2F–G)

PeakSTA=[]; peakvec=[]; TroughSTA=[]; troughvec=[]; PeakTrough_FR=[];
nTauPeak=[1000 1000];
PeakTroughMat=[]; PeakSTASpClass=[]; TroughSTASpClass=[];
for f=[foi]
    PeakSTA{f}=[];
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
        else
            PeakSTA{f}(:,:,n)=NaN(nROI,sum(nTauPeak)+1);
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
%% Peak triggering average kymo (Figure 2F)
figure(132); clf; tiledlayout(2,1);
showf=[23]; cax=[-0.01 0.035];
[nROI, nTime]=size(PCresult{showf}.Subthreshold);
showInd=[3 8 39];
showSTA=PeakSTA{showf}-median(PeakSTA{showf}(:,1:1000,:),2);
for i=[1 3]
    nexttile([1 1]);
    [dendaxissorted, dsort]=sort(PCresult{f}.dendaxis);
    imagesc([-nTauPeak(1):nTauPeak(2)],[1:nROI],showSTA(dsort,:,showInd(i)),cax);
    set_kymoYtick(dendaxissorted);
    xlabel('Time (ms)'); %ylabel('Distance from soma (µm)');
    xlim([-500 500]);
    cb=colorbar;
    cb.Label.String='∆F/F';
end
colormap(turbo);
set_fontsize(18); set_font('Arial');
set_figsize(120,220);

SWCpoints=zeros(size(PCresult{f}.swc,1),4);
SWCpoints(:,1:3)=PCresult{f}.swc;
SWCpoints(:,4)=SWCpoints(:,4)+6; SWCpoints(1,4)=30;
ftprints=PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist);
[H, W, T] = size(ftprints);
t2show=[-120:60:120]; shift_x=180;

figure(133); clf;
tiledlayout(2,1,'TileSpacing','tight');
for i=[1 3]
    nexttile([1 1]);
    [dendaxissorted, dsort]=sort(PCresult{f}.dendaxis);
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
        drawScaleBar(100*Pixelsize(showf),'vertical','color','k','position',[shift_x*(g-0.1) 50],'Linewidth',1.5);
    end
end
set_figsize(230,200);

%% Peak triggering average by ROIs (Figure 2G)
PeakTraceROI=[];
cmap_errorbar=[0.5 0 1; 0.3 0.3 1; 1 0.4 0.1; 1 0 0];
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
%dClass2plot=[[{1},{1}];[{[3 5]} {[3 5]}];[{[4]} {[4]}]];
dClass2plot=[[{[1 2 3 4 5]},{[1 2 3 4 5]}]]; % Trigger, Read
pair_str=[]; Length_PeakMat=[]; %neuronInd, Trigger Pair i, Trigger Label, Measure Label, Length C, R^2, Legnth C dX, R^2 dX
cmap_corrscatter=[0.5 0.8 1; 1 0.5 0.5];
cmap_neuron=hsv(max(foi));
M=[]; S=[]; M1d=[]; S1d=[]; dendbin=[0 15:45:450 500 600];
PeakDecayData_Pair=[]; PeakDecayData_Cat=[];
for f=foi
    DistMat=PCresult{f}.interDendDist(PCresult{f}.sortDist,PCresult{f}.sortDist);
    DistMat=max(cat(3,DistMat,DistMat'),[],3);
    DistMat1d = abs(PCresult{f}.dendaxis1d - PCresult{f}.dendaxis1d');
    PeakZeroMat=squeeze(PeakSTA{f}(:,nTauPeak(1)+1,:)); %row: Read, column: Trigger
    ft=PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist);
    ftcoord=get_coord(ft);

    for dc=1:size(dClass2plot,1) %region class to plot

        TrigROIs=find(ismember(PCresult{f}.labelvec,dClass2plot{dc,1}));
        ReadROIs=find(ismember(PCresult{f}.labelvec,dClass2plot{dc,2}));

        for tr=TrigROIs
            dax=DistMat(ReadROIs,tr);
            dax1d=DistMat1d(ReadROIs,tr);
            Volt=PeakZeroMat(ReadROIs,tr);
            [~, trInd_inReadROI]=min(dax);            
            [grad_mag, grad_dir, grad_vec] = get_gradient(Volt, ftcoord(ReadROIs,:), 100);
            PeakDecayData_Cat=[PeakDecayData_Cat; {f, dClass2plot{dc,1}, dClass2plot{dc,2}, tr, ReadROIs, dax, dax1d, Volt, grad_mag(trInd_inReadROI),grad_dir(trInd_inReadROI), grad_vec(trInd_inReadROI,:), ftcoord(tr,:)}];
        end

        uppertri_ind=triu(ones(size(DistMat)))';
        label_ind=ismember(PCresult{f}.labelvec,dClass2plot{dc,2})'*ismember(PCresult{f}.labelvec,dClass2plot{dc,1});
        daxpair=tovec(DistMat(uppertri_ind & label_ind));
        dax1dpair=tovec(DistMat1d(uppertri_ind & label_ind));
        peakdat2fit=tovec(PeakZeroMat(uppertri_ind & label_ind));
        peakVoltage=mean(peakdat2fit(daxpair==0),'omitnan');
        N_points=sum(~isnan(peakdat2fit));
        prc_triggeringpoint=sum(peakdat2fit<peakVoltage)./N_points;

        if N_points>5
            [fitted_exp params R2]=expfit_wBd(daxpair,peakdat2fit,[0:500],[peakVoltage -500],[peakVoltage*0.7 -max(daxpair)*2],[peakVoltage*1.2 -100]);
            [fitted_exp1d params1d R21d]=expfit_wBd(dax1dpair,peakdat2fit,[0:500],[peakVoltage -500],[peakVoltage*0.7 -max(dax1dpair)*2],[peakVoltage*1.5 -100]);
            PeakDecayData_Pair=[PeakDecayData_Pair; {f, dClass2plot{dc,1}, dClass2plot{dc,2}, TrigROIs, ReadROIs, -params(2), R2, -params1d(2), R21d, N_points, ...
                prc_triggeringpoint,daxpair, dax1dpair,peakdat2fit, fitted_exp, fitted_exp1d}];
        end
    end
end

fieldName={'NeuronID','TriggerClass','ReadClass','TriggerROI','ReadROIs','dax','dax1d','TriggerVoltage','Gradient_mag','Gradient_dir','Gradient_vec','Trigger_coord'};
PeakDecayData_Cat=array2table(PeakDecayData_Cat,'VariableNames',fieldName);

fieldName={'NeuronID','TriggerClass','ReadClass','TriggerROI','ReadROIs','LengthConstant','Rsquare','LengthConstant1d','Rsquare1d','NFitpoints','Percentile','dax','dax1d','TriggerVoltage','Expfit','Expfit1d'};
PeakDecayData_Pair=array2table(PeakDecayData_Pair,'VariableNames',fieldName);
%%
figure(131); clf; tiledlayout(1,3);
f2show=18;
showInd=cell2mat(PeakDecayData_Cat.NeuronID)==f2show;
nexttile([1 1]);
for j=find(showInd)'
    d2show=PeakDecayData_Cat.dax{j};
    v2show=PeakDecayData_Cat.TriggerVoltage{j};
    [d2show, dsort]=sort(d2show);
    v2show=v2show(dsort);
    %v2show=v2show./v2show(1);
    plot(d2show,v2show,'color',[0.5 0.5 0.5]); hold all;
end

nexttile([1 1]);
for j=find(showInd)'
    d2show=PeakDecayData_Cat.dax1d{j};
    v2show=PeakDecayData_Cat.TriggerVoltage{j};
    [d2show, dsort]=sort(d2show);
    v2show=v2show(dsort);
    %v2show=v2show./v2show(1);
    plot(d2show,v2show,'color',[0.5 0.5 0.5]); hold all;
end

%%
f2show=20;
ftcoord=get_coord(PCresult{f2show}.Ftprnts(:,:,PCresult{f2show}.sortDist));
showInd=find(cell2mat(PeakDecayData_Cat.NeuronID)==f2show);
figure(142); clf;
for j=showInd(1:3:30)'
    nexttile([1 1]);
Coord2show=PeakDecayData_Cat.Trigger_coord{j};
Gradvec2show=PeakDecayData_Cat.Gradient_vec{j}*700000;
Volt2show=PeakDecayData_Cat.TriggerVoltage{j};
% Arrow plot — scale arrow length by magnitude
imshow2(PCresult{f2show}.ref_im,[]); hold all;
scatter(ftcoord(PeakDecayData_Cat.ReadROIs{j},1), ftcoord(PeakDecayData_Cat.ReadROIs{j},2), 40, vec2cmap(Volt2show,'turbo'),'filled')
quiver(Coord2show(:,1), Coord2show(:,2), Gradvec2show(:,1), Gradvec2show(:,2), 0.5, 'w', 'LineWidth', 1);
end



%%

f=26;
DistMat=PCresult{f}.interDendDist(PCresult{f}.sortDist,PCresult{f}.sortDist);
DistMat=max(cat(3,DistMat,DistMat'),[],3);
DistMat1d = abs(PCresult{f}.dendaxis1d - PCresult{f}.dendaxis1d');
PeakZeroMat=squeeze(PeakSTA{f}(:,nTauPeak(1)+1,:)); %row: Read, column: Trigger
ft=PCresult{f}.Ftprnts(:,:,PCresult{f}.sortDist);
ftcoord=get_coord(ft);
dc=1;
TrigROIs=find(ismember(PCresult{f}.labelvec,dClass2plot{dc,1}));
ReadROIs=find(ismember(PCresult{f}.labelvec,dClass2plot{dc,2}));
tr=1;
dax=DistMat(ReadROIs,tr);
dax1d=DistMat1d(ReadROIs,tr);
Volt=PeakZeroMat(ReadROIs,tr);

x=ftcoord(ReadROIs,1); y=ftcoord(ReadROIs,2);
d = [0; cumsum(sqrt(diff(x).^2 + diff(y).^2))];

f2show=26;
figure;
imshow2(PCresult{f2show}.ref_im,[]); hold all;
scatter(x, y, 40, vec2cmap(Volt,'turbo'),'filled')

% Activity at one time point
A = Volt;

% Create regular grid
[xq,yq] = meshgrid(linspace(min(x),max(x),100), ...
                   linspace(min(y),max(y),100));

% Interpolate activity
Agrid = griddata(x,y,A,xq,yq,'natural');

% Compute gradient
[Ax,Ay] = gradient(Agrid);

% Gradient magnitude
grad_mag = sqrt(Ax.^2 + Ay.^2);

% Plot
figure
imagesc(grad_mag)
colorbar




%%
%linkaxes([ax1, ax2]);
xlim([0 400])
%NeuronInd, Pair class, Length Cst, R^2, Legnth Cst dX, R^2 dX, N fit points, range of distance, range dx, peak zscore


% Plot local stimulation data
LocalStimResult=importdata(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism/LocalStimulationResult_20260609.mat']);
t2show=[-5:5:60];
% tInd=Lstaall.OnsetTime>=20 & Lstaall.OnsetTime<=40;
% path2show=find(Lstaall.StimROI_percentile>0.5 & Lstaall.Fit_R2>0.52 & tInd & Lstaall.Fit_length<Lstaall.Max_distance*0.5); % Percentile, R^2
% Stim2show=unique(Lstaall.StimInd(path2show));
% Lsta2show=ismember(Lstaall.StimInd,Stim2show) & Lstaall.Fit_length<240 & Lstaall.Fit_length>20;

tInd2show=find(t2show>=20 & t2show<=40);
tInd=LocalStimResult.Lstaall.OnsetTime>=20 & LocalStimResult.Lstaall.OnsetTime<=40;
path2show=find(LocalStimResult.Lstaall.StimROI_percentile>0.5 & LocalStimResult.Lstaall.Fit_R2>0.52 & tInd); % Percentile, R^2
Stim2show=unique(LocalStimResult.Lstaall.StimInd(path2show));
Lsta2show=ismember(LocalStimResult.Lstaall.StimInd,Stim2show) & LocalStimResult.Lstaall.Fit_length<240 & LocalStimResult.Lstaall.Fit_length>20;

data2show=LocalStimResult.STADistmat(Stim2show,tInd2show);
[deltax, dsort]=cellfun(@(x) sort(x(:,1)),data2show(:,1),'UniformOutput',false);
[~, dsort1d]=cellfun(@(x) sort(x(:,3)),data2show1d(:,1),'UniformOutput',false);
data2show_ave=[]; data2show_ave1d=[];
volt2show=cellfun(@(x) x(:,2),data2show,'UniformOutput',false);
for ss=1:size(data2show,1)
data2show_ave{ss,1}=[data2show{ss,1}(:,1) mean(cell2mat(volt2show(ss,:)),2) data2show{ss,1}(:,3) ];
end

data2show_ave2=cellfun(@(x,y) [x(y,1) x(y,2)/x(x(:,1)==0,2)],data2show_ave,dsort,'UniformOutput',false); %normalize
data2show_ave1d=cellfun(@(x,y) [(x(y,3)) x(y,2)/x(x(:,1)==0,2)],data2show_ave,dsort1d,'UniformOutput',false);
LocalstimVoltage_bin=binning_data(data2show_ave2,dendbin);

l(4)=errorbar_shade(LocalstimVoltage_bin.centers,LocalstimVoltage_bin.mean,LocalstimVoltage_bin.sem,cmap_errorbar(4,:));
xlim([0 200]); box off;

nexttile([1 2])
cmap=[0 0.6 1;0 0 0; 0.7 0 0.9];
%SponLength2show=find(Length_PeakMat.Rsquare>0.2 & Length_PeakMat.peak_zscore>1);
SponLength2show=find(Length_PeakMat.MaxDistance>0);
LengthC2show={LocalStimResult.Lstaall.Fit_length(path2show), Length_PeakMat.LengthConstant(SponLength2show)};
p=Boxplot_wPoints2(LengthC2show,cmap([1 2],:)); axis tight; hold all
drawPValueLines(p,300,'TextYOffset',0.1,'StepHeight',0.1)
set(gca,'yscale','log','xtick',[1 2],'xticklabel',{'Opto. Stim.','Spontaneous'},'ytick',[10 100 1000])
ylabel('Decay length (\mum)'); ylim([10 1700]); xlim([0.5 2.5]); box off;
