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
%%
ISItotal=[]; ISIfile=[]; ISIburstorder=[]; ISIburstorderfile=[];
for f=foi
    ISItotal=[ISItotal diff(find(PCresult{f}.SpikeMat(1,:)>0))];
    ISIfile =[ISIfile repmat(f,1,numel(diff(find(PCresult{f}.SpikeMat(1,:)))))];
    [~, SporderMat_CS cs_frame]=get_STA(PCresult{f}.SpikeOrder,PCresult{f}.SpikeClassMat(2,:),1,1000);
    for bs=1:5
        ISIburstorder{f,bs}=[]; ISIburstorderfile{f,bs}=[];
        for cs=1:size(SporderMat_CS,2)
            if any(SporderMat_CS(1,cs,:)==(bs+1))
                ISIburstorder{f,bs}=[ISIburstorder{f,bs}; find(SporderMat_CS(1,cs,:)==(bs+1),1)-find(SporderMat_CS(1,cs,:)==bs,1)];
                ISIburstorderfile{f,bs}=[ISIburstorderfile{f,bs}; [cs_frame(cs)]];
            end
        end
    end
end

figure(300); clf;
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile([1 1]);
histogram(ISItotal,[1:100]); box off;
xlabel('Inter-spike interval (ms)'); ylabel('Counts');
title('ISI distribution of total SomAP')

countstr=counting_string(1:5);
for bs=1
nexttile([1 1]);
histogram(cell2mat(ISIburstorder(foi,bs)),[1:100]); box off;
xlabel('Inter-spike interval (ms)'); ylabel('Counts');
title(['ISI distribution from ' countstr{bs} ' to ' countstr{bs+1} ' spike'])
end
%%
figure(301); clf;
f2show=[ [18 2]; [18 5]; [20 4]; [23 2] ; [23 3]; [25 11]];
for show_idx=1:size(f2show,1)
    [orderidx]=find(ISIburstorder{f2show(show_idx,1),1}>16);
    [~, CS2show]=get_STA(PCresult{f2show(show_idx,1)}.NormalizedTrace_dirt,ISIburstorderfile{f2show(show_idx,1),1}(orderidx(f2show(show_idx,2))),150,100);
    [dax_sorted, sorted_dax_idx]=sort(PCresult{f2show(show_idx,1)}.dendaxis,'ascend');
    nexttile([1 1]);
    CS2show=CS2show-median(CS2show(:,:,80:100),3);
imagesc([-50:100],[1:size(CS2show,1)],permute(CS2show(sorted_dax_idx,:,100:end),[1 3 2]))
set_kymoYtick(dax_sorted);     
colormap(turbo);
end
