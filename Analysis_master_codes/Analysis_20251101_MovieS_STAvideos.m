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
nTauDelay=[-3:4]; % Short time window to estimate delay.
nTauAlign=[20 15];
dclass_str={'Basal','Soma','Trunk','Oblique','Distal'};
dSpikePropsMat=[]; FiringRate=[];
nTauSTA=[150 150];

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
    CS_trace_new(CS_fistSp-1)=0;    %do not connect two CSs

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

    PCresult{f}.BlueStim = Result.Blue;
    PCresult{f}.CStrace = CS_trace;
    PCresult{f}.BStrace = BS_trace;
    PCresult{f}.avgImg = Result.ref_im;
    PCresult{f}.Ftprnts = Result.ftprnt;
    PCresult{f}.VR = Result.VR;

    ActivityTr=max([Result.spike(1,:); BS_trace; CS_trace],[],1)>0;
    PCresult{f}.Subthreshold_larger = get_subthreshold(NormalizedTrace_dirt, ActivityTr, 40, 500);
    PCresult{f}.Subthreshold_larger(isnan(NormalizedTrace_dirt)) = NaN;

    %dSpikes
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
%% Save STA movies
time_segment=15000; bound=5; overlap=200; Dirtthreshold=0.3;
for f=[4 6 11 10 18 20 23]
    f
    cd(fpath{f})
    SpikeLabelMat=[]; MovSTAinfo=[];
    [nROI nTime]=size(PCresult{f}.Subthreshold);
    load(fullfile(fpath{f},'PC_Result.mat'),'Result');
    load([fpath{f} '/output_data.mat'])
    load([fpath{f} '/MovSTAinfo.mat'])
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
    Mov_SpikeTA=zeros(sz(2)*sz(1),sum(nTauSTA)+1,4);
    Movfilt_SpikeTA=zeros(sz(2)*sz(1),sum(nTauSTA)+1,4);
    N_added=zeros(sz(2)*sz(1),sum(nTauSTA)+1,4);
    SpClassVec=PCresult{f}.SpikeClassMat;

    DendOnlySpikes=find(dSpikePropsMat{f}.Branch_Index~=1 & ismember(dSpikePropsMat{f}.dClass,[3 5]) &...
        ~dSpikePropsMat{f}.IsbAP & ~dSpikePropsMat{f}.IsCS & dSpikePropsMat{f}.dt_pre_somAP>=4 & dSpikePropsMat{f}.dt_post_somAP>=2.5);
    dSpikeFrame=unique(dSpikePropsMat{f}.Spike_frame(DendOnlySpikes));
    dSpikeFrame(find(diff(dSpikeFrame)==1)+1)=[];
    ApicaldSpikeFrame=setdiff(dSpikeFrame,FiringRate{f}.DendIntSpFrame);
    SpClassVec(4,:)=ind2vec(nTime,ApicaldSpikeFrame,1);

    for j=1:length(f_seg)-1

        [nInd fInd]=find(SpClassVec(:,[f_seg_real(j):f_seg_real(j+1)]));

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
                mov_res_sub=toimg(get_subthreshold(tovec(mov_res),ind2vec(size(mov_res,3),fInd,1),11,25),sz(2),sz(1));
                mov_res_sub=maskEdge(mov_res_sub,5,0);
                mov_STA_tmp=get_STA(tovec(mov_res),fInd,3,3);
                mov_STA_tmp=toimg(mov_STA_tmp,sz(2),sz(1));
                %[~,~,Vsub]=get_eigvector(tovec(imgaussfilt(mov_res_sub,1.5))',10);
                %[~,~,Vsub]=get_eigvector(tovec(mov_res_sub)',20);
                [~,~,Vsub]=get_eigvector(tovec(imgaussfilt(mov_res_sub(:,:,1:end),1))',20);
                imshow2_patch(toimg(Vsub,sz(2),sz(1)))
                sub2keep=input('sub component to keep \n');
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
    MovTAname=['/SpikeTriggeredMatMovie_Type.bin'];
    Mov_SpikeTA.transpose.savebin([fpath{f} MovTAname])

    Movfilt_SpikeTA=-Movfilt_SpikeTA./N_added;
    Range_MovFiltSTA=[min(Movfilt_SpikeTA(:)) max(Movfilt_SpikeTA(:))];
    Movfilt_SpikeTA=vm(mat2gray(double(Movfilt_SpikeTA))*10000);
    MovTAfiltname=['/SpikeTriggeredMatMovieFilt_Type.bin'];
    Movfilt_SpikeTA.transpose.savebin([fpath{f} MovTAfiltname])

    fieldName={'StackedFrame','StackedStype'};
    SpikeLabelMat=array2table(SpikeLabelMat,'VariableNames',fieldName);
    MovSTAinfo.StackedSpike=SpikeLabelMat;
    MovSTAinfo.Range_MovFiltSTA=Range_MovFiltSTA;
    MovSTAinfo.Range_MovSTA=Range_MovSTA;
    MovSTAinfo.N_added=N_added;
    MovSTAinfo.template_basis=template_basis;
    save([fpath{f} '/MovSTAinfo.mat'],'MovSTAinfo','-v7.3')
end
%% Convert to dF/F and color
for f=[4 6 10 11 18 20 23]
    cd(fpath{f}); F0img=[];
    load(fullfile(fpath{f},'PC_Result.mat'),'Result');
    load([fpath{f} '/output_data.mat'])
    load([fpath{f} '/MovSTAinfo.mat'])
    load(fullfile(fpath{f},'coloredSTA.mat'));

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

    j=2;
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
    MovSTA=double(readBinMov('SpikeTriggeredMatMovie_Type.bin',sz(2)*sz(1),301));
    MovSTA=reshape(MovSTA,sz(2),sz(1),301,[]);
    MovSTAfilt=double(readBinMov('SpikeTriggeredMatMovieFilt_Type.bin',sz(2)*sz(1),301));
    MovSTAfilt=reshape(MovSTAfilt,sz(2),sz(1),301,[]);

    oldRange = [min(MovSTA(:)) max(MovSTA(:))];
    NewRange = [min(MovSTAinfo.Range_MovSTA(:)) max(MovSTAinfo.Range_MovSTA(:))];
    MovSTA = (MovSTA - oldRange(1)) * diff(NewRange) / diff(oldRange) + NewRange(1);

    oldRange = [min(MovSTAfilt(:)) max(MovSTAfilt(:))];
    NewRange = [min(MovSTAinfo.Range_MovFiltSTA(:)) max(MovSTAinfo.Range_MovFiltSTA(:))];
    MovSTAfilt = (MovSTAfilt - oldRange(1)) * diff(NewRange) / diff(oldRange) + NewRange(1);

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
    strImg_bin2=strImg_bin2(:,:,1); strImg_bin2(strImg_bin2<0.05)=0;
    ax1=[ax1 nexttile([1 1])];
    imshow2(strImg_bin2,[])
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
    dFFmov_filt=MovSTAfilt./F0img.*strImg_bin.*strImg_bin2;
    dFFmov=MovSTA./F0img.*strImg_bin.*strImg_bin2;;
    dFFmov_tmp=dFFmov; dFFmov_tmp(isnan(dFFmov_tmp) | isinf(dFFmov_tmp))=0;
    SpikeHeight_dFFmov=tovec(squeeze(dFFmov_tmp(:,:,nTauSTA(1)+1,:)))'*tovec(Result.ftprnt(:,:,1)>0)./sum(tovec(Result.ftprnt(:,:,1)>0));
    dFFmov=dFFmov./mean(SpikeHeight_dFFmov(1:3));
    dFF_range=[-0.3 1.5];

    dFFmov2show=[];
    dFFfiltmov2show=[];
    colored_dFFmov=[];
    colored_dFFfiltmov=[];
    cmap=gen_colormap([gen_colormap([0 0.2 1; 0 0 0; 1 0 0],5); gen_colormap([1 0 0; 1 1 0],6)],256);

    tic;
    for stype=1:4
        for t=1:size(dFFmov,3);
            dFF_tmp=dFFmov(:,:,t,stype);
            dFF_tmp(max(Result.bvMask,[],3)>0)=NaN;
            dFF_tmp=maskEdge(dFF_tmp,6,0);
            dFF_tmp=medfilt2nan(dFF_tmp,[8 8]);
            dFFmov2show(:,:,t,stype) = imgaussfiltnan(dFF_tmp, 2).*strImg_bin(:,:,1);
            dFFmov2show(:,:,t,stype)= fillNaNmovie(dFFmov2show(:,:,t,stype));

            dFF_tmp=dFFmov_filt(:,:,t,stype);
            dFF_tmp(max(Result.bvMask,[],3)>0)=NaN;
            dFF_tmp=maskEdge(dFF_tmp,6,0);
            dFF_tmp=medfilt2nan(dFF_tmp,[8 8]);
            dFFfiltmov2show(:,:,t,stype) = imgaussfiltnan(dFF_tmp, 2).*strImg_bin(:,:,1);
            dFFfiltmov2show(:,:,t,stype)= fillNaNmovie(dFFfiltmov2show(:,:,t,stype));

            colored_dFFmov(:,:,:,t,stype) = grs2rgb(double(dFFmov2show(:,:,t,stype)), cmap ,dFF_range(1),dFF_range(2)).*strImg_bin(:,:,1).*strImg_bin2;
            colored_dFFmov(:,:,:,t,stype) = grs2rgb(double(strImg), colormap("gray")).*strImg_bin + colored_dFFmov(:,:,:,t,stype);

            colored_dFFfiltmov(:,:,:,t,stype) = grs2rgb(double(dFFfiltmov2show(:,:,t,stype)), cmap ,dFF_range(1),dFF_range(2)).*strImg_bin(:,:,1).*strImg_bin2;
            colored_dFFfiltmov(:,:,:,t,stype) = grs2rgb(double(strImg), colormap("gray")).*strImg_bin + colored_dFFfiltmov(:,:,:,t,stype);
        end
    end
    toc;

    save(fullfile(fpath{f},'coloredSTA.mat'),'colored_dFFmov','colored_dFFfiltmov','dFFmov2show','dFFfiltmov2show','F0img','-v7.3');
end
%% Adjust color axis
dFF_range=[-0.5 2];
cmap=gen_colormap([gen_colormap([0 0.2 1; 0 0 0; 1 0 0],5); gen_colormap([1 0 0; 1 1 0],6)],256);
for f=[4 10 11 18 20 23]
    load([fpath{f} '/coloredSTA.mat'])
    load(fullfile(fpath{f},'PC_Result.mat'),'Result');

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
    strImg_bin2=strImg_bin2(:,:,1); strImg_bin2(strImg_bin2<0.05)=0;
    ax1=[ax1 nexttile([1 1])];
    imshow2(strImg_bin2,[])
    linkaxes(ax1)

    if isfield(Result,'Structure')
        strImg=Result.Structure.*point2img(Result.SWC(:,[1 2]),7,size(Result.ref_im));
        strImg_bin=strImg>0;
    else
        strImg=avgImg;
        strImg_bin=strImg_bin2;
    end

    for stype=1:4
        for t=1:size(dFFmov2show,3)
            colored_dFFmov(:,:,:,t,stype) = grs2rgb(double(dFFmov2show(:,:,t,stype)), cmap ,dFF_range(1),dFF_range(2)).*strImg_bin(:,:,1).*strImg_bin2;
            colored_dFFmov(:,:,:,t,stype) = grs2rgb(double(strImg), colormap("gray")).*strImg_bin + colored_dFFmov(:,:,:,t,stype);

            colored_dFFfiltmov(:,:,:,t,stype) = grs2rgb(double(dFFfiltmov2show(:,:,t,stype)), cmap ,dFF_range(1),dFF_range(2)).*strImg_bin(:,:,1).*strImg_bin2;
            colored_dFFfiltmov(:,:,:,t,stype) = grs2rgb(double(strImg), colormap("gray")).*strImg_bin + colored_dFFfiltmov(:,:,:,t,stype);
        end
    end
    save(fullfile(fpath{f},'coloredSTA.mat'),'colored_dFFmov','colored_dFFfiltmov','dFFmov2show','dFFfiltmov2show','F0img','-v7.3');
    disp([num2str(f) ' is done.'])
end

%% Generate Movie
%[4 6 10 11 18 20 23]
Neuron2show=[4 10 11 18 20 23];
%Neuron2show=[4];
nTauSTA=[150 150];
colored_dFFmov2show=[]; STAsKymos2show=[]; STAsKymo2show_line=[];
cmap2=gen_colormap([gen_colormap([0 0.2 1; 1 1 1; 1 0 0],5); gen_colormap([1 0 0; 1 1 0],6)],256);
for f=Neuron2show
    load(fullfile(fpath{f},'coloredSTA.mat'));
    load([fpath{f} '/MovSTAinfo.mat'])
    [nROI nTime]=size(PCresult{f}.NormalizedTrace_dirt);

    SpClassVec=PCresult{f}.SpikeClassMat;
    DendOnlySpikes=find(dSpikePropsMat{f}.Branch_Index~=1 & ismember(dSpikePropsMat{f}.dClass,[3 5]) &...
        ~dSpikePropsMat{f}.IsbAP & ~dSpikePropsMat{f}.IsCS & dSpikePropsMat{f}.dt_pre_somAP>=4 & dSpikePropsMat{f}.dt_post_somAP>=2.5);
    dSpikeFrame=unique(dSpikePropsMat{f}.Spike_frame(DendOnlySpikes));
    dSpikeFrame(find(diff(dSpikeFrame)==1)+1)=[];
    ApicaldSpikeFrame=setdiff(dSpikeFrame,FiringRate{f}.DendIntSpFrame);
    SpClassVec(4,:)=ind2vec(nTime,ApicaldSpikeFrame,1);

    colored_dFFmov2show{f}=colored_dFFmov;
    for sptype=1:4
        sptime=find(SpClassVec(sptype,:));
        [STAsKymos2show{f}(:,:,sptype)]=get_STA(PCresult{f}.NormalizedTrace_dirt,sptime,nTauSTA(1),nTauSTA(2));
        STAsKymos2show{f}(:,:,sptype)=STAsKymos2show{f}(:,:,sptype)./PCresult{f}.SpikeHeight;
        if ~isempty(sptime)
            STAsKymos2show{f}(:,:,sptype)=STAsKymos2show{f}(:,:,sptype)-prctile(STAsKymos2show{f}(:,1:nTauSTA(1)/2,sptype),50,2);
        end

        for dclass=[1 2 3 4 5]
            dclassroi=PCresult{f}.roi_dClass==dclass;
            if ~isempty(find(dclassroi))
                STAsKymo2show_line{f}(dclass,:,sptype)=mean(STAsKymos2show{f}(dclassroi,:,sptype),1,'omitnan');
            else
                STAsKymo2show_line{f}(dclass,:,sptype)=NaN(1,sum(nTauSTA)+1);
            end
        end
    end

end

%%
figure(20); clf;
saveMov2='/Volumes/cohen_lab/Lab/Papers/2025 Voltron Optopatch prism dendrites in vivo/MovieSs';
time2show=[-15:20]+nTauSTA(1);
stypes={'Simple spike','Complex spike',['Spike train without plateau'],'dSpike'};
v = VideoWriter([saveMov2 '/MovieS8_STAmovies'],'MPEG-4');
v.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
v.Quality= 100;
open(v);
xlimend=[400 400 300]; cmapROI=gen_colormap(Plasma,4);
for j = time2show
    clf;
    set(gca,'units','pixels','position',[200 0 1000 800]);
    tiledlayout(3*2,3,'Padding','none', 'TileSpacing', 'none');
    %tiledlayout(2,4,'Padding','compact');
    ax2=[]; g=1;
    for f=Neuron2show([1 3 4])
        s1=1;
        for stype=[1 2 4]
            ax1=nexttile(6*(g-1)+s1,[1 1]);
            if f<20
                imshow2(colored_dFFmov2show{f}(:,1:xlimend(g),:,j,stype),[0 1]);
            else
                imshow2(fliplr(colored_dFFmov2show{f}(:,1:xlimend(g),:,j,stype)),[0 1]);
            end
            pbaspect([size(double(colored_dFFmov2show{f}(:,:,:,j,stype)),2) size(double(colored_dFFmov2show{f}(:,:,:,j,stype)),1) 1]),colormap(gray); hold all;
            colormap(ax1,cmap2)
            title(stypes{stype})
            axis off
            if stype==1
                text(7,12,[num2str((j-nTauSTA(1))) ' ms'], 'FontSize', 10, 'color', [0.99 0.99 0.99])
            end

            ax2=[ax2 nexttile(6*(g-1)+3+s1,[1 1])];
            l=plot(time2show-nTauSTA(1),STAsKymo2show_line{f}([1 2 3 5],j+[-15:20],stype),'linewidth',1.5); hold all
            arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmapROI,2));
            plot(repmat(0,1,2),[-0.3 1],'k--'); axis off;
            ylim([-0.6 1.2])
            s1=s1+1;
        end
        arrayfun(@(x) colormap(x,turbo),ax1);
        nexttile(6*(g-1)+s1-1,[1 1]);
        cb=colorbar;  cb.Label.String='Z score';
        cb.Ticks=[0 1]; cb.TickLabels=[-0.3 1.5];
        drawScaleBar(100/1.17,'horizontal','color',[1 1 1],'Linewidth',3);
        text(xlimend(g)*0.7,145,'100 \mum','color','w','Fontsize',10);
        g=g+1;
    end
    nexttile(3,[1 1]); colormap(cmap2);
    nexttile(9,[1 1]); colormap(cmap2);
    arrayfun(@(x) colormap(x,cmap2),ax1);
    nexttile(18,[1 1]);
    drawScaleBar(10,'horizontal','color',[0 0 0],'Linewidth',3,'Position',[19 0.9]);
    drawScaleBar(1,'vertical','color',[0 0 0],'Linewidth',3,'Position',[19 -0.1]);
    text(11,1.1,'10 ms','color','k','Fontsize',8);
    text(20.5,-0.3,'Somatic spike height','color','k','Fontsize',8,'Rotation',90);
    
    set_fontsize(11);
    nexttile(4,[1 1]); legend({'Basal','Soma','Trunk','Distal'},"Box","off",'NumColumns',2,'Location','northwest','fontsize',7)

    pause(0.1);
    set(gcf,'color','w')    % Sets background to white
    frame = getframe(gcf);
    writeVideo(v,frame);
    pause(0.1);
end
close(v);

