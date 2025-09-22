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
%%
nTau={[-200:50],[-200:50],[-200:50]}; %SS, CS, Brst
clear SpikeInd MatSpike MatSTA MatBlue MatCStrace MatSub SpikeList NormalizedTrace_ch NormalizedTrace_dirt SpikeIndBlueOff Dist_order allSpikeMat noi interDendDist noi_dist derivSub LapSubSilent
clear Subthreshold dendaxis BrstOrder roisD roisD_order LapSpclassVec
clear SI MI

for f=foi(1:end)
    f
    load(fullfile(fpath{f},'PC_Result.mat'),'Result') %load traces
    if ~isempty(find(ismember(Result.spike(1,:),Result.SpClass(1,:))==0)) | ~isempty(find(ismember(Result.spike(1,:),Result.SpClass(2,:))==0))
        error('Spclass and spike mat is not matching');
    end

    Result.spike=Result.spike>0; Result.SpClass=Result.SpClass>0;
    Dist_order{f}=Result.BrancDist_order;
    interDendDist{f}=Result.interDendDist*Pixelsize(f);
    Dsign=ones(1,size(interDendDist{f},1));
    Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
    perisomaROI=setdiff(find(interDendDist{f}(1,:)<40),BadROI{f}); % ROIs < 40 um from soma
    noi=setdiff([1:size(Result.ftprnt,3)],BadROI{f});
    noi_dist{f}=ismember(Dist_order{f},noi);
    roisD(f,:)={basal_ROI{f},PeriSoma_ROI{f},apical_ROI{f},oblique_ROI{f},distal_ROI{f}}; %set the ROIs
    for dClass=1:size(roisD,2)
        g=1;
        if ~isnan(roisD{f,dClass})
            for d=roisD{f,dClass}
                dind=setdiff(find(Result.BranchLabel==d),BadROI{f});
                roisD_order{f}{dClass,g}=ismember(Dist_order{f}(noi_dist{f}),dind)*d;
                g=g+1;
            end
        end
    end

    dendaxis{f}=interDendDist{f}(1,:).*Dsign;
    dendaxis{f}=dendaxis{f}(Dist_order{f}(noi_dist{f}));

    NormalizedTrace=(Result.normTraces)./Result.F0_PCA;
    bAP_STA=get_STA(NormalizedTrace, Result.spike(1,:), 30, 20);
    bAP_STA=bAP_STA-prctile(bAP_STA,20,2);
    SpikeHeight=max(mean(bAP_STA(perisomaROI,:),1),[],2);

    NormalizedTrace=NormalizedTrace/SpikeHeight;
    NormalizedTrace_dirt{f,1}=NormalizedTrace;
    NormalizedTrace_dirt{f,1}(:,Result.motionReject>0)=NaN;
    NormalizedTrace_ch(f,:)=cellfun(@(x) x./Result.F0_PCA,Result.norm_trace_check,'UniformOutput',false);
    NormalizedTrace_ch{f,1}(:,Result.motionReject>0)=NaN; NormalizedTrace_ch{f,2}(:,Result.motionReject>0)=NaN;
    if ifdirtRemov(f)
        NormalizedTrace_dirt{f,1}(Result.dirtTrace>0)=NaN;
        NormalizedTrace_ch{f,1}(Result.dirtTrace>0)=NaN; NormalizedTrace_ch{f,2}(Result.dirtTrace>0)=NaN;
    end
    NormalizedTrace_dirt{f,1}=NormalizedTrace_dirt{f,1}(Dist_order{f}(noi_dist{f},:),:);
    NormalizedTrace_ch{f,1}=NormalizedTrace_ch{f,1}(Dist_order{f}(noi_dist{f},:),:);
    NormalizedTrace_ch{f,2}=NormalizedTrace_ch{f,2}(Dist_order{f}(noi_dist{f},:),:);

    Subthreshold{f}=get_subthreshold(NormalizedTrace_dirt{f},max(Result.spike(1,:),[],1)>0,7,17);
    Subthreshold{f}(isnan(NormalizedTrace_dirt{f}))=NaN;
    [nROI nTime]=size(Subthreshold{f});

    ShuffleV.Pvalue{f}=Result.ShuffleV_percentile;
    ShuffleV.Zscore{f}=Result.ShuffleV_zscore;
    ShuffleV.PCPvalue{f}=Result.ShuffleEigTr_percentile;
    ShuffleV.PCZscore{f}=Result.ShuffleEigTr_zscore;
    ShuffleV.FRPvalue{f}=Result.ShuffleFR_percentile;
    ShuffleV.FRZscore{f}=Result.ShuffleFR_zscore;

    SI{f}.FRreal=Result.SI_FRreal;
    SI{f}.FRnull=Result.SI_FRnull;
    SI{f}.FRclassReal=Result.SI_FRClassReal;
    SI{f}.FRclassNull=Result.SI_FRClassnull;

    MI{f}.EigTrPos=Result.MI_EigTrReal;
    MI{f}.EigTrFR=Result.MI_EigTrFR;
    MI{f}.EigTrFRdClass=Result.MI_EigTrFRdClass;
    MI{f}.SubVFR=Result.MI_SubVFR;
    MI{f}.SubVFRdClass=Result.MI_SubVFRdClass;
    MI{f}.SubVPos=Result.MI_SubVReal;

    % spike class, SS:1, CS:2, dSP:3, BS:4
    ss_time=find(Result.SpClass(1,:)); % BS is subclass of SS
    brst=bwlabel((ss_time(2:end)-ss_time(1:end-1))<=20); % SSs that have an ISI shorter than 20 ms are BS.
    SpClass=Result.SpClass; BS_trace=zeros(1,size(Result.traces,2));
    for b=1:max(brst)
        bwn=find(brst==b);
        SpClass(1,ss_time([bwn bwn(end)+1]))=0;
        SpClass(4,ss_time([bwn(1)]))=1;
        BS_trace(1,[ss_time(bwn(1)): ss_time(bwn(end)+1)])=b;
    end
    SpClass=SpClass([1 2 4 3],:);
    Classvec = get_Class2index(SpClass);
    SpikeClassvec=Classvec.*Result.spike(1,:);

    BrstOrder = get_BurstOrder(Result.spike(1,:), 20) - SpClass(1,:);
    BrstOrder(find(SpClass(3,:)))=1;
    ComplexSpikeOrder=get_spikeOrder(Result.CStrace,Result.spike(1,:));

    for stype=1:3
        s_list=find(SpClass(stype,:)>0);
        [~, MatSpClass{f,stype}]=get_STA(SpikeClassvec,SpClass(stype,:),-nTau{stype}(1),nTau{stype}(end));
        [~, MatSpike{f,stype}]=get_STA(Result.spike(1,:),SpClass(stype,:),-nTau{stype}(1),nTau{stype}(end));
        [~, MatBlue{f,stype}]=get_STA(Result.Blue,SpClass(stype,:),-nTau{stype}(1),nTau{stype}(end));
        [~, MatCStrace{f,stype}]=get_STA(Result.CStrace,SpClass(stype,:),-nTau{stype}(1),nTau{stype}(end));
        [~, MatSub{f,stype}]=get_STA(Subthreshold{f},SpClass(stype,:),-nTau{stype}(1),nTau{stype}(end));
        [~, MatOrder{f,stype}]=get_STA(BrstOrder,SpClass(stype,:),-nTau{stype}(1),nTau{stype}(end));

        switch stype
            case 1
                [~, MatSTA{f,stype,1} sptime]=get_STA(NormalizedTrace_dirt{f,1},SpClass(stype,:),-nTau{stype}(1),nTau{stype}(end));
                SpikeList{f,stype,1}=ind2vec(size(NormalizedTrace_dirt{f},2),sptime,1);
            case 2
                for ns=1:5
                    [~, MatSTA{f,stype,ns} sptime]=get_STA(NormalizedTrace_dirt{f,1},(ComplexSpikeOrder)==ns,-nTau{stype}(1),nTau{stype}(end));
                    SpikeList{f,stype,ns}=ind2vec(size(NormalizedTrace_dirt{f},2),sptime,1);
                end
            case 3
                for ns=1:5
                    [~, MatSTA{f,stype,ns} sptime]=get_STA(NormalizedTrace_dirt{f,1},(BrstOrder.*(1-Result.CStrace))==ns,-nTau{stype}(1),nTau{stype}(end));
                    SpikeList{f,stype,ns}=ind2vec(size(NormalizedTrace_dirt{f},2),sptime,1);
                end
        end
    end

    allSpikeMat{f}=double(Result.spike);
    allSpikeMat{f}(:,Result.motionReject>0)=NaN;
    if ifdirtRemov(f)
        allSpikeMat{f}(Result.dirtTrace>0)=NaN;
    end
    allSpikeClassMat{f}=SpClass;
    allSpikeClassVecMat{f}=SpikeClassvec;

    BlueStim{f}=Result.Blue;
    VRtrack{f}=Result.VR;
    CStrace{f}=Result.CStrace;
    Ftprnts{f}=Result.ftprnt;
    AvgImg{f}=Result.ref_im;

    if ~isempty(PlaceFieldList{f}) % in place field
        binTrack=(ceil(VRtrack{f}(5,:)/((115)/150)));
        PFvec{f}=zeros(1,nTime);
        for p=1:length(PlaceFieldBin{f})/2
            if PlaceFieldBin{f}(2*(p-1)+1)>PlaceFieldBin{f}(2*(p-1)+2)
                Pvec=~(binTrack<(PlaceFieldBin{f}(2*(p-1)+1)) & binTrack>(PlaceFieldBin{f}(2*(p-1)+2)));
            else
                Pvec=(binTrack>(PlaceFieldBin{f}(2*(p-1)+1)) & binTrack<(PlaceFieldBin{f}(2*(p-1)+2)));
            end
            Lapvec=(VRtrack{f}(8,:)>PlaceFieldList{f}(2*(p-1)+1) & VRtrack{f}(8,:)<PlaceFieldList{f}(2*(p-1)+2));
            PFvec{f}=PFvec{f}| (Lapvec & Pvec);
        end
    end

    runVec{f}= double(imdilate(movmean((VRtrack{f}(end,:)),200)>0.002,[ones(1,2001)]));


    %LapFR{f}=PlaceTrigger_average(double(allSpikeMat{f}(1,:)==1),150,VRtrack{f},0.002,115); %total trace
    % LapV{f}=PlaceTrigger_average(NormalizedTrace_dirt{f},150,VRtrack{f},0.002,115); %total trace
    %LapSub{f}=PlaceTrigger_average(Subthreshold{f},150,VRtrack{f},0.002,115); %total trace
    % LapSpclass{f}=PlaceTrigger_average(double(allSpikeClassMat{f}>0),150,VRtrack{f},-0.002,115); %total trace
    % LapSpclassVec{f}=PlaceTrigger_average(double(allSpikeClassVecMat{f}>0),150,VRtrack{f},-0.002,115); %total trace
    %
    % subthreshold_silent=Subthreshold{f};
    % subthreshold_silent(:,unique(get_perispikeIndex(allSpikeMat{f}(1,:),[-8:30])))=NaN;
    % LapSubSilent{f}=PlaceTrigger_average(subthreshold_silent,150,VRtrack{f},0.002,115); %total trace
end

MatSTA=cellfun(@(x) permute(x,[1 3 2]),MatSTA,'UniformOutput',false);
MatSpClass=cellfun(@(x) permute(x,[1 3 2]),MatSpClass,'UniformOutput',false);
MatSpike=cellfun(@(x) permute(x,[1 3 2]),MatSpike,'UniformOutput',false);
MatBlue=cellfun(@(x) permute(x,[1 3 2]),MatBlue,'UniformOutput',false);
MatSub=cellfun(@(x) permute(x,[1 3 2]),MatSub,'UniformOutput',false);
MatOrder=cellfun(@(x) permute(x,[1 3 2]),MatOrder,'UniformOutput',false);
MatCStrace=cellfun(@(x) permute(x,[1 3 2]),MatCStrace,'UniformOutput',false);

%%
PeakSTA=[]; peakvec=[];
nTauPeak=[150 150];
for f=foi
StopFreq=[12 200];
[nROI nTime]=size(Subthreshold{f});
Subthreshold_int=interpolateNaN(Subthreshold{f});
FilterTr=[];
peakvec{f}=zeros(nROI,nTime);
figure(f+400); clf;
for n=1:nROI
    %[PhaseTrace BasalSubFilt BasalthetaPower] = get_phase(Subthreshold_int(n,:), 1000, FilterFreq);
perispikefrm=unique(find(allSpikeMat{f}(1,:))'+[-5:30]);
perispikefrm(perispikefrm<1 | perispikefrm>nTime)=[];
    FilterTr(n,:)=get_bandstop(Subthreshold_int(n,:),1000,StopFreq);
FilterTr(n,:)=FilterTr(n,:)-movmedian(FilterTr(n,:),300,2);
[pks, locs] = findpeaks(FilterTr(n,:),'MinPeakHeight', 0.5,'MinPeakDistance',50, ...
                                  'MinPeakProminence', 0.2);  % Prominence can also be tuned
peakvec{f}(n,:)=ind2vec(nTime,locs,1);
peakvec{f}(n,BlueStim{f}>0)=0;
peakvec{f}(n,perispikefrm)=0;
[PeakSTA{f}(:,:,n), peakMat, peakloc]=get_STA(Subthreshold{f},peakvec{f}(n,:),nTauPeak(1),nTauPeak(2));
peakvec{f}(n,:)=ind2vec(nTime,peakloc(sum(isnan(peakMat),[1 3])==0),1);
nexttile([1 1]);
imagesc(PeakSTA{f}(:,:,n)); 
colormap(turbo);
end
end
%%
time_segment=15000; bound=5; overlap=200; nTauPeak=[150 150];
for f=[foi(13)]%length(fpath)
    f
    [nROI nTime]=size(Subthreshold{f});
    load(fullfile(fpath{f},'PC_Result.mat'),'Result');
    load([fpath{f} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4])); blueDMDcontour=[];
    frm_end=EndFrame(f); BlueOffFrame=Result.Blue==0;
    
    f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
    f_seg_real=[f_seg(1:end-1)' f_seg(2:end)'];
    f_seg_real(1:end-1,2)=f_seg_real(1:end-1,2)+overlap;
    f_seg_real(2:end,1)=f_seg_real(2:end,1)+overlap+1;

    take_window=repmat([1 time_segment],length(f_seg)-1,1);
    take_window(2:end,1)=take_window(2:end,1)+overlap; take_window(1:end-1,2)=take_window(1:end-1,2)+overlap;
    take_window(end)=mod(f_seg(end),time_segment);
    take_window(take_window==0)=time_segment;

    perispike_time=unique(find(allSpikeMat{f}(1,:))'+[-5:30]); perispike_time(perispike_time<1 | perispike_time>nTime)=[];
    periblue_time=unique(find(BlueStim{f}>0)'+[-5:30]); periblue_time(periblue_time<1 | periblue_time>nTime)=[];
    t_fit= (ind2vec(size(Result.traces_bvMask,2),periblue_time,1)==0) & (ind2vec(size(Result.traces_bvMask,2),perispike_time,1)==0);
    [bleaching_fit] = expfitDM_2(find(t_fit)',-mean(Result.traces_bvMask(:,t_fit))',[1:size(Result.traces_bvMask,2)]',[100000 10000]);

    SpikeHeight=Result.SpikeHeight_fit;
    Mov_PeakTA=zeros(sz(2)*sz(1),sum(nTauPeak)+1,nROI);
    N_added=zeros(1,nROI);
    F0img=[];
    ftp=Ftprnts{f}(:,:,Dist_order{f}(noi_dist{f}));

    for j=1:length(f_seg)-1
        
        [nInd fInd]=find(peakvec{f}(:,[f_seg_real(j):f_seg_real(j+1)]));

        if ~isempty(nInd)
        mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
        load([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat']);

        mov_mc=mov_mc(:,:,[take_window(j,1):take_window(j,2)]);
        %mc= movmean(mcTrace.xymean-movmedian(mcTrace.xymean,500,1),3,1);
        mc=mcTrace.xymean([take_window(j,1):take_window(j,2)],:);

        bkg = zeros(1, size(mov_mc,3));
        bkg(1,:) = bleaching_fit(f_seg_real(j,1):f_seg_real(j,2));  % linear term

        mov_res= mov_mc-median(mov_mc,3);
        mov_res = SeeResiduals(mov_res,mc);
        mov_res = SeeResiduals(mov_res,mc.^2);
        mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
        mov_res = SeeResiduals(mov_res,bkg,1);
        mov_res=tovec(mov_res);
        mov_res= mov_res./SpikeHeight(f_seg_real(j,1):f_seg_real(j,2));
        if isempty(F0img)
        F0img=get_F0img(toimg(mov_res,sz(2),sz(1)));
        end
        mov_res=mov_res./tovec(F0img);
        
        for k=1:length(fInd)
            f2add=fInd(k);
            n2add=nInd(k);
            if f2add-nTauPeak(1)>1 & f2add+nTauPeak(2)<size(mov_mc,3)
        Mov_PeakTA(:,:,n2add)=Mov_PeakTA(:,:,n2add)+mov_res(:,f2add+[-nTauPeak(1):nTauPeak(2)]);
        N_added(1,n2add)=N_added(1,n2add)+1;
            end
        end
        disp([num2str(j) ' th segment is stacked']);
        else
            disp([num2str(j) ' has no valid index']);
        end
    end
    Mov_PeakTA=-Mov_PeakTA./reshape(N_added,1,1,[]);
% 
% for n=find(N_added>0)
%     PeakMov=toimg(Mov_PeakTA(:,:,n),sz(2),sz(1));
%     PeakMov=imgaussfilt(PeakMov,2);
%     cax_sub=[prctile(PeakMov(:),5) prctile(PeakMov(:),97)];
%     colorSTA=grs2rgb(tovec(PeakMov),colormap(turbo),cax_sub(1),cax_sub(2));
%     colorSTA=reshape(colorSTA,sz(2),sz(1),sum(nTauPeak)+1,[]);
%     colorSTA=permute(colorSTA,[1 2 4 3]);
%     if isfield(Result,'Structure')
%         PeakMov_sub_Struc=colorSTA.*mat2gray(Result.Structure)*3;
%     else
%         PeakMov_sub_Struc=colorSTA.*mat2gray(Result.ref_im-100)*2;
%     end
% 
%     ftp_boundary=bwboundaries(ftp(:,:,n));
%     figure(161); clf;
%     writeMov4d(fullfile(fpath{f},['PeakTriggered_movie_' num2str(n,2)]), ...
%         PeakMov_sub_Struc,[-nTauPeak(1):nTauPeak(2)],10,1,cax_sub,ftp_boundary{1}(:,[2 1]))
% end
    Mov_PeakTA=vm(mat2gray(double(Mov_PeakTA))*10000);
    Mov_PeakTA.transpose.savebin([fpath{f} '/PeakTriggered_movie.bin'])
end

%%




[BasalPhase BasalSubFilt BasalthetaPower] = get_phase(BasalSub_intNan, 1000, FilterFreq);
[ApicalPhase ApicalSubFilt ApicalthetaPower] = get_phase(ApicalSub_intNan, 1000, FilterFreq);
Run_thetapower(f,:,1)=[mean(BasalthetaPower(runtime>0)) mean(BasalthetaPower(runtime==0))];
Run_thetapower(f,:,2)=[mean(ApicalthetaPower(runtime>0)) mean(ApicalthetaPower(runtime==0))];

    crossingVecApical = [((ApicalPhase(1:end-1) < phaseTarget) & (ApicalPhase(2:end) > phaseTarget))' 0];
    crossingVecBasal = [((BasalPhase(1:end-1) < phaseTarget) & (BasalPhase(2:end) > phaseTarget))' 0];
    ApicalphaseVec=(crossingVecApical & zscore(ApicalthetaPower)>1.5);
    ApicalphaseVec_silent=(ApicalphaseVec & silenttime_vec);
    BasalphaseVec=(crossingVecBasal & zscore(BasalthetaPower)>1.5);
    BasalphaseVec_silent=(BasalphaseVec & silenttime_vec);