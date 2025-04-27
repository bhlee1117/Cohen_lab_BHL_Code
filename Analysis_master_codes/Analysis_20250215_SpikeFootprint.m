clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:Z31');

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
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
PlaceFieldList=cellfun(@(x) (str2num(num2str(x))),raw(:,21),'UniformOutput',false);
PlaceFieldBin=cellfun(@(x) (str2num(num2str(x))),raw(:,22),'UniformOutput',false);
set(0,'DefaultFigureWindowStyle','docked')
foi=[1 4 5 6 8 10 11 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27];
%foi=[1 4 5 6 8 10 11 15 16 17 18 19 20 21 22 23 24 25 26 27];
%foi=18;
%%
nTau={[-200:200],[-200:200],[-200:200]}; %SS, CS, Brst
clear SpikeInd MatSpike MatSTA MatBlue MatCStrace MatSub SpikeList NormalizedTrace_ch NormalizedTrace_dirt SpikeIndBlueOff Dist_order allSpikeMat noi interDendDist noi_dist derivSub LapSubSilent
clear Subthreshold dendaxis BrstOrder roisD roisD_order

for f=foi
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
    NormalizedTrace_ch(f,:)=cellfun(@(x) x./Result.F_ref,Result.norm_trace_check,'UniformOutput',false);
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

    % spike class, SS:1, CS:2, dSP:3, BS:4
    ss_time=find(Result.SpClass(1,:)); % BS is subclass of SS
    brst=bwlabel((ss_time(2:end)-ss_time(1:end-1))<=20); % SSs that have an ISI shorter than 20 ms are BS.
    SpClass=Result.SpClass; BS_trace=zeros(1,size(Result.traces,2));
    for b=1:max(brst)
        bwn=find(brst==b);
        SpClass(1,ss_time([bwn bwn(end)+1]))=0;
        SpClass(4,ss_time([bwn(1)]))=1;
        BS_trace(1,[ss_time(bwn): ss_time(bwn(end)+1)])=b;
    end
    SpClass=SpClass([1 2 4],:);
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

    BlueStim{f}=Result.Blue;
    VRtrack{f}=Result.VR;
    CStrace{f}=Result.CStrace;
    Ftprnts{f}=Result.ftprnt;
    AvgImg{f}=Result.ref_im;

    % LapFR{f}=PlaceTrigger_average(double(allSpikeMat{f}(1,:)==1),150,VRtrack{f},0.002,115); %total trace
    % LapV{f}=PlaceTrigger_average(NormalizedTrace_dirt{f},150,VRtrack{f},0.002,115); %total trace
    % LapSub{f}=PlaceTrigger_average(Subthreshold{f},150,VRtrack{f},0.002,115); %total trace
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
f=26;
%CStr=reshape(MatSTA{f,2,1}(:,-nTau{stype}(1)-5:end-105,:),size(MatSTA{f,2,1},1),[]);
%tr=movmean(CStr,3,2,'omitnan');
tr=movmean(NormalizedTrace_dirt{f},5,2,'omitnan');
eigTrace=NaN(size(tr,1),size(tr,2));
validFrame=find(sum(isnan(tr),1)==0);
tr=tr(:,validFrame);
[V D T]=get_eigvector(tr,size(tr,1));
eigTrace(:,validFrame)=T';

figure(10); clf; tiledlayout(4,1);
ax1=nexttile([1 1]);
imagesc(NormalizedTrace_dirt{f},[-0.5 1.5]); colormap(turbo);
title('dF/F')
ax2=nexttile([1 1]);
filteredTr=[];
filteredTr(:,validFrame)=pcafilterTrace(tr,[1 2 3]);
imagesc(filteredTr,[-0.5 1.5]); colormap(turbo);
title('PC#1-3')
ax3=nexttile([1 1]);
ResidualTr=[];
ResidualTr(:,validFrame)=pcafilterTrace(tr,[4:20]);
imagesc(ResidualTr,[-0.5 1.5]); colormap(turbo);
title('Residual')
ax4=nexttile([1 1]);
%plot(mean(NormalizedTrace_dirt{f}(find(roisD_order{f}{2, 1}>0),:),1,'omitnan'));
%plot(mean(CStr(find(roisD_order{f}{2, 1}>0),validFrame),1,'omitnan'));
plot(mean(ResidualTr.^2,1,'omitnan'));
title('Mean square')
linkaxes([ax1 ax4 ax2 ax3],'x')
%%
read_frm=354923+[-2000:3000];
mov=readBinMov_BHL_multiple(fpath{f},3,read_frm,15000);
mov_res=mov-median(mov,3);
mov_res=SeeResiduals(mov_res,Result.mcTrace(read_frm,:));
mov_res=SeeResiduals(mov_res,Result.mcTrace(read_frm,:).^2);
mov_res=SeeResiduals(mov_res,Result.mcTrace(read_frm,1).*Result.mcTrace(read_frm,end));
%%