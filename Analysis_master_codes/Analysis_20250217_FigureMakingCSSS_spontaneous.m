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
%foi=[1 4 5 6 8 10 11 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27];
foi=[1 4 5 6 8 10 11 15 16 17 18 19 20 21 22 23 24 25 26 27];
%foi=23;
%%
nTau={[-200:50],[-200:50],[-200:50]}; %SS, CS, Brst
clear SpikeInd MatSpike MatSTA MatBlue MatCStrace MatSub SpikeList NormalizedTrace_ch NormalizedTrace_dirt SpikeIndBlueOff Dist_order allSpikeMat noi interDendDist noi_dist derivSub LapSubSilent
clear Subthreshold dendaxis BrstOrder roisD roisD_order LapSpclassVec

for f=23
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

    %NormalizedTrace=NormalizedTrace/SpikeHeight;
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
    SpikeClassvec=Classvec.*max([Result.spike(1,:); SpClass(4,:)]);

    BrstOrder{f} = get_BurstOrder(Result.spike(1,:), 20) - SpClass(1,:);
    BrstOrder{f}(find(SpClass(3,:)))=1;
    ComplexSpikeOrder=get_spikeOrder(Result.CStrace,Result.spike(1,:));
    

    for stype=1:3
        s_list=find(SpClass(stype,:)>0);
        [~, MatSpClass{f,stype}]=get_STA(SpikeClassvec,SpClass(stype,:),-nTau{stype}(1),nTau{stype}(end));
        [~, MatSpike{f,stype}]=get_STA(Result.spike(1,:),SpClass(stype,:),-nTau{stype}(1),nTau{stype}(end));
        [~, MatBlue{f,stype}]=get_STA(Result.Blue,SpClass(stype,:),-nTau{stype}(1),nTau{stype}(end));
        [~, MatCStrace{f,stype}]=get_STA(Result.CStrace,SpClass(stype,:),-nTau{stype}(1),nTau{stype}(end));
        [~, MatSub{f,stype}]=get_STA(Subthreshold{f},SpClass(stype,:),-nTau{stype}(1),nTau{stype}(end));
        [~, MatOrder{f,stype}]=get_STA(BrstOrder{f},SpClass(stype,:),-nTau{stype}(1),nTau{stype}(end));

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
                    [~, MatSTA{f,stype,ns} sptime]=get_STA(NormalizedTrace_dirt{f,1},(BrstOrder{f}.*(1-Result.CStrace))==ns,-nTau{stype}(1),nTau{stype}(end));
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

    % LapFR{f}=PlaceTrigger_average(double(allSpikeMat{f}(1,:)==1),150,VRtrack{f},0.002,115); %total trace
    % LapV{f}=PlaceTrigger_average(NormalizedTrace_dirt{f},150,VRtrack{f},0.002,115); %total trace
    % LapSub{f}=PlaceTrigger_average(Subthreshold{f},150,VRtrack{f},-0.002,115); %total trace
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

%% Figure 3
f=20; 
load(fullfile(fpath{f},'PC_Result.mat'),'Result');
% CS_STA=read_tiff(fullfile(fpath{f},'CS_1.tiff'));
% STAmovieCS=-mean(reshape(CS_STA,size(CS_STA,1),size(CS_STA,2),101,[]),4);
% figure(4); clf;
% imshow2(imgaussfilt(max(STAmovieCS,[],3),0.5),[0 230]);

Result.spike=Result.spike>0; Result.SpClass=Result.SpClass>0;
Dist_order=Result.BrancDist_order;
interDendDist=Result.interDendDist*Pixelsize(f);
Dsign=ones(1,size(interDendDist,1));
Dsign(Dist_order(1:find(Dist_order==1)-1))=-1;
perisomaROI=setdiff(find(interDendDist(1,:)<60),BadROI{f}); % ROIs < 40 um from soma
noi=setdiff([1:size(Result.ftprnt,3)],BadROI{f});
%noi=setdiff([1:size(Result.ftprnt,3)],[]);
noi_dist=ismember(Dist_order,noi);

roisD(f,:)={basal_ROI{f},PeriSoma_ROI{f},apical_ROI{f},oblique_ROI{f},distal_ROI{f}}; %set the ROIs
for dClass=1:size(roisD,2)
    g=1;
    if ~isnan(roisD{f,dClass})
        for d=roisD{f,dClass}
            dind=setdiff(find(Result.BranchLabel==d),BadROI{f});
            roisD_order{f}{dClass,g}=ismember(Dist_order(noi_dist),dind);
            g=g+1;
        end
    end
end

roisD_order_ind=cellfun(@find,roisD_order{f},'UniformOutput',false);

NormalizedTrace=(Result.normTraces)./Result.F_ref;
bAP_STA=get_STA(NormalizedTrace, Result.spike(1,:).*double(Result.Blue==0), 30, 20);
bAP_STA=bAP_STA-prctile(bAP_STA,20,2);
SpikeHeight=max(mean(bAP_STA(perisomaROI,:),1),[],2);
NormalizedTrace=NormalizedTrace(Dist_order(noi_dist),:)/SpikeHeight;

ss_time=find(Result.SpClass(1,:)); % BS is subclass of SS
brst=bwlabel((ss_time(2:end)-ss_time(1:end-1))<=20); % SSs that have an ISI shorter than 20 ms are BS.
SpClass=Result.SpClass; BS_trace=zeros(1,size(Result.traces,2));
for b=1:max(brst)
    bwn=find(brst==b);
    SpClass(1,ss_time([bwn bwn(end)+1]))=0;
    SpClass(4,ss_time([bwn(1)]))=1;
    BS_trace(1,[ss_time(bwn(1)): ss_time(bwn(end)+1)])=b;
end
SpClass=SpClass([1 2 4],:);
Classvec = get_Class2index(SpClass);
SpikeClassvec=Classvec.*Result.spike(1,:);

figure(301); clf; tiledlayout(1,3); ax1=[]; cmap=distinguishable_colors(3);
tshow=[87454 22019 19084]; nTau=[-50:80]; show_noi=setdiff([1:size(NormalizedTrace,1)],[2 39 40 38 15 17]);
tr_show=mean(NormalizedTrace([25 10 9],:),1,'omitnan');
%tr_show=pcafilterTrace(NormalizedTrace,5);
sub_show=get_subthreshold(mean(NormalizedTrace([25 10 9],:),1,'omitnan'),Result.spike(1,:),7,15);
for t=1:length(tshow)
    ax1=[ax1 nexttile([1 1])];
    tr_show_tmp=tr_show(tshow(t)+nTau);
    %tr_show_tmp=tr_show_tmp-median(tr_show_tmp(:,1:50),2);
    l=plot(nTau,tr_show_tmp','k','linewidth',1.5);
    %arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(gen_colormap([0 0 0; 0.5 0 0; 1 0 0],length(show_noi)),2));
    hold all
    plot(nTau,sub_show(tshow(t)+nTau),'color',[1 0.5 0.5],'linewidth',2)
    hold all
    plot(0,1.6,'.','color',cmap(t,:),'markersize',20); hold all 
    axis off
end
linkaxes(ax1,'xy');


%% Figure 3, show STA kymo
f=23; 
figure(101); clf; tiledlayout(2,6)
stype_str={'SS','CS','BS','dSP'}; ax1=[]; ax2=[]; FR_bin_Frame=12; ax3=[];
nTau=[200 70]; cmap=hsv(4);
nTime=size(NormalizedTrace_dirt{f},2);
%NormalizedTraceFiltered=pcafilterTrace(NormalizedTrace_dirt{f},[1:15]);
NormalizedTraceFiltered=NormalizedTrace_dirt{f};
Subthreshold=get_subthreshold(NormalizedTrace_dirt{f},allSpikeMat{f}(1,:),7,17); cax=[-0.001 0.01];

roisD_order_ind=cellfun(@find,roisD_order{f},'UniformOutput',false);
if isempty(cell2mat(roisD_order_ind(1,:)'))
    basalind=cell2mat(roisD_order_ind(2,:)'); %if there is no basal, use soma
else
    basalind=cell2mat(roisD_order_ind(1,:)');
end
apicalind=cell2mat(roisD_order_ind(5,:)'); %apical = distal dend
somaind=cell2mat(roisD_order_ind(2,:)');

Basal_Trace= mean(NormalizedTraceFiltered(basalind,:),1,'omitnan');
Apical_Trace= mean(NormalizedTraceFiltered(apicalind,:),1,'omitnan');

silenttime=setdiff([1:nTime],unique(find(max([allSpikeMat{f}(1,:); allSpikeClassMat{f}; CStrace{f}; BlueStim{f}])>0)'+[-5:30]));
silenttime_vec=ind2vec(nTime,silenttime,1);

BASub =[Basal_Trace; Apical_Trace];
BASub_silent=[Basal_Trace; Apical_Trace];
BASub_silent=get_subthreshold(BASub_silent,allSpikeMat{f}(1,:),7,17);
BASub_silent(:,silenttime_vec==0)=NaN;

STApreSpikeBA=[];
for sclass=1:3
    TriggerSpike=allSpikeClassMat{f}(sclass,:).*double(BlueStim{f}==0);
    [STATrace STATraceMat Triggerspike_used]=get_STA(NormalizedTrace_dirt{f},TriggerSpike,nTau(1),nTau(2));
    [~, NMat]=get_STA(max([allSpikeMat{f}(1,:); CStrace{f}]),ind2vec(nTime,Triggerspike_used,1),nTau(1),-1);
    [~, SpMat]=get_STA(allSpikeClassVecMat{f},ind2vec(nTime,Triggerspike_used,1),nTau(1),nTau(2));
    Isolate_ind=find(sum(squeeze(sum(NMat,1,'omitnan')),2,'omitnan')==0);
    %Isolate_ind=[1:size(STATraceMat,2)];
    %length(Isolate_ind)
    STATrace=mean(STATraceMat(:,Isolate_ind,:),2,'omitnan');
    [~, STASpikeClassMat]=get_STA(allSpikeClassVecMat{f},TriggerSpike,nTau(1),nTau(2));
    STASpikeClass=mean(STASpikeClassMat(:,Isolate_ind,:),2,'omitnan');
    %STASpikeClass(:,nTau(1):end)=0;
    STAsub=get_STA(Subthreshold,TriggerSpike,nTau(1),nTau(2));
    [~, STApreSpikeBA{sclass}] = get_STA(BASub,TriggerSpike,nTau(1),nTau(2));

    [STASpikeClass_bin binT]= Bin_Vector(STASpikeClass, [-nTau(1):nTau(2)], [-nTau(1):FR_bin_Frame:-1 0:1 2:FR_bin_Frame:nTau(2)]);
    %STASpikeClass_bin=STASpikeClass;
    [~, dind]=sort(dendaxis{f},'ascend');
    [X, Y] = meshgrid([-nTau(1):nTau(2)], min(dendaxis{f}):15:max(dendaxis{f}));
    STAtraceq = interp2([-nTau(1):nTau(2)], dendaxis{f}(dind), STATrace(dind,:), X, Y, 'linear');
    STAsubq = interp2([-nTau(1):nTau(2)], dendaxis{f}(dind), STAsub(dind,:), X, Y, 'linear');

    ax1=[ax1 nexttile(sclass,[1 1])];
   %imagesc([-nTau(1):nTau(2)],Y(:,1),STAtraceq,cax) 
   imagesc([-nTau(1):nTau(2)],dendaxis{f}(dind),STATrace(dind,:),cax) 
   ylabel('Distance from soma (\mum)')
colormap(turbo(256))

    ax3=[ax3 nexttile(sclass+4,[1 1])];
    STAtraceMatB=mean(STATraceMat(basalind,Isolate_ind,:),1,'omitnan'); 
    STAtraceMatA=mean(STATraceMat(apicalind,Isolate_ind,:),1,'omitnan');
    STAtraceMatS=mean(STATraceMat(somaind,Isolate_ind,:),1,'omitnan');
    
    BasalM=squeeze(mean(STAtraceMatB,[2],'omitnan')); BasalS=squeeze(std(STAtraceMatB,0,2,'omitnan'));
    ApicalM=squeeze(mean(STAtraceMatA,2,'omitnan')); ApicalS=squeeze(std(STAtraceMatA,0,2,'omitnan'));
    SomaM=squeeze(mean(STAtraceMatS,[2],'omitnan')); SomaS=squeeze(std(STAtraceMatS,0,2,'omitnan'));
    h=[];
    h(1)=plot([-nTau(1):nTau(2)],BasalM,'color',[0.1 0.5 1]); hold all
    h(2)=plot([-nTau(1):nTau(2)],ApicalM,'color',[1 0 0]); hold all
    h(3)=plot([-nTau(1):nTau(2)],SomaM,'color',[1 0.7 0]); hold all
    title(['N = ' num2str(length(Isolate_ind))])

ax2=[ax2 nexttile(sclass+8,[1 1])];
% l=plot(binT,STASpikeClass_bin');
% arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap,2));
%imagesc(permute(SpMat,[2 3 1]));
for spclass=1:3
[r c]=find(permute(SpMat(spclass,:,:)>0,[2 3 1]));
scatter(c-nTau(1),r,15,cmap(spclass,:),'filled'); hold all
xlabel('Peri-spike time (ms)')
set(gca,'YDir','reverse')
axis tight
end

if sclass==4
mean(sum(SpMat(:,:,nTau(1):nTau(1)+3),3,'omitnan'),2)
end

ylabel([stype_str{sclass} ' #'])
end
legend(stype_str)
ax1=[ax1 nexttile(3,[1 1])];
%colorbar;
nexttile(6,[1 1]);
legend(h,{'Basal','Distal','Soma'})
linkaxes(ax1,'xy');
linkaxes(ax2,'x');
linkaxes(ax3,'xy');
linkaxes([ax1 ax2 ax3],'x');
xlim([-nTau(1) nTau(2)])
%% SS, CS, BS, dSP examples
f=23;
STATraceMat_show=[];
STATrace_show=[];
for sclass=1:3
    TriggerSpike=allSpikeClassMat{f}(sclass,:).*double(BlueStim{f}==0);
    [STATrace_show{sclass} STATraceMat_show{sclass} Triggerspike_used]=get_STA(NormalizedTrace_dirt{f},TriggerSpike,30,30);
end

SS_show=[2 3 4 9 10 12 33 39 40 42];
CS_show=[1 6 7 17 9 16 20 29 36 54];
figure(30); clf; tiledlayout(2,5);
for sp=SS_show
    nexttile([1 1]);
imagesc([-30:30],[1:size(NormalizedTrace_dirt{f},1)],permute(STATraceMat_show{1}(:,sp,:),[1 3 2]),cax)
Dsign=ones(1,size(interDendDist{f},1));
Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
if find(Dist_order{f}(noi_dist{f})==1)~=1
    set(gca,'YTick',[1 find(Dist_order{f}(noi_dist{f})==1) sum(noi_dist{f})-1],'YTickLabel',num2str([min(interDendDist{f}(1,:).*Dsign) 0 max((interDendDist{f}(1,:).*Dsign))]',3))
else
    set(gca,'YTick',[1 sum(noi_dist{f})-1],'YTickLabel',num2str([0 max((interDendDist{f}(1,:)))]',3))
end
end
colormap(turbo)
figure(31); clf; tiledlayout(2,5)
for sp=CS_show
    nexttile([1 1]);
imagesc([-30:30],[1:size(NormalizedTrace_dirt{f},1)],permute(STATraceMat_show{2}(:,sp,:),[1 3 2]),cax)
Dsign=ones(1,size(interDendDist{f},1));
Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
if find(Dist_order{f}(noi_dist{f})==1)~=1
    set(gca,'YTick',[1 find(Dist_order{f}(noi_dist{f})==1) sum(noi_dist{f})-1],'YTickLabel',num2str([min(interDendDist{f}(1,:).*Dsign) 0 max((interDendDist{f}(1,:).*Dsign))]',3))
else
    set(gca,'YTick',[1 sum(noi_dist{f})-1],'YTickLabel',num2str([0 max((interDendDist{f}(1,:)))]',3))
end
end

colormap(turbo)

f=20;
STATraceMat_show=[];
STATrace_show=[];
for sclass=4
    TriggerSpike=allSpikeClassMat{f}(sclass,:).*double(BlueStim{f}==0);
    [STATrace_show{sclass} STATraceMat_show{sclass} Triggerspike_used]=get_STA(NormalizedTrace_dirt{f},TriggerSpike,50,50);
end
%% plateau example
PlateauShowIndex=[4 86400 87200 16 -0.004 0.013;18 25470 25870 10 -0.004 0.009]; % Neuron & CS index
figure(32); clf; tiledlayout(2,1);
for p=[1 2]
nexttile([1 1]);
show_tr=NormalizedTrace_dirt{PlateauShowIndex(p,1)}(:,PlateauShowIndex(p,2):PlateauShowIndex(p,3));
imagesc(show_tr(setdiff([1:size(show_tr,1)],PlateauShowIndex(p,4)),:),[PlateauShowIndex(p,5) PlateauShowIndex(p,6)])
f=PlateauShowIndex(p,1);
Dsign=ones(1,size(interDendDist{f},1));
Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
if find(Dist_order{f}(noi_dist{f})==1)~=1
    set(gca,'YTick',[1 find(Dist_order{f}(noi_dist{f})==1) sum(noi_dist{f})-1],'YTickLabel',num2str([min(interDendDist{f}(1,:).*Dsign) 0 max((interDendDist{f}(1,:).*Dsign))]',3))
else
    set(gca,'YTick',[1 sum(noi_dist{f})-1],'YTickLabel',num2str([0 max((interDendDist{f}(1,:)))]',3))
end
colorbar;
ylabel(['Distance from' newline 'soma (\mum)'])
end
colormap(turbo)
%% Blue on vs off

nTau2=[30 50];
STAmat_blueOnoff=[]; SPmat_blueOnoff=[];
SlopeOff=[]; SlopeOn=[];
figure(15); clf; tiledlayout(5,3); cax=[-0.002 0.015]*50;
f=26;

roisD_order_ind=cellfun(@find,roisD_order{f},'UniformOutput',false);
if isempty(cell2mat(roisD_order_ind(1,:)'))
    basalind=cell2mat(roisD_order_ind(2,:)'); %if there is no basal, use soma
else
    basalind=cell2mat(roisD_order_ind(1,:)');
end
apicalind=cell2mat(roisD_order_ind(5,:)'); %apical = distal dend

for spclass=1:3
[~, STAmat_ch]=get_STA(NormalizedTrace_dirt{f},allSpikeClassMat{f}(spclass,:),nTau2(1),nTau2(2));
[~, SPmat]=get_STA(max([allSpikeMat{f}(1,:); CStrace{f}]),allSpikeClassMat{f}(spclass,:),nTau2(1),nTau2(2));
[~, BlueMat]=get_STA(BlueStim{f}>0,allSpikeClassMat{f}(spclass,:),nTau2(1),nTau2(2));

Blueoffind=squeeze(sum(BlueMat,3))==0;
Blueonind=squeeze(sum(BlueMat,3))==(sum(nTau2)+1);

Isolatedind=squeeze(sum(SPmat(1,:,1:nTau2(1)),3))==0;

STAmat_blueOnoff{f,spclass,1}=STAmat_ch(:,Blueoffind & Isolatedind,:);
STAmat_blueOnoff{f,spclass,2}=STAmat_ch(:,Blueonind & Isolatedind,:);

% STAmat_blueOnoff{f,spclass,1}=STAmat(:,Blueoffind,:);
% STAmat_blueOnoff{f,spclass,2}=STAmat(:,Blueonind,:);

SPmat_blueOnoff{f,spclass,1}=squeeze(SPmat(:,Blueoffind,:)); %off
SPmat_blueOnoff{f,spclass,2}=squeeze(SPmat(:,Blueonind,:)); %on

for s=1:size(STAmat_blueOnoff{f,spclass,1},2)
SlopeOff{spclass}(s)=mean(STAmat_blueOnoff{f,spclass,1}(apicalind,s,nTau2(1)-10:nTau2(1)-5),[1 3],'omitnan') - mean(STAmat_blueOnoff{f,spclass,1}(basalind,s,nTau2(1)-10:nTau2(1)-5),[1 3],'omitnan');
end
for s=1:size(STAmat_blueOnoff{f,spclass,2},2)
SlopeOn{spclass}(s)=mean(STAmat_blueOnoff{f,spclass,2}(apicalind,s,nTau2(1)-10:nTau2(1)-5),[1 3],'omitnan') - mean(STAmat_blueOnoff{f,spclass,2}(basalind,s,nTau2(1)-10:nTau2(1)-5),[1 3],'omitnan');
end

STA_off=squeeze(mean(STAmat_blueOnoff{f,spclass,1},2,'omitnan'));
STA_on=squeeze(mean(STAmat_blueOnoff{f,spclass,2},2,'omitnan'));
STA_diff=STA_on-STA_off;
nexttile(spclass,[1 1])
imagesc(STA_off,cax)
colorbar;
title(['N = ' num2str(size(STAmat_blueOnoff{f,spclass,1},2))])
nexttile(spclass+3,[1 1])
imagesc(STA_on,cax)
colorbar;
title(['N = ' num2str(size(STAmat_blueOnoff{f,spclass,2},2))])
nexttile(spclass+6,[1 1])
imagesc(STA_diff,cax/2)
colorbar;
nexttile(spclass+9,[1 1])
% hs=histogram(SlopeOff{spclass},[-1.5:0.1:1.5],'Normalization','probability'); hold all
% histogram(SlopeOn{spclass},[-1.5:0.1:1.5],'Normalization','probability'); hold all
scatter(mean(STAmat_blueOnoff{f,spclass,1}(basalind,:,nTau2(1)-10:nTau2(1)-2),[1 3],'omitnan'),mean(STAmat_blueOnoff{f,spclass,1}(apicalind,:,nTau2(1)-10:nTau2(1)-5),[1 3],'omitnan'),10,'red','filled'); hold all
scatter(mean(STAmat_blueOnoff{f,spclass,2}(basalind,:,nTau2(1)-10:nTau2(1)-2),[1 3],'omitnan'),mean(STAmat_blueOnoff{f,spclass,2}(apicalind,:,nTau2(1)-10:nTau2(1)-5),[1 3],'omitnan'),10,'cyan','filled')
xlabel('Basal');
ylabel('Apical')
box on

nexttile(spclass+12,[1 1])
hs=histogram(SlopeOff{spclass},[-1.5:0.1:1.5],'Normalization','probability'); hold all
histogram(SlopeOn{spclass},[-1.5:0.1:1.5],'Normalization','probability'); hold all
xlabel('Apical - Basal')
ylabel('Probability')
end
legend({'Blue Off','Blue On'})
nexttile(12,[1 1])
legend({'Blue Off','Blue On'})
colormap(turbo);

figure(16); clf; ax3=[];
ax3=[ax3 nexttile([1 1])];
p_values = Boxplot_wPoints2(SlopeOff, cmap(1:3,:));
set(gca,'XTick',[1 2 3],'XTickLabel',{'SS','CS','BS'})
ylabel('Apical - Basal')
title('Blue Off')
box on

ax3=[ax3 nexttile([1 1])];
p_values = Boxplot_wPoints2(SlopeOn, cmap(1:3,:));
set(gca,'XTick',[1 2 3],'XTickLabel',{'SS','CS','BS'})
ylabel('Apical - Basal')
title('Blue On')
box on

linkaxes(ax3)

%% Rest vs Run
RestTime=[];
for f=foi
bwrun=bwlabel(VRtrack{f}(end,:)<0.002);
runarea=regionprops('table',bwrun,'Area');
runarea=table2array(runarea);
RestTime{f}=ismember(bwrun,find(runarea>500));
end

RestFiringRate=cellfun(@(x,y,z) mean(x(:,(y>0 & z==0)),2,'omitnan')*1000,allSpikeClassMat,RestTime,BlueStim,'Uniformoutput',false);
RunFiringRate=cellfun(@(x,y,z) mean(x(:,(y==0 & z==0)),2,'omitnan')*1000,allSpikeClassMat,RestTime,BlueStim,'Uniformoutput',false);

RestFiringRateMat=cell2mat(RestFiringRate(foi));
RunFiringRateMat=cell2mat(RunFiringRate(foi));

figure(16); clf; tiledlayout(1,4); ax1=[];
for sp=1:4
    ax1=[ax1 nexttile([1 1])];
    Boxplot_wPoints(RestFiringRateMat(sp,:),RunFiringRateMat(sp,:),cmap(sp,:))
    set(gca,'XTick',[1 2],'XTickLabel',{'Rest','Run'})
    ylabel('Firing rate (Hz)')
    title(stype_str{sp})
end
linkaxes(ax1,'x')


%% Blue on vs off spike rate
foi_blue=[20 21 22 23 24 25 26 27];
RestTime=[];

BlueOnFiringRate=cellfun(@(x,y) mean(x(:,(y>0)),2,'omitnan')*1000,allSpikeClassMat,BlueStim,'Uniformoutput',false);
BlueOffFiringRate=cellfun(@(x,y) mean(x(:,(y==0)),2,'omitnan')*1000,allSpikeClassMat,BlueStim,'Uniformoutput',false);

BlueOnFiringRateMat=cell2mat(BlueOnFiringRate(foi_blue));
BlueOffFiringRateMat=cell2mat(BlueOffFiringRate(foi_blue));

figure(16); clf;
FRratioMat=BlueOnFiringRateMat./BlueOffFiringRateMat;
for sp=1:4
    FRratio{sp}=FRratioMat(sp,:)';
end

ax1=[ax1 nexttile([1 1])];
Boxplot_wPoints2(FRratio,cmap)
set(gca,'XTick',[1 2 3 4],'XTickLabel',stype_str)
ylabel('On/Off ratio')

%% Rank of bAP


nTau2=[1 3];
STAmat_blueOnoff=[]; SPmat_blueOnoff=[];
SlopeOff=[]; SlopeOn=[];
figure(15); clf; cax=[-0.002 0.015]*50;

for f=23;
    roisD_order_ind=cellfun(@find,roisD_order{f},'UniformOutput',false);
    if isempty(cell2mat(roisD_order_ind(1,:)'))
        basalind=cell2mat(roisD_order_ind(2,:)'); %if there is no basal, use soma
    else
        basalind=cell2mat(roisD_order_ind(1,:)');
    end
    apicalind=cell2mat(roisD_order_ind(5,:)'); %apical = distal dend
    somaind=cell2mat(roisD_order_ind(2,:)');

    [~, dsort]=sort(dendaxis{f});

    STAmat_ch=[]; STAmat_Amp_dend=[];
    %sptrace=allSpikeMat{f}(1,:);
    sptrace=max(allSpikeClassMat{f}(1:3,:),[],1)>0;

    [~, STAmat sptime]=get_STA(NormalizedTrace_dirt{f},sptrace,nTau2(1),nTau2(2));
    for ch=1:2
        NormTr_hi=NormalizedTrace_ch{f,ch};%-movmedian(NormalizedTrace_ch{f},3000,2);
        [~, STAmat_ch{ch} sptime]=get_STA(NormTr_hi,sptrace,nTau2(1),nTau2(2));
        %STAmat_Amp_dend{ch}=max(STAmat_ch{ch},[],3);
    end
    brst_trace=[];
    for b=1:max(BrstOrder{f})
        brst_trace(b,:)=BrstOrder{f}==b;
    end
    [~, SPmat]=get_STA(allSpikeClassVecMat{f}(1:3,:),sptrace,0,0);
    [~, BstMat]=get_STA(brst_trace,sptrace,0,0);
    [~, BlueMat]=get_STA(BlueStim{f}>0,sptrace,nTau2(1),nTau2(2));

    Blueoffsp=sum(BlueMat,3,'omitnan')==0;
    nonnansp=sum(isnan(STAmat),[1 3],'omitnan')==0;

    STAmat_Amp_ch=cellfun(@(x) max(x(:,Blueoffsp & nonnansp,:),[],3),STAmat_ch,'UniformOutput',false);
    sptime=sptime(Blueoffsp & nonnansp);
    STAmat_Amp=max(STAmat(:,Blueoffsp & nonnansp,:),[],3);
    
    STAmat_AmpNorm=STAmat_Amp./mean(STAmat_Amp(somaind,:),1,'omitnan');
    [V D eigSpike]=get_eigvector(STAmat_AmpNorm,10);
    [~, eigSort]=sort(eigSpike(:,1));
    
    if find(Dist_order{f}(noi_dist{f})==1)~=1
    set(gca,'YTick',[1 find(Dist_order{f}(noi_dist{f})==1) sum(noi_dist{f})-1],'YTickLabel',num2str([min(interDendDist{f}(1,:).*Dsign) 0 max((interDendDist{f}(1,:).*Dsign))]',3))
else
    set(gca,'YTick',[1 sum(noi_dist{f})-1],'YTickLabel',num2str([0 max((interDendDist{f}(1,:)))]',3))
end

    %set label index
    nTime=size(NormalizedTrace_dirt{f},2);
    roisD_order_ind=cellfun(@find,roisD_order{f},'UniformOutput',false);
    labelvec=NaN(1,size(NormalizedTrace_dirt{f},1));
    labelclass=[];
    STAmat_Amp_dend{1}=[]; STAmat_Amp_dend{2}=[];
    for dClass=1:5
        labelvec(cell2mat(roisD_order_ind(dClass,:)'))=dClass;
        for dend=1:size(roisD_order_ind,2)
            if ~isempty(roisD_order_ind{dClass,dend})
                labelclass=[labelclass dClass];
                STAmat_Amp_dend{1}=[STAmat_Amp_dend{1}; mean(STAmat_Amp_ch{1}(roisD_order_ind{dClass,dend},:),1,'omitnan')];
                STAmat_Amp_dend{2}=[STAmat_Amp_dend{2}; mean(STAmat_Amp_ch{2}(roisD_order_ind{dClass,dend},:),1,'omitnan')];
            end
        end
    end

    
    corrMat=get_corrMat(STAmat_Amp_ch{1}./mean(STAmat_Amp(somaind,:),1,'omitnan'),STAmat_Amp_ch{2}./mean(STAmat_Amp(somaind,:),1,'omitnan'));

    corrMat_dend=get_corrMat(STAmat_Amp_dend{1},STAmat_Amp_dend{2});
    Corrcoeff{f}=NaN(5,5);
   
    intdistMat=max(cat(3,interDendDist{f},interDendDist{f}'),[],3);
    intdistMat=intdistMat(Dist_order{f}(noi_dist{f}),Dist_order{f}(noi_dist{f}));
    
    figure(39); clf; tiledlayout(1,3); cmap=hsv(3);
    nexttile([1 2]);
    imagesc(STAmat_AmpNorm(dsort,eigSort),[-0.2 1.5])
    colormap(turbo)
    xlabel('Spike ID')
    % ax2=nexttile([1 1]);
    % corrMat(corrMat==0)=NaN;
    % triuInd=(triu(true(size(corrMat,1)),1));
    % imagesc(corrMat(dsort,dsort))
    % colormap(ax2,hot)
    % colorbar
    nexttile([1 1]);
    
    basalOnMat=(Dsign(Dist_order{f}(noi_dist{f}))<0).*(Dsign(Dist_order{f}(noi_dist{f}))<0)';
    apicalOnMat=(Dsign(Dist_order{f}(noi_dist{f}))>0).*(Dsign(Dist_order{f}(noi_dist{f}))>0)';
    baApOnMat=(Dsign(Dist_order{f}(noi_dist{f}))>0).*(Dsign(Dist_order{f}(noi_dist{f}))<0)';

    scatter(tovec(intdistMat(triuInd & basalOnMat)),tovec(corrMat(triuInd & basalOnMat)),10,cmap(1,:),'o','filled'); hold all
    scatter(tovec(intdistMat(triuInd & apicalOnMat)),tovec(corrMat(triuInd & apicalOnMat)),10,cmap(2,:),'o','filled'); hold all
    scatter(tovec(intdistMat(triuInd & baApOnMat)),tovec(corrMat(triuInd & baApOnMat)),10,cmap(3,:),'o','filled'); hold all
    xlabel('Pairwise distance (\mum)')
    ylabel('Correlation')
    legend({'Among basal','Among apical','Basal and apical'})

    for dClassI=1:5
        for dClassJ=dClassI:5
            if dClassI==dClassJ
                if length(find(labelclass==dClassI))>1
                    corrmat_tmp=corrMat_dend(find(labelclass==dClassI),find(labelclass==dClassJ));
                    [a]=triu(ones(size(corrMat_dend(find(labelclass==dClassI),find(labelclass==dClassI)))),1);
                    Corrcoeff{f}(dClassI,dClassJ)=mean(tovec(corrmat_tmp(a>0)),'omitnan');
                else
                    Corrcoeff{f}(dClassI,dClassJ)=NaN;
                end
            else
                Corrcoeff{f}(dClassI,dClassJ)=mean(tovec(corrMat_dend(find(labelclass==dClassI),find(labelclass==dClassJ))),'omitnan');
            end
        end
    end

    cmap_label=hsv(5);
    cmap_label=cmap_label([4 2 3 5 1],:);
    distMatrix = sqrt(2*(1 - corrMat));
    distMatrix(logical(eye(size(distMatrix)))) = 0;
    Z_ref = linkage(squareform(distMatrix), 'average');
    leafOrder = optimalleaforder(Z_ref,squareform(distMatrix));
    Cluster_ref= switchlabel(cluster(Z_ref, 'maxclust', 5));
    [~, order_cluster]=sort(Cluster_ref,'ascend');

    l=imshow_label(corrMat(leafOrder,leafOrder),labelvec(leafOrder),cmap_label,{'Basal','Soma','Trunk','Oblique','Apical (distal)'});
    axis equal tight off
    if f==foi(end)
    else
        legend off
    end
    caxis([-0.3 1])
    title(num2str(f))
    drawnow;
    colormap(hot(256))
end

figure(119); clf;
label_str={'B','S','T','O','D'}; pair_str=[];
Corr_to_plot=[1 1; 3 3; 4 4; 5 5; 1 2; 1 3; 1 4; 1 5; 2 3; 2 4; 2 5; 3 4; 3 5; 4 5;]; g=1; show_c=[];
for f=foi
    for p=1:size(Corr_to_plot,1)
        show_c{p}(g)=Corrcoeff{f}(Corr_to_plot(p,1),Corr_to_plot(p,2));
        pair_str{p}=[label_str{Corr_to_plot(p,1)} '-' label_str{Corr_to_plot(p,2)}];
    end
    g=g+1;
end



Boxplot_wPoints2(show_c,hsv(size(Corr_to_plot,1)))

%bar([1:size(Corr_to_plot,1)],mean(show_c,1,'omitnan'),'FaceColor',[0.5 .5 .5]); hold all
%plot(1:size(Corr_to_plot,1),show_c','marker','o','color',[0.7 0.7 0.7]); hold all
% x = repelem(1:size(show_c, 2), size(show_c, 1))'; % Group identifiers
% y = show_c(:); % Flatten data for boxchart
% b=boxchart(x, y,'BoxFaceColor', [0.2 0.5 0.9]); hold all
% b.MarkerStyle = 'none'; b.BoxFaceAlpha=0.6;
% b.WhiskerLineColor= [0.1 0.3 0.5]; b.LineWidth= 1.5;
%errorbar(1:size(Corr_to_plot,1),mean(show_c,1,'omitnan'),std(show_c,0,1,'omitnan')./sqrt(sum(~isnan(show_c))),'color','k','linewidth',2,'LineStyle','none','Marker','+')

set(gca,'XTick',[1:size(Corr_to_plot,1)],'XTickLabel',pair_str);
xlim([0.5 size(Corr_to_plot,1)+0.5]);
grid on
ylabel('Correlation coefficient')


%% scatter plot of dSpikes and soma amplitude
f=23;

normTr_show=NormalizedTrace_dirt{f};%-movmedian(NormalizedTrace_dirt{f},3000,2,'omitnan');

roisD_order_ind=cellfun(@find,roisD_order{f},'UniformOutput',false);
if isempty(cell2mat(roisD_order_ind(1,:)'))
    basalind=cell2mat(roisD_order_ind(2,:)'); %if there is no basal, use soma
else
    basalind=cell2mat(roisD_order_ind(1,:)');
end
apicalind=setdiff(cell2mat(roisD_order_ind(5,:)'),38); %apical = distal dend
somaind=cell2mat(roisD_order_ind(2,:)'); %apical = distal dend

SomaTr=mean(normTr_show(somaind,:),1,'omitnan');
DendTr=mean(normTr_show(apicalind,:),1,'omitnan');

[allsta]=get_STA([SomaTr; DendTr],allSpikeMat{f}(1,:),2,3);
[~, dt_ref]=max(allsta,[],2);
dt_ref=dt_ref(2)-dt_ref(1);

DendTr_hi=DendTr-movmedian(DendTr,50,2);
DendTr_peak=find_spike_bh(DendTr_hi,0.002,0.0005).*double(BlueStim{f}==0);
ShiftedTrace=[SomaTr; [DendTr(dt_ref+1:end) zeros(1,dt_ref)]];
[~, peridSpMAt]=get_STA(ShiftedTrace,DendTr_peak,3,3);
%[~, periSpMAt]=get_STA([SomaTr; DendTr],max(allSpikeClassMat{f}).*double(BlueStim{f}==0),1,4);
%[~, STAMat]=get_STA(normTr_show,max(allSpikeClassMat{f}).*double(BlueStim{f}==0),35,15);
[~, periSpMAt]=get_STA(ShiftedTrace,allSpikeMat{f}(1,:).*double(BlueStim{f}==0),3,3);
[~, STAMat sptime]=get_STA(normTr_show,allSpikeMat{f}(1,:).*double(BlueStim{f}==0),50,35);
[~, STAdSPMat dsptime]=get_STA(normTr_show,DendTr_peak,50,35);
[~, STAdSP_SomSpMat]=get_STA(allSpikeMat{f}(1,:),DendTr_peak,2,3);
[~, STASP_SomSpMat]=get_STA(allSpikeMat{f}(1,:),allSpikeMat{f}(1,:).*double(BlueStim{f}==0),50,35);
nonnanSpike=sum(isnan(STAMat),[1 3])==0;
nonnandSpike=sum(isnan(STAdSPMat),[1 3])==0;

[peridSpPeak dt]=max(peridSpMAt,[],3);
[periSpPeak dt2]=max(periSpMAt,[],3);

dt=dt(:,nonnandSpike);
dt2=dt2(:,nonnanSpike);
dsptime=dsptime(nonnandSpike);
sptime=sptime(nonnanSpike);
periSpPeak=periSpPeak(:,nonnanSpike);
peridSpPeak=peridSpPeak(:,nonnandSpike);

figure(33); clf;
%[Count cx cy]=scatter_heatmap2(periSpPeak(2,:),periSpPeak(1,:),[-0.01:0.003:0.025],[-0.01:0.003:0.025]);
%[dCount dcx dcy]=scatter_heatmap2(peridSpPeak(1,:),peridSpPeak(2,:),[-0.01:0.001:0.025],[-0.01:0.001:0.025]); hold all
scatter(peridSpPeak(2,:),peridSpPeak(1,:),10,[0.5 0.5 1],'o','filled'); hold all
scatter(periSpPeak(2,:),periSpPeak(1,:),35,[1 0.7 0],'o','MarkerEdgeAlpha',0.7); hold all

% [X, Y] = meshgrid(cx, cy);
% [c, h]=contour(Y, X, Count,[-0.001 0 0.001 0.005:0.015:0.1],'LineWidth',1,'LabelSpacing', 200);
% colormap(turbo)
% clabel(c,h,'color','k','Fontsize',12,'LabelSpacing', 800)
xlabel('Apical dendrite peak amplitude (\DeltaF/F)')
ylabel('Soma amplitude (\DeltaF/F)')
legend({'Dendritic peak','bAP'})

interactive_scatter_matrix_viewer(STAMat(:,nonnanSpike,:), STAdSPMat(:,nonnandSpike,:), periSpPeak([2 1],:), peridSpPeak([2 1],:))
%interactive_multi_select_viewer(STAMat, STAdSPMat, periSpPeak([2 1],:), peridSpPeak([2 1],:))

figure(34); clf;
histogram(dt2(2,:)-1,[-1:4],'FaceColor',[1 0.7 0]);
xlabel('Dendrite peak time - soma peak time (ms)')
ylabel('Number of spikes')

        figure(35); clf; tiledlayout(1,3);
        nexttile([1 1])
        imagesc(squeeze(STAMat(:,353,:)),cax)
        
        nexttile([1 1])
        plot([-50:35],squeeze(mean(STAdSPMat(somaind,2742,:),1,'omitnan'))); hold all
        plot([-50:35],squeeze(mean(STAdSPMat(basalind,2742,:),1,'omitnan'))); hold all
        plot([-50:35],squeeze(mean(STAdSPMat(apicalind,2742,:),1,'omitnan'))); hold all
        
        nexttile([1 1])
        plot([-50:35],squeeze(mean(STAdSPMat(somaind,7194,:),1,'omitnan'))); hold all
        plot([-50:35],squeeze(mean(STAdSPMat(basalind,7194,:),1,'omitnan'))); hold all
        plot([-50:35],squeeze(mean(STAdSPMat(apicalind,7194,:),1,'omitnan'))); hold all

% figure(34); clf;
% % [Count cx cy]=scatter_heatmap2(periSpPeak(2,:),dt2(2,:)-2,[-0.01:0.003:0.025],[-3:5]);
% scatter(peridSpPeak(2,nonnandSpike),4-dt(1,nonnandSpike),10,[0.5 0.5 1],'o','filled'); hold all
% scatter(periSpPeak(2,nonnanSpike),dt2(2,nonnanSpike)-2,50,[1 0 0],'o','AlphaData',0.5);
% % [X, Y] = meshgrid(cx, cy);
% % [c, h]=contour(X', Y', Count,[-0.001 0 0.001 0.005:0.02:0.1],'LineWidth',1,'LabelSpacing', 200);
% % colormap(turbo)
% % clabel(c,h,'color','k','Fontsize',12,'LabelSpacing', 800)
% xlabel('Apical dendrite peak amplitude (\DeltaF/F)')
% ylabel('Dendrite peak time - soma peak time (ms)')
% ylim([-2.5 3.5])







