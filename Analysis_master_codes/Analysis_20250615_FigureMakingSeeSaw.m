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
%% Concatenate data
nTau={[-200:50],[-200:50],[-200:50]}; %SS, CS, Brst
clear SpikeInd MatSpike MatSTA MatBlue MatCStrace MatSub SpikeList NormalizedTrace_ch NormalizedTrace_dirt SpikeIndBlueOff Dist_order allSpikeMat noi interDendDist noi_dist derivSub LapSubSilent
clear Subthreshold dendaxis BrstOrder roisD roisD_order LapSpclassVec dendaxis1d
clear SI MI

% Load local stimulation data
saveto='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism';
load([saveto '/LocalStimulationResult.mat']);

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
    nROI=length(noi); nTime=size(Result.normTraces,2);
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
    SWC{f}=Result.SWC;

    roisD_order_ind=cellfun(@find,roisD_order{f},'UniformOutput',false);
    labelvec{f}=NaN(1,nROI);
    for dClass=1:5
        labelvec{f}(cell2mat(roisD_order_ind(dClass,:)'))=dClass;
    end

    dendaxis{f}=interDendDist{f}(1,:).*Dsign;
    dendaxis{f}=dendaxis{f}(Dist_order{f}(noi_dist{f}));

    ftcoord=get_coord(Result.ftprnt);
    dendaxis1d{f} = projectTrunkaxis(ftcoord);
    dendaxis1d{f} = dendaxis1d{f}(Dist_order{f}(noi_dist{f}));

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

    % SI{f}.FRreal=Result.SI_FRreal;
    % SI{f}.FRnull=Result.SI_FRnull;
    % SI{f}.FRclassReal=Result.SI_FRClassReal;
    % SI{f}.FRclassNull=Result.SI_FRClassnull;
    % 
    % MI{f}.EigTrPos=Result.MI_EigTrReal;
    % MI{f}.EigTrFR=Result.MI_EigTrFR;
    % MI{f}.EigTrFRdClass=Result.MI_EigTrFRdClass;
    % MI{f}.SubVFR=Result.MI_SubVFR;
    % MI{f}.SubVFRdClass=Result.MI_SubVFRdClass;
    % MI{f}.SubVPos=Result.MI_SubVReal;

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

    % LapFR{f}=PlaceTrigger_average(double(allSpikeMat{f}(1,:)==1),150,VRtrack{f},0.002,115); %total trace
    % LapV{f}=PlaceTrigger_average(NormalizedTrace_dirt{f},150,VRtrack{f},0.002,115); %total trace
    % LapSub{f}=PlaceTrigger_average(Subthreshold{f},150,VRtrack{f},0.002,115); %total trace
    % LapSpclass{f}=PlaceTrigger_average(double(allSpikeClassMat{f}>0),150,VRtrack{f},-0.002,115); %total trace
    % LapSpclassVec{f}=PlaceTrigger_average(double(allSpikeClassVecMat{f}>0),150,VRtrack{f},-0.002,115); %total trace

    subthreshold_silent=Subthreshold{f};
    subthreshold_silent(:,unique(get_perispikeIndex(allSpikeMat{f}(1,:),[-8:30])))=NaN;
    %LapSubSilent{f}=PlaceTrigger_average(subthreshold_silent,150,VRtrack{f},0.002,115); %total trace
end

MatSTA=cellfun(@(x) permute(x,[1 3 2]),MatSTA,'UniformOutput',false);
MatSpClass=cellfun(@(x) permute(x,[1 3 2]),MatSpClass,'UniformOutput',false);
MatSpike=cellfun(@(x) permute(x,[1 3 2]),MatSpike,'UniformOutput',false);
MatBlue=cellfun(@(x) permute(x,[1 3 2]),MatBlue,'UniformOutput',false);
MatSub=cellfun(@(x) permute(x,[1 3 2]),MatSub,'UniformOutput',false);
MatOrder=cellfun(@(x) permute(x,[1 3 2]),MatOrder,'UniformOutput',false);
MatCStrace=cellfun(@(x) permute(x,[1 3 2]),MatCStrace,'UniformOutput',false);

%% show STA kymo
figure(101); clf; tiledlayout(9,9)
stype_str={'SS','CS','BS'};
sub_time=[-19:-5]; %from spike
ax1=[];
for f=foi
    [~, BlueOffSpike]=find_wCondition(SpikeList{f,1,1},BlueStim{f}==0);
    refkymo=mean(MatSTA{f,1,1}(:,:,BlueOffSpike)-median(MatSTA{f,1,1}(:,-nTau{1}(1)+sub_time,BlueOffSpike),2,'omitnan'),3,'omitnan');
    cax=[prctile(refkymo(:),0.5) prctile(refkymo(:),99.5)];

    Dsign=ones(1,size(interDendDist{f},1));
    Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
    for stype=1:3
        [~, BlueOffSpike]=find_wCondition(SpikeList{f,stype,1},BlueStim{f}==0);
        if ~isempty(BlueOffSpike) % If there is no spike
            imshow=mean(MatSTA{f,stype,1}(:,:,BlueOffSpike)-median(MatSTA{f,stype,1}(:,-nTau{stype}(1)+sub_time,BlueOffSpike),2,'omitnan'),3,'omitnan');
        else
            imshow=zeros(sum(noi_dist{f}),length(nTau{stype}));
        end
        %imshow=mean(STAmat{f,stype}(:,:,SpikeIndBlueOff{f,sty  pe}),3,'omitnan');
        if isnan(cax); cax=[0 1]; end;
        ax1=[ax1 nexttile([1 1])];
        imagesc(imshow,cax)
        title(['Cell# ' num2str(f) ', ' stype_str{stype} ', N=' num2str(length(BlueOffSpike))])
        if find(Dist_order{f}(noi_dist{f})==1)~=1
            set(gca,'YTick',[1 find(Dist_order{f}(noi_dist{f})==1) sum(noi_dist{f})],'YTickLabel',num2str([min(interDendDist{f}(1,:).*Dsign) 0 max((interDendDist{f}(1,:).*Dsign))]',3))
        else
            set(gca,'YTick',[1 sum(noi_dist{f})],'YTickLabel',num2str([0 max((interDendDist{f}(1,:)))]',3))
        end
        set(gca,'XTick',[1 -nTau{stype}(1)+1 length(nTau{stype})],'XTickLabel',num2str([nTau{stype}(1) 0 nTau{stype}(end)]',3))
    end
end
linkaxes(ax1,'x')
colormap(turbo)

figure(100); clf;
for f=foi
    nexttile([1 1])
    imshow2(AvgImg{f},[])
    title(num2str(f))
end
%% See-saw and Global power during blue on and off

f=27;
[nROI nTime]=size(Subthreshold{f});
[SubV SubD subTrace]=get_eigvector(Subthreshold{f}(:,sum(isnan(Subthreshold{f}),1)==0),nROI);
subTrace_onFrame=NaN(nROI,nTime);
subTrace_onFrame(:,sum(isnan(Subthreshold{f}),1)==0)=subTrace';
perispike_time=[-2:20]; 
perispike_frame=ind2vec(nTime,unique([tovec(find(double(allSpikeMat{f}(1,:)==1))'+perispike_time); find(CStrace{f})']),1);

Blue_on_frame=imdilate(BlueStim{f}>0, [ones(1, 1), 1, ones(1, 200)])>0;
Blue_off_frame=imdilate(BlueStim{f}==0, [1, ones(1, 50)])>0;

roisD_order_ind=cellfun(@find,roisD_order{f},'UniformOutput',false);
if isempty(cell2mat(roisD_order_ind(1,:)'))
basalind=cell2mat(roisD_order_ind(2,:)'); %if there is no basal, use soma
else
basalind=cell2mat(roisD_order_ind(1,:)');
end
apicalind=cell2mat(roisD_order_ind(5,:)'); %apical = distal dend 

SubBasalApical=[mean(Subthreshold{f}(basalind,:),1,'omitnan'); mean(Subthreshold{f}(apicalind,:),1,'omitnan')];

figure(301); clf; tiledlayout(2,3);
for pc=1:3
nexttile([1 1])
h=histogram(subTrace_onFrame(pc,Blue_off_frame & ~perispike_frame),'Normalization','probability'); hold all
h2=histogram(subTrace_onFrame(pc,Blue_on_frame & ~perispike_frame),h.BinEdges,'Normalization','probability'); hold all
h.FaceColor=[0.2 0.2 0.2]; h2.FaceColor=[0 0.5 1];
h.EdgeAlpha=0; h2.EdgeAlpha=0;
xlabel('Projection onto eigvector')
ylabel('Probability')
title(['PC #' num2str(pc)])
end

ax1=nexttile([1 1]);
[EigTrHistoOff, x_edges, y_edges] = scatter_heatmap(subTrace_onFrame(1,Blue_off_frame & ~perispike_frame),subTrace_onFrame(2,Blue_off_frame & ~perispike_frame),50); hold all
colormap(ax1,'hot')
ylabel('See-saw axis')
xlabel('Global axis')
title('Blue off')
ax2=nexttile([1 1]);
[EigTrHistoOn, x_center, y_center] = scatter_heatmap2(subTrace_onFrame(1,Blue_on_frame & ~perispike_frame),subTrace_onFrame(2,Blue_on_frame & ~perispike_frame),x_edges,y_edges);
colormap(ax2,'parula')
ylabel('See-saw axis')
xlabel('Global axis')
title('Blue on')
nexttile([1 1])
[X, Y]=meshgrid(x_center,y_center);
contour_cmap(Y,X,EigTrHistoOff, 10, hot(10)); hold all
contour_cmap(Y,X,EigTrHistoOn, 10, parula(10));
ylabel('See-saw axis')
xlabel('Global axis')

%% Show Global, See-saw power during blue On/off of multiple neurons
figure(306); clf;
for f=[20 21 22 25 26 27]
[nROI nTime]=size(Subthreshold{f});
[SubV SubD subTrace]=get_eigvector(Subthreshold{f}(:,sum(isnan(Subthreshold{f}),1)==0)',nROI);
subTrace_onFrame=NaN(nROI,nTime);
subTrace_onFrame(:,sum(isnan(Subthreshold{f}),1)==0)=subTrace';
perispike_time=[-2:20]; 
perispike_frame=ind2vec(nTime,unique([tovec(find(double(allSpikeMat{f}(1,:)==1))'+perispike_time); find(CStrace{f})']),1);

Blue_on_frame=imdilate(BlueStim{f}>0, [ones(1, 1), 1, ones(1, 200)])>0;
Blue_off_frame=imdilate(BlueStim{f}==0, [1, ones(1, 50)])>0;

[EigTrHistoOff, x_edges, y_edges] = hist_heatmap(subTrace_onFrame(1,Blue_off_frame & ~perispike_frame),subTrace_onFrame(2,Blue_off_frame & ~perispike_frame),50); hold all
[EigTrHistoOn, x_center, y_center] = hist_heatmap2(subTrace_onFrame(1,Blue_on_frame & ~perispike_frame),subTrace_onFrame(2,Blue_on_frame & ~perispike_frame),x_edges,y_edges);
nexttile([1 1])
[X, Y]=meshgrid(x_center,y_center);
contour_cmap(Y,X,EigTrHistoOff, 10, hot(10)); hold all
contour_cmap(Y,X,EigTrHistoOn, 10, parula(10));
ylabel('See-saw axis')
xlabel('Global axis')
title(['Neuron ID# ' num2str(f)])
end

%% Scatter plot, global-seesaw, B+A/B-A/S
figure(307); clf;  tiledlayout(4,3); g=1;
for f=[20 23 25 26]

[nROI nTime]=size(Subthreshold{f}); perispike_time=[-2:20]; 
[SubV SubD subTrace]=get_eigvector(Subthreshold{f}(:,sum(isnan(Subthreshold{f}),1)==0)',nROI);
subTrace_onFrame=NaN(nROI,nTime);
subTrace_onFrame(:,sum(isnan(Subthreshold{f}),1)==0)=subTrace';

perispike_frame=ind2vec(nTime,unique([tovec(find(double(allSpikeMat{f}(1,:)==1))'+perispike_time); find(CStrace{f})']),1);
Blue_on_frame=imdilate(BlueStim{f}>0, [ones(1, 1), 1, ones(1, 200)])>0;
Blue_off_frame=imdilate(BlueStim{f}==0, [1, ones(1, 50)])>0;

roisD_order_ind=cellfun(@find,roisD_order{f},'UniformOutput',false);
if isempty(cell2mat(roisD_order_ind(1,:)'))
basalind=cell2mat(roisD_order_ind(2,:)'); %if there is no basal, use soma
else
basalind=cell2mat(roisD_order_ind(1,:)');
end
apicalind=cell2mat(tovec(roisD_order_ind([3 5],:))); %apical = distal dend 

SubBasalApical=[mean(Subthreshold{f}(basalind,:),1,'omitnan'); mean(Subthreshold{f}(apicalind,:),1,'omitnan')];

    ax1=[]; Prespike_t=5;
nexttile([1 1])
[hist_subV_GlSS]=scatter_heatmap(subTrace_onFrame(1,Blue_off_frame),subTrace_onFrame(2,Blue_off_frame),50);
[~, subTraceSSmat]=get_STA(subTrace_onFrame,allSpikeClassMat{f}(1,:).*Blue_off_frame,Prespike_t,-Prespike_t);
[~, subTraceCSmat]=get_STA(subTrace_onFrame,allSpikeClassMat{f}(2,:).*Blue_off_frame,Prespike_t,-Prespike_t);
hold all
scatter(subTraceSSmat(1,:),subTraceSSmat(2,:),10,[1 0 1],'filled');
scatter(subTraceCSmat(1,:),subTraceCSmat(2,:),10,[0 1 1],'filled');
xlabel('Global'); ylabel('See-saw'); legend({'SS','CS'})

ax1=[ax1 nexttile([1 1])];
BAsumdiffTrace=[sum(SubBasalApical,1); diff(SubBasalApical,1,1)];
[hist_subV_ApBa]=scatter_heatmap(BAsumdiffTrace(1,Blue_off_frame),BAsumdiffTrace(2,Blue_off_frame),50);
[~, subTraceSSmat]=get_STA(BAsumdiffTrace,allSpikeClassMat{f}(1,:).*Blue_off_frame,Prespike_t,-Prespike_t);
[~, subTraceCSmat]=get_STA(BAsumdiffTrace,allSpikeClassMat{f}(2,:).*Blue_off_frame,Prespike_t,-Prespike_t);
hold all
scatter(subTraceSSmat(1,:),subTraceSSmat(2,:),10,[1 0 1],'filled');
scatter(subTraceCSmat(1,:),subTraceCSmat(2,:),10,[0 1 1],'filled')
xlabel('Basal + apical'); ylabel('Apical - basal'); legend({'SS','CS'})

nexttile([1 1])
[hist_subV_A_B]=scatter_heatmap(SubBasalApical(1,Blue_off_frame),SubBasalApical(2,Blue_off_frame),50);
[~, subTraceSSmat]=get_STA(SubBasalApical,allSpikeClassMat{f}(1,:).*Blue_off_frame,Prespike_t,-Prespike_t);
[~, subTraceCSmat]=get_STA(SubBasalApical,allSpikeClassMat{f}(2,:).*Blue_off_frame,Prespike_t,-Prespike_t);
hold all
scatter(subTraceSSmat(1,:),subTraceSSmat(2,:),10,[1 0 1],'filled');
scatter(subTraceCSmat(1,:),subTraceCSmat(2,:),10,[0 1 1],'filled')
xlabel('Basal'); ylabel('Apical'); legend({'SS','CS'})
colormap('bone')

title(ax1,['Neuron ID# ' num2str(f)])
end

%% Mutual information
SI_zscore=NaN(max(foi),1); SI_SpClass_zscore=NaN(max(foi),4); SI_value=NaN(max(foi),5);
MI_EigTrPos=NaN(max(foi),10); MI_EigTrFR=NaN(10,5,max(foi)); FR=NaN(max(foi),5);
MI_SubVPos=[]; MI_SubVFR=[];
for f=foi(1:end-2)
    SI_value(f,:)=[SI{f}.FRreal SI{f}.FRclassReal'];
    SI_zscore(f)=(SI{f}.FRreal-mean(SI{f}.FRnull))/std(SI{f}.FRnull);
    SI_SpClass_zscore(f,:)=(SI{f}.FRclassReal-mean(SI{f}.FRclassNull,2))'./std(SI{f}.FRclassNull,0,2)';
    MI_EigTrPos(f,:)=MI{f}.EigTrPos;
    MI_EigTrFR(:,:,f)=[MI{f}.EigTrFR MI{f}.EigTrFRdClass'];
    MI_SubVPos{f}=[MI{f}.SubVPos];
    MI_SubVFR{f}=[MI{f}.SubVFR MI{f}.SubVFRdClass'];
    MI_SubVPos{f}((MI_SubVPos{f}<0))=NaN; MI_SubVFR{f}(MI_SubVFR{f}<0)=NaN;
    MI_SubVPos{f}=[dendaxis{f}' MI{f}.SubVPos];
    MI_SubVFR{f}=[dendaxis{f}' MI{f}.SubVFR MI{f}.SubVFRdClass'];
    BlueOffFrame=BlueStim{f}==0;
    RunTime=VRtrack{f}(end,:)>0.001;
    StimOnLap=unique(VRtrack{f}(8,:).*(BlueStim{f}>0));
    TrainingEnvLap=unique(VRtrack{f}(8,:).*(VRtrack{f}(2,:)==1));
    StimOfflap=setdiff(unique(VRtrack{f}(8,:)),[StimOnLap TrainingEnvLap]);
    FrameOfInterest= BlueOffFrame & RunTime & ismember(VRtrack{f}(8,:),StimOfflap);
    FR(f,:)=[mean(allSpikeMat{f}(1,FrameOfInterest),'omitnan') mean(allSpikeClassMat{f}(:,FrameOfInterest),2,'omitnan')']*1000;
end
MI_EigTrPos((MI_EigTrPos<0))=NaN; MI_EigTrFR(MI_EigTrFR<0)=NaN;


figure(308); clf;
nexttile([1 1]);
p1=Boxplot_wPoints2(SI_value,hsv(5));
set(gca,'XTick',[1:5],'XTickLabel',{'AP','SS','CS','BS','dSP'});
ylabel('Spatial information (bits/spike)')
nexttile([1 1]);
p2=Boxplot_wPoints2([SI_zscore(foi) SI_SpClass_zscore(foi,:)],hsv(5)); hold all
plot([0.5 5.5],[1.96 1.96],'k--')
%drawPValueLines(p,-2,'StepHeight',3)
set(gca,'XTick',[1:5],'XTickLabel',{'AP','SS','CS','BS','dSP'});
ylabel('Z-scored spatial information')

figure(309); clf;
nexttile([1 1])
p1=Boxplot_wPoints2(MI_EigTrPos(foi,:),hsv(10));
drawPValueLines(p1,0,'TextYOffset',0.01,'StepHeight',0.05)
set(gca,'XTick',[1:10],'XTickLabel',[1:10]);
xlabel('Principle components')
ylabel('Mutual information (bits)')

figure(311); clf;
nexttile([1 1])
MIsubVPosMat=cell2mat(MI_SubVPos(foi)');
scatter(MIsubVPosMat(:,1),MIsubVPosMat(:,2),10,[0.6 0.6 0.6],'filled'); hold all
%Norm_SubVPos=cellfun(@(x) [x(:,1) x(:,2)/mean(x(x(:,1)>-40 & x(:,1)<40,2),'omitnan')],MI_SubVPos(foi),'UniformOutput',false);
[mean_amplitudes std_amplitudes x_bin_centers indicies]=binning_data(MI_SubVPos(foi),[-300 -150 -50 50 150 250 400 550]);
errorbar(x_bin_centers,mean_amplitudes,std_amplitudes./length(foi),'color',[1 0 0]);
xlabel('Principle components')
ylabel('Mutual information (bits)')
title('Subthreshold vs VR position')
ylim([0 0.15])

nexttile([1 1]); 
cmap=hsv(5);
for dclass=1:5
    show_dat=cellfun(@(x) x(:,[1 dclass+1]),MI_SubVFR(foi),'UniformOutput',false);
    [mean_amplitudes std_amplitudes x_bin_centers indicies]=binning_data(show_dat,[-300 -150 -50 50 150 250 400 550]);
    errorbar(x_bin_centers,mean_amplitudes,std_amplitudes./length(foi),'color',cmap(dclass,:)); hold all
end
legend({'AP','SS','CS','BS','dSP'})
xlabel('Distance from soma')
ylabel('Mutual information (bits)')
title('Subthreshold vs Firing rate')

figure(310); clf; tiledlayout(1,4);
dclass_str={'AP','SS','CS','BS','dSP'};
for dclass=1:4
nexttile([1 1])
showfoi=FR(foi,dclass)>0;
p1=Boxplot_wPoints2(squeeze(MI_EigTrFR([1:10],dclass,foi(showfoi)))',hsv(10));
drawPValueLines(p1,0,'TextYOffset',0.01,'StepHeight',0.05)
set(gca,'XTick',[1:10],'XTickLabel',[1:10]);
xlabel('PC component #')
set(gca,'yscale','log'); box off;
title(dclass_str{dclass})
ylabel('Mutual information (bits)')
xlim([0.5 10.5])
end

%% PCA analysis, plot variance and eigenvector map

%f=10;
cumsumD=[]; SubV=[]; SubD=[];
figure(217); clf;
for f=foi
    f
    nTime=size(Subthreshold{f},2);
    nROI=size(Subthreshold{f},1);
    [SubV{f} SubD{f} subTrace]=get_eigvector(Subthreshold{f}(:,sum(isnan(Subthreshold{f}),1)==0),nROI);
    cumsumD{f}=[[1:nROI]' [cumsum(SubD{f})/sum(SubD{f})]];  
    plot(cumsumD{f}(:,2),'color',[0.7 0.7 0.7]); hold all
    drawnow;
end
[M S x]=binning_data(cumsumD,[0.5:1:10.5]);
errorbar(x,M,S,'color','k');
xlabel('PC component #'); box off;
ylabel('Fraction of variance')

show_f=[1 4 15 18 20 23 25 26 27];
figure(216); clf; tiledlayout(length(show_f),5);
for f=show_f
for pc=1:5%size(subTrace_onFrame,1)
    
    nexttile([1 1])
    if f>=20
    ScaleImage=showScaleImage(fliplr(Ftprnts{f}(:,:,Dist_order{f}(noi_dist{f}))>0), (SubV{f}(:,pc)),gen_colormap([0 0.5 1; 1 1 1; 1 0 0]),[prctile(SubV{f}(:),1) prctile(SubV{f}(:),99)]);    
    else
    ScaleImage=showScaleImage(Ftprnts{f}(:,:,Dist_order{f}(noi_dist{f}))>0, (SubV{f}(:,pc)),gen_colormap([0 0.5 1; 1 1 1; 1 0 0]),[prctile(SubV{f}(:),1) prctile(SubV{f}(:),99)]);
    end
    %ScaleImage(ScaleImage==0)=1;
    imshow2(ScaleImage,[])
    se = strel('disk', 5); % Define a structuring element

    FrprntBW=max(Ftprnts{f}(:,:,Dist_order{f}(noi_dist{f})),[],3)>0;
    se90 = strel('line',2,90);
    se0 = strel('line',2,0);
    FrprntBW = imdilate(FrprntBW,se);

    b=bwboundaries(FrprntBW);
    hold all
    b=cell2mat(b);
    %plot(b(:,2),b(:,1),'color',[0.5 0.5 0.5])
    axis equal tight off
    if pc==0
        title(['N# ' num2str(f) ', Sub., ' num2str(SubD{f}(pc)/sum(SubD{f})*100,2),' %'])
    else
        title(['Neuron ID: #' num2str(f) ', PC#' num2str(pc) newline 'fraction: ' num2str(SubD{f}(pc)/sum(SubD{f})*100,2),' %'])
    end
end
cb=colorbar;
colormap(gen_colormap([0 0.5 1; 1 1 1; 1 0 0]));
cb.Ticks=[0 1];
cb.TickLabels=num2str([prctile(SubV{f}(:),1); prctile(SubV{f}(:),99)],2);
end
% LapSub=PlaceTrigger_average(subTrace_onFrame,150,VRtrack{f},-0.002,115)*1000; %total trace
% figure(217); clf;
% imshow_patch(ringmovMean(LapSub,3))
%% Subthreshold map vs p-value map
nBin=100;
figure(237); clf; tiledlayout(13,length(foi)/2); g=1;
for f=foi(11:20)
    StimOnLap=unique(VRtrack{f}(8,:).*(BlueStim{f}>0));
    TrainingEnvLap=unique(VRtrack{f}(8,:).*(VRtrack{f}(2,:)==1));
    StimOfflap=setdiff(unique(VRtrack{f}(8,:)),[StimOnLap TrainingEnvLap]);
    ax1=nexttile(g,[5 1]);
    imagesc(permute(mean(LapSub{f}(StimOfflap,:,:),1,'omitnan'),[3 2 1]))
    colormap(ax1,turbo)
    ax2=nexttile(g+length(foi)/2*5,[5 1]);
    imagesc(ShuffleV.Pvalue{f},[0 1])
    colormap(ax2,"parula")
    title(['Neuron #' num2str(f)])
    ax3=nexttile(g+length(foi)/2*10,[3 1]);
    plot(mean(LapFR{f}(StimOfflap,:),1,'omitnan')*1000,'r')
    g=g+1;
end

% FR, Sub, EigTrace p-value map
figure(236); clf; tiledlayout(22,length(foi)/2); g=1;
for f=foi(11:20)
    nexttile(g,[4 1]);
    imagesc(ShuffleV.FRPvalue{f}',[0 1]); axis on;
    nexttile(g+length(foi)/2*4,[7 1]);
    imagesc(ShuffleV.PCPvalue{f},[0 1]); axis off;
    nexttile(g+length(foi)/2*11,[11 1]);
    imagesc(ShuffleV.Pvalue{f},[0 1]); axis off;
    title(['Neuron #' num2str(f)])
    g=g+1;
end

figure(237); clf;
p_threshold=0.05;
Fieldsize_threshold=5;
for f=foi
    StimOnLap=unique(VRtrack{f}(8,:).*(BlueStim{f}>0));
    TrainingEnvLap=unique(VRtrack{f}(8,:).*(VRtrack{f}(2,:)==1));
    StimOfflap=setdiff(unique(VRtrack{f}(8,:)),[StimOnLap TrainingEnvLap]);

FRint=interp1([1:150]/150*2,1000*mean(LapFR{f}(StimOfflap,:),1,'omitnan'),[1:nBin]/nBin*2);
FRpmat=ShuffleV.FRPvalue{f}';
FRpmat(end)=mean(FRpmat([1 end-1]));

[Fsize Fpos]=get_CircularField(FRpmat<p_threshold,3);
[FsizeD FposD FmatD]=get_CircularField(ShuffleV.Pvalue{f}<p_threshold,3);
[~, ~, FmatD_hyp]=get_CircularField(ShuffleV.Pvalue{f}>1-p_threshold,3);

meanFRcand=[]; cand_bin=[];
for fcand=1:size(Fpos{1},1)
    if Fpos{1}(fcand,1)<=Fpos{1}(fcand,2)
    cand_bin=[cand_bin [Fpos{1}(fcand,1):Fpos{1}(fcand,2)]];
    meanFRcand(fcand)=mean(FRint(Fpos{1}(fcand,1):Fpos{1}(fcand,2)));
    else
    meanFRcand(fcand)=mean(FRint(setdiff([1:nBin],[Fpos{1}(fcand,2)+1:Fpos{1}(fcand,1)-1])));    
    cand_bin=[cand_bin setdiff([1:nBin],[Fpos{1}(fcand,2)+1:Fpos{1}(fcand,1)-1])];    
    end
end
meanFRcand=meanFRcand/mean(FRint(setdiff([1:nBin],cand_bin)));

if ~isempty(find(Fsize{1}>=Fieldsize_threshold))
    Fieldind=find((Fsize{1}>=Fieldsize_threshold) & (meanFRcand'>4));

    for field=Fieldind'
if Fpos{1}(field,1)>Fpos{1}(field,2)
FieldPos=[Fpos{1}(field,1) Fpos{1}(field,2)+nBin];
else
FieldPos=[Fpos{1}(field,1) Fpos{1}(field,2)];
end
PFcenter=round(mean(FieldPos));
Fieldwidth=diff(FieldPos);

pvalMat=repmat(FmatD>0,1,3);
pvalMat=pvalMat(:,PFcenter+nBin+[-nBin/2+1:nBin/2]);

pvalMat_hyp=repmat(FmatD_hyp>0,1,3);
pvalMat_hyp=pvalMat_hyp(:,PFcenter+nBin+[-nBin/2+1:nBin/2]);

FR=repmat(FRint,1,3);
FR=FR(PFcenter+nBin+[-nBin/2+1:nBin/2]);
nexttile([1 1])
imagesc([ind2vec(100,50+[-round(Fieldwidth/2)+1:round(Fieldwidth/2)],2);(pvalMat + pvalMat_hyp*3)])
hold all
yyaxis right
plot(FR,'r-')
title(num2str(f))
    end

end
end
%% Show p-value maps

figure(230); clf;
for f=foi
    StimOnLap=unique(VRtrack{f}(8,:).*(BlueStim{f}>0));
    TrainingEnvLap=unique(VRtrack{f}(8,:).*(VRtrack{f}(2,:)==1));
    StimOfflap=setdiff(unique(VRtrack{f}(8,:)),[StimOnLap TrainingEnvLap]);

    nexttile([1 1])
imagesc([1:100]/100*2,[1:size(ShuffleV.Pvalue{f},1)],ShuffleV.Pvalue{f}); hold all
yyaxis right
plot([1:150]/150*2,mean(LapFR{f}(StimOfflap,:),1,'omitnan'),'r-')
title([num2str(f) ', p-value Map'])
caxis([0 1])
end

figure(234);clf;
for f=foi
    StimOnLap=unique(VRtrack{f}(8,:).*(BlueStim{f}>0));
    TrainingEnvLap=unique(VRtrack{f}(8,:).*(VRtrack{f}(2,:)==1));
    StimOfflap=setdiff(unique(VRtrack{f}(8,:)),[StimOnLap TrainingEnvLap]);
    nexttile([1 1])
    plot([1:150]/150*2,mean(LapFR{f}(StimOfflap,:),1,'omitnan'),'r-'); hold all
    yyaxis right
    plot([1:150]/150*2,mean(LapSub{f}(StimOfflap,:,:),[1 3],'omitnan'),'k-')
end

figure(231); clf;
for f=foi
    StimOnLap=unique(VRtrack{f}(8,:).*(BlueStim{f}>0));
    TrainingEnvLap=unique(VRtrack{f}(8,:).*(VRtrack{f}(2,:)==1));
    StimOfflap=setdiff(unique(VRtrack{f}(8,:)),[StimOnLap TrainingEnvLap]);

    nexttile([1 1])
imagesc([1:100]/100*2,[1:size(LapSub{f},3)],permute(mean(LapSub{f}(StimOfflap,:,:),1,'omitnan'),[3 2 1])); hold all
yyaxis right
plot([1:150]/150*2,mean(LapFR{f}(StimOfflap,:),1,'omitnan'),'r-')
plot((1:150)/150*2,mean(LapSub{f}(StimOfflap,:),1,'omitnan'),'k-')
end

figure(232); clf;
for f=foi
    corrMat_pval=get_corrMat(ShuffleV.Zscore{f},ShuffleV.Zscore{f});
    corrMat_Sub=get_corrMat(Subthreshold{f},Subthreshold{f});
    nexttile([1 1]);
    imagesc([corrMat_pval corrMat_Sub],[-0.2 1]);
    colormap(gen_colormap([0 0.5 1; 1 1 1; 1 0 0]))
    title(['Neuron #' num2str(f)])
end
cb=colorbar;
cb.Label.String=['Correlation'];

figure(235); clf;
CorrPvalFR=[]; CorrPvalFRMat=[];
cmap_label=hsv(5); cmap_label=cmap_label([4 2 3 5 1],:);
for f=foi
    StimOnLap=unique(VRtrack{f}(8,:).*(BlueStim{f}>0));
    TrainingEnvLap=unique(VRtrack{f}(8,:).*(VRtrack{f}(2,:)==1));
    StimOfflap=setdiff(unique(VRtrack{f}(8,:)),[StimOnLap TrainingEnvLap]);

    TC_Sub=permute(mean(LapSub{f}(StimOfflap,:,:),1,'omitnan'),[3 2 1]);
    TC_FR=mean(LapFR{f}(StimOfflap,:),1,'omitnan');
    TC_FR_int=interp1([1:150]/150*2,TC_FR,[1:100]/100*2,'linear','extrap');
    %CorrMatPvalFR=corr(ShuffleV.Pvalue{f}',TC_FR_int');
    CorrMatPvalFR=corr(TC_Sub',TC_FR');

    roisD_order_ind=cellfun(@find,roisD_order{f},'UniformOutput',false);
    for r_ind=1:5
    CorrPvalFR{f,r_ind}=CorrMatPvalFR(cell2mat(roisD_order_ind(r_ind,:)'));
    end
end
for r=1:5
CorrPvalFRMat(:,r)=cell2mat(cellfun(@nanmean, CorrPvalFR(foi,r),'UniformOutput',false));
end
pval=Boxplot_wPoints2(CorrPvalFRMat,cmap_label);
set(gca,'XTick',[1:5],'XTickLabel',{'Basal','Soma','Trunk','Oblique','Distal'});
drawPValueLines(pval,-0.1,'TextYOffset',0.07);
ylabel('Corr. between Subth. TC & firing TC')

%% Correlation matrix between ROIs
figure(118); clf
perispike_time=[-3:10];
Corrcoeff=[]; corrMat=[]; labelvec=[];
for f=foi
    nTime=size(Subthreshold{f},2);
    sub_ch=[];
    for ch=1:2
        sub_ch{ch}=get_subthreshold(NormalizedTrace_ch{f,ch},max(allSpikeMat{f}(1,:),[],1)>0,7,17);
    end
    perispike_frame=unique([tovec(find(double(allSpikeMat{f}(1,:)==1))'+perispike_time); find(CStrace{f})']);
    perispike_frame(perispike_frame<=0 | perispike_frame>nTime)=[];
    silenttime_vec=ind2vec(nTime,perispike_frame,0,1);
    nonvalid_frame=find(sum(isnan(Subthreshold{f}),1)>0);

    Blue_on_frame=find(imdilate(BlueStim{f}>0, [ones(1, 1), 1, ones(1, 200)]));
    Badframe=unique([Blue_on_frame nonvalid_frame]);
    Goodframe=setdiff([1:nTime],Badframe);

    roisD_order_ind=cellfun(@find,roisD_order{f},'UniformOutput',false);

    nexttile([1 1]);
    PFinframe=find(ind2vec(nTime, Goodframe,1)>0 & PFvec{f}); PFoutframe=find(ind2vec(nTime, Goodframe,1)>0 & (1-PFvec{f}));
    Runframe=find(ind2vec(nTime, Goodframe,1)>0 & runVec{f}); Restframe=find(ind2vec(nTime, Goodframe,1)>0 & (1-runVec{f}));
    silentGoodframe=find(ind2vec(nTime, Goodframe,1)>0 & silenttime_vec);
    corrMat{f,1}=get_corrMat(sub_ch{1},sub_ch{2},Goodframe);
    corrMat{f,2}=get_corrMat(sub_ch{1},sub_ch{2},silentGoodframe); %silent
    corrMatPF{f,1}=get_corrMat(sub_ch{1},sub_ch{2},PFinframe);
    corrMatPF{f,2}=get_corrMat(sub_ch{1},sub_ch{2},PFoutframe);
    corrMatRun{f,1}=get_corrMat(sub_ch{1},sub_ch{2},Runframe);
    corrMatRun{f,2}=get_corrMat(sub_ch{1},sub_ch{2},Restframe);
    
    for silenceind=1:2
    corrMat2=corrMat{f,silenceind};
    corrMat2(eye(size(corrMat2,1)) == 1) = NaN;
    Corrcoeff{f,silenceind}=NaN(5,5);

    for dClassI=1:5 %Basal, Soma, Trunk, Oblique, Distal
        for dClassJ=dClassI:5
            Corrcoeff{f,silenceind}(dClassI,dClassJ)=mean(corrMat2(labelvec{f}==dClassI,labelvec{f}==dClassJ),[1 2],'omitnan');
        end
    end
    end

    cmap_label=hsv(5);
    cmap_label=cmap_label([4 2 3 5 1],:);
    distMatrix = sqrt(2*(1 - corrMat{f,1}));
    distMatrix(logical(eye(size(distMatrix)))) = 0;
    Z_ref = linkage(squareform(distMatrix), 'average');
    leafOrder = optimalleaforder(Z_ref,squareform(distMatrix));
    Cluster_ref= switchlabel(cluster(Z_ref, 'maxclust', 5));
    [~, order_cluster]=sort(Cluster_ref,'ascend');

    l=imshow_label(corrMat{f,1}(:,:),labelvec{f},cmap_label,{'Basal','Soma','Trunk','Oblique','Apical (distal)'});
    axis equal tight off
    if f==foi(end)
    else
        legend off
    end
    caxis([-0.3 1])
    title(num2str(f))
    drawnow;
end
colormap(turbo(256))
%% Correlation coefficient between regions
figure(119); clf;
label_str={'Basal','Soma','Trunk','Oblique','Distal'}; pair_str=[];
Corr_to_plot=[1 1;5 5;3 3;4 4; 1 5]; g=1; show_c=[];
for f=foi
    for p=1:size(Corr_to_plot,1)
        show_c(g,p)=Corrcoeff{f,1}(Corr_to_plot(p,1),Corr_to_plot(p,2));
        pair_str{p}=[label_str{Corr_to_plot(p,1)} ' & ' label_str{Corr_to_plot(p,2)}];
    end
    g=g+1;
end
p=Boxplot_wPoints2(show_c,distinguishable_colors(size(Corr_to_plot,1)));
%drawPValueLines(p,0,'TextYOffset',0.035,'StepHeight',0.08)
set(gca,'XTick',[1:size(Corr_to_plot,1)],'XTickLabel',pair_str);
xlim([0.5 size(Corr_to_plot,1)+0.5]);
ylabel('Correlation coefficient')
ylim([-0.4 1.1])
%% Correlation over distance
figure(121); clf; figure(221); clf; 
pair_str=[]; t_consts_linear=[]; 
l_const=[]; Rsquare=[]; N_points=[];
l_const1d=[]; Rsquare1d=[];
dClass2plot=[[{1},{1}];[{[3 5]} {[3 5]}];[{[4]} {[4]}]];
corr_str={'B & B','A & A','O & O'};
cmap_corrscatter=[0.5 0.8 1; 1 0.5 0.5]; 
cmap_errorbar=[0 0.5 1; 1 0 0; 1 0.8 0];
cmap_neuron=hsv(max(foi));
M=[]; S=[]; M1d=[]; S1d=[]; dendbin=[0:30:300 450 500 600];
for dc=1:size(dClass2plot,1) %region class
    dax2=[]; cax2=[]; nax2=[];
    for f=foi
        DistMat=interDendDist{f}(Dist_order{f}(noi_dist{f}),Dist_order{f}(noi_dist{f}));
        DistMat=max(cat(3,DistMat,DistMat'),[],3);
        DistMat1d = abs(dendaxis1d{f} - dendaxis1d{f}');
        
        uppertri_ind=triu(ones(size(DistMat)));
        label_ind=ismember(labelvec{f},dClass2plot{dc,1})'*ismember(labelvec{f},dClass2plot{dc,2});
        dax=tovec(DistMat(uppertri_ind & label_ind));
        cax=tovec(corrMat{f,1}(uppertri_ind & label_ind));
        dax1d=tovec(DistMat1d(uppertri_ind & label_ind));
        dax2=[dax2; dax]; cax2=[cax2; cax]; nax2=[nax2; ones(length(cax),1)*f];
        N_points(f,dc)=length(dax);
        if N_points(f,dc)>5
            [M(f,:,dc) S(f,:,dc) corr_xbin indicies]=binning_data({[dax cax]},dendbin);
            [M1d(f,:,dc) S1d(f,:,dc) corr_xbin1d indicies]=binning_data({[dax1d cax]},dendbin);
            [y_fit, params, Rsquare(f,dc)] = fitExpDecay(dax,cax,[0:500]');
            [y_fit1d, params1d, Rsquare1d(f,dc)] = fitExpDecay(dax1d,cax,[0:500]');
            l_const(f,[dc dc+size(dClass2plot,1)])=[params(2) range(dax)];
            l_const1d(f,[dc dc+size(dClass2plot,1)])=[params1d(2) range(dax1d)];
            % figure(221); nexttile([1 1]);
            % plot(dax,cax,'o',[0:500],y_fit,'r')
            % title([num2str(f) ', ' corr_str{dc} ', l_c=' num2str(l_const(f,dc))])
        else
            N_points(f,dc)=NaN; l_const(f,dc)=NaN; l_const1d(f,dc)=NaN;
            M(f,1:(length(dendbin)-1),dc)=NaN; S(f,1:(length(dendbin)-1),dc)=NaN;
            M1d(f,1:(length(dendbin)-1),dc)=NaN; S1d(f,1:(length(dendbin)-1),dc)=NaN;
        end
    end
    % scatter(dax2,cax2,15,cmap_corrscatter(dc,:),'filled'); hold all
    % hold all
    M_total=mean(M(foi,:,dc),1,'omitnan'); S_total=sem(M(foi,:,dc),1); N_total=sum(~isnan(M(foi,:,dc)),1);
    M_total(N_total<4)=NaN; S_total(N_total<4)=NaN;

    M_total1d=mean(M1d(foi,:,dc),1,'omitnan'); S_total1d=sem(M1d(foi,:,dc),1); N_total=sum(~isnan(M1d(foi,:,dc)),1);
    M_total1d(N_total<4)=NaN; S_total1d(N_total<4)=NaN;

    figure(121);
    ax1=[nexttile(1,[1 1])];
    errorbar(corr_xbin,M_total,S_total,'linewidth',2,'color',cmap_errorbar(dc,:)); hold all
    xlabel('Pairwise contour distance (\mum)'); ylabel('Correlation');
    ax2=[nexttile(2,[1 1])];
    errorbar(corr_xbin1d,M_total1d,S_total1d,'linewidth',2,'color',cmap_errorbar(dc,:)); hold all
    xlabel('\Deltax (\mum)'); ylabel('Correlation');
    % nexttile([1 1])
    % scatter_heatmap(dax2,cax2,10); hold all
    
    %scatter(dax2,cax2,10,cmap_neuron(nax2,:),'filled'); hold all
    pair_str{dc}=[label_str{dClass2plot{dc,1}} ' & ' label_str{dClass2plot{dc,1}}];
end
legend([{'Basal & Basal'},{'Apical & Apical'},{'Oblique & Oblique'}])
linkaxes([ax1, ax2]);
xlim([0 400])
% Correlation length constant
figure(122); clf;
nexttile([1 1])
l_const_show=l_const(foi,1:size(dClass2plot,1));
l_const_show(Rsquare(foi,:)<0)=NaN;
M=mean(l_const_show,1,'omitnan'); S=std(l_const_show,0,1,'omitnan');
p=Boxplot_wPoints2(l_const_show,cmap_errorbar);
drawPValueLines(p,0,'TextYOffset',50,'StepHeight',100)
set(gca,'XTick',[1:size(l_const_show,1)],'XTickLabel',{'Basal','Apical','Oblique'});
ylabel('Correlation length (\mum)')
nexttile([1 1])

l_const_show=l_const1d(foi,1:size(dClass2plot,1));
l_const_show(Rsquare1d(foi,1:size(dClass2plot,1))<0 | l_const1d(foi,size(dClass2plot,1)+1:2*size(dClass2plot,1))<100)=NaN;
M=mean(l_const_show,1,'omitnan'); S=std(l_const_show,0,1,'omitnan');
p=Boxplot_wPoints2(l_const_show,cmap_errorbar);
drawPValueLines(p,0,'TextYOffset',50,'StepHeight',100)
set(gca,'XTick',[1:size(l_const_show,1)],'XTickLabel',{'Basal','Apical','Oblique'});
ylabel('Correlation length (\mum)')

%% correlation length during run/rest, PF in/out
pair_str=[]; l_const_Run=[]; l_const_PF=[];
RsquareP=[]; RsquareR=[];
dClass2plot=[[{1},{1}];[{[3 5]} {[3 5]}]; [{[4]} {[4]}]];
dcStr={'Basal','Apical','Oblique'}; MP=[]; SP=[]; MR=[]; SR=[];
cmap_corrscatter=distinguishable_colors(size(dClass2plot,1));
figure(123); clf; %tiledlayout(2,size(dClass2plot,1));
for dc=1:size(dClass2plot,1)    
    for onoff=1:2
    dax2=[]; cax2=[]; cax3=[];
    for f=foi
        DistMat=interDendDist{f}(Dist_order{f}(noi_dist{f}),Dist_order{f}(noi_dist{f}));
        DistMat=max(cat(3,DistMat,DistMat'),[],3);
        DistMat1d = abs(dendaxis1d{f} - dendaxis1d{f}');
        uppertri_ind=triu(ones(size(DistMat1d)));
        label_ind=ismember(labelvec{f},dClass2plot{dc,1})'*ismember(labelvec{f},dClass2plot{dc,2});
        dax=tovec(DistMat1d(uppertri_ind & label_ind));
        cax=tovec(corrMatPF{f,onoff}(uppertri_ind & label_ind));        
        caxR=tovec(corrMatRun{f,onoff}(uppertri_ind & label_ind));        
        dax2=[dax2; dax]; cax2=[cax2; cax]; cax3=[cax3; caxR];
        if length(dax)>10
        [MP(f,:,dc,onoff) SP(f,:,dc,onoff) x_bin_centers indicies]=binning_data({[dax cax]},dendbin);            
        [y_fit2, params, RsquareP(f,dc,onoff)] = fitExpDecay(dax(dax>-1),cax(dax>-1),[0:500]);
        % nexttile([1 1]);
        % plot(dax(dax>-1),cax(dax>-1),'.',[0:500],y_fit2)
        l_const_PF(f,dc,onoff)=params(2);        
        [MR(f,:,dc,onoff) SR(f,:,dc,onoff) x_bin_centers indicies]=binning_data({[dax caxR]},dendbin);            
        [y_fit2, paramsR, RsquareR(f,dc,onoff)] = fitExpDecay(dax(dax>-1),caxR(dax>-1),[0:500]);
        % nexttile([1 1]);
        % plot(dax(dax>-1),caxR(dax>-1),'.',[0:500],y_fit2)
        l_const_Run(f,dc,onoff)=paramsR(2);        
        else
        l_const_Run(f,dc,onoff)=NaN;    
        l_const_PF(f,dc,onoff)=NaN;
        end
    end
    nexttile(dc,[1 1]);
    MP_total=mean(MP(foi,:,dc,onoff),1,'omitnan'); SP_total=sem(MP(foi,:,dc,onoff),1); NP_total=sum(~isnan(MP(foi,:,dc,onoff)),1);
    MP_total(NP_total<4)=NaN; SP_total(NP_total<4)=NaN;
    errorbar(x_bin_centers,MP_total,SP_total,'linewidth',2,'color',cmap_errorbar(dc,:)/onoff); hold all
    %xlabel('Pairwise distance (\mum)'); ylabel('Correlation');
    xlabel('\Deltax (\mum)'); ylabel('Correlation');
    title([dcStr{dc} ' PF on/off']); legend({'In','Out'});

    nexttile(dc+size(dClass2plot,1),[1 1]);
    MR_total=mean(MR(foi,:,dc,onoff),1,'omitnan'); SR_total=sem(MR(foi,:,dc,onoff),1); NR_total=sum(~isnan(MR(foi,:,dc,onoff)),1);
    MR_total(NR_total<4)=NaN; SR_total(NR_total<4)=NaN;
    errorbar(x_bin_centers,MR_total,SR_total,'linewidth',2,'color',cmap_errorbar(dc,:)/onoff); hold all
    xlabel('\Deltax (\mum)'); ylabel('Correlation'); legend({'Run','Rest'});
    title([dcStr{dc} ' run/rest'])
    end
end

figure(124); clf; tiledlayout(2,3);
for dc=1:3
    nexttile([1 1])
l_const_show=permute(l_const_PF(foi,dc,:),[1 3 2]);
p=Boxplot_wPoints2(l_const_show,cmap_errorbar(dc,:));
drawPValueLines(p,0,'TextYOffset',50,'StepHeight',100)
set(gca,'XTick',[1:2],'XTickLabel',{'PF in','PF out'});
ylabel('Correlation length (\mum)'); xlim([0.5 2.5]);
    nexttile(3+dc,[1 1])
l_const_show=permute(l_const_Run(foi,dc,:),[1 3 2]);
p=Boxplot_wPoints2(l_const_show,cmap_errorbar(dc,:));
drawPValueLines(p,0,'TextYOffset',50,'StepHeight',100)
set(gca,'XTick',[1:2],'XTickLabel',{'Run','Rest'});
ylabel('Correlation length (\mum)'); xlim([0.5 2.5]);
end

f=20;
figure(125);
clf;
for j=[2 7 20 33];%1:size(Subthreshold{f},1)
    nexttile([1 1]);
    Ft=Ftprnts{f}(:,:,Dist_order{f}(noi_dist{f}));
    bd=cell2mat(bwboundaries(imdilate(Ft(:,:,j)>0,strel('diamond',5))));
    showScaleImageSmooth(permute(Ft>0,[1 2 3]), corrMat{f,2}(j,:),colormap(gen_colormap([0 0.5 1; 1 1 1; 1 0 0])),[-0.25 1]); hold all
    plot(bd(:,2),bd(:,1),'color',[1 1 1])
end
cb=colorbar;
colormap(gen_colormap([0 0.5 1; 1 1 1; 1 0 0]));
cb.Ticks=[0 1];
cb.TickLabels=[-0.4 1];
%% SWC image (Figure 2a)
figure(100); clf;
f=20;
cmap_ExTr=gen_colormap(Plasma,10);
SWCpoints=SWC{f};
SWCpoints(:,3)=SWCpoints(:,3)+5;
SWCpoints(1,3)=50;
SomaROI=[13 14 15]; ApicalROI=[28 29 30]; BasalROI=[8 9 10];
nROI=size(NormalizedTrace_dirt{f},1);
ftprint=Ftprnts{f}(:,:,Dist_order{f}(noi_dist{f}));
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

%% Peak triggered average (Figure 2)
PeakSTA=[]; peakvec=[]; TroughSTA=[]; troughvec=[]; PeakTrough_FR=[]; 
nTauPeak=[500 500];
PeakTroughMat=[]; PeakSTASpClass=[]; TroughSTASpClass=[]; 
for f=foi
    StopFreq=[12 200]; %filter high frequency
    [nROI nTime]=size(Subthreshold{f});
    Subthreshold_int=interpolateNaN(Subthreshold{f});
    FilterTr=[]; 
    for n=1:nROI
        %[PhaseTrace BasalSubFilt BasalthetaPower] = get_phase(Subthreshold_int(n,:), 1000, FilterFreq);
        perispikefrm=unique([tovec(find(allSpikeMat{f}(1,:))'+[-10:40]); find(CStrace{f})']);
        perispikefrm(perispikefrm<1 | perispikefrm>nTime)=[];
        FilterTr(n,:)=get_bandstop(Subthreshold_int(n,:),1000,StopFreq);
        FilterTr(n,:)=FilterTr(n,:)-movmedian(FilterTr(n,:),300,2);
        [pks, locs] = findpeaks(FilterTr(n,:),'MinPeakHeight', 0.35,'MinPeakDistance',50, ...
            'MinPeakProminence', 0.15);  % Prominence can also be tuned
        [trough, loc_trgh] = findpeaks(-FilterTr(n,:),'MinPeakHeight', 0.4,'MinPeakDistance',50, ...
            'MinPeakProminence', 0.15);  % Prominence can also be tuned
        % Turn to 0 during blue Stim and peri-spike frame
        peakvec{f}(n,:)=ind2vec(nTime,locs,1); 
        peakvec{f}(n,ind2vec(nTime,perispikefrm,1)>0)=0;
        troughvec{f}(n,:)=ind2vec(nTime,loc_trgh,1); 
        troughvec{f}(n,ind2vec(nTime,perispikefrm,1)>0)=0;

        [PeakSTA{f}(:,:,n), peakMat, peakloc]=get_STA(Subthreshold{f},peakvec{f}(n,:),nTauPeak(1),nTauPeak(2));
        [TroughSTA{f}(:,:,n), troughMat, troughloc]=get_STA(Subthreshold{f},troughvec{f}(n,:),nTauPeak(1),nTauPeak(2));
        [~, PeakSTASpClass{f,n}]=get_STA(allSpikeClassVecMat{f},peakvec{f}(n,:),nTauPeak(1),nTauPeak(2));
        [~, TroughSTASpClass{f,n}]=get_STA(allSpikeClassVecMat{f},troughvec{f}(n,:),nTauPeak(1),nTauPeak(2));

        % peakvec{f}(n,:)=ind2vec(nTime,peakloc(sum(isnan(peakMat{f}),[1 3])==0),1);
        % troughvec{f}(n,:)=ind2vec(nTime,troughloc(sum(isnan(troughMat),[1 3])==0),1);
        [~, post_spikeFrame, ~, postPeakdist] = find_nearestSpike(peakloc, max(allSpikeClassMat{f}(1:3,:),[],1)>0);
        nanfrm=find(isnan(post_spikeFrame)); post_spikeFrame(nanfrm)=1;
        [postPeakSpclass, ~]=find(allSpikeClassVecMat{f}(1:3,post_spikeFrame));
        postPeakSpclass(nanfrm)=NaN;

        [~, post_spikeFrame, ~, postTroughdist] = find_nearestSpike(troughloc, max(allSpikeClassMat{f}(1:3,:),[],1)>0);
        nanfrm=find(isnan(post_spikeFrame)); post_spikeFrame(nanfrm)=1;
        [postTroughSpclass, ~]=find(allSpikeClassVecMat{f}(1:3,post_spikeFrame));
        postTroughSpclass(nanfrm)=NaN;

        [~, ~, ~, postCSPeakdist] = find_nearestSpike(peakloc, max(allSpikeClassMat{f}(2,:),[],1)>0);
        [~, ~, ~, postCSTroughdist] = find_nearestSpike(troughloc, max(allSpikeClassMat{f}(2,:),[],1)>0);

        [~, AUC_peak, ~, Amp_peak] = get_AUC(permute(peakMat,[1 3 2]), repmat(nTauPeak(1)+1,nROI,1,length(peakloc)), 10, 10);
        [~, AUC_trough, ~, Amp_trough] = get_AUC(permute(troughMat,[1 3 2]), repmat(nTauPeak(1)+1,nROI,1,length(troughloc)), 10, 10);
        mat2cat_peak=[ones(length(peakloc),1) [1:length(peakloc)]' peakloc' repmat([f n labelvec{f}(n) dendaxis{f}(n)],length(peakloc),1)...
                 Amp_peak(n,:)' AUC_peak(n,:)' PFvec{f}(peakloc)' BlueStim{f}(peakloc)' VRtrack{f}(5,peakloc)' postPeakSpclass postPeakdist' postCSPeakdist'...
                 mean(AUC_peak(labelvec{f}==2,:),1,'omitnan')' mean(AUC_peak(labelvec{f}==1,:),1,'omitnan')' mean(AUC_peak(labelvec{f}==3,:),1,'omitnan')' mean(AUC_peak(labelvec{f}==5,:),1,'omitnan')'];

        mat2cat_trough=[zeros(length(troughloc),1) [1:length(troughloc)]' troughloc' repmat([f n labelvec{f}(n) dendaxis{f}(n)],length(troughloc),1)...
                 Amp_trough(n,:)' AUC_trough(n,:)' PFvec{f}(troughloc)' BlueStim{f}(troughloc)' VRtrack{f}(5,troughloc)' postTroughSpclass postTroughdist' postCSTroughdist'...
                 mean(AUC_trough(labelvec{f}==2,:),1,'omitnan')' mean(AUC_trough(labelvec{f}==1,:),1,'omitnan')' mean(AUC_trough(labelvec{f}==3,:),1,'omitnan')' mean(AUC_trough(labelvec{f}==5,:),1,'omitnan')'];
        PeakTroughMat=[PeakTroughMat; mat2cat_peak; mat2cat_trough];
    end
    disp(['Peak triggered average finished, file:' num2str(f)])
    % Peak/trough, Peak/trough ID, frame, Neuron, ROI, dendrite class, distance from soma, amplitude, AUC ,PF in, Blue, VR pos, AUC soma, AUC basal, AUC trunk, AUC distal
    % tempFR=PlaceTrigger_average(double([peakvec{f}; troughvec{f}]),100,VRtrack{f},0.002,115); %total trace
    % PeakTrough_FR{f,1}=tempFR(:,:,1:nROI);
    % PeakTrough_FR{f,2}=tempFR(:,:,nROI+1:2*nROI);
end
fieldName={'IsPeak', 'PeakID', 'Frame', 'FileInd', 'ROI', 'dClass', 'Distance', 'Amp', 'AUC' , 'InPF', 'IsBlue', 'VRpos','NearestSpClass','NearestSpdist','NearestCSdist' ...
           'AUC_soma', 'AUC_basal', 'AUC_trunk', 'AUC_distal'};
PeakTroughMat=array2table(PeakTroughMat,'VariableNames',fieldName);
%% Temporal distance between peak and CS/SS
Spdist_bin=10.^[0:0.3:5];
Spdist_c=mean([Spdist_bin(2:end); Spdist_bin(1:end-1)],1);
cmap=hsv(3);
figure(141); clf; 
for dc=[1 5]
    nexttile([1 1]);
showInd=(PeakTroughMat.dClass==dc & PeakTroughMat.IsPeak==1); 
SpDistCount=histcounts(PeakTroughMat.NearestSpdist(showInd),Spdist_bin);
CSDistCount=histcounts(PeakTroughMat.NearestCSdist(showInd),Spdist_bin,'Normalization','probability');
CSRatioDist=CSDistCount./SpDistCount;
plot(Spdist_c,CSDistCount)
hold all;
showInd=(PeakTroughMat.dClass==dc & PeakTroughMat.IsPeak==0);
SpDistCount=histcounts(PeakTroughMat.NearestSpdist(showInd),Spdist_bin);
CSDistCount=histcounts(PeakTroughMat.NearestCSdist(showInd),Spdist_bin,'Normalization','probability');
CSRatioDist=CSDistCount./SpDistCount;
plot(Spdist_c,CSDistCount)
set(gca,'xscale','log'); %ylim([0 1]);
end
%% Temporal distance between peaks and troughs
figure(142); clf;
cmap=hsv(4); noi=32; t_edge=[0:30:1000]; %8 9 /13 /21 22/ 28 29 32
nexttile([1 1]);
PeakMatsub=permute(PeakSTASpClass{f,noi},[2 3 1])>0;
[E, T ,Cl] = ind2sub(size(PeakMatsub), find(PeakMatsub));
scatter(T-500,E,20,cmap(Cl,:),'filled'); hold all
for dc=[1:4]
[y x]=histcounts(T(Cl==dc),t_edge);
plot(mean([x(2:end); x(1:end-1)])-500,y,'color',cmap(dc,:));
end

nexttile([1 1]);
TroughMatsub=permute(TroughSTASpClass{f,noi},[2 3 1])>0;
[E, T ,Cl] = ind2sub(size(TroughMatsub), find(TroughMatsub));
scatter(T-500,E,20,cmap(Cl,:),'filled'); hold all
for dc=[1:4]
[y x]=histcounts(T(Cl==dc),t_edge);
plot(mean([x(2:end); x(1:end-1)])-500,y,'color',cmap(dc,:));
end

%% plot Peak triggered average
figure(132); clf; tiledlayout(1,2);
showf=[23]; cax=[-0.2 0.6];
[nROI, nTime]=size(Subthreshold{showf});
[~,showInd(1)]=min(dendaxis{showf}); [showInd(2)]=find(dendaxis{showf}==0); [~,showInd(3)]=max(dendaxis{showf});
for i=[1 3]
nexttile([1 1]);
showSTA=PeakSTA{showf}-mean(PeakSTA{showf}(:,1:10,:),2);
imagesc([-nTauPeak(1):nTauPeak(2)],[1:nROI],showSTA(:,:,showInd(i)),cax);
set(gca,'YTick',[1 find(Dist_order{showf}(noi_dist{showf})==1) sum(noi_dist{showf})],'YTickLabel',num2str([min(dendaxis{showf}) 0 max((dendaxis{showf}))]',3))
xlabel('Time (ms)')
ylabel('Distance from soma (\mum)')
% nexttile([1 1]);
% showSTA=TroughSTA{showf}-mean(TroughSTA{showf}(:,1:10,:),2);
% imagesc([-nTauPeak(1):nTauPeak(2)],[1:nROI],showSTA(:,:,showInd(i)),cax);
% set(gca,'YTick',[1 find(Dist_order{showf}(noi_dist{showf})==1) sum(noi_dist{showf})],'YTickLabel',num2str([min(dendaxis{showf}) 0 max((dendaxis{showf}))]',3))
% xlabel('Time (ms)')
end
colormap(turbo);

%% Peak correlation over distance (figure 2)
figure(131); clf; tiledlayout(1,5);
dClass2plot=[[{1},{1}];[{[3 5]} {[3 5]}];[{[4]} {[4]}]];
dClass2temporal=[[{1},{1}];[{[3 5]} {[3 5]}];[{[1]} {[5]}];[{[5]} {[1]}]]; %Read & Triggering
pair_str=[]; Length_PeakMat=[]; %neuronInd, Trigger Pair i, Trigger Label, Measure Label, Length C, R^2, Legnth C dX, R^2 dX
cmap_corrscatter=[0.5 0.8 1; 1 0.5 0.5]; 
cmap_errorbar=[0.5 0 1; 0.3 0.3 1; 1 0.4 0.1; 1 0 0];
cmap_neuron=hsv(max(foi));
M=[]; S=[]; M1d=[]; S1d=[]; dendbin=[0 15:45:450 500 600];
PeakTraceROI=[];
for dc=1:size(dClass2plot,1) %region class to plot
    dax2=[]; cax2=[]; nax2=[];
    for f=foi
        DistMat=interDendDist{f}(Dist_order{f}(noi_dist{f}),Dist_order{f}(noi_dist{f}));
        DistMat=max(cat(3,DistMat,DistMat'),[],3);
        DistMat1d = abs(dendaxis1d{f} - dendaxis1d{f}');
        PeakZeroMat=squeeze(PeakSTA{f}(:,nTauPeak(1)+1,:));
        uppertri_ind=triu(ones(size(DistMat)));
        label_ind=ismember(labelvec{f},dClass2plot{dc,1})'*ismember(labelvec{f},dClass2plot{dc,2});
        dax=tovec(DistMat(uppertri_ind & label_ind));
        peakdat2fit=tovec(PeakZeroMat(uppertri_ind & label_ind));
        peakdat2fit=peakdat2fit./mean(peakdat2fit(dax==0));
        dax1d=tovec(DistMat1d(uppertri_ind & label_ind));
        dax2=[dax2; dax]; cax2=[cax2; peakdat2fit]; nax2=[nax2; ones(length(peakdat2fit),1)*f];
        N_points=sum(~isnan(peakdat2fit));
        if N_points>5 % number of fitting point is more than 5
            [M(f,:,dc) S(f,:,dc) dendbin_c]=binning_data({[dax peakdat2fit]},dendbin);
            [M1d(f,:,dc) S1d(f,:,dc) dendbin1d_c]=binning_data({[dax1d peakdat2fit]},dendbin);
            [y_fit, params, Rsquare] = fitExpDecay(dax,peakdat2fit,[0:500]');
            [y_fit1d, params1d, Rsquare1d] = fitExpDecay(dax1d,peakdat2fit,[0:500]');
            l_const=[params(2)];  l_const1d=[params1d(2)];
            MaxDistance=[range(dax) range(dax1d)];
            %PeakPercentile=sum(mean(peakdat2fit(dax<10))>peakdat2fit)/length(peakdat2fit);
            peakdatzscore=zscore(peakdat2fit);
            % figure(221); nexttile([1 1]);
            % plot(dax,cax,'o',[0:500],y_fit,'r')
            % title([num2str(f) ', ' corr_str{dc} ', l_c=' num2str(l_const(f,dc))])
            Length_PeakMat=[Length_PeakMat; [f, dc, l_const, Rsquare, l_const1d, Rsquare1d, N_points, MaxDistance, mean(peakdatzscore(dax<30))]]; 
        else
            l_const=NaN; l_const1d=NaN; MaxDistance=NaN(1,2); Rsquare=NaN; Rsquare1d=NaN; peakdatzscore=NaN;
            M(f,1:(length(dendbin)-1),dc)=NaN; S(f,1:(length(dendbin)-1),dc)=NaN;
            M1d(f,1:(length(dendbin)-1),dc)=NaN; S1d(f,1:(length(dendbin)-1),dc)=NaN;
        end
        %NeuronInd, Pair class, Length Cst, R^2, Legnth Cst dX, R^2 dX, N fit points, range of distance, range dx, peak zscore
    end
    M_total=mean(M(foi,:,dc),1,'omitnan'); S_total=sem(M(foi,:,dc),1); N_total=sum(~isnan(M(foi,:,dc)),1);
    M_total(N_total<4)=NaN; S_total(N_total<4)=NaN;

    M_total1d=mean(M1d(foi,:,dc),1,'omitnan'); S_total1d=sem(M1d(foi,:,dc),1); N_total=sum(~isnan(M1d(foi,:,dc)),1);
    M_total1d(N_total<4)=NaN; S_total1d(N_total<4)=NaN;

    ax1=[nexttile(1,[1 3])];
    l(dc)=errorbar_shade(dendbin_c,M_total,S_total,cmap_errorbar(dc,:)); hold all
    xlabel('Distance from peak ROI (\mum)'); ylabel('Normalized \DeltaF/F');
    % ax2=[nexttile(2,[1 1])];
    % errorbar(dendbin1d_c,M_total1d,S_total1d,'linewidth',2,'color',cmap_errorbar(dc,:)); hold all
    % xlabel('\Deltax (\mum)'); ylabel('Normalized \DeltaF/F');
    % pair_str{dc}=[label_str{dClass2plot{dc,1}} ' & ' label_str{dClass2plot{dc,1}}];
end
%linkaxes([ax1, ax2]);
xlim([0 400])

% Plot local stimulation data
saveto='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism';
load([saveto '/LocalStimulationResult.mat']);
path2show=find(Lstaall(:,2)>0.5 & Lstaall(:,3)>0.4); % Percentile, R^2
data2show=STADistmat(path2show); data2show1d=STADistmat(path2show);
[deltax, dsort]=cellfun(@(x) sort(x(:,1)),data2show,'UniformOutput',false);
[~, dsort1d]=cellfun(@(x) sort(x(:,3)),data2show1d,'UniformOutput',false);
data2show=cellfun(@(x,y) [x(y,1) x(y,2)/x(x(:,1)==0,2)],data2show,dsort,'UniformOutput',false); %normalize
data2show1d=cellfun(@(x,y) [(x(y,3)) x(y,2)/x(x(:,1)==0,2)],data2show1d,dsort1d,'UniformOutput',false);
%cellfun(@(x) plot(x(:,1),x(:,2),'color',[0.7 0.7 0.7]),data2show)
[M S x_bin N]=binning_data(data2show,dendbin);
l(4)=errorbar_shade(x_bin,M,S./sqrt(sum(cellfun(@sum,N),2))',cmap_errorbar(4,:));
xlim([0 200]); box off;
legend(l,[{'Basal & Basal'},{'Apical & Apical'},{'Oblique & Oblique'},{'Local Stim.'}])

nexttile([1 2])
cmap=[0 0.6 1;0 0 0; 0.7 0 0.9];
SponLength2show=find(Length_PeakMat(:,4)>0.2 & Length_PeakMat(:,10)>1);
LengthC2show={Lstaall(path2show,1), Length_PeakMat(SponLength2show,3)};
p=Boxplot_wPoints2(LengthC2show,cmap([1 2],:)); axis tight; hold all
drawPValueLines(p,300,'TextYOffset',0.1,'StepHeight',0.1)
set(gca,'yscale','log','xtick',[1 2],'xticklabel',{'Opto. Stim.','Spontaneous'},'ytick',[10 100 1000])
ylabel('Decay length (\mum)'); ylim([10 1700]); xlim([0.5 2.5]); box off;

figure(134); clf; tiledlayout(2,2); cmap_temp=gen_colormap(Plasma,4);
title_str={'Basal peak, Basal voltage','Apical peak, Apical voltage','Apical peak, Basal voltage','Basal peak, Apical voltage'};
for dc=[1 4 3 2]
    for f=foi
    ROI2plot=ismember(labelvec{f},dClass2temporal{dc,2})'*ismember(labelvec{f},dClass2temporal{dc,1});
    % uppertri_ind=1-triu(ones(length(labelvec{f})));
    % ROI2plot=ROI2plot & uppertri_ind;
    STAsub=permute(PeakSTA{f},[3 1 2]);
    PeakTraceROI{f}(:,dc)=tovec(STAsub)'*tovec(ROI2plot)/sum(ROI2plot(:));
    end

    nexttile([1 1]);
    PeakTrace2show=cell2mat(cellfun(@(x) x(:,dc),PeakTraceROI(foi),'UniformOutput',false));
    plot([-nTauPeak(1):nTauPeak(2)],PeakTrace2show,'color',[0.7 0.7 0.7]); hold all;
    errorbar_shade([-nTauPeak(1):nTauPeak(2)],mean(PeakTrace2show,2,'omitnan'),std(PeakTrace2show,0,2,'omitnan'),cmap_temp(round(dc/4-0.1)*3+1,:));
    xlabel('Time (ms)'); ylabel('Subthreshold (spike height)');
    title(title_str{dc})
    box off;
end

figure(135); clf; tiledlayout(1,2); cmap_temp=gen_colormap(Plasma,4);
ddd={[1 4],[3 2]};
title_str={'Basal peak-TA','Apical peak-TA'};
for dc2=1:length(ddd)
    ldc=[]; g=1;
for dc=ddd{dc2}
    nexttile(dc2,[1 1]);
    PeakTrace2show=cell2mat(cellfun(@(x) x(:,dc),PeakTraceROI(foi),'UniformOutput',false));
    %plot([-nTauPeak(1):nTauPeak(2)],PeakTrace2show,'color',[0.7 0.7 0.7]); hold all;
    ldc(g)=errorbar_shade([-nTauPeak(1):nTauPeak(2)],mean(PeakTrace2show,2,'omitnan'),std(PeakTrace2show,0,2,'omitnan'),cmap_errorbar(dc,:));
    xlabel('Time (ms)'); ylabel('Subthreshold (spike height)');
    title(title_str{dc2})
    g=g+1;
end
legend(ldc,{'Basal','Apical'})
end

%% Correlation length constant
figure(122); clf;
nexttile([1 1])
l_const_show=l_const(foi,1:size(dClass2plot,1));
l_const_show(Rsquare(foi,:)<0)=NaN;
M=mean(l_const_show,1,'omitnan'); S=std(l_const_show,0,1,'omitnan');
p=Boxplot_wPoints2(l_const_show,cmap_errorbar);
drawPValueLines(p,0,'TextYOffset',50,'StepHeight',100)
set(gca,'XTick',[1:size(l_const_show,1)],'XTickLabel',{'Basal','Apical','Oblique'});
ylabel('Correlation length (\mum)')
nexttile([1 1])




DistMat=interDendDist{f}(Dist_order{f}(noi_dist{f}),Dist_order{f}(noi_dist{f}));
        DistMat=max(cat(3,DistMat,DistMat'),[],3);
        DistMat1d = abs(dendaxis1d{f} - dendaxis1d{f}');
        uppertri_ind=triu(ones(size(DistMat1d)));
        label_ind=ismember(labelvec{f},dClass2plot{dc,1})'*ismember(labelvec{f},dClass2plot{dc,2});
        dax=tovec(DistMat1d(uppertri_ind & label_ind));
        cax=tovec(corrMatPF{f,onoff}(uppertri_ind & label_ind));        
        caxR=tovec(corrMatRun{f,onoff}(uppertri_ind & label_ind));    



figure(125); clf;
nexttile([1 1])
l_const_show=l_const(foi,1:size(dClass2plot,1));
l_const_show(Rsquare(foi,:)<0)=NaN;
M=mean(l_const_show,1,'omitnan'); S=std(l_const_show,0,1,'omitnan');
p=Boxplot_wPoints2(l_const_show,cmap_errorbar);
drawPValueLines(p,0,'TextYOffset',50,'StepHeight',100)
set(gca,'XTick',[1:size(l_const_show,1)],'XTickLabel',{'Basal','Apical','Oblique'});
ylabel('Correlation length (\mum)')
nexttile([1 1])

l_const_show=l_const1d(foi,1:size(dClass2plot,1));
l_const_show(Rsquare1d(foi,1:size(dClass2plot,1))<0 | l_const1d(foi,size(dClass2plot,1)+1:2*size(dClass2plot,1))<100)=NaN;
M=mean(l_const_show,1,'omitnan'); S=std(l_const_show,0,1,'omitnan');
p=Boxplot_wPoints2(l_const_show,cmap_errorbar);
drawPValueLines(p,0,'TextYOffset',50,'StepHeight',100)
set(gca,'XTick',[1:size(l_const_show,1)],'XTickLabel',{'Basal','Apical','Oblique'});
ylabel('Correlation length (\mum)')


figure(121); clf; figure(221); clf; 
pair_str=[]; t_consts_linear=[]; 
l_const=[]; Rsquare=[]; N_points=[];
l_const1d=[]; Rsquare1d=[];
dClass2plot=[[{1},{1}];[{[3 5]} {[3 5]}];[{[4]} {[4]}]];
corr_str={'B & B','A & A','O & O'};
cmap_corrscatter=[0.5 0.8 1; 1 0.5 0.5]; 
cmap_errorbar=[0 0.5 1; 1 0 0; 1 0.8 0];
cmap_neuron=hsv(max(foi));
M=[]; S=[]; M1d=[]; S1d=[]; dendbin=[0:30:300 450 500 600];
for dc=1:size(dClass2plot,1) %region class
    dax2=[]; cax2=[]; nax2=[];
    for f=foi
        DistMat=interDendDist{f}(Dist_order{f}(noi_dist{f}),Dist_order{f}(noi_dist{f}));
        DistMat=max(cat(3,DistMat,DistMat'),[],3);
        DistMat1d = abs(dendaxis1d{f} - dendaxis1d{f}');
        
        uppertri_ind=triu(ones(size(DistMat)));
        label_ind=ismember(labelvec{f},dClass2plot{dc,1})'*ismember(labelvec{f},dClass2plot{dc,2});
        dax=tovec(DistMat(uppertri_ind & label_ind));
        cax=tovec(corrMat{f,1}(uppertri_ind & label_ind));
        dax1d=tovec(DistMat1d(uppertri_ind & label_ind));
        dax2=[dax2; dax]; cax2=[cax2; cax]; nax2=[nax2; ones(length(cax),1)*f];
        N_points(f,dc)=length(dax);
        if N_points(f,dc)>5
            [M(f,:,dc) S(f,:,dc) corr_xbin indicies]=binning_data({[dax cax]},dendbin);
            [M1d(f,:,dc) S1d(f,:,dc) corr_xbin1d indicies]=binning_data({[dax1d cax]},dendbin);
            [y_fit, params, Rsquare(f,dc)] = fitExpDecay(dax,cax,[0:500]');
            [y_fit1d, params1d, Rsquare1d(f,dc)] = fitExpDecay(dax1d,cax,[0:500]');
            l_const(f,[dc dc+size(dClass2plot,1)])=[params(2) range(dax)];
            l_const1d(f,[dc dc+size(dClass2plot,1)])=[params1d(2) range(dax1d)];
            % figure(221); nexttile([1 1]);
            % plot(dax,cax,'o',[0:500],y_fit,'r')
            % title([num2str(f) ', ' corr_str{dc} ', l_c=' num2str(l_const(f,dc))])
        else
            N_points(f,dc)=NaN; l_const(f,dc)=NaN; l_const1d(f,dc)=NaN;
            M(f,1:(length(dendbin)-1),dc)=NaN; S(f,1:(length(dendbin)-1),dc)=NaN;
            M1d(f,1:(length(dendbin)-1),dc)=NaN; S1d(f,1:(length(dendbin)-1),dc)=NaN;
        end
    end
    % scatter(dax2,cax2,15,cmap_corrscatter(dc,:),'filled'); hold all
    % hold all
    M_total=mean(M(foi,:,dc),1,'omitnan'); S_total=sem(M(foi,:,dc),1); N_total=sum(~isnan(M(foi,:,dc)),1);
    M_total(N_total<4)=NaN; S_total(N_total<4)=NaN;

    M_total1d=mean(M1d(foi,:,dc),1,'omitnan'); S_total1d=sem(M1d(foi,:,dc),1); N_total=sum(~isnan(M1d(foi,:,dc)),1);
    M_total1d(N_total<4)=NaN; S_total1d(N_total<4)=NaN;

    figure(121);
    ax1=[nexttile(1,[1 1])];
    errorbar(corr_xbin,M_total,S_total,'linewidth',2,'color',cmap_errorbar(dc,:)); hold all
    xlabel('Pairwise contour distance (\mum)'); ylabel('Correlation');
    ax2=[nexttile(2,[1 1])];
    errorbar(corr_xbin1d,M_total1d,S_total1d,'linewidth',2,'color',cmap_errorbar(dc,:)); hold all
    xlabel('\Deltax (\mum)'); ylabel('Correlation');
    % nexttile([1 1])
    % scatter_heatmap(dax2,cax2,10); hold all
    
    %scatter(dax2,cax2,10,cmap_neuron(nax2,:),'filled'); hold all
    pair_str{dc}=[label_str{dClass2plot{dc,1}} ' & ' label_str{dClass2plot{dc,1}}];
end
legend([{'Basal & Basal'},{'Apical & Apical'},{'Oblique & Oblique'}])
linkaxes([ax1, ax2]);
xlim([0 400])
% Correlation length constant
figure(122); clf;
nexttile([1 1])
l_const_show=l_const(foi,1:size(dClass2plot,1));
l_const_show(Rsquare(foi,:)<0)=NaN;
M=mean(l_const_show,1,'omitnan'); S=std(l_const_show,0,1,'omitnan');
p=Boxplot_wPoints2(l_const_show,cmap_errorbar);
drawPValueLines(p,0,'TextYOffset',50,'StepHeight',100)
set(gca,'XTick',[1:size(l_const_show,1)],'XTickLabel',{'Basal','Apical','Oblique'});
ylabel('Correlation length (\mum)')
nexttile([1 1])

l_const_show=l_const1d(foi,1:size(dClass2plot,1));
l_const_show(Rsquare1d(foi,1:size(dClass2plot,1))<0 | l_const1d(foi,size(dClass2plot,1)+1:2*size(dClass2plot,1))<100)=NaN;
M=mean(l_const_show,1,'omitnan'); S=std(l_const_show,0,1,'omitnan');
p=Boxplot_wPoints2(l_const_show,cmap_errorbar);
drawPValueLines(p,0,'TextYOffset',50,'StepHeight',100)
set(gca,'XTick',[1:size(l_const_show,1)],'XTickLabel',{'Basal','Apical','Oblique'});
ylabel('Correlation length (\mum)')

%% Basal/Apical/Soma subthreshold dynamics

run_threshold=0.002; max_xcorrT=500;
PSD_sub=[];
preCS_SubBA=[]; preSS_SubBA=[]; ThetaPow=[]; PhaseTr=[];
STmat_SubthOrth=[]; STmat_SubthOrth2=[]; %Orth: A+B, B-A, S; Orth2: B, A, S;
OrthAxis=[1 0 1 0 1; 1 0 -1 0 -1; 0 1 0 0 0]; %basal, soma, trunk, oblique, distal
OrthAxis2=[1 0 0 0 0; 0 0 1 0 1; 0 0 0 0 1]; %basal, soma, trunk, oblique, distal
xcorr_silentMat=NaN(max(foi),max_xcorrT*2+1); BA=NaN(max(foi),max_xcorrT*2+1);
xcorrMat=NaN(max(foi),max_xcorrT*2+1);
VecPower=NaN(max(foi),size(OrthAxis,1)); VecPowerPF=NaN(max(foi),size(OrthAxis,1),2); VecPowerRun=NaN(max(foi),size(OrthAxis,1),2);
VecthetaPower=NaN(max(foi),size(OrthAxis,1)); VecthetaPowerPF=NaN(max(foi),size(OrthAxis,1),2); VecthetaPowerRun=NaN(max(foi),size(OrthAxis,1),2);
OrthaxisXcorr_silent=NaN(size(OrthAxis2,1),size(OrthAxis2,1),max_xcorrT*2+1,max(foi));
OrthaxisXcorr=NaN(size(OrthAxis2,1),size(OrthAxis2,1),max_xcorrT*2+1,max(foi));
FilterFreq=[5 12]; psd_norm=[];

for f=foi
    f
    [nROI nTime]=size(Subthreshold{f});
    filteredSubthres_orthAxis=[]; filteredSubthres_orthAxis2=[];

    % set ROI indices
    roisD_order_ind=cellfun(@find,roisD_order{f},'UniformOutput',false);
    OrthAxis_vec=zeros(size(OrthAxis,1),nROI); OrthAxis_vec2=zeros(size(OrthAxis,1),nROI);
    for d=1:5
        dclassinds=cell2mat(roisD_order_ind(d,:)');
        for v=1:size(OrthAxis,1)
            if abs(OrthAxis(v,d))>0 & isempty(dclassinds)
                OrthAxis_vec(v,:)=NaN;
            else
                suborthvec=ind2vec(nROI,dclassinds,OrthAxis(v,d)/length(dclassinds));
                OrthAxis_vec(v,abs(suborthvec)>0)=suborthvec(abs(suborthvec)>0);
            end
        end
        for v=1:size(OrthAxis2,1)
            if abs(OrthAxis2(v,d))>0 & isempty(dclassinds)
                OrthAxis_vec2(v,:)=NaN;
            else
                suborthvec2=ind2vec(nROI,dclassinds,OrthAxis2(v,d)/length(dclassinds));
                OrthAxis_vec2(v,abs(suborthvec2)>0)=suborthvec2(abs(suborthvec2)>0);
            end
        end
    end
    OrthAxis_vec=OrthAxis_vec(:,cell2mat(roisD_order_ind(:)));
    OrthAxis_vec2=OrthAxis_vec2(:,cell2mat(roisD_order_ind(:)));

    % set ROI indices
    [~, D]=get_eigvector(Subthreshold{f}(cell2mat(roisD_order_ind(:)),:));
    sub_ch=[];
    for ch=1:2
        sub_ch{ch}=get_subthreshold(NormalizedTrace_ch{f,ch},max(allSpikeMat{f}(1,:),[],1)>0,7,17);
        sub_ch{ch}=sub_ch{ch}(cell2mat(roisD_order_ind(:)),:);
        sub_ch{ch} = pcafilterTrace(sub_ch{ch},1:find(cumsum(D)/sum(D)>0.95,1));
        sub_ch{ch}(:,BlueStim{f}>0)=NaN;
        filteredSubthres_orthAxis2(:,:,ch)=(sub_ch{ch}'*OrthAxis_vec2')';
    end
    filteredSubthres_orthAxis = pcafilterTrace(Subthreshold{f}(cell2mat(roisD_order_ind(:)),:),1:find(cumsum(D)/sum(D)>0.95,1));
    filteredSubthres_orthAxis(:,BlueStim{f}>0)=NaN;
    filteredSubthres_orthAxis=(filteredSubthres_orthAxis'*OrthAxis_vec')';

    for dclass=1:4
        [~, STmat_SubthOrth{f,dclass}]=get_STA(filteredSubthres_orthAxis, allSpikeClassMat{f}(dclass,:).*(BlueStim{f}==0),30,30); %pre complex spike
    end

    silenttime=setdiff([1:nTime],unique(find(max([allSpikeMat{f}(1,:); allSpikeClassMat{f}; CStrace{f}; BlueStim{f}])>0)'+[-5:30]));
    silenttime_vec=ind2vec(nTime,silenttime,1);
    filteredSubthres_orthAxis_silent=filteredSubthres_orthAxis;
    filteredSubthres_orthAxis_silent(:,~silenttime_vec)=NaN;
    filteredSubthres_orthAxis_silent2=filteredSubthres_orthAxis2;
    filteredSubthres_orthAxis_silent2(:,~silenttime_vec,:)=NaN;

    for oi=1:size(OrthAxis2,1)
        for oj=oi:size(OrthAxis2,1)
           OrthaxisXcorr_silent(oi,oj,:,f)=nanXCorr(filteredSubthres_orthAxis_silent2(oi,:,1),filteredSubthres_orthAxis_silent2(oj,:,2),max_xcorrT,1);
           OrthaxisXcorr(oi,oj,:,f)=nanXCorr(filteredSubthres_orthAxis2(oi,:,1),filteredSubthres_orthAxis2(oj,:,2),max_xcorrT,1);
        end
        filteredSubthres_orthAxis_silent_interpNaN = interpolateNaN(filteredSubthres_orthAxis_silent);
        if any(~isnan(filteredSubthres_orthAxis_silent_interpNaN(oi,:)))
        [PhaseTr{f}(oi,:),~,ThetaPow{f}(oi,:)] = get_phase(filteredSubthres_orthAxis_silent_interpNaN(oi,:), 1000, FilterFreq);
        [frq, pw] = nanPSD(filteredSubthres_orthAxis_silent2(oi,:,1), 1000, 3000);
        psd_norm{f}=[frq pw];
        end
    end
end

% Basal & apical & soma autocorr and crosscorr
figure(127); clf; t_lag=[-max_xcorrT:max_xcorrT]; ax2=[];
subplot_loc=[1 3 2 9 6 5]; g=1; ROI_str={'Basal','Apical','Distal'};
cmap_label=jet(9);
for oi=1:size(OrthAxis2,1)
    for oj=oi:size(OrthAxis2,1)
        nexttile([1 1])
        plot(t_lag,squeeze(OrthaxisXcorr_silent(oi,oj,:,foi)),'color',[0.8 0.8 0.8]); hold all
        errorbar_shade(t_lag,mean(squeeze(OrthaxisXcorr_silent(oi,oj,:,foi))',1,'omitnan'),std(squeeze(OrthaxisXcorr_silent(oi,oj,:,foi))',0,1,'omitnan'),cmap_label(g,:));
        title([ROI_str{oi} ' & ' ROI_str{oj}])
        g=g+1;
        xlabel('\tau (ms)')
ylabel('R(\tau)')
    end
end


figure(125); clf; cmap=lines(3); tiledlayout(2,3);
% Basal & apical during run or rest
nexttile([1 1]);
for r=1:2
    plot([1 2]+(r-1)*2,Run_BA(foi,[1 3]+r-1),'marker','o','color',[0.6 0.6 0.6]); hold all
    M=mean(Run_BA(foi,[1 3]+r-1),1,'omitnan'); S=std(Run_BA(foi,[1 3]+r-1),0,1,'omitnan')./sqrt(length(foi));
    errorbar([1 2]+(r-1)*2,M,S,'color',cmap(2,:));
    [p]=signrank(Run_BA(foi,[1]+r-1),Run_BA(foi,[3]+r-1));
    mm=max(tovec(Run_BA(foi,[1 3]+r-1)));
    line([1 2]+(r-1)*2,[mm*1.1 mm*1.1],'color',[0 0 0]);
    text(1.5+(r-1)*2,mm*1.1+0.01,['p = ' num2str(p,2)],'HorizontalAlignment','center')
end
xlim([0.5 4.5]);
ylabel('Voltage (spike height')
set(gca,'XTick',[1:4],'XTickLabel',{'Basal, Run','Basal, Rest','Apical, Run','Apical, Rest'})

% Basal & apical in PF or out PF
nexttile([1 1]);
for r=1:2
    plot([1 2]+(r-1)*2,PlaceField_BA(foi,[1 3]+r-1),'marker','o','color',[0.6 0.6 0.6]); hold all
    M=mean(PlaceField_BA(foi,[1 3]+r-1),1,'omitnan'); S=std(PlaceField_BA(foi,[1 3]+r-1),0,1,'omitnan')./sqrt(length(foi));
    errorbar([1 2]+(r-1)*2,M,S,'color',cmap(2,:));
    [p]=signrank(PlaceField_BA(foi,[1]+r-1),PlaceField_BA(foi,[3]+r-1));
    mm=max(tovec(PlaceField_BA(foi,[1 3]+r-1)));
    line([1 2]+(r-1)*2,[mm*1.1 mm*1.1],'color',[0 0 0]);
    text(1.5+(r-1)*2,mm*1.1+0.01,['p = ' num2str(p,2)],'HorizontalAlignment','center')
end
xlim([0.5 4.5]);
ylabel('Voltage (spike height')
set(gca,'XTick',[1:4],'XTickLabel',{'Basal, in PF','Basal, out PF','Apical, in PF','Apical, out PF'})

% Basal - apical during run or rest
nexttile([1 1]);
plot([1 2],Run_BAunbal(foi,3:4),'marker','o','color',[0.6 0.6 0.6]); hold all
M=mean(Run_BAunbal(foi,3:4),1,'omitnan'); S=std(Run_BAunbal(foi,3:4),0,1,'omitnan')./sqrt(length(foi));
errorbar([1 2],M,S,'color',cmap(1,:));
[p]=signrank(Run_BAunbal(foi,3),Run_BAunbal(foi,4));
xlim([0.5 2.5])
set(gca,'XTick',[1 2],'XTickLabel',{'Run','Rest'})
title(['p = ' num2str(p,2)])
ylabel('Basal- Apical (Spike height)')

% Basal - Apical in PF or Out PF
nexttile([1 1]);
plot([1 2],PlaceField_BAunbal(foi,3:4),'marker','o','color',[0.6 0.6 0.6]); hold all
M=mean(PlaceField_BAunbal(foi,3:4),1,'omitnan'); S=std(PlaceField_BAunbal(foi,3:4),0,1,'omitnan')./sqrt(length(foi));
errorbar([1 2],M,S,'color',cmap(3,:));
[p]=signrank(PlaceField_BAunbal(foi,3),PlaceField_BAunbal(foi,4));
xlim([0.5 2.5])
ylabel('Basal- Apical (Spike height)')
set(gca,'XTick',[1 2],'XTickLabel',{'In PF','Out PF'})
title(['p = ' num2str(p,2)])

% Theta power during run or rest
nexttile([1 1]);
for r=1:2
    plot([1 2]+(r-1)*2,Run_thetapower(foi,:,r),'marker','o','color',[0.6 0.6 0.6]); hold all
    M=mean(Run_thetapower(foi,:,r),1,'omitnan'); S=std(Run_thetapower(foi,:,r),0,1,'omitnan')./sqrt(length(foi));
    errorbar([1 2]+(r-1)*2,M,S,'color',cmap(2,:));
    [p]=signrank(Run_thetapower(foi,1,r),Run_thetapower(foi,2,r));
    mm=max(tovec(Run_thetapower(foi,:,r)));
    line([1 2]+(r-1)*2,[mm*1.1 mm*1.1],'color',[0 0 0]);
    text(1.5+(r-1)*2,mm*1.1+0.01,['p = ' num2str(p,2)],'HorizontalAlignment','center')
end
xlim([0.5 4.5]);
ylabel('Theta-band power')
set(gca,'XTick',[1:4],'XTickLabel',{'Basal, Run','Basal, Rest','Apical, Run','Apical, Rest'})

% Theta power in PF or out PF
nexttile([1 1]);
for r=1:2
    plot([1 2]+(r-1)*2,PlaceField_thetapower(foi,:,r),'marker','o','color',[0.6 0.6 0.6]); hold all
    M=mean(PlaceField_thetapower(foi,:,r),1,'omitnan'); S=std(PlaceField_thetapower(foi,:,r),0,1,'omitnan')./sqrt(length(foi));
    errorbar([1 2]+(r-1)*2,M,S,'color',cmap(2,:));
    [p]=signrank(PlaceField_thetapower(foi,1,r),PlaceField_thetapower(foi,2,r));
    mm=max(tovec(PlaceField_thetapower(foi,:,r)));
    line([1 2]+(r-1)*2,[mm*1.1 mm*1.1],'color',[0 0 0]);
    text(1.5+(r-1)*2,mm*1.1+0.01,['p = ' num2str(p,2)],'HorizontalAlignment','center')
end
xlim([0.5 4.5]);
ylabel('Theta-band power')
set(gca,'XTick',[1:4],'XTickLabel',{'Basal, in PF','Basal, out PF','Apical, in PF','Apical, out PF'})
    
% nexttile([1 1]);
% plot(t_lag,BAxcorr(foi,:)','color',[0.8 0.8 0.8]); hold all
% errorbar_shade(t_lag,mean(BAxcorr(foi,:),1,'omitnan'),std(BAxcorr(foi,:),0,1,'omitnan'));
% xlabel('\tau (ms)')
% ylabel('R(\tau)')
% nexttile([1 1]);
% plot(t_lag,BBxcorr(foi,:)','color',[0.8 0.8 0.8]); hold all
% errorbar_shade(t_lag,mean(BBxcorr(foi,:),1,'omitnan'),std(BBxcorr(foi,:),0,1,'omitnan'));
% xlabel('\tau (ms)')
% ylabel('R(\tau)')
% nexttile([1 1]);
% plot(t_lag,AAxcorr(foi,:)','color',[0.8 0.8 0.8]); hold all
% errorbar_shade(t_lag,mean(AAxcorr(foi,:),1,'omitnan'),std(AAxcorr(foi,:),0,1,'omitnan'));
% xlabel('\tau (ms)')
% ylabel('R(\tau)')

figure(128); clf; cmap=distinguishable_colors(2);
[mean_PSDbs std_PSDbs b_cent]=binning_data(PSD_sub(foi,1),[0:1:20]);
[mean_PSDap std_PSDap b_cent]=binning_data(PSD_sub(foi,2),[0:1:20]);
errorbar_shade(b_cent,mean_PSDbs,std_PSDbs,cmap(1,:)); hold all
errorbar_shade(b_cent,mean_PSDap,std_PSDap,cmap(2,:));

%
% sample_gridSS=linspace(-1,2,50);
% sample_gridCS=linspace(-1,2,30);
% interp_grid=linspace(-1,2,100);
% [Xq,Yq] = meshgrid(interp_grid,interp_grid);
%
% clf;
% nexttile([1 1])
% [Prob_preCSsubBA xc yc]=scatter_heatmap2(preCS_SubBA(1,:),preCS_SubBA(2,:),sample_gridCS,sample_gridCS);
% shading faceted
% [Xmesh, Ymesh] = meshgrid(xc,yc);
% Prob_preCSsubBA=interp2(Xmesh,Ymesh,Prob_preCSsubBA,Xq,Yq);
%
% nexttile([1 1])
% [Prob_preSSsubBA xc yc]=scatter_heatmap2(preSS_SubBA(1,:),preSS_SubBA(2,:),sample_gridSS,sample_gridSS);
% [Xmesh, Ymesh] = meshgrid(xc,yc);
% Prob_preSSsubBA=interp2(Xmesh,Ymesh,Prob_preSSsubBA,Xq,Yq);
% colormap(turbo)
% ax = gca; % Get current axes
% ax.GridColor = 'w'; % Set grid color (e.g., black)
% ax.GridAlpha = 0.5; % Set transparency of the grid lines

%% bAP amplitude in place fields
nTau_short=[1 3];
nTau_show=[15 5];

placefield_bAPAMP=[];
placefield_SubAmp=[];
placeTfield_bAPAMP=[];
placefield_spClassVec=[];
PF_center=[];
Placefield_LapFR=[];

t_bin=[-200 200:400:4000];

for f=foi
    nTime=size(Subthreshold{f},2);
    roisD_order_ind=cellfun(@find,roisD_order{f},'UniformOutput',false);

    if isempty(cell2mat(roisD_order_ind(1,:)'))
        basalind=cell2mat(roisD_order_ind(2,:)'); %if there is no basal, use soma
    else
        basalind=cell2mat(roisD_order_ind(1,:)');
    end
    apicalind=cell2mat(roisD_order_ind(5,:)'); %apical = distal dend
    somaind=cell2mat(roisD_order_ind(2,:)');

    %[~, bAP_MAT sp_time]=get_STA(NormalizedTrace_dirt{f},max(allSpikeClassMat{f}),1,3);
    [bAP_STA, bAP_MAT sp_time]=get_STA(NormalizedTrace_dirt{f},allSpikeMat{f}(1,:).*double(BlueStim{f}==0),nTau_short(1),nTau_short(2));
    [bAP_STA_Bon, bAP_MAT_Bon sp_time_Bon]=get_STA(NormalizedTrace_dirt{f},allSpikeMat{f}(1,:).*double(BlueStim{f}>0),nTau_short(1),nTau_short(2));
    [bAP_Amp, shift]=max(bAP_MAT,[],3);
    %bAP_Amp = bAP_Amp-mean(bAP_Amp(somaind,:),1,'omitnan');

    %     apicalIntSp=find(max(shift(apicalind,:),[],1)<3 & mean(bAP_STA));
    %     basalIntSp=find(max(shift(basalind,:),[],1)<3);
    %
    % apicalIntSp_vec=ind2vec(nTime,sp_time(apicalIntSp),1);
    % basalIntSp_vec=ind2vec(nTime,sp_time(basalIntSp),1);

    % [~, AP_INT_MAT]=get_STA(NormalizedTrace_dirt{f},apicalIntSp_vec,nTau_show(1),nTau_show(2));
    % [~, BS_INT_MAT]=get_STA(NormalizedTrace_dirt{f},basalIntSp_vec,nTau_show(1),nTau_show(2));

    subMat=tovec(permute(bAP_Amp,[1 3 2]));

    if ~isempty(PlaceFieldList{f}) % in place field
        binTrack=(ceil(VRtrack{f}(5,:)/((115)/150)));
        PFvec=zeros(1,nTime);
        for p=1:length(PlaceFieldBin{f})/2
            if PlaceFieldBin{f}(2*(p-1)+1)>PlaceFieldBin{f}(2*(p-1)+2) %place field includes teleport
                Pvec=~(binTrack<(PlaceFieldBin{f}(2*(p-1)+1)) & binTrack>(PlaceFieldBin{f}(2*(p-1)+2)));
            else
                Pvec=(binTrack>(PlaceFieldBin{f}(2*(p-1)+1)) & binTrack<(PlaceFieldBin{f}(2*(p-1)+2)));
            end
            Lapvec=(VRtrack{f}(8,:)>PlaceFieldList{f}(2*(p-1)+1) & VRtrack{f}(8,:)<PlaceFieldList{f}(2*(p-1)+2));
            PFvec=PFvec| (Lapvec & Pvec);

            bin_dist=ceil(VRtrack{f}(5,:)/(115/50));
            PF_rescale=[floor(PlaceFieldBin{f}(2*(p-1)+[1])/150*50) ceil(PlaceFieldBin{f}(2*(p-1)+[2])/150*50)];

            bin_STA=[];
            for b=1:max(bin_dist)
                bin_STA(:,:,b)=get_STA(NormalizedTrace_dirt{f},max(allSpikeMat{f}).*double(BlueStim{f}==0).*double(bin_dist==b).*double(PFvec),nTau_short(1),nTau_short(2));
            end

            tbin_STA=[];
            first_PFspike_time=find(get_spikeOrder(PFvec,max(allSpikeMat{f}))==1);
            for t=1:length(t_bin)-1
                t_vec=ind2vec(nTime,tovec(first_PFspike_time'+[t_bin(t):t_bin(t+1)-1]),1);
                tbin_STA(:,:,t)=get_STA(NormalizedTrace_dirt{f},max(allSpikeMat{f}).*double(BlueStim{f}==0).*double(t_vec),nTau_short(1),nTau_short(2));
            end

            bin_STA_Amp=squeeze(max(bin_STA,[],2));
            tbin_STA_Amp=squeeze(max(tbin_STA,[],2));
            tmp=repmat(bin_STA_Amp,1,2); tmp2=repmat(LapSub{f},1,2); tmp3=repmat(LapSpclassVec{f},1,2);

            Placefield_LapFR{f,p}=sum(LapSpclassVec{f}([PlaceFieldList{f}(2*(p-1)+1):PlaceFieldList{f}(2*(p-1)+2)],:,:),3,'omitnan');
            if PlaceFieldBin{f}(2*(p-1)+1)>PlaceFieldBin{f}(2*(p-1)+2) %place field includes teleport
                bin_rescale=([PF_rescale(1):PF_rescale(2)+50]-PF_rescale(1))/50*2;
                placefield_bAPAMP{f,p}=[[NaN bin_rescale]; [dendaxis{f}' tmp(:,PF_rescale(1):PF_rescale(2)+50)]];
                placefield_SubAmp{f,p}=permute(mean(tmp2([PlaceFieldList{f}(2*(p-1)+1):PlaceFieldList{f}(2*(p-1)+2)],[PlaceFieldBin{f}(2*(p-1)+1):PlaceFieldBin{f}(2*(p-1)+2)+150],:),1,'omitnan'),[3 2 1]);
                placefield_SubAmp{f,p}=[[NaN [PlaceFieldBin{f}(2*(p-1)+1):PlaceFieldBin{f}(2*(p-1)+2)+150]-PlaceFieldBin{f}(2*(p-1)+1)]; [dendaxis{f}' placefield_SubAmp{f,p}]];
                placefield_spClassVec{f,p}=tmp3([PlaceFieldList{f}(2*(p-1)+1):PlaceFieldList{f}(2*(p-1)+2)],[PlaceFieldBin{f}(2*(p-1)+1):PlaceFieldBin{f}(2*(p-1)+2)+150],:);
                PF_center(f,p)=mod(mean([PlaceFieldBin{f}(2*(p-1)+1) PlaceFieldBin{f}(2*(p-1)+2)+150]),150)/150*2;
            else
                bin_rescale=([PF_rescale(1):PF_rescale(2)]-PF_rescale(1))/50*2;
                placefield_bAPAMP{f,p}=[[NaN bin_rescale]; [dendaxis{f}' tmp(:,PF_rescale(1):PF_rescale(2))]];
                placefield_SubAmp{f,p}=permute(mean(tmp2([PlaceFieldList{f}(2*(p-1)+1):PlaceFieldList{f}(2*(p-1)+2)],[PlaceFieldBin{f}(2*(p-1)+1):PlaceFieldBin{f}(2*(p-1)+2)],:),1,'omitnan'),[3 2 1]);
                placefield_SubAmp{f,p}=[[NaN [PlaceFieldBin{f}(2*(p-1)+1):PlaceFieldBin{f}(2*(p-1)+2)]-PlaceFieldBin{f}(2*(p-1)+1)]; [dendaxis{f}' placefield_SubAmp{f,p}]];
                %placefield_spClassVec{f,p}=permute(mean(tmp3([PlaceFieldList{f}(2*(p-1)+1):PlaceFieldList{f}(2*(p-1)+2)],[PlaceFieldBin{f}(2*(p-1)+1):PlaceFieldBin{f}(2*(p-1)+2)],:),1,'omitnan'),[3 2 1]);
                placefield_spClassVec{f,p}=tmp3([PlaceFieldList{f}(2*(p-1)+1):PlaceFieldList{f}(2*(p-1)+2)],[PlaceFieldBin{f}(2*(p-1)+1):PlaceFieldBin{f}(2*(p-1)+2)],:);
                PF_center(f,p)=mod(mean([PlaceFieldBin{f}(2*(p-1)+1) PlaceFieldBin{f}(2*(p-1)+2)]),150)/150*2;
            end
            %placeTfield_bAPAMP{f,p}=[[NaN mean([t_bin(2:end); t_bin(1:end-1)])]; [dendaxis{f}' tbin_STA_Amp]];
            %NormTr=[mean(NormalizedTrace_dirt{f}(somaind,:),1,'omitnan'); mean(NormalizedTrace_dirt{f}(apicalind,:),1,'omitnan')];
            meanAmp= max(bAP_STA,[],2);
            Normcst= [mean(meanAmp(somaind,:),1,'omitnan'); mean(meanAmp(apicalind,:),1,'omitnan')];
            placeTfield_bAPAMP{f,p}=[mean(tbin_STA_Amp(somaind,:),1,'omitnan'); mean(tbin_STA_Amp(apicalind,:),1,'omitnan')]./Normcst;
        end
    end
end
PF_center(PF_center==0)=NaN;

v=cellfun(@(x) ~isempty(x),placefield_bAPAMP);
PFbAPdat=placefield_bAPAMP(v);

PFbAPdat_normPF=[]; dbin=80;
for p=1:length(PFbAPdat)
    PFbAPdat_normPF{p}=PFbAPdat{p};
    PFbAPdat_normPF{p}(1,2:end)=rescale(PFbAPdat_normPF{p}(1,2:end));
end

figure(31); clf;
imshow_patch(PFbAPdat_normPF,[0.5 2]);
colormap(turbo)
figure(36); clf;
PF_show=[1 2 3 4 7 8 9 10 11 12 13 14];
show3Dbinning(PFbAPdat_normPF(PF_show),[0:0.07:1],[-dbin*3/2:dbin:400],'image')
shading flat; set(gca,'YDir','reverse')
colormap(turbo)
caxis([0.6 1.3])
xlabel('Normalized place field')
ylabel('Distance from soma (\mum)')
cb=colorbar;
cb.Label.String='bAP amplitude';

figure(32); clf;
imshow_patch(placefield_SubAmp(v),[-0.2 0.7]);
colormap(turbo)
PFSubdat=placefield_SubAmp(v);

figure(37); clf; PFSubdat_normPF=[];
for p=1:length(PFSubdat)
    PFSubdat_normPF{p}=PFSubdat{p};
    PFSubdat_normPF{p}(1,2:end)=rescale(PFSubdat_normPF{p}(1,2:end));
end
show3Dbinning(PFSubdat_normPF(PF_show),[0 0.03:0.05:1],[-dbin*3/2:dbin:400],'image')
shading flat; set(gca,'YDir','reverse')
colormap(turbo)
caxis([-0.05 0.25])
xlabel('Normalized place field')
ylabel('Distance from soma (\mum)')
cb=colorbar;
cb.Label.String='Subthreshold';

figure(35); clf;
%imshow_patch(cellfun(@(x) im_merge(ringmovMean(x(:,:,[1 2 3]),3),[1 1 1;1 0 0;1 1 1]), placefield_spClassVec(v),'UniformOutput',false));
tmppf=placefield_spClassVec(v);
for p=1:size(tmppf,1)
    nexttile([1 1])
    l=plot(squeeze(mean(tmppf{p},1,'omitnan')));
    arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell([0 0 0;1 0 0;0 0.5 1],2));
    xlabel('VR position')
    ylabel('Firing rate');
end
legend({'SS','CS','BS'})


figure(33); clf;
% [binnedZ, centerX, centerY]=show3Dbinning(placeTfield_bAPAMP(v),t_bin,[-100:70:450],'image');
% colormap(turbo); caxis([0.5 2]);
% shading flat
% set(gca,'YDir','reverse')
somaTfield=cell2mat(cellfun(@(x) x(1,:),placeTfield_bAPAMP(v),'UniformOutput',false));
apicalTfield=cell2mat(cellfun(@(x) x(2,:),placeTfield_bAPAMP(v),'UniformOutput',false));

dat=apicalTfield;%./max(apicalTfield,[],2);
dat_som=somaTfield;%./somaTfield(:,1);
T=mean([t_bin(2:end); t_bin(1:end-1)])/1000;
M=mean(dat,1,'omitnan'); S=std(dat,0,1,'omitnan'); N=sum(~isnan(dat),1);
Ms=mean(dat_som,1,'omitnan'); Ss=std(dat_som,0,1,'omitnan'); Ns=sum(~isnan(dat_som),1);
plot(T,dat','color',[1 0.7 0.7]);
hold all
plot(T,dat_som','color',[0.7 0.7 1]);
errorbar_shade(T,M,S/sqrt(N),[1 0 0])
errorbar_shade(T,Ms,Ss/sqrt(Ns),[0 0 1])

%% Somatic spike vs dendritic spike, how many spikes are preceding? -> No spikes
nTau_short=[2 3];
FractionMat=[];
figure(38); clf;
for f=foi
    nTime=size(Subthreshold{f},2);
    roisD_order_ind=cellfun(@find,roisD_order{f},'UniformOutput',false);

    if isempty(cell2mat(roisD_order_ind(1,:)'))
        basalind=cell2mat(roisD_order_ind(2,:)'); %if there is no basal, use soma
    else
        basalind=cell2mat(roisD_order_ind(1,:)');
    end
    apicalind=cell2mat(roisD_order_ind(5,:)'); %apical = distal dend
    somaind=cell2mat(roisD_order_ind(2,:)');

    [bAP_STA, bAP_MAT sp_time]=get_STA(NormalizedTrace_dirt{f},max(allSpikeClassMat{f}),nTau_short(1),nTau_short(2));
    [~, Sp_STA]=get_STA(allSpikeMat{f},max(allSpikeClassMat{f}),nTau_short(1),nTau_short(2));
    Sp_STA=max(Sp_STA(Dist_order{f}(noi_dist{f}),:,1:nTau_short(1)),[],3);
    %[bAP_STA, bAP_MAT sp_time]=get_STA(NormalizedTrace_dirt{f},allSpikeMat{f}(1,:).*double(BlueStim{f}==0),nTau_short(1),nTau_short(2));
    [bAP_Amp, shift]=max(bAP_MAT,[],3);
    bAP_STA_Amp=max(bAP_STA,[],2);
    %bAP_Amp = bAP_Amp-mean(bAP_Amp(somaind,:),1,'omitnan');

    %average_bAPamp=[max(bAP_STA(basalind,:)) max(mean(bAP_STA(apicalind,:),1,'omitnan'))];
    % apicalIntSp=find(min(shift(apicalind,:),[],1)<(nTau_short(1)+1) & sum(bAP_Amp(apicalind,:) > bAP_STA_Amp(apicalind)*3)>1);
    % basalIntSp=find(min(shift(basalind,:),[],1)<(nTau_short(1)+1) & sum(bAP_Amp(basalind,:) > bAP_STA_Amp(basalind)*3)>1);
    apicalIntSp=find(max(Sp_STA(apicalind,:),[],1)>0);
    basalIntSp=find(max(Sp_STA(basalind,:),[],1)>0);

    apicalIntSp_vec=ind2vec(nTime,sp_time(apicalIntSp),1);
    basalIntSp_vec=ind2vec(nTime,sp_time(basalIntSp),1);

    [bAP_STA_show ]=get_STA(NormalizedTrace_dirt{f},max(allSpikeClassMat{f}),nTau_show(1),nTau_show(2));
    [~, AP_INT_MAT]=get_STA(NormalizedTrace_dirt{f},apicalIntSp_vec,nTau_show(1),nTau_show(2));
    [~, BS_INT_MAT]=get_STA(NormalizedTrace_dirt{f},basalIntSp_vec,nTau_show(1),nTau_show(2));

    FractionMat(f,:)=[sum(basalIntSp_vec), sum(apicalIntSp_vec), sum(basalIntSp_vec & apicalIntSp_vec)]/sum(max(allSpikeClassMat{f}));

    nexttile([1 1]);
    imagesc([squeeze(mean(BS_INT_MAT,2,'omitnan')) squeeze(mean(AP_INT_MAT,2,'omitnan')) bAP_STA_show])
    title(num2str(f))
end
colormap(turbo)

%%