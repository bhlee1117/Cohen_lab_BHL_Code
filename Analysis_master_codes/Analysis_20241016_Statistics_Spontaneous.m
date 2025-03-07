clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:W31');

ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,10),'UniformOutput',false);
basal_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,11),'UniformOutput',false);
apical_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false);
fpath=raw(:,1)';
StructureData=raw(:,10);
BadROI=cellfun(@(x) (str2num(num2str(x))),raw(:,15),'UniformOutput',false);
EndFrame=cell2mat(raw(:,13));
ifmotionReject=cell2mat(raw(:,14));
ifdirtRemov=cell2mat(raw(:,16));
Pixelsize=cell2mat(raw(:,6));
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
PlaceFieldList=cellfun(@(x) (str2num(num2str(x))),raw(:,19),'UniformOutput',false);
PlaceFieldBin=cellfun(@(x) (str2num(num2str(x))),raw(:,20),'UniformOutput',false);
set(0,'DefaultFigureWindowStyle','docked')
foi=[1 4 5 6 8 10 11 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27];
%foi=[1 13 17];%7 9 17
%foi=[6 7 12 13 14 22];
% %%
% for f=8
%     load(fullfile(fpath{f},'PC_Result.mat'),'Result')
%     interDendDist=[];
%     SkelDend = Skeletonize_dendrite(Result.ref_im,6,0.02,10);
%     nROI=size(Result.ftprnt,3);
%     for i=1:nROI
%         i
%         for j=1:nROI
%             [interDendDist(i,j), path]=geodesic_distance(SkelDend,get_coord(Result.ftprnt(:,:,i)),get_coord(Result.ftprnt(:,:,j)));
%         end
%     end
%     Result.interDendDist=interDendDist;
%
%     coord_1d=dim_reduce(get_coord(Result.ftprnt));
%     som_coord=get_coord(Result.ftprnt(:,:,1));
%     if som_coord(1)<size(Result.ref_im,2)/2
%         [~, Result.dist_order]=sort(coord_1d,'ascend');
%     else
%         [~, Result.dist_order]=sort(coord_1d,'descend');
%     end
%
%     save(fullfile(fpath{f},'PC_Result.mat'),'Result')
% end
%
% %%
% for f=23:25%length(fpath)
%     load(fullfile(fpath{f},'PC_Result.mat'),'Result');
%     figure(f+100); clf;
%     nexttile([1 1]);
%     rois={basal_ROI{f},apical_ROI{f},ref_ROI{f}};
%     ft_roi=[];
%     for r=1:3
%         if ~isnan(rois{r})
%             ft_roi(:,:,r)=max(Result.ftprnt(:,:,rois{r}),[],3);
%         end
%     end
%     show_footprnt_contour(ft_roi,Result.ref_im)
%     nexttile([1 1]);
%     plot(mean(Result.traces(ref_ROI{f},:),1));
%     drawnow;
% end
%%
nTau={[-70:50],[-70:50],[-70:50]}; %SS, CS, Brst
clear SpikeInd SpikeMat STAmat NormalizedTrace_ch NormalizedTrace_dirt SpikeIndBlueOff Dist_order allSpikeMat noi interDendDist noi_dist derivSub LapSubSilent
clear Subthreshold

for f=23
    f
    load(fullfile(fpath{f},'PC_Result.mat'),'Result')
    rois={basal_ROI{f},apical_ROI{f},ref_ROI{f}};
    for r=1:length(rois)
        if isnan(rois{r}); rois{r}=[]; end;
    end
    Dist_order{f}=Result.dist_order;
    nROI=size(Result.normTraces,1);
    nTau_bAP=[-30:20];

    StimOn_Lap=[];
    stimlap=Result.VR(8,:).*double(Result.Blue>0);
    for b=1:max(Result.VR(8,:))
        if sum(stimlap==b)>100
            StimOn_Lap=[StimOn_Lap b];
        end
    end
    StimOff_Lap=setdiff(unique(Result.VR(8,:)),StimOn_Lap);

    % spike class, SS:1, CS:2, dSP:3, BS:4
    som_spike=find(Result.spike(1,:)); ss_time=find(Result.SpClass(1,:)); % BS is subclass of SS
    brst=bwlabel((ss_time(2:end)-ss_time(1:end-1))<15); % SSs that have an ISI shorter than 15 ms are BS.
    SpClass=Result.SpClass; BS_trace=zeros(1,size(Result.traces,2));
    for b=1:max(brst)
        bwn=find(brst==b);
        SpClass(1,ss_time([bwn bwn(end)+1]))=0;
        SpClass(4,ss_time([bwn(1)]))=1;
        BS_trace(1,[ss_time(bwn): ss_time(bwn(end)+1)])=b;
    end

    bAP_ref=[]; % This is to calculate F_ref
    for s=som_spike
        if sum((s+nTau_bAP)<0 | (s+nTau_bAP)>size(Result.traces,2),2)==0
            isnearby=sum(ismember(s+nTau_bAP,som_spike))>1;
            isnearbyCS=sum(ismember(s+nTau_bAP,find(Result.CStrace)))>1;
            ispartCS=sum(Result.CStrace(s+nTau_bAP))>0;
            if ~isnearby & ~isnan(s) & ~isnearbyCS & ~ispartCS
                bAP_ref=[bAP_ref s];
            end
        end
    end
    sp_na=sum((bAP_ref'+nTau_bAP)<0 | (bAP_ref'+nTau_bAP)>size(Result.traces,2),2)==0;
    bAP_ref=bAP_ref(sp_na);

    % Isolated Somatic spike
    SS_s=[];
    for s=som_spike
        if sum((s+nTau{1})<0 | (s+nTau{1})>size(Result.traces,2),2)==0
            isnearby=sum(ismember(s+nTau{1},som_spike))>1;
            isnearbyCS=sum(ismember(s+nTau{1},find(Result.CStrace)))>1;
            ispartCS=sum(Result.CStrace(s+nTau{1}))>0;
            if ~isnearby & ~isnan(s) & ~isnearbyCS & ~ispartCS
                SS_s=[SS_s s];
            end
        end
    end
    sp_na=sum((SS_s'+nTau{1})<0 | (SS_s'+nTau{1})>size(Result.traces,2),2)==0;
    SS_s=SS_s(sp_na);

    % Isolated Complex spike
    CS_s=[];
    C_list=find(Result.SpClass(2,:));
    CS_label=max(Result.spike,[],1).*bwlabel(Result.CStrace);
    CS_list=[];
    for g=1:length(C_list)
        s=C_list(g);
        if sum((s+nTau{2})<0 | (s+nTau{2})>size(Result.traces,2),2)==0
            s_tmp=Result.spike(1,:);
            s_tmp(find(CS_label==g))=0;
            isnearbyCS=max(bwlabel(Result.CStrace(s+nTau{2})))>1;
            isnearbyS=sum(ismember(s+nTau{2},find(s_tmp)))>0;
            if ~isnearbyCS & ~isnearbyS
                CS_s=[CS_s s];
                CS_list=[CS_list g];
            end
        end
    end
    if ~isempty(CS_s)
        sp_na=sum((CS_s'+nTau{2})<0 | (CS_s'+nTau{2})>size(Result.traces,2),2)==0;
        CS_s=CS_s(sp_na);
    end

    % Isolated Burst spike
    BS_s=[];
    B_list=find(SpClass(4,:));
    for g=1:length(B_list)
        s=B_list(g);
        isnearby=sum(ismember(s+nTau{3},find(BS_trace~=BS_trace(s) & Result.spike(1,:)==1)))>0;
        isnearbyCS=sum(ismember(s+nTau{3},find(Result.CStrace)))>0;
        if ~isnearbyCS & ~isnearby
            BS_s=[BS_s s];
        end
    end
    sp_na=sum((BS_s'+nTau{3})<0 | (BS_s'+nTau{3})>size(Result.traces,2),2)==0;
    BS_s=BS_s(sp_na);

    STA_SS=squeeze(mean(reshape(Result.normTraces(:,bAP_ref'+nTau_bAP),nROI,[],length(nTau_bAP)),2,'omitnan'));
    STA_SS= STA_SS - prctile(STA_SS,10,2);
    %STA_SS= STA_SS - min(STA_SS(:,1:10),2,'omitnan');

    %F_ref=max(STA_SS,[],2)./Result.dFF_slope';

    %plot(mean(STA_SS(:,-nTau_bAP(1)+[7:11]),2),max(STA_SS,[],2)./Result.dFF_slope','.','MarkerSize',20); hold all
    %F_ref=mean(STA_SS(:,-nTau_bAP(1)+[7:10]),2);
    %F_ref=(tovec(imgaussfilt(Result.ref_im,1))'*tovec(Result.ftprnt)/Result.SpikeHeight_fit(1))';
    %F_ref=-median(Result.traces_bvMask,2,'omitnan');
    F_ref=Result.F_ref;

    SilentPeriod=ones(1,size(Result.traces,2));
    sp_time=find(max(Result.spike,[],1))';
    sp_na=sum((find(max(Result.spike,[],1))'+[-10:150])<0 | (find(max(Result.spike,[],1))'+[-10:150])>size(Result.traces,2),2)==0;
    SilentPeriod(sp_time(sp_na)+[-10:150])=NaN;
    t_fit=find(~isnan(SilentPeriod) & Result.Blue==0);

    NormalizedTrace=(Result.normTraces)./F_ref;

    NormalizedTrace_dirt{f,1}=NormalizedTrace;
    NormalizedTrace_dirt{f,1}(:,Result.motionReject>0)=NaN;
    NormalizedTrace_ch(f,:)=cellfun(@(x) x./F_ref,Result.norm_trace_check,'UniformOutput',false);
    NormalizedTrace_ch{f,1}(:,Result.motionReject>0)=NaN; NormalizedTrace_ch{f,2}(:,Result.motionReject>0)=NaN;
    if ifdirtRemov(f)
        NormalizedTrace_dirt{f,1}(Result.dirtTrace>0)=NaN;
        NormalizedTrace_ch{f,1}(Result.dirtTrace>0)=NaN; NormalizedTrace_ch{f,2}(Result.dirtTrace>0)=NaN;
    end

    S_list=[];
    S_list={SS_s,CS_s,BS_s};

    SpikeInd(f,:)=S_list';
    for stype=1:3
        if ~isempty(S_list{stype})
            STAmat{f,stype}=reshape(NormalizedTrace_dirt{f,1}(:,S_list{stype}'+nTau{stype}),nROI,[],length(nTau{stype}));
            SpikeMat{f,stype}=permute(reshape(Result.spike(1,S_list{stype}'+nTau{stype}),1,[],length(nTau{stype})),[2 3 1]);
        else
            STAmat{f,stype}=NaN(nROI,1,length(nTau{stype}));
            SpikeMat{f,stype}=NaN(nROI,length(nTau{stype}));
        end
        SpikeIndBlueOff{f,stype}=find(Result.Blue(S_list{stype})==0);
    end

    allSpikeMat{f}=Result.spike;
    allSpikeMat{f}(:,Result.motionReject>0)=NaN;
    if ifdirtRemov(f)
        allSpikeMat{f}(Result.dirtTrace>0)=NaN;
    end

    allSpikeClassMat{f}=SpClass;

    interDendDist{f}=Result.interDendDist*Pixelsize(f);
    Subthreshold{f}=get_subthreshold(NormalizedTrace_dirt{f},max(allSpikeMat{f}(1,:),[],1)>0,7,17);
    Subthreshold{f}(isnan(NormalizedTrace_dirt{f}))=NaN;

    BlueStim{f}=Result.Blue;
    VRtrack{f}=Result.VR;
    CStrace{f}=Result.CStrace;
    Ftprnts{f}=Result.ftprnt;

    %noi=unique([1 setdiff([1:size(NormalizedTrace_dirt{f},1)],find(Result.dFF_slopefit<0.4 | Result.dFF_slope<0))]); %Filterting based on dFF fit result
    %noi=setdiff(unique([1 cell2mat(rois)]),BadROI{f});
    noi=setdiff([1:size(NormalizedTrace_dirt{f},1)],BadROI{f});
    noi_dist{f}=ismember(Dist_order{f},noi);

    AvgImg{f}=Result.ref_im;

    % for n=1:size(Subthreshold{f},1)
    %     n
    % derivSub{f}(n,:)=get_slope(NormalizedTrace_dirt{f}(n,:),7);
    % end

    LapFR{f}=PlaceTrigger_average(double(allSpikeMat{f}(1,:)==1),150,VRtrack{f},0.002,115); %total trace
    LapV{f}=PlaceTrigger_average(NormalizedTrace_dirt{f},150,VRtrack{f},0.002,115); %total trace
    LapSub{f}=PlaceTrigger_average(Subthreshold{f},150,VRtrack{f},0.002,115); %total trace

    subthreshold_silent=Subthreshold{f};
    subthreshold_silent(:,unique(get_perispikeIndex(allSpikeMat{f}(1,:),[-8:30])))=NaN;
    LapSubSilent{f}=PlaceTrigger_average(subthreshold_silent,150,VRtrack{f},0.002,115); %total trace
end

STAmat=cellfun(@(x) permute(x,[1 3 2]),STAmat,'UniformOutput',false);
%% show STA kymo
figure(101); clf; tiledlayout(9,9)
stype_str={'SS','CS','BS'};
sub_time=[-69:-40]; %from spike
ax1=[];
for f=23
    refkymo=mean(STAmat{f,1}(Dist_order{f}(noi_dist{f}),:,SpikeIndBlueOff{f,1})-median(STAmat{f,1}(Dist_order{f}(noi_dist{f}),-nTau{1}(1)+sub_time,SpikeIndBlueOff{f,1}),2,'omitnan'),3,'omitnan');
    cax=[prctile(refkymo(:),0.5) prctile(refkymo(:),99.5)];

    Dsign=ones(1,size(interDendDist{f},1));
    Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
    for stype=1:3
        %sub_time=-nTau{stype}(1)+[-10:-2];
        if ~isempty(SpikeIndBlueOff{f,stype}) % If there is no spike
            imshow=mean(STAmat{f,stype}(Dist_order{f}(noi_dist{f}),:,SpikeIndBlueOff{f,stype})-median(STAmat{f,stype}(Dist_order{f}(noi_dist{f}),-nTau{stype}(1)+sub_time,SpikeIndBlueOff{f,stype}),2,'omitnan'),3,'omitnan');
        else
            imshow=zeros(sum(noi_dist{f}),length(nTau{stype}));
        end
        %imshow=mean(STAmat{f,stype}(:,:,SpikeIndBlueOff{f,sty  pe}),3,'omitnan');
        if isnan(cax); cax=[0 1]; end;
        ax1=[ax1 nexttile([1 1])];
        imagesc(imshow,cax)
        title(['Cell# ' num2str(f) ', ' stype_str{stype} ', N=' num2str(length(SpikeIndBlueOff{f,stype}))])
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

%% Measure the Spike attenuation rate
clear SpikeAmp SpikeAmp_kink STAmat_crop STAmat_crop_kink SpikeDelay SpikeOrder AllSpikeAmp AllSpikeAmp_kink preSub AllSpikepreSub Allactivity preSubprox
pre_spike_window=[-5:-1];
pre_spike_window2=[-20:-1];
Amp_window=[0];
for f=foi
    nROI=size(NormalizedTrace_dirt{f},1);
    spikes=find(max(allSpikeMat{f}(1,:),[],1)>0 & BlueStim{f}==0);
    sp_na=sum((spikes'+[-5:4])<0 | (spikes'+[-5:4])>size(NormalizedTrace_dirt{f},2),2)==0;
    spikes=spikes(sp_na);

    catMat=[];
    for stype=1:3
        if ~isempty(STAmat{f,stype})
            catMat=cat(3,catMat,STAmat{f,stype});
        end
    end
    AllSTA=mean(catMat(Dist_order{f}(noi_dist{f}),:,:),3,'omitnan');
    AllSTA=AllSTA-mean(AllSTA(:,1:50),2,'omitnan');
    AllSTA=AllSTA(:,1:-nTau{1}(1)+5);
    [~, maxAmptime]=max(AllSTA,[],2);
    maxAmptime=maxAmptime+nTau{1}(1);

    Allspmat=reshape(NormalizedTrace_dirt{f}(Dist_order{f}(noi_dist{f}),spikes'+[-2:4]),sum(noi_dist{f}),[],7);
    for n=1:size(AllSTA,1)
        AllSpikeAmp{f}(n,:)=mean(Allspmat(n,:,maxAmptime(n)+2+Amp_window),3,'omitnan');
    end
    AllSpikepreSub{f}=mean(reshape(NormalizedTrace_dirt{f}(Dist_order{f}(noi_dist{f}),spikes'+[-5:-1]),sum(noi_dist{f}),[],5),3,'omitnan');
    AllSpikeAmp_kink{f}=AllSpikeAmp{f}-AllSpikepreSub{f};
    Allactivity{f}=max([double(allSpikeMat{f}==1); CStrace{f}],[],1);
    Allactivity{f}(1,isnan(allSpikeMat{f}(1,:)))=NaN;

    for stype=1:3
        preSubprox{f,stype}=mean(STAmat{f,stype}(Dist_order{f}(noi_dist{f}),-nTau{stype}(1)+pre_spike_window,SpikeIndBlueOff{f,stype}),2,'omitnan');
        preSub{f,stype}=mean(STAmat{f,stype}(Dist_order{f}(noi_dist{f}),-nTau{stype}(1)+pre_spike_window2,SpikeIndBlueOff{f,stype}),2,'omitnan');
        SubSpmat=SpikeMat{f,stype}(SpikeIndBlueOff{f,stype},:); % During Blue off

        if ~isempty(SubSpmat)

            SubSpmat(:,end-4:end)=0; % avoid overlap
            spind=find(SubSpmat'); % find all the spikes

            STAmat_perm=STAmat{f,stype}(:,:,SpikeIndBlueOff{f,stype});%-median(STAmat{f,stype}(:,-nTau{stype}(1)+[-69:-20],SpikeIndBlueOff{f,stype}),2,'omitnan');
            STAmat_perm_kink=STAmat{f,stype}(Dist_order{f}(noi_dist{f}),:,SpikeIndBlueOff{f,stype})-preSubprox{f,stype};

            STAmat_perm_crop=reshape(STAmat_perm(Dist_order{f}(noi_dist{f}),spind'+[-2:max(maxAmptime)+1]'),sum(noi_dist{f}),max(maxAmptime)+4,[]);
            STAmat_perm_crop_kink=reshape(STAmat_perm_kink(:,spind'+[-2:max(maxAmptime)+1]'),sum(noi_dist{f}),max(maxAmptime)+4,[]);

            for n=1:size(STAmat_perm_crop,1)
                SpikeAmp{f,stype}(n,:)=mean(STAmat_perm_crop(n,maxAmptime(n)+2+Amp_window,:),2,'omitnan');
                SpikeAmp_kink{f,stype}(n,:)=mean(STAmat_perm_crop_kink(n,maxAmptime(n)+2+Amp_window,:),2,'omitnan');
            end

            SpOrderMat=NaN(length(SpikeIndBlueOff{f,stype}),max(sum(SubSpmat,2)));
            for eventS=1:size(SpOrderMat,1) % N th (SS/CS/BS)
                SpOrderMat(eventS,1:sum(SubSpmat(eventS,:),2))=1;
            end
            [sp_order, sp_event]=ind2sub([size(SpOrderMat,2) size(SpOrderMat,1)],find(SpOrderMat'==1));

            SpikeOrder{f,stype}=sp_order;

            for eventS=1:size(SpOrderMat,2)
                SpikeDelay{f,stype}(:,eventS)=get_delay(squeeze(mean(STAmat_perm_crop(:,:,find(sp_order==eventS)),3,'omitnan')),100);
                STAmat_crop{f,stype}(:,eventS)=tovec(squeeze(mean(STAmat_perm_crop(:,:,find(sp_order==eventS)),3,'omitnan')));
                STAmat_crop_kink{f,stype}(:,eventS)=tovec(squeeze(mean(STAmat_perm_crop_kink(:,:,find(sp_order==eventS)),3,'omitnan')));
            end
            SpikeDelay{f,stype}=SpikeDelay{f,stype}-SpikeDelay{f,stype}(find(Dist_order{f}(noi_dist{f})==1),:);
        end
    end
end

%% plots

% Show triggering average traces
figure(102); clf;
tiledlayout(7,26); g2=1;
for f=foi

    Dsign=ones(1,size(interDendDist{f},1));
    Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
    cax=[prctile(STAmat_crop{f,1}(:,1),1) prctile(STAmat_crop{f,1}(:,1),99)];
    g=1;
    for stype=1:3
        for s=1:min([size(STAmat_crop{f,stype},2) 3])
            nexttile(g2+26*(g-1),[1 1]);
            imagesc(reshape(STAmat_crop{f,stype}(:,s),sum(noi_dist{f}),[]),cax)
            if g==1;
                title(num2str(f))
            end
            g=g+1;

            if find(Dist_order{f}(noi_dist{f})==1)~=1
                set(gca,'YTick',[1 find(Dist_order{f}(noi_dist{f})==1) sum(noi_dist{f})],'YTickLabel',num2str([min(interDendDist{f}(1,:).*Dsign) 0 max((interDendDist{f}(1,:).*Dsign))]',3))
            else
                set(gca,'YTick',[1 sum(noi_dist{f})],'YTickLabel',num2str([0 max((interDendDist{f}(1,:)))]',3))
            end
        end
    end

    g2=g2+1;
end
colormap(turbo)

%% Show conduction speed
figure(103); clf;
conduction_Speed=NaN(max(foi),7);
conduction_Speed_basal=NaN(max(foi),7);
numsp2show=[1 3 3]; %SS, CS, BS
for f=foi
    Dsign=ones(1,size(interDendDist{f},1));
    Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
    dendaxis=interDendDist{f}(1,:).*Dsign;
    %nexttile([1 1]);
    g=1;
    for stype=1:3
        for s=1:numsp2show(stype)

            xaxis=dendaxis(Dist_order{f}(noi_dist{f}));
            ap=find(xaxis>=0);
            bs=find(xaxis<=0);

            if ~isempty(STAmat_crop{f,stype}) && s<=size(STAmat_crop{f,stype},2) && max(dendaxis)>240

                [p S] = polyfit(xaxis(ap), SpikeDelay{f,stype}(ap,s), 1);
                y_fit = polyval(p, xaxis(ap));

                SS_res = sum((SpikeDelay{f,stype}(ap,s)' - y_fit).^2);
                SS_tot = sum((SpikeDelay{f,stype}(ap,s)' - mean(SpikeDelay{f,stype}(ap,s)')).^2);
                R_squared = 1 - (SS_res / SS_tot);

                if length(bs)>2
                    [p_bs S] = polyfit(-xaxis(bs), SpikeDelay{f,stype}(bs,s), 1);
                    y_fit_bs = polyval(p_bs, -xaxis(bs));
                    SS_res = sum((SpikeDelay{f,stype}(bs,s)' - y_fit_bs).^2);
                    SS_tot = sum((SpikeDelay{f,stype}(bs,s)' - mean(SpikeDelay{f,stype}(bs,s)')).^2);
                    R_squared_bs = 1 - (SS_res / SS_tot);
                else
                    p_bs=NaN;
                    R_squared_bs=0;
                end

                if R_squared>0.5
                    conduction_Speed(f,g)=1/p(1);
                else
                    conduction_Speed(f,g)=NaN;
                end

                if R_squared>0.7 && max(abs(xaxis(bs)))>80
                    conduction_Speed_basal(f,g)=1/p_bs(1);
                else
                    conduction_Speed_basal(f,g)=NaN;
                end

                g=g+1;
            else
                conduction_Speed(f,g)=NaN;
                conduction_Speed_basal(f,g)=NaN;
                g=g+1;
            end
        end
    end
end
%conduction_Speed(12,4)=NaN;
cmap=hsv(7);
%boxplot(reshape([conduction_Speed; conduction_Speed_basal],[],14), 'PlotStyle', 'traditional', 'Colors', 'k'); hold all
boxplot(conduction_Speed, 'PlotStyle', 'traditional', 'Colors', 'k'); hold all
for j = 1:size(conduction_Speed,2)
    % Scatter plot points for each category with a slight horizontal offset for visibility
    jitter = (rand(size(conduction_Speed,1), 1) - 0.5) * 0.1; % Small random jitter for horizontal spacing
    scatter(j + jitter, conduction_Speed(:, j), 35, cmap(j,:), 'filled', 'MarkerFaceAlpha', 0.6);

    %jitter2 = (rand(size(conduction_Speed,1), 1) - 0.5) * 0.1; % Small random jitter for horizontal spacing
    %scatter(2*j + jitter2, conduction_Speed_basal(:, j), 25, cmap(j,:), 'filled', 'MarkerFaceAlpha', 0.6);

    %plot([2*j-1 + jitter 2*j + jitter2]',[conduction_Speed(:,j) conduction_Speed_basal(:,j)]','marker','o','color',[cmap(j,:) 0.5])
end
%set(gca,'XTick',[1:size(conduction_Speed,2)*2],'XTickLabel',{'SS ap','SS bs','CS 1st ap','CS 1st bs','CS 2nd ap','CS 2nd bs','CS 3rd ap','CS 3rd bs','BS 1st ap','BS 1st bs','BS 2nd ap','BS 2nd bs','BS 3rd ap','BS 3rd bs'})
set(gca,'XTick',[1:size(conduction_Speed,2)*2],'XTickLabel',{'SS','CS 1st','CS 2nd','CS 3rd','BS 1st','BS 2nd','BS 3rd'})
ylabel('Conduction speed (\mum/ms)')

%% Show attenuation of spike
dendriteaxis_bin=[-360:100:550];
ylimval=[0 1.4];
figure(104); clf; tiledlayout(2,3)
NormAmp_X=[]; perisoma_distance=30;
for f=foi
    nexttile(1,[1 1])
    Dsign=ones(1,size(interDendDist{f},1));
    Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
    dendaxis=interDendDist{f}(1,:).*Dsign;
    dendaxis=dendaxis(Dist_order{f}(noi_dist{f}));

    [~, maxROI]=max(mean(SpikeAmp{f,1},2,'omitnan'));
    perisomaROI=find(abs(dendaxis-dendaxis(maxROI))<perisoma_distance);
    Normconst=mean(SpikeAmp{f,1}(perisomaROI,:),[1 2],'omitnan');

    [dendaxis indaxis]=sort(dendaxis,'ascend');
    som_ROI=find(Dist_order{f}(noi_dist{f})==1);
    NormAmp=mean(SpikeAmp{f,1},2,'omitnan');
    %NormAmp=NormAmp./NormAmp(som_ROI);
    %NormAmp=NormAmp./max(NormAmp);
    NormAmp=NormAmp/Normconst;
    plot(dendaxis,NormAmp(indaxis),'.','color',[0.6 0.6 0.6],'MarkerSize',15); hold all
    NormAmp_X{f}=[dendaxis' NormAmp(indaxis)];
end
[mean_amplitudes std_amplitudes dendBin_center ind]=binning_data(NormAmp_X(foi),dendriteaxis_bin);
N_neuron=sum((cellfun(@sum,ind)),2)';
errorbar(dendBin_center, mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',[1 0 0],'LineWidth',2);
xlabel('Distance from Soma (\mum)')
ylabel('Normalized bAP amplitude of SS')
ylim(ylimval);

NormAmp_X=[];
for f=foi
    nexttile(4,[1 1])
    Dsign=ones(1,size(interDendDist{f},1));
    Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
    dendaxis=interDendDist{f}(1,:).*Dsign;
    dendaxis=dendaxis(Dist_order{f}(noi_dist{f}));

    [~, maxROI]=max(mean(SpikeAmp{f,1},2,'omitnan'));
    perisomaROI=find(abs(dendaxis-dendaxis(maxROI))<perisoma_distance);
    Normconst=mean(SpikeAmp{f,1}(perisomaROI,:),[1 2],'omitnan');

    [dendaxis indaxis]=sort(dendaxis,'ascend');
    som_ROI=find(Dist_order{f}(noi_dist{f})==1);
    NormAmp=mean(SpikeAmp_kink{f,1},2,'omitnan');
    %NormAmp=NormAmp./NormAmp(som_ROI);
    %NormAmp=NormAmp./max(NormAmp);
    NormAmp=NormAmp/Normconst;
    plot(dendaxis,NormAmp(indaxis),'.','color',[0.6 0.6 0.6],'MarkerSize',15); hold all
    NormAmp_X{f}=[dendaxis' NormAmp(indaxis)];
end
[mean_amplitudes std_amplitudes dendBin_center ind]=binning_data(NormAmp_X(foi),dendriteaxis_bin);
N_neuron=sum((cellfun(@sum,ind)),2)';
errorbar(dendBin_center, mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',[1 0 0],'LineWidth',2);
xlabel('Distance from Soma (\mum)')
ylabel('Normalized bAP amplitude of SS (kink)')
ylim(ylimval);

cmap=distinguishable_colors(6); cmap=cmap([1 2 3 4 6 5],:);
for s=1:4
    NormAmp_X=[];
    for f=foi
        ref_ROIdist=ismember(Dist_order{f}(noi_dist{f}),ref_ROI{f});
        nexttile(2,[1 1])
        if ~isempty(SpikeAmp_kink{f,2})
            %Normconst=mean(SpikeAmp{f,2}(ref_ROIdist,find(SpikeOrder{f,2}==s)),[1 2],'omitnan');
            Dsign=ones(1,size(interDendDist{f},1));
            Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
            dendaxis=interDendDist{f}(1,:).*Dsign;
            dendaxis=dendaxis(Dist_order{f}(noi_dist{f}));

            [~, maxROI]=max(mean(SpikeAmp{f,1},2,'omitnan'));
            perisomaROI=find(abs(dendaxis-dendaxis(maxROI))<perisoma_distance);
            Normconst=mean(SpikeAmp{f,2}(perisomaROI,find(SpikeOrder{f,2}==s)),[1 2],'omitnan');

            [dendaxis indaxis]=sort(dendaxis,'ascend');
            som_ROI=find(Dist_order{f}(noi_dist{f})==1);
            NormAmp=mean(SpikeAmp{f,2}(:,find(SpikeOrder{f,2}==s)),2,'omitnan');
            %NormAmp=NormAmp./NormAmp(som_ROI);
            %NormAmp=NormAmp./max(NormAmp);
            NormAmp=NormAmp/Normconst;
            NormAmp_X{f}=[dendaxis' NormAmp(indaxis)];
        else
            NormAmp_X{f}=[NaN NaN];
        end
    end
    [mean_amplitudes std_amplitudes dendBin_center ind]=binning_data(NormAmp_X(foi),dendriteaxis_bin);
    N_neuron=sum((cellfun(@sum,ind)),2)';
    errorbar(dendBin_center, mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(s,:),'LineWidth',2); hold all
end
xlabel('Distance from Soma (\mum)')
ylabel('Normalized bAP amplitude of CS')
legend(counting_string([1:4]))
ylim(ylimval);

cmap=distinguishable_colors(6); cmap=cmap([1 2 3 4 6 5],:);
for s=1:4
    NormAmp_X=[];
    for f=foi
        ref_ROIdist=ismember(Dist_order{f}(noi_dist{f}),ref_ROI{f});
        nexttile(5,[1 1])
        if ~isempty(SpikeAmp_kink{f,2})
            %Normconst=mean(SpikeAmp{f,2}(ref_ROIdist,find(SpikeOrder{f,2}==s)),[1 2],'omitnan');
            Dsign=ones(1,size(interDendDist{f},1));
            Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
            dendaxis=interDendDist{f}(1,:).*Dsign;
            dendaxis=dendaxis(Dist_order{f}(noi_dist{f}));

            [~, maxROI]=max(mean(SpikeAmp{f,1},2,'omitnan'));
            perisomaROI=find(abs(dendaxis-dendaxis(maxROI))<perisoma_distance);
            Normconst=mean(SpikeAmp{f,2}(perisomaROI,find(SpikeOrder{f,2}==s)),[1 2],'omitnan');

            [dendaxis indaxis]=sort(dendaxis,'ascend');
            som_ROI=find(Dist_order{f}(noi_dist{f})==1);
            NormAmp=mean(SpikeAmp_kink{f,2}(:,find(SpikeOrder{f,2}==s)),2,'omitnan');
            %NormAmp=NormAmp./NormAmp(som_ROI);
            %NormAmp=NormAmp./max(NormAmp);
            NormAmp=NormAmp/Normconst;
            NormAmp_X{f}=[dendaxis' NormAmp(indaxis)];
        else
            NormAmp_X{f}=[NaN NaN];
        end
    end
    [mean_amplitudes std_amplitudes dendBin_center ind]=binning_data(NormAmp_X(foi),dendriteaxis_bin);
    N_neuron=sum((cellfun(@sum,ind)),2)';
    errorbar(dendBin_center, mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(s,:),'LineWidth',2); hold all
end
xlabel('Distance from Soma (\mum)')
ylabel('Normalized bAP amplitude of CS (kink)')
legend(counting_string([1:4]))
ylim(ylimval);

cmap=distinguishable_colors(6); cmap=cmap([1 2 3 4 6 5],:);
for s=1:3
    NormAmp_X=[];
    for f=foi
        ref_ROIdist=ismember(Dist_order{f}(noi_dist{f}),ref_ROI{f});
        nexttile(3,[1 1])
        if ~isempty(SpikeAmp_kink{f,3})
            % Normconst=mean(SpikeAmp{f,3}(ref_ROIdist,find(SpikeOrder{f,3}==s)),[1 2],'omitnan');
            Dsign=ones(1,size(interDendDist{f},1));
            Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
            dendaxis=interDendDist{f}(1,:).*Dsign;
            dendaxis=dendaxis(Dist_order{f}(noi_dist{f}));

            [~, maxROI]=max(mean(SpikeAmp{f,1},2,'omitnan'));
            perisomaROI=find(abs(dendaxis-dendaxis(maxROI))<perisoma_distance);
            Normconst=mean(SpikeAmp{f,3}(perisomaROI,find(SpikeOrder{f,3}==s)),[1 2],'omitnan');

            [dendaxis indaxis]=sort(dendaxis,'ascend');
            som_ROI=find(Dist_order{f}(noi_dist{f})==1);
            NormAmp=mean(SpikeAmp{f,3}(:,find(SpikeOrder{f,3}==s)),2,'omitnan');
            %NormAmp=NormAmp./NormAmp(som_ROI);
            %NormAmp=NormAmp./max(NormAmp);
            NormAmp=NormAmp/Normconst;
            NormAmp_X{f}=[dendaxis' NormAmp(indaxis)];
        else
            NormAmp_X{f}=[NaN NaN];
        end
    end
    [mean_amplitudes std_amplitudes dendBin_center, ind, bindat]=binning_data(NormAmp_X(foi),dendriteaxis_bin);
    N_neuron=sum((cellfun(@sum,ind)),2)';
    errorbar(dendBin_center, mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(s,:),'LineWidth',2); hold all
end
xlabel('Distance from Soma (\mum)')
ylabel('Normalized bAP amplitude of BS')
legend(counting_string([1:4]))
ylim(ylimval);

cmap=distinguishable_colors(6); cmap=cmap([1 2 3 4 6 5],:);
for s=1:3
    NormAmp_X=[];
    for f=foi
        ref_ROIdist=ismember(Dist_order{f}(noi_dist{f}),ref_ROI{f});
        nexttile(6,[1 1])
        if ~isempty(SpikeAmp_kink{f,3})
            %Normconst=mean(SpikeAmp{f,3}(ref_ROIdist,find(SpikeOrder{f,3}==s)),[1 2],'omitnan');
            Dsign=ones(1,size(interDendDist{f},1));
            Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
            dendaxis=interDendDist{f}(1,:).*Dsign;
            dendaxis=dendaxis(Dist_order{f}(noi_dist{f}));

            [~, maxROI]=max(mean(SpikeAmp{f,1},2,'omitnan'));
            perisomaROI=find(abs(dendaxis-dendaxis(maxROI))<perisoma_distance);
            Normconst=mean(SpikeAmp{f,3}(perisomaROI,find(SpikeOrder{f,3}==s)),[1 2],'omitnan');

            [dendaxis indaxis]=sort(dendaxis,'ascend');
            som_ROI=find(Dist_order{f}(noi_dist{f})==1);
            NormAmp=mean(SpikeAmp_kink{f,3}(:,find(SpikeOrder{f,3}==s)),2,'omitnan');
            %NormAmp=NormAmp./NormAmp(som_ROI);
            %NormAmp=NormAmp./max(NormAmp);
            NormAmp=NormAmp/Normconst;
            NormAmp_X{f}=[dendaxis' NormAmp(indaxis)];
        else
            NormAmp_X{f}=[NaN NaN];
        end
    end
    [mean_amplitudes std_amplitudes dendBin_center ind]=binning_data(NormAmp_X(foi),dendriteaxis_bin);
    N_neuron=sum((cellfun(@sum,ind)),2)';
    errorbar(dendBin_center, mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(s,:),'LineWidth',2); hold all
end
xlabel('Distance from Soma (\mum)')
ylabel('Normalized bAP amplitude of BS (kink)')
legend(counting_string([1:4]))
ylim(ylimval);

figure(105); clf;
cmap=distinguishable_colors(6);
nexttile([1 1])
for stype=1:3
    NormAmp_X=[];
    for f=foi
        ref_ROIdist=ismember(Dist_order{f}(noi_dist{f}),ref_ROI{f});
        %Normconst=mean(SpikeAmp{f,1}(ref_ROIdist,:),[1 2],'omitnan');
        if ~isempty(SpikeAmp_kink{f,stype})
            Dsign=ones(1,size(interDendDist{f},1));
            Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
            dendaxis=interDendDist{f}(1,:).*Dsign;
            dendaxis=dendaxis(Dist_order{f}(noi_dist{f}));

            [~, maxROI]=max(mean(SpikeAmp{f,1},2,'omitnan'));
            perisomaROI=find(abs(dendaxis-dendaxis(maxROI))<perisoma_distance);
            Normconst=mean(SpikeAmp{f,stype}(perisomaROI,find(SpikeOrder{f,stype}==1)),[1 2],'omitnan');

            [dendaxis indaxis]=sort(dendaxis,'ascend');
            som_ROI=find(Dist_order{f}(noi_dist{f})==1);
            NormAmp=mean(SpikeAmp{f,stype}(:,find(SpikeOrder{f,stype}==1)),2,'omitnan');
            %NormAmp=NormAmp./NormAmp(som_ROI);
            %NormAmp=NormAmp./max(NormAmp);
            NormAmp=NormAmp/Normconst;
            NormAmp_X{f}=[dendaxis' NormAmp(indaxis)];
        else
            NormAmp_X{f}=[NaN NaN];
        end
    end
    [mean_amplitudes std_amplitudes dendBin_center ind]=binning_data(NormAmp_X(foi),dendriteaxis_bin);
    N_neuron=sum((cellfun(@sum,ind)),2)';
    errorbar(dendBin_center, mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(stype,:),'LineWidth',2); hold all
end
xlabel('Distance from Soma (\mum)')
ylabel('Normalized bAP amplitude')
legend({'SS','CS 1st','BS 1st'})
ylim(ylimval);

nexttile([1 1])
for stype=1:3
    NormAmp_X=[];
    for f=foi
        ref_ROIdist=ismember(Dist_order{f}(noi_dist{f}),ref_ROI{f});
        %Normconst=mean(SpikeAmp{f,1}(ref_ROIdist,:),[1 2],'omitnan');
        if ~isempty(SpikeAmp_kink{f,stype})
            Dsign=ones(1,size(interDendDist{f},1));
            Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
            dendaxis=interDendDist{f}(1,:).*Dsign;
            dendaxis=dendaxis(Dist_order{f}(noi_dist{f}));

            [~, maxROI]=max(mean(SpikeAmp{f,1},2,'omitnan'));
            perisomaROI=find(abs(dendaxis-dendaxis(maxROI))<perisoma_distance);
            Normconst=mean(SpikeAmp{f,1}(perisomaROI,:),[1 2],'omitnan');
            %Normconst=mean(SpikeAmp{f,stype}(perisomaROI,find(SpikeOrder{f,stype}==1)),[1 2],'omitnan');

            [dendaxis indaxis]=sort(dendaxis,'ascend');
            som_ROI=find(Dist_order{f}(noi_dist{f})==1);
            NormAmp=mean(SpikeAmp_kink{f,stype}(:,find(SpikeOrder{f,stype}==1)),2,'omitnan');
            %NormAmp=NormAmp./NormAmp(som_ROI);
            %NormAmp=NormAmp./max(NormAmp);
            NormAmp=NormAmp/Normconst;
            NormAmp_X{f}=[dendaxis' NormAmp(indaxis)];
        else
            NormAmp_X{f}=[NaN NaN];
        end
    end
    [mean_amplitudes std_amplitudes dendBin_center ind]=binning_data(NormAmp_X(foi),dendriteaxis_bin);
    N_neuron=sum((cellfun(@sum,ind)),2)';
    errorbar(dendBin_center, mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(stype,:),'LineWidth',2); hold all
end
xlabel('Distance from Soma (\mum)')
ylabel('Normalized bAP amplitude (kink)')
legend({'SS','CS 1st','BS 1st'})
ylim(ylimval);

figure(106); clf;
cmap=distinguishable_colors(6);
nexttile([1 1])
for stype=1:3
    NormAmp_X=[];
    for f=foi
        ref_ROIdist=ismember(Dist_order{f}(noi_dist{f}),ref_ROI{f});
        %Normconst=mean(SpikeAmp{f,1}(ref_ROIdist,:),[1 2],'omitnan');
        if ~isempty(SpikeAmp_kink{f,stype})
            Dsign=ones(1,size(interDendDist{f},1));
            Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
            dendaxis=interDendDist{f}(1,:).*Dsign;
            dendaxis=dendaxis(Dist_order{f}(noi_dist{f}));

            [~, maxROI]=max(mean(SpikeAmp{f,1},2,'omitnan'));
            perisomaROI=find(abs(dendaxis-dendaxis(maxROI))<perisoma_distance);
            Normconst=mean(SpikeAmp{f,1}(perisomaROI,:),[1 2],'omitnan');

            [dendaxis indaxis]=sort(dendaxis,'ascend');
            som_ROI=find(Dist_order{f}(noi_dist{f})==1);
            NormAmp=mean(squeeze(preSub{f,stype}),2,'omitnan');
            %NormAmp=NormAmp./NormAmp(som_ROI);
            NormAmp=NormAmp/Normconst;
            NormAmp_X{f}=[dendaxis' NormAmp(indaxis)];
        else
            NormAmp_X{f}=[NaN NaN];
        end
    end
    [mean_amplitudes std_amplitudes dendBin_center ind]=binning_data(NormAmp_X(foi),dendriteaxis_bin);
    N_neuron=sum((cellfun(@sum,ind)),2)';
    errorbar(dendBin_center, mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(stype,:),'LineWidth',2); hold all
end
xlabel('Distance from Soma (\mum)')
ylabel('Mean Pre-spike amplitude (-20 ms to -1 ms)')
legend({'SS','CS','BS'})

nexttile([1 1])
for stype=1:3
    NormAmp_X=[];
    for f=foi
        if ~isempty(SpikeAmp_kink{f,stype})
            Dsign=ones(1,size(interDendDist{f},1));
            Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
            dendaxis=interDendDist{f}(1,:).*Dsign;
            dendaxis=dendaxis(Dist_order{f}(noi_dist{f}));

            [~, maxROI]=max(mean(SpikeAmp{f,1},2,'omitnan'));
            perisomaROI=find(abs(dendaxis-dendaxis(maxROI))<perisoma_distance);
            Normconst=mean(SpikeAmp{f,1}(perisomaROI,:),[1 2],'omitnan');

            [dendaxis indaxis]=sort(dendaxis,'ascend');
            som_ROI=find(Dist_order{f}(noi_dist{f})==1);
            NormAmp=mean(squeeze(preSubprox{f,stype}),2,'omitnan');
            %NormAmp=NormAmp./NormAmp(som_ROI);
            NormAmp=NormAmp/Normconst;
            NormAmp_X{f}=[dendaxis' NormAmp(indaxis)];
        else
            NormAmp_X{f}=[NaN NaN];
        end
    end
    [mean_amplitudes std_amplitudes dendBin_center ind]=binning_data(NormAmp_X(foi),dendriteaxis_bin);
    N_neuron=sum((cellfun(@sum,ind)),2)';
    errorbar(dendBin_center, mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(stype,:),'LineWidth',2); hold all
end
xlabel('Distance from Soma (\mum)')
ylabel('Mean Pre-spike amplitude (-5 ms to -1 ms)')
legend({'SS','CS','BS'})

figure(108); clf; clear StdSub_X StdSubprox_X
cmap=distinguishable_colors(6);

NormAmp_X=[];
for f=foi
    %ref_ROIdist=ismember(Dist_order{f}(noi_dist{f}),ref_ROI{f});
    %Normconst=mean(SpikeAmp{f,1}(ref_ROIdist,:),[1 2],'omitnan');

    Dsign=ones(1,size(interDendDist{f},1));
    Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
    dendaxis=interDendDist{f}(1,:).*Dsign;
    dendaxis=dendaxis(Dist_order{f}(noi_dist{f}));

    [~, maxROI]=max(mean(SpikeAmp{f,1},2,'omitnan'));
    perisomaROI=find(abs(dendaxis-dendaxis(maxROI))<perisoma_distance);
    Normconst=mean(SpikeAmp{f,1}(perisomaROI,:),[1 2],'omitnan');

    [dendaxis indaxis]=sort(dendaxis,'ascend');
    Presub_cat=cell2mat(cellfun(@squeeze,preSub(f,:),'UniformOutput',false));
    Presubprox_cat=cell2mat(cellfun(@squeeze,preSubprox(f,:),'UniformOutput',false));
    StdSub=std(Presub_cat/Normconst,0,2,'omitnan');
    StdSubprox=std(Presubprox_cat/Normconst,0,2,'omitnan');
    StdSub_X{f}=[dendaxis' StdSub(indaxis)];
    StdSubprox_X{f}=[dendaxis' StdSubprox(indaxis)];
end
[mean_amplitudes, std_amplitudes, dendBin_center, ind, binnedAmp]=binning_data(StdSub_X(foi),dendriteaxis_bin);
[mean_amplitudes2, std_amplitudes2, dendBin_center2, ~, binnedAmp2]=binning_data(StdSubprox_X(foi),dendriteaxis_bin);
N_neuron=sum((cellfun(@sum,ind)),2)';
errorbar(dendBin_center, mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(1,:),'LineWidth',2); hold all
errorbar(dendBin_center2, mean_amplitudes2, std_amplitudes2./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(2,:),'LineWidth',2); hold all
xlabel('Distance from Soma (\mum)')
ylabel(['Std of Pre-spike amplitude'])
legend('-20 ms to -1 ms','-5 ms to -1 ms')
%%
% show bimodality
clear BimCoeff
figure(106); clf; tiledlayout(2,1);
cmap=distinguishable_colors(6);
for f=foi
    BC=[]; SD=[];
    nexttile(1,[1 1])
    Dsign=ones(1,size(interDendDist{f},1));
    Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
    dendaxis=interDendDist{f}(1,:).*Dsign;
    dendaxis=dendaxis(Dist_order{f}(noi_dist{f}));
    [dendaxis indaxis]=sort(dendaxis,'ascend');
    som_ROI=find(ismember(Dist_order{f}(noi_dist{f}),ref_ROI{f}));
    SpAmp_tmp=SpikeAmp{f,1}(:,find(SpikeOrder{f,1}==1)); % SS
    SpAmp_tmp=[SpAmp_tmp SpikeAmp{f,2}(:,find(SpikeOrder{f,2}==1))]; % CS 1st
    SpAmp_tmp=[SpAmp_tmp SpikeAmp{f,3}(:,find(SpikeOrder{f,3}==1))]; % BS 1st

    SpAmp_tmp=SpAmp_tmp./max(SpAmp_tmp,[],1);
    %SpAmp_tmp=(SpAmp_tmp-mean(SpAmp_tmp(som_ROI,:),1,'omitnan'));
    for n=1:size(SpAmp_tmp,1)
        BC(n,1)=bimodality_coefficient(SpAmp_tmp(n,:));
    end
    SD=std(SpAmp_tmp,0,2,'omitnan');
    BimCoeff{f}=[dendaxis' BC(indaxis)];
    Width{f}=[dendaxis' SD(indaxis)];
    plot(dendaxis',BC(indaxis),'.-','color',[0.5 0.5 0.5]); hold all

    nexttile(2,[1 1])
    plot(dendaxis',SD(indaxis),'.-','color',[0.5 0.5 0.5]); hold all
end
[mean_BC std_BC Bin_BC]=binning_data(BimCoeff(foi),[-300:60:500]);
[mean_SD std_SD Bin_SD]=binning_data(Width(foi),[-300:60:500]);
nexttile(1,[1 1])
errorbar(Bin_BC, mean_BC, std_BC, 'o-', 'CapSize', 5,'color',[1 0 0],'LineWidth',2);
xlabel('Distance from Soma (\mum)')
ylabel('Bimodality Coefficient')
nexttile(2,[1 1])
errorbar(Bin_SD, mean_SD, std_SD, 'o-', 'CapSize', 5,'color',[1 0 0],'LineWidth',2);
xlabel('Distance from Soma (\mum)')
ylabel('Standard deviation')

% %%
% f=22;
% STAmovieSS=readBinMov(fullfile(fpath{f},'STAmovieSS.bin'),size(AvgImg{f},1),size(AvgImg{f},2));
% STAmovieCS=readBinMov(fullfile(fpath{f},'STAmovieCS.bin'),size(AvgImg{f},1),size(AvgImg{f},2));
%
%
%
%
% %%
% f=23;
%
% Dsign=ones(1,size(interDendDist{f},1));
% Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
% dendaxis=interDendDist{f}(1,:).*Dsign;
% dendaxis=dendaxis(Dist_order{f}(noi_dist{f}));
% som_ROI=find(ismember(Dist_order{f}(noi_dist{f}),ref_ROI{f}));
% SpAmp_tmp=SpikeAmp_kink{f,1}(:,find(SpikeOrder{f,1}==1)); % SS
% SpAmp_tmp=[SpAmp_tmp SpikeAmp_kink{f,2}(:,find(SpikeOrder{f,2}==1))]; % CS 1st
% SpAmp_tmp=[SpAmp_tmp SpikeAmp_kink{f,3}(:,find(SpikeOrder{f,3}==1))]; % BS 1st
%
% [~, maxroi]=max(mean(SpAmp_tmp,2,'omitnan'));
% dendaxis_adj=abs(dendaxis-dendaxis(maxroi));
% figure(10); clf;
% for j=[100:50:450]
%     d=[j j+50];
%     doi=[maxroi+[-1:1] find(dendaxis_adj>=d(1) & dendaxis_adj<d(2))]; % ratio
%     %doi=[maxroi+[-1:1] find(dendaxis_adj>=0 & dendaxis_adj<d(2))]; % fit ROI slope
%
%     decay_fit=[]; fitparm=[]; expfitparm=[];
%     for sp=1:size(SpAmp_tmp,2)
%         doi_notnan=find(~isnan(SpAmp_tmp(doi,sp)));
%         initial_parm=[3 200 3]; %-> [amp decay offset]
%         lowerbd=[0 100 0]; upperbd=[100 1000 50];
%         %[expfitparm(sp,:) ,decay_fit(sp,:)]=expfit_wBd(dendaxis_adj(doi(doi_notnan))',SpAmp_tmp(doi(doi_notnan),sp),[0:d],initial_parm,lowerbd,upperbd);
%         %fitparm(sp,:) = polyfit(dendaxis_adj(doi(doi_notnan)), SpAmp_tmp(doi(doi_notnan),sp), 1); % slope
%         fitparm(sp,:)=mean(SpAmp_tmp(doi(4:end),sp),1,'omitnan')./mean(SpAmp_tmp(doi(1:3),sp),1,'omitnan'); % ratio
%     end
%     nonnan=find(~isnan(fitparm(:,1)) & ~isnan(mean(SpAmp_tmp(end-4:end-3,:),1,'omitnan'))');
%     nexttile([1 1])
%     scatter(fitparm(:,1),mean(SpAmp_tmp(end-4:end-3,:),1,'omitnan'),15,'filled','MarkerFaceAlpha', 0.8)
%     %title(['Fit from soma to ' num2str(d(2)+dendaxis(maxroi)) '\mum, Corr: ', num2str(corr(-fitparm(nonnan,1),mean(SpAmp_tmp(end-4:end-3,nonnan),1,'omitnan')'),2)]) % slope
%     title(['Amp. from ' num2str(d(1)+dendaxis(maxroi)) ' to ' num2str(d(2)+dendaxis(maxroi)) '\mum' newline 'divided by soma, Corr: ', num2str(corr(fitparm(nonnan,1),mean(SpAmp_tmp(end-4:end-3,nonnan),1,'omitnan')'),2)]) % ratio
%     %xlabel('Slope: Attenuation per \mum'); ylabel('Distal dendrite bAP Amp.')
%     xlabel('Attenuation ratio'); ylabel('Distal dendrite bAP Amp.')
%     axis tight
%     box on
%     xlim([0 1.5])
% end

%% Relation between pre-spike and Spike amplitude
DistBinEdge=[-360:180:650];

figure(107); clf; ax1=[];
cmap=distinguishable_colors(max(foi));
for d=1:length(DistBinEdge)-1
    dbin=[DistBinEdge(d) DistBinEdge(d+1)];
    preAmp=[]; SpAmp=[];
    for f=foi
        Dsign=ones(1,size(interDendDist{f},1));
        Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
        dendaxis=interDendDist{f}(1,:).*Dsign;
        dendaxis=dendaxis(Dist_order{f}(noi_dist{f}));
        dbin_ind=find(dendaxis>=dbin(1) & dendaxis<=dbin(2));

        [~, maxROI]=max(mean(SpikeAmp{f,1},2,'omitnan'));
        perisomaROI=find(abs(dendaxis-dendaxis(maxROI))<perisoma_distance);
        Normconst=mean(SpikeAmp{f,1}(perisomaROI,:),[1 2],'omitnan');

        for stype=1:3
            if ~isempty(dbin_ind) && ~isempty(SpikeAmp_kink{f,stype})
                showx=tovec(squeeze(preSubprox{f,stype}(dbin_ind,:,:)))/Normconst;
                showy=tovec(SpikeAmp{f,stype}(dbin_ind,find(SpikeOrder{f,stype}==1)))/Normconst;
                % sptmp=SpikeAmp{f,stype}(:,find(SpikeOrder{f,stype}==1))./max(SpikeAmp{f,stype}(:,find(SpikeOrder{f,stype}==1)),[],1);
                % showy=tovec(sptmp(dbin_ind,:));
                nonnan=find(~isnan(showx) & ~isnan(showy));
                preAmp=[preAmp; showx(nonnan)];
                SpAmp=[SpAmp; showy(nonnan)];
                %scatter(showx(nonnan),showy(nonnan),10,cmap(f,:),'filled','MarkerFaceAlpha', 0.8); hold all
            end
        end
    end
    if ~isempty(preAmp)
        ax1=[ax1 nexttile([1 1])];
        scatter_heatmap2(preAmp,SpAmp,linspace(-1,1,40),linspace(-1,2,40));
        xlabel(['Mean Pre-spike subthreshold' newline '(from -5 ms to -1 ms, scaled by spike height)'])
        ylabel('bAP amplitude (spike height)')
        title(['From ' num2str(dbin(1)) ' \mum to ' num2str(dbin(2)) ' \mum'])
    end
end
colormap(Aurora);
linkaxes(ax1,'xy')

%% Preceding ISI
max_precedingtime=350;
%time_bin=[1:50:max_precedingtime];
time_bin=[1 20 40 75 100 200 max_precedingtime];
num_bins = length(time_bin) - 1; % Number of time bins
time_bin_strings = strings(1, num_bins);

for t = 1:num_bins
    time_bin_strings(t) = sprintf('Preceding spike is within %d-%d ms', time_bin(t), time_bin(t+1));
end
clear SpAmp_silent
for f=foi
    spikes=find(max(allSpikeMat{f}(1,:),[],1)>0 & BlueStim{f}==0);
    sp_na=sum((spikes'+[-5:4])<0 | (spikes'+[-5:4])>size(NormalizedTrace_dirt{f},2),2)==0;
    spikes=spikes(sp_na);

    sind=find(spikes>max_precedingtime);
    actmat=Allactivity{f}(spikes(sind)'+[-max_precedingtime:-1]);

    %Normconst=mean(SpikeAmp{f,1}(perisomaROI,:),[1 2],'omitnan');

    [~, recent_activitiy] = max(flip(actmat, 2), [], 2);
    recent_activitiy(recent_activitiy==1)=NaN;

    SpAmp_silent{f}=[recent_activitiy AllSpikeAmp{f}(:,sind)'];
end

[~, ~, time_bin_centered, ind]=binning_data(SpAmp_silent,time_bin);

figure(109); clf; clear MeanSpAmp
cmap=distinguishable_colors(length(time_bin_centered));
for b=1:length(time_bin_centered)
    for f=foi
        Dsign=ones(1,size(interDendDist{f},1));
        Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
        dendaxis=interDendDist{f}(1,:).*Dsign;
        dendaxis=dendaxis(Dist_order{f}(noi_dist{f}));

        [~, maxROI]=max(mean(SpikeAmp{f,1},2,'omitnan'));
        perisomaROI=find(abs(dendaxis-dendaxis(maxROI))<perisoma_distance);

        Normconst=mean(SpAmp_silent{f}(ind{b,f},perisomaROI+1),[1 2],'omitnan');

        MeanSpAmp{f}=[dendaxis' mean(SpAmp_silent{f}(ind{b,f},2:end),1,'omitnan')'/Normconst];
    end
    [mean_amplitudes std_amplitudes dendBin_center ind2]=binning_data(MeanSpAmp(foi),dendriteaxis_bin);
    N_neuron=sum((cellfun(@sum,ind2)),2)';
    errorbar(dendBin_center, mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(b,:),'LineWidth',2); hold all
end
xlabel('Distance from Soma (\mum)')
ylabel('bAP Amplitude (spike height)')
legend(time_bin_strings)

%% Subthreshold dynamics
figure(110); clf;
tiledlayout(8,7);
for f=[27];
    area_threshold=150; perispike_time=[-5:20];
    EXC_kymo=Subthreshold{f}(Dist_order{f}(noi_dist{f}),:) > prctile(Subthreshold{f}(Dist_order{f}(noi_dist{f}),BlueStim{f}==0),90,2);
    INH_kymo=Subthreshold{f}(Dist_order{f}(noi_dist{f}),:)<prctile(Subthreshold{f}(Dist_order{f}(noi_dist{f}),BlueStim{f}>0),10,2);

    %derivSub{f}

    EXC_kymobw=bwlabel(EXC_kymo);
    excbwprops=struct2array(regionprops(EXC_kymobw,'Area'));
    INH_kymobw=bwlabel(INH_kymo);
    inhbwprops=struct2array(regionprops(INH_kymobw,'Area'));

    EXC_kymobw(ismember(EXC_kymobw,find(excbwprops<area_threshold)))=0;
    INH_kymobw(ismember(INH_kymobw,find(inhbwprops<area_threshold)))=0;

    perispike_frame=unique([tovec(find(double(allSpikeMat{f}(1,:)==1))'+perispike_time); find(CStrace{f})']);
    perispike_frame(perispike_frame<=0 | perispike_frame>size(NormalizedTrace_dirt{f},2))=[];
    nonvalid_frame=find(sum(isnan(NormalizedTrace_dirt{f}),1)>0);
    Blue_on_frame=find(imdilate(BlueStim{f}>0, [ones(1, 1), 1, ones(1, 200)]));
    Blue_off_frame=find(imdilate(BlueStim{f}==0, [1, ones(1, 50)]));

    Badframe=unique([perispike_frame' Blue_on_frame nonvalid_frame]);
    Badframe_inh=unique([perispike_frame' Blue_off_frame nonvalid_frame]);
    %Badframe=unique([Blue_on_frame nonvalid_frame]);
    %Badframe_inh=unique([Blue_off_frame nonvalid_frame]);
    Goodframe=setdiff([1:size(NormalizedTrace_dirt{f},2)],Badframe);
    Goodframe_inh=setdiff([1:size(NormalizedTrace_dirt{f},2)],Badframe_inh);

    EXC_kymobw(:,Badframe)=NaN;
    INH_kymobw(:,Badframe_inh)=NaN;

    ExcCandFrame=bwlabel(sum(EXC_kymobw,1)>0);
    InhCandFrame=bwlabel(sum(INH_kymobw,1)>0);

    ExcCandFrame = bwlabel(imdilate(ExcCandFrame>0, [ones(1,5) 1 ones(1,5)]));
    InhCandFrame = bwlabel(imdilate(InhCandFrame>0, [ones(1,5) 1 ones(1,5)]));

    ExcSub=Subthreshold{f}(Dist_order{f}(noi_dist{f}),:).*(ExcCandFrame>0);
    InhSub=Subthreshold{f}(Dist_order{f}(noi_dist{f}),:).*(InhCandFrame>0);

    ExcPatt=[];
    for b=1:max(ExcCandFrame)
        bseg=find(ExcCandFrame==b);
        [~, maxfrm]=max(mean(Subthreshold{f}(Dist_order{f}(noi_dist{f}),bseg),1,'omitnan'));
        if bseg(maxfrm)+3 < size(Subthreshold{f},2)
            ExcPatt=[ExcPatt mean(Subthreshold{f}(Dist_order{f}(noi_dist{f}),bseg(maxfrm)+[-3:3]),2,'omitnan')];
        end
    end
    ExcCandFrame(ismember(ExcCandFrame,find(sum(isnan(ExcPatt),1)>0)))=0;
    ExcCandFrame=bwlabel(ExcCandFrame>0);
    ExcPatt(:,sum(isnan(ExcPatt),1)>0)=[];

    InhPatt=[];
    for b=1:max(InhCandFrame)
        bseg=find(InhCandFrame==b);
        [~, minfrm]=min(mean(Subthreshold{f}(Dist_order{f}(noi_dist{f}),bseg),1,'omitnan'));
        InhPatt=[InhPatt mean(Subthreshold{f}(Dist_order{f}(noi_dist{f}),bseg(minfrm)+[-3:3]),2,'omitnan')];
    end
    InhCandFrame(ismember(InhCandFrame,find(sum(isnan(InhPatt),1)>0)))=0;
    InhCandFrame=bwlabel(InhCandFrame>0);
    InhPatt(:,sum(isnan(InhPatt),1)>0)=[];

    figure(3); clf;
    [ExcPatt_reduce NclusterExc Excindmat]= hierachyCluster_BH(ExcPatt,100,1);
    ExcPatt_reduce=ExcPatt_reduce-min(ExcPatt_reduce(:));

    figure(4); clf;
    [InhPatt_reduce NclusterInh Inhindmat]= hierachyCluster_BH(InhPatt,80,1);
    InhPatt_reduce=InhPatt_reduce-min(InhPatt_reduce(:));
    figure(113); clf;
    for n=1:NclusterExc %Silent during blue off
        nexttile([1 1])
        showim=max(double(Ftprnts{f}(:,:,Dist_order{f}(noi_dist{f}))>0).*reshape((rescale(ExcPatt_reduce(:,n))+0.1),1,1,[]),[],3);
        imshow2(showim,[0 1])
        title(['Exc. Patt. #' num2str(n), ', Fraction :', num2str(sum(Excindmat(:,n)/size(Excindmat,1)),2)])
    end
    colormap(turbo)

    figure(114); clf;
    for n=1:NclusterInh %Silent during blue off
        nexttile([1 1])
        showim=max(double(Ftprnts{f}(:,:,Dist_order{f}(noi_dist{f}))>0).*reshape(((InhPatt_reduce(:,n))),1,1,[]),[],3);
        imshow2(showim,[0 10])
        title(['Inh. Patt. #' num2str(n), ', Fraction :', num2str(sum(Inhindmat(:,n)/size(Inhindmat,1)),2)])
    end
    colormap(turbo)

    [matchingPatt, ~] = find(Inhindmat'>0);
    mathcingPattTr=zeros(1,size(InhCandFrame,2));
    for b=1:max(InhCandFrame)
        mathcingPattTr(find(InhCandFrame==b))=matchingPatt(b);
    end
    %ExcSub=Subthreshold{f}(Dist_order{f}(noi_dist{f}),:).*(EXC_kymobw>0);
    %InhSub=Subthreshold{f}(Dist_order{f}(noi_dist{f}),:).*(INH_kymobw>0);

    [ExcV ExcD ExcTrace]=get_eigvector(ExcSub(:,Goodframe)',20);
    [InhV InhD InhTrace]=get_eigvector(InhSub(:,Goodframe_inh)',20);
    [SubV SubD subTrace]=get_eigvector(Subthreshold{f}(Dist_order{f}(noi_dist{f}),sum(isnan(Subthreshold{f}),1)==0)',20);

    % [icsExc, ~, sepmat]=sorted_ica(ExcTrace(:,1:10),10);
    % ExcI=ExcSub(:,Goodframe)*icsExc;
    %
    % [icsInh, ~, sepmat]=sorted_ica(InhTrace(:,1:10),10);
    % InhI=InhSub(:,Goodframe_inh)*icsInh;

    figure(115); clf;
    for pc=1:7
        nexttile([1 1])
        showim=max(double(Ftprnts{f}(:,:,Dist_order{f}(noi_dist{f}))>0).*reshape((rescale(ExcV(:,pc))+0.1),1,1,[]),[],3);
        %showim=max(double(Ftprnts{f}(:,:,Dist_order{f}(noi_dist{f}))>0).*reshape((rescale(ExcI(:,pc))+0.1),1,1,[]),[],3);
        imshow2(showim,[])
        title(['Exc., ' num2str(ExcD(pc)/sum(ExcD)*100,2),' %'])
        %title(['Exc., ica comp. #' num2str(pc)])
    end

    for pc=1:7
        nexttile([1 1])
        showim=max(double(Ftprnts{f}(:,:,Dist_order{f}(noi_dist{f}))>0).*reshape((rescale(InhV(:,pc))+0.1),1,1,[]),[],3);
        %showim=max(double(Ftprnts{f}(:,:,Dist_order{f}(noi_dist{f}))>0).*reshape((rescale(InhI(:,pc))+0.1),1,1,[]),[],3);
        imshow2(showim,[])
        title(['Inh., ' num2str(InhD(pc)/sum(InhD)*100,2),' %'])
        %title(['Inh., ica comp. #' num2str(pc)])
    end
    colormap(Aurora)

end

%%
figure;
for f=foi
    nexttile([1 1])
    try
        imagesc(squeeze(mean(LapSub{f}(3:11,:,Dist_order{f}(noi_dist{f})),1,'omitnan'))',[-1 4]);
    catch
        imagesc(squeeze(mean(LapSub{f}(3:end,:,Dist_order{f}(noi_dist{f})),1,'omitnan'))',[-1 4])
    end
    nexttile([1 1])
    imagesc(LapFR{f});
    title(num2str(f))
end
colormap(turbo)

%%

f=26;
Lap_Nat_Pla=[4];
figure(120); clf; tiledlayout(4,1)
caxval=[-1 3];
avgLap=4;
nTauSilent=[-5:20];
%preLap=[Lap_Nat_Pla(1)-avgLap:Lap_Nat_Pla(1)-1];
preLap=[Lap_Nat_Pla(1)-1:Lap_Nat_Pla(1)-1];
postLap=[Lap_Nat_Pla(end)+1:Lap_Nat_Pla(end)+avgLap];
LasSub_prePla=squeeze(mean(ringmovMean(LapSub{f}(preLap,:,Dist_order{f}(noi_dist{f})),3),1,'omitnan'));
LasSub_postPla=squeeze(mean(ringmovMean(LapSub{f}(postLap,:,Dist_order{f}(noi_dist{f})),3),1,'omitnan'));
nexttile([1 1])
imagesc(LasSub_prePla',caxval); title('Average of 5 laps before plateau');
nexttile([1 1])
imagesc(LasSub_postPla',caxval); title('Average of 5 laps after plateau');
nexttile([1 1])
imagesc(LasSub_postPla'-LasSub_prePla',caxval); title('After plateau - before plateau');
colormap(turbo)
nexttile([1 1])
dFR=mean(LapFR{f}(postLap,:),1,'omitnan')-mean(LapFR{f}(preLap,:),1,'omitnan');
plot(dFR*1000); title('\Delta Firing rate (Post - pre)');

subthreshold_silent=Subthreshold{f};
subthreshold_silent(:,unique(get_perispikeIndex(allSpikeMat{f}(1,:),[-5:20])))=NaN;

%%
figure(121); clf;
subcell_boundary=[-40 40; -300 -60; 150 500];
for f=foi
    if ~isnan(PlaceFieldList{f})

        Dsign=ones(1,size(interDendDist{f},1));
        Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
        dendaxis=interDendDist{f}(1,:).*Dsign;

        ref_ROIdist=find(dendaxis>subcell_boundary(1,1) & dendaxis<subcell_boundary(1,2));
        bas_ROIdist=find(dendaxis>subcell_boundary(2,1) & dendaxis<subcell_boundary(2,2));
        api_ROIdist=find(dendaxis>subcell_boundary(3,1) & dendaxis<subcell_boundary(3,2));

        for pf=1:size(PlaceFieldList{f},2)/2
            nexttile([1 1]);
            pflap=PlaceFieldList{f}([2*pf-1 2*pf]);

            plot(squeeze(mean(ringmovMean(LapSubSilent{f}(pflap(1):pflap(2),:,ref_ROIdist),5),[1 3],'omitnan'))); hold all
            plot(squeeze(mean(ringmovMean(LapSubSilent{f}(pflap(1):pflap(2),:,bas_ROIdist),5),[1 3],'omitnan'))); hold all
            plot(squeeze(mean(ringmovMean(LapSubSilent{f}(pflap(1):pflap(2),:,api_ROIdist),5),[1 3],'omitnan'))); hold all
            grid on
            ylabel('Subthreshold')
            yyaxis right
            plot(mean(ringmovMean(LapFR{f}(pflap(1):pflap(2),:),3),1,'omitnan')*1000,'k')
            legend({'Peri-soma','Basal','Apical','Firing rate'})
            xlabel('VR position (bin)')
            ylabel('Firing rate')
            title(num2str(f))
        end

    end

end

PlaceFieldcat=[]; PlaceTimecat=[];
nBin=[-30:30]; nTauPF=[-5000:5000];
g=1; dendytick=[];
for f=[4 5 6 10 11 14 15 16 17 18 19 20 26 27]
    if ~isnan(PlaceFieldList{f})

        Dsign=ones(1,size(interDendDist{f},1));
        Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
        dendaxis=interDendDist{f}(1,:).*Dsign;
        dendaxis=dendaxis(Dist_order{f}(noi_dist{f}));


        for pf=1:size(PlaceFieldList{f},2)/2
            pflap=PlaceFieldList{f}([2*pf-1 2*pf]);
            pfbin=PlaceFieldBin{f}([2*pf-1 2*pf])+size(LapSub{f},2);
            if pfbin(1)>pfbin(2)
                pfbin(2)=pfbin(2)+size(LapSub{f},2);
            end
            pfpeakbin=mod(mean(pfbin),size(LapSub{f},2));

            repLap=squeeze(mean(repmat(ringmovMean(LapSub{f}(pflap(1):pflap(2),:,Dist_order{f}(noi_dist{f})),3),1,3),1,'omitnan'))';
            repLapFR=repmat(ringmovMean(LapFR{f},3),1,3);
            pftmp=repLap(:,round(mean(pfbin))+nBin);

            [PosTriggerTrace]=get_PositionAlignTrace(Subthreshold{f}(Dist_order{f}(noi_dist{f}),:),pfpeakbin/size(LapSub{f},2)*115,nTauPF,VRtrack{f});
            pttmp=squeeze(mean(movmean(PosTriggerTrace,30,2,'omitnan'),1,'omitnan'))';

            [posSpike]=get_PositionAlignTrace(allSpikeMat{f}(1,:),pfpeakbin/size(LapSub{f},2)*115,nTauPF,VRtrack{f});
            pttmp_sp=squeeze(sum(posSpike,1,'omitnan'))';

            PlaceFieldcat{g}=[[dendaxis' pftmp]; [NaN mean(repLapFR(pflap(1):pflap(2),round(mean(pfbin))+nBin),1,'omitnan')]];
            PlaceTimecat{g}=[[dendaxis' pttmp]; [NaN pttmp_sp']];

            if find(Dist_order{f}(noi_dist{f})==1)~=1
                dendytick{g}=[[1 find(Dist_order{f}(noi_dist{f})==1) sum(noi_dist{f})]' [min(interDendDist{f}(1,:).*Dsign) 0 max((interDendDist{f}(1,:).*Dsign))]'];
            else
                dendytick{g}=[[1 sum(noi_dist{f})]' [0 max((interDendDist{f}(1,:)))]'];
            end
            g=g+1;
        end
    end
end

figure(122); clf; tiledlayout(3,6);
for g=1:length(PlaceFieldcat)
nexttile([1 1])
imagesc(nBin*2/150,[1:size(PlaceFieldcat{g},1)-1],PlaceFieldcat{g}(1:end-1,2:end)); hold all
set(gca,'YTick',dendytick{g}(:,1),'YTicklabel',num2str(dendytick{g}(:,2),3))
yyaxis right
plot(nBin*2/150,PlaceFieldcat{g}(end,2:end)*1000,'k')
xlabel('Peri-place field (m)')
ylabel('Mean firing rate (Hz)')
end
colormap(turbo)


figure(123); clf; tiledlayout(9,6)
for g=1:length(PlaceTimecat)
    cax=[prctile(tovec(PlaceTimecat{g}(1:end-1,2:end)),0.5) prctile(tovec(PlaceTimecat{g}(1:end-1,2:end)),99.5)];
nexttile(mod(g,6)+6*(mod(g,6)==0)+18*floor(g/6.1),[2 1])
imagesc(nTauPF,[1:size(PlaceTimecat{g},1)-1],PlaceTimecat{g}(1:end-1,2:end)); hold all
set(gca,'YTick',dendytick{g}(:,1),'YTicklabel',num2str(dendytick{g}(:,2),3))
nexttile(mod(g,6)+6*(mod(g,6)==0)+18*floor(g/6.1)+12,[1 1])
plot(nTauPF,movmean(PlaceTimecat{g}(end,2:end),200)*1000,'k')
ylabel('Number of spike')
end
colormap(turbo)























