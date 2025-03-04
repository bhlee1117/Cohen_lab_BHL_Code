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
clear SpikeInd SpikeMat STAmat NormalizedTrace_ch NormalizedTrace_dirt SpikeIndBlueOff Dist_order allSpikeMat noi interDendDist noi_dist derivSub
clear Subthreshold

f=22
   
    load(fullfile(fpath{f},'PC_Result.mat'),'Result')
    rois={basal_ROI{f},apical_ROI{f},ref_ROI{f}};
    for r=1:length(rois)
        if isnan(rois{r}); rois{r}=[]; end;
    end
    Dist_order{f}=Result.dist_order;
    nROI=size(Result.normTraces,1);
    nTau_bAP=[-60:50];

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

    STA_SS=squeeze(mean(reshape(Result.normTraces(:,bAP_ref'+nTau_bAP),nROI,[],length(nTau_bAP)),2,'omitnan'));
    STA_SS= STA_SS - prctile(STA_SS,10,2);
    %STA_SS= STA_SS - prctile(STA_SS(:,1:30),30,2);
    %STA_SS= STA_SS - mean(STA_SS(:,1:37),2,'omitnan');

    %F_ref=max(STA_SS,[],2)./Result.dFF_slope';

    %plot(mean(STA_SS(:,-nTau_bAP(1)+[7:11]),2),max(STA_SS,[],2)./Result.dFF_slope','.','MarkerSize',20); hold all
    F_ref=mean(STA_SS(:,-nTau_bAP(1)+[8:10]),2);

    clf; tiledlayout(4,1);
nexttile([1 1])
plot(STA_SS(Result.dist_order,:)'); hold all
%imagesc(rescale2(STA_SS(Result.dist_order,:),2)); hold all
ylim([-mean(F_ref) mean(F_ref)*3]); plot(-nTau_bAP(1)+[7 10],[0 0],'r');
nexttile([1 1])
plot(F_ref(Result.dist_order),max(STA_SS(Result.dist_order,:),[],2),'.');
nexttile([1 1])
imagesc(STA_SS(Result.dist_order,:)./F_ref(Result.dist_order))
colormap(turbo)
nexttile([1 1])
plot(Result.interDendDist(1,:),max(STA_SS./Result.F_ref,[],2),'.')
%%

Result.F_ref=F_ref;
save(fullfile(fpath{f},'PC_Result.mat'),'Result','-v7.3')
load(fullfile(fpath{f},'PC_Result.mat'),'Result')
backupServer(fpath{f},'BHL18TB_D2','cohen_lab/Lab/Labmembers/Byung Hun Lee/Data')