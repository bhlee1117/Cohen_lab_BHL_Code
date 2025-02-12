
clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'B5:P175');

fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);
StimROI=raw(:,6);
StimWfn=raw(:,7);
isGoodCell=cell2mat(raw(:,11));
PixelSize=cell2mat(cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false));
refROI=cellfun(@(x) (str2num(num2str(x))),raw(:,14),'UniformOutput',false);
maintrunkROI=cellfun(@(x) (str2num(num2str(x))),raw(:,15),'UniformOutput',false);
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
[~, unqInd] = unique([Mouse NeuronInd] ,'row');
set(0,'DefaultFigureWindowStyle','docked')


%%
% for f=[168]
% load(fullfile(fpath{f},'OP_Result.mat'),'Result');
%     interDendDist=[];
%     SkelDend = Skeletonize_dendrite(Result.ref_im,10,0.02,10);
%     nROI=size(Result.ftprnt,3);
%     for i=1%:nROI
%         i
%         for j=1:nROI
%             [interDendDist(i,j), path]=geodesic_distance(SkelDend,get_coord(Result.ftprnt(:,:,i)),get_coord(Result.ftprnt(:,:,j)));
%         end
%     end
%     Result.interDendDist=interDendDist;
%     save(fullfile(fpath{f},'OP_Result.mat'),'Result','-v7.3');
% end
%
% %%
% [~, unqInd] = unique([Mouse NeuronInd] ,'row');
% for i=unqInd([37])'
% i
%     cat_trace=[];
%     cat_spike=[];
%     SameCellInd=find(Mouse==Mouse(i) & NeuronInd==NeuronInd(i) & isGoodCell==1);
%     for j=SameCellInd'
% load(fullfile(fpath{j},'OP_Result.mat'),'Result');
%     cat_trace=[cat_trace Result.normTraces];
%     cat_spike=[cat_spike Result.spike];
%     end
%
%     cat_trace=cat_trace-movmedian(cat_trace,300,2);
%     [V D eigTrace]=get_eigvector(cat_trace',20);
%     ind2use=find(cumsum(D)/sum(D)>0.80,1)
%     F0_PCA=sqrt(sum((V(:,[1:ind2use]).*sqrt(D([1:ind2use]))').^2,2));
%
%     for j=SameCellInd'
%         j
%         load(fullfile(fpath{j},'OP_Result.mat'),'Result');
%         Result.F0_PCA=F0_PCA;
%         save(fullfile(fpath{j},'OP_Result.mat'),'Result','-v7.3');
%     end
% end
%% Figure 1. Representative image
f=125;
load(fullfile(fpath{f},'OP_Result.mat'),'Result');
mov_mc=readBinMov_BHL(fpath{f});
load(fullfile(fpath{f},'mcTrace01.mat'))
mov_res= mov_mc-mean(mov_mc,3);
mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,:));
mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,:).^2);
mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,1).*mcTrace.xymean(:,end));

STAmov=get_STA(tovec(mov_res),Result.SpClass(1,:),30,20);
STAmov=reshape(STAmov,size(Result.ref_im,1),size(Result.ref_im,2),[]);
STAmov=(-STAmov).*(max(Result.bvMask,[],3)==0);
STAmov=pcafilt(STAmov,5);
STAmov=STAmov-mean(STAmov(:,:,1:10),3);

SkelDend = Skeletonize_dendrite(Result.ref_im,7,0.02,10);
RobustdFF_const=get_robustdFF(STAmov,Result.ftprnt,(Result.ref_im-100).*(max(Result.bvMask,[],3)==0));
%normTrace = Result.traces_bvMask./RobustdFF_const';
STAtmp=get_STA(Result.normTraces,Result.spike(1,:),20,20);
STAtmp= STAtmp-median(STAtmp(:,1:10),2);
F_ref=mean(STAtmp(:,20+[9:10]),2);
Subthreshold=get_subthreshold(Result.normTraces,Result.spike(1,:),7,17);
F0_PCA=get_F0PCA(Subthreshold,3);

normTrace = Result.normTraces./F0_PCA;
%normTrace = Result.normTraces./F_ref;
normTrace = normTrace - prctile(normTrace,20,2);
normTrace = normTrace(Result.dist_order,:);
STAtr=get_STA(normTrace,Result.spike(1,:),20,20);
normTrace_filt = pcafilterTrace(normTrace, 3); %filter subthreshold
%normTrace_filt = normTrace; %filter subthreshold
%normTrace_filt_bd=get_bandstop(normTrace_filt,1000,[45 60]);

ftprnt=Result.ftprnt(:,:,Result.dist_order).*SkelDend;

%%
ROIs = {[2 5],[7],[11 14 16],[21],[29]}; %Basal, Soma, trunk, Oblique, tuft
ROIaxis = [2 5 7 11 14 20 26 29 30 31]; %Basal, Soma, trunk, Oblique, tuft
ROIaxis2 = [2 5 7 11 14 20 26 29];

% ROIs = {[6 9 10],[11],[15 17 20],[14 16],[22 24 27]}; %Basal, Soma, trunk, Oblique, tuft
% ROIaxis = [6 9 10 11 15 17 20 22 24 27]; %Basal, Soma, trunk, Oblique, tuft

Show_tr=[]; cmap=[0 0.5 1; 0 0 0; 1 0.6 0.3; 1 0.4 0.4 ; 1 0 0];
%cmap=distinguishable_colors(5);
Dsign=ones(1,size(Result.interDendDist,2));
Dsign(Result.dist_order(1:find(Result.dist_order==1)-1))=-1;
contourdist=Result.interDendDist(1,Result.dist_order).*Dsign(Result.dist_order)*PixelSize(f);
contourdist_show=contourdist(ROIaxis);
contourdist_show2=contourdist(ROIaxis2);

figure(1); clf; tiledlayout(6,1); scale=0.047;
ax1=nexttile([6 1]);
ftprnt_mask=[];
for r=[1:5]
    Show_tr(r,:)=mean(normTrace_filt(ROIs{r},:),1,'omitnan');
    Show_tr(r,:)=Show_tr(r,:)-prctile(Show_tr(r,Result.Blue==0),50);
    ftprnt_mask(:,:,r)=max(ftprnt(:,:,ROIs{r}),[],3);
    plot(Show_tr(r,:)-scale*r,'color',cmap(r,:)); hold all
end
axis off;

plot(Result.Blue/100-scale*6,'color',[0 0.6 1]); axis off;
%linkaxes([ax1 ax2],'x')
xlim([10 9999])

cmap_tr=hot(62);
cmap_tr=cmap_tr(1:31,:);
figure(2); clf; l=[];
tiledlayout(8,1);
ax1=nexttile([6 1]);
g=1;
% l=plot(normTrace_filt');
% arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap_tr,2)); hold all
for r=[2 5]
    l(g)=plot(Show_tr(r,:)-mean(Show_tr(r,3900:3950),2,'omitnan'),'color',cmap(r,:)); hold all
    g=g+1;
end
plot(Result.Blue/100-scale/2,'color',[0 0.6 1]);
%legend(l,{'Soma','Distal'})
axis tight off;

ax3=nexttile([2 1]);
[~, dind]=sort(contourdist_show2,'ascend');
[X, Y] = meshgrid([1:size(Result.normTraces,2)], min(contourdist_show2):40:max(contourdist_show2));
normTrace_showq = interp2([1:size(Result.normTraces,2)], contourdist_show2(dind), normTrace_filt(ROIaxis2(dind),:), X, Y, 'linear');
imagesc([1:size(Result.normTraces,2)],Y(:,1),normTrace_showq)
set(gca,'YDir','reverse');
caxis([-0.01 0.055])
colormap(turbo(256));

linkaxes([ax1 ax3],'x')
xlim([3950 4300])
set(gca,'XTick',[4000 4250],'XTickLabel',[0 250])
%set(gca,'YTick',[contourdist_show(1) 0 contourdist_show(end)],'YTickLabel',num2str([contourdist_show(1) 0 contourdist_show(end)]',3))
xlabel('Time (ms)')
ylabel(['Distance from' newline 'soma (\mum)'])

% figure(2); clf;
% show_footprnt_contour(ftprnt_mask,Result.Structure_bin,cmap)

figure(3); clf; tiledlayout(1,2);
h=nexttile([1 1]);
STAtr_show=STAtr-mean(STAtr(:,1:10),2);
l=plot([-20:20],STAtr_show');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(gen_colormap([0 0 0; 0.5 0 0; 1 0 0],31),2));
xlim([-4 10])
set(gca,'XTick',[-4 0 4 8])
xlabel('Time (ms)')
ylabel(['\DeltaF/F'])
cb = colorbar;
cb.Ticks = []; % Disable all ticks
cb.Label.String = 'Basal to distal';
colormap(h, gen_colormap([0 0 0; 0.5 0 0; 1 0 0]));

h2=nexttile([1 1]);
[~, dind]=sort(contourdist_show2,'ascend');
[X, Y] = meshgrid([-20:0.5:20], min(contourdist_show2):15:max(contourdist_show2));
STAtr_showq = interp2([-20:20], contourdist_show(dind), STAtr_show(ROIaxis(dind),:), X, Y, 'linear');
imagesc(X(1,:),Y(:,1),STAtr_showq,[0 0.025]);
colormap(h2,"turbo")
set(gca,'YDir','reverse')
xlim([-4 10])
xlabel('Time (ms)')
ylabel(['Distance from' newline 'soma (\mum)'])
cb = colorbar;
cb.Label.String = '\DeltaF/F';

%%
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
%% Figure 2
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
dendaxis=interDendDist(1,:).*Dsign;
dendaxis=dendaxis(Dist_order(noi_dist));

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

basalind=[8 9 10];
apicalind=[28 29 30];
somaind=[13 14 15];

NormalizedTrace=(Result.normTraces)./Result.F_ref;
NormalizedTrace_noNaN=NormalizedTrace;
%NormalizedTrace(:,Result.motionReject>0)=NaN;
% if ifdirtRemov(f)
%     NormalizedTrace(Result.dirtTrace>0)=NaN;
% end

bAP_STA=get_STA(NormalizedTrace, Result.spike(1,:).*double(Result.Blue==0), 30, 20);
bAP_STA=bAP_STA-prctile(bAP_STA,20,2);
SpikeHeight=max(mean(bAP_STA(perisomaROI,:),1),[],2);

NormalizedTrace=NormalizedTrace/SpikeHeight;
NormalizedTrace_noNaN=NormalizedTrace_noNaN/SpikeHeight;
NormalizedTrace=NormalizedTrace(Dist_order(noi_dist),:);
NormalizedTrace_noNaN=NormalizedTrace_noNaN(Dist_order(noi_dist),:);
filteredNormTr = pcafilterTrace(NormalizedTrace,5);
filteredNormTr2 = pcafilterTrace(NormalizedTrace_noNaN, 5);
Subthreshold=get_subthreshold(filteredNormTr2,max(Result.spike(1,:),[],1)>0,7,17);
%filteredNormTr = NormalizedTrace;
ftprnt=Result.ftprnt(:,:,Dist_order(noi_dist));

Show_tr=[mean(filteredNormTr2(basalind,:),1,'omitnan'); mean(filteredNormTr2(somaind,:),1,'omitnan'); mean(filteredNormTr2(apicalind,:),1,'omitnan')];
Show_tr=interpolateNaN(Show_tr);
Show_subtr=get_subthreshold(Show_tr,max(Result.spike,[],1)>0,7,17);
Show_ftprnt=cat(3,max(ftprnt(:,:,basalind)>0,[],3),max(ftprnt(:,:,somaind)>0,[],3),max(ftprnt(:,:,apicalind)>0,[],3));
%roi_show=setdiff([1:sum(noi_dist)],[2 38 21]);
roi_show=setdiff([1:sum(noi_dist)],[17]);

figure(4); clf;
show_footprnt_contour(Show_ftprnt,Result.ref_im,[0 0.4 1; 1 0.7 0; 1 0 0])

figure(5); clf; tiledlayout(4,1);
tax=[1:size(Result.normTraces,2)]/1000;
% ax1=nexttile([2 1]);
% imagesc(tax,dendaxis(roi_show),filteredNormTr(roi_show,:),[-0.2 1.5]);
% colormap(ax1,turbo(256));

ax3=nexttile([3 1]);
l=plot(tax,Show_tr'-[1:3]*2);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell([0 0.4 1; 1 0.7 0; 1 0 0],2));
hold all;
l2=plot(tax,Show_subtr'-[1:3]*2);
arrayfun(@(l,c) set(l,'Color',c{:}),l2,num2cell([0 0.4 1; 1 0.7 0; 1 0 0]/2,2));
legend(l,{'Basal','Soma','Distal'})
axis off

ax2=nexttile([1 1]);
plot(tax,Result.VR(5,:),'color',[0.1 0.9 0.2],'linewidth',2);
axis off
linkaxes([ax2 ax3],'x')
xlim([45 160])

figure(6); clf; tiledlayout(5,1); %  zoom in version of figure(5);
tax=[1:size(Result.normTraces,2)]/1000;
ax1=nexttile([2 1]);
imagesc(tax,dendaxis(roi_show),filteredNormTr(roi_show,:),[-0.2 1]);
%imagesc(NormalizedTrace(roi_show,:),[-0.2 1.5]);
colormap("turbo")
set(gca,'XTick',[77 78 79],'XTickLabel',[0 1 2])

ax3=nexttile([3 1]);
l=plot(tax,Show_tr'-[1:3]*1.5);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell([0 0.4 1; 1 0.7 0; 1 0 0],2));
hold all;
l=plot(tax,Show_subtr'-[1:3]*1.5);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell([0 0.4 1; 1 0.7 0; 1 0 0]/2,2));
axis off
linkaxes([ax1 ax3],'x')
xlim([76.9 79.1])

%% Figure 3
f=23; 
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
dendaxis=interDendDist(1,:).*Dsign;
dendaxis=dendaxis(Dist_order(noi_dist));

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
    BS_trace(1,[ss_time(bwn): ss_time(bwn(end)+1)])=b;
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
figure(101); clf; tiledlayout(3,3)
stype_str={'SS','CS','BS'}; ax1=[]; ax2=[]; FR_bin_Frame=12; ax3=[];
nTau=[200 100]; cmap=distinguishable_colors(3);
NormalizedTraceFiltered=pcafilterTrace(NormalizedTrace,5);
%NormalizedTraceFiltered=NormalizedTrace;
Subthreshold=get_subthreshold(NormalizedTrace,Result.spike(1,:),7,17); cax=[-0.05 1];
basalind=[1 2 3 4 5 6]; apicalind=[36 37 38 39 40];
Basal_Trace= mean(NormalizedTraceFiltered(basalind,:),1,'omitnan');
Apical_Trace= mean(NormalizedTraceFiltered(apicalind,:),1,'omitnan');

silenttime=setdiff([1:nTime],unique(find(max([allSpikeMat{f}(1,:); allSpikeClassMat{f}; CStrace{f}; BlueStim{f}])>0)'+[-5:30]));
silenttime_vec=ind2vec(nTime,silenttime,1);

BASub =[Basal_Trace; Apical_Trace];
BASub_silent=[Basal_Trace; Apical_Trace];
BASub_silent=get_subthreshold(BASub_silent,Result.spike(1,:),7,17);
BASub_silent(:,silenttime_vec==0)=NaN;

STApreSpikeBA=[];
for sclass=1:3
    TriggerSpike=SpClass(sclass,:).*double(Result.Blue==0);
    [STATrace STATraceMat]=get_STA(NormalizedTrace,TriggerSpike,nTau(1),nTau(2));
    STASpikeClass=get_STA(SpikeClassvec,TriggerSpike,nTau(1),nTau(2));
    %STASpikeClass(:,nTau(1):end)=0;
    STAsub=get_STA(Subthreshold,TriggerSpike,nTau(1),nTau(2));
    [~, STApreSpikeBA{sclass}] = get_STA(BASub,TriggerSpike,nTau(1),nTau(2));

    [STASpikeClass_bin binT]= Bin_Vector(STASpikeClass, [-nTau(1):nTau(2)], [-nTau(1):FR_bin_Frame:-1 0:1 2:FR_bin_Frame:nTau(2)]);
    %STASpikeClass_bin=STASpikeClass;
    [~, dind]=sort(dendaxis,'ascend');
    [X, Y] = meshgrid([-nTau(1):nTau(2)], min(dendaxis):36:max(dendaxis));
    STAtraceq = interp2([-nTau(1):nTau(2)], dendaxis(dind), STATrace(dind,:), X, Y, 'linear');
    STAsubq = interp2([-nTau(1):nTau(2)], dendaxis(dind), STAsub(dind,:), X, Y, 'linear');

    ax1=[ax1 nexttile(sclass,[1 1])];
   imagesc([-nTau(1):nTau(2)],Y(:,1),STAtraceq,cax) 
   ylabel('Distance from soma (\mum)')
colormap(turbo(256))

    ax3=[ax3 nexttile(sclass+3,[1 1])];
    STAtraceMatB=mean(STATraceMat(basalind,:,:),1,'omitnan'); STAtraceMatA=mean(STATraceMat(apicalind,:,:),1,'omitnan');
    BasalM=squeeze(mean(STAtraceMatB,[2],'omitnan')); BasalS=squeeze(std(STAtraceMatB,0,2,'omitnan'));
    ApicalM=squeeze(mean(STAtraceMatA,2,'omitnan')); ApicalS=squeeze(std(STAtraceMatA,0,2,'omitnan'));
    h(1)=errorbar_shade([-nTau(1):nTau(2)],BasalM,BasalS./sqrt(size(STATraceMat,2)),[0.1 0.5 1]); hold all
    h(2)=errorbar_shade([-nTau(1):nTau(2)],ApicalM,ApicalS./sqrt(size(STATraceMat,2)),[1 0 0]); hold all

ax2=[ax2 nexttile(sclass+6,[1 1])];
l=plot(binT,STASpikeClass_bin');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap,2));
xlabel('Peri-spike time (ms)')
ylabel(['Spike rate (per ' stype_str{sclass} ')'])
ylim([0 0.02])
end
legend(l,stype_str)
ax1=[ax1 nexttile(3,[1 1])];
colorbar;
nexttile(6,[1 1]);
legend(h,{'Basal','Distal'})
linkaxes(ax1,'xy');
linkaxes(ax2,'xy');
linkaxes(ax3,'xy');

%[V D]=get_eigvector(BASub_silent(:,sum(isnan(BASub_silent))==0),2);
figure(302); clf; %tiledlayout(2,2);
ax1=nexttile([1 1]); l=[];
marker_class={'>','o','diamond'};
[~,~,~,h2]=scatter_heatmap2(BASub_silent(1,:),BASub_silent(2,:),linspace(-1,1,50),linspace(-1.5,1.5,50)); hold all
colormap(h2.Parent,bone);
caxis([0.0001 0.005])
colormap_list={'spring','hot','winter'};
for sclass=1:3
    preSub=mean(STApreSpikeBA{sclass}(:,:,nTau(1)+[-6:0]),3,'omitnan');
    %l(sclass)=scatter(preSub(1,:),preSub(2,:),40,'filled','Marker',marker_class{sclass},'MarkerFaceColor',cmap(sclass,:),'MarkerFaceAlpha',0.7); hold all
    [count xc yc h_heatmap]=scatter_heatmap2(preSub(1,:),preSub(2,:),linspace(-1,1,10),linspace(-1.5,1.5,10));
    delete(h_heatmap)
    [~, h{sclass}]=contour(xc,yc,count,10); hold all
    colormap(h{sclass}.Parent,colormap_list{sclass});
end
plot([-2 3],[0 0],'color',[0.7 0.7 0.7]); plot([0 0],[-2 3],'color',[0.7 0.7 0.7]);
legend(l,{'SS','CS','BS'})

ax3=nexttile([1 1]); l=[];
scatter_heatmap2(BASub_silent(1,:),BASub_silent(2,:),linspace(-2,3,100),linspace(-2,3,100)); hold all
plot([-2 3],[0 0],'color',[0.7 0.7 0.7]); plot([0 0],[-2 3],'color',[0.7 0.7 0.7]);

quiver(mean(BASub_silent(1,:),'omitnan'), mean(BASub_silent(2,:),'omitnan'), V(1,1)*sqrt(D(1)), V(2,1)*sqrt(D(1)), 'r', 'LineWidth', 1.5, 'DisplayName', 'Eigenvector 1');
quiver(mean(BASub_silent(1,:),'omitnan'), mean(BASub_silent(2,:),'omitnan'), V(1,2)*sqrt(D(2)), V(2,2)*sqrt(D(2)), 'g', 'LineWidth', 1.5, 'DisplayName', 'Eigenvector 2');
colormap(turbo); xlabel('Basal'); ylabel('Apical'); legend(l,{'pre CS','pre SS'})
ax4=nexttile([1 1]); l=[];
scatter_heatmap2(BASub_silent(1,:),BASub_silent(2,:),linspace(-2,3,100),linspace(-2,3,100)); hold all
plot([-2 3],[0 0],'color',[0.7 0.7 0.7]); plot([0 0],[-2 3],'color',[0.7 0.7 0.7]);
l(1)=scatter(BApostCS1st(1,:),BApostCS1st(2,:),40,'filled','Marker','>','MarkerFaceColor',[1 0 0],'MarkerFaceAlpha',0.7); hold all
l(2)=scatter(BApostSS(1,:),BApostSS(2,:),40,'filled','Marker','o','MarkerFaceColor',[0 0 0],'MarkerFaceAlpha',0.7);
quiver(mean(BASub_silent(1,:),'omitnan'), mean(BASub_silent(2,:),'omitnan'), V(1,1)*sqrt(D(1)), V(2,1)*sqrt(D(1)), 'r', 'LineWidth', 1.5, 'DisplayName', 'Eigenvector 1');
quiver(mean(BASub_silent(1,:),'omitnan'), mean(BASub_silent(2,:),'omitnan'), V(1,2)*sqrt(D(2)), V(2,2)*sqrt(D(2)), 'g', 'LineWidth', 1.5, 'DisplayName', 'Eigenvector 2');
colormap(turbo); xlabel('Basal'); ylabel('Apical'); legend(l,{'3-8 ms after 1st CS','3-8 ms after 1st SS'})
linkaxes([ax3 ax4],'xy')







