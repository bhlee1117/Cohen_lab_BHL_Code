
clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'C5:Q175');

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
%% Figure 1. Representative image
f=125;
load(fullfile(fpath{f},'OP_Result.mat'),'Result');
load(fullfile(fpath{f},'RobustdFFfit.mat'));
% mov_mc=readBinMov_BHL(fpath{f});
% load(fullfile(fpath{f},'mcTrace01.mat'))
% mov_res= mov_mc-mean(mov_mc,3);
% mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,:));
% mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,:).^2);
% mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,1).*mcTrace.xymean(:,end));
% 
% STAmov=get_STA(tovec(mov_res),Result.SpClass(1,:),30,20);
% STAmov=reshape(STAmov,size(Result.ref_im,1),size(Result.ref_im,2),[]);
% STAmov=(-STAmov).*(max(Result.bvMask,[],3)==0);
% STAmov=pcafilt(STAmov,5);
% STAmov=STAmov-mean(STAmov(:,:,1:10),3);

SkelDend = Skeletonize_dendrite(Result.ref_im,7,0.02,10,0);
%RobustdFF_const=get_robustdFF(STAmov,Result.ftprnt,(Result.ref_im-100).*(max(Result.bvMask,[],3)==0));
%normTrace = Result.traces_bvMask./RobustdFF_const';
STAtmp=get_STA(Result.normTraces,Result.spike(1,:),20,20);
STAtmp= STAtmp-median(STAtmp(:,1:7),2);
F_ref=mean(STAtmp(:,20+[9:10]),2);
Subthreshold=get_subthreshold(Result.normTraces,Result.spike(1,:),7,17);
%F0_PCA=get_F0PCA(Subthreshold,3);

ROIvec=Result.maintrunk;
normTrace = Result.normTraces./RobustdFFfit.ScaleFactor;
%normTrace = Result.normTraces./F_ref;
normTrace = normTrace - prctile(normTrace,20,2);
normTrace_main= normTrace(ROIvec,:);
normTrace = normTrace(Result.dist_order,:);
STAtr=get_STA(normTrace,Result.spike(1,:),20,20);
STAtr_main=get_STA(normTrace_main,Result.spike(1,:),20,20);
%normTrace_filt = normTrace; %filter subthreshold
normTrace_filt = pcafilterTrace(normTrace, [1:10]); %filter subthreshold
normTrace_filt_main=normTrace_filt(find_index_bh(Result.dist_order,ROIvec),:);
ftprnt=Result.ftprnt(:,:,Result.dist_order).*SkelDend;

ROIs = {[2 5],[7],[11 14 16],[21],[29]}; %Basal, Soma, trunk, Oblique, tuft
ROIaxis = [2 5 7 11 14 20 26 29 30 31]; %Basal, Soma, trunk, Oblique, tuft
ROIaxis2 = [2 5 7 11 14 20 26 29];

% ROIs = {[6 9 10],[11],[15 17 20],[14 16],[22 24 27]}; %Basal, Soma, trunk, Oblique, tuft
% ROIaxis = [6 9 10 11 15 17 20 22 24 27]; %Basal, Soma, trunk, Oblique, tuft
Show_tr=[]; cmap=gen_colormap(Plasma,6);
%cmap=distinguishable_colors(5);
Dsign=ones(1,size(Result.interDendDist,2));
Dsign(Result.dist_order(1:find(Result.dist_order==1)-1))=-1;
contourdist=Result.interDendDist(1,ROIvec).*Dsign(ROIvec)*PixelSize(f);
%contourdist_show2=contourdist(ROIaxis2);

figure(1); clf; tiledlayout(6,1); scale=0.2;
ax1=nexttile([6 1]);
ftprnt_mask=[];
for r=[1:5]
    Show_tr(r,:)=mean(normTrace_filt(ROIs{r},:),1,'omitnan');
    Show_tr(r,:)=Show_tr(r,:)-prctile(Show_tr(r,Result.Blue==0),50);
    ftprnt_mask(:,:,r)=max(ftprnt(:,:,ROIs{r}),[],3);
    plot(Show_tr(r,:)-scale*r,'color',cmap(r,:),'linewidth',1); hold all
end
axis off;
drawScaleBar(0.1,'vertical','position',[10050 -0.93],'color','k','linewidth',2.5);
drawScaleBar(1000,'horizontal','position',[10050 -0.93],'color','k','linewidth',2.5);
plot(Result.Blue/20-scale*5.7,'color',[0 0.6 1]); axis off;
drawScaleBar(0.0578,'vertical','position',[10050 -1.23],'color',[0 0.6 1],'linewidth',2.5);
%linkaxes([ax1 ax2],'x')
xlim([10 10055])
set_figsize(200,150);
%% Zoom in trace
figure(2); clf; l=[];
tiledlayout(8,1,'Padding','tight');
ax1=nexttile([5 1]);
g=1;
% l=plot(normTrace_filt');
% arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap_tr,2)); hold all
for r=[2 5]
    l(g)=plot(Show_tr(r,:)-mean(Show_tr(r,3900:3950),2,'omitnan'),'color',cmap(r,:),'linewidth',1.5); hold all
    g=g+1;
end
plot(Result.Blue/30-scale/2,'color',[0 0.6 1],'linewidth',1.5);
axis tight off;
drawScaleBar(0.1,'vertical','position',[4000 0],'color','k','linewidth',2.5);

ax3=nexttile([3 1]);
[~, dind]=sort(contourdist,'ascend');
[X, Y] = meshgrid([1:size(normTrace_main,2)], min(contourdist):35:max(contourdist));
normTrace_showq = interp2([1:size(normTrace_main,2)], contourdist, normTrace_filt_main, X, Y, 'linear');
imagesc([1:size(normTrace_main,2)],Y(:,1),normTrace_showq)
%imagesc([1:size(normTrace_main,2)],contourdist,normTrace_filt_main)
cb=colorbar;
cb.Label.String = '∆F/F'; cb.LineWidth=1;
set(gca,'YDir','reverse');
caxis([-0.005 0.2])
colormap(turbo(256));

linkaxes([ax1 ax3],'x')
xlim([3950 4300])
set(gca,'XTick',[4000 4250],'XTickLabel',[0 250])
%set(gca,'YTick',[contourdist_show(1) 0 contourdist_show(end)],'YTickLabel',num2str([contourdist_show(1) 0 contourdist_show(end)]',3))
xlabel('Time (ms)')
ylabel(['Distance from' newline 'soma (µm)'])

set_font('Arial'); set_fontsize(18); 
set(gca,'LineWidth',1);
set_figsize(270,170);
% figure(2); clf;
% show_footprnt_contour(ftprnt_mask,Result.Structure_bin,cmap)
%% Kymo figure (Figure 1G)
figure(3); clf; tiledlayout(1,2);
h=nexttile([1 1]);
STAtr_show=STAtr-mean(STAtr(:,1:20),2);
STAtr_show_main=STAtr_main-mean(STAtr_main(:,1:20),2);
l=plot([-20:20],STAtr_show','linewidth',1);
set(gca, 'LineWidth', 2)
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(gen_colormap(Plasma,31),2));
xlim([-4 10])
set(gca,'XTick',[-4 0 4 8])
xlabel('Time (ms)')
ylabel(['∆F/F']); box off;
cb = colorbar;
cb.Ticks = []; % Disable all ticks
cb.Label.String = 'Basal to distal'; cb.Label.FontSize=16;
colormap(h, gen_colormap(Plasma));

h2=nexttile([1 1]);
[~, dind]=sort(contourdist,'ascend');
[X, Y] = meshgrid([-20:0.5:20], min(contourdist):15:max(contourdist));
STAtr_showq = interp2([-20:20], contourdist, STAtr_show_main, X, Y, 'linear');
imagesc(X(1,:),Y(:,1),STAtr_showq,[-0.005 0.1]);
colormap(h2,turbo(256))
set(gca,'YDir','reverse')
xlim([-4 10])
xlabel('Time (ms)')
ylabel(['Distance from soma (µm)'])
cb = colorbar;
cb.Label.String = '∆F/F';
set_font('Arial'); set_fontsize(22); set_figsize(330,150);
%% SNAPT

[superlocColormov dtimg dFF]=generate_SNAPTmov(-STAmov(:,:,26:36),Result.Structure,[],[]);
%[yR xR zR]=size(Result.Structure);
    %bluePatt = bwboundaries(imwarp(Result.BlueDMDimg,Result.tform,'OutputView', imref2d([yR xR])));
    bluePatt= bwboundaries(Result.BlueDMDimg);
    figure(20); clf;
    v = VideoWriter([fpath{f} '/SNAPT_movie'],'MPEG-4');
    %v = VideoWriter([fpath{f} '/SNAPT_movie'],'Uncompressed AVI');

    open(v);
    subframeT = 0.025; % ms
    initialT = -2; % ms
    finalT = 2; % ms
    times = initialT:subframeT:finalT;
    SNAPTmov=superlocColormov;%.*Result.Structure;
    SNAPTmov=flipud(fliplr(SNAPTmov));
    SNAPTmov=SNAPTmov(:,125:640,:,:);

    for j = 1:length(times)
        clf;
        %set(gca,'units','pixels','position',[200 0 1000 800])
        imshow(SNAPTmov(:,:,:,j),[])
        pbaspect([size(double(SNAPTmov(:,:,:,j)),2) size(double(SNAPTmov(:,:,:,j)),1) 1]),colormap(gray)
        hold all
        %plot(bluePatt{1}(:,2),bluePatt{1}(:,1),'color',[0 0.6 1],'linewidth',2)
        axis off
        text(5,17,[num2str(times(j)+0.75,2) 'ms'], 'FontSize', 20, 'color', [0.99 0.99 0.99])% the value 1. is to adjust timing by eyes
        pause(0.1)
        set(gcf,'color','w')    % Sets background to white
        frame = getframe(gcf);
        writeVideo(v,frame);
        pause(0.1);
    end;
    close(v);

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
f=20; cmap=gen_colormap(Plasma,3);
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

NormalizedTrace=(Result.normTraces)./Result.F0_PCA;
NormalizedTrace_noNaN=NormalizedTrace;
%NormalizedTrace(:,Result.motionReject>0)=NaN;
% if ifdirtRemov(f)
%     NormalizedTrace(Result.dirtTrace>0)=NaN;
% end

bAP_STA=get_STA(NormalizedTrace, Result.spike(1,:).*double(Result.Blue==0), 30, 20);
bAP_STA=bAP_STA-prctile(bAP_STA,20,2);
SpikeHeight=max(mean(bAP_STA(perisomaROI,:),1),[],2);

%NormalizedTrace=NormalizedTrace/SpikeHeight;
%NormalizedTrace_noNaN=NormalizedTrace_noNaN/SpikeHeight;
NormalizedTrace=NormalizedTrace(Dist_order(noi_dist),:);
NormalizedTrace_noNaN=NormalizedTrace_noNaN(Dist_order(noi_dist),:);
filteredNormTr = pcafilterTrace(NormalizedTrace,[1:5]);
filteredNormTr2 = pcafilterTrace(NormalizedTrace_noNaN, [1:5]);
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
show_footprnt_contour(Show_ftprnt,Result.ref_im,cmap)

figure(5); clf; tiledlayout(4,1);
tax=[1:size(Result.normTraces,2)]/1000; lscale=[0.02];
% ax1=nexttile([2 1]);
% imagesc(tax,dendaxis(roi_show),filteredNormTr(roi_show,:),[-0.2 1.5]);
% colormap(ax1,turbo(256));

ax3=nexttile([3 1]);
l=plot(tax(45000:160000),Show_tr(:,45000:160000)'-[1:3]*lscale);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap,2));
hold all;
l2=plot(tax(45000:160000),Show_subtr(:,45000:160000)'-[1:3]*lscale);
arrayfun(@(l,c) set(l,'Color',c{:}),l2,num2cell(cmap/2,2));
legend(l,{'Basal','Soma','Distal'})
axis off

ax2=nexttile([1 1]);
plot(tax(45000:160000),Result.VR(5,45000:160000),'color',[0.1 0.9 0.2],'linewidth',2);
axis off
linkaxes([ax2 ax3],'x')
xlim([45 160])

figure(6); clf; tiledlayout(5,1); %  zoom in version of figure(5);
tax=[1:size(Result.normTraces,2)]/1000;
ax1=nexttile([2 1]);
imagesc(tax,dendaxis(roi_show),filteredNormTr(roi_show,:),[-0.005 0.017]);
%imagesc(NormalizedTrace(roi_show,:),[-0.2 1.5]);
colormap("turbo")
set(gca,'XTick',[77 78 79],'XTickLabel',[0 1 2])

ax3=nexttile([3 1]);
l=plot(tax,Show_tr'-[1:3]*0.017);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap,2));
hold all;
l=plot(tax,Show_subtr'-[1:3]*0.017);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap/2,2));
axis off
linkaxes([ax1 ax3],'x')
xlim([76.9 79.1])







