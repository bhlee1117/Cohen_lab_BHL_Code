clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/FromBackup/PP72_PlaceCellResults';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:P23');

% [~, ~, NeuronsToUse]=xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
%     'PlaceCellData_Arrangement.xlsx'], 'Sheet1', 'L8:M46');
%
% NeuronsToUse=cellfun(@(x) (str2num(num2str(x))),NeuronsToUse,'UniformOutput',false);
ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,10),'UniformOutput',false);
basal_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,11),'UniformOutput',false);
apical_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false);
fpath=raw(:,1)';
StructureData=raw(:,10);
EndFrame=cell2mat(raw(:,13));
ifmotionReject=cell2mat(raw(:,14));
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
save_at='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Updates/2024/20240408_Movs_Figs/';
%%
f=18; load(fullfile(fpath{f},'PC_Result.mat'),'Result')
nROI=size(Result.normTraces,1);
nTau={[-100:50],[-100:200],[-20:20]}; %SS, CS, dSP
nTau_bAP=[-20:20];
% Isolated Somatic spike
som_spike=find(Result.spike(1,:));
tr_ref=Result.normTraces(ref_ROI{f},:);
tr_sub=mean(tr_ref,1)-movprc(mean(tr_ref,1),200,20,2);
tr_sub=get_subthreshold(tr_sub,Result.spike(1,:),5,10);
[trans tr_trace]=detect_transient2(tr_sub,[5 1.5],Result.spike(1,:),15);

% Calculate distance order
SkelDend = Skeletonize_dendrite(Result.ref_im,4,0.02,25);
interDendDist=[];
for i=1:size(Result.normTraces,1)
    for j=1:size(Result.normTraces,1)
        [interDendDist(i,j), ~]=geodesic_distance(SkelDend,get_coord(Result.ftprnt(:,:,i)),get_coord(Result.ftprnt(:,:,j)));
    end
end
coord_1d=dim_reduce(get_coord(Result.ftprnt));
[~, dist_order]=sort(coord_1d,'ascend');
som_roi=find(dist_order==1);
geodist=interDendDist(1,:)'.*sign(coord_1d-coord_1d(1));
show_footprnt_contour(Result.ftprnt(:,:,dist_order),Result.ref_im)

% Isolated Somatic spike
bAP_s=[];
for s=som_spike
    isnearby=sum(ismember(s+nTau{1},som_spike))>1;
    isnearbyCS=sum(ismember(s+nTau{1},find(Result.CStrace)))>1;
    ispartCS=tr_trace(s)>0;
    if ~isnearby & ~isnan(s) & ~isnearbyCS & ~ispartCS
        bAP_s=[bAP_s s];
    end
end

% Reference bAPs for F_ref
bAP_ref=[];
for s=som_spike
    isnearby=sum(ismember(s+nTau_bAP,som_spike))>1;
    isnearbyCS=sum(ismember(s+nTau_bAP,find(Result.CStrace)))>1;
    ispartCS=tr_trace(s)>0;
    if ~isnearby & ~isnan(s) & ~isnearbyCS & ~ispartCS
        bAP_ref=[bAP_ref s];
    end
end

% Isolated Complex spike
CS_s=[];
C_list=find(Result.SpClass(2,:));
CS_label=max(Result.spike,[],1).*bwlabel(Result.CStrace);
CS_list=[];
for g=1:length(C_list)
    s=C_list(g);
    s_tmp=Result.spike(1,:);
    s_tmp(find(CS_label==g))=0;
    isnearbyCS=max(bwlabel(Result.CStrace(s+nTau{2})))>1;
    isnearbyS=sum(ismember(s+nTau{2},find(s_tmp)))>0;
    if ~isnearbyCS & ~isnearbyS
        CS_s=[CS_s s];
        CS_list=[CS_list g];
    end
end

prc_normTr=Result.normTraces;
%prc_normTr=Result.normTraces-movprc(Result.normTraces,500,35,2);
STA_SSmat=reshape(prc_normTr(:,bAP_s'+nTau{1}),nROI,[],length(nTau{1}));
STA_SS=squeeze(mean(reshape(prc_normTr(:,bAP_ref'+nTau_bAP),nROI,[],length(nTau_bAP)),2));
F_ref=mean(STA_SS(:,-nTau_bAP(1)+[10:14]),2);
STA_SSmat=STA_SSmat./F_ref;
prc_normTrCS=Result.normTraces;
%prc_normTrCS=Result.normTraces-movprc(Result.normTraces,500,35,2);
STA_CSmat=reshape(prc_normTrCS(:,CS_s'+nTau{2}),nROI,[],length(nTau{2}));
STA_CSmat=STA_CSmat./F_ref;
CSpike=Result.spike.*Result.CStrace;
STA_CSpikemat=reshape(CSpike(:,CS_s'+nTau{2}),nROI,[],length(nTau{2}));

figure(2); clf;
cax=[-0.5 1.5];
tiledlayout(2,2)
nexttile([1 1])
imagesc(squeeze(mean(STA_SSmat(dist_order,:,:),2)),cax); hold all
plot(1,som_roi,'marker','>','MarkerFaceColor','r')
ylabel('ROI (distance order)')
title('STA of SS')
nexttile([1 1])
imagesc(squeeze(mean(STA_CSmat(dist_order,:,:),2)),cax); hold all
plot(1,som_roi,'marker','>','MarkerFaceColor','r')
ylabel('ROI (distance order)')
title('STA of CS')
colormap(turbo)
nexttile([1 1])
l=plot(nTau{1},squeeze(mean(STA_SSmat(dist_order,:,:),2))');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2))
xlabel('Peri-spike time (ms)')
nexttile([1 1])
l=plot(nTau{2},squeeze(mean(STA_CSmat(dist_order,:,:),2))');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2))
xlabel('Peri-spike time (ms)')
for n=1:nROI
    roi_l{n}=['ROI#' num2str(n)];
end
legend(roi_l,'FontSize',8)

figure(3); clf; cmap=distinguishable_colors(6); ax1=[];
rois={basal_ROI{f},apical_ROI{f}};
tiledlayout(2,2)
ax1=[ax1 nexttile([1 1])];
errorbar_shade(nTau{1},squeeze(mean(mean(STA_SSmat(rois{1},:,:),1),2))',squeeze(std(mean(STA_SSmat(rois{1},:,:),1),0,2))',cmap(1,:)); hold all
errorbar_shade(nTau{1},squeeze(mean(mean(STA_SSmat(rois{2},:,:),1),2))',squeeze(std(mean(STA_SSmat(rois{2},:,:),1),0,2))',cmap(2,:));
xlabel('Peri-spike time (ms)')
ax1=[ax1 nexttile([1 1])];
errorbar_shade(nTau{2},squeeze(mean(mean(STA_CSmat(rois{1},:,:),1),2))',squeeze(std(mean(STA_CSmat(rois{1},:,:),1),0,2))',cmap(1,:)); hold all
errorbar_shade(nTau{2},squeeze(mean(mean(STA_CSmat(rois{2},:,:),1),2))',squeeze(std(mean(STA_CSmat(rois{2},:,:),1),0,2))',cmap(2,:));
xlabel('Peri-spike time (ms)')
ax1=[ax1 nexttile([1 1])];
errorbar_shade(nTau{1},squeeze(mean(mean(STA_SSmat(rois{1},:,:),1),2))',squeeze(std(mean(STA_SSmat(rois{1},:,:),1),0,2))',cmap(3,:)); hold all
errorbar_shade(nTau{2},squeeze(mean(mean(STA_CSmat(rois{1},:,:),1),2))',squeeze(std(mean(STA_CSmat(rois{1},:,:),1),0,2))',cmap(6,:));
xlabel('Peri-spike time (ms)')
ax1=[ax1 nexttile([1 1])];
errorbar_shade(nTau{1},squeeze(mean(mean(STA_SSmat(rois{2},:,:),1),2))',squeeze(std(mean(STA_SSmat(rois{2},:,:),1),0,2))',cmap(3,:)); hold all
errorbar_shade(nTau{2},squeeze(mean(mean(STA_CSmat(rois{2},:,:),1),2))',squeeze(std(mean(STA_CSmat(rois{2},:,:),1),0,2))',cmap(6,:));
xlabel('Peri-spike time (ms)')
linkaxes(ax1,'xy')

figure(4); clf;
tiledlayout(ceil(length(CS_s)/20/2),4)
for t=1:floor(length(CS_s)/20)
    nexttile([1 1])
    imagesc(squeeze(mean(STA_CSmat(dist_order,(t-1)*20+1:t*20,:),2)),cax)
    colormap('turbo')
    nexttile([1 1])
    errorbar_shade(nTau{2},squeeze(mean(mean(STA_CSmat(rois{1},(t-1)*20+1:t*20,:),1),2))',squeeze(std(mean(STA_CSmat(rois{1},(t-1)*20+1:t*20,:),1),0,2))',cmap(1,:)); hold all
    errorbar_shade(nTau{2},squeeze(mean(mean(STA_CSmat(rois{2},(t-1)*20+1:t*20,:),1),2))',squeeze(std(mean(STA_CSmat(rois{2},(t-1)*20+1:t*20,:),1),0,2))',cmap(2,:));
    ylim([-1 4])
end

% figure(12); clf;
% for t=10:10:120
% tfit=[-t:-2];
% dvdt_SSmat=[]; dvdt_CSmat=[];
% tss=-nTau{1}(1)+tfit+1;
% tcs=-nTau{2}(1)+tfit+1;
% ft = fittype('poly1'); % 'poly1' specifies a first-degree polynomial (linear fit)
% for s=1:size(STA_SSmat,2)
% for n=1:size(STA_SSmat,1)
%     [fitresult, gof] = fit(tfit', squeeze(STA_SSmat(n,s,tss)), ft);
%     dvdt_SSmat(n,s)=fitresult.p1;
% end
% end
% for s=1:size(STA_CSmat,2)
%     for n=1:size(STA_SSmat,1)
% [fitresult, gof] = fit(tfit', squeeze(STA_CSmat(n,s,tcs)), ft);
% dvdt_CSmat(n,s)=fitresult.p1;
%     end
% end
% nexttile([1 1])
% plot(mean(dvdt_SSmat(rois{1},:),1),mean(dvdt_SSmat(rois{2},:),1),'.','markersize',10); hold all
% plot(mean(dvdt_CSmat(rois{1},:),1),mean(dvdt_CSmat(rois{2},:),1),'.','markersize',10); hold all
% xlabel('Basal dV/dt')
% ylabel('Apical dV/dt')
% legend({'SS','CS'})
% title(['From ' num2str(-tfit(1)) 'ms to ' num2str(-tfit(end)) 'ms before spike'])
% end
%
% figure(13); clf;
% dvdt_t_SSmat=[]; dvdt_t_CSmat=[]; ax1=[];
% myVideo = VideoWriter(['dvdtV_basalvsapical'],"MPEG-4"); %open video file
% myVideo.FrameRate = 5;  %can adjust this, 5 - 10 works well for me
% myVideo.Quality= 100; open(myVideo);
% tseq=[130:-5:0];
% for t=1:length(tseq)
% tfit=[-(tseq(t)+2):-(tseq(t)-8)];
% tss=-nTau{1}(1)+tfit+1;
% tcs=-nTau{2}(1)+tfit+1;
% ft = fittype('poly1'); % 'poly1' specifies a first-degree polynomial (linear fit)
% for s=1:size(STA_SSmat,2)
% for n=1:size(STA_SSmat,1)
%     [fitresult, gof] = fit(tfit', squeeze(STA_SSmat(n,s,tss)), ft);
%     dvdt_t_SSmat(n,s,t)=fitresult.p1;
% end
% end
% for s=1:size(STA_CSmat,2)
%     for n=1:size(STA_SSmat,1)
% [fitresult, gof] = fit(tfit', squeeze(STA_CSmat(n,s,tcs)), ft);
% dvdt_t_CSmat(n,s,t)=fitresult.p1;
%     end
% end
% end
% %ax1=[ax1 nexttile([1 1])];
% for t=1:length(tseq)
% clf;
% tfit=[-(tseq(t)+2):-(tseq(t)-8)];
% tss=-nTau{1}(1)+tfit+1;
% tcs=-nTau{2}(1)+tfit+1;
% %plot(mean(dvdt_t_SSmat(rois{1},:,t),1),mean(dvdt_t_SSmat(rois{2},:,t),1),'.','markersize',10); hold all
% %plot(mean(dvdt_t_CSmat(rois{1},:,t),1),mean(dvdt_t_CSmat(rois{2},:,t),1),'.','markersize',10); hold all
% plot(mean(dvdt_t_SSmat(rois{1},:,t),1),mean(STA_SSmat(rois{2},:,tss),[1 3]),'.','markersize',10); hold all
% plot(mean(dvdt_t_CSmat(rois{1},:,t),1),mean(STA_CSmat(rois{2},:,tcs),[1 3]),'.','markersize',10); hold all
% xlabel('Basal dV/dt (F/ms)'); %ylabel('Apical dV/dt (F/ms)');
% ylabel('Distal V');
% xlim([-1 1]); ylim([-1 7]);
% legend({'SS','CS'})
% title(['From ' num2str(-tfit(1)) 'ms to ' num2str(-tfit(end)) 'ms before spike'])
% frame = getframe(gcf); %get frame
% writeVideo(myVideo, frame);
% end
% close(myVideo);

%% calculate SS/CS property
bwCS=bwlabel(Result.CStrace);
CSstats = regionprops(bwCS, 'all');
CSstats = CSstats(CS_list);
STA_CSmat_sub=STA_CSmat;%-prctile(STA_CSmat,30,3);
dt=5; %ms

subthreshold_time=25; %ms
CSprop=[];
for c=1:size(STA_CSmat_sub,2)
    CSprop.area(:,c)=sum(STA_CSmat_sub(:,c,-nTau{2}(1)+1:-nTau{2}(1)+CSstats(c).Area),3);
    CSprop.density(:,c)=mean(STA_CSmat_sub(:,c,-nTau{2}(1)+1:-nTau{2}(1)+CSstats(c).Area),3);
    CSprop.preArea(:,c)=sum(STA_CSmat_sub(:,c,-nTau{2}(1)-subthreshold_time:-nTau{2}(1)),3);
    CSprop.preAmp(:,c)=max(movmedian(STA_CSmat_sub(:,c,-nTau{2}(1)-subthreshold_time:-nTau{2}(1)),5,3),[],3);
    CSprop.length(c)=CSstats(c).Area;
    ft = fittype('poly1'); % 'poly1' specifies a first-degree polynomial (linear fit)
    [fitresult, gof] = fit(geodist, CSprop.preArea(:,c), ft);
    CSprop.dVdx(c)=fitresult.p1;
    CSprop.intercept_preArea(c)=fitresult.p2;

    CSprop.nSpike(c)=max(bwlabel(max(Result.spike(:,CS_s(c):CS_s(c)+CSstats(c).Area),[],1)));
    sp_t=find(STA_CSpikemat(1,c,-nTau{2}(1)+1:-nTau{2}(1)+CSstats(c).Area))-nTau{2}(1);
    CSprop.SpikeAmp{c}=squeeze(max(permute(reshape(STA_CSmat(:,c,sp_t+[0:3]),size(STA_CSmat,1),length(sp_t),4),[1 3 2]),[],2)); %roi, time, spike -> roi, spike amp
    for t=1:-nTau{2}(1)/dt
        for n=1:size(STA_CSmat,1)
    tfit=[t*dt:(t+1)*dt]-2;            
    [fitresult, gof] = fit(tfit', squeeze(STA_CSmat(n,c,tfit)),ft);
    CSprop.dVdt(n,c,t)=fitresult.p1;
        end
    end
    
    for t=1:4
    for n=1:size(STA_CSmat,1)
    tfit_post=-nTau{2}(1)+[1:3]+t;
    [fitresult, gof] = fit([tfit_post]', squeeze(STA_CSmat(n,c,tfit_post)),ft);
    CSprop.dVdt_post(n,c,t)=fitresult.p1;
    end
    end

    subth=[];
    for n=1:size(STA_CSmat_sub,1)
        subth(n,:)=get_subthreshold(squeeze(STA_CSmat_sub(n,c,:)),max(Result.spike(:,CS_s(c)+nTau{2}),[],1),7,20);
    end
    CSprop.ADPamp(:,c)=max(subth(:,-nTau{2}(1)+1:-nTau{2}(1)+CSstats(c).Area),[],2);
    CS_subspike=bwlabel(max(Result.spike(:,CS_s(c):CS_s(c)+CSstats(c).Area),[],1));
    [Nsp, sp_onset]=unique(CS_subspike);
    sp_onset=sp_onset(find(Nsp>0));
    CSprop.Sp_time{c}=sp_onset;
    CSprop.meanISI(c)=mean(sp_onset(2:end)-sp_onset(1:end-1));
    CSprop.nSpike(c)=length(Nsp)-1;
end

SSprop=[];
for c=1:size(STA_SSmat,2)
    SSprop.area(:,c)=sum(STA_SSmat(:,c,-nTau{1}(1)+1:-nTau{1}(1)+10),3);
    SSprop.preArea(:,c)=sum(STA_SSmat(:,c,-nTau{1}(1)-subthreshold_time:-nTau{1}(1)),3);
    SSprop.preAmp(:,c)=max(movmedian(STA_SSmat(:,c,-nTau{1}(1)-subthreshold_time:-nTau{1}(1)),5,3),[],3);
    ft = fittype('poly1'); % 'poly1' specifies a first-degree polynomial (linear fit)
    [fitresult, gof] = fit(geodist, SSprop.preArea(:,c), ft);
    SSprop.dVdx(c)=fitresult.p1;
    SSprop.intercept_preArea(c)=fitresult.p2;
    SSprop.SpikeAmp{c}=squeeze(max(permute(reshape(STA_SSmat(:,c,-nTau{1}(1)+[0:3]),size(STA_CSmat,1),1,4),[1 3 2]),[],2)); %roi, time, spike -> roi, spike amp
     
    for t=1:-nTau{1}(1)/dt
        for n=1:size(STA_SSmat,1)
    tfit=[t*dt:(t+1)*dt]-2;
    [fitresult, gof] = fit(tfit', squeeze(STA_SSmat(n,c,tfit)),ft);
    SSprop.dVdt(n,c,t)=fitresult.p1;
        end
    end
    
    for t=1:4
    for n=1:size(STA_SSmat,1)
    tfit_post=-nTau{1}(1)+[1:3]+t;
    [fitresult, gof] = fit([tfit_post]', squeeze(STA_SSmat(n,c,tfit_post)),ft);
    SSprop.dVdt_post(n,c,t)=fitresult.p1;
    end
    end

    subth=[];
    for n=1:size(STA_SSmat,1)
        subth(n,:)=get_subthreshold(squeeze(STA_SSmat(n,c,:)),max(Result.spike(:,bAP_s(c)+nTau{1}),[],1),5,20);
    end
    SSprop.ADPamp(:,c)=max(subth(:,-nTau{1}(1)+1:-nTau{1}(1)+10),[],2);
end
%%
figure(13); clf;
nexttile([1 1])
plot(sum(SSprop.preArea,1),SSprop.dVdx,'.','markersize',12); hold all
plot(sum(CSprop.preArea,1),CSprop.dVdx,'.','markersize',12); hold all
% h=histogram(SSprop.poly_preArea,30,'Normalization','probability'); hold all
% histogram(CSprop.poly_preArea,h.BinEdges,'Normalization','probability')
legend({'SS','CS'});
ylabel('dV/dx')
xlabel('Sum of Pre-area')
title(num2str(subthreshold_time))

nexttile([1 1])
h=histogram(squeeze(mean(SSprop.preArea(rois{1},:),1)),30,'Normalization','probability'); hold all
histogram(squeeze(mean(CSprop.preArea(rois{1},:),1)),h.BinEdges,'Normalization','probability')
legend({'SS','CS'});
xlabel('Pre-area (Basal)')

nexttile([1 1])
h=histogram(squeeze(mean(SSprop.preArea(rois{2},:),1)),30,'Normalization','probability'); hold all
histogram(squeeze(mean(CSprop.preArea(rois{2},:),1)),h.BinEdges,'Normalization','probability')
legend({'SS','CS'});
xlabel('Pre-area (Distal)')

nexttile([1 1])
h=histogram(squeeze(mean(SSprop.dVdt_post(rois{1},:),1)),30,'Normalization','probability'); hold all
histogram(squeeze(mean(CSprop.dVdt_post(rois{1},:),1)),h.BinEdges,'Normalization','probability')
legend({'SS','CS'});
xlabel('Post-spike dV/dt (Basal)')

nexttile([1 1])
h=histogram(squeeze(mean(SSprop.dVdt_post(rois{2},:),1)),30,'Normalization','probability'); hold all
histogram(squeeze(mean(CSprop.dVdt_post(rois{2},:),1)),h.BinEdges,'Normalization','probability')
legend({'SS','CS'});
xlabel('Post-spike dV/dt (distal)')

nexttile([1 1])
h=histogram(squeeze(mean(SSprop.dVdt(rois{1},:,end-1),1)),30,'Normalization','probability'); hold all
histogram(squeeze(mean(CSprop.dVdt(rois{1},:,end-1),1)),h.BinEdges,'Normalization','probability')
legend({'SS','CS'});
xlabel('dV/dt (Basal)')

nexttile([1 1])
h=histogram(squeeze(mean(SSprop.dVdt(rois{2},:,end-1),1)),30,'Normalization','probability'); hold all
histogram(squeeze(mean(CSprop.dVdt(rois{2},:,end-1),1)),h.BinEdges,'Normalization','probability')
legend({'SS','CS'});
xlabel('dV/dt (Distal)')

nexttile([1 1])
h=histogram(SSprop.dVdx,30,'Normalization','probability'); hold all
histogram(CSprop.dVdx,h.BinEdges,'Normalization','probability')
legend({'SS','CS'});
xlabel('dV/dx')

nexttile([1 1])
h=histogram(mean(SSprop.preAmp(rois{1},:),1),30,'Normalization','probability'); hold all
histogram(mean(CSprop.preAmp(rois{1},:),1),h.BinEdges,'Normalization','probability')
legend({'SS','CS'});
xlabel('Sum of Pre-Amp. (Basal)')

nexttile([1 1])
h=histogram(mean(SSprop.preAmp(rois{2},:),1),30,'Normalization','probability'); hold all
histogram(mean(CSprop.preAmp(rois{2},:),1),h.BinEdges,'Normalization','probability')
legend({'SS','CS'});
xlabel('Sum of Pre-Amp. (Distal)')

nexttile([1 1])
h=histogram(mean(CSprop.ADPamp(rois{1},:),1),30,'Normalization','probability'); hold all
histogram(mean(CSprop.ADPamp(rois{2},:),1),h.BinEdges,'Normalization','probability')
legend({'Basal','Distal'});
xlabel('Complex spike amplitude')

ISImat=[];
for c=1:size(CSprop.Sp_time,2)
    ISImat(c,1:length(CSprop.Sp_time{c})-1)=CSprop.Sp_time{c}(2:end)-CSprop.Sp_time{c}(1:end-1);
end
ISImat(ISImat==0)=NaN;
nexttile([1 1])
plot([1:5],ISImat(:,1:5),'color',[0.6 0.6 0.6]); hold all
errorbar([1:5],mean(ISImat(:,1:5),1,'omitnan'),std(ISImat(:,1:5),0,1,'omitnan'),'color',[1 0 0]); hold all
xlabel('N^t^h Spike'); ylabel('ISI (ms)');
xlim([0.5 5.5]); ylim([0 30])

nexttile([1 1])
plot(ISImat(:,1),sum(CSprop.ADPamp(rois{2},:),1),'.','markersize',12); hold all
plot(ISImat(:,2),sum(CSprop.ADPamp(rois{2},:),1),'.','markersize',12)
xlabel('ISI (ms)'); ylabel('Complex spike amplitude (Distal)');
legend({'2nd - 1st','3rd - 2nd'});
xlim([3 30])

nexttile([1 1]);
plot(mean(SSprop.dVdt_post(rois{1},:,2),1),squeeze(mean(STA_SSmat(rois{1},:,-nTau{1}(1)+[3:5]),[1 3])),'.','markersize',12); hold all
plot(mean(CSprop.dVdt_post(rois{1},:,2),1),squeeze(mean(STA_CSmat(rois{1},:,-nTau{2}(1)+[3:5]),[1 3])),'.','markersize',12);
xlabel('Post-spike dV/dt (Basal)'); ylabel('Post-spike V (Basal)');
title(['2-4 ms after spike'])
nexttile([1 1]);
plot(mean(SSprop.dVdt_post(rois{2},:,2),1),squeeze(mean(STA_SSmat(rois{2},:,-nTau{1}(1)+[3:5]),[1 3])),'.','markersize',12); hold all
plot(mean(CSprop.dVdt_post(rois{2},:,2),1),squeeze(mean(STA_CSmat(rois{2},:,-nTau{2}(1)+[3:5]),[1 3])),'.','markersize',12);
xlabel('Post-spike dV/dt (Distal)'); ylabel('Post-spike V (Distal)');
title(['2-4 ms after spike'])

nexttile([1 1]);
plot(mean(SSprop.dVdt_post(rois{1},:,3),1),squeeze(mean(STA_SSmat(rois{1},:,-nTau{1}(1)+[4:6]),[1 3])),'.','markersize',12); hold all
plot(mean(CSprop.dVdt_post(rois{1},:,3),1),squeeze(mean(STA_CSmat(rois{1},:,-nTau{2}(1)+[4:6]),[1 3])),'.','markersize',12);
xlabel('Post-spike dV/dt (Basal)'); ylabel('Post-spike V (Basal)');
title(['3-5 ms after spike'])
nexttile([1 1]);
plot(mean(SSprop.dVdt_post(rois{2},:,3),1),squeeze(mean(STA_SSmat(rois{2},:,-nTau{1}(1)+[4:6]),[1 3])),'.','markersize',12); hold all
plot(mean(CSprop.dVdt_post(rois{2},:,3),1),squeeze(mean(STA_CSmat(rois{2},:,-nTau{2}(1)+[4:6]),[1 3])),'.','markersize',12);
xlabel('Post-spike dV/dt (Distal)'); ylabel('Post-spike V (Distal)');
title(['3-5 ms after spike'])

nexttile([1 1]);
plot(mean(SSprop.dVdt_post(rois{1},:,4),1),squeeze(mean(STA_SSmat(rois{1},:,-nTau{1}(1)+[5:7]),[1 3])),'.','markersize',12); hold all
plot(mean(CSprop.dVdt_post(rois{1},:,4),1),squeeze(mean(STA_CSmat(rois{1},:,-nTau{2}(1)+[5:7]),[1 3])),'.','markersize',12);
xlabel('Post-spike dV/dt (Basal)'); ylabel('Post-spike V (Basal)');
title(['4-6 ms after spike'])
nexttile([1 1]);
plot(mean(SSprop.dVdt_post(rois{2},:,4),1),squeeze(mean(STA_SSmat(rois{2},:,-nTau{1}(1)+[5:7]),[1 3])),'.','markersize',12); hold all
plot(mean(CSprop.dVdt_post(rois{2},:,4),1),squeeze(mean(STA_CSmat(rois{2},:,-nTau{2}(1)+[5:7]),[1 3])),'.','markersize',12);
xlabel('Post-spike dV/dt (Distal)'); ylabel('Post-spike V (Distal)');
title(['4-6 ms after spike'])
%%
figure(14); clf;
for t=1:10
    nexttile([1 1])
plot(mean(SSprop.preArea(rois{1},:),1),squeeze(mean(SSprop.dVdt(rois{1},:,end-2-t),1)),'.'); hold all
plot(mean(CSprop.preArea(rois{1},:),1),squeeze(mean(CSprop.dVdt(rois{1},:,end-2-t),1)),'.')
title('Basal')
    nexttile([1 1])
plot(mean(SSprop.preArea(rois{2},:),1),squeeze(mean(SSprop.dVdt(rois{2},:,end-2-t),1)),'.'); hold all
plot(mean(CSprop.preArea(rois{2},:),1),squeeze(mean(CSprop.dVdt(rois{2},:,end-2-t),1)),'.')
title('Distal')
end


%% Calculate Phase
FFT_SSmat=[]; FFT_CSmat=[]; Phase_SSmat=[]; Phase_CSmat=[]; PhaseMag_SSmat=[]; PhaseMag_CSmat=[];
Hilbert_SSmat=[]; Hilbert_CSmat=[];
theta_pass_frq=[7 35]; 
freq_lowhigh=theta_pass_frq/(1000/2);
[b4, a4] = butter(4, freq_lowhigh, 'bandpass');

for s=1:size(STA_SSmat,2)
    for n=1:size(STA_SSmat,1)
        [FFT_SSmat(n,s,:) freqSS] = fft_simple(squeeze(STA_SSmat(n,s,1:-nTau{1}(1))), 1000);
        x=[1:-nTau{1}(1)]/1000;
        trTofit=squeeze(STA_SSmat(n,s,1:-nTau{1}(1)))';
        trTofit_filt=filtfilt(b4, a4, trTofit);
        Hilbert_SSmat(n,s,:) = hilbert(trTofit_filt);
        Phase_SSmat(n,s,:) = angle(hilbert(trTofit_filt));
        PhaseMag_SSmat(n,s,:)=abs(hilbert(trTofit_filt));
    end
end

for s=1:size(STA_CSmat,2)
    for n=1:size(STA_CSmat,1)
        [FFT_CSmat(n,s,:) freqSS] = fft_simple(squeeze(STA_CSmat(n,s,1:-nTau{1}(1))), 1000);
        x=[1:-nTau{1}(1)]/1000;
        trTofit=squeeze(STA_CSmat(n,s,1:-nTau{1}(1)))';
        trTofit_filt=filtfilt(b4, a4, trTofit);
        Hilbert_CSmat(n,s,:) = hilbert(trTofit_filt);
        Phase_CSmat(n,s,:) = angle(hilbert(trTofit_filt));
        PhaseMag_CSmat(n,s,:)=abs(hilbert(trTofit_filt));
    end
end

figure(10); clf;
ax1=[]; ax2=[];
for s=1:10
    nexttile([1 1])
    for n=5
        trTofit=squeeze(STA_CSmat(n,s,1:-nTau{1}(1)))';
        trTofit_filt=filtfilt(b4, a4, trTofit);
        plot(trTofit); hold all
        plot(trTofit_filt);
        plot(angle(squeeze(Hilbert_CSmat(n,s,:))))
    end
end

figure(11); clf;
title_str={'Distal D1','Distal D2','Basal D1','Basal D2','Basal D3'};
for t=-5:1:-1
     g=1;
for n=[5 6 12 14 15]
    nexttile([1 1])
polarhistogram(squeeze(Phase_SSmat(n,:,end+t)),30,'Normalization','probability'); hold all
polarhistogram(squeeze(Phase_CSmat(n,:,end+t)),30,'Normalization','probability');
title([title_str{g} ', ' num2str(-t) ' ms before spike']); g=g+1;
end
end
legend({'SS','CS'});

figure(13); clf;
t=[-5:1:0];
nexttile([1 1]);
plot3(t,squeeze(Phase_SSmat([14],:,end+t)),squeeze(Phase_SSmat([5],:,end+t)),'color','r'); hold all
plot3(t,squeeze(Phase_CSmat([14],:,end+t)),squeeze(Phase_CSmat([5],:,end+t)),'color','b')
polar(squeeze(Phase_SSmat([5],:,end+t))',squeeze(PhaseMag_SSmat([5],:,end+t))','r'); hold all
polar(squeeze(Phase_CSmat([5],:,end+t))',squeeze(PhaseMag_CSmat([5],:,end+t))','b'); hold all

figure(12); clf; ax1=[];
for c=1:40
ax1=[ax1 nexttile([1 1])];
%imagesc(squeeze(Phase_SSmat(dist_order,c,:))); hold all
imagesc([squeeze(Phase_SSmat(dist_order,c,:)) squeeze(STA_SSmat(dist_order,c,-nTau{1}(1)+1:end))],2*[-pi pi]); hold all
%plot(squeeze(STA_SSmat(1,c,:)))
axis tight
title('SS')
end
for c=1:40
ax1=[ax1 nexttile([1 1])];
%imagesc(squeeze(Phase_CSmat(dist_order,c,:))); hold all
imagesc([squeeze(Phase_CSmat(dist_order,c,:)) squeeze(STA_CSmat(dist_order,c,-nTau{2}(1)+1:end))],2*[-pi pi]); hold all
%plot(squeeze(STA_CSmat(1,c,:)))
axis tight
title('CS')
end
colormap('turbo')
linkaxes(ax1','x')
%%
figure; clf; cmap=distinguishable_colors(2);
ax1=[];
for t=1:22
    ax1=[ax1 nexttile([1 1])];
    plot(squeeze(mean(STA_SSmat(rois{1},:,-nTau{1}(1)-20+t),1)),squeeze(mean(STA_SSmat(rois{2},:,-nTau{1}(1)-20+t),1)),'.','color',cmap(1,:));
    hold all
    plot(squeeze(mean(STA_CSmat_sub(rois{1},:,-nTau{2}(1)-20+t),1)),squeeze(mean(STA_CSmat_sub(rois{2},:,-nTau{2}(1)-20+t),1)),'.','color',cmap(2,:));
    xlabel('Basal fluorescence')
    ylabel('Distal fluorescence')
    title([num2str(-(-20+t)) 'ms before spike'])
end
linkaxes(ax1,'xy')
%%
figure(1); clf; cmap=distinguishable_colors(2); ax1=[];
for t=5:10:95
    preAreaCS=mean(STA_CSmat_sub(:,:,-nTau{2}(1)-t:-nTau{2}(1)-1),3);
    preAreaSS=mean(STA_SSmat(:,:,-nTau{1}(1)-t:-nTau{1}(1)-1),3);
    preAreaCSbs=mean(preAreaCS(rois{1},:)); preAreaCSap=mean(preAreaCS(rois{2},:));
    preAreaSSbs=mean(preAreaSS(rois{1},:)); preAreaSSap=mean(preAreaSS(rois{2},:));
    [pap]=ranksum(preAreaSSap,preAreaCSap);
    [pbs]=ranksum(preAreaSSbs,preAreaCSbs);
    ax1=[ax1 nexttile([1 1])];
    plot(preAreaSSbs,preAreaSSap,'.','color',cmap(1,:),'markersize',10); hold all
    plot(preAreaCSbs,preAreaCSap,'.','color',cmap(2,:),'markersize',10)
    xlabel('Basal fluorescence')
    ylabel('Apical fluorescence')
    legend({'SS','CS'})
    title(['Average from ' num2str(t+1) ' ms to ' num2str(2) ' ms before spike' ...
        '\newline p-basal : ' num2str(pbs,2) ', p-apical : ' num2str(pap,2)])
end
linkaxes(ax1,'xy')
%%
figure(5); clf;
PreBasalApical=[mean(CSprop.preArea(rois{1},:))' mean(CSprop.preArea(rois{2},:))'];

nexttile([1 1])
plot(mean(CSprop.preArea(rois{1},:)),mean(CSprop.area(rois{1},:)),'.','markersize',10); hold all
plot(mean(CSprop.preArea(rois{2},:)),mean(CSprop.area(rois{1},:)),'.','markersize',10)
[b bp]=corr(mean(CSprop.preArea(rois{1},:))',mean(CSprop.area(rois{1},:))');
[a ap]=corr(mean(CSprop.preArea(rois{2},:))',mean(CSprop.area(rois{1},:))');
xlabel('Area of Pre-spiking subthreshold')
ylabel('Area of Complex spike (basal)')
legend({['Basal, corr:' num2str(b,2), ', p:' num2str(bp,2)], ...
    ['Apical, corr:' num2str(a,2), ', p:' num2str(ap,2)]})

nexttile([1 1])
plot(mean(CSprop.preArea(rois{1},:)),mean(CSprop.area(rois{2},:)),'.','markersize',10); hold all
plot(mean(CSprop.preArea(rois{2},:)),mean(CSprop.area(rois{2},:)),'.','markersize',10)
[b bp]=corr(mean(CSprop.preArea(rois{1},:))',mean(CSprop.area(rois{2},:))');
[a ap]=corr(mean(CSprop.preArea(rois{2},:))',mean(CSprop.area(rois{2},:))');
xlabel('Area of Pre-spiking subthreshold')
ylabel('Area of Complex spike (apical)')
legend({['Basal, corr:' num2str(b,2), ', p:' num2str(bp,2)], ...
    ['Apical, corr:' num2str(a,2), ', p:' num2str(ap,2)]})

nexttile([1 1])
plot(mean(CSprop.preArea(rois{1},:)),CSprop.area(1,:),'.','markersize',10); hold all
plot(mean(CSprop.preArea(rois{2},:)),CSprop.area(1,:),'.','markersize',10)
[b bp]=corr(mean(CSprop.preArea(rois{1},:))',CSprop.area(1,:)');
[a ap]=corr(mean(CSprop.preArea(rois{2},:))',CSprop.area(1,:)');
xlabel('Area of Pre-spiking subthreshold')
ylabel('Area of Complex spike (soma)')
legend({['Basal, corr:' num2str(b,2), ', p:' num2str(bp,2)], ...
    ['Apical, corr:' num2str(a,2), ', p:' num2str(ap,2)]})

nexttile([1 1])
plot(mean(CSprop.preAmp(rois{1},:)),mean(CSprop.area(rois{1},:)),'.','markersize',10); hold all
plot(mean(CSprop.preAmp(rois{2},:)),mean(CSprop.area(rois{1},:)),'.','markersize',10)
[b bp]=corr(mean(CSprop.preAmp(rois{1},:))',mean(CSprop.area(rois{1},:))');
[a ap]=corr(mean(CSprop.preAmp(rois{2},:))',mean(CSprop.area(rois{1},:))');
xlabel('Amplitude of Pre-spiking subthreshold')
ylabel('Area of Complex spike (basal)')
legend({['Basal, corr:' num2str(b,2), ', p:' num2str(bp,2)], ...
    ['Apical, corr:' num2str(a,2), ', p:' num2str(ap,2)]})

nexttile([1 1])
plot(mean(CSprop.preAmp(rois{1},:)),mean(CSprop.area(rois{2},:)),'.','markersize',10); hold all
plot(mean(CSprop.preAmp(rois{2},:)),mean(CSprop.area(rois{2},:)),'.','markersize',10)
[b bp]=corr(mean(CSprop.preAmp(rois{1},:))',mean(CSprop.area(rois{2},:))');
[a ap]=corr(mean(CSprop.preAmp(rois{2},:))',mean(CSprop.area(rois{2},:))');
xlabel('Amplitude of Pre-spiking subthreshold')
ylabel('Area of Complex spike (apical)')
legend({['Basal, corr:' num2str(b,2), ', p:' num2str(bp,2)], ...
    ['Apical, corr:' num2str(a,2), ', p:' num2str(ap,2)]})

nexttile([1 1])
plot(mean(CSprop.preAmp(rois{1},:)),CSprop.area(1,:),'.','markersize',10); hold all
plot(mean(CSprop.preAmp(rois{2},:)),CSprop.area(1,:),'.','markersize',10)
[b bp]=corr(mean(CSprop.preAmp(rois{1},:))',CSprop.area(1,:)');
[a ap]=corr(mean(CSprop.preAmp(rois{2},:))',CSprop.area(1,:)');
xlabel('Amplitude of Pre-spiking subthreshold')
ylabel('Area of Complex spike (soma)')
legend({['Basal, corr:' num2str(b,2), ', p:' num2str(bp,2)], ...
    ['Apical, corr:' num2str(a,2), ', p:' num2str(ap,2)]})


%%
figure(6); clf;
tiledlayout(3,4)
dvdx_range=[prctile(CSprop.dVdx,30) prctile(CSprop.dVdx,70)]; ax1=[];
c_show=find(CSprop.dVdx<dvdx_range(1)); c_show2=find(CSprop.dVdx>dvdx_range(2));
s_show=find(SSprop.dVdx<dvdx_range(1)); s_show2=find(SSprop.dVdx>dvdx_range(2));
d_label={'B','B','B','B','B','S','A','A','A','A','A','A','D','D','D'};
ax1=[ax1 nexttile(1,[1 1])];
imagesc(squeeze(mean(STA_CSmat(dist_order,c_show2,:),2)),[-0.5 3])
colormap(turbo)
title(['High dV/dx'])
ax1=[ax1 nexttile(2,[1 1])];
imagesc(squeeze(mean(STA_CSmat(dist_order,c_show,:),2)),[-0.5 3])
colormap(turbo)
title(['Low dV/dx'])
ax1=[ax1 nexttile(5,[1 1])];
l=plot(squeeze(mean(STA_CSmat(dist_order,c_show2,:),2))');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2))
ylim([-2 10])
ax1=[ax1 nexttile(6,[1 1])];
l=plot(squeeze(mean(STA_CSmat(dist_order,c_show,:),2))');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2))
ylim([-2 10])

ax1=[ax1 nexttile(9,[1 1])];
M=squeeze(mean(STA_CSmat(rois{1},c_show2,:),[1 2]));
S=squeeze(std(mean(STA_CSmat(rois{1},c_show2,:),1),0,2));
errorbar_shade([1:size(STA_CSmat,3)],M,S,cmap(1,:)); hold all
M=squeeze(mean(STA_CSmat(rois{1},c_show,:),[1 2]));
S=squeeze(std(mean(STA_CSmat(rois{1},c_show,:),1),0,2));
errorbar_shade([1:size(STA_CSmat,3)],M,S,cmap(2,:)); hold all
title(['Basal'])

ax1=[ax1 nexttile(10,[1 1])];
M=squeeze(mean(STA_CSmat(rois{2},c_show2,:),[1 2]));
S=squeeze(std(mean(STA_CSmat(rois{2},c_show2,:),1),0,2));
errorbar_shade([1:size(STA_CSmat,3)],M,S,cmap(1,:)); hold all
M=squeeze(mean(STA_CSmat(rois{2},c_show,:),[1 2]));
S=squeeze(std(mean(STA_CSmat(rois{2},c_show,:),1),0,2));
errorbar_shade([1:size(STA_CSmat,3)],M,S,cmap(2,:)); hold all
title(['Distal'])

ax1=[ax1 nexttile(3,[1 1])];
imagesc(squeeze(mean(STA_SSmat(dist_order,s_show2,:),2)),[-0.5 3])
colormap(turbo); title(['High dV/dx'])
ax1=[ax1 nexttile(4,[1 1])];
imagesc(squeeze(mean(STA_SSmat(dist_order,s_show,:),2)),[-0.5 3])
colormap(turbo); title(['Low dV/dx'])

ax1=[ax1 nexttile(7,[1 1])];
l=plot(squeeze(mean(STA_SSmat(dist_order,s_show2,:),2))');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2))
ylim([-2 10]); 
ax1=[ax1 nexttile(8,[1 1])];
l=plot(squeeze(mean(STA_SSmat(dist_order,s_show,:),2))');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2))
ylim([-2 10])
legend(d_label)

ax1=[ax1 nexttile(11,[1 1])];
M=squeeze(mean(STA_SSmat(rois{1},s_show2,:),[1 2]));
S=squeeze(std(mean(STA_SSmat(rois{1},s_show2,:),1),0,2));
errorbar_shade([1:size(STA_SSmat,3)],M,S,cmap(1,:)); hold all
M=squeeze(mean(STA_SSmat(rois{1},s_show,:),[1 2]));
S=squeeze(std(mean(STA_SSmat(rois{1},s_show,:),1),0,2));
errorbar_shade([1:size(STA_SSmat,3)],M,S,cmap(2,:)); hold all
title(['Basal'])

ax1=[ax1 nexttile(12,[1 1])];
M=squeeze(mean(STA_SSmat(rois{2},s_show2,:),[1 2]));
S=squeeze(std(mean(STA_SSmat(rois{2},s_show2,:),1),0,2));
errorbar_shade([1:size(STA_SSmat,3)],M,S,cmap(1,:)); hold all
M=squeeze(mean(STA_SSmat(rois{2},s_show,:),[1 2]));
S=squeeze(std(mean(STA_SSmat(rois{2},s_show,:),1),0,2));
errorbar_shade([1:size(STA_SSmat,3)],M,S,cmap(2,:)); hold all
title(['Distal'])
linkaxes(ax1,'x')
%%
dvdx_range=[prctile(CSprop.dVdx,30) prctile(CSprop.dVdx,70)];
c_show=find(CSprop.dVdx<dvdx_range(1)); c_show2=find(CSprop.dVdx>dvdx_range(2));
s_show=find(SSprop.dVdx<dvdx_range(1)); s_show2=find(SSprop.dVdx>dvdx_range(2));
figure(7); clf;
ax1=nexttile([1 1]);
imagesc(squeeze(mean(STA_CSmat(dist_order,c_show2,:),2)),[-0.5 4]);
colormap('turbo')

ax1=nexttile([1 1]);
imagesc(squeeze(mean(STA_CSmat(dist_order,c_show,:),2)),[-0.5 4]);
colormap('turbo')

ax2=nexttile([1 1]);
errorbar_shade(nTau{2},squeeze(mean(STA_CSmat(rois{1},c_show,:),[1 2])),squeeze(std(mean(STA_CSmat(rois{1},c_show,:),1),0,2)),cmap(1,:))
errorbar_shade(nTau{2},squeeze(mean(STA_CSmat(rois{1},c_show2,:),[1 2])),squeeze(std(mean(STA_CSmat(rois{1},c_show2,:),1),0,2)),cmap(2,:))

ax2=nexttile([1 1]);
errorbar_shade(nTau{2},squeeze(mean(STA_CSmat(rois{2},c_show,:),[1 2])),squeeze(std(mean(STA_CSmat(rois{2},c_show,:),1),0,2)),cmap(1,:))
errorbar_shade(nTau{2},squeeze(mean(STA_CSmat(rois{2},c_show2,:),[1 2])),squeeze(std(mean(STA_CSmat(rois{2},c_show2,:),1),0,2)),cmap(2,:))

ax1=nexttile([1 1]);
imagesc(squeeze(mean(STA_SSmat(dist_order,s_show2,:),2)),[-0.5 4]);
colormap('turbo')

ax1=nexttile([1 1]);
imagesc(squeeze(mean(STA_SSmat(dist_order,s_show,:),2)),[-0.5 4]);
colormap('turbo')

ax2=nexttile([1 1]);
errorbar_shade(nTau{1},squeeze(mean(STA_SSmat(rois{1},s_show,:),[1 2])),squeeze(std(mean(STA_SSmat(rois{1},s_show,:),1),0,2)),cmap(1,:))
errorbar_shade(nTau{1},squeeze(mean(STA_SSmat(rois{1},s_show2,:),[1 2])),squeeze(std(mean(STA_SSmat(rois{1},s_show2,:),1),0,2)),cmap(2,:))

ax2=nexttile([1 1]);
errorbar_shade(nTau{1},squeeze(mean(STA_SSmat(rois{2},s_show,:),[1 2])),squeeze(std(mean(STA_SSmat(rois{2},s_show,:),1),0,2)),cmap(1,:))
errorbar_shade(nTau{1},squeeze(mean(STA_SSmat(rois{2},s_show2,:),[1 2])),squeeze(std(mean(STA_SSmat(rois{2},s_show2,:),1),0,2)),cmap(2,:))

%%
figure(9); clf; cmap=distinguishable_colors(2);
bshigh=find(CSprop.dVdx<prctile(CSprop.dVdx,30));
aphigh=find(CSprop.dVdx>prctile(CSprop.dVdx,70));
concatAmpbs=[]; concatAmpap=[];
for c=1:length(bshigh)
    concatAmpbs(1:nROI,1:size(CSprop.SpikeAmp{bshigh(c)},2),c)=CSprop.SpikeAmp{bshigh(c)};
end
for c=1:length(aphigh)
    concatAmpap(1:nROI,1:size(CSprop.SpikeAmp{aphigh(c)},2),c)=CSprop.SpikeAmp{aphigh(c)};
end
concatAmpap(concatAmpap==0)=NaN; concatAmpbs(concatAmpbs==0)=NaN;

Mba=[mean(mean(concatAmpbs(rois{1},:,:),1),3,'omitnan'); mean(mean(concatAmpbs(rois{2},:,:),1),3,'omitnan')];
Map=[mean(mean(concatAmpap(rois{1},:,:),1),3,'omitnan'); mean(mean(concatAmpap(rois{2},:,:),1),3,'omitnan')];
Sba=[std(mean(concatAmpbs(rois{1},:,:),1),0,3,'omitnan'); std(mean(concatAmpbs(rois{2},:,:),1),0,3,'omitnan')];
Sap=[std(mean(concatAmpap(rois{1},:,:),1),0,3,'omitnan'); std(mean(concatAmpap(rois{2},:,:),1),0,3,'omitnan')];
tiledlayout(1,2)
nexttile([1 1])
errorbar_shade([1:size(concatAmpbs,2)],Mba(1,:),Sba(1,:),cmap(1,:)); hold all
errorbar_shade([1:size(concatAmpap,2)],Map(1,:),Sap(1,:),cmap(2,:))
xlabel('N^t^h spike in CS')
ylabel('bAP amplitude (at Basal)')
nexttile([1 1])
errorbar_shade([1:size(concatAmpbs,2)],Mba(2,:),Sba(2,:),cmap(1,:)); hold all
errorbar_shade([1:size(concatAmpap,2)],Map(2,:),Sap(2,:),cmap(2,:))
xlabel('N^t^h spike in CS')
ylabel('bAP amplitude (at Apical)')

ApAmpap=squeeze(mean(concatAmpap(rois{2},:,:),1,'omitnan'));
ApAmpbs=squeeze(mean(concatAmpbs(rois{2},:,:),1,'omitnan'));
for s=1:min(size(concatAmpbs,2),size(concatAmpap,2))
    [~, pval(s)]=ttest2(ApAmpap(s,:),ApAmpbs(s,:));
end
title(num2str(pval,2))


figure(10); clf;
pACS=squeeze(mean(mean(STA_CSmat(rois{2},:,100:298),1),3));
pASS=squeeze(mean(mean(STA_SSmat(rois{2},:,100:298),1),3));
plot3(pASS,sum(SSprop.area,1),SSprop.dVdx,'.'); hold all
plot3(pACS,sum(CSprop.area,1),CSprop.dVdx,'.');
grid on
xlabel('Pre-Area (apical)')
ylabel('Area')
zlabel('dVdx')

%%
figure(15); clf;
tiledlayout(1,4)
[~, dvdxarg]= sort(CSprop.dVdx);
MaxSpAmp=cell2mat(cellfun(@(x) max(x,[],2),CSprop.SpikeAmp,'UniformOutput',false));
nexttile([1 1])
imagesc(CSprop.preArea(dist_order,dvdxarg)')
title(['Pre spike Area'])
nexttile([1 1])
imagesc(CSprop.area(dist_order,dvdxarg)')
title(['Complex spike Area'])
nexttile([1 1])
imagesc(CSprop.ADPamp(dist_order,dvdxarg)')
title(['Pre spike Amplitude'])
nexttile([1 1])
imagesc(CSprop.density(dist_order,dvdxarg)'./MaxSpAmp(dist_order,dvdxarg)')
title(['Complex spike density'])
colormap(turbo)
%% Show CSs
figure(9); clf; cax=[-2 5]; ax1=[];
tiledlayout(10,5)
%c_show=find(CSprop.dVdx<prctile(CSprop.dVdx,30));
%c_show=[1:50];
[~, c_show]=sort(ISImat(:,1));
for c=c_show(1:2:100)'
    ax1=[ax1 nexttile([1 1])];
    imagesc(squeeze(STA_CSmat(dist_order,c,:)),cax)
    colormap("turbo")
    title([num2str(c) ' , ' num2str(CSprop.dVdx(c),2),' , ' num2str(sum(CSprop.ADPamp(rois{2},c)),2)])
end
linkaxes(ax1,'xy')

%% Show SSs
figure(9); clf; cax=[-2 5]; ax1=[];
tiledlayout(10,5)
%c_show=find(CSprop.dVdx<prctile(CSprop.dVdx,30));
c_show=[1:50];
for c=c_show
    ax1=[ax1 nexttile([1 1])];
    imagesc(squeeze(STA_SSmat(dist_order,c,:)),cax)
    colormap("turbo")
    title([num2str(c) ' , ' num2str(SSprop.dVdx(c),2),' , ' num2str(sum(SSprop.ADPamp(rois{2},c)),2)])
end
linkaxes(ax1,'xy')
%%
dt=30; %ms
bin=15;%ms
bin_edge=[bin:bin:10000];
dTrace=Result.normTraces./F_ref;

coord_1d=dim_reduce(get_coord(Result.ftprnt));
geodist=interDendDist(1,:)'.*sign(coord_1d-coord_1d(1));

for t=1:length(bin_edge)-1
toi=[bin_edge(t):bin_edge(t+1)];
dArea=sum(dTrace(:,toi),2);

[fitresult, gof] = fit(geodist, dArea, ft);
dVdxTrace(t)=fitresult.p1;

end