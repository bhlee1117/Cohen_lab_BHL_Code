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

%%
f=18; load(fullfile(fpath{f},'PC_Result.mat'),'Result')
rois={basal_ROI{f},apical_ROI{f}};
nROI=size(Result.normTraces,1);
nTau_bAP=[-20:20];
nTau={[-55:15],[-150:200],[-55:20]}; %SS, CS, dSP
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

bAP_ref=[];
for s=som_spike
    isnearby=sum(ismember(s+nTau_bAP,som_spike))>1;
    isnearbyCS=sum(ismember(s+nTau_bAP,find(Result.CStrace)))>1;
    ispartCS=tr_trace(s)>0;
    if ~isnearby & ~isnan(s) & ~isnearbyCS & ~ispartCS
        bAP_ref=[bAP_ref s];
    end
end

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

dSP_s=[];
for s=find(Result.SpClass(3,:))
    isnearby=sum(ismember(s+nTau{3},som_spike))>1;
    isnearbyCS=sum(ismember(s+nTau{3},find(Result.CStrace)))>1;
    ispartCS=tr_trace(s)>0;
    if ~isnearby & ~isnan(s) & ~isnearbyCS & ~ispartCS
        dSP_s=[dSP_s s];
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

STA_SS=squeeze(mean(reshape(Result.normTraces(:,bAP_s'+nTau_bAP),nROI,[],length(nTau_bAP)),2));
STA_SS= STA_SS;% - mean(STA_SS(:,1:10),2);
F_ref=mean(STA_SS(:,-nTau_bAP(1)+[7:14]),2);

NomalizedTrace=(Result.normTraces)./F_ref;

STA_CSmat=reshape(NomalizedTrace(:,CS_s'+nTau{2}),nROI,[],length(nTau{2}));
STA_SSmat=reshape(NomalizedTrace(:,bAP_s'+nTau{1}),nROI,[],length(nTau{1}));
STA_dSPmat=reshape(NomalizedTrace(:,dSP_s'+nTau{3}),nROI,[],length(nTau{3}));
%%
figure(5); clf; cmap=distinguishable_colors(6); ax1=[];
rois={basal_ROI{f},apical_ROI{f}}; cax=[0 4];
tiledlayout(3,2);
ax3=nexttile([1 1]);
imagesc(nTau{1},[1:nROI],squeeze(mean(STA_SSmat(dist_order,:,:),2)),cax)
title('Simple spike')
ax4=nexttile([1 1]);
imagesc(nTau{2},[1:nROI],squeeze(mean(STA_CSmat(dist_order,:,:),2)),cax)
title('Complex spike')
colormap(turbo); linkaxes([ax3 ax4],'xy')
ax1=[ax1 nexttile([1 1])];
h1=errorbar_shade(nTau{1},squeeze(mean(mean(STA_SSmat(rois{1},:,:),1),2))',squeeze(std(mean(STA_SSmat(rois{1},:,:),1),0,2))',cmap(1,:)); hold all
h2=errorbar_shade(nTau{1},squeeze(mean(mean(STA_SSmat(rois{2},:,:),1),2))',squeeze(std(mean(STA_SSmat(rois{2},:,:),1),0,2))',cmap(2,:));
legend([h1 h2],{'Basal','Apical'})
xlabel('Peri-spike time (ms)')
title('Simple spike')
ax1=[ax1 nexttile([1 1])];
h1=errorbar_shade(nTau{2},squeeze(mean(mean(STA_CSmat(rois{1},:,:),1),2))',squeeze(std(mean(STA_CSmat(rois{1},:,:),1),0,2))',cmap(1,:)); hold all
h2=errorbar_shade(nTau{2},squeeze(mean(mean(STA_CSmat(rois{2},:,:),1),2))',squeeze(std(mean(STA_CSmat(rois{2},:,:),1),0,2))',cmap(2,:));
legend([h1 h2],{'Basal','Apical'})
xlabel('Peri-spike time (ms)')
title('Complex spike')
ax1=[ax1 nexttile([1 1])];
h1=errorbar_shade(nTau{1},squeeze(mean(mean(STA_SSmat(rois{1},:,:),1),2))',squeeze(std(mean(STA_SSmat(rois{1},:,:),1),0,2))',cmap(3,:)); hold all
h2=errorbar_shade(nTau{2},squeeze(mean(mean(STA_CSmat(rois{1},:,:),1),2))',squeeze(std(mean(STA_CSmat(rois{1},:,:),1),0,2))',cmap(6,:));
legend([h1 h2],{'SS','CS'})
xlabel('Peri-spike time (ms)')
title('Basal')
ax1=[ax1 nexttile([1 1])];
h=errorbar_shade(nTau{1},squeeze(mean(mean(STA_SSmat(rois{2},:,:),1),2))',squeeze(std(mean(STA_SSmat(rois{2},:,:),1),0,2))',cmap(3,:)); hold all
h2=errorbar_shade(nTau{2},squeeze(mean(mean(STA_CSmat(rois{2},:,:),1),2))',squeeze(std(mean(STA_CSmat(rois{2},:,:),1),0,2))',cmap(6,:));
xlabel('Peri-spike time (ms)')
legend([h h2],{'SS','CS'})
linkaxes(ax1,'xy')
title('Apical')
%% SS, CS vs Non Spike dynamics
nTau_silent=[-50:200]; presp_time=[15 2]; Nbin=70;
all_spike=find(max(Result.spike,[],1));
SilentTrace=NomalizedTrace;
SilentTrace(:,all_spike'+nTau_silent)=NaN;
SilentTrace(:,Result.motionReject)=NaN;
SilentTrace=movmean(SilentTrace,range(presp_time),2);

figure(2); clf; cmap=distinguishable_colors(6);
tiledlayout(6,2)
nexttile([3 1])
plot(SilentTrace(1,:),mean(SilentTrace(rois{2},:),1),'.','color',[0.6 0.6 0.6]); hold all
SubSS_soma=squeeze(mean(STA_SSmat(1,:,-nTau{1}(1)-presp_time(1):-nTau{1}(1)-presp_time(2)),3));
SubSS_basal=squeeze(mean(STA_SSmat(rois{1},:,-nTau{1}(1)-presp_time(1):-nTau{1}(1)-presp_time(2)),[1 3]));
SubSS_apical=squeeze(mean(STA_SSmat(rois{2},:,-nTau{1}(1)-presp_time(1):-nTau{1}(1)-presp_time(2)),[1 3]));
plot(SubSS_soma,SubSS_apical,'.','color',cmap(1,:),'MarkerSize',15)
SubCS_soma=squeeze(mean(STA_CSmat(1,:,-nTau{2}(1)-presp_time(1):-nTau{2}(1)-presp_time(2)),3));
SubCS_basal= squeeze(mean(STA_CSmat(rois{1},:,-nTau{2}(1)-presp_time(1):-nTau{2}(1)-presp_time(2)),[1 3]));
SubCS_apical=squeeze(mean(STA_CSmat(rois{2},:,-nTau{2}(1)-presp_time(1):-nTau{2}(1)-presp_time(2)),[1 3]));
plot(SubCS_soma,SubCS_apical,'.','color',cmap(2,:),'MarkerSize',15)
SubdSP_soma=squeeze(mean(STA_dSPmat(1,:,-nTau{3}(1)-presp_time(1):-nTau{3}(1)-presp_time(2)),3));
SubdSP_basal=squeeze(mean(STA_dSPmat(rois{1},:,-nTau{3}(1)-presp_time(1):-nTau{3}(1)-presp_time(2)),[1 3]));
SubdSP_apical=squeeze(mean(STA_dSPmat(rois{2},:,-nTau{3}(1)-presp_time(1):-nTau{3}(1)-presp_time(2)),[1 3]));
plot(SubdSP_soma,SubdSP_apical,'.','color',cmap(3,:),'MarkerSize',15)
xlabel('V Soma')
ylabel('V Apical')
legend({'Non-Spiking','SS','CS','dSpike'})
nexttile([1 1])
h=histogram(SilentTrace(1,:),Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
histogram(SubSS_soma,h.BinEdges,'Normalization','probability','FaceColor',cmap(1,:))
xlabel('V Soma'); ylabel('Probability');
legend({'Non spike','SS'})
nexttile([1 1])
h=histogram(SilentTrace(1,:),Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
histogram(SubCS_soma,h.BinEdges,'Normalization','probability','FaceColor',cmap(2,:))
xlabel('V Soma'); ylabel('Probability');
legend({'Non spike','CS'})
nexttile([1 1])
h=histogram(SilentTrace(1,:),Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
histogram(SubdSP_soma,h.BinEdges,'Normalization','probability','FaceColor',cmap(3,:))
xlabel('V Soma'); ylabel('Probability')
legend({'Non spike','dSP'})

nexttile([1 1])
h=histogram(mean(SilentTrace(rois{1},:),1),Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
histogram(SubSS_basal,h.BinEdges,'Normalization','probability','FaceColor',cmap(1,:))
xlabel('V Basal'); ylabel('Probability');
legend({'Non spike','SS'})
nexttile(9,[1 1])
h=histogram(mean(SilentTrace(rois{1},:),1),Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
histogram(SubCS_basal,h.BinEdges,'Normalization','probability','FaceColor',cmap(2,:))
xlabel('V Basal'); ylabel('Probability');
legend({'Non spike','CS'})
nexttile(11,[1 1])
h=histogram(mean(SilentTrace(rois{1},:),1),Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
histogram(SubdSP_basal,h.BinEdges,'Normalization','probability','FaceColor',cmap(3,:))
xlabel('V Basal'); ylabel('Probability');
legend({'Non spike','dSP'})

nexttile([1 1])
h=histogram(mean(SilentTrace(rois{2},:),1),Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
histogram(SubSS_apical,h.BinEdges,'Normalization','probability','FaceColor',cmap(1,:))
xlabel('V Apical'); ylabel('Probability');
legend({'Non spike','SS'})
nexttile([1 1])
h=histogram(mean(SilentTrace(rois{2},:),1),Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
histogram(SubCS_apical,h.BinEdges,'Normalization','probability','FaceColor',cmap(2,:))
xlabel('V Apical'); ylabel('Probability');
legend({'Non spike','CS'})
nexttile([1 1])
h=histogram(mean(SilentTrace(rois{2},:),1),Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
histogram(SubdSP_apical,h.BinEdges,'Normalization','probability','FaceColor',cmap(3,:))
xlabel('V Apical'); ylabel('Probability');
legend({'Non spike','dSP'})

figure(6); clf;
nexttile([1 1])
ddat={SubSS_basal,SubCS_basal};
x=[ones(1,length(SubSS_basal))+randn(1,length(SubSS_basal))*0.05 ones(1,length(SubCS_basal))+randn(1,length(SubCS_basal))*0.05+1];
M=cellfun(@mean,ddat); S=cellfun(@std,ddat);
violin(ddat,'x',[1 2],'facecolor',cmap,'edgecolor',[],'mc',[],'medc',[],'bw',0.2,'plotlegend',0);
hold all
scatter(x,cell2mat(ddat),3,cmap([ones(1,length(SubSS_basal)) ones(1,length(SubCS_basal))+1],:))
errorbar([1 2],M,S,'LineWidth',2,'linestyle','none',...
    'color','k','Capsize',10,'marker','+','markersize',10)
xlim([0.5 2.5]); 
set(gca,'XTick',[1 2],'XTickLabel',{'SS','CS'})
ylabel('Subthreshold ')
title(['Basal, p = ' num2str(ranksum(SubSS_basal,SubCS_basal))]);

ddat={SubSS_apical,SubCS_apical};
nexttile([1 1])
x=[ones(1,length(SubSS_apical))+randn(1,length(SubSS_apical))*0.05 ones(1,length(SubCS_apical))+randn(1,length(SubCS_apical))*0.05+1];
M=cellfun(@mean,ddat); S=cellfun(@std,ddat);
violin(ddat,'x',[1 2],'facecolor',cmap,'edgecolor',[],'mc',[],'medc',[],'bw',0.2,'plotlegend',0);
hold all
scatter(x,cell2mat(ddat),3,cmap([ones(1,length(SubSS_apical)) ones(1,length(SubCS_apical))+1],:))
errorbar([1 2],M,S,'LineWidth',2,'linestyle','none',...
    'color','k','Capsize',10,'marker','+','markersize',10)
xlim([0.5 2.5]); 
set(gca,'XTick',[1 2],'XTickLabel',{'SS','CS'})
ylabel('Subthreshold')
title(['Apical, p = ' num2str(ranksum(SubSS_apical,SubCS_apical))]);

figure(3); clf; cmap=distinguishable_colors(6);
APSS_soma=squeeze(mean(STA_SSmat(1,:,-nTau{1}(1)+1),3));
APSS_apical=squeeze(mean(max(STA_SSmat(rois{2},:,-nTau{1}(1)+[1:3]),[],3),[1])); hold all
plot(APSS_soma,APSS_apical,'o','color',cmap(1,:),'MarkerSize',6)
APCS_soma=squeeze(mean(STA_CSmat(1,:,-nTau{2}(1)+1),3));
APCS_apical=squeeze(mean(max(STA_CSmat(rois{2},:,-nTau{2}(1)+[1:3]),[],3),[1])); hold all
plot(APCS_soma,APCS_apical,'o','color',cmap(2,:),'MarkerSize',6)
APdSP_soma=squeeze(mean(STA_dSPmat(1,:,-nTau{3}(1)+1),3));
APdSP_apical=squeeze(mean(STA_dSPmat(rois{2},:,-nTau{3}(1)+1),[1 3]));
plot(APdSP_soma,APdSP_apical,'o','color',cmap(3,:),'MarkerSize',6)
legend({'SS','CS','dSpike'})

%%
dt=11; dVdt=[]; nFrame=size(NomalizedTrace,2);
for t=1:nFrame
    for n=1:nROI
        if t<(dt/2) || nFrame-t < dt/2
            dVdt(n,t)=NaN;
        else
            p=polyfit([1:dt],NomalizedTrace(n,t+[-floor(dt/2):floor(dt/2)]),1);
            dVdt(n,t)=p(1);
        end
    end
end
dVdt_soma=[]; dVdt_basal=[]; dVdt_apical=[];
for t=1:nFrame
        if t<(dt/2) || nFrame-t < dt/2
            dVdt_soma(1,t)=NaN;
            dVdt_basal(1,t)=NaN;
            dVdt_apical(1,t)=NaN;
        else
            p1=polyfit([1:dt],NomalizedTrace(1,t+[-floor(dt/2):floor(dt/2)]),1);
            dVdt_soma(1,t)=p1(1);
            p2=polyfit([1:dt],mean(NomalizedTrace(rois{1},t+[-floor(dt/2):floor(dt/2)])),1);
            dVdt_basal(1,t)=p2(1);
            p3=polyfit([1:dt],mean(NomalizedTrace(rois{2},t+[-floor(dt/2):floor(dt/2)])),1);
            dVdt_apical(1,t)=p3(1);
        end
end
dvdt_region=[dVdt_soma;dVdt_basal;dVdt_apical];
%%
dVdt_CSmat=reshape(dvdt_region(:,CS_s'+nTau{2}),3,[],length(nTau{2}));
dVdt_SSmat=reshape(dvdt_region(:,bAP_s'+nTau{1}),3,[],length(nTau{1}));
dVdt_dSPmat=reshape(dvdt_region(:,dSP_s'+nTau{3}),3,[],length(nTau{3}));

nTau_silent=[-50:200]; presp_time=[25 10]; Nbin=70;
all_spike=find(max(Result.spike,[],1));
SilentTrace=dvdt_region;
SilentTrace(:,all_spike'+nTau_silent)=NaN;
SilentTrace(:,Result.motionReject)=NaN;
SilentTrace=movmean(SilentTrace,range(presp_time),2);

figure(2); clf; cmap=distinguishable_colors(6);
tiledlayout(6,2)
nexttile([3 1])
plot(SilentTrace(1,:),mean(SilentTrace(3,:),1),'.','color',[0.6 0.6 0.6]); hold all
SubdVdtSS_soma=squeeze(mean(dVdt_SSmat(1,:,-nTau{1}(1)-presp_time(1):-nTau{1}(1)-presp_time(2)),3));
SubdVdtSS_basal=squeeze(mean(dVdt_SSmat(2,:,-nTau{1}(1)-presp_time(1):-nTau{1}(1)-presp_time(2)),[1 3]));
SubdVdtSS_apical=squeeze(mean(dVdt_SSmat(3,:,-nTau{1}(1)-presp_time(1):-nTau{1}(1)-presp_time(2)),[1 3]));
plot(SubdVdtSS_soma,SubdVdtSS_apical,'.','color',cmap(1,:),'MarkerSize',15)

SubdVdtCS_soma=squeeze(mean(dVdt_CSmat(1,:,-nTau{2}(1)-presp_time(1):-nTau{2}(1)-presp_time(2)),3));
SubdVdtCS_basal= squeeze(mean(dVdt_CSmat(2,:,-nTau{2}(1)-presp_time(1):-nTau{2}(1)-presp_time(2)),[1 3]));
SubdVdtCS_apical=squeeze(mean(dVdt_CSmat(3,:,-nTau{2}(1)-presp_time(1):-nTau{2}(1)-presp_time(2)),[1 3]));
plot(SubdVdtCS_soma,SubdVdtCS_apical,'.','color',cmap(2,:),'MarkerSize',15)

SubdVdtdSP_soma=squeeze(mean(dVdt_dSPmat(1,:,-nTau{3}(1)-presp_time(1):-nTau{3}(1)-presp_time(2)),3));
SubdVdtdSP_basal=squeeze(mean(dVdt_dSPmat(2,:,-nTau{3}(1)-presp_time(1):-nTau{3}(1)-presp_time(2)),[1 3]));
SubdVdtdSP_apical=squeeze(mean(dVdt_dSPmat(3,:,-nTau{3}(1)-presp_time(1):-nTau{3}(1)-presp_time(2)),[1 3]));
plot(SubdVdtdSP_soma,SubdVdtdSP_apical,'.','color',cmap(3,:),'MarkerSize',15)
xlabel('dV/dt Soma')
ylabel('dV/dt Apical')
legend({'Non-Spiking','SS','CS','dSpike'})

nexttile([1 1])
h=histogram(SilentTrace(1,:),Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
histogram(SubdVdtSS_soma,h.BinEdges,'Normalization','probability','FaceColor',cmap(1,:))
xlabel('dV/dt Soma'); ylabel('Probability');
legend({'Non spike','SS'})
nexttile([1 1])
h=histogram(SilentTrace(1,:),Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
histogram(SubdVdtCS_soma,h.BinEdges,'Normalization','probability','FaceColor',cmap(2,:))
xlabel('dV/dt Soma'); ylabel('Probability');
legend({'Non spike','CS'})
nexttile([1 1])
h=histogram(SilentTrace(1,:),Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
histogram(SubdVdtdSP_soma,h.BinEdges,'Normalization','probability','FaceColor',cmap(3,:))
xlabel('dV/dt Soma'); ylabel('Probability')
legend({'Non spike','dSP'})

nexttile([1 1])
h=histogram(mean(SilentTrace(2,:),1),Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
histogram(SubdVdtSS_basal,h.BinEdges,'Normalization','probability','FaceColor',cmap(1,:))
xlabel('dV/dt Basal'); ylabel('Probability');
legend({'Non spike','SS'})
nexttile(9,[1 1])
h=histogram(mean(SilentTrace(2,:),1),Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
histogram(SubdVdtCS_basal,h.BinEdges,'Normalization','probability','FaceColor',cmap(2,:))
xlabel('dV/dt Basal'); ylabel('Probability');
legend({'Non spike','CS'})
nexttile(11,[1 1])
h=histogram(mean(SilentTrace(2,:),1),Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
histogram(SubdVdtdSP_basal,h.BinEdges,'Normalization','probability','FaceColor',cmap(3,:))
xlabel('dV/dt Basal'); ylabel('Probability');
legend({'Non spike','dSP'})

nexttile([1 1])
h=histogram(mean(SilentTrace(3,:),1),Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
histogram(SubdVdtSS_apical,h.BinEdges,'Normalization','probability','FaceColor',cmap(1,:))
xlabel('dV/dt Apical'); ylabel('Probability');
legend({'Non spike','SS'})
nexttile([1 1])
h=histogram(mean(SilentTrace(3,:),1),Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
histogram(SubdVdtCS_apical,h.BinEdges,'Normalization','probability','FaceColor',cmap(2,:))
xlabel('dV/dt Apical'); ylabel('Probability');
legend({'Non spike','CS'})
nexttile([1 1])
h=histogram(mean(SilentTrace(3,:),1),Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
histogram(SubdVdtdSP_apical,h.BinEdges,'Normalization','probability','FaceColor',cmap(3,:))
xlabel('dV/dt Apical'); ylabel('Probability');
legend({'Non spike','dSP'})

figure(6); clf;
nexttile([1 1])
ddat={SubdVdtSS_basal,SubdVdtCS_basal};
x=[ones(1,length(SubSS_basal))+randn(1,length(SubSS_basal))*0.05 ones(1,length(SubCS_basal))+randn(1,length(SubCS_basal))*0.05+1];
M=cellfun(@mean,ddat); S=cellfun(@std,ddat);
violin(ddat,'x',[1 2],'facecolor',cmap,'edgecolor',[],'mc',[],'medc',[],'bw',0.02,'plotlegend',0);
hold all
scatter(x,cell2mat(ddat),3,cmap([ones(1,length(SubSS_basal)) ones(1,length(SubCS_basal))+1],:))
errorbar([1 2],M,S,'LineWidth',2,'linestyle','none',...
    'color','k','Capsize',10,'marker','+','markersize',10)
xlim([0.5 2.5]); 
set(gca,'XTick',[1 2],'XTickLabel',{'SS','CS'})
ylabel('Subthreshold ')
title(['Basal, p = ' num2str(ranksum(SubdVdtSS_basal,SubdVdtCS_basal))]);

ddat={SubdVdtSS_apical,SubdVdtCS_apical};
nexttile([1 1])
x=[ones(1,length(SubSS_apical))+randn(1,length(SubSS_apical))*0.05 ones(1,length(SubCS_apical))+randn(1,length(SubCS_apical))*0.05+1];
M=cellfun(@mean,ddat); S=cellfun(@std,ddat);
violin(ddat,'x',[1 2],'facecolor',cmap,'edgecolor',[],'mc',[],'medc',[],'bw',0.02,'plotlegend',0);
hold all
scatter(x,cell2mat(ddat),3,cmap([ones(1,length(SubSS_apical)) ones(1,length(SubCS_apical))+1],:))
errorbar([1 2],M,S,'LineWidth',2,'linestyle','none',...
    'color','k','Capsize',10,'marker','+','markersize',10)
xlim([0.5 2.5]); 
set(gca,'XTick',[1 2],'XTickLabel',{'SS','CS'})
ylabel('dV/dt')
title(['Apical, p = ' num2str(ranksum(SubdVdtSS_apical,SubdVdtCS_apical))]);

%% Plateau vs CSs
plateauCS=[8 9 17];
%RestCS=setdiff([1:size(STA_CSmat,2)],plateauCS);
RestCS=setdiff([10:60],plateauCS);
STA_rest=squeeze(mean(STA_CSmat(:,RestCS,:),2));
% 
% figure(7); clf;
% t_show=[-nTau{2}(1)-100:-nTau{2}(1)+200];
% colorline(STA_rest(1,t_show)',mean(STA_rest(rois{2},t_show))','winter','o'); hold all
% colorline(squeeze(STA_CSmat(1,9,t_show)),squeeze(mean(STA_CSmat(rois{2},9,t_show))),'turbo','*'); hold all

f=figure(7); cmap=distinguishable_colors(6);
set(f, 'Color', 'w');
myVideo = VideoWriter(['20240613_PlateauPhase'],"MPEG-4"); %open video file
myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
myVideo.Quality= 100;
open(myVideo)
f.Position = [100, 100, 600, 700];
for t=51:250
    clf;
colorline(STA_rest(1,50:t)',mean(STA_rest(rois{2},50:t))','winter','o'); hold all    
xlabel('V Soma'); ylabel('V apical');
title(['Average of CSs, ' num2str(t+nTau{2}(1)) ' ms'])
xlim([-4 8.5]); ylim([-4 7.5]);
frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    pause(0.05)
end

for t=2:350
    clf;
plot(STA_rest(1,51:250)',mean(STA_rest(rois{2},51:250))','color',[0.5 0.5 0.5]); hold all
colorline(squeeze(STA_CSmat(1,9,1:t)),squeeze(mean(STA_CSmat(rois{2},9,1:t))),'turbo','*'); hold all
xlabel('V Soma'); ylabel('V apical');
xlim([-4 8.5]); ylim([-4 7.5]);
title(['Plateau potential, ' num2str(t) ' ms'])
frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    pause(0.05)
end
close(myVideo);
close(f);
%%
nTau_silent=[-50:200]; presp_time=[30 2]; Nbin=50;
SubCS_soma=squeeze(mean(STA_CSmat(1,:,-nTau{2}(1)-presp_time(1):-nTau{2}(1)-presp_time(2)),3));
SubCS_basal= squeeze(mean(STA_CSmat(rois{1},:,-nTau{2}(1)-presp_time(1):-nTau{2}(1)-presp_time(2)),[1 3]));
SubCS_apical=squeeze(mean(STA_CSmat(rois{2},:,-nTau{2}(1)-presp_time(1):-nTau{2}(1)-presp_time(2)),[1 3]));
SubdVdtCS_soma=squeeze(mean(dVdt_CSmat(1,:,-nTau{2}(1)-presp_time(1):-nTau{2}(1)-presp_time(2)),3));
SubdVdtCS_basal= squeeze(mean(dVdt_CSmat(2,:,-nTau{2}(1)-presp_time(1):-nTau{2}(1)-presp_time(2)),[1 3]));
SubdVdtCS_apical=squeeze(mean(dVdt_CSmat(3,:,-nTau{2}(1)-presp_time(1):-nTau{2}(1)-presp_time(2)),[1 3]));

Subplateau= [mean(STA_CSmat(1,9,-nTau{2}(1)-presp_time(1):-nTau{2}(1)-presp_time(2)),3) ...
    squeeze(mean(STA_CSmat(rois{1},9,-nTau{2}(1)-presp_time(1):-nTau{2}(1)-presp_time(2)),[1 3])) ...
    squeeze(mean(STA_CSmat(rois{2},9,-nTau{2}(1)-presp_time(1):-nTau{2}(1)-presp_time(2)),[1 3]))];

SubdVdtplateau= [squeeze(mean(dVdt_CSmat(1,:,-nTau{2}(1)-presp_time(1):-nTau{2}(1)-presp_time(2)),3)) ...
    squeeze(mean(dVdt_CSmat(2,:,-nTau{2}(1)-presp_time(1):-nTau{2}(1)-presp_time(2)),[1 3])) ...
    squeeze(mean(dVdt_CSmat(3,:,-nTau{2}(1)-presp_time(1):-nTau{2}(1)-presp_time(2)),[1 3]))];

figure(8); clf;
tiledlayout(3,2)
nexttile(1,[1 1])
h=histogram(SubCS_soma,Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
plot([Subplateau(1) Subplateau(1)],[0 0.1],'r')
xlabel('V soma'); ylabel('Probability');
nexttile(3,[1 1])
h=histogram(SubCS_basal,Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
plot([Subplateau(2) Subplateau(2)],[0 0.1],'r')
xlabel('V basal'); ylabel('Probability');
nexttile(5,[1 1])
h=histogram(SubCS_apical,Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
plot([Subplateau(3) Subplateau(3)],[0 0.1],'r')
xlabel('V apical'); ylabel('Probability');

nexttile([1 1])
h=histogram(SubdVdtCS_soma,Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
plot([SubdVdtplateau(1) SubdVdtplateau(1)],[0 0.1],'r')
xlabel('dV/dt Soma'); ylabel('Probability');

nexttile([1 1])
h=histogram(SubdVdtCS_basal,Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
plot([SubdVdtplateau(2) SubdVdtplateau(2)],[0 0.1],'r')
xlabel('dV/dt basal'); ylabel('Probability');

nexttile([1 1])
h=histogram(SubdVdtCS_apical,Nbin,'Normalization','probability','FaceColor',[0.6 0.6 0.6]); hold all
plot([SubdVdtplateau(3) SubdVdtplateau(3)],[0 0.1],'r')
xlabel('dV/dt apical'); ylabel('Probability');

%% plateau movie

load(fullfile(fpath{f},"output_data.mat"))
load([fpath{f} '/mcTrace' num2str(2,'%02d') '.mat']);
sz=double(Device_Data{1, 3}.ROI([2 4]));
mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(2,'%02d') '.bin'],sz(2),sz(1)));

mc=mcTrace.xymean;
meanF=squeeze(mean(mov_mc,[1 2]));
meanF2=-(meanF-movmedian(meanF,100));
sp=find_spike_bh(meanF2'./get_threshold(meanF2',1),4,2);
SilentPeriod=ones(1,size(mov_mc,3));
SilentPeriod(find(max(sp,[],1))'+[-10:150])=NaN;
t_fit=find(~isnan(SilentPeriod));
[y_fit t_consts coeffY]  = expfitDM_2(t_fit',meanF2(t_fit),[1:size(mov_mc,3)]',1000);

bkg = zeros(1, size(mov_mc,3));
bkg(1,:) = y_fit;  % linear term

mov_res= mov_mc-median(mov_mc,3);
mov_res2 = SeeResiduals(mov_res,mc);
mov_res2 = SeeResiduals(mov_res2,mc.^2);
mov_res2 = SeeResiduals(mov_res2,mc(:,1).*mc(:,end));
mov_res2= SeeResiduals(mov_res2,bkg,1);

t_show=10595+[-500:1000];
bound=7;
mov_filt=pcafilt(imgaussfilt3(-mov_res2(bound:end-bound,bound:end-bound,t_show),[2 2 0.1]),10);

Rfixed = imref2d([size(Result.ref_im,1) size(Result.ref_im,2)]);
inverseTform = invert(Result.tform);
revertedStruct = imwarp(Result.Struckmax, inverseTform,'OutputView',Rfixed);
revertedStruct(revertedStruct==0)=prctile(revertedStruct(:),30);
revertedStruct=mat2gray(revertedStruct);
revertedStruct=revertedStruct(bound:end-bound,bound:end-bound);
revertedStruct_filt=revertedStruct-imgaussfilt(revertedStruct,15);
revertedStruct(revertedStruct==1)=0; revertedStruct(revertedStruct<0.05)=0;

colormovie=grs2rgb(tovec(mov_filt.*revertedStruct),colormap('jet'),-3,6);
colormovie=reshape(colormovie,size(mov_filt,1),size(mov_filt,2),size(mov_filt,3),[]);
colormovie=permute(colormovie,[1 2 4 3]);
colormovie=colormovie.*revertedStruct*5;%.*SkelDend(bound:end-bound,bound:end-bound);

    mov_show=imrotate(colormovie,-14);
    mov_show=mov_show(56:180,22:350,:,:);
    writeMov4d(['Plateau'],mov_show,[1:size(mov_filt,3)],10,1,[-5 10])
%%
f=figure(10); cmap=distinguishable_colors(6);
set(f, 'Color', 'w');
myVideo = VideoWriter(['20240604_dvdtTrace'],"MPEG-4"); %open video file
myVideo.FrameRate = 50;  %can adjust this, 5 - 10 works well for me
myVideo.Quality= 100;
open(myVideo)
T_seg=200; Tau=[-50 50];
clf;
f.Position = [100, 100, 600, 700];
tiledlayout(5,3);
ax11=nexttile([2 3]);
ax12=[nexttile([1 3])];
linkaxes([ax11 ax12 ax13],'x')
ax3=nexttile([1 1]);
ax4=nexttile([1 1]);
ax5=nexttile([1 1]);
ax6=nexttile([1 1]);
ax7=nexttile([1 1]);
ax8=nexttile([1 1]);

for i=1:nFrame-max(T_seg,Tau(2))
    if i<T_seg
        t_show=[1:i];
    else
        t_show=[i-T_seg+1:i];
    end

    if i+Tau(1)<1
        t_show2=[1:i+Tau(2)-1];
        t_x=t_show2-i;
    else
        t_show2=i+[Tau(1):Tau(2)];
        t_x=t_show2-median(t_show2);
    end
    axes(ax11);
    plot(t_x,NomalizedTrace(1,t_show2)); hold all
    plot(t_x,mean(NomalizedTrace(rois{1},t_show2),1))
    plot(t_x,mean(NomalizedTrace(rois{2},t_show2),1))
    plot([0 0],[-2 10],'color',[1 0 0],'LineStyle','--'); hold off
    title(num2str(i/1000))
    set([ax11],'YLim',[-4 12])
    legend({'Soma','Basal','Apical'},'Location','northeastoutside')
    
    axes(ax12);
    plot(t_x,Result.SpClass(:,t_show2)');
    ylim([0 1.2])
    yyaxis right
    plot(t_x,Result.VR(5,t_show2)); hold all
    plot([0 0],[0 100],'color',[1 0 0],'LineStyle','--'); hold off
    ylabel('VR position')
    set([ax12],'XLim',Tau,'YLim',[0 115])
    legend({'SS','CS','dSp'},'Location','northeastoutside')
    axes(ax3);
    plot(mean(NomalizedTrace(rois{1},t_show)),mean(dVdt(rois{2},t_show))); hold all
    plot(mean(NomalizedTrace(rois{1},t_show(end))),mean(dVdt(rois{2},t_show(end))),'ro'); hold off
    xlabel('Basal V')
    ylabel('Apical dV/dt')
    axes(ax4);
    plot(mean(NomalizedTrace(rois{1},t_show)),mean(dVdt(rois{1},t_show))); hold all
    plot(mean(NomalizedTrace(rois{1},t_show(end))),mean(dVdt(rois{1},t_show(end))),'ro'); hold off
    xlabel('Basal V')
    ylabel('Basal dV/dt')
    axes(ax5);
    plot(mean(NomalizedTrace(rois{2},t_show)),mean(dVdt(rois{2},t_show))); hold all
    plot(mean(NomalizedTrace(rois{2},t_show(end))),mean(dVdt(rois{2},t_show(end))),'ro'); hold off
    xlabel('Apical V')
    ylabel('Apical dV/dt')
    axes(ax6);
    plot(mean(NomalizedTrace(rois{2},t_show)),mean(dVdt(rois{1},t_show))); hold all
    plot(mean(NomalizedTrace(rois{2},t_show(end))),mean(dVdt(rois{1},t_show(end))),'ro'); hold off
    xlabel('Apical V')
    ylabel('Basal dV/dt')
    axes(ax7);
    plot(mean(dVdt(rois{1},t_show)),mean(dVdt(rois{2},t_show))); hold all
    plot(mean(dVdt(rois{1},t_show(end))),mean(dVdt(rois{2},t_show(end))),'ro'); hold off
    xlabel('Basal dV/dt')
    ylabel('Apical dV/dt')
    axes(ax8);
    plot(mean(NomalizedTrace(rois{1},t_show)),mean(NomalizedTrace(rois{2},t_show))); hold all
    plot(mean(NomalizedTrace(rois{1},t_show(end))),mean(NomalizedTrace(rois{2},t_show(end))),'ro'); hold off
    xlabel('Basal V')
    ylabel('Apical V')

    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    pause(0.05)
end
close(myVideo);
close(f);


%%
Nbin=50;
STACSbasal=mean(STA_CSmat(rois{1},:,[1:-nTau{2}(1)]),[1 2]); STACSapical=mean(STA_CSmat(rois{2},:,[1:-nTau{2}(1)]),[1 2]);
edgeX=[min(STACSbasal):range(STACSbasal)/Nbin:max(STACSbasal)]*4;
edgeY=[min(STACSapical):range(STACSapical)/Nbin:max(STACSapical)]*4;
[CSsta_subCounts] = histcounts2(STACSbasal,STACSapical,edgeX,edgeY);
[SSsta_subCounts] = histcounts2(mean(STA_SSmat(rois{1},:,[1:-nTau{1}(1)]),[1 2]),mean(STA_SSmat(rois{2},:,[1:-nTau{1}(1)]),[1 2]),edgeX,edgeY);
SS_subCounts=[]; CS_subCounts=[];
for ss=1:size(STA_SSmat,2)
    SS_subCounts(:,:,ss) = histcounts2(squeeze(mean(STA_SSmat(rois{1},ss,[1:-nTau{1}(1)]),1)),squeeze(mean(STA_SSmat(rois{2},ss,[1:-nTau{1}(1)]),1)),edgeX,edgeY);
end
for cs=1:size(STA_CSmat,2)
    CS_subCounts(:,:,cs) = histcounts2(squeeze(mean(STA_CSmat(rois{1},cs,[1:-nTau{2}(1)]),1)),squeeze(mean(STA_CSmat(rois{2},cs,[1:-nTau{2}(1)]),1)),edgeX,edgeY);
end

figure(3); clf;
nexttile([1 1])
colorline(movmean(squeeze(mean(STA_CSmat(rois{1},:,[1:-nTau{2}(1)+1]),[1 2])),1),movmean(squeeze(mean(STA_CSmat(rois{2},:,[1:-nTau{2}(1)+1]),[1 2])),1),'jet','*')
hold all
colorline(movmean(squeeze(mean(STA_SSmat(rois{1},:,[1:-nTau{1}(1)+1]),[1 2])),1),movmean(squeeze(mean(STA_SSmat(rois{2},:,[1:-nTau{1}(1)+1]),[1 2])),1),'jet','o')
xlabel('Basal V')
ylabel('Apical V')
nexttile([1 1])

