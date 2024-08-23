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
rois={basal_ROI{f},apical_ROI{f}};
nROI=size(Result.normTraces,1);
nTau={[-55:15],[-60:200],[-20:20]}; %SS, CS, dSP
nTau_db=[-60:200];
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

dt=5; %ms
% CS_s=[];
% C_list=find(Result.SpClass(2,:));
% CS_label=max(Result.spike,[],1).*bwlabel(Result.CStrace);
% CS_list=[];
% for g=1:length(C_list)
%     s=C_list(g);
%     s_tmp=Result.spike(1,:);
%     s_tmp(find(CS_label==g))=0;
%     CS_s=[CS_s s];
%     CS_list=[CS_list g];
% end

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

doublet_ISI_threshold= 20; %ms
somSpike=find(Result.SpClass(1,:));
somSpike_label=zeros(1,length(somSpike));
ISI_label=bwlabel((somSpike(2:end)-somSpike(1:end-1))<doublet_ISI_threshold);
Sp_Class=zeros(1,size(Result.normTraces,2));
Sp_Class_doubletTrace=zeros(1,size(Result.normTraces,2));
for b=1:max(ISI_label)
    spikeind=find(ISI_label==b);
    Sp_Class(somSpike(spikeind(1)))=2; %doublet
    Sp_Class_doubletTrace(somSpike([spikeind spikeind(end)+1]))=b;
    somSpike_label([spikeind spikeind(end)+1])=b;
    doublet_length(b)=length(spikeind)+1;
end
Sp_Class(somSpike(find(somSpike_label==0)))=1; %siglet
Sp_Class(find(Result.SpClass(2,:)))=3; %Complex

dblet_S=find(Sp_Class==2); %doublet
%dblet_S=dblet_S(find(doublet_length==3));
rmv_list=[];
siglet_spike=find(Sp_Class==1);
for g=1:length(dblet_S)
    s=dblet_S(g);
    isnearby=sum(ismember(s+nTau_db,siglet_spike))>1;
    isnearbyCS=sum(ismember(s+nTau_db,find(Result.CStrace)))>1;
    ispartCS=tr_trace(s)>0;
    if ~isnearby & ~isnan(s) & ~isnearbyCS & ~ispartCS
        rmv_list=[rmv_list g];
    end
end
dblet_S(rmv_list)=[];

bwCS=bwlabel(Result.CStrace);
CSstats = regionprops(bwCS, 'all');
CSstats = CSstats(CS_list);
STA_SS=squeeze(mean(reshape(Result.normTraces(:,bAP_ref'+nTau_bAP),nROI,[],length(nTau_bAP)),2));
STA_SS= STA_SS - mean(STA_SS(:,1:10),2);
F_ref=mean(STA_SS(:,-nTau_bAP(1)+[7:14]),2);
CSpike=Result.spike.*Result.CStrace;

%NomalizedTrace=(Result.normTraces-movprc(Result.normTraces,2000,35,2))./F_ref;
NomalizedTrace=(Result.normTraces)./F_ref;
NomalizedTrace_omitSpike=NomalizedTrace;
NomalizedTrace_omitSpike(:,find(max(Result.spike,[],1))'+[-5:5])=NaN;

STA_CSmat=reshape(NomalizedTrace(:,CS_s'+nTau{2}),nROI,[],length(nTau{2}));
STA_CSmat_omitSpike=reshape(NomalizedTrace_omitSpike(:,CS_s'+nTau{2}),nROI,[],length(nTau{2}));
STA_CSmat_sub=STA_CSmat;%-median(STA_CSmat_omitSpike,3,'omitnan');
STA_CSpikemat=reshape(CSpike(:,CS_s'+nTau{2}),nROI,[],length(nTau{2}));

STA_dbSmat=reshape(NomalizedTrace(:,dblet_S'+nTau_db),nROI,[],length(nTau_db));
STA_dbSmat_omitNaN=reshape(NomalizedTrace_omitSpike(:,dblet_S'+nTau_db),nROI,[],length(nTau_db));
STA_dbSmat=STA_dbSmat;%-median(STA_dbSmat_omitNaN,3,'omitnan');

STA_SSmat=reshape(NomalizedTrace(:,bAP_s'+nTau{1}),nROI,[],length(nTau{1}));

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
    CSprop.SpikeAmp{c}=squeeze(max(permute(reshape(STA_CSmat_sub(:,c,sp_t+[-1:3]),size(STA_CSmat,1),length(sp_t),5),[1 3 2]),[],2)); %roi, time, spike -> roi, spike amp
    % for t=1:-nTau{2}(1)/dt
    %     for n=1:size(STA_CSmat,1)
    %         tfit=[t*dt:(t+1)*dt]-2;
    %         [fitresult, gof] = fit(tfit', squeeze(STA_CSmat_sub(n,c,tfit)),ft);
    %         CSprop.dVdt(n,c,t)=fitresult.p1;
    %     end
    % end

    % for t=1:4
    %     for n=1:size(STA_CSmat,1)
    %         tfit_post=-nTau{2}(1)+[1:3]+t;
    %         [fitresult, gof] = fit([tfit_post]', squeeze(STA_CSmat_sub(n,c,tfit_post)),ft);
    %         CSprop.dVdt_post(n,c,t)=fitresult.p1;
    %     end
    % end

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

dbSprop=[];
for c=1:length(dblet_S)
    spikeind=find(Sp_Class_doubletTrace==c);
    dbSprop.nSpike(c)=length(spikeind);
    dbSprop.SpikeAmp{c}=squeeze(max(permute(reshape(NomalizedTrace(:,spikeind'+[-1:3]),nROI,[],5),[1 3 2]),[],2)); %roi, time, spike -> roi, spike amp
end
%% Show STAs
figure(2); clf; axh=[];
cax=[-0.5 2.5];
tiledlayout(2,3)
axh=[axh nexttile([1 1])];
imagesc(nTau{1},[1:nROI],squeeze(mean(STA_SSmat(dist_order,:,:),2)),cax); hold all
plot(nTau{1}(1)+1,som_roi,'marker','>','MarkerFaceColor','r')
ylabel('ROI (distance order)')
title('STA of SS')
axh=[axh nexttile([1 1])];
imagesc(nTau_db,[1:nROI],squeeze(mean(STA_dbSmat(dist_order,:,:),2)),cax); hold all
plot(nTau_db(1)+1,som_roi,'marker','>','MarkerFaceColor','r')
ylabel('ROI (distance order)')
title('STA of Burst')
axh=[axh nexttile([1 1])];
imagesc(nTau{2},[1:nROI],squeeze(mean(STA_CSmat(dist_order,:,:),2)),cax); hold all
plot(nTau{2}(1)+1,som_roi,'marker','>','MarkerFaceColor','r')
ylabel('ROI (distance order)')
title('STA of CS')
colormap(turbo)
linkaxes(axh,'x')
nexttile([1 1])
l=plot(nTau{1},squeeze(mean(STA_SSmat(dist_order,:,:),2))');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2))
xlabel('Peri-spike time (ms)')
nexttile([1 1])
l=plot(nTau_db,squeeze(mean(STA_dbSmat(dist_order,:,:),2))');
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
tiledlayout(1,3)
ax1=[ax1 nexttile([1 1])];
errorbar_shade(nTau{1},squeeze(mean(mean(STA_SSmat(rois{1},:,:),1),2))',squeeze(std(mean(STA_SSmat(rois{1},:,:),1),0,2))',cmap(1,:)); hold all
errorbar_shade(nTau{1},squeeze(mean(mean(STA_SSmat(rois{2},:,:),1),2))',squeeze(std(mean(STA_SSmat(rois{2},:,:),1),0,2))',cmap(2,:));
xlabel('Peri-spike time (ms)')
ax1=[ax1 nexttile([1 1])];
errorbar_shade(nTau_db,squeeze(mean(mean(STA_dbSmat(rois{1},:,:),1),2))',squeeze(std(mean(STA_CSmat(rois{1},:,:),1),0,2))',cmap(1,:)); hold all
errorbar_shade(nTau_db,squeeze(mean(mean(STA_dbSmat(rois{2},:,:),1),2))',squeeze(std(mean(STA_CSmat(rois{2},:,:),1),0,2))',cmap(2,:));
xlabel('Peri-spike time (ms)')
ax1=[ax1 nexttile([1 1])];
errorbar_shade(nTau{2},squeeze(mean(mean(STA_CSmat(rois{1},:,:),1),2))',squeeze(std(mean(STA_CSmat(rois{1},:,:),1),0,2))',cmap(1,:)); hold all
errorbar_shade(nTau{2},squeeze(mean(mean(STA_CSmat(rois{2},:,:),1),2))',squeeze(std(mean(STA_CSmat(rois{2},:,:),1),0,2))',cmap(2,:));
xlabel('Peri-spike time (ms)')
linkaxes(ax1,'xy')
%%

figure(2); clf;
length_dist=[ones(1,sum(somSpike_label==0)) doublet_length];
histogram(length_dist,[0.5:1:10]);
hold all
histogram(CSprop.nSpike,[0.5:1:10]);
xlabel('Number of spikes')
ylabel('Event #')

figure(3); clf; cax=[-3 15]; ax1=[]; ax2=[];
for s=1:40
    ax1=[ax1 nexttile([1 1])];
    imagesc(nTau_db,[1:nROI],squeeze(STA_dbSmat(dist_order,s,:)),cax)
end
colormap(turbo)

figure(4); clf;
for s=1:3:120
    ax2=[ax2 nexttile([1 1])];
    imagesc(nTau{2},[1:nROI],squeeze(STA_CSmat_sub(dist_order,s,:)),cax)
end
linkaxes(ax1,'x')
linkaxes(ax2,'x')
colormap(turbo)

%% Show phase plot
dblet_S=find(Sp_Class==2); %doublet
STA_CSmat_sub2=STA_CSmat_sub;
STA_dbSmat=reshape(NomalizedTrace(:,dblet_S'+nTau_db),nROI,[],length(nTau_db));
STA_dbSmat=STA_dbSmat-median(reshape(NomalizedTrace_omitSpike(:,dblet_S'+nTau_db),nROI,[],length(nTau_db)),3,'omitnan');
STA_dbSmat2=STA_dbSmat;
show_cs=[13 15 18 20 24 36 50 84 88]; ax_c=[]; ax1=[]; cax=[-3 10];
show_ds=[4 14 32 42 44 46 50 53 58]; ax_d=[]; cmap=distinguishable_colors(6);
figure(4); clf;
tiledlayout(3,6)
for s=show_cs
    ax_c=[ax_c nexttile([1 1])];
sp_tr=max(Result.spike(:,CS_s(s)+nTau{2}),[],1);
sp_tr(1:-nTau{2}(1))=0;
sub_spike_ind=bwlabel_firstind(sp_tr);
if sub_spike_ind(end)+20 < length(nTau{2})
    sub_spike_ind_seg=[[[1; sub_spike_ind(1:end-1)] sub_spike_ind]; [sub_spike_ind(end) sub_spike_ind(end)+20]];
else
    sub_spike_ind_seg=[[[1; sub_spike_ind(1:end-1)] sub_spike_ind]; [sub_spike_ind(end) length(nTau{2})]];
end
cmap_sub=turbo(size(sub_spike_ind_seg,1));
lg_label=[];
for s_sub=1:size(sub_spike_ind_seg,1)
    show_seg=[sub_spike_ind_seg(s_sub,1):sub_spike_ind_seg(s_sub,2)];
    soma_seg=squeeze(STA_CSmat_sub2(1,s,show_seg));
    apical_seg=squeeze(mean(STA_CSmat_sub2(rois{2},s,show_seg),1));
    plot(soma_seg,apical_seg,'color',cmap_sub(s_sub,:),'Linewidth',1.5); hold all
    lg_label{s_sub}=['Spike#' num2str(s_sub-1)];
end
plot(squeeze(STA_CSmat_sub2(1,s,sub_spike_ind)),squeeze(mean(STA_CSmat_sub2(rois{2},s,sub_spike_ind))),'r*')
text(squeeze(STA_CSmat_sub2(1,s,sub_spike_ind))+1,squeeze(mean(STA_CSmat_sub2(rois{2},s,sub_spike_ind))),num2str([1:size(sub_spike_ind_seg,1)-1]'))
lg_label(1)={'Pre spike'};
lg_label(end)={'Post spike'};
legend(lg_label,'Location','northwest','FontSize',8);
xlabel('V Soma')
ylabel('V Apical')

ax1=[ax1 nexttile([1 1])];
%imagesc(nTau{2},[1:nROI],squeeze(STA_CSmat_sub2(dist_order,s,:)),cax);
plot(nTau{2},squeeze(STA_CSmat_sub2(1,s,:)),'color',cmap(1,:)); hold all
plot(nTau{2},squeeze(mean(STA_CSmat_sub2(rois{2},s,:),1)),'color',cmap(2,:)); hold all
xlabel('Time (ms)')
end

figure(5); clf;
tiledlayout(3,6)

for s=show_ds
    ax_d=[ax_d nexttile([1 1])];
sp_tr=max(Result.spike(:,dblet_S(s)+nTau_db),[],1);
sp_tr(1:-nTau_db(1))=0;
sub_spike_ind=bwlabel_firstind(sp_tr);
sub_spike_ind=sub_spike_ind(1:doublet_length(s));
if sub_spike_ind(end)+15 < length(nTau_db)
    sub_spike_ind_seg=[[[1; sub_spike_ind(1:end-1)] sub_spike_ind]; [sub_spike_ind(end) sub_spike_ind(end)+20]];
else
    sub_spike_ind_seg=[[[1; sub_spike_ind(1:end-1)] sub_spike_ind]; [sub_spike_ind(end) length(nTau_db)]];
end
cmap_sub=turbo(size(sub_spike_ind_seg,1));
lg_label=[];
for s_sub=1:size(sub_spike_ind_seg,1)
    show_seg=[sub_spike_ind_seg(s_sub,1):sub_spike_ind_seg(s_sub,2)];
    soma_seg=squeeze(STA_dbSmat2(1,s,show_seg));
    apical_seg=squeeze(mean(STA_dbSmat2(rois{2},s,show_seg),1));
    plot(soma_seg,apical_seg,'color',cmap_sub(s_sub,:),'Linewidth',1.5); hold all
    lg_label{s_sub}=['Spike#' num2str(s_sub-1)];
end
plot(squeeze(STA_dbSmat2(1,s,sub_spike_ind)),squeeze(mean(STA_dbSmat2(rois{2},s,sub_spike_ind))),'r*')
text(squeeze(STA_dbSmat2(1,s,sub_spike_ind))+1,squeeze(mean(STA_dbSmat2(rois{2},s,sub_spike_ind))),num2str([1:size(sub_spike_ind_seg,1)-1]'))
lg_label(1)={'Pre spike'};
lg_label(end)={'Post spike'};
legend(lg_label,'Location','northwest','FontSize',8);
xlabel('V Soma')
ylabel('V Apical')

ax1=[ax1 nexttile([1 1])];
%imagesc(nTau_db,[1:nROI],squeeze(STA_dbSmat2(dist_order,s,:)),cax);
plot(nTau_db,squeeze(STA_dbSmat2(1,s,:)),'color',cmap(1,:)); hold all
plot(nTau_db,squeeze(mean(STA_dbSmat2(rois{2},s,:),1)),'color',cmap(2,:)); hold all
xlabel('Time (ms)')
end
colormap(turbo)
linkaxes([ax_c ax_d],'xy');
linkaxes([ax1],'xy');

%% Measure Attenuation
nTau_Amp=[-1:2]; cmap=distinguishable_colors(6); ax1=[];
[~, dist_order2]=sort(geodist,'ascend');
singleMat=reshape(NomalizedTrace(:,find(Sp_Class==1)'+nTau_Amp),nROI,[],length(nTau_Amp));
SSSpikeAmp=max(singleMat,[],3);
%SpikeAmp=SpikeAmp./SpikeAmp(1,:);
Attenuation_complex=[]; Attenuation_simple=[]; Attenuation_dbS=[];
CSSpikeAmp=[]; dbSSpikeAmp=[];
for s=1:size(SSSpikeAmp,2)
    [fitresult, gof] = fit(geodist, SSSpikeAmp(:,s), ft);
    Attenuation_simple(s,1)=fitresult.p1;
end

for sp_order=1:5
    for s=1:size(CSprop.SpikeAmp,2)
        if size(CSprop.SpikeAmp{s},2)<sp_order
            Attenuation_complex(s,sp_order)=NaN; CSSpikeAmp(:,s,sp_order)=NaN;
        else
            CSSpikeAmp(:,s,sp_order)=CSprop.SpikeAmp{s}(:,sp_order);
            [fitresult, gof] = fit(geodist, squeeze(CSSpikeAmp(:,s,sp_order)), ft);
            Attenuation_complex(s,sp_order)=fitresult.p1;
        end
    end
end

for sp_order=1:3
    for s=1:size(dbSprop.SpikeAmp,2)
        if size(dbSprop.SpikeAmp{s},2)<sp_order
            Attenuation_dbS(s,sp_order)=NaN; dbSSpikeAmp(:,s,sp_order)=NaN;
        else
            dbSSpikeAmp(:,s,sp_order)=dbSprop.SpikeAmp{s}(:,sp_order);
            [fitresult, gof] = fit(geodist, squeeze(dbSSpikeAmp(:,s,sp_order)), ft);
            Attenuation_dbS(s,sp_order)=fitresult.p1;
        end
    end
end

figure(4); clf; ax1=[];
ax1=[ax1 nexttile([1 1])];
h=histogram(SSSpikeAmp(1,:),[0:0.5:40],'Normalization','probability'); hold all
histogram(mean(SSSpikeAmp(rois{2},:),1),h.BinEdges,'Normalization','probability');
title(['Simple Spike'])
xlabel('Intensity (A.U.)'); legend({'Soma','Apical'});
ax1=[ax1 nexttile([1 1])];
h=histogram(dbSSpikeAmp(1,:,1),[0:0.5:40],'Normalization','probability'); hold all
histogram(mean(dbSSpikeAmp(rois{2},:,1),1),h.BinEdges,'Normalization','probability');
title(['Burst spike #1'])
ax1=[ax1 nexttile([1 1])];
h=histogram(dbSSpikeAmp(1,:,2),[0:0.5:40],'Normalization','probability'); hold all
histogram(mean(dbSSpikeAmp(rois{2},:,2),1),h.BinEdges,'Normalization','probability');
title(['Burst spike #2'])
ax1=[ax1 nexttile([1 1])];
h=histogram(dbSSpikeAmp(1,:,3),[0:0.5:40],'Normalization','probability'); hold all
histogram(mean(dbSSpikeAmp(rois{2},:,3),1),h.BinEdges,'Normalization','probability');
title(['Burst spike #3'])
ax1=[ax1 nexttile([1 1])];
h=histogram(CSSpikeAmp(1,:,1),[0:0.5:40],'Normalization','probability'); hold all
histogram(mean(CSSpikeAmp(rois{2},:,1),1),h.BinEdges,'Normalization','probability');
title(['Complex spike #1'])
ax1=[ax1 nexttile([1 1])];
h=histogram(CSSpikeAmp(1,:,2),[0:0.5:40],'Normalization','probability'); hold all
histogram(mean(CSSpikeAmp(rois{2},:,2),1),h.BinEdges,'Normalization','probability');
title(['Complex spike #2'])
ax1=[ax1 nexttile([1 1])];
h=histogram(CSSpikeAmp(1,:,3),[0:0.5:40],'Normalization','probability'); hold all
histogram(mean(CSSpikeAmp(rois{2},:,3),1),h.BinEdges,'Normalization','probability');
title(['Complex spike #3'])
ax1=[ax1 nexttile([1 1])];
h=histogram(CSSpikeAmp(1,:,4),[0:0.5:40],'Normalization','probability'); hold all
histogram(mean(CSSpikeAmp(rois{2},:,4),1),h.BinEdges,'Normalization','probability');
title(['Complex spike #4'])

figure(5); clf;
M_data={mean(SSSpikeAmp(rois{2},:),1),mean(dbSSpikeAmp(rois{2},:,1),1),mean(dbSSpikeAmp(rois{2},:,2),1),mean(dbSSpikeAmp(rois{2},:,3),1)...
    ,mean(CSSpikeAmp(rois{2},:,1),1),mean(CSSpikeAmp(rois{2},:,2),1),mean(CSSpikeAmp(rois{2},:,3),1),mean(CSSpikeAmp(rois{2},:,4),1)};
M_dataS={mean(SSSpikeAmp(1,:),1),mean(dbSSpikeAmp(1,:,1),1),mean(dbSSpikeAmp(1,:,2),1),mean(dbSSpikeAmp(1,:,3),1),...
    mean(CSSpikeAmp(1,:,1),1),mean(CSSpikeAmp(1,:,2),1),mean(CSSpikeAmp(1,:,3),1),mean(CSSpikeAmp(1,:,4),1)};
M=cell2mat(cellfun(@(x) mean(x,2,'omitnan'),M_data,'UniformOutput',false));
S=cell2mat(cellfun(@(x) std(x,0,2,'omitnan'),M_data,'UniformOutput',false));
Ms=cell2mat(cellfun(@(x) mean(x,2,'omitnan'),M_dataS,'UniformOutput',false));
Ss=cell2mat(cellfun(@(x) std(x,0,2,'omitnan'),M_dataS,'UniformOutput',false));
errorbar([1:8],M,S,'LineStyle','none','Marker','+','color',cmap(2,:)); hold all;
errorbar([1:8],Ms,Ss,'LineStyle','none','Marker','+','color',cmap(1,:)); xlim([0.5 8.5]); ylim([4 28]); 
set(gca,'xtick',[1:8],'xticklabel',{'SS','Brst #1','Brst #2','Brst #3','CS#1','CS#2','CS#3','CS#4'})
legend({'Apical','Soma'})
ylabel('')
