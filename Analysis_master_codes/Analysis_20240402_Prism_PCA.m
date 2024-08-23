clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/FromBackup/PP72_PlaceCellResults';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:P23');

% [~, ~, NeuronsToUse]=xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
%     'PlaceCellData_Arrangement.xlsx'], 'Sheet1', 'L8:M46');
%
% NeuronsToUse=cellfun(@(x) (str2num(num2str(x))),NeuronsToUse,'UniformOutput',false);
ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,11),'UniformOutput',false);
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
f=11; load(fullfile(fpath{f},'PC_Result.mat'),'Result')
nROI=size(Result.normTraces,1);
nTau={[-20:20],[-60:200],[-20:20]}; %SS, CS, dSP
spclass_ind=1;
%load aligned movie somatic spike
alignmovlist=dir(fullfile(fpath{f},[alignedMovFN{spclass_ind} '*.tiff']));
AlignMov=[];
for l=1:length(alignmovlist)
    l
    AlignMov=cat(3,AlignMov,readtiff(fullfile(fpath{f},alignmovlist(l).name)));
end
sz_align=size(AlignMov);
AlignMov=double(reshape(AlignMov,sz_align(1),sz_align(2),length(nTau{spclass_ind}),[]));
%AlignMov=AlignMov-median(AlignMov(:,:,:,:),3);
AlignMov=AlignMov-mean(AlignMov(:,:,1:6,:),3);
%AlignMov=AlignMov-mean(maxk(AlignMov,7,3),3);
%AlignMov=AlignMov-prctile(AlignMov,80,3);

noi=setdiff([1:nROI],[7 15 8 9 10]); %f11
%noi=setdiff([1:nROI],[20 18 6 7 8]); %f14
SkelDend = Skeletonize_dendrite(Result.ref_im,4,0.02,25);
interDendDist=[];
for i=1:size(Result.normTraces,1)
    for j=1:size(Result.normTraces,1)
        [interDendDist(i,j), ~]=geodesic_distance(SkelDend,get_coord(Result.ftprnt(:,:,i)),get_coord(Result.ftprnt(:,:,j)));
    end
end
[~, dist_order_noi]=sort(interDendDist(1,noi),'ascend');
[~, dist_order]=sort(interDendDist(1,:),'ascend');
dist_order_inv_noi=mod(find([1:length(noi)]==dist_order_noi'),length(noi));
dist_order_inv_noi(dist_order_inv_noi==0)=length(noi);
dist_order_inv=mod(find([1:nROI]==dist_order'),nROI);
dist_order_inv(dist_order_inv==0)=nROI;

toi_heatmap=[18:19]; toi_pca=[13:19]; toi_movie=[10:25]; 

% Isolated Somatic spike
som_spike=Result.StackedSpike{1}(2,:); bAP_s=[]; 
dSpike=Result.StackedSpike{3}(2,:);
for s=som_spike
    isnearby=sum(ismember(s+nTau{1},som_spike))>1;
    isnearby=isnearby | sum(ismember(s+nTau{1},dSpike))>1;
    if ~isnearby & ~isnan(s)
        bAP_s=[bAP_s s];
    end
end
som_list=ismember(Result.StackedSpike{1}(2,:),bAP_s);

% % dSpike inducing Soma spike
% dSP_s=find(Result.SpClass(3,:)>0 & sum(Result.spike(2:end,:),1)>0);
% dSP_s_IndSoma=[];
% som_spike=find(Result.spike(1,:)>0);
% for s=dSP_s
%     isnearby=sum(ismember(s+[0 1 2],som_spike))>0;
%     if isnearby
%         [~, shift]=max(Result.normTraces(1,s+[0 1 2]));
%         dSP_s_IndSoma=[dSP_s_IndSoma s+shift-1];
%     end
% end

%PCA on Soma spike trace
nKeep = 20;
somRoi=tovec(Result.ftprnt(bound:end-bound,bound:end-bound,1)>0);
AlignMov_sub=reshape(AlignMov(:,:,:,som_list),sz_align(1)*sz_align(2)*length(nTau{1}),[]);
somF=squeeze(min(sum(reshape(AlignMov_sub.*repmat(somRoi,length(nTau{1}),1),sz_align(1)*sz_align(2),length(nTau{1}),[]),1),[],2));
AlignMov_sub=AlignMov_sub./somF';
AlignMov_sub=reshape(AlignMov_sub,sz_align(1),sz_align(2),length(nTau{1}),[]);

prc_normTr=Result.normTraces-movprc(Result.normTraces,60,35,2);
STA_SSmat=reshape(prc_normTr(:,bAP_s'+nTau{1}),nROI,[],length(nTau{1}));
%STA_SSmat=STA_SSmat_large(:,:,1-nTau_large(1)+nTau2);
%STA_SSmat=STA_SSmat-median(STA_SSmat(:,:,1:12),3);
%STA_SSmat=STA_SSmat-median(mink(STA_SSmat,10,3),3);
STA_SS=squeeze(mean(STA_SSmat,2));
F_ref=mean(STA_SS(:,ceil(size(STA_SSmat,3)/2)+[9:12]),2);
STA_SSmat=STA_SSmat./F_ref;

STA_SSmatVec=tovec(permute(STA_SSmat(noi,:,toi_pca),[1 3 2]));
%STA_SSmatVec=tovec(permute(mean(STA_SSmat(noi,:,toi_pca),3),[1 3 2]));
covMat=STA_SSmatVec*STA_SSmatVec';
[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;

figure(8); clf;
nexttile([1 1])
plot(cumsum(D./sum(D))); hold all
bar([1:length(D)],D./sum(D)); xlim([0.5 nKeep+0.5]); ylim([0 1])

% STA_SSeigImg = (V'*STA_SSmatVec);
% STA_SSeigImg=toimg(STA_SSeigImg',length(noi),length(toi_pca));
% STA_SSeigTr=tovec(STA_SSeigImg)'*STA_SSmatVec';
% 
STA_SSeigTr = (V'*STA_SSmatVec);
STA_SSeigImg=STA_SSeigTr(1:nKeep,:)*STA_SSmatVec';
[STA_SScoV, STA_SScoV_P]=corr(STA_SSeigImg',STA_SSmatVec);
STA_SSeigImg=toimg(STA_SSeigImg',length(noi),length(toi_pca));
figure(12); clf;
for i=1:nKeep
    nexttile([1 1])
    imagesc(STA_SSeigImg(dist_order_noi,:,i))
    colormap('turbo')
    title(num2str(i))
end

figure(13); clf;
pos_bin=150; pos_range=[113 120]; 
for i=[4]
PCs=[1 i];
nexttile([1 1])
cmap_pos=turbo(pos_range(2)-pos_range(1)+1);
cmap_pos=[repmat(cmap_pos(1,:),pos_range(1)-1,1); cmap_pos; repmat(cmap_pos(end,:),pos_bin-pos_range(2),1)];
Sp_PosBin=ceil(Result.VR(5,bAP_s)/115*pos_bin);
sp_PosRange = Sp_PosBin>pos_range(1) & Sp_PosBin<pos_range(2);
scatter(STA_SSeigTr(PCs(1),sp_PosRange),STA_SSeigTr(PCs(2),sp_PosRange),[],cmap_pos(Sp_PosBin(sp_PosRange),:),'filled'); hold all
%scatter(STA_SSeigTr(PCs(1),~sp_PosRange),STA_SSeigTr(PCs(2),~sp_PosRange),[],'k','filled');
xlabel(['PC #' num2str(PCs(1))])
ylabel(['PC #' num2str(PCs(2))])
end

Ls=[];
for i=[4]%:19

    gridN=4; PCs=[1 i];
    PCmovie=zeros(sz_align(1),sz_align(2),length(nTau{1}),gridN^2);
    PCSSmat=zeros(size(STA_SSmat,1),size(STA_SSmat,3),gridN^2);
    P1_grid=[min(STA_SSeigTr(PCs(1),:)): range(STA_SSeigTr(PCs(1),:))/gridN :max(STA_SSeigTr(PCs(1),:))];
    P2_grid=[min(STA_SSeigTr(PCs(2),:)): range(STA_SSeigTr(PCs(2),:))/gridN :max(STA_SSeigTr(PCs(2),:))];
    PCmovietile=zeros(sz_align(1)*gridN,sz_align(2)*gridN,length(nTau{1}));
    PC_nspike=zeros(gridN^2,1);
    for g1=1:gridN
        for g2=1:gridN
            g_list= P1_grid(g1)<STA_SSeigTr(PCs(1),:) & P1_grid(g1+1)>STA_SSeigTr(PCs(1),:) & P2_grid(g2)<STA_SSeigTr(PCs(2),:) & P2_grid(g2+1)>STA_SSeigTr(PCs(2),:);
            PC_nspike(g1+(g2-1)*gridN)=sum(g_list);
            if ~isempty(find(g_list))
                PCSSmat(:,:,g1+(g2-1)*gridN)=mean(STA_SSmat(:,g_list,:),2);
                PCmovie(:,:,:,g1+(g2-1)*gridN)=mean(AlignMov_sub(:,:,:,g_list),4);
            else
                PCSSmat(:,:,g1+(g2-1)*gridN)=zeros(size(STA_SSmat,1),size(STA_SSmat,3));
                PCmovie(:,:,:,g1+(g2-1)*gridN)=zeros(sz_align(1),sz_align(2),length(nTau{1}));
            end
            PCmovietile(sz_align(1)*(g2-1)+1:sz_align(1)*g2,sz_align(2)*(g1-1)+1:sz_align(2)*g1,:)=PCmovie(:,:,:,g1+(g2-1)*gridN);

            s=zeros(1,size(Result.normTraces,2));
            s(1,bAP_s(find(g_list)))=1;
            Ls(:,:,g1+(g2-1)*gridN,i)=PlaceTrigger_average(s,300,Result.VR,0,115,'sum');
        end
    end

    figure(15); clf;
    tiledlayout(gridN,gridN*2)
    nexttile([gridN/2 gridN])
    imshow2([mean(PCmovietile(:,:,[16:19]),3)],[-0.2 0.5]*0.001)
    nexttile(gridN^2+1,[gridN/2 gridN])
    imshow2([imgaussfilt(mean(PCmovietile(:,:,[16:19]),3),1)],[-0.2 0.5]*0.001)
    for g=1:gridN^2
        nexttile([1 1])
        imagesc(PCSSmat(dist_order,:,g),[-7 4])
        xlim([12 22])
        title([num2str(g) ', N = ', num2str(PC_nspike(g))])
    end
    colormap(turbo)
%    saveas(gca,[save_at 'F14PCImg#1#' num2str(i) '.fig'])
 %   writeMov([save_at 'F14PCMovie#1#' num2str(i)],imgaussfilt3(PCmovietile,[1 1 0.1]),[],5,1,[-0.2 0.6]*0.001)
end
figure; moviefixsc(imgaussfilt3(PCmovietile,[1 1 0.1]),[-0.2 0.6]*0.001)

figure(17); clf;
LS_mov=ringmovMean(Ls,3);
ax3=[];
for i=2
    for g=1:16
ax3=[ax3 nexttile([1 1])];
imagesc(Ls(:,:,g,i))
    end
end
linkaxes(ax3,'xy')
colormap('turbo')

figure(17); clf;
LsM=[];
rois={[2 3 4 7 8],[10 11],[12 13 14 15]};
for i=1:length(rois)
    LsM(i,:)=squeeze(mean(sum(Ls(:,:,rois{i},2),3,'omitnan'),1,'omitnan'));
end
LsM=ringmovMean(LsM,3);
plot(LsM')

%%
pos_bin=300;
%pcratio=rescale(STA_SSeigTr(1,:))./rescale(STA_SSeigTr(2,:));
%pcratio(abs(pcratio)>200)=0;
pcs=[1 2];

s_bin=zeros(2,size(Result.normTraces,2));
for i=[1:2]
s_bin(i,bAP_s)=rescale(STA_SSeigTr(pcs(i),:));
Lap_PCratio(:,:,i) = PlaceTrigger_average(s_bin(i,:),pos_bin,Result.VR,0,115,'sum');
end
Lap_Spike = PlaceTrigger_average(double(s_bin(1,:)~=0),pos_bin,Result.VR,0,115,'sum');
figure(16); clf;
Lap_PCratio(Lap_PCratio==0)=NaN;
Lap_Spike(Lap_Spike==0)=NaN;

ax1=nexttile([1 1]);
imagesc(Lap_PCratio(:,:,1)./Lap_Spike)
ax2=nexttile([1 1]);
imagesc(Lap_PCratio(:,:,2)./Lap_Spike)
colormap(turbo)
linkaxes([ax1 ax2],'xy')
%%
pos_bin=300;
SomSp_VRbin=ceil(Result.VR(5,bAP_s)/115*pos_bin);
STAmovie_pc=[]; STA_SSmat_pc=[];
for p=1:pos_bin
    if ~isempty(find(SomSp_VRbin==p))
    STAmovie_pc(:,:,:,p)=mean(AlignMov_sub(:,:,:,find(SomSp_VRbin==p)),4);
    STA_SSmat_pc(:,:,p)=squeeze(mean(STA_SSmat(:,find(SomSp_VRbin==p),:),2));
    else
    STAmovie_pc(:,:,:,p)=zeros(sz_align(1),sz_align(2),length(nTau{1}));
    STA_SSmat_pc(:,:,p)=zeros(size(STA_SSmat,1),size(STA_SSmat,3));
    end
end
figure(17); clf;
for p=[238:248]
    nexttile([1 1])
    imagesc(STA_SSmat_pc(dist_order,:,p),[0.5 1.5])
    colormap(turbo)
end


%% correlation between branches

[power_tr freq_traces]=fft_simple(Result.normTraces,1000);
for n=1:size(Result.normTraces,1)
ACF(n,:)=xcorr(Result.normTraces(n,:));
end
figure(31); clf;        
%plot(freq_traces',power_tr')
plot(rescale2(ACF,2)')




%%
nTau2=[-40:40];
% STA_SSmat_large=reshape(Result.normTraces(:,bAP_s'+nTau_large),nROI,[],length(nTau_large));
% tic;
% STA_SSmat_large=STA_SSmat_large - reshape(movprc(tovec(STA_SSmat_large),60,25,2),nROI,[],length(nTau_large));
% toc;
prc_normTr=Result.normTraces-movprc(Result.normTraces,60,35,2);
STA_SSmat=reshape(prc_normTr(:,bAP_s'+nTau2),nROI,[],length(nTau2));
%STA_SSmat=STA_SSmat_large(:,:,1-nTau_large(1)+nTau2);
%STA_SSmat=STA_SSmat-median(STA_SSmat(:,:,1:12),3);
%STA_SSmat=STA_SSmat-median(mink(STA_SSmat,10,3),3);
STA_SS=squeeze(mean(STA_SSmat,2));
F_ref=mean(STA_SS(:,ceil(size(STA_SSmat,3)/2)+[9:12]),2);
STA_SSmat=STA_SSmat./F_ref;

figure(9); clf;
nexttile([1 1])
l=plot(squeeze(mean(STA_SSmat(dist_order,:,:),2))');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2))
nexttile([1 1])
l=plot(squeeze(std(STA_SSmat(dist_order,:,:),0,2))');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2))

t_sub=[17:19];
SomSp_VRbin=ceil(Result.VR(5,bAP_s)/115*pos_bin);
rois={[12 13 14],[3 4 5 6]}; %basal & apical f11
%rois={[21 17 18 19 15],[10 11 12 4 5 9]}; %basal & apical f14
STA_SS2=squeeze(mean(STA_SSmat,2));

basal_SSmat=squeeze(mean(STA_SSmat(rois{1},:,:),1));
apical_SSmat=squeeze(mean(STA_SSmat(rois{2},:,:),1));
SomSp_VRbin=Result.VR(5,bAP_s)/115*pos_bin;
SomSP_Lap=Result.VR(8,bAP_s);

figure(89); clf; scale=5;
tiledlayout(1,2)
ax1=nexttile([1 1]);
plot(basal_SSmat([1:40:700],:)'+[1:18]*scale)
hold all;
plot([0 size(STA_SSmat,3)],repmat([0; 0],1,18)+[1:18]*scale,'k')
%basal_SSmat=basal_SSmat-max(basal_SSmat,[],2)+(SpikeH(1)-mean(F_ref(rois{1})));
%basal_SSmat=basal_SSmat-mean(basal_SSmat(:,21+[5:8]),2);
ax2=nexttile([1 1]);
plot(basal_SSmat([1:40:700],:)'+[1:18]*scale)
hold all;
plot([0 size(STA_SSmat,3)],repmat([0; 0],1,18)+[1:18]*scale,'k')
linkaxes([ax1 ax2],'xy')

basal_SStrace=NaN(1,size(Result.normTraces,2)); Apical_SStrace=NaN(1,size(Result.normTraces,2));
basal_SStrace(1,tovec(bAP_s+nTau2'))=tovec(basal_SSmat');
Apical_SStrace(1,tovec(bAP_s+nTau2'))=tovec(apical_SSmat');
figure(99); clf;
plot(mean(Result.normTraces(rois{1},:),1)); hold all
plot(mean(Result.normTraces(rois{2},:),1)); 
plot(basal_SStrace+10)
plot(Apical_SStrace+10)
plot([1 599999],[0 0],'color',[0.7 0.7 0.7])
plot([1 599999],[10 10],'color',[0.7 0.7 0.7])
yyaxis right
plot(Result.mcTrace(:,[1 end]),'k')

pos_bin=150; pf_bin=7;  pf=[108 122];
pf_edge=[pf(1):range(pf)/(pf_bin-1):pf(2)];

figure(18); clf;
ax1=[]; cmap=turbo(length(pf_edge-1));

for i=1:length(pf_edge)-1
%ax1=[ax1 nexttile([1 1])];
Sp_PF=find(SomSp_VRbin>=pf_edge(i) & SomSp_VRbin<pf_edge(i+1));
length(Sp_PF)
%l=plot(basal_SSmat(Sp_PF,:)',apical_SSmat(Sp_PF,:)');
plot(mean(basal_SSmat(Sp_PF,:),1)',mean(apical_SSmat(Sp_PF,:),1)','color',cmap(i,:)); hold all
%arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(length(Sp_PF)),2))
xlabel('Basal')
ylabel('Apical')
end
%linkaxes(ax1,'xy')
PFBinMovietile=[];
for i=1:length(pf_edge)-1
Sp_PF=find(SomSp_VRbin>=pf_edge(i) & SomSp_VRbin<pf_edge(i+1));
PFBinMovietile=[PFBinMovietile; mean(AlignMov_sub(:,:,:,Sp_PF),4)];
end
figure; moviefixsc(imgaussfilt3(PFBinMovietile,[2 2 0.1]),[-0.25 0.5]*0.001)
%colormap(turbo)
colormap(gen_colormap([0 0.5 1; 1 1 1; 1 0 0]))


figure(19); clf;
ax2=[];
for l=unique(SomSP_Lap)
    ax2=[ax2 nexttile([1 1])];
Sp_PF=find(SomSp_VRbin>=pf_edge(1) & SomSp_VRbin<pf_edge(end) & SomSP_Lap==l);
%l=plot(basal_SSmat(Sp_PF,t_sub)',apical_SSmat(Sp_PF,t_sub)')
plot(mean(basal_SSmat(Sp_PF,t_sub),2)',mean(apical_SSmat(Sp_PF,t_sub),2)','k'); hold all
scatter(mean(basal_SSmat(Sp_PF,t_sub),2)',mean(apical_SSmat(Sp_PF,t_sub),2)',[],turbo(length(Sp_PF)),'Filled'); 
plot([-10 10],[-10 10],'--','color',[0.7 0.7 0.7])
%arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(length(Sp_PF)),2))
%arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(length(Sp_PF)),2))
xlabel('Basal')
ylabel('Apical')
title(['Lap #' num2str(l)])
end
linkaxes(ax2,'xy')
xlim([-4 4]); ylim([-4 4])

figure(20); clf;
ax2=[]; ax3=[];
cmap=turbo(100);
bin_tick=pf_edge;
for l=unique(SomSP_Lap)
    ax2=[ax2 nexttile([1 1])];
Sp_PF=find(SomSp_VRbin>=pf_edge(1) & SomSp_VRbin<pf_edge(end) & SomSP_Lap==l);
%l=plot(basal_SSmat(Sp_PF,t_sub)',apical_SSmat(Sp_PF,t_sub)');
sub_position=ceil((SomSp_VRbin(Sp_PF)-pf_edge(1))/range(pf_edge)*100);
sub_position(sub_position==0)=1;
PF_frm=find(Result.VR(5,:)/115*pos_bin>pf_edge(1) & Result.VR(5,:)/115*pos_bin<pf_edge(end) & Result.VR(8,:)==l);
[~, t_tick]=min(abs((Result.VR(5,PF_frm)/115*pos_bin)'-bin_tick),[],1);
t_tick=t_tick*0.001;
tick_subpos=ceil((bin_tick-pf_edge(1))/range(pf_edge)*100);
tick_subpos(tick_subpos==0)=1;

scatter(mean(basal_SSmat(Sp_PF,t_sub),2)',mean(apical_SSmat(Sp_PF,t_sub),2)',12,cmap(sub_position,:),'Filled'); hold all
plot([-10 10],[-10 10],'--','color',[0.7 0.7 0.7])
%arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(length(Sp_PF)),2))
%arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(length(Sp_PF)),2))
xlabel('Basal')
ylabel('Apical')
title(['Lap #' num2str(l)])
xlim([-3 3])
ylim([-3 3])
ax3=[ax3 nexttile([1 1])];
%plot(SomSp_VRbin(Sp_PF),mean(basal_SSmat(Sp_PF,t_sub),2)-mean(apical_SSmat(Sp_PF,t_sub),2),'marker','.','markersize',12)
if ~isempty(Sp_PF)
    ft = fittype('poly1'); % 'poly1' specifies a first-degree polynomial (linear fit)
t=(bAP_s(Sp_PF)-bAP_s(Sp_PF(1)))*0.001;
dif=mean(basal_SSmat(Sp_PF,t_sub),2)-mean(apical_SSmat(Sp_PF,t_sub),2);
[fitresult, gof] = fit(t', dif, ft);
    plot((bAP_s(Sp_PF)-bAP_s(Sp_PF(1)))*0.001,mean(basal_SSmat(Sp_PF,t_sub),2)-mean(apical_SSmat(Sp_PF,t_sub),2),'marker','.','markersize',12);
    hold all
    for tt=1:length(t_tick)
    plot([t_tick(tt) t_tick(tt)],[-3 3],'color',cmap(tick_subpos(tt),:))
    end
    %plot(SomSp_VRbin(Sp_PF),mean(basal_SSmat(Sp_PF,t_sub),2)-mean(apical_SSmat(Sp_PF,t_sub),2),'marker','.','markersize',12)
ylim([-6 6])
xlabel('Time (s)')
%xlabel('VR Position (bin)')
ylabel('Basal - Apical');
plot(t,fitresult(t),'r','LineStyle','--')
title([num2str(fitresult.p1,2) 't+' num2str(fitresult.p2,2)])
end
end
linkaxes(ax2,'xy');

figure(21); clf;
Sp_PF=find(SomSp_VRbin>=pf_edge(1) & SomSp_VRbin<120);
plot(mean(basal_SSmat(Sp_PF,t_sub),2)',mean(apical_SSmat(Sp_PF,t_sub),2)','.','color',[0 0.5 1],'markersize',14); hold all
Sp_PF=find(SomSp_VRbin>=120 & SomSp_VRbin<pf_edge(end));
plot(mean(basal_SSmat(Sp_PF,t_sub),2)',mean(apical_SSmat(Sp_PF,t_sub),2)','.','color',[1 0 0],'markersize',14)
plot([-10 10],[-10 10],'--','color',[0.7 0.7 0.7])
xlabel('Basal')
ylabel('Apical')
legend('Before Reward','After Reward')

figure(22); clf;
pos_edgeN=150;
pos_edge=[0:(pos_bin)/pos_edgeN:pos_bin];
cmap2=turbo(length(pos_edge)-1);
ax4=[]; ax5=[];
%tiledlayout(2,pos_edgeN)
g=1;
for p=1:length(pos_edge)-1
Sp_PF=find(SomSp_VRbin>=pos_edge(p) & SomSp_VRbin<pos_edge(p+1));
if length(Sp_PF)>5
ax4=[ax4 nexttile(2*g-1,[1 1])];
M_B=mean(basal_SSmat(Sp_PF,:)); M_A=mean(apical_SSmat(Sp_PF,:));
M_B=M_B-mean(M_B(ceil(length(nTau2)/2)+[7:12])); M_A=M_A-mean(M_A(ceil(length(nTau2)/2)+[7:12]));
plot(M_B,M_A,'color',cmap2(p,:)); hold all
ax5=[ax5 nexttile(2*g,[1 1])];
errorbar_shade(nTau2,M_B,std(basal_SSmat(Sp_PF,:),0,1),[1 0 0]); hold all
errorbar_shade(nTau2,M_A,std(apical_SSmat(Sp_PF,:),0,1),[0 0.5 1]);
g1=mod(p,4)+4*(mod(p,4)==0);
g2=floor((p-1)/4)+1;
title(['N = ', num2str(length(Sp_PF)), ' ,Bin = ', num2str(p)])
g=g+1;
end
end
linkaxes(ax4,'xy'); linkaxes(ax5,'xy');

figure(25); clf;
PosBinMovietile=[]; g1=1;
tile_to_show=[110 116 120 121 123 125];
for p=1:length(tile_to_show)
Sp_PF=find(SomSp_VRbin>=tile_to_show(p)-1 & SomSp_VRbin<tile_to_show(p));
if isempty(Sp_PF)
    PosBinMovietile(sz_align(1)*(g1-1)+1:sz_align(1)*g1,1:sz_align(2),:)=zeros(sz_align(1),sz_align(2),length(nTau{1}));
else
    PosBinMovietile(sz_align(1)*(g1-1)+1:sz_align(1)*g1,1:sz_align(2),:)=mean(AlignMov_sub(:,:,:,Sp_PF),4);
end
g1=g1+1;
end
moviefixsc(imgaussfilt3(PosBinMovietile,[1 1 0.1]),[-0.1 0.4]*0.001)
colormap(gen_colormap([0 0.5 1; 1 1 1; 1 0 0]))

figure(26); clf;
imagesc(imgaussfilt(mean(PFBinMovietile(:,:,17:19),3),2))
axis equal tight off
colormap(gen_colormap([0 0.5 1; 1 1 1; 1 0 0]))
caxis([-1 1]*0.0001)

figure(23); clf;
nexttile([1 1]);
imagesc(Result.Lap_FR)
colormap('turbo')
hold all
for p=1:length(pos_edge)
plot([pos_edge(p) pos_edge(p)],[0.5 size(Result.Lap_FR,1)+0.5],'r'); hold all
end


deltaSubBA=mean(basal_SSmat(:,t_sub),2)-mean(apical_SSmat(:,t_sub),2);
del_trace=zeros(1,size(Result.normTraces,2));
del_trace(bAP_s)=deltaSubBA;
del_trace_b=del_trace; del_trace_b(del_trace<0)=0;
del_trace_a=del_trace; del_trace_a(del_trace>0)=0;
Lap_SumSpike=PlaceTrigger_average(double(del_trace~=0),pos_bin*2,Result.VR,0,115,'sum');
Lap_deltaSubBA=PlaceTrigger_average(del_trace,pos_bin*2,Result.VR,0,115,'sum');
Lap_deltaSubBA_b=PlaceTrigger_average(del_trace_b,pos_bin*2,Result.VR,0,115,'sum');
Lap_SumSpike_b=PlaceTrigger_average(double(del_trace_b~=0),pos_bin*2,Result.VR,0,115,'sum');
Lap_deltaSubBA_a=PlaceTrigger_average(del_trace_a,pos_bin*2,Result.VR,0,115,'sum');
Lap_SumSpike_a=PlaceTrigger_average(double(del_trace_a~=0),pos_bin*2,Result.VR,0,115,'sum');

figure(24); clf;
tiledlayout(3,1); ax7=[];
ax7=[ax7 nexttile([1 1])];
Lap_deltaSubBA(isnan(Lap_deltaSubBA))=0;
imagesc(Lap_deltaSubBA,[-5 5])
colormap(gen_colormap([0 0.5 1; 1 1 1; 1 0 0]))
xlabel('VR position (bin)')
ylabel('Laps')

ax7=[ax7 nexttile([1 1])];
Lap_show=Lap_deltaSubBA./Lap_SumSpike; Lap_show(isnan(Lap_show))=0;
imagesc(Lap_show,[-3 3])
colormap(gen_colormap([0 0.5 1; 1 1 1; 1 0 0]))
xlabel('VR position (bin)')
ylabel('Laps')

ax7=[ax7 nexttile([1 1])];
M_b=mean(Lap_deltaSubBA_b./Lap_SumSpike_b,1,'omitnan'); M_b(isnan(M_b))=0;
M_a=mean(Lap_deltaSubBA_a./Lap_SumSpike_a,1,'omitnan'); M_a(isnan(M_a))=0;
plot([1:pos_bin*2],M_b,'color',[1 0 0]); hold all
plot([1:pos_bin*2],-M_a,'color',[0 0.5 1]);
xlabel('VR position (bin)')
ylabel('Mean pre-spike subthreshold')
linkaxes(ax7,'x')


basal_trace=zeros(1,size(Result.normTraces,2));
basal_trace(bAP_s)=mean(basal_SSmat(:,t_sub),2);
apical_trace=zeros(1,size(Result.normTraces,2));
apical_trace(bAP_s)=mean(apical_SSmat(:,t_sub),2);
Lap_Basal=PlaceTrigger_average(basal_trace,pos_bin*2,Result.VR,0,115,'sum');
Lap_Apical=PlaceTrigger_average(apical_trace,pos_bin*2,Result.VR,0,115,'sum');
figure(27); clf;
tiledlayout(6,1); ax7=[];
ax7=[ax7 nexttile([2 1])];
Lap_Basal(isnan(Lap_Basal))=0;
Lap_Basal_perSpike=Lap_Basal./Lap_SumSpike;
Mbin_SubB=mean(Lap_Basal_perSpike,'omitnan');
Mbin_SubB(isnan(Mbin_SubB))=0;
Lap_Basal_perSpike(isnan(Lap_Basal_perSpike))=0;
imagesc(Lap_Basal_perSpike,[-4 4])
colormap(gen_colormap([0 0.5 1; 1 1 1; 1 0 0]))
xlabel('VR position (bin)')
ylabel('Laps')

ax7=[ax7 nexttile([2 1])];
Lap_Apical(isnan(Lap_Apical))=0;
Lap_Apical_perSpike=Lap_Apical./Lap_SumSpike;
Mbin_SubA=mean(Lap_Apical_perSpike,'omitnan');
Mbin_SubA(isnan(Mbin_SubA))=0;
Lap_Apical_perSpike(isnan(Lap_Apical_perSpike))=0;
imagesc(Lap_Apical_perSpike,[-2 2])
colormap(gen_colormap([0 0.5 1; 1 1 1; 1 0 0]))
xlabel('VR position (bin)')
ylabel('Laps')
nexttile([2 1])
plot(movmean(Mbin_SubA,5,'omitnan'),'r'); hold all
plot(movmean(Mbin_SubB,5,'omitnan'),'color',[0 0.5 1])
colormap(jet)

%%
rois={[12 13 14],[3 4 5 6]}; %basal & apical f11
%rois={[21 17 18 19 15],[10 11 12 4 5 9]}; %basal & apical f14

pos_bin=300;
prc_normTr=Result.normTraces-movprc(Result.normTraces,60,35,2);
STA_SSmat=reshape(prc_normTr(:,bAP_s'+nTau2),nROI,[],length(nTau2));
STA_SS=squeeze(mean(STA_SSmat,2));
F_ref=mean(STA_SS(:,ceil(size(STA_SSmat,3)/2)+[9:12]),2);
STA_SSmat=STA_SSmat./F_ref;
tr_dff=Result.normTraces./F_ref;

t=[1:size(Result.normTraces,2)]/1000;
voltage_b=mean(tr_dff(rois{1},:),1);
voltage_a=mean(tr_dff(rois{2},:),1);
subth_b=get_subthreshold(voltage_b,Result.spike(1,:),5,7);
subth_a=get_subthreshold(voltage_a,Result.spike(1,:),5,7);
subth_b(Result.motionReject)=NaN; subth_a(Result.motionReject)=NaN;

figure;
plot(t,voltage_b,t,subth_b); hold all
plot(t,voltage_a+20,t,subth_a+20);

Lap_SubB=PlaceTrigger_average(subth_b,pos_bin,Result.VR,0,115,'rate');
Lap_SubA=PlaceTrigger_average(subth_a,pos_bin,Result.VR,0,115,'rate');
Lap_SubB(isnan(Lap_SubB))=0; Lap_SubA(isnan(Lap_SubA))=0;
figure(25); clf;
tiledlayout(4,1)
ax2=nexttile([1 1]);
imagesc(ringmovMean(Lap_SubB,3),[-2000 2000]);      
colormap(ax2,gen_colormap([0 0.5 1; 1 1 1; 1 0 0]))
title('Basal')
ax3=nexttile([1 1]);
imagesc(ringmovMean(Lap_SubA,3),[-2000 2000]);      
colormap(ax3,gen_colormap([0 0.5 1; 1 1 1; 1 0 0]))
title('Apical')
ax4=nexttile([1 1]);
imagesc(ringmovMean(Lap_SubB-Lap_SubA,3),[-2000 2000]);      
colormap(ax4,gen_colormap([0 0.5 1; 1 1 1; 1 0 0]))
title('Basal - Apical')
ax5=nexttile([1 1]);
plot(movmean(mean(Lap_SubB,1,'omitnan'),3,'omitnan'),'color',[1 0 0]); hold all
plot(movmean(mean(Lap_SubA,1,'omitnan'),3,'omitnan'),'color',[0 0.5 1])
plot(movmean(mean(Lap_SubB-Lap_SubA,1,'omitnan'),3,'omitnan'),'k')
linkaxes([ax2 ax3 ax4],'xy')
linkaxes([ax2 ax5],'x')
xlabel('VR position (bin)')

CS_spike_Trace=max(Result.spike([1 rois{1}],:),[],1).*Result.CStrace;
Lap_CS=PlaceTrigger_average(CS_spike_Trace,pos_bin,Result.VR,0,115,'rate');
Lap_dSP=PlaceTrigger_average(Result.SpClass(3,:),pos_bin,Result.VR,0,115,'rate');
figure(29); clf;
tiledlayout(3,1)
nexttile([1 1])
imagesc(im_merge(cat(3,ringmovMean(Lap_SubA,3),ringmovMean(Lap_CS,3)),[1 1 1; 1 0 0]))
nexttile([1 1])
imagesc(im_merge(cat(3,ringmovMean(Lap_SubA,3),ringmovMean(Lap_dSP,3)),[1 1 1; 1 0 0]))
nexttile([1 1])
plot(rescale(movmean(mean(Lap_SubB,1,'omitnan'),3,'omitnan'))); hold all
plot(rescale(movmean(mean(Lap_dSP,1,'omitnan'),3,'omitnan')))
plot(rescale(movmean(mean(Lap_CS,1,'omitnan'),3,'omitnan')))
%%
% imagesc(abs(STA_SScoV.*(STA_SScoV_P<0.01)));
% colormap('turbo')

STA_SScoVSig=STA_SScoV.*(STA_SScoV_P<0.01);
CorrEvents=abs(STA_SScoVSig)>0.5;

Corr_movie=[]; Eigmov_corr=[];
AlignMov_sub=reshape(AlignMov(:,:,toi_movie,som_list),sz_align(1)*sz_align(2)*length(toi_movie),[]);
AlignMov_sub=AlignMov_sub./min(AlignMov_sub,[],1);
%AlignMov_sub=AlignMov_sub./prctile(AlignMov_sub,0.01,1);
for i=1:size(CorrEvents,1)
Eigmov_corr(:,:,:,i)=reshape(mean(AlignMov_sub(:,find(CorrEvents(i,:))),2),sz_align(1),sz_align(2),[]);
end
mean_STA=reshape(mean(AlignMov_sub,2),sz_align(1),sz_align(2),[]);
for j=1:5
    for i=1:2   
        Corr_movie(sz_align(1)*(i-1)+1:sz_align(1)*i,sz_align(2)*(j-1)+1:sz_align(2)*j,:)=Eigmov_corr(:,:,:,j+(i-1)*5);%-mean_STA;
    end
end
figure(7); clf;
nexttile([1 1])
%moviefixsc([Corr_movie; imgaussfilt3(Corr_movie,[2 2 0.1])],[-0.01 0.02])
moviefixsc([Corr_movie; imgaussfilt3(Corr_movie,[2 2 0.1])],[-0.01 0.05])

% for i=1:nKeep
%     %plot(rescale(STA_SSeigTr(i,:)-movmedian(STA_SSeigTr(i,:),30))+i-0.5)
%     hold all
% end


figure(9); clf;
for i=1:12
    nexttile([1 1])
    %imagesc(STA_SSeigImg(dist_order,:,i))
    STA_eigCorrImg(:,:,i)=squeeze(mean(STA_SSmat(:,find(CorrEvents(i,:)),:),2));
    imagesc(STA_eigCorrImg(dist_order,[7:25],i),[0.5 1.5]);
    title(['Pattern#' num2str(i) ', N=' num2str(sum(CorrEvents(i,:)))])
    %STA_eigImgSub=mean(STA_SSeigImg(:,toi_heatmap-toi_pca(1)+1,i),2,'omitnan');
    ftmap=double(Result.ftprnt(:,:,dist_order)>0.07).*reshape(mean(STA_eigCorrImg(dist_order,toi_heatmap,i),2),1,1,[]);
    ftprnt_XY=get_coord(Result.ftprnt(:,:,dist_order));
    ftmap(ftmap==0)=NaN;
    STA_SSeigImg_ftprnt(:,:,i)=max(ftmap,[],3);
    nexttile([1 1])
    imagesc(STA_SSeigImg_ftprnt(:,:,i),[0.5 1.5])
    text(ftprnt_XY(:,1),ftprnt_XY(:,2),num2str([1:nROI]'),'color','w')
    colormap('turbo')
end

pcoi={[1],[2],[3]};
pcoi_vennmat=[];
for p=1:length(pcoi)
    pcoi_vennmat(:,p)=max(double(CorrEvents(pcoi{p},:)),[],1)+1;
end
V=venn_data(pcoi_vennmat);
figure(10); clf;
[~, S]=venn(V,'FaceColor',{'r','y','b'},'FaceAlpha',{0.5,0.6,0.7},'EdgeColor','black');
text(S.ZoneCentroid(:,1),S.ZoneCentroid(:,2),num2str(V'),'HorizontalAlignment','center');
text(S.Position(:,1),S.Position(:,2)+S.Radius'/2,cellfun(@num2str,pcoi,'UniformOutput',false),'color','r','HorizontalAlignment','center')
axis equal off

sp_corr_trace=[];
for i=1:10
sp_corr_trace(i,:)=zeros(1,size(Result.normTraces,2));
sp_corr_trace(i,bAP_s(find(CorrEvents(i,:))))=1;
Lap_sp_Corr(:,:,i)=PlaceTrigger_average(sp_corr_trace(i,:),150,Result.VR,0,115);  
end

figure(11); clf;
for i=1:8
    nexttile([1 1])
imagesc(ringmovMean(Lap_sp_Corr(:,:,i),5))
colormap('turbo')
colorbar
title(num2str(i))
end


figure(12); clf;
rois={[13 8 9],[16 19 17]}; g=1;
cmap=distinguishable_colors(4);
cmap_light=distinguishable_colors(4)+0.3; cmap_light(cmap_light>1)=1;
for i=[1 2 3 7]
CorrContour=permute(STA_SSmat(dist_order,find(CorrEvents(i,:)),:),[1 3 2]);
MContour=mean(CorrContour,3);
CorrContour=reshape(CorrContour,nROI,[]);
%plot(mean(CorrContour(rois{1},:),1),mean(CorrContour(rois{2},:),1),'color',cmap_light(g,:)); hold all
plot(mean(MContour(rois{1},:),1),mean(MContour(rois{2},:),1),'marker','.','markersize',8,'color',cmap(g,:)); hold all
xlabel(['ROI ' num2str(rois{1})])
ylabel(['ROI ' num2str(rois{2})])
g=g+1;
end

%%
%rois={[12 13 14],[3 4 5 6]}; %basal & apical f11
%rois={[1 2],[5 6]}; f1
%rois={[21 17 18 19 15],[10 11 12 4 5 9]}; %basal & apical f14
%rois={[1 10],[6 7 8]}; %basal & apical f4
rois={[9 10],[4 5 6]}; %basal & apical f19
som_spike=Result.SpClass(1,:); bAP_s=[]; 
for s=find(som_spike)
    isnearby=sum(ismember(s+nTau{1},som_spike))>1;
    %isnearby=isnearby | sum(ismember(s+nTau{1},dSpike))>1;
    if ~isnearby & ~isnan(s)
        bAP_s=[bAP_s s];
    end
end

CS_spike=Result.SpClass(2,:); CS_s=find(CS_spike); 

prc_normTr=Result.normTraces;%-movprc(Result.normTraces,120,35,2);
STA_SSmat=reshape(prc_normTr(:,bAP_s'+nTau{1}),nROI,[],length(nTau{1}));
STA_SSmat=STA_SSmat-mean(STA_SSmat(:,:,1:6),3);
STA_SS=squeeze(mean(STA_SSmat,2));
F_ref=mean(STA_SS(:,ceil(size(STA_SSmat,3)/2)+[9:12]),2);
STA_SSmat=STA_SSmat./F_ref;

STA_CSmat=reshape(prc_normTr(:,CS_s'+nTau{2}),nROI,[],length(nTau{2}));
STA_CSmat=STA_CSmat-median(STA_CSmat(:,:,1:30),3);
STA_CSmat=STA_CSmat./F_ref;
STA_CS=squeeze(mean(STA_CSmat,2));


cmap=distinguishable_colors(length(rois));
figure(33); clf;
for i=1:length(rois)
    nexttile([1 1])
M_ss=squeeze(mean(mean(STA_SSmat(rois{i},:,:),1),2));
S_ss=squeeze(std(mean(STA_SSmat(rois{i},:,:),1),0,2));

M_cs=squeeze(mean(mean(STA_CSmat(rois{i},:,:),1),2));
S_cs=squeeze(std(mean(STA_CSmat(rois{i},:,:),1),0,2));

errorbar_shade(nTau{1},M_ss,S_ss,cmap(1,:))
hold all
errorbar_shade(nTau{2},M_cs,S_cs,cmap(2,:))
xlim([-40 150])
end

figure(34); clf;
for i=1:length(rois)
    nexttile([1 1])
l=plot(nTau{2}',squeeze(mean(STA_CSmat(rois{i},:,:),1))','color',cmap(2,:));
hold all
plot(nTau{1},mean(STA_SS(rois{i},:),1),'color',[0 0 0])
%arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(40),2))
xlim([-40 150])
end

