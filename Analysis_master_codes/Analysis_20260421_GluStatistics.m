%% Load file path
clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Prism_V2+Glu_Data_Arrangement.xlsx'], 'Sheet1', 'C5:AA31');
load('/Volumes/BHL18TB_D2/20260203_SD_V2+iGluSNFR4/25X_transformationMatrix.mat');
fpath=raw(:,1)';
V2moviemaxTime=15000;
GlumoviemaxTime=5000;
Endframe=cell2mat(raw(:,5));
foi=[3 5 6 8 9 10 12 13];
StructureData=raw(:,6);
pixelsize=0.468; %µm
%% Show stitched FOVs
SortStruct= unique_cell(StructureData);
FOV=4; maxmatchingdist=5;
StructureStack=read_tiff(StructureData{find(SortStruct.groupIdx==FOV,1)});
Savefilepath=fileparts(StructureData{find(SortStruct.groupIdx==FOV,1)});
swcfiles=dir(fullfile(StructureData{find(SortStruct.groupIdx==FOV,1)}, 'Tracing*.swc'));
isSynDistExist=isfile(fullfile(Savefilepath, 'SynDist_Neuron.mat'));
if isSynDistExist
    load(fullfile(fileparts(StructureData{find(SortStruct.groupIdx==FOV,1)}), 'SynDist_Neuron.mat'));
    disp('SynDist matrix loaded');
else
    SynDist_Neuron=[];
end
Z_step_all=[];
GluResultAll=[]; VoltResultAll=[]; Device_Data_all=[]; g=1;
for f=find(SortStruct.groupIdx==FOV)'
    GluResultAll{g}=importdata(fullfile(fpath{f},"Glu_Result.mat"));
    VoltResultAll{g}=importdata(fullfile(fpath{f},"Volt_Result.mat"));
    load(fullfile(fpath{f},"output_data.mat"));
    Device_Data_all{g}=Device_Data;
    Z_step_all=[Z_step_all VoltResultAll{g}.Z_step];
    g=g+1;
end
try
    StructureImg=max(StructureStack(:,:,[min(Z_step_all)-2:max(Z_step_all)+2]),[],3);
catch
    StructureImg=max(StructureStack(:,:,[min(Z_step_all):max(Z_step_all)]),[],3);
end

swcfiles=dir(fullfile(fileparts(StructureData{f}), 'Tracing*.swc'));
dendritePoints=[];
for si=1:length(swcfiles)
    swcname=fullfile(swcfiles(si).folder, swcfiles(si).name);
    tree = load_tree(swcname);
    dendritePoints{si}=[tree.X, tree.Y];
end
dendritePointsAll=cell2mat(dendritePoints');
DendriteMask=point2img(cell2mat(dendritePoints'),3,size(StructureImg));

VoltFtprnt_bdry=[]; GluFtprnt=[]; ProxyGlu=[]; distance2dendrite=[];
for fi=1:length(GluResultAll)
    GluCoord=get_coord(GluResultAll{fi}.S_glu); GluCoord_str=[];
    [xGluWorld, yGluWorld] = intrinsicToWorld(VoltResultAll{fi}.Rglu, GluCoord(:,1), GluCoord(:,2));
    [xStrWorld, yStrWorld] = transformPointsInverse(VoltResultAll{fi}.tform_Str2Glu, xGluWorld, yGluWorld);
    Rstr = imref2d(size(StructureImg));
    [GluCoord_str(:,1), GluCoord_str(:,2)] = worldToIntrinsic(Rstr, xStrWorld, yStrWorld);
    AlignedFtprnt=transformCamera_O2B(Device_Data_all{fi},tformReg,VoltResultAll{fi}.ftprnt,GluResultAll{fi}.AvgGluImg);
    AlignedFtprnt = imwarp(AlignedFtprnt, VoltResultAll{fi}.Rglu, invert(VoltResultAll{fi}.tform_Str2Glu), 'OutputView', Rstr);
    GluFtprnt{fi}=GluCoord_str;
    distance2dendrite{fi} = matchSglu2SWC(dendritePoints, GluCoord_str);
    ProxyGlu{fi}=distance2dendrite{fi}.neuronID>0 & distance2dendrite{fi}.minDist<maxmatchingdist;
    for roi=1:size(AlignedFtprnt,3)
        bdry=bwboundaries(AlignedFtprnt(:,:,roi)>0);
        bdry_sz=cellfun(@(x) size(x,1),bdry);
        [~, m]=max(bdry_sz);
        VoltFtprnt_bdry{fi}{roi}=(bdry{m});
    end
end

figure(FOV); clf;
cmap_fi=hsv(length(GluResultAll));
cmap_voltfi=lines(length(GluResultAll));
cmap_N=autumn(length(dendritePoints));
imshow2(StructureImg,[]); hold all;
for N=1:length(dendritePoints)
    scatter(dendritePoints{N}(:,1),dendritePoints{N}(:,2),0.1,cmap_N(N,:)*0.6);
end
for fi=1:length(GluResultAll)
    scatter(GluFtprnt{fi}(ProxyGlu{fi},1),GluFtprnt{fi}(ProxyGlu{fi},2),5,cmap_fi(fi,:),'filled');
    for r=1:length(VoltFtprnt_bdry{fi})
        plot(VoltFtprnt_bdry{fi}{r}(:,2),VoltFtprnt_bdry{fi}{r}(:,1),'color',cmap_voltfi(fi,:));
    end
end
%% Calculate glutamate event rate and clusterness
GluSpike=cellfun(@(x) find_spike_bh(zscore(x.dFF_glu-movmedian(x.dFF_glu,25,2),0,2),3,1,5),GluResultAll,'UniformOutput',false);
GluTrace_hi=cellfun(@(x) x.dFF_glu-movmedian(x.dFF_glu,25,2),GluResultAll,'UniformOutput',false);
figure(FOV*3); clf;
show_traces_spikes(GluTrace_hi{1},GluSpike{1},sum(GluSpike{1}));
GluRate=cellfun(@(x,y) sum(x,2)./range(y.t_ax),GluSpike,GluResultAll,'UniformOutput',false);
GluEvents_cluster=[]; GluEvents_stat=[];
ClusterThreshold=2; %event size;
ClusterDistThreshold=10; %µm
GluEvents_stat.GluRate=[];
GluEvents_stat.clusterness=[];
GluEvents_stat.coord=[];
GluEvents_stat.NeuronID=[];
GluEvents_stat.FOVID=[];

for N=1:length(swcfiles)
    SkelDend=point2img(dendritePoints{N},3,size(StructureImg));
    N
    for fi=1:length(distance2dendrite)
        Nsub=find(distance2dendrite{fi}.neuronID==N & distance2dendrite{fi}.minDist<maxmatchingdist);
        if ~isempty(Nsub)
            nSyn=sum(length(Nsub));
            GluCoordSub=GluFtprnt{fi}(Nsub,:);
            if ~isSynDistExist
                SynDist_Neuron{N,fi}=[];
                for i=1:nSyn
                    for j=i:nSyn
                        [SynDist_Neuron{N,fi}(i,j), p]=geodesic_distance(SkelDend,GluCoordSub(i,:),GluCoordSub(j,:));
                    end
                    if mod(i, 20) == 0
                        fprintf('  Frame %d / %d  (%.1f%%)\n', i, nSyn, 100*i/nSyn);
                    end
                end
            end
            GluEvents_cluster{N,fi}= find_spatiotemporal_events(GluSpike{fi}(Nsub,:), SynDist_Neuron{N,fi}, GluCoordSub,'MaxDist',ClusterDistThreshold/pixelsize,'Tolerance',1);
            roi_event_idx=[]; roi_event_size=[];
            for n = 1:nSyn
                roi_event_idx{n}  = find(arrayfun(@(e) ismember(n, e.rois),  GluEvents_cluster{N,fi}));
                roi_event_size{n} = [ GluEvents_cluster{N,fi}(roi_event_idx{n}).size];
            end
            Clusterness=cell2mat(cellfun(@(x) sum(x>=ClusterThreshold)./length(x),roi_event_size,'UniformOutput',false));
            GluEvents_stat.GluRate=[GluEvents_stat.GluRate; GluRate{fi}(Nsub)];
            GluEvents_stat.clusterness=[GluEvents_stat.clusterness; Clusterness'];
            GluEvents_stat.coord=[GluEvents_stat.coord; GluCoordSub];
            GluEvents_stat.NeuronID=[GluEvents_stat.NeuronID; ones(nSyn,1)*N];
            GluEvents_stat.FOVID=[GluEvents_stat.FOVID; ones(nSyn,1)*fi];
        end
    end
end
if ~isSynDistExist
    save(fullfile(Savefilepath, 'SynDist_Neuron.mat'),'SynDist_Neuron');
    disp('Inter-synapse distance matrix is saved.')
end

%% Plot rate and clusterness
figure(2*FOV); clf; tiledlayout(2,1,'Padding','tight');
cmap_N=gen_colormap([[0.3 0.3 0.3]; [0.7 0.7 0.7]],length(dendritePoints));
rate_rng=[0 0.3]; clst_rng=[0 0.2];
RateCmap=vec2cmap(GluEvents_stat.GluRate,turbo(256),rate_rng);
clusternessCmap=vec2cmap(GluEvents_stat.clusterness,turbo(256),clst_rng);

ax1=nexttile([1 1]);
imshow2(StructureImg,[]); hold all;
for N=1:length(dendritePoints)
    scatter(dendritePoints{N}(:,1),dendritePoints{N}(:,2),0.1,cmap_N(N,:)*0.6);
end
scatter(GluEvents_stat.coord(:,1),GluEvents_stat.coord(:,2),10,RateCmap,"filled");
title('Glutamate event rate (Hz)'); 
cb = colorbar; cb.Colormap = turbo;
cb.Ticks = linspace(cb.Limits(1), cb.Limits(2), 2);   % or however many ticks you want
cb.TickLabels = {num2str(rate_rng(1)), num2str(rate_rng(end))};
cb.Label.String="Event rate (Hz)";

nexttile([1 1]);
imshow2(StructureImg,[]); hold all;
for N=1:length(dendritePoints)
    scatter(dendritePoints{N}(:,1),dendritePoints{N}(:,2),0.1,cmap_N(N,:)*0.6);
end
scatter(GluEvents_stat.coord(:,1),GluEvents_stat.coord(:,2),10,clusternessCmap,"filled");
title(sprintf('Coactive probability with synapses within %2.0d µm',ClusterDistThreshold));
cb = colorbar; cb.Colormap = turbo;
cb.Ticks = linspace(cb.Limits(1), cb.Limits(2), 2);   % or however many ticks you want
cb.TickLabels = {num2str(clst_rng(1)), num2str(clst_rng(end))};
cb.Label.String="Cluster probability";

saveas(gcf,fullfile(fileparts(StructureData{find(SortStruct.groupIdx==FOV,1)}), 'Rate_cluster_map.fig'))

%% Plateau vs glutamate analysis
FOV_i=1; VoltROI2look=3;
nTau=[-1000:1000];  SkelDend=point2img(dendritePointsAll,3,size(StructureImg));
Nsub=distance2dendrite{FOV_i}.neuronID==1 & distance2dendrite{FOV_i}.minDist<maxmatchingdist;
distance2soma=[];
for si=1:size(GluResultAll{FOV_i}.S_glu,3)
    if Nsub(si)
[distance2soma(si), p]=geodesic_distance(SkelDend,dendritePointsAll(1,:),GluFtprnt{FOV_i}(si,:));
    else
        distance2soma(si)=inf;
    end
end
% [coord1D principal_axis]= projectTrunkaxis(dendritePointsAll);
% GluCoord_1D=GluFtprnt{FOV_i}*principal_axis;
% [~, Synsort]=sort(GluCoord_1D,'ascend');
[distsorted Synsort]=sort(distance2soma,'ascend');
Nsub=Nsub(Synsort);
GluCoord_sub=GluFtprnt{FOV_i}(Synsort,:);

t_volt=VoltResultAll{FOV_i}.t_ax(1:size(VoltResultAll{FOV_i}.normTraces,2));
t_glu=GluResultAll{FOV_i}.t_ax(1:size(GluResultAll{FOV_i}.dFF_glu,2));
GluTr=GluTrace_hi{FOV_i}(Synsort,:);
GluTr2look= interp1(t_glu, GluTr', t_volt, 'linear', 'extrap')';  % M x T1

GluSp2look= interp1(t_glu, GluSpike{FOV_i}', t_volt, 'linear', 'extrap')';  % M x T1
GluSp2look=double(GluSp2look(Synsort,:)>0.9);

VoltageTr=VoltResultAll{FOV_i}.normTraces;
VoltageTr=zscore(VoltageTr,0,2);
VoltageTr2look=VoltageTr(VoltROI2look,:);
%VoltageTr2look=VoltageTr2look-movprc(VoltageTr2look,300,30,2);
putative_CSframe=[37870 384877 40659 49165 49304 49396 55004 55343 55466 75502 76087 77300 77988 80881 85060 87763 89005 91912 93100 93222 95573 103364 105892 108874 112564 113733 ...
                  114528 115073 115506 116052 117522 118517 118646 119634 124907 151327 155595 161380 166895 209768 226145 244723];

[STAglu STAgluMat taken_CSframe]=get_STA(GluTr2look,putative_CSframe,-nTau(1),nTau(end));
[STAVolt STAVoltMat]=get_STA(VoltageTr,putative_CSframe,-nTau(1),nTau(end));
[STAgluSp STAgluSPMat]=get_STA(GluSp2look,putative_CSframe,-nTau(1),nTau(end));

%% Show STA traces
dendbin=[0:100:750];
[~, ~, ~, distbin_ind]=binning_data({[distsorted' distsorted']*pixelsize},dendbin);

figure(101); clf;
tiledlayout(7,2);
ax1=nexttile([3 1]);
imagesc(nTau,[1:sum(Nsub)],STAglu(Nsub,:),[-0.2 0.5]/10); 
title('CS-triggered average glutamate (∆F/F)');
ylabel('Synapse ID (sorted by distance from soma)'); 

ax1=[ax1 nexttile([3 1])];
imagesc(nTau,[1:sum(Nsub)],squeeze(mean(movsum(STAgluSPMat(Nsub,:,:),36,3),2))); 
title('CS-triggered average glutamate (spike)');

ax2=nexttile([2 1]);
bin_lg=[]; li=[]; g=1;
for b=1:length(distbin_ind)
    if sum(double(distbin_ind{b}))>5
        meanSTAglu=mean(STAglu(distbin_ind{b} & Nsub,:),1,'omitnan');
    li(g)=plot(nTau,meanSTAglu); hold all;
    bin_lg{g}=sprintf('From %2.0d µm to %2.0d µm',dendbin(b),dendbin(b+1));
    g=g+1;
    end
end
arrayfun(@(l,c) set(l,'Color',c{:}),li',num2cell(hsv(length(li)),2)); box off;
ylim([-0.005 0.02]); xlabel('Time (ms)'); ylabel('∆F/F');
legend(bin_lg);

ax4=nexttile([2 1]);
bin_lg=[]; li=[]; g=1;
STAgluSpTol=squeeze(mean(movsum(STAgluSPMat,36,3),2));
for b=1:length(distbin_ind)
    if sum(double(distbin_ind{b}))>6
        meanSTAglu=mean(STAgluSpTol(distbin_ind{b} & Nsub,:),1,'omitnan');
        li(g)=plot(nTau,meanSTAglu); hold all;
        bin_lg{g}=sprintf('From %2.0d µm to %2.0d µm',dendbin(b),dendbin(b+1));
        g=g+1;
    end
end
arrayfun(@(l,c) set(l,'Color',c{:}),li',num2cell(hsv(length(li)),2)); box off;
ylim([0 0.1]); xlabel('Time (ms)'); ylabel('Mean glutamate activity (Hz)');

ax3=nexttile([2 2]);
plot(nTau,STAVolt([3 6],:)');
box off; ylabel('Voltage (Z score)'); xlabel('Time (ms)');
legend('Trunk','Distal');
linkaxes([ax1 ax2 ax3 ax4],'x');
xlim([-400 300]);

%% Show each CS case
figure(100); clf;
tiledlayout(6,7,'padding','tight');
for frm=1:length(taken_CSframe)
    nexttile([1 1])
    %arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(Nsub)),2)); hold all; 
    plot(nTau,squeeze(STAVoltMat(VoltROI2look,frm,:)),'k-','linewidth',1.5); box off;
    ylim([-5 12]); hold all; yyaxis right;
    [rr cc]=find(squeeze(STAgluSPMat(Nsub,frm,:)));
    plot(cc+nTau(1)-1,rr,'r|','linewidth',1.5); 
    xlim([-1000 1000]); ylim([0 sum(Nsub)+0.5])
    title(num2str(frm))
end
%% Glutamate triggered average
nTauGlu=[-300:300];
GluInterest=[26 65 66 70]; VoltageInterest=[3 6];
cmap_GluInterest=lines(length(GluInterest));
cmap_VfprntInterest=hsv(length(VoltageInterest));

GluTr=GluResultAll{FOV_i}.dFF_glu(Synsort,:);
GluTr_sorted=GluTr(Nsub,:);
GluCoord_sub_sorted=GluCoord_sub(Nsub,:);

GluSpike_interest=GluSpike{FOV_i}(Synsort,:);
GluSpike_interest=GluSpike_interest(Nsub,:);
GluSpike_interest=GluSpike_interest(GluInterest,:);

STAgluMat_sorted=STAgluMat(Nsub,:,:);
STAgluSPMat_sorted=STAgluSPMat(Nsub,:,:);

GluSTAvoltage=[];
for glu2show=1:length(GluInterest)
gsp=find(GluSpike_interest(glu2show,:)>0);
[glusp2voltax, ~] = match_nearest(t_glu(gsp), t_volt);
GluSTAvoltage(:,:,glu2show)=get_STA(VoltageTr,glusp2voltax,-nTauGlu(1),nTauGlu(end));
end

figure(103); clf; ax1=[];
tiledlayout(6,length(GluInterest),'padding','tight')
nexttile([2 length(GluInterest)]);
imshow2(StructureImg,[800 13000]); hold all;
scatter(GluCoord_sub_sorted(GluInterest,1),GluCoord_sub_sorted(GluInterest,2),30,cmap_GluInterest,'filled');
g=1;
for vbd=VoltageInterest   
plot(VoltFtprnt_bdry{FOV_i}{vbd}(:,2),VoltFtprnt_bdry{FOV_i}{vbd}(:,1),'color',cmap_VfprntInterest(g,:),'linewidth',1);
g=g+1;
end
xlim([400 1549]); ylim([0 500]);
drawScaleBar(100/pixelsize,'horizontal','color','w')

ax1=[ax1 nexttile([1 length(GluInterest)])];
plot(t_glu,GluTr_sorted(GluInterest,:)'+[1:length(GluInterest)]*1); box off;
title('Glutamate (∆F/F)');

ax1=[ax1 nexttile([1 length(GluInterest)])];
l=plot(t_volt,VoltageTr(VoltageInterest,:)'+[1:length(VoltageInterest)]*15); box off;
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap_VfprntInterest,2));
title('Voltage (Z-score)'); xlabel('Time (s)');
legend('Trunk','Distal');
linkaxes(ax1,'x');

for gin=1:length(GluInterest)
nexttile([1 1]);
l2=plot(nTauGlu,GluSTAvoltage(VoltageInterest,:,gin)');
arrayfun(@(l,c) set(l,'Color',c{:}),l2,num2cell(cmap_VfprntInterest,2));
title(sprintf('Glu-triggered average Syn. #%1.f',gin));
xlabel('Time (ms)'); ylabel('Voltage (Z score)'); 
ylim([-0.3 1]); box off;
xlim([-300 100]);
end


for gin=1:length(GluInterest)
nexttile([1 1]);
plot(nTau,squeeze(STAgluMat_sorted(GluInterest(gin),:,:))','color',[0.5 0.5 0.5]); hold all;
plot(nTau,squeeze(mean(STAgluMat_sorted(GluInterest(gin),:,:),2)),'k')
title(sprintf('CS-triggered Syn. #%1.f',gin));
xlabel('Time (ms)'); ylabel('Glutamate (∆F/F)'); 
%ylim([-0.3 1]); box off;
xlim([-600 200]);
end


%% Generate shuffled STAs
time_bin=[-500:50:300]; N_shuffle=1000;
STAglu_perm=zeros([sum(Nsub),length(time_bin)-1, N_shuffle, 3]);
for perm=1:N_shuffle
    randpoint=randi(size(GluTr2look,2),length(taken_CSframe),1);
    [~, tmp]=get_STA(GluTr2look(Nsub,:),randpoint,-nTau(1),nTau(end));
    [~, tmpSp]=get_STA(GluSp2look(Nsub,:),randpoint,-nTau(1),nTau(end));
    tmpmean=[]; tmpmax=[]; tmpSpsum=[];
    for t=1:length(time_bin)-1
        tmpmax(:,:,t)=max(tmp(:,:,[time_bin(t)+1:time_bin(t+1)]-nTau(1)),[],3);
        tmpmean(:,:,t)=mean(tmp(:,:,[time_bin(t)+1:time_bin(t+1)]-nTau(1)),3);
        tmpSpsum(:,:,t)=sum(tmpSp(:,:,[time_bin(t)+1:time_bin(t+1)]-nTau(1)),3);
    end
    STAglu_perm(:,:,perm,1)=mean(tmpmax,2);
    STAglu_perm(:,:,perm,2)=mean(tmpmean,2);
    STAglu_perm(:,:,perm,3)=mean(tmpSpsum,2);
    if mod(perm, 100) == 0
        fprintf('  Permutation %d / %d  (%.1f%%)\n', perm, N_shuffle, 100*perm/N_shuffle);
    end
end

STA_glu_Real=zeros([sum(Nsub),length(time_bin)-1, 3]);
tmpmean=[]; tmpmax=[]; tmpSpsum=[];
for t=1:length(time_bin)-1
    tmpmax(:,:,t)=max(STAgluMat_sorted(:,:,[time_bin(t)+1:time_bin(t+1)]-nTau(1)),[],3);
    tmpmean(:,:,t)=mean(STAgluMat_sorted(:,:,[time_bin(t)+1:time_bin(t+1)]-nTau(1)),3);
    tmpSpsum(:,:,t)=sum(STAgluSPMat_sorted(:,:,[time_bin(t)+1:time_bin(t+1)]-nTau(1)),3);    
end
STA_glu_Real=cat(4,squeeze(mean(tmpmax,2)),squeeze(mean(tmpmean,2)),squeeze(mean(tmpSpsum,2)));

count = sum(STAglu_perm <= STA_glu_Real, 3);
STA_glu_pval = (1 + count) ./ (1 + N_shuffle);

%% Make p-value movie
pvalSTAglu=[];
output_filename=fullfile(Savefilepath,'GluPvalueMovie.mp4');
v = VideoWriter(output_filename, 'MPEG-4');
v.FrameRate = 10;
v.Quality   = 95;
open(v);
fig = figure('Color','k', 'units','pixels', 'Position',[200 0 1000 800]);
dffclim=[0.05 0.2];
statname={'Amplitude','Mean','# of glu. spike'};

% ===================== RENDER LOOP ======================================
fprintf('Rendering frames...\n');
tiledlayout('padding','tight');
for t=1:length(time_bin)-1
    clf(fig);
    set(fig, 'Color','k');
    tl = tiledlayout(fig, 2, 2, 'TileSpacing','tight', 'Padding','tight');

    sgtitle(sprintf('From %3.0f ms to %3.0f ms',time_bin(t)+1, time_bin(t+1)));

    for stat=1:3
        nexttile([1 1]);
    imshow2(StructureImg,[]); hold all;
    pvalcmap=vec2cmap(STA_glu_pval(:,t,:,stat),'hot',[0.95 1]);
    scatter(GluCoord_sub(Nsub,1),GluCoord_sub(Nsub,2),5,pvalcmap,'filled');
    
    xlim([400 1549]); ylim([0 500]);
    cb = colorbar; cb.Colormap = hot;
    cb.Ticks = linspace(cb.Limits(1), cb.Limits(2), 2);   % or however many ticks you want
    cb.TickLabels = {num2str(0.95), num2str(1)};
    cb.Label.String="Percentile";
    title(statname{stat});
    end

    nexttile([1 1]);
    imshow2(StructureImg,[]); hold all;
    dffcmap=vec2cmap(STA_glu_Real(:,t,:,1),'hot',dffclim);
    scatter(GluCoord_sub(Nsub,1),GluCoord_sub(Nsub,2),5,dffcmap,'filled');
    
    xlim([400 1549]); ylim([0 500]);
    cb = colorbar; cb.Colormap = hot;
    cb.Ticks = linspace(cb.Limits(1), cb.Limits(2), 2);   % or however many ticks you want
    cb.TickLabels = {num2str(dffclim(1)), num2str(dffclim(2))};
    cb.Label.String="∆F/F amplitude";
    title('∆F/F amplitude')

    drawnow;
    writeVideo(v, getframe(fig));
end

close(v);
close(fig);
fprintf('Done. Movie saved to:\n  %s\n', output_filename);




