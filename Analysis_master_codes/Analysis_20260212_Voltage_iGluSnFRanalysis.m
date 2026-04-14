%% Load file path
clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Prism_V2+Glu_Data_Arrangement.xlsx'], 'Sheet1', 'C5:AA31');
load('/Volumes/BHL18TB_D2/20260203_SD_V2+iGluSNFR4/25X_transformationMatrix.mat');
fpath=raw(:,1)';
V2moviemaxTime=15000;
GlumoviemaxTime=5000;
Endframe=cell2mat(raw(:,5));
foi=[3 5 6 8 9 10 11 12 13];
StructureData=raw(:,6);
%% Map Glu Footprint onto morphology
% Voltage image >> Structure image >> Glu

f=3;
GluResult=importdata(fullfile(fpath{f},"Glu_Result.mat"));
VoltResult=importdata(fullfile(fpath{f},"Volt_Result.mat"));
load(fullfile(fpath{f},"output_data.mat"));
swcfiles=dir(fullfile(fpath{f}, 'Tracing*.swc'));
if isempty(swcfiles)
swcfiles=dir(fullfile(fileparts(StructureData{f}), 'Tracing*.swc'));
end
sz=size(GluResult.AvgGluImg);

cam1_vsyn = Device_Data{1,2}.Counter_Inputs(1,1).data;
start_idx = find(cam1_vsyn == max(cam1_vsyn), 1);
end_idx   = numel(Device_Data{1,2}.buffered_tasks(1,2).channels(1,1).data);
segment_size = 103 * floor(Device_Data{3}.exposuretime/0.001);
n_to_add = end_idx - start_idx + 1;
vals = repelem((cam1_vsyn(start_idx-1)+1 : cam1_vsyn(start_idx-1)+ceil(n_to_add/segment_size))', segment_size);
cam1_vsyn(start_idx:end_idx) = vals(1:n_to_add);

cam2_vsyn = Device_Data{1, 2}.Counter_Inputs(1, 2).data;
cam2_vsyn_trig = find (diff(cam2_vsyn) == 1)+1;
segment_size2 = cam2_vsyn_trig (10) - cam2_vsyn_trig(9);
start_idx = min(find(cam2_vsyn ==max (cam2_vsyn)));
end_idx = length(Device_Data{1, 2}.buffered_tasks(1, 2).channels(1, 2).data);
last_val = cam2_vsyn(start_idx - 1);
n_to_add = end_idx - start_idx + 1;
n_segments = ceil(n_to_add / segment_size2);
added_part = repelem((last_val + 1 : last_val + n_segments)', segment_size2);
added_part = added_part(1:n_to_add);
cam2_vsyn(start_idx:end_idx) = added_part;

EncoderTick=double(diff(Device_Data{2}.buffered_tasks(3).channels(2).data(1:segment_size*max(cam1_vsyn)+1)>2))>0;
Speed=sum(reshape(EncoderTick,segment_size,[]),1);
Speed=Speed(cam1_vsyn(200001)+2:end);

EncoderTick2=double(diff(Device_Data{2}.buffered_tasks(3).channels(2).data(1:segment_size2*max(cam2_vsyn)+1)>2))>0;
Speed2=sum(reshape(EncoderTick2,segment_size2,[]),1);
Speed2=Speed2(cam2_vsyn(200001)+3:end);

TracingCoord=[]; TracingCoord_structure=[];
for si=1:length(swcfiles)
swcname=fullfile(swcfiles(si).folder, swcfiles(si).name);
tree = load_tree(swcname);
if ~isnan(StructureData{f})
    [xGlu, yGlu] = transformPointsForward(VoltResult.tform_Str2Glu, tree.X, tree.Y);
    [xPlot, yPlot] = worldToIntrinsic(VoltResult.Rglu, xGlu, yGlu);
    TracingCoord{si}=[xPlot, yPlot];
    TracingCoord_structure{si}=TracingCoord{si};
else
    TracingCoord{si}=[tree.X, tree.Y];
    TracingCoord_structure{si}=TracingCoord{si};
end
end
fprintf('%2.0f swc files are loaded\n', length(swcfiles))

    alignedVolt = transformCamera_O2B(Device_Data, tformReg, VoltResult.ref_im, GluResult.AvgGluImg);
    alignedVoltftprnt =[];
    for ft=1:size(VoltResult.ftprnt,3)
        alignedVoltftprnt(:,:,ft) = transformCamera_O2B(Device_Data, tformReg, VoltResult.ftprnt(:,:,ft), GluResult.AvgGluImg);
    end
    TracingCoord_Glu=TracingCoord_structure;
GluCoord=get_coord(GluResult.S_glu);
GluMatched2Neuron= matchSglu2SWC(TracingCoord_Glu, GluCoord, 'MaxDist', 5);
DendCmap=hsv(length(swcfiles)); GluColored=zeros(sz(1),sz(2),3);
for si=1:length(swcfiles)
GluColored=GluColored+repmat(max(GluResult.S_glu(:,:,GluMatched2Neuron.neuronID==si),[],3),1,1,3).*reshape(DendCmap(si,:),1,1,[]);
end
GluSpike=find_spike_bh(zscore(GluResult.dFF_glu,0,2),3,1);

%% Plot Glu and voltage traces
clf; ax1=[];
tiledlayout(4,1);
ax1=[ax1 nexttile([2 1])];
GluSubset=GluMatched2Neuron.neuronID==2;
plot(GluResult.t_ax,rescale2(GluResult.dFF_glu(GluSubset,:),2)'+[1:sum(GluSubset)])
ax1=[ax1 nexttile([1 1])];
plot(GluResult.t_ax,sum(GluSpike(GluSubset,:),1));
ax1=[ax1 nexttile([1 1])];
Vtr=VoltResult.traces;
Vtr=Vtr-movmedian(Vtr,500,2);
plot(VoltResult.tax(1:end-1),rescale2(Vtr,2)'+[1:size(Vtr,1)]*0.5)
linkaxes(ax1,'x')

figure(11); clf;
nexttile([1 1])
imshow2(GluColored+grs2rgb(alignedVolt,gray(256)),[]);

nexttile([1 1])
for si=1:length(swcfiles)
scatter(TracingCoord_Glu{si}(:,1),TracingCoord_Glu{si}(:,2),6,DendCmap(si,:),'filled'); hold all
axis equal tight off
set(gca,'YDir','reverse')
end

nexttile([1 1])
SumTglu=[];
for si=1:length(swcfiles)
SumTglu(si,:)=mean(GluResult.T_glu(GluMatched2Neuron.neuronID==si,:),1);    
end
l=plot(GluResult.t_ax,SumTglu);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(DendCmap,2))

show_traces_spikes(GluResult.dFF_glu,GluSpike,SumTglu)
%%
load(fullfile(fpath{f},"output_data.mat"));
alignedVolt_ftprnt = transformCamera_O2B(Device_Data, tformReg, VoltResult.ftprnt, GluResult.AvgGluImg);
vTax=VoltResult.tax(1:size(VoltResult.normTraces,2));
interactive2chViewer(GluResult.dFF_glu, GluResult.S_glu, GluResult.t_ax, VoltResult.normTraces, alignedVolt_ftprnt, vTax)

%% Generate movie to show Glu and voltage together



%% Plot representative time points
%time2show=[242 414 2033]; 
%time2show=[1258 2596 2599]; 
%time2show=[31 45 65 111 129 142 242 268 306 352 362 414 418 421 506 563 638 701 769 982 993 1514 2033 2107 2464 3561]; 
time2show=find(GluSpike(402,:)>0);
GluWindow=[-1:1]; DendCmap=lines(length(TracingCoord_Glu));
GluSpikeN_dendrite=zeros(length(TracingCoord_Glu),length(time2show));
figure(100); clf; hold all; tiledlayout('TileSpacing','tight');
set(gcf, 'Color', 'k');
for t=1:length(time2show)
    nexttile([1 1]); hold all
    ActiveGlu=max(GluSpike(:,time2show(t)+GluWindow),[],2)>0;
    GluSpikeN_dendrite(1,t)=sum(ActiveGlu & GluMatched2Neuron.neuronID==0);
    for si=1:length(TracingCoord_Glu)
    scatter(TracingCoord_Glu{si}(:,1),TracingCoord_Glu{si}(:,2),6,DendCmap(si,:),'filled','MarkerFaceAlpha',0.6); 
    GluSpikeN_dendrite(si+1,t)=sum(ActiveGlu & GluMatched2Neuron.neuronID==si);
    end
    scatter(GluCoord(:,1),GluCoord(:,2),15,[0.6 0.6 0.6],'filled'); 
    scatter(GluCoord(ActiveGlu,1),GluCoord(ActiveGlu,2),20,[1 0 0],'filled'); 
    title(['Frame #', num2str(time2show(t))]);
    axis equal tight off
end

figure(101); clf;
h=bar([1:length(time2show)],GluSpikeN_dendrite,'stacked');
arrayfun(@(l,c) set(l,'FaceColor',c{:}),h',num2cell([[0.2 0.2 0.2]; DendCmap],2))

%% Calculate the Joint Probability and population coactivity
Tuplesize=2; tau=1; % tolerance for synchronization
JointProb_pair = get_JointProbability(GluSpike, Tuplesize, tau);
Jmat=zeros(size(GluSpike,1),size(GluSpike,1));
for i=1:size(JointProb_pair.combos,1)
Jmat(JointProb_pair.combos(i,1),JointProb_pair.combos(i,2))=JointProb_pair.P_joint(i)-JointProb_pair.P_ind(i);
end
kernel = ones(1, 2*tau + 1);
GluSpike_D = conv2(double(GluSpike), kernel, 'same') > 0;  % [M x T] logical
CoactiveN=sum(GluSpike_D,1)/size(GluSpike_D,1);

figure(22); clf;
tiledlayout(2,2,'TileSpacing','tight');
nexttile([1 2]);
plot(GluResult.t_ax,CoactiveN); box off;
xlabel('Time (s)'); ylabel('Coactive probability');

nexttile([1 1]);
histogram(sum(GluSpike_D,1),35);
set(gca,'yscale','log');
xlabel('# of coactive synapses'); ylabel('Probability'); box off;

nexttile([1 1]);
SameMat=GluMatched2Neuron.neuronID == GluMatched2Neuron.neuronID'; SameMat_vec=[];
for i=1:size(JointProb_pair.combos,1)
    SameMat_vec(i)=SameMat(JointProb_pair.combos(i,1),JointProb_pair.combos(i,2));
end
scatter(JointProb_pair.P_ind(SameMat_vec==0),JointProb_pair.P_joint(SameMat_vec==0),6,'k','filled'); hold all
scatter(JointProb_pair.P_ind(SameMat_vec>0),JointProb_pair.P_joint(SameMat_vec>0),6,'r','filled');
mi=max([JointProb_pair.P_ind; JointProb_pair.P_joint]);
plot([0:0.0001:mi],[0:0.0001:mi],'r');
set(gca,'XScale','log','YScale','log'); 
xlabel('Independent probability (P_i * P_j)');
ylabel('Joint probability (P_{ij})')
%% Find clusters
K = 10; % number of clusters
tau=1; % tolerance for synchronization
kernel = ones(1, 2*tau + 1);
GluSpike_D = conv2(double(GluSpike), kernel, 'same') > 0;  % [M x T] logical
GluSpike_demean = GluSpike_D- mean(GluSpike_D,1);
C = corr(double(GluSpike_demean)');
Z = linkage(squareform(1-C), 'average');
T = cluster(Z, 'maxclust', K);
[H,~,perm] = dendrogram(Z, 0); % perm is leaf order
N = size(C,1);
leafPos = zeros(N,1); leafPos(perm) = 1:N;
orderTbl = [(1:N)' T leafPos];
orderTbl = sortrows(orderTbl, [2 3]);   % by cluster, then leaf position
ord = orderTbl(:,1);
Csort = C(ord, ord);  Tsort = T(ord);
cols = hsv(K);

figure(23); clf;
tiledlayout(1,3,'TileSpacing','tight');
axM=nexttile([1 1]);
imagesc(Csort); caxis([-0.05 0.1]);
axis square tight; colorbar;
title('Correlation matrix (sorted by cluster)');
xlabel('Synapses (sorted)'); ylabel('Synapses (sorted)');

nexttile([1 2]);
Cluster2show=5; Bplot = GluSpike(ord,:);
for k = 1:K
    rows = find(Tsort==k); [r,c] = find(Bplot(rows,:));
    plot(GluResult.t_ax(c), rows(r), '|', 'Color', cols(k,:), 'MarkerSize',4,'linewidth',2); hold all;
end
xlabel('Time (s)'); ylabel('Synapse ID'); axis tight off;

pos = axM.Position;
axTop = axes('Position',[pos(1), pos(2)+pos(4), pos(3), 0.03]);
imagesc(Tsort');               % 1 x N
colormap(axTop, cols);
set(axTop,'XTick',[],'YTick',[]);
axis tight; axTop.XAxisLocation = 'top';

figure(24); clf; hold all;
%imshow2(alignedVolt,[]); hold all;
for si=1:length(TracingCoord_Glu)
    scatter(TracingCoord_Glu{si}(:,1),TracingCoord_Glu{si}(:,2),3,DendCmap(si,:)/2,'filled')
end
markers = {'o','s','^','d','v','>','<','d','o','^'};  % define enough markers
for k = [1:10]
    scatter(GluCoord(T==k,1),GluCoord(T==k,2),30,cols(k,:),markers{k},'filled');
end
axis equal tight off;

%% Cross-correlation of synapses on different dendrite
tau=1; % tolerance for synchronization
kernel = ones(1, 2*tau + 1);
GluSpike_D = conv2(double(GluSpike), kernel, 'same') > 0;  % [M x T] logical
%CorrMat=corr(GluSpike_D');
CorrMat=Jmat;
dendIdx= GluMatched2Neuron.neuronID;
clustIdx= T;
Cluster2show=[7 8 9 10];

N = size(CorrMat,1); 
catA = [];  % same dendrite, same cluster
catB = [];  % same dendrite, diff cluster
catC = [];  % diff dendrite, same cluster
catD = [];  % diff dendrite, diff cluster

for i = 1:N-1
    for j = i+1:N
        sameD = dendIdx(i)==dendIdx(j);
        sameC = clustIdx(i)==clustIdx(j);
        isclusterinterested=all(ismember([clustIdx(j) clustIdx(i)],Cluster2show));

        if sameD && sameC && isclusterinterested
            catA(end+1) = CorrMat(i,j);
        elseif sameD && ~sameC && isclusterinterested
            catB(end+1) = CorrMat(i,j);
        elseif ~sameD && sameC && isclusterinterested
            catC(end+1) = CorrMat(i,j);
        elseif isclusterinterested
            catD(end+1) = CorrMat(i,j);
        end        
    end
end

figure(26); clf; tiledlayout(1,4);
nexttile([1 1]); cmap=bone(2);
Boxplot_wPoints2({[catA catB],[catC catD]},cmap);
pmat=get_pValue({[catA catB],[catC catD]},0,'Method','ttest2');
ylim([-0.005 0.015]); drawPValueLines(pmat,0.001,'StepHeight',0.001,'TextYOffset', 0.0005);
ylim([-0.005 0.02])
set(gca,'XTick',[1:2],'XTickLabel',{'Same D','Different D'}); ylabel('P_{ij} - P_i*P_j'); box off;
nexttile([1 1]); cmap=bone(2);
Boxplot_wPoints2({[catA catC],[catB catD]},cmap);
pmat2=get_pValue({[catA catC],[catB catD]},0,'Method','ttest2');
ylim([-0.005 0.015]); drawPValueLines(pmat2,0.001,'StepHeight',0.001,'TextYOffset', 0.0005);
ylim([-0.005 0.02])
set(gca,'XTick',[1:2],'XTickLabel',{'Same C','Different C'}); ylabel('P_{ij} - P_i*P_j'); box off;
nexttile([1 2]); cmap=bone(4);
Boxplot_wPoints2({catA,catB,catC,catD},cmap);
pmat3=get_pValue({catA,catB,catC,catD},0,'Method','ttest2');
ylim([-0.005 0.015]); drawPValueLines(pmat3,0.001,'StepHeight',0.001,'TextYOffset', 0.0005);
ylim([-0.005 0.02])
set(gca,'XTick',[1:4],'XTickLabel',{'Same D & Same C','Same D & Different C','Different D & Same C','Different D & Different C'})
ylabel('P_{ij} - P_i*P_j'); box off;
set_fontsize(10);

%% Glutamate triggered average voltage
GluTA=[];
VoltTr=VoltResult.traces-movmedian(VoltResult.traces,1000,2);
VoltTr=VoltTr./get_threshold(VoltTr,1);
for sy=1:size(GluSpike,1)
GluOntime=GluResult.t_ax(find(GluSpike(sy,:)>0));
D = abs(bsxfun(@minus, VoltResult.tax(:), GluOntime(:).'));
[dt, idx] = min(D, [], 1);
idx=idx(dt<0.1);
closestVals = VoltResult.tax(idx);
GluTA(:,:,sy)=get_STA(VoltTr,idx,100,100);
end

%% Load movie
nCol2 = double(Device_Data{4}.ROI([2]));  % ROI on the camera
nRow2 = double(Device_Data{4}.ROI([4]));  % ROI on the camera
nTime = size(GluResult.t_ax,2);
try
mov2_mc=double(readBinMov_times([fpath{f} '/mc2' num2str(1,'%02d') '.bin'], nRow2, nCol2,[1:nTime]));
catch
mov2_mc=double(readBinMov([fpath{f} '/mc2' num2str(1,'%02d') '.bin'], nRow2, nCol2));
end
load([fpath{f} '/mc2Trace' num2str(1,'%02d') '.mat'])
meanF=squeeze(mean(mov2_mc,[1 2]));
y_fit=expfitDM_2([1:size(mov2_mc,3)]',meanF,[1:size(mov2_mc,3)]',[10000 100]);
mov2_mc_filt=mov2_mc./reshape(y_fit,1,1,[]);
mov2_mc_filt_lw=movmedian(mov2_mc_filt,round(5/Device_Data{4}.exposuretime),3);
mov_sub=mov2_mc_filt-mov2_mc_filt_lw;

[V, D, u_sub] = get_eigvector(tovec(imgaussfilt(mov_sub,1))',40);
u_sub=reshape(u_sub,nRow2(1),nCol2(1),[]);

figure(2); clf; imshow2_patch(u_sub); drawnow;
n=input("PCs to regress out\n");
mov_sub=SeeResiduals(mov_sub,V(:,n));
mov_sub=SeeResiduals(mov_sub,mcTrace2.xymean([1:nTime],:));
mov2_mc_filt2=mov_sub+mov2_mc_filt_lw;
mov_sub=mov_sub.*(max(GluResult.bvMask,[],3)==0);
mov_sub_DoG=DoG_filter(mov_sub,1.6,27);
se = strel('sphere', 3); % Define a structuring element
GluMask = imdilate(max(GluResult.S_glu(:,:,GluMatched2Neuron.neuronID>0),[],3)>0,[se]);
mov_sub_DoG=mov_sub_DoG.*GluMask;
StdImg=std(mov_sub_DoG,0,3);
mov_sub_zscore=mov_sub_DoG./StdImg;
mov_sub_zscore(isnan(mov_sub_zscore) | isinf(mov_sub_zscore))=median(mov_sub_zscore(:),'omitnan');
%% Gen Movie
fig=figure(100);
set(fig, 'Color', 'w');
maskdendtrace=point2img(cell2mat(TracingCoord_Glu'),2,sz);
imshowrange=[0.1 0.25];
myVideo = VideoWriter(fullfile(fpath{f},'GluMovie_wTracing'),"MPEG-4"); %open video file 
myVideo.FrameRate = 20;  %can adjust this, 5 - 10 works well for me
myVideo.Quality= 100;
pixelsize=6.5/25*180/100;
open(myVideo)
cax=[0 1];
for i=[1:nTime]
    clf; tiledlayout(2,1,'TileSpacing','compact');
    nexttile([1 1]);
    pbaspect([size(double(mov_sub_DoG),2) size(double(mov_sub_DoG),1) 1]*2)
    mov_sub_DoG_colored=grs2rgb(mov_sub_DoG(:,:,i),gen_colormap([0 0 0;1 0 0],256),imshowrange(1),imshowrange(2));
    imshow2(mov_sub_DoG_colored*3+repmat(mat2gray(alignedVolt),1,1,3)*2.*reshape([0 0.6 1],1,1,[]),cax); hold all

    nexttile([1 1]);
    DendCmap=lines(length(TracingCoord_Glu));
    imshow2(mov_sub_DoG_colored*3,cax); hold all
    for si=1:length(TracingCoord_Glu)
    scatter(TracingCoord_Glu{si}(:,1),TracingCoord_Glu{si}(:,2),6,DendCmap(si,:),'filled','MarkerFaceAlpha',0.1);
    end
    axis tight off equal
    title([num2str(GluResult.t_ax(i),3) ' s'])
    drawScaleBar(100/pixelsize,'horizontal','color','w');
    drawnow;
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end
close(myVideo)

%%
f=5;
Devicedata_filename=fullfile(fpath{f},'output_data.mat');
load(Devicedata_filename);
% set orange cam
sz1 = double(Device_Data{3}.ROI([2 4]));  % ROI on the camera
exposuretime1 = Device_Data{3}.exposuretime;

cam1_vsyn = Device_Data{1, 2}.Counter_Inputs(1, 1).data;
start_idx = min(find(cam1_vsyn ==max (cam1_vsyn)));
end_idx = length(Device_Data{1, 2}.buffered_tasks(1, 2).channels(1, 1).data);
segment_size = 103*floor(exposuretime1/0.001);
last_val = cam1_vsyn(start_idx - 1);
n_to_add = end_idx - start_idx + 1;
n_segments = ceil(n_to_add / segment_size);
added_part = repelem((last_val + 1 : last_val + n_segments)', segment_size);
added_part = added_part(1:n_to_add);
cam1_vsyn(start_idx:end_idx) = added_part;
VmovTimesegments=[(cam1_vsyn(o_laser)+2):V2moviemaxTime:cam1_vsyn(end)];
nFrame2analyze=VmovTimesegments(end)-VmovTimesegments(1)+1;
t_vol = (0:nFrame2analyze-1)*exposuretime1;

cam2_vsyn = Device_Data{1, 2}.Counter_Inputs(1, 2).data;
cam2_vsyn_trig = find (diff(cam2_vsyn) == 1)+1;
segment_size = cam2_vsyn_trig (10) - cam2_vsyn_trig(9);
start_idx = min(find(cam2_vsyn ==max (cam2_vsyn)));
end_idx = length(Device_Data{1, 2}.buffered_tasks(1, 2).channels(1, 2).data);
last_val = cam2_vsyn(start_idx - 1);
n_to_add = end_idx - start_idx + 1;
n_segments = ceil(n_to_add / segment_size);
added_part = repelem((last_val + 1 : last_val + n_segments)', segment_size);
added_part = added_part(1:n_to_add);
cam2_vsyn(start_idx:end_idx) = added_part;

mod488 = 200001;
mod607 = 200001;
cam2_trig = find (diff(cam2_vsyn) ==1)+1;
cam2_trig = cam1_vsyn (cam2_trig);
cam2_trig = cam2_trig (cam2_vsyn(mod488)+2:cam2_vsyn(end))'- (cam1_vsyn(mod607)+2);
if cam2_vsyn(end)<GlumoviemaxTime
    GluemovTimesegments=[(cam2_vsyn(mod488)+2) cam2_vsyn(end)];
else
    GluemovTimesegments=[(cam2_vsyn(mod488)+2):GlumoviemaxTime:cam2_vsyn(end)];
end

Blue_Nframe = Device_Data{4}.frames_requested;
nCol2 = double(Device_Data{4}.ROI([2]));  % ROI on the camera
nRow2 = double(Device_Data{4}.ROI([4]));  % ROI on the camera
exposuretime2 = Device_Data{4}.exposuretime;
t_cal = (cam2_trig-1)*exposuretime2;