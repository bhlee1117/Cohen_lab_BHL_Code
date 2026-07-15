%% Section 0 : Setup
clear
clc;

addpath('D:\Labmember\Data\ByungHun\VU_analysiscode');
addpath(genpath('C:\Users\Lab Member\Documents\GitHub\Cohen_lab_BHL_Code'));

% Remove any shadow of the built-in 'graph' (old cvx/sdpt3 defines graph.m,
% which breaks the centroid clustering in extractVoltronST / extractGluSNFR3).
gshadow = which('graph','-all');
for gi = 1:numel(gshadow)
    if ~contains(gshadow{gi}, fullfile('toolbox','matlab'))   % keep MATLAB's built-in
        rmpath(fileparts(gshadow{gi}));
        fprintf('Removed shadowing graph.m folder from path: %s\n', fileparts(gshadow{gi}));
    end
end
iswindow=0;
try
    [~, ~, raw] = xlsread('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Prism_V2+Glu_Data_Arrangement.xlsx', 'Sheet1', 'C5:AA31');
catch
    raw = readcell(macToWindowsPath('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Prism_V2+Glu_Data_Arrangement.xlsx'),...
        'Sheet', 'Sheet1', 'Range', 'C5:AA31');
    iswindow=1;
end

fpath           = raw(:,1)';
fpath_valid=cell2mat(cellfun(@(x) any(~ismissing(x)),fpath,'UniformOutput',false));
if iswindow
    for f=find(fpath_valid>0)
        fpath{f}=macToWindowsPath(fpath{f});
    end
end
V2moviemaxTime  = 15000;   % voltage movie chunk length (frames)
GlumoviemaxTime = 3500;    % glutamate movie chunk length (frames)
o_laser         = 200001;  % 607/488 modulation onset (DAQ sample index, ~2 s after acq start)
mod488          = 200001;
foi             = [14:19];   % files of interest
set(0,'DefaultFigureWindowStyle','docked')

% --- extra columns + transform used by the f=3/f=5 sections below (cross-platform) ---
Endframe      = cell2mat(raw(:,5));
StructureData = raw(:,6);            % NB: paths in here also need macToWindowsPath on Windows
tformFile = '/Volumes/BHL18TB_D2/20260203_SD_V2+iGluSNFR4/25X_transformationMatrix.mat';
if iswindow, tformFile = macToWindowsPath(tformFile); end
load(tformFile);                     % tformReg
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
mov2_mc=double(readBinMov_times(fullfile(fpath{f}, ['mc2' num2str(1,'%02d') '.bin']), nRow2, nCol2,[1:nTime]));
catch
mov2_mc=double(readBinMov(fullfile(fpath{f}, ['mc2' num2str(1,'%02d') '.bin']), nRow2, nCol2));
end
load(fullfile(fpath{f}, ['mc2Trace' num2str(1,'%02d') '.mat']))
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

%% =====================================================================
%  SPIKE-TRIGGERED-AVERAGE + PLACE-TUNING PIPELINE  (files f = 14:19)
%  Soma voltage (Voltron-ST) + dendritic glutamate (iGluSNFR) + treadmill VR.
%
%  Per file, in order, run the cells below:
%    S1  detect simple (SS) + complex (CS) spikes on the soma voltage
%    S2  STA voltage trace  (SS- and CS-triggered)
%    S3  STA voltage movie  (streamed from mc*.bin voltage chunks)
%    S4  STA glutamate movie + per-ROI STA glutamate trace (voltage-spike-triggered)
%    S5  glutamate place tuning vs treadmill position (PlaceTrigger_average)
%
%  Inputs expected in each fpath{f} (produced by Analysis_20260713_SD_VoltronST_plus_Glu.m):
%    Volt_Result.mat  -> Result: .traces (nROI x T), .t_ax, .ref_im, .ftprnt,
%                                 .VmovTimesegments, .mc
%    Glu_Result.mat   -> fields : .S_glu, .dFF_glu (Nsyn x Tg), .t_ax,
%                                 .AvgGluImg, .bvMask, .RegressComponentTile
%    output_data.mat  -> Device_Data
%  Voltage bins  : mc%02d.bin  + mcTrace%02d.mat   (orange cam, Device_Data{3})
%  Glutamate bins: mc2%02d.bin + mc2Trace%02d.mat  (blue   cam, Device_Data{4})
%
%  Results are collected in STA_Result and saved as STA_Glu_Volt_Result.mat per file.
% =====================================================================
foi_STA = foi;   % from Section 0 (14:19)

% -- shared parameters --------------------------------------------------
% NB: o_laser, mod488 and GlumoviemaxTime are defined in Section 0 (cross-platform
%     setup) and are used as-is here — do NOT redefine them, or the glutamate tile
%     boundaries reconstructed in S4/S5 will not match how the movies were chunked.
staTauV   = [100 100];   % STA voltage window  [pre post] in voltage frames (~ms)
staTauG   = [15 15];     % STA glutamate window [pre post] in glutamate frames
matchTol  = 1;           % max |dt| (glu frames) allowed when mapping a V-spike to a glu frame

SS_thr    = [4 2];       % find_spike_bh [height prominence] (z-scored high-pass)
CS_onthr  = [5 1.5];     % detect_transient2 [on off] thresholds for complex-spike bursts
CS_minSp  = 3;           % >= this many spikes in a burst to call it a complex spike
CS_maxISI = 20;          % mean of first ISIs (frames) below this to call it a complex spike

% -- treadmill / place-tuning parameters (S5) ---------------------------
place_bin     = 100;     % # of position bins along the belt
vel_thresh    = 1;       % running threshold, in encoder ticks / glutamate-frame
% NOTE: beltLength_ticks is the belt circumference expressed in encoder ticks
%       (one lap).  It is rig-specific and NOT stored in the data, so SET IT here.
%       If a Virmen .data file is the real position source instead, tell me and
%       this section can be swapped to virmen_interpolate().
beltLength_ticks = 3000; % <-- CONFIRM/SET this per your treadmill belt

%% S1 : Detect simple + complex spikes on the soma voltage  (f = 14:19)
for f = foi_STA
    fprintf('=== S1 spike detection, file #%d ===\n', f);
    VoltResult = importdata(fullfile(fpath{f},'Volt_Result.mat'));
    [nROI, nT] = size(VoltResult.traces);

    % z-scored high-pass trace (spikes are positive-going in .traces)
    tr      = VoltResult.traces;
    tr_ref  = tr ./ get_threshold(tr,1);            % noise ~ 1 sigma
    tr_hi   = tr_ref - movmedian(tr_ref,50,2);

    spike   = zeros(nROI,nT);
    CStrace = zeros(nROI,nT);
    SpClass = zeros(2,nT,nROI);                     % row1 = SS(bAP), row2 = CS onset
    for n = 1:nROI
        % --- simple spikes ---
        sp = find_spike_bh(tr_hi(n,:), SS_thr(1), SS_thr(2));

        % --- complex-spike bursts (slow depolarizing envelope + >2 spikes) ---
        tr_sub = movmean(tr_ref(n,:),20,2);
        tr_sub = tr_sub - movmedian(tr_sub,300,2);
        [trans, tr_trace] = detect_transient2(tr_sub, CS_onthr, sp, 20);
        if isempty(trans.amp)
            CS_tr = zeros(1,nT);
        else
            hasISI   = cellfun(@(x) numel(x)>1, trans.ISI);
            meanISI  = inf(1,numel(trans.length));
            meanISI(hasISI) = cellfun(@(x) mean(x(1:2)), trans.ISI(hasISI));
            CS_ind   = find(trans.spike_number > CS_minSp-1 & meanISI < CS_maxISI);
            CS_tr    = ismember(tr_trace, CS_ind);
        end

        % trim each burst so it starts at its first spike (matches lab convention)
        bwCS = bwlabel(CS_tr);
        CS_spk = sp .* bwCS;
        [~, CS_spk_t] = unique(CS_spk); CS_spk_t(CS_spk_t<6) = [];
        for b = 1:max(bwCS)
            bfrm = find(bwCS==b);
            CS_tr(bfrm(1) : CS_spk_t(bwCS(CS_spk_t)==b)-5) = 0;
        end

        spike(n,:)     = sp;
        CStrace(n,:)   = CS_tr;
        SpClass(1,:,n) = sp .* (1-CS_tr);           % simple spikes (bAPs)
        SpClass(2,CS_spk_t(1:end),n) = 1;           % complex-spike onsets
        fprintf('  ROI %d: %d SS, %d CS\n', n, sum(SpClass(1,:,n)), sum(SpClass(2,:,n)));
    end

    VoltResult.spike   = spike;
    VoltResult.CStrace = CStrace;
    VoltResult.SpClass = SpClass;
    Result = VoltResult;  save(fullfile(fpath{f},'Volt_Result.mat'),'Result','-v7.3');

    % quick look (first ROI)
    figure(300+f); clf;
    plot(VoltResult.t_ax(1:nT), tr_ref(1,:)); hold on;
    onlyCS = nan(1,nT); onlyCS(CStrace(1,:)>0) = tr_ref(1,CStrace(1,:)>0);
    plot(VoltResult.t_ax(1:nT), onlyCS, 'r');
    plot(VoltResult.t_ax(find(SpClass(1,:,1))), tr_ref(1,SpClass(1,:,1)>0), 'k.');
    plot(VoltResult.t_ax(find(SpClass(2,:,1))), tr_ref(1,SpClass(2,:,1)>0), 'ro');
    title(sprintf('File %d  black=SS  red=CS',f)); xlabel('Time (s)'); axis tight; drawnow;
end

%% S2 : STA voltage trace  (SS- and CS-triggered, f = 14:19)
for f = foi_STA
    fprintf('=== S2 STA voltage trace, file #%d ===\n', f);
    VoltResult = importdata(fullfile(fpath{f},'Volt_Result.mat'));
    [nROI, nT] = size(VoltResult.traces);
    tr = VoltResult.traces - movmedian(VoltResult.traces,1000,2);
    tr = tr ./ get_threshold(tr,1);

    tV = (-staTauV(1):staTauV(2)) / 1000;           % s, assuming ~1 kHz voltage
    STA_V_SS = nan(nROI,sum(staTauV)+1);
    STA_V_CS = nan(nROI,sum(staTauV)+1);
    for n = 1:nROI
        ssIdx = find(VoltResult.SpClass(1,:,n) > 0);
        csIdx = find(VoltResult.SpClass(2,:,n) > 0);
        s = get_STA(tr(n,:), ssIdx, staTauV(1), staTauV(2)); if ~isempty(s), STA_V_SS(n,:) = s; end
        c = get_STA(tr(n,:), csIdx, staTauV(1), staTauV(2)); if ~isempty(c), STA_V_CS(n,:) = c; end
    end

    STA_Result = struct();
    STA_Result.tV = tV;
    STA_Result.STA_V_SS = STA_V_SS;
    STA_Result.STA_V_CS = STA_V_CS;
    save(fullfile(fpath{f},'STA_Glu_Volt_Result.mat'),'-struct','STA_Result','-v7.3');

    figure(320+f); clf; hold on;
    plot(tV, STA_V_SS','k'); plot(tV, STA_V_CS','r');
    xlabel('Time from spike (s)'); ylabel('Voltage (\sigma)');
    title(sprintf('File %d STA voltage: black=SS red=CS',f)); axis tight; drawnow;
end

%% S3 : STA voltage movie  (streamed from mc*.bin, f = 14:19)
for f = foi_STA
    fprintf('=== S3 STA voltage movie, file #%d ===\n', f);
    VoltResult = importdata(fullfile(fpath{f},'Volt_Result.mat'));
    load(fullfile(fpath{f},'output_data.mat'));                 % Device_Data
    sz         = double(Device_Data{3}.ROI([2 4]));            % [nCol nRow]
    readFrames = diff(VoltResult.VmovTimesegments);
    nSeg       = numel(readFrames);
    wV         = (-staTauV(1):staTauV(2));

    ssVec = VoltResult.SpClass(1,:,1) > 0;                     % ROI 1 = soma
    csVec = VoltResult.SpClass(2,:,1) > 0;
    H = sz(2); W = sz(1); nss = 0; ncs = 0;
    STA_ss = zeros(H,W,numel(wV)); STA_cs = zeros(H,W,numel(wV));
    off = 0;
    for j = 1:nSeg
        gFrames = off + (1:readFrames(j));                     % global frames of this chunk
        if ~any(ssVec(gFrames)) && ~any(csVec(gFrames))
            off = off + readFrames(j); continue;
        end
        movR = loadVoltChunk(fpath{f}, j, sz, readFrames(j));  % H x W x nT, spike-positive
        nT   = size(movR,3);
        locSS = find(ssVec(off+(1:nT)));  locCS = find(csVec(off+(1:nT)));
        for k = locSS(:)'
            if k+wV(1)>=1 && k+wV(end)<=nT
                STA_ss = STA_ss + movR(:,:,k+wV); nss = nss+1;
            end
        end
        for k = locCS(:)'
            if k+wV(1)>=1 && k+wV(end)<=nT
                STA_cs = STA_cs + movR(:,:,k+wV); ncs = ncs+1;
            end
        end
        off = off + nT;
        fprintf('  chunk %2d/%2d  (SS %d, CS %d)\n', j, nSeg, nss, ncs);
    end
    if nss>0, STA_ss = STA_ss/nss; end
    if ncs>0, STA_cs = STA_cs/ncs; end

    STA_Result = importdata(fullfile(fpath{f},'STA_Glu_Volt_Result.mat'));
    STA_Result.STAmov_V_SS = STA_ss;  STA_Result.STAmov_V_SS_n = nss;
    STA_Result.STAmov_V_CS = STA_cs;  STA_Result.STAmov_V_CS_n = ncs;
    save(fullfile(fpath{f},'STA_Glu_Volt_Result.mat'),'-struct','STA_Result','-v7.3');

    figure(340+f); clf; tiledlayout(1,2);
    nexttile; imshow2(max(STA_ss - median(STA_ss(:,:,1:15),3),[],3),[]); title(sprintf('File %d STA V SS (n=%d)',f,nss));
    nexttile; imshow2(max(STA_cs - median(STA_cs(:,:,1:15),3),[],3),[]); title(sprintf('CS (n=%d)',ncs)); drawnow;
end

%% S4 : STA glutamate movie + per-ROI STA glutamate trace (V-spike-triggered, f = 14:19)
for f = foi_STA
    fprintf('=== S4 STA glutamate, file #%d ===\n', f);
    VoltResult = importdata(fullfile(fpath{f},'Volt_Result.mat'));
    GluResult  = importdata(fullfile(fpath{f},'Glu_Result.mat'));
    load(fullfile(fpath{f},'output_data.mat'));                 % Device_Data
    nCol2 = double(Device_Data{4}.ROI(2)); nRow2 = double(Device_Data{4}.ROI(4));
    exposuretime2 = Device_Data{4}.exposuretime;

    % ---- reconstruct glutamate (blue) camera clock -> tile boundaries ----
    cam2_vsyn = Device_Data{1,2}.Counter_Inputs(1,2).data;
    cam2_vsyn_trig = find(diff(cam2_vsyn)==1)+1;
    segment_size2  = cam2_vsyn_trig(10) - cam2_vsyn_trig(9);
    start_idx = min(find(cam2_vsyn==max(cam2_vsyn)));
    end_idx   = length(Device_Data{1,2}.buffered_tasks(1,2).channels(1,2).data);
    last_val  = cam2_vsyn(start_idx-1); n_to_add = end_idx-start_idx+1;
    added = repelem((last_val+1 : last_val+ceil(n_to_add/segment_size2))', segment_size2);
    cam2_vsyn(start_idx:end_idx) = added(1:n_to_add);
    if cam2_vsyn(end) < GlumoviemaxTime
        GluSeg = [(cam2_vsyn(mod488)+2) cam2_vsyn(end)];
    else
        GluSeg = (cam2_vsyn(mod488)+2) : GlumoviemaxTime : cam2_vsyn(end);
        GluSeg(end+1) = cam2_vsyn(end);
    end
    Ntile   = numel(GluSeg)-1;
    tileLen = diff(GluSeg);
    nTg     = size(GluResult.dFF_glu,2);

    % ---- map soma V-spikes -> nearest glutamate frame ----
    tgt = VoltResult.t_ax(:); gax = GluResult.t_ax(1:nTg)';
    ssFrames = find(VoltResult.SpClass(1,:,1) > 0);   % soma SS
    csFrames = find(VoltResult.SpClass(2,:,1) > 0);   % soma CS
    map2glu = @(vf) local_nearestFrame(tgt(vf), gax, matchTol);
    gluIdx_SS = map2glu(ssFrames);
    gluIdx_CS = map2glu(csFrames);

    % ---- per-ROI STA glutamate TRACE (dF/F) ----
    Nsyn = size(GluResult.dFF_glu,1);
    tG   = (-staTauG(1):staTauG(2)) * exposuretime2;
    STA_glu_SS = get_STA(GluResult.dFF_glu, gluIdx_SS, staTauG(1), staTauG(2));
    STA_glu_CS = get_STA(GluResult.dFF_glu, gluIdx_CS, staTauG(1), staTauG(2));

    % ---- STA glutamate MOVIE (stream tiles, reprocess with stored PCs) ----
    wG = (-staTauG(1):staTauG(2));
    STAmov_SS = zeros(nRow2,nCol2,numel(wG)); nSS = 0;
    STAmov_CS = zeros(nRow2,nCol2,numel(wG)); nCS = 0;
    bvMask = max(GluResult.bvMask,[],3)==0;
    off = 0;
    for k = 1:Ntile
        gRange = off + (1:tileLen(k));
        inSS = gluIdx_SS(gluIdx_SS>=off+1 & gluIdx_SS<=off+tileLen(k)) - off;
        inCS = gluIdx_CS(gluIdx_CS>=off+1 & gluIdx_CS<=off+tileLen(k)) - off;
        if isempty(inSS) && isempty(inCS), off = off+tileLen(k); continue; end

        mov2 = double(readBinMov_times(fullfile(fpath{f}, ['mc2' num2str(k,'%02d') '.bin']), nRow2, nCol2, 1:tileLen(k)));
        load(fullfile(fpath{f}, ['mc2Trace' num2str(k,'%02d') '.mat']));   % mcTrace2
        meanF = squeeze(mean(mov2,[1 2]));
        y_fit = expfitDM_2((1:size(mov2,3))', meanF, (1:size(mov2,3))', [10000 100]);
        mov_f = mov2 ./ reshape(y_fit,1,1,[]);
        mov_lw= movmedian(mov_f, round(5/exposuretime2), 3);
        mov_s = mov_f - mov_lw;
        [V,~,~] = get_eigvector(tovec(imgaussfilt(mov_s,1))', 40);
        mov_s = SeeResiduals(mov_s, V(:,GluResult.RegressComponentTile{k}));
        mov_s = SeeResiduals(mov_s, mcTrace2.xymean(1:tileLen(k),:));
        mov_s = mov_s .* bvMask;

        for c = inSS(:)'
            if c+wG(1)>=1 && c+wG(end)<=tileLen(k), STAmov_SS = STAmov_SS + mov_s(:,:,c+wG); nSS = nSS+1; end
        end
        for c = inCS(:)'
            if c+wG(1)>=1 && c+wG(end)<=tileLen(k), STAmov_CS = STAmov_CS + mov_s(:,:,c+wG); nCS = nCS+1; end
        end
        off = off + tileLen(k);
        fprintf('  glu tile %2d/%2d  (SS %d, CS %d)\n', k, Ntile, nSS, nCS);
    end
    if nSS>0, STAmov_SS = STAmov_SS/nSS; end
    if nCS>0, STAmov_CS = STAmov_CS/nCS; end

    STA_Result = importdata(fullfile(fpath{f},'STA_Glu_Volt_Result.mat'));
    STA_Result.tG          = tG;
    STA_Result.STA_glu_SS  = STA_glu_SS;   % Nsyn x window, dF/F, SS-triggered
    STA_Result.STA_glu_CS  = STA_glu_CS;   % Nsyn x window, dF/F, CS-triggered
    STA_Result.STAmov_glu_SS = STAmov_SS;  STA_Result.STAmov_glu_SS_n = nSS;
    STA_Result.STAmov_glu_CS = STAmov_CS;  STA_Result.STAmov_glu_CS_n = nCS;
    save(fullfile(fpath{f},'STA_Glu_Volt_Result.mat'),'-struct','STA_Result','-v7.3');

    figure(360+f); clf; tiledlayout(2,2);
    nexttile; imagesc(tG, 1:Nsyn, STA_glu_SS); title(sprintf('File %d STA glu dF/F (SS)',f)); xlabel('t from spike (s)'); ylabel('synapse');
    nexttile; imagesc(tG, 1:Nsyn, STA_glu_CS); title('STA glu dF/F (CS)'); xlabel('t from spike (s)');
    nexttile; imshow2(max(STAmov_SS,[],3),[]); title(sprintf('STA glu movie SS (n=%d)',nSS));
    nexttile; imshow2(max(STAmov_CS,[],3),[]); title(sprintf('STA glu movie CS (n=%d)',nCS)); drawnow;
end

%% S5 : Glutamate place tuning vs treadmill position  (PlaceTrigger_average, f = 14:19)
for f = foi_STA
    fprintf('=== S5 glutamate place tuning, file #%d ===\n', f);
    GluResult = importdata(fullfile(fpath{f},'Glu_Result.mat'));
    load(fullfile(fpath{f},'output_data.mat'));                 % Device_Data
    nTg = size(GluResult.dFF_glu,2);

    % ---- reconstruct blue-cam clock (glutamate frame base) ----
    cam2_vsyn = Device_Data{1,2}.Counter_Inputs(1,2).data;
    cam2_vsyn_trig = find(diff(cam2_vsyn)==1)+1;
    segment_size2  = cam2_vsyn_trig(10) - cam2_vsyn_trig(9);
    start_idx = min(find(cam2_vsyn==max(cam2_vsyn)));
    end_idx   = length(Device_Data{1,2}.buffered_tasks(1,2).channels(1,2).data);
    last_val  = cam2_vsyn(start_idx-1); n_to_add = end_idx-start_idx+1;
    added = repelem((last_val+1 : last_val+ceil(n_to_add/segment_size2))', segment_size2);
    cam2_vsyn(start_idx:end_idx) = added(1:n_to_add);

    % ---- encoder ticks -> per-glutamate-frame speed ----
    enc = Device_Data{2}.buffered_tasks(3).channels(2).data;
    nUse = min(segment_size2*max(cam2_vsyn)+1, numel(enc));
    EncoderTick = double(diff(enc(1:nUse) > 2)) > 0;
    EncoderTick = EncoderTick(1 : segment_size2*floor(numel(EncoderTick)/segment_size2));
    Speed = sum(reshape(EncoderTick, segment_size2, []), 1);
    Speed = Speed(cam2_vsyn(mod488)+2 : end);                   % align onset to first glu frame
    if numel(Speed) >= nTg, Speed = Speed(1:nTg); else, Speed(end+1:nTg) = 0; end

    % ---- cumulative distance -> position (row5), lap (row8), velocity (last row) ----
    cumDist = cumsum(Speed);
    pos     = mod(cumDist, beltLength_ticks);
    lap     = floor(cumDist / beltLength_ticks) + 1;
    Virmen  = zeros(9, nTg);
    Virmen(5,:) = pos;
    Virmen(8,:) = lap;
    Virmen(9,:) = Speed;                                        % velocity, ticks/frame

    % ---- glutamate events (per synapse) ----
    GluSpike = find_spike_bh(zscore(GluResult.dFF_glu,0,2), 3, 1);

    % ---- lap x position maps: event rate, and mean dF/F ----
    [Lap_FR, Lap_SpkN, Lap_Nvalid] = PlaceTrigger_average(GluSpike,          place_bin, Virmen, vel_thresh, beltLength_ticks);
    [Lap_dFF]                      = PlaceTrigger_average(GluResult.dFF_glu,  place_bin, Virmen, vel_thresh, beltLength_ticks);

    % ---- spatial information per synapse (running frames only) ----
    binTr = ceil(Virmen(5,:)/(beltLength_ticks/place_bin));
    binTr(Virmen(9,:) < vel_thresh) = NaN;                     % drop non-running
    Nsyn  = size(GluResult.dFF_glu,1);
    SI = nan(Nsyn,1);
    for s = 1:Nsyn
        SI(s) = SpatialInfo(GluSpike(s,:), binTr);
    end

    STA_Result = importdata(fullfile(fpath{f},'STA_Glu_Volt_Result.mat'));
    STA_Result.Virmen     = Virmen;
    STA_Result.Lap_FR     = Lap_FR;      % nLap x place_bin x Nsyn (glu event rate)
    STA_Result.Lap_dFF    = Lap_dFF;     % nLap x place_bin x Nsyn (mean dF/F)
    STA_Result.Lap_SpkN   = Lap_SpkN;
    STA_Result.Lap_Nvalid = Lap_Nvalid;
    STA_Result.SpatialInfo = SI;
    STA_Result.place_bin  = place_bin;
    STA_Result.beltLength_ticks = beltLength_ticks;
    save(fullfile(fpath{f},'STA_Glu_Volt_Result.mat'),'-struct','STA_Result','-v7.3');

    % population place map (mean over laps) sorted by peak position
    popMap = squeeze(mean(Lap_FR,1,'omitnan'))';               % Nsyn x place_bin
    [~,pk] = max(ringmovMean(popMap,5),[],2); [~,ord] = sort(pk);
    figure(380+f); clf; tiledlayout(1,2);
    nexttile; imagesc(zscore_binning(popMap(ord,:),2)); xlabel('Position bin'); ylabel('synapse (sorted)');
    title(sprintf('File %d glu place map (peak-sorted)',f)); colorbar;
    nexttile; histogram(SI,30); xlabel('Spatial information (bits/event)'); ylabel('# synapses');
    title('Glu tuning strength'); drawnow;
end

%% ===================== Local functions (STA pipeline) =====================
function gluIdx = local_nearestFrame(spikeTimes, gluAxis, tol)
% Map each spike time (s) to the nearest glutamate-frame index; drop matches
% whose gap exceeds tol (in glutamate frames).
if isempty(spikeTimes), gluIdx = []; return; end
D = abs(gluAxis(:) - spikeTimes(:)');      % nGlu x nSpike
[mn, idx] = min(D,[],1);
dtFrame = mn ./ median(diff(gluAxis));
gluIdx  = idx(dtFrame < tol);
end