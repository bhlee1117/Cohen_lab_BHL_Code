cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20230917'
clear
[fpath] = uigetfile_n_dir(); %only Treadmill data
%% Motion correction

for i=1:length(fpath)
    disp(['Motion correction processing on ' fpath{i}])
    load([fpath{i} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    mov=double(readBinMov([fpath{i} '/frames1.bin'],sz(2),sz(1)));
    CamTrigger=find(Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1));
    mov_test=mov(:,:,250:350);

    try mov_test = single(mov_test)./single(max(mov_test.data(:)));
    catch disp('change to vm')
    mov_test=vm(mov_test); mov_test = single(mov_test)./single(max(mov_test.data(:))); end

    mov_test = movmean(mov_test,10,3);
    mov_ref = squeeze(median(mov_test,3));
    [mov_mc,xyField]=optical_flow_motion_correction_LBH(vm(mov(:,:,1:length(CamTrigger))),mov_ref, ...
        'normcorre');
    mov_mc=vm(mov_mc);
    mov_mc.transpose.savebin([fpath{i} '/mc.bin'])
    mcTrace=xyField;

    save([fpath{i} '/mcTrace.mat'],'mcTrace')
end

%% Soma Pulses Stim
i=1; 

    load([fpath{i} '/output_data.mat'])
    load([fpath{i} '/mcTrace.mat'],'mcTrace')
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    mov_mc=double(readBinMov([fpath{i} '/mc.bin'],sz(2),sz(1)));
    Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    CamTrigger=find(Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1));
    Blue=Blue(CamTrigger);

    nFrames=size(mov_mc,3);
    T_mean=squeeze(mean(mov_mc,[1 2]));
    avgImg=mean(mov_mc,3);
    bkTrace=get_blueoffTrace(T_mean,Blue,80);

    bkg = zeros(2, nFrames);
    bkg(1,:) = linspace(-1, 1, nFrames);  % linear term
    bkg(2,:) = linspace(-1, 1, nFrames).^2;  % quadratic term
    %bkg(3,:) = bkTrace;

    mov_res= SeeResiduals(mov_mc,bkg,1);
    mov_res = SeeResiduals(mov_res,mcTrace);
    mov_res = -SeeResiduals(mov_res,mcTrace.^2);
    mov_res=mov_res(:,:,10:end-10);
    [ROIs, F]=clicky(mov_res,avgImg);

    %%
nFrames2=size(mov_res,3);
datDS = imresize(mov_res(10:end-10,10:end-10,:), 0.3, 'bilinear', 'Antialiasing',true);
% datDS = imfilter(dat2, fspecial('average', [3, 3]), 'replicate');
% datDS = datDS(3:end-3,22:289,:);
% datDS = datDS(1:3:end,1:3:end,:);
datDS = datDS - medfilt1(datDS, 20, [], 3);

MovVec = tovec(datDS);
covMat = MovVec*MovVec';
[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;
eigTraces = V'*MovVec;
figure(8); clf
for i=1:10
plot(rescale(eigTraces(i,:)')+i-0.5)
hold all
end

nKeep = 10;
eigImgs = zeros(sz(2), sz(1), nKeep);
for j = 1:nKeep;
    eigImgs(:,:,j) = mean(mov_res.*reshape(eigTraces(j,:), [1, 1, nFrames2]),3);
end;
figure; clf;
for j = 1:nKeep;
    nexttile([1 1]);
    imshow2(eigImgs(3:end-3,3:end-3,j), []);
    title(num2str(j))
end;


keep_ind=[1:8];
[ics, mixmat, sepmat] = sorted_ica(eigTraces(keep_ind,:)',length(keep_ind));
figure(9); clf
stackplot(ics)

eigImgsVec = tovec(eigImgs(:,:,keep_ind));
footPrintsVec = eigImgsVec*sepmat';
footPrints = toimg(footPrintsVec, [sz(2), sz(1)]);
figure(10); clf
for j = 1:length(keep_ind);
    nexttile([1 1]);
    imshow2(footPrints(3:end-3,3:end-3,j), []);
    title(num2str(j))
end;

keep_ind_ics=[1 3];
ReconMov=toimg(tovec(mov_res).*mean(footPrintsVec(:,keep_ind_ics),2),sz(2),sz(1));

%% StimTriggered average
window=[-10:60]; bound=3;
Blue_crop=Blue(10:end-10)>0;
triggerFrame=find((Blue_crop(2:end)-Blue_crop(1:end-1)==1));
StimTA_mov=zeros(size(ReconMov,1)-(2*bound+1),size(ReconMov,2)-(2*bound+1),length(window));
g=1;
for i=1:length(triggerFrame)
    try
StimTA_mov=StimTA_mov+ReconMov(bound+1:end-bound-1,bound+1:end-bound-1,triggerFrame(i)+window);    
g=g+1;
    end
end
StimTA_mov=StimTA_mov./g;        
StimTA_mov=imgaussfilt3(StimTA_mov,[0.5 0.5 1]);
writeMov_wTrace([fpath{1} ,'/StimTA'],StimTA_mov,[],10,1,[-0.15 100],[],Blue_crop(triggerFrame(1)+window))

%% SpikeTriggered average
window=[-10:60];
[ROI, reconInt]=clicky(ReconMov,avgImg);
Blue_crop=Blue(10:end-10)>0;
ref_ind=1;
reconInt=reconInt(:,ref_ind)'./get_threshold(reconInt(:,ref_ind)',1);
spikes=find_spike_bh(reconInt,3,2);
figure; plot(reconInt)
hold all
plot(find(spikes),reconInt(find(spikes)),'r.')
Blue_bw=bwlabel(Blue_crop);
SpikeTA_mov=zeros(size(ReconMov,1),size(ReconMov,2),length(window));
STA=zeros(1,length(window));
g=0;
for i=1:max(Blue_bw)
bwSp=find((Blue_bw==i).*spikes);
if ~isempty(bwSp) && bwSp(1)+window(end)<size(ReconMov,3)
SpikeTA_mov=SpikeTA_mov+ReconMov(:,:,bwSp(1)+window);    
STA=STA+reconInt(bwSp(1)+window);
g=g+1;
end
end
SpikeTA_mov=SpikeTA_mov./g;        
STA=STA./g;

interpolationFactor = 10; % Doubling the resolution in the third dimension

originalSize = size(SpikeTA_mov);
interpolatedSize = [originalSize(1:2), originalSize(3) * interpolationFactor];

[Xq, Yq, Zq] = meshgrid(1:originalSize(2), 1:originalSize(1), linspace(1, originalSize(3), interpolatedSize(3)));

interpSpikeTA_mov = interp3(SpikeTA_mov, Xq, Yq, Zq, 'linear');
interpSpikeTA_mov=imgaussfilt3(interpSpikeTA_mov,[0.5 0.5 1]);
interSTA=interp1([1:length(STA)],STA,[1:originalSize(3) * interpolationFactor]/interpolationFactor);
writeMov_wTrace([fpath{1} ,'/SpikeTA_2'],interpSpikeTA_mov,[],50,10,[-0.10 200],[],interSTA)

%%
i=3; % dendrite

    load([fpath{i} '/output_data.mat'])
    load([fpath{i} '/mcTrace.mat'],'mcTrace')
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    mov_mc=double(readBinMov([fpath{i} '/mc.bin'],sz(2),sz(1)));
    Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    CamTrigger=find(Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1));
    Blue=Blue(CamTrigger);

    nFrames=size(mov_mc,3);
    T_mean=squeeze(mean(mov_mc,[1 2]));
    avgImg=mean(mov_mc,3);
    bkTrace=get_blueoffTrace(T_mean,Blue,80);

    bkg = zeros(2, nFrames);
    bkg(1,:) = linspace(-1, 1, nFrames);  % linear term
    bkg(2,:) = linspace(-1, 1, nFrames).^2;  % quadratic term
    %bkg(3,:) = bkTrace;

    mov_res= SeeResiduals(mov_mc,bkg,1);
    mov_res = SeeResiduals(mov_res,mcTrace);
    mov_res = -SeeResiduals(mov_res,mcTrace.^2);
    mov_res=mov_res(:,:,10:end-10);
    [ROIs, F]=clicky(mov_res,avgImg);

    %%
nFrames2=size(mov_res,3);
datDS = imresize(mov_res(10:end-10,10:end-10,:), 0.3, 'bilinear', 'Antialiasing',true);
% datDS = imfilter(dat2, fspecial('average', [3, 3]), 'replicate');
% datDS = datDS(3:end-3,22:289,:);
% datDS = datDS(1:3:end,1:3:end,:);
datDS = datDS - medfilt1(datDS, 20, [], 3);

MovVec = tovec(datDS);
covMat = MovVec*MovVec';
[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;
eigTraces = V'*MovVec;
figure(8); clf
nKeep = 20;

for i=1:nKeep
plot(rescale(eigTraces(i,:)')+i-0.5)
hold all
end


eigImgs = zeros(sz(2), sz(1), nKeep);
for j = 1:nKeep;
    eigImgs(:,:,j) = mean(mov_res.*reshape(eigTraces(j,:), [1, 1, nFrames2]),3);
end;
figure; clf;
for j = 1:nKeep;
    nexttile([1 1]);
    imshow2(eigImgs(3:end-3,3:end-3,j), []);
    title(num2str(j))
end;


keep_ind=[1 3 5 7 8];
[ics, mixmat, sepmat] = sorted_ica(eigTraces(keep_ind,:)',length(keep_ind));
figure(9); clf
stackplot(ics)

eigImgsVec = tovec(eigImgs(:,:,keep_ind));
footPrintsVec = eigImgsVec*sepmat';
footPrints = toimg(footPrintsVec, [sz(2), sz(1)]);
figure(10); clf
for j = 1:length(keep_ind);
    nexttile([1 1]);
    imshow2(footPrints(3:end-3,3:end-3,j), []);
    title(num2str(j))
end;

keep_ind_ics=[1 2 4 5];
ReconMov=toimg(tovec(mov_res).*mean(footPrintsVec(:,keep_ind_ics),2),sz(2),sz(1));

%% StimTriggered average
window=[-10:60]; bound=3;
Blue_crop=Blue(10:end-10)>0;
triggerFrame=find((Blue_crop(2:end)-Blue_crop(1:end-1)==1));
StimTA_mov=zeros(size(ReconMov,1)-(2*bound+1),size(ReconMov,2)-(2*bound+1),length(window));
g=1;
for i=1:length(triggerFrame)
    try
StimTA_mov=StimTA_mov+ReconMov(bound+1:end-bound-1,bound+1:end-bound-1,triggerFrame(i)+window);    
g=g+1;
    end
end
StimTA_mov=StimTA_mov./g;        
StimTA_mov=imgaussfilt3(StimTA_mov,[0.5 0.5 1]);
writeMov_wTrace([fpath{3} ,'/StimTA_dendriteStim'],StimTA_mov,[],10,1,[-0.1 15],[],Blue_crop(triggerFrame(1)+window))

%% SpikeTriggered average
window=[-10:60];
[ROI, reconInt]=clicky(mov_res,avgImg);
Blue_crop=Blue(10:end-10)>0;
Blue_crop_dilate = imdilate(Blue_crop, [ones(1,10) ones(1,20)]);

ref_ind=1;
reconInt=reconInt(:,ref_ind)'./get_threshold(reconInt(:,ref_ind)',1);
spikes=find_spike_bh(reconInt,5,2);
figure; plot(reconInt)
hold all
plot(find(spikes),reconInt(find(spikes)),'r.')


Blue_bw=bwlabel(Blue_crop_dilate);
SpikeTA_movCS=zeros(size(ReconMov,1),size(ReconMov,2),length(window));
SpikeTA_movSS=zeros(size(ReconMov,1),size(ReconMov,2),length(window));

STACS=zeros(1,length(window));
STASS=zeros(1,length(window));
triggerSS=[]; triggerCS=[];
g=0; gg=0;
for i=1:max(Blue_bw)
bwSp=find((Blue_bw==i).*spikes);
if ~isempty(bwSp) && bwSp(1)+window(end)<size(ReconMov,3)
    if length(bwSp)>1 

triggerCS=[triggerCS; bwSp(1)];
SpikeTA_movCS=SpikeTA_movCS+mov_res(:,:,bwSp(1)+window);    
STACS=STACS+reconInt(bwSp(1)+window);
g=g+1;
    else
triggerSS=[triggerSS; bwSp(1)];        
SpikeTA_movSS=SpikeTA_movSS+mov_res(:,:,bwSp(1)+window);    
STASS=STASS+reconInt(bwSp(1)+window);        
gg=gg+1;
    end
end
end

SpikeTA_movCS=SpikeTA_movCS./g;        
STACS=STACS./g;

SpikeTA_movSS=SpikeTA_movSS./gg;        
STASS=STASS./gg;

interpolationFactor = 5; % Doubling the resolution in the third dimension

originalSize = size(SpikeTA_movCS);
interpolatedSize = [originalSize(1:2), originalSize(3) * interpolationFactor];

[Xq, Yq, Zq] = meshgrid(1:originalSize(2), 1:originalSize(1), linspace(1, originalSize(3), interpolatedSize(3)));

interpSpikeTACS_mov = interp3(SpikeTA_movCS, Xq, Yq, Zq, 'linear');
interpSpikeTACS_mov=imgaussfilt3(interpSpikeTACS_mov,[0.5 0.5 1]);
interSTACS=interp1([1:length(STACS)],STACS,[1:originalSize(3) * interpolationFactor]/interpolationFactor);
figure;
writeMov_wTrace([fpath{3} ,'/SpikeTA_dendriteStimCS'],interpSpikeTACS_mov,[],50,5,[-1 30],[],interSTACS)

interpSpikeTASS_mov = interp3(SpikeTA_movSS, Xq, Yq, Zq, 'linear');
interpSpikeTASS_mov=imgaussfilt3(interpSpikeTASS_mov,[0.5 0.5 1]);
interSTASS=interp1([1:length(STASS)],STASS,[1:originalSize(3) * interpolationFactor]/interpolationFactor);
figure;
writeMov_wTrace([fpath{3} ,'/SpikeTA_dendriteStimSS'],interpSpikeTASS_mov,[],50,5,[-1 30],[],interSTASS)


%% interneurons
foi=10; % dendrite

    load([fpath{foi} '/output_data.mat'])
    load([fpath{foi} '/mcTrace.mat'],'mcTrace')
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    mov_mc=double(readBinMov([fpath{foi} '/mc.bin'],sz(2),sz(1)));
    Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    CamTrigger=find(Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1));
    Blue=Blue(CamTrigger);

    nFrames=size(mov_mc,3);
    T_mean=squeeze(mean(mov_mc,[1 2]));
    avgImg=mean(mov_mc,3);
    bkTrace=get_blueoffTrace(T_mean,Blue,80);

    bkg = zeros(2, nFrames);
    bkg(1,:) = linspace(-1, 1, nFrames);  % linear term
    bkg(2,:) = linspace(-1, 1, nFrames).^2;  % quadratic term
    %bkg(3,:) = bkTrace;

    mov_res= SeeResiduals(mov_mc,bkg,1);
    mov_res = SeeResiduals(mov_res,mcTrace);
    mov_res = -SeeResiduals(mov_res,mcTrace.^2);
    mov_res=mov_res(:,:,10:end-10);
    [ROIs, F]=clicky(mov_res,avgImg);

    %%
nFrames2=size(mov_res,3);
datDS = imresize(mov_res(10:end-10,10:end-10,:), 0.3, 'bilinear', 'Antialiasing',true);
% datDS = imfilter(dat2, fspecial('average', [3, 3]), 'replicate');
% datDS = datDS(3:end-3,22:289,:);
% datDS = datDS(1:3:end,1:3:end,:);
datDS = datDS - medfilt1(datDS, 20, [], 3);

MovVec = tovec(datDS);
covMat = MovVec*MovVec';
[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;
eigTraces = V'*MovVec;
figure(8); clf
nKeep = 20;

for i=1:nKeep
plot(rescale(eigTraces(i,:)')+i-0.5)
hold all
end


eigImgs = zeros(sz(2), sz(1), nKeep);
for j = 1:nKeep;
    eigImgs(:,:,j) = mean(mov_res.*reshape(eigTraces(j,:), [1, 1, nFrames2]),3);
end;
figure; clf;
for j = 1:nKeep;
    nexttile([1 1]);
    imshow2(eigImgs(3:end-3,3:end-3,j), []);
    title(num2str(j))
end;


keep_ind=[1:5];
[ics, mixmat, sepmat] = sorted_ica(eigTraces(keep_ind,:)',length(keep_ind));
figure(9); clf
stackplot(ics)

eigImgsVec = tovec(eigImgs(:,:,keep_ind));
footPrintsVec = eigImgsVec*sepmat';
footPrints = toimg(footPrintsVec, [sz(2), sz(1)]);
figure(10); clf
for j = 1:length(keep_ind);
    nexttile([1 1]);
    imshow2(footPrints(3:end-3,3:end-3,j), []);
    title(num2str(j))
end;

keep_ind_ics=[1 2];
ReconMov=toimg(tovec(mov_res).*mean(footPrintsVec(:,keep_ind_ics),2),sz(2),sz(1));

%% StimTriggered average
window=[-10:60]; bound=3;
Blue_crop=Blue(10:end-10)>0;
triggerFrame=find((Blue_crop(2:end)-Blue_crop(1:end-1)==1));
StimTA_mov=zeros(size(ReconMov,1)-(2*bound+1),size(ReconMov,2)-(2*bound+1),length(window));
g=1;
for i=1:length(triggerFrame)
    try
StimTA_mov=StimTA_mov+ReconMov(bound+1:end-bound-1,bound+1:end-bound-1,triggerFrame(i)+window);    
g=g+1;
    end
end
StimTA_mov=StimTA_mov./g;        
StimTA_mov=imgaussfilt3(StimTA_mov,[0.5 0.5 1]);
writeMov_wTrace([fpath{foi} ,'/StimTA_dendriteStim'],StimTA_mov,[],10,1,[-10 50],[],Blue_crop(triggerFrame(1)+window))

%% SpikeTriggered average
window=[-10:60];
[ROI, reconInt]=clicky(ReconMov,avgImg);
Blue_crop=Blue(10:end-10)>0;
ref_ind=1;
reconInt=reconInt(:,ref_ind)'./get_threshold(reconInt(:,ref_ind)',1);
spikes=find_spike_bh(reconInt,3,2);
figure; plot(reconInt)
hold all
plot(find(spikes),reconInt(find(spikes)),'r.')


Blue_bw=bwlabel(Blue_crop);
SpikeTA_mov=zeros(size(ReconMov,1),size(ReconMov,2),length(window));
STA=zeros(1,length(window));
g=0;
for i=1:max(Blue_bw)
bwSp=find((Blue_bw==i).*spikes);
if ~isempty(bwSp) && bwSp(1)+window(end)<size(ReconMov,3)
SpikeTA_mov=SpikeTA_mov+ReconMov(:,:,bwSp(1)+window);    
%SpikeTA_mov=SpikeTA_mov+mov_res(:,:,bwSp(1)+window);    
STA=STA+reconInt(bwSp(1)+window);
g=g+1;
end
end
SpikeTA_mov=SpikeTA_mov./g;        
STA=STA./g;

interpolationFactor = 10; % Doubling the resolution in the third dimension

originalSize = size(SpikeTA_mov);
interpolatedSize = [originalSize(1:2), originalSize(3) * interpolationFactor];

[Xq, Yq, Zq] = meshgrid(1:originalSize(2), 1:originalSize(1), linspace(1, originalSize(3), interpolatedSize(3)));

interpSpikeTA_mov = interp3(SpikeTA_mov, Xq, Yq, Zq, 'linear');
interpSpikeTA_mov=imgaussfilt3(interpSpikeTA_mov,[0.5 0.5 1]);
interSTA=interp1([1:length(STA)],STA,[1:originalSize(3) * interpolationFactor]/interpolationFactor);
 writeMov_wTrace([fpath{foi} ,'/SpikeTA_2'],interpSpikeTA_mov,[],50,10,[-10 30],[],interSTA)


