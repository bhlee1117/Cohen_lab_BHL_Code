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
%% Load movie and find the footprints
f=13;
Devicedata_filename=fullfile(fpath{f},'output_data.mat');
Vdata_filename=fullfile(fpath{f},'frames1.bin');
Gludata_filename=fullfile(fpath{f},'frames2.bin');
load(Devicedata_filename);

% set orange cam
nCol1 = double(Device_Data{3}.ROI([2]));  % ROI on the camera
nRow1 = double(Device_Data{3}.ROI([4]));  % ROI on the camera
exposuretime1 = Device_Data{3}.exposuretime;
% ---------- Correct spinning disk artifact ----------
% Discard first some second of frames without 607 ex light, drop last 1s to avoid sync issues
o_laser = 200001; % all the recording started 2s after the acquisition start

% voltage camera clock
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
segment_size2 = cam2_vsyn_trig (10) - cam2_vsyn_trig(9);
start_idx = min(find(cam2_vsyn ==max (cam2_vsyn)));
end_idx = length(Device_Data{1, 2}.buffered_tasks(1, 2).channels(1, 2).data);
last_val = cam2_vsyn(start_idx - 1);
n_to_add = end_idx - start_idx + 1;
n_segments = ceil(n_to_add / segment_size2);
added_part = repelem((last_val + 1 : last_val + n_segments)', segment_size2);
added_part = added_part(1:n_to_add);
cam2_vsyn(start_idx:end_idx) = added_part;

mod488 = 200001;
mod607 = 200001;

cam2_trig = find (diff(cam2_vsyn) ==1)+1; %DAQ frame
%cam2_trig = cam1_vsyn (cam2_trig); %Cam1 frame
t_start= cam2_trig (cam2_vsyn(mod488)+2)-mod488+1;
if cam2_vsyn(end)<GlumoviemaxTime
    GluemovTimesegments=[(cam2_vsyn(mod488)+2) cam2_vsyn(end)];
else
    GluemovTimesegments=[(cam2_vsyn(mod488)+2):GlumoviemaxTime:cam2_vsyn(end)];
    GluemovTimesegments(end+1)=cam2_vsyn(end);
end
nFrame2analyze2=GluemovTimesegments(end)-GluemovTimesegments(1);
t_glu=[t_start:segment_size2:(t_start+segment_size2*nFrame2analyze2)]/DAQ_rate;

Blue_Nframe = Device_Data{4}.frames_requested;
nCol2 = double(Device_Data{4}.ROI([2]));  % ROI on the camera
nRow2 = double(Device_Data{4}.ROI([4]));  % ROI on the camera

Vmov_read_tmp=double(readBinMov_times([fpath{f} '/mc' num2str(1,'%02d') '.bin'], nRow1, nCol1,[2000:3000]));
AvgVoltageImg=mean(Vmov_read_tmp,3);

Ntile=length(GluemovTimesegments)-1;
opts = struct();
opts.doPlot = true;
T=GluemovTimesegments(end)-GluemovTimesegments(1);
GluResult = struct();

%% Extract GluResult
% ---------- PASS 1: detect per chunk ----------
tileResults = cell(Ntile,1);
RegressComponentTile=cell(Ntile,1);
mcTrace_glu=[];
for k = 1:Ntile
    fprintf('Reading %1.0f out of %1.0f movie...\n', k,length(GluemovTimesegments)-1)

    mov2_mc=double(readBinMov_times([fpath{f} '/mc2' num2str(k,'%02d') '.bin'], nRow2, nCol2,[1:(GluemovTimesegments(k +1)-GluemovTimesegments(k))]));
    load([fpath{f} '/mc2Trace' num2str(k,'%02d') '.mat'])
    meanF=squeeze(mean(mov2_mc,[1 2]));
    y_fit=expfitDM_2([1:size(mov2_mc,3)]',meanF,[1:size(mov2_mc,3)]',[10000 100]);
    mov2_mc_filt=mov2_mc./reshape(y_fit,1,1,[]);
    mov2_mc_filt_lw=movmedian(mov2_mc_filt,round(5/exposuretime2),3);
    mov_sub=mov2_mc_filt-mov2_mc_filt_lw;
    
    [V, D, u_sub] = get_eigvector(tovec(imgaussfilt(mov_sub,1))',40);
    u_sub=reshape(u_sub,nRow2(1),nCol2(1),[]);

    if ~isfield(GluResult,'bvMask')
    %[~, GluResult.bvMask]=get_ROI(max(abs(u_sub),[],3));
    [~, GluResult.bvMask]=get_ROI(mean(mov2_mc,3));
    end
    figure(2); clf; imshow2_patch(u_sub); drawnow;
    n=input("PCs to regress out\n");
    RegressComponentTile{k}=n;
    GluResult.RegressComponentTile{k}=RegressComponentTile{k};

    mov_sub=SeeResiduals(mov_sub,V(:,n));
    mov_sub=SeeResiduals(mov_sub,mcTrace2.xymean([1:(GluemovTimesegments(k +1)-GluemovTimesegments(k))],:));

    mcTrace_glu=[mcTrace_glu; mcTrace2.xymean([1:(GluemovTimesegments(k +1)-GluemovTimesegments(k))],:)];
    mov2_mc_filt2=mov_sub+mov2_mc_filt_lw;
    mov_sub=mov_sub.*(max(GluResult.bvMask,[],3)==0);

    if ~isfield(tileResults{1},'opts') %first trial
        opts.diagSaveDir       = fpath{f};  % where PNGs are saved
        opts.diagPrefix        = 'exp01';                 % filename prefix
        opts.diagZoom_pad      = 15;                      % px around centroid for zoomed plots
        opts.diagMaxSyn        = 16;
        opts.eventMaxArea    = 30;   % px² threshold for NMF split
        opts.nmfReps         = 5;
        opts.nmfDiagMaxEvent = 6;
        opts.showNMFDiag     = true;
        
    else
        opts=tileResults{1}.opts;
        opts.uiClusterSize=0;
        opts.plotCentroids=0;
        opts.uiSynFilter=1;
        opts.uiSeedThresh=1;
    end

    tileResults{k} = extractGluSNFR3(mov_sub, mov2_mc_filt2, exposuretime2, opts);  % detect mode default
    drawnow;
end

if Ntile>1
Dict = mergeSynapseDict(tileResults, 3);
S_full = Dict.S_glu;     % H×W×Nsyn_global

opts2 = struct();
opts2.mode = 'project';

for k = 1:Ntile
    fprintf('Reading %1.0f out of %1.0f movie...\n', k,length(GluemovTimesegments)-1)
    try
        mov2_mc=double(readBinMov([fpath{f} '/mc2' num2str(k,'%02d') '.bin'], nRow2, nCol2));
    catch
        mov2_mc=double(readBinMov_times([fpath{f} '/mc2' num2str(k,'%02d') '.bin'], nRow2, nCol2,[1:(GluemovTimesegments(k+1)-GluemovTimesegments(k))]));
    end
    meanF=squeeze(mean(mov2_mc,[1 2]));
    y_fit=expfitDM_2([1:size(mov2_mc,3)]',meanF,[1:size(mov2_mc,3)]',[10000 100]);
    mov2_mc_filt=mov2_mc./reshape(y_fit,1,1,[]);
    mov2_mc_filt_lw=movmedian(mov2_mc_filt,100,3);
    mov_sub=mov2_mc_filt-mov2_mc_filt_lw;
    [V, D, u_sub] = get_eigvector(tovec(imgaussfilt(mov_sub,1))',40);
    mov_sub=SeeResiduals(mov_sub,V(:,RegressComponentTile{k}));
    mov2_mc_filt2=mov_sub+mov2_mc_filt_lw;
    Rproj{k} = extractGluSNFR3(mov_sub, mov2_mc_filt2, exposuretime2, opts2, S_full);
end
else
    Rproj=tileResults;
    S_full =Rproj{1,1}.S_glu;
end

frameRate = 1/exposuretime2;

% Find first non-empty
first = find(~cellfun(@isempty, Rproj), 1, 'first');
if isempty(first)
    GluResult = struct('S_glu', S_full, 'T_glu', [], 'F0_glu', [], 'dFF_glu', [], 'frameRate', frameRate);
    return
end

Nsyn = size(S_full,3);

% Concatenate along time
T_all = [];
F0_all = [];
dFF_all = [];

for k = 1:numel(Rproj)
    R = Rproj{k};
    if isempty(R) || isempty(R.T_glu)
        continue
    end

    % sanity: synapse count matches
    assert(size(R.T_glu,1) == Nsyn, 'Chunk %d has wrong Nsyn.', k);

    T_all  = [T_all,  R.T_glu]; %#ok<AGROW>

    if isfield(R,'F0_glu') && ~isempty(R.F0_glu)
        F0_all = [F0_all, R.F0_glu]; %#ok<AGROW>
    end

    if isfield(R,'dFF_glu') && ~isempty(R.dFF_glu)
        dFF_all = [dFF_all, R.dFF_glu]; %#ok<AGROW>
    end
end

GluResult.mc= mcTrace_glu;
GluResult.S_glu = S_full;
GluResult.T_glu = T_all;
GluResult.F0_glu = F0_all;
GluResult.dFF_glu = dFF_all;
GluResult.frameRate = frameRate;
GluResult.AvgGluImg = mean(mov2_mc,3);
GluResult.t_ax=t_glu;

% Save and show
if exist('fpath','var') && exist('f','var')
    save(fullfile(fpath{f}, 'Glu_Result.mat'), '-struct', 'GluResult', '-v7.3');
    fprintf('Saved Glu_Result.mat to %s\n', fpath{f});
end

% Quick plots for sanity
tax = (0:T-1)/frameRate;
glucoord=round(get_coord(GluResult.S_glu));
[alignedGlu, Rglu, Rvolt] = transformCamera(Device_Data, tformReg, AvgVoltageImg, max(S_full,[],3));
alignedVolt = transformCamera_O2B(Device_Data, tformReg, AvgVoltageImg, GluResult.AvgGluImg);

validSynapse=zeros(1,size(glucoord,1));
for s=1:size(glucoord,1)
validSynapse(s)=alignedVolt(glucoord(s,2),glucoord(s,1))>0;
end

figure(99); clf;
line_color=hsv(sum(validSynapse));
l=plot(GluResult.t_ax,dFF_all(validSynapse>0,:)'+[1:sum(validSynapse)]*1);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(line_color,2))
title('dF/F');
xlabel('Time (s)'); ylabel('Spine ID'); axis tight;
saveas(gca,fullfile(fpath{f},'Glu_trace.fig'));

figure(100); clf;
offsetscale=0.7;
coloredS_glu=showScaleImage(S_full(:,:,validSynapse>0),[1:sum(validSynapse)],'hsv');

nexttile([1 1]);
imshow2(coloredS_glu+ grs2rgb(alignedVolt,gray(256))*offsetscale,[])

nexttile([1 1]);
imshow2(im_merge(cat(3,mat2gray(alignedVolt),mat2gray(max(S_full,[],3))),[0 0.6 1; 1 0 0]),[])

nexttile([1 1]);
imshow2(im_merge(cat(3,mat2gray(GluResult.AvgGluImg),mat2gray(max(S_full,[],3))),[0 0.6 1; 1 0 0]),[])
saveas(gca,fullfile(fpath{f},'Glu_footprint.png'));

% %% Show Orange and blue average image
% for f=2:7
% 
%     Devicedata_filename=fullfile(fpath{f},'output_data.mat');
% Vdata_filename=fullfile(fpath{f},'frames1.bin');
% Gludata_filename=fullfile(fpath{f},'frames2.bin');
% load(Devicedata_filename);
% 
% % set orange cam
% nCol1 = double(Device_Data{3}.ROI([2]));  % ROI on the camera
% nRow1 = double(Device_Data{3}.ROI([4]));  % ROI on the camera
% exposuretime1 = Device_Data{3}.exposuretime;
% % ---------- Correct spinning disk artifact ----------
% % Discard first some second of frames without 607 ex light, drop last 1s to avoid sync issues
% o_laser = 200001; % all the recording started 2s after the acquisition start
% 
% % voltage camera clock
% cam1_vsyn = Device_Data{1, 2}.Counter_Inputs(1, 1).data;
% start_idx = min(find(cam1_vsyn ==max (cam1_vsyn)));
% end_idx = length(Device_Data{1, 2}.buffered_tasks(1, 2).channels(1, 1).data);
% segment_size = 103*floor(exposuretime1/0.001);
% last_val = cam1_vsyn(start_idx - 1);
% n_to_add = end_idx - start_idx + 1;
% n_segments = ceil(n_to_add / segment_size);
% added_part = repelem((last_val + 1 : last_val + n_segments)', segment_size);
% added_part = added_part(1:n_to_add);
% cam1_vsyn(start_idx:end_idx) = added_part;
% VmovTimesegments=[(cam1_vsyn(o_laser)+2):V2moviemaxTime:cam1_vsyn(end)];
% nFrame2analyze=VmovTimesegments(end)-VmovTimesegments(1)+1;
% t_vol = (0:nFrame2analyze-1)*exposuretime1;
% 
% cam2_vsyn = Device_Data{1, 2}.Counter_Inputs(1, 2).data;
% cam2_vsyn_trig = find (diff(cam2_vsyn) == 1)+1;
% segment_size = cam2_vsyn_trig (10) - cam2_vsyn_trig(9);
% start_idx = min(find(cam2_vsyn ==max (cam2_vsyn)));
% end_idx = length(Device_Data{1, 2}.buffered_tasks(1, 2).channels(1, 2).data);
% last_val = cam2_vsyn(start_idx - 1);
% n_to_add = end_idx - start_idx + 1;
% n_segments = ceil(n_to_add / segment_size);
% added_part = repelem((last_val + 1 : last_val + n_segments)', segment_size);
% added_part = added_part(1:n_to_add);
% cam2_vsyn(start_idx:end_idx) = added_part;
% 
% mod488 = 200001;
% mod607 = 200001;
% cam2_trig = find (diff(cam2_vsyn) ==1)+1;
% cam2_trig = cam1_vsyn (cam2_trig);
% cam2_trig = cam2_trig (cam2_vsyn(mod488)+2:cam2_vsyn(end))'- (cam1_vsyn(mod607)+2);
% if cam2_vsyn(end)<GlumoviemaxTime
%     GluemovTimesegments=[(cam2_vsyn(mod488)+2) cam2_vsyn(end)];
% else
%     GluemovTimesegments=[(cam2_vsyn(mod488)+2):GlumoviemaxTime:cam2_vsyn(end)];
% end
% 
% Blue_Nframe = Device_Data{4}.frames_requested;
% nCol2 = double(Device_Data{4}.ROI([2]));  % ROI on the camera
% nRow2 = double(Device_Data{4}.ROI([4]));  % ROI on the camera
% exposuretime2 = Device_Data{4}.exposuretime;
% t_cal = (cam2_trig-1)*exposuretime2;
% 
% Vmov_read_tmp=double(readBinMov_times([fpath{f} '/mc' num2str(1,'%02d') '.bin'], nRow1, nCol1,[2000:3000]));
% AvgVoltageImg=mean(Vmov_read_tmp,3);    
% 
% Glumov_read_tmp=double(readBinMov_times([fpath{f} '/mc2' num2str(1,'%02d') '.bin'], nRow2, nCol2,[500:500+round(5/exposuretime2)]));
% AvgGluImg=sum(Glumov_read_tmp,3);
% AvgGluImg_trf=transformCamera(Device_Data,tformReg,AvgVoltageImg,AvgGluImg);
% 
% figure(f); clf;
% nexttile([1 1])
% imshow2(AvgVoltageImg,[5 70])
% nexttile([1 1])
% imshow2(AvgGluImg_trf,[2000 20000])
% end