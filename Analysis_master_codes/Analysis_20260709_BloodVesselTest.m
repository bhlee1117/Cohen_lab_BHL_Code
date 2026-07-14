%% Load file path
clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:AA31');

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
NSeesawComponent=cell2mat(raw(:,25));
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
PlaceFieldList=cellfun(@(x) (str2num(num2str(x))),raw(:,21),'UniformOutput',false);
PlaceFieldBin=cellfun(@(x) (str2num(num2str(x))),raw(:,22),'UniformOutput',false);
set(0,'DefaultFigureWindowStyle','docked')
%foi=[1 4 5 6 8 10 11 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27];
foi=[1 4 5 6 8 10 11 15 16 17 18 19 20 21 22 23 24 25 26 27];

%%
actualframe=[78000:83000];
fpath2readmov='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20240807/201405BHLm141_N1_VR_LowStim/';
f=20;
load(fullfile(fpath{f},'PC_Result.mat'));
load(fullfile(fpath2readmov,'mcTrace_seg.mat'))
load(fullfile(fpath2readmov,'output_data.mat'))
sz=double(Device_Data{1, 3}.ROI([2 4]));

mov_mc=double(readBinMov(fullfile(fpath2readmov,'mov_mc_segment.bin'),sz(2),sz(1)));
bv_trace=tovec(mov_mc)'*tovec(Result.bvMask);
V_trace=tovec(mov_mc)'*tovec(Result.ftprnt);
bv_trace=SeeResiduals(permute(bv_trace,[2 3 1]),V_trace);
bv_trace=squeeze(bv_trace)';

mov_res= mov_mc-mean(mov_mc,3);
mov_res = SeeResiduals(mov_res,mcTrace);
mov_res = SeeResiduals(mov_res,mcTrace.^2);
mov_res = SeeResiduals(mov_res,mcTrace(:,1).*mcTrace(:,2));
mov_res2 = SeeResiduals(mov_res,bv_trace);
% mov_res = SeeResiduals(mov_res,bv_trace);
% F0img=get_F0img_PCA(movmean(mov_res,20,3));

% [~, mov_res_perispike]=get_STA(tovec(mov_res),Result.spike(1,78000+[0:size(mov_res,3)-1]),3,3);
% mov_res_perispike=reshape(permute(mov_res_perispike,[1 3 2]),sz(2),sz(1),[]);
% mov_res_sub=movmean(mov_res,20,3);
% [u,s,v] = svds(tovec(mov_res_perispike),20);
% reshape_u=reshape(u,sz(2),sz(1),[]);
% figure(22); clf; imshow2_patch(reshape_u);
% [u_sub,s,v] = svds(tovec(mov_res_sub),20);
% reshape_u_sub=reshape(u_sub,sz(2),sz(1),[]);
% figure(23); clf; imshow2_patch(reshape_u_sub);
% pc2use=[1:10];
% [mov_resfilt] = pcafilt_template(mov_res, cat(3,reshape_u(:,:,pc2use),reshape_u_sub(:,:,[4 5])));

% ---- 당신 무비 준비 ----
% mov: [Ly, Lx, T] (세로, 가로, 프레임)
% 일부만 테스트하려면 프레임을 자르세요 (최소 1024프레임 이상 권장! 아래 주의 참고)
mov_mc_sub=mov_mc(:,:,1:2000)-mean(mov_mc(:,:,1:2000),3);

%%
[~, mov_res_perispike]=get_STA(tovec(mov_res),Result.spike(1,78000+[0:size(mov_res,3)-1]),3,3);
mov_res_perispike=reshape(permute(mov_res_perispike,[1 3 2]),sz(2),sz(1),[]);
mov_res_sub=movmean(mov_res,20,3);
[u,s,v] = svds(tovec(mov_res_perispike),20);
reshape_u=reshape(u,sz(2),sz(1),[]);
figure(22); clf; imshow2_patch(reshape_u);
[u_sub,s,v] = svds(tovec(mov_res_sub),20);
reshape_u_sub=reshape(u_sub,sz(2),sz(1),[]);
figure(23); clf; imshow2_patch(reshape_u_sub);
pc2use=[1:10];
[mov_resfilt] = pcafilt_template(mov_res2(:,:,1:2000), cat(3,reshape_u(:,:,pc2use),reshape_u_sub(:,:,[4 5])));

%% Run blood vessel correction

venvpy  = '/Users/bhlee1117/Documents/GitHub/bleach-dr/.venv/bin/python';
script  = '/Users/bhlee1117/Documents/GitHub/bleach-dr/run_correction.py';
in_h5  = '/Users/bhlee1117/Documents/GitHub/bleach-dr/data/my_movie.h5';
out_h5 = '/Users/bhlee1117/Documents/GitHub/bleach-dr/data/my_movie_corrected.h5';

% --- 1. write movie to HDF5 (mov: [Ly,Lx,T]) ---
mov = single(-mov_res(:, :, 1:2000));          % your movie, >=1024 frames
mov_out = permute(mov, [2 1 3]);          % -> Python reads (T, Ly, Lx)
if exist(in_h5,'file'), delete(in_h5); end
h5create(in_h5, '/mov', size(mov_out), 'Datatype', 'single');
h5write(in_h5,  '/mov', mov_out);

% --- 2. run the Python correction ---
cmd = sprintf('"%s" "%s" --input "%s" --output "%s" --frame-rate 1000 --rank1-method skip', ...
              venvpy, script, in_h5, out_h5);
[status, out] = system(cmd, '-echo');     % prints live progress
assert(status == 0, 'Python failed:\n%s', out);

% --- 3. read the result back (undo the permute) ---
corrected = permute(h5read(out_h5, '/corrected'), [2 1 3]);
blood     = permute(h5read(out_h5, '/blood'),     [2 1 3]);
%%
figure(3);
mov2write= [[double(imgaussfilt(mov,1)) -mov_resfilt]; corrected blood];
moviefixsc(mov2write)
writefile='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20240807/201405BHLm141_N1_VR_LowStim/BloodVesselCorrectionCompare.mp4';
writeMov(writefile,mov2write(:,:,1:1000),[1:1000],10,1,[-70 100])
%%
tr_corrected=(tovec(double(corrected))'*tovec(Result.ftprnt))';
figure(2); clf; ax1=[];
ax1=[ax1 nexttile([1 1])];
imagesc(rescale2(Result.normTraces(Result.dist_order,actualframe(1:2000)).*(Result.dirtTrace(Result.dist_order,actualframe(1:2000))==0),2))
title('Signal from -∆F movie')
ax1=[ax1 nexttile([1 1])];
imagesc(rescale2(tr_corrected(Result.dist_order,:).*(Result.dirtTrace(Result.dist_order,actualframe(1:2000))==0),2));
title('Signal from blood vessel corrected movie')
colormap(turbo);
linkaxes(ax1,'xy')
xlim([0 1000]);
set_figsize(280,215)
