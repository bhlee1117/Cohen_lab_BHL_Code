clear;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20240224_BHLm112_BHLm115_VR/Bloodvessel_Correction'
load(['output_data.mat'])
load("mcTrace01.mat")
sz=double(Device_Data{1, 3}.ROI([2 4]));
mov_mc=double(readBinMov(['mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1)));
%%
bkg = zeros(2, size(mov_mc,3));
bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term

mc=mcTrace.xymean;
mov_res=mov_mc-mean(mov_mc,3);
mov_res = SeeResiduals(mov_res,mc);
mov_res = SeeResiduals(mov_res,mc.^2);
mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
mov_res= SeeResiduals(mov_res,bkg,1);

%%

mov_lag=mov_res(:,:,1:end-1).*mov_res(:,:,2:end);

%%
silent_time=[4661:5612];
mov_lag_silent=mov_res(:,:,silent_time).*mov_res(:,:,silent_time+1);
imshow2(std(mov_lag_silent,0,3),[])