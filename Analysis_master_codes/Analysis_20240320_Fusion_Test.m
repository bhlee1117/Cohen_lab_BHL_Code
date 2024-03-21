clear; clc;
fpath='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20240224_BHLm112_BHLm115_VR/125644BHLm115_VR2';
cd(fpath);
fname=fullfile(fpath,'125644BHLm115_VR2/frames1.bin');

%%

toi=[13500:15000];

load(fullfile(fpath,"output_data.mat"))
sz=double(Device_Data{1, 3}.ROI([2 4]));
ref_time=[9000:10000];
mov=double(readBinMov_times([fpath '/frames1.bin'],sz(2),sz(1),toi));
mov_roll=rollingShutter_correction(mov,Device_Data{1, 3}.exposuretime,'fusion');

mov_roll=mov_roll(:,:,2:end);

[mov_roll_mc,xyField_roll]=optical_flow_motion_correction_LBH(mov_roll,mean(mov,3),'normcorre');
[mov_mc,xyField]=optical_flow_motion_correction_LBH(mov,mean(mov,3),'normcorre');
avgImg=mean(mov_mc,3);

mc_roll=xyField_roll.xymean;
mc=xyField.xymean;

mov_roll_res= mov_roll_mc-mean(mov_roll_mc,3);
bkg = zeros(2, size(mov_roll_mc,3));
bkg(1,:) = linspace(-1, 1, size(mov_roll_mc,3));  % linear term
bkg(2,:) = linspace(-1, 1, size(mov_roll_mc,3)).^2;  % quadratic term
mov_roll_res = SeeResiduals(mov_roll_res,mc_roll);
mov_roll_res = SeeResiduals(mov_roll_res,mc_roll.^2);
mov_roll_res = SeeResiduals(mov_roll_res,mc_roll(:,1).*mc_roll(:,end));
mov_roll_res= SeeResiduals(mov_roll_res,bkg,1);

mov_res= mov_mc-mean(mov_mc,3);
bkg = zeros(2, size(mov_mc,3));
bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
mov_res = SeeResiduals(mov_res,mc);
mov_res = SeeResiduals(mov_res,mc.^2);
mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
mov_res= SeeResiduals(mov_res,bkg,1);
%%
[roi, tr_roll]=clicky(mov_roll_res,avgImg);
[tr_mc]=apply_clicky(roi,mov_res(:,:,2:end));
[tr_raw]=apply_clicky(roi,mov(:,:,2:end));
[tr_raw_roll]=apply_clicky(roi,mov_roll);

%%
figure(10); clf; roiN=[1 3]; t_ref=[1100:1170];
nexttile([1 1])
imshow2(std(mov(:,:,2:end)-mov_roll,0,3),[]);
nexttile([1 1])
plot([1:size(mov,1)],(mean(std(mov(:,300:end,2:end)-mov_roll(:,300:end,:),0,3),2)));
nexttile([1 1])
plot(-rescale2(tr_roll(t_ref,roiN),1)+[1 2]);
title('Rolling Shutter Corrected, Motion Corrected')

nexttile([1 1])
plot(-rescale2(tr_mc(t_ref,roiN),1)+[1 2]);
title('Motion Corrected')

nexttile([1 1])
plot(-rescale2(tr_raw_roll(t_ref,roiN),1)+[1 2]);
title('Rolling Shutter Corrected')

nexttile([1 1])
plot(-rescale2(tr_raw(t_ref,roiN),1)+[1 2]);
title('Raw')
