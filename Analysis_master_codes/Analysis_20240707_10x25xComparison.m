clear; clc;
fpth25x='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20240705_Optopatch_Prism/20240704/205405BHLm144_N4_SomRP_50mm_25x_Or3V_Blue3mW';
fpth10x='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20240705_Optopatch_Prism/20240704/220049BHLm144_N4_SomRP_180mm_10x_Or5V_Blue7mW';

%% MC 25x movie
load([fpth25x '/output_data.mat']);
sz=double(Device_Data{1, 3}.ROI([2 4]));
mov=double(readBinMov([fpth25x '/frames1.bin'],sz(2),sz(1)));

DAQ_rate=Device_Data{1, 2}.buffered_tasks(1, 1).rate;
CamDAQ_rate=Device_Data{1, 2}.Counter_Inputs.rate;
CamTrig=Device_Data{1, 2}.Counter_Inputs.data;
CamTrig2=find(CamTrig(2:end)-CamTrig(1:end-1)>0);
Frm_rate=(CamTrig2(2)-CamTrig2(1))/CamDAQ_rate;

mov=rollingShutter_correction(mov,1/Frm_rate,'fusion');
mov=mov(:,:,2:end);
[mov_mc_25x,xyField]=optical_flow_motion_correction_LBH(mov,mean(mov,3),'normcorre');
avg_im=mean(mov_mc_25x,3);
mov_mc_25x=vm(mov_mc_25x);
mov_mc_25x.transpose.savebin([fpth25x '/mc_ShutterReg' '.bin'])
mcTrace=xyField; % Normcorre

save([fpth25x '/mcTrace.mat'],'mcTrace','avg_im')

%% MC 10x movie

load([fpth10x '/output_data.mat']);
sz=double(Device_Data{1, 3}.ROI([2 4]));
mov=double(readBinMov([fpth10x '/frames1.bin'],sz(2),sz(1)));

DAQ_rate=Device_Data{1, 2}.buffered_tasks(1, 1).rate;
CamDAQ_rate=Device_Data{1, 2}.Counter_Inputs.rate;
CamTrig=Device_Data{1, 2}.Counter_Inputs.data;
CamTrig2=find(CamTrig(2:end)-CamTrig(1:end-1)>0);
Frm_rate=(CamTrig2(2)-CamTrig2(1))/CamDAQ_rate;

mov=rollingShutter_correction(mov,1/Frm_rate,'fusion');
mov=mov(:,:,2:end);
[mov_mc_10x,xyField]=optical_flow_motion_correction_LBH(mov,mean(mov,3),'normcorre');
avg_im=mean(mov_mc_10x,3);
mov_mc_10x=vm(mov_mc_10x);
mov_mc_10x.transpose.savebin([fpth10x '/mc_ShutterReg' '.bin'])
mcTrace=xyField; % Normcorre
save([fpth10x '/mcTrace.mat'],'mcTrace','avg_im')

%% Load MC movie

load([fpth25x '/output_data.mat']);
sz=double(Device_Data{1, 3}.ROI([2 4]));
mov_25x=double(readBinMov([fpth25x '/mc_ShutterReg.bin'],sz(2),sz(1)));

load([fpth10x '/output_data.mat']);
sz=double(Device_Data{1, 3}.ROI([2 4]));
mov_10x=double(readBinMov([fpth10x '/mc_ShutterReg.bin'],sz(2),sz(1)));

%%

[roi_25x, tr_25x]=clicky(mov_25x);
[roi_10x, tr_10x]=clicky(mov_10x);

% tr_25x=tr_25x-movmedian(tr_25x,50);
% tr_10x=tr_10x-movmedian(tr_10x,50);
%%
% load([fpth25x '/mcTrace.mat'])
% bkg = zeros(2, size(mov_25x,3));
% bkg(1,:) = linspace(-1, 1, size(mov_25x,3));  % linear term
% bkg(2,:) = linspace(-1, 1, size(mov_25x,3)).^2;  % quadratic term
% %bkg(3,:) = tr_25x(:,2);  % quadratic term
% 
% mc=mcTrace.xymean;
% mov_25x_res= mov_25x-median(mov_25x,3);
% mov_25x_res = SeeResiduals(mov_25x_res,mc);
% mov_25x_res = SeeResiduals(mov_25x_res,mc.^2);
% mov_25x_res = SeeResiduals(mov_25x_res,mc(:,1).*mc(:,end));
% mov_25x_res= SeeResiduals(mov_25x_res,bkg,1);
% 
% tr_25x=apply_clicky(roi_25x,mov_25x_res);
% 
% load([fpth10x '/mcTrace.mat'])
% bkg = zeros(2, size(mov_10x,3));
% bkg(1,:) = linspace(-1, 1, size(mov_10x,3));  % linear term
% bkg(2,:) = linspace(-1, 1, size(mov_10x,3)).^2;  % quadratic term
% %bkg(3,:) = mov_10x(:,2);  % quadratic term
% 
% mc=mcTrace.xymean;
% mov_10x_res= mov_10x-median(mov_10x,3);
% mov_10x_res = SeeResiduals(mov_10x_res,mc);
% mov_10x_res = SeeResiduals(mov_10x_res,mc.^2);
% mov_10x_res = SeeResiduals(mov_10x_res,mc(:,1).*mc(:,end));
% mov_10x_res= SeeResiduals(mov_10x_res,bkg,1);
% 
% tr_10x=apply_clicky(roi_10x,mov_10x_res);
%%
tr_25x_norm=-(tr_25x-movmedian(tr_25x,30))';
tr_25x_norm=tr_25x_norm./get_threshold(tr_25x_norm,1);
sp_25x=find_spike_bh(tr_25x_norm,5,3);

tr_10x_norm=-(tr_10x-movmedian(tr_10x,30))';
tr_10x_norm=tr_10x_norm./get_threshold(tr_10x_norm,1);
sp_10x=find_spike_bh(tr_10x_norm,5,3);

figure(2); clf; nTau=[-30:10];
tiledlayout(2,6)
nexttile([1 5])
plot(-tr_25x(1:end-1)); hold all
plot(find(sp_25x),-tr_25x(find(sp_25x)),'ro')
nexttile([1 1])
sp_25x([1:-nTau(1) end-nTau(2):end])=0;
STA_25x=-tr_25x(find(sp_25x)'+nTau);
errorbar_shade(nTau,mean(STA_25x,1),std(STA_25x,0,1))
title(['\DeltaF/F = ' num2str(-(max(mean(STA_25x,1))-prctile(mean(STA_25x,1),20))/prctile(mean(STA_25x,1),20)*100,2) '%'])
nexttile([1 5])
plot(-tr_10x(1:end-1)); hold all
plot(find(sp_10x),-tr_10x(find(sp_10x)),'ro')
nexttile([1 1])
sp_10x([1:-nTau(1) end-nTau(2):end])=0;
STA_10x=-tr_10x(find(sp_10x)'+nTau);
errorbar_shade(nTau,mean(STA_10x,1),std(STA_10x,0,1))
title(['\DeltaF/F = ' num2str(-(max(mean(STA_10x,1))-prctile(mean(STA_10x,1),20))/prctile(mean(STA_10x,1),20)*100,2) '%'])
%%
STA_mov_25x=squeeze(mean(reshape(mov_25x(:,:,find(sp_25x)'+nTau),size(mov_25x,1),size(mov_25x,2),[],length(nTau)),3));
STA_mov_25x=STA_mov_25x-median(STA_mov_25x,3);
STA_mov_10x=squeeze(mean(reshape(mov_10x(:,:,find(sp_10x)'+nTau),size(mov_10x,1),size(mov_10x,2),[],length(nTau)),3));
STA_mov_10x=STA_mov_10x-median(STA_mov_10x,3);

%%

writeMov([fpth25x '\STAmov'],-STA_mov_25x,nTau,10,1,[-50 200])
writeMov([fpth10x '\STAmov'],-STA_mov_10x,nTau,10,1,[-7 15])