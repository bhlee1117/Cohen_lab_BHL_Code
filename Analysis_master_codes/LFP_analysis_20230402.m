%% Data and Parameters
% Number of rows and columns for the image
close all; clear all;
cd /Volumes/cohen_lab/Lab/Labmembers/Pojeong Park/Data/230328_LFP_invivo/154100invivo_LFP_Test1_gain_10_10kHz
cameradata=load(fullfile('output_data.mat'));
nrow=cameradata.Device_Data{1, 4}.ROI(4);ncol=cameradata.Device_Data{1, 4}.ROI(2);

% Frame rate of camera in seconds

dt=cameradata.Device_Data{1, 4}.exposuretime; 


%% Patch data loading_CC

amp_data=cameradata.Device_Data{1, 2}.buffered_tasks.channels(1, 1).data  ; 


dt_A = 1e-5;
tsteps = dt_A:dt_A:length(amp_data)*dt_A;

figure;plot(tsteps,amp_data*1e1);ylabel('Voltage (mV)');xlabel('Time (s)'); 
% gain 10 = 100 mV/mV (value = 1e1); gain 100 =  1000 mV/mV (value = 1)
%ylim([-80 -10]);
saveas(gcf,'MP_patchdata.fig');
saveas(gcf,'MP_patchdata.png');
save('amp_data.mat','amp_data');


%% CC filter

dt_A = 1e-5;
tsteps = dt_A:dt_A:length(amp_data)*dt_A;

figure;plot(tsteps,smooth(amp_data*1e1,100)-mean(amp_data(1:1000)*1e1));ylabel('Voltage (mV)');xlabel('Time (s)');
%ylim([-8000 1000]);
% xlim([0 38]);

saveas(gcf,'MP_fitlered.fig');
saveas(gcf,'MP_filtered.png');


%%

  fpass=[1 7500]; 
  %[b, a] = butter(10, fpass/((1/dt_A)/2), 'bandpass');
  %filt_trace=filtfilt(b, a, amp_data*1e1);
  %filt_trace=bandpass(amp_data*1e1,fpass,1/dt_A);
  filt_trace=amp_data;
  filt_theta=bandpass(filt_trace,[3 12],1/dt_A);
  filt_gamma=bandpass(filt_trace,[30 90],1/dt_A);
  filt_SWR=bandpass(filt_trace,[150 300],1/dt_A);



figure; 
tiledlayout(4,1)
ax1=[];
ax1=[ax1 nexttile([1 1])];
plot(tsteps,amp_data*1e1-mean(amp_data(1:1000)*1e1))
title('Raw trace')
ax1=[ax1 nexttile([1 1])];
plot(tsteps,filt_theta)
title('3-12 Hz filtered')
ax1=[ax1 nexttile([1 1])];
plot(tsteps,filt_gamma)
title('30-90 Hz filtered')
ax1=[ax1 nexttile([1 1])];
plot(tsteps,filt_SWR)
title('150-300 Hz filtered')
linkaxes(ax1,'x')
%     for n=1:size(traces,1)
%     Result{i}.traces_res_hi_filtered(n,:) = filtfilt(b, a, double(Result{1}.traces_res_hi(n,:)));
%     end



