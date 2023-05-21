function allWaves = dragonflySimpleWaveform(nCells,recordFrames,maxBlueIntensity,centeredROI,exposureTime)
%
%   2018-2019 Vicente Parot
%   Cohen Lab - Harvard university
%
if ~exist('nCells','var')
    nCells = [];
end
if ~exist('recordFrames','var')
    recordFrames = 10000;
end
if ~exist('maxBlueIntensity','var')
    maxBlueIntensity = 400;
end
maxBlueAmp = blueLUT_mWcm2_2_V(maxBlueIntensity);
if ~exist('centeredROI','var')
    centeredROI = [];
end
if ~exist('exposureTime','var')
    exposureTime = 1;
end
fprintf('using max blue intensity %g mW/cm2, amp %g V\n',maxBlueIntensity,maxBlueAmp)

dt                  =   1e-5; % s
t_pulse_duration    =   4e-4; % s % use dt max precision
t_cam_frame         =  exposureTime*1e-3; % s
n_frames            = recordFrames;
t_total             = n_frames*t_cam_frame+.1;   % s

wp0 = baseWave(dt,t_pulse_duration); % 0-valued constant wave of pulse duration
wp1 = wp0+1; % 1-valued constant wave of pulse duration

w_cam_frame = wp1.extend(t_cam_frame);
w_cam = w_cam_frame.repeat(n_frames);
w_cam = extend([wp0 w_cam wp1],t_total);

w_stim = flexRampTrainWave(dt,t_total);
w_stim.tBefore = 0; % duration off
w_stim.tRamp = 1; % duration on
w_stim.period = 1; % seconds

w_lime = flexRampTrainWave(dt,t_total);
w_lime.tBefore = 0.1; % duration off
w_lime.tRamp = 0.2; % duration on
w_lime.period = 0.5; % seconds

w_dmd = extend([wp0 wp1],t_total);
w_blue = w_stim.*maxBlueAmp;
w_red = w_cam*0+1;
% w_lime = w_cam*0;
allWaves = [
    w_cam                 
    w_dmd
    w_blue
    w_red
    w_lime
    ];
figure
plot(allWaves+[1:5]'*0.3)
legend cam dmd blue red lime
end

