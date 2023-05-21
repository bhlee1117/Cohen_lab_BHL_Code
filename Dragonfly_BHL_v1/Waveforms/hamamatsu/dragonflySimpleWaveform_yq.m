function allWaves = dragonflySimpleWaveform_yq(df_app)%nCells,recordFrames,maxBlueIntensity,centeredROI,exposureTime)
%
%   2018-2019 Vicente Parot, 2021 Yitong Qi
%   Cohen Lab - Harvard university
%

% Camera in sychronous readout trigger mode
%   trigger to exposure delay: 165.58 us
%   trigger to exposure delay jitter: 9.74 us

% Collect app params
dt                  =   1/df_app.ni6343.rateConfigData; % DAQ clock rate
t_pulse_duration    =   dt;  % use dt max precision
t_cam_frame         =  df_app.ExposureTimeMsEditField_waveform.Value*1e-3; % s
n_frames            = df_app.camFramesInEditField.Value;
t_total             = (n_frames+2)*t_cam_frame;   % extra trigger to initiate recording, 
                                                  % one extra frame time
v_blue              = df_app.BlueAmpVEditField.Value;
% v_blue_2            = df_app.AOTFAmpVEditField.Value;
v_orange            = df_app.OrangeAmpVEditField.Value;

daq_channels        = df_app.ni6343.outChannelsConfigData(:,1);
% Waveform definition start


wp0 = baseWave(dt,t_pulse_duration); % 0-valued constant wave of pulse duration
wp1 = wp0+1; % 1-valued constant wave of pulse duration

w_cam_frame = wp1.extend(t_cam_frame);

w_cam = extend([wp0 w_cam_frame.repeat(n_frames+1)],t_total);  % extra trigger to initiate recording


%----constant stim sequence------
% w_stim = flexRampTrainWave(dt,t_total);
% w_stim.tBefore = 13; % duration off 
% w_stim.tRamp = 1; % duration on
% w_stim.period = t_total; % seconds
%------------end---------------------

%----optopatch stim sequence---
w_stim = baseWave(dt,.5);
for amp_factor = [.1 .2 .5 1]
w_stim = [w_stim extend(baseWave(dt,.5)+1*amp_factor,1)];
end
w_ramp = rampWave(dt,5); w_ramp.tRamp = 5;
w_stim = [w_stim w_ramp];
w_stim = w_stim.extend(t_total);
%---------- end --------------

%----optopatch stim sequence---

% w_ramp = rampWave(dt,5); w_ramp.tRamp = 5;
% w_stim = w_ramp;
% w_stim = w_stim.extend(t_total);
%---------- end --------------

%----sine stim sequence------
% w_pstim = flexSineWave(dt,t_total);
% w_pstim.tBefore = 1; % duration off 
% w_pstim.tOn = 2; % duration on
% w_pstim.period = 1/10; 
%------------end---------------------

%----stair sequence---
% w_stim = [];
% w_wait = baseWave(dt,.5); 
% for amp_factor = logspace(-1.92,0,10) %0.1-5 V
%     w_stim = [w_stim w_wait baseWave(dt,1.5)+1*amp_factor];
% end
% w_stim = w_stim.extend(t_total);
% w_dmd = extend([wp0],t_total);
%---------- end --------------

%----pulse sequence test---
% w_stim = [];
% w_wait = baseWave(dt,.5); 
% 
% for pulse_width = (0.1:0.1:1)*1e-3 %0.1-5 V
%     w_stim_p = flexRampTrainWave(dt,0.5);
%     w_stim_p.tBefore = 0; % duration off 
%     w_stim_p.tRamp = pulse_width-5e-5; % duration on
%     w_stim_p.period = pulse_width; % seconds
%     w_stim = [w_stim w_wait w_stim_p];
%     pulse_width;
% end
% w_stim = w_stim.extend(t_total);
% % w_dmd = extend([wp0],t_total);
%---------- end --------------


w_blue = w_stim * v_blue;
w_orange = (wp0.extend(t_total)+1) * v_orange;
% w_orange = w_stim * v_orange;
%%

w_out{1,1} = 'cam';
w_out{1,2} = w_cam;

w_out{3,1} = 'blueIntensity';
w_out{3,2} = w_blue;

w_out{2,1} = 'orangeIntensity';
w_out{2,2} = w_orange;

% w_out{3,1} = 'AOTFVoltage';
% w_out{3,2} = w_stim * 3;

% w_out{4,1} = 'piezoTrigger';
% w_out{4,2} = w_pstim*10;
% 
% w_out{4,1} = 'pololu';
% w_out{4,2} = wp0.extend(t_total)+0;

% w_out{5,1} = 'dmd';
% w_out{5,2} = extend([wp0 repeat(wp1.extend(.5),7)],t_total);
%%
idx_daq_ch = [];
for ii=1:size(w_out,1)
    idx_daq_ch = [idx_daq_ch find(cellfun(@(x)...
                    contains(x,w_out{ii,1},'ignorecase',1),daq_channels))];

end

allWaves = [];
for ii=1:length(daq_channels)
    is_used = find(idx_daq_ch==ii);
    if is_used
        allWaves = [allWaves; w_out{is_used,2}];
    else
        allWaves = [allWaves; wp0.extend(t_total)];
    end
    
end

figure
plot(allWaves+[0:size(allWaves.amplitude,1)-1]'*2)
legend(daq_channels)
end