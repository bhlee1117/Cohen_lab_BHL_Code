function allWaves = dragonflyPiezoTestWaveform_yq(df_app)%nCells,recordFrames,maxBlueIntensity,centeredROI,exposureTime)
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

daq_channels        = df_app.ni6343.outChannelsConfigData(:,1);
% Waveform definition start


wp0 = baseWave(dt,t_pulse_duration); % 0-valued constant wave of pulse duration
wp1 = wp0+1; % 1-valued constant wave of pulse duration

w_cam_frame = wp1.extend(t_cam_frame);

w_cam = extend([wp0 w_cam_frame.repeat(n_frames+1)],t_total);  % extra trigger to initiate recording


%----sine stim sequence------
w_stim = flexSineWave(dt,t_total);
w_stim.tBefore = 4; % duration off 
w_stim.tOn = 4; % duration on
w_stim.period = 1/10; 

% w_env = flexRampTrainWave(dt,t_total);
% w_env.tBefore = 0; % duration off 
% w_env.tRamp = .1; % duration on
% w_env.period = 1/2; % seconds
% 
% w_stim = w_env*w_stim;
%------------end---------------------

%----constant stim sequence------
% w_stim = flexRampTrainWave(dt,t_total);
% w_stim.tBefore = 0; % duration off 
% w_stim.tRamp = .05; % duration on
% w_stim.period = 1/3; % seconds
%------------end---------------------


%%

w_out{1,1} = 'cam';
w_out{1,2} = w_cam;

w_out{2,1} = 'piezo';
w_out{2,2} = w_stim*9;

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
line(allWaves.time,allWaves.amplitude+[0:size(allWaves.amplitude,1)-1]'*2)
legend cam dmd blue red lime orange
end