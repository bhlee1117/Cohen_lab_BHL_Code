function allWaves = dragonflyMultiPatWaveform_sw_yq(df_app)%nCells,recordFrames,maxBlueIntensity,centeredROI,exposureTime)
%
%   2018-2019 Vicente Parot, 2021 Yitong Qi
%   Cohen Lab - Harvard university
%
%---------- Usage Notes: -----------------
% 1. load dmd sequence of length, n_pat
% 2. specify waveform for each pattern
%-----------------------------------------

% Hardware notes:
% Camera (Hamamatsu Flash) in sychronous readout trigger mode:
%   trigger to exposure delay: 165.58 us
%   trigger to exposure delay jitter: 9.74 us

% DMD max switch rate: 22727 Hz
%% parameter definition section

t_jit = 1e-5; %s trigger to exposure delay jitter
t_exp_delay = 16e-5; %s trigger to exposure delay
% Collect app params
dt                  =   1/df_app.ni6343.rateConfigData; % DAQ clock rate
t_pulse_duration    =   dt;  % use dt max precision
t_cam_frame         =  df_app.ExposureTimeMsEditField_waveform.Value*1e-3; % s
n_frames            = df_app.camFramesInEditField.Value;
t_total             = (n_frames+2)*t_cam_frame;   % extra trigger to initiate recording, 
                                                  % one extra frame time
                                                  
n_pat               = size(df_app.dmdApp.mask_sequence,3); % # of dmd patterns; 
orange_duty_cycle   = 0.8;
blue_duty_cycle     = 1 - orange_duty_cycle;
t_dmd_frame_orange  = floor((t_cam_frame-t_jit)*orange_duty_cycle/dt)*dt; % expo time orange
t_dmd_frame_blue    = floor((t_cam_frame-t_jit)*blue_duty_cycle/(n_pat-1)/dt)*dt; % expo time blue

daq_channels        = df_app.ni6343.outChannelsConfigData(:,1);

v_blue              = df_app.BlueAmpVEditField.Value;
v_orange            = df_app.OrangeAmpVEditField.Value;

% Waveform definition start

wp0 = baseWave(dt,t_pulse_duration); % 0-valued constant wave of pulse duration
wp1 = wp0+1; % 1-valued constant wave of pulse duration

w_stim_frame_orange = [wp0 wp0.extend(t_dmd_frame_orange-dt)+1];
w_stim_frame_blue = [wp0 wp0.extend(t_dmd_frame_blue-dt)+1];

w_cam_frame = wp1.extend(t_cam_frame);
w_dmd_frame = extend([wp0.extend(t_jit) wp1.extend(t_dmd_frame_orange) ...
                    repeat(wp1.extend(t_dmd_frame_blue),n_pat-1)],...
                    t_cam_frame);


w_pat = [wp0, wp0.extend(t_exp_delay),... 
        repeat(extend(...
        [wp0.extend(t_jit) w_stim_frame_orange],...
        t_cam_frame),n_frames)];
for ii=2:n_pat
    w_pat = [w_pat;
        wp0, wp0.extend(t_exp_delay),... 
        repeat(extend(...
        [wp0.extend(t_jit) wp0.extend(t_dmd_frame_orange) ...
        wp0.extend(t_dmd_frame_blue*(ii-2)) w_stim_frame_blue],...
        t_cam_frame),n_frames)
        ];
end
w_cam = extend([wp0 w_cam_frame.repeat(n_frames+1)],t_total);
w_dmd = extend([wp0 wp0.extend(t_exp_delay) w_dmd_frame.repeat(n_frames+1)],t_total);

%% waveform definition section
w_orange = flexRampTrainWave(dt,t_total);
w_orange.tBefore = 0;
w_orange.tRamp = t_total;
w_orange.period = t_total;

intens = [v_orange;1;1];
t_on_1 = 0.05;
t_on = 1;
%----2 constant stim sequence------
% w_stim1 = flexRampTrainWave(dt,t_total);
% w_stim1.tBefore = 0.5; % duration off 
% w_stim1.tRamp = .1; % duration on
% w_stim1.period = t_total; % seconds
% 
% w_stim2 = flexRampTrainWave(dt,t_total);
% w_stim2.tBefore = 0.6; % duration off 
% w_stim2.tRamp = 0.4; % duration on
% w_stim2.period = t_total; % seconds
% % w_stim2 = [baseWave(dt,0.35) w_stim2 baseWave(dt,0.25)];
% 
% w_stim = [w_stim1; w_stim2]*intens;
%------------end---------------------

%----2 constant stim sequence------
w_stim1 = flexRampTrainWave(dt,t_total);
w_stim1.tBefore = 0.25; % duration off 
w_stim1.tRamp = t_on_1; % duration on
w_stim1.period = t_total; % seconds

w_stim2 = flexRampTrainWave(dt,t_total);
w_stim2.tBefore = .25+t_on_1; % duration off 
w_stim2.tRamp = t_on-t_on_1; % duration on
w_stim2.period = t_total; % seconds
% % w_stim2 = [baseWave(dt,0.35) w_stim2 baseWave(dt,0.25)];
% 
% w_stim3 = flexRampTrainWave(dt,t_total);
% w_stim3.tBefore = 0.55; % duration off 
% w_stim3.tRamp = 0.4; % duration on
% w_stim3.period = t_total; % seconds

w_stim = [w_orange; w_stim1; w_stim2;]*intens;
%------------end---------------------

%----sequential step stim sequence------
% w_stim1 = flexRampTrainWave(dt,t_total);
% w_stim1.tBefore = .5; % duration off 
% w_stim1.tRamp = 2; % duration on
% w_stim1.period = t_total; % seconds
% 
% w_stim2 = flexRampTrainWave(dt,t_total);
% w_stim2.tBefore = .5; % duration off 
% w_stim2.tRamp = 2; % duration on
% w_stim2.period = t_total; % seconds
% 
% 
% w_stim = [w_stim1; w_stim2]*intens;
%------------end---------------------

%---- 1 ramp + 2 const sequence------
% w_const = flexRampTrainWave(dt,t_total);
% w_const.tBefore = 0.5; % duration off
% w_const.tRamp = 0.5; % duration on
% w_const.period = t_total; % seconds
% 
% w_ramp = rampWave(dt,t_total);
% w_ramp.tBefore = 0.6;
% w_ramp.tRamp = 0.4;
% 
% w_stim = [w_const; w_ramp]*intens;
%------------end---------------------

%---- 1 ramp + 2 const + 3 ramp down sequence------
% w_const = flexRampTrainWave(dt,t_total);
% w_const.tBefore = 0; % duration off
% w_const.tRamp = 1; % duration on
% w_const.period = t_total; % seconds
% 
% w_ramp_up = rampWave(dt,t_total); 
% w_ramp_up.tRamp = t_total;
% 
% w_ramp_down = rampWave(dt,t_total); 
% w_ramp_down.tRamp = t_total;
% w_ramp_down = 1-w_ramp_down;
% 
% w_stim = [w_const; w_ramp_up;w_ramp_down]*intens;
%------------end---------------------

%----optopatch stim sequence---
% w_stim = baseWave(dt,.5);
% for amp_factor = [.1 .2 .5 1]
% w_stim = [w_stim extend(baseWave(dt,.5)+1*amp_factor,1)];
% end
% w_ramp = rampWave(dt,5); w_ramp.tRamp = 5;
% w_stim = [w_stim w_ramp];
% w_stim = w_stim.extend(t_total);
% w_dmd = extend([wp0],t_total);
%---------- end --------------

%----1 ramp sequences---
% w_ramp1 = rampWave(dt,t_total); 
% w_ramp1.tBefore = 0;
% w_ramp1.tRamp = 5;
% w_stim = w_ramp1*intens;
%---------- end --------------

%----2 ramp sequences---
% w_ramp1 = rampWave(dt,t_total); 
% w_ramp2 = rampWave(dt,t_total); 
% 
% w_ramp1.tBefore = 0;
% w_ramp1.tRamp = 5;
% 
% w_ramp2.tBefore = 5;
% w_ramp2.tRamp = 5;
% 
% w_stim = [w_ramp1; w_ramp2]*intens;
%---------- end --------------

%----3 ramp sequences---
% w_ramp1 = rampWave(dt,t_total); 
% w_ramp2 = rampWave(dt,t_total); 
% w_ramp3 = rampWave(dt,t_total); 
% 
% w_ramp1.tBefore = 0;
% w_ramp1.tRamp = 5;
% 
% w_ramp2.tBefore = 5;
% w_ramp2.tRamp = 5;
% 
% w_ramp3.tBefore = 10;
% w_ramp3.tRamp = 5;
% 
% w_stim = [w_ramp1; w_ramp2; w_ramp3]*intens;
%---------- end --------------

w_stim_all = extend(w_pat,t_total).*w_stim;
w_blue = baseWave(dt,t_total);
w_orange = baseWave(dt,t_total);

w_orange = w_orange+double(w_stim_all.amplitude(1,:));
for ii=2:n_pat
    w_blue = w_blue+double(w_stim_all.amplitude(ii,:));
end
%% DAQ waveform definition

w_out{1,1} = 'cam';
w_out{1,2} = w_cam;

w_out{2,1} = 'blue';
w_out{2,2} = w_blue;

w_out{3,1} = 'dmd';
w_out{3,2} = w_dmd;

w_out{4,1} = 'orangeIntensity';
w_out{4,2} = w_orange;
% w_out{4,1} = 'orangeShutter';
% w_out{4,2} = wp0.extend(t_total)+0;
% 
% w_out{5,1} = 'pololu';
% w_out{5,2} = wp0.extend(t_total)+1;


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
legend cam dmd blue red lime orange
end