function allWaves = dragonflyMultiPatWaveform_2color_yq(df_app)%nCells,recordFrames,maxBlueIntensity,centeredROI,exposureTime)
%
%   2018-2019 Vicente Parot, 2021 Yitong Qi
%   Cohen Lab - Harvard university
%
%---------- Usage Notes: -----------------
% 1. load dmd sequence of length, n_pat
% 2. specify waveform for each pattern
%-----------------------------------------

% Hardware notes:
% Camera (Teledyne Kinetix) in sychronous readout trigger mode:
%   trigger to exposure delay: 3.75 us
%   trigger to exposure delay jitter: 3.75 us

% DMD max switch rate: 22727 Hz
%% parameter definition section

% Collect app params
dt                  =   1/df_app.ni6343.rateConfigData; % DAQ clock rate
t_pulse_duration    =   dt;  % use dt max precision
t_cam_frame         =  df_app.ExposureTimeMsEditField_waveform.Value*1e-3; % s
t_gs_frame          =  df_app.GSExposureTimemsEditField.Value*1e-3; % s
n_frames            = df_app.camFramesInEditField.Value;
n_gs_frames         = floor(n_frames*t_cam_frame/t_gs_frame);
t_total             = (n_frames)*t_cam_frame;   % extra trigger to initiate recording, 
                                                  % one extra frame time
                                                  
n_pat               = size(df_app.dmdApp.mask_sequence,3); % # of dmd patterns; 
orange_duty_cycle   = 0.8;
blue_duty_cycle     = 1 - orange_duty_cycle;
t_dmd_frame_orange  = floor((t_cam_frame)*orange_duty_cycle/dt)*dt-dt; % expo time orange
t_dmd_frame_blue    = floor((t_cam_frame)*blue_duty_cycle/(n_pat-1)/dt)*dt - dt; % expo time blue

daq_channels        = df_app.ni6343.outChannelsConfigData(:,1);

v_blue              = df_app.BlueAmpVEditField.Value;
v_ir                = df_app.IRAmpVEditField.Value;
v_orange            = df_app.OrangeAmpVEditField.Value;

% Waveform definition start

wp0 = baseWave(dt,t_pulse_duration); % 0-valued constant wave of pulse duration
wp1 = wp0+1; % 1-valued constant wave of pulse duration

w_stim_frame_orange = [wp0 wp0.extend(t_dmd_frame_orange-dt*2)+1 wp0];
w_stim_frame_blue = [extend(wp0,dt*2) wp0.extend(t_dmd_frame_blue-dt*2)+1 wp0];

w_cam_frame = wp1.extend(t_cam_frame);
w_gs_frame = wp1.extend(t_gs_frame);
w_dmd_frame = extend([wp0.extend(dt) wp1.extend(t_dmd_frame_orange) ...
                    repeat(wp1.extend(t_dmd_frame_blue),n_pat-1)],...
                    t_cam_frame);


w_pat = [wp0, wp0.extend(dt),... 
        repeat(extend(...
        [wp0.extend(dt) w_stim_frame_orange],...
        t_cam_frame),n_frames)];
for ii=2:n_pat
    w_pat = [w_pat;
        wp0, wp0.extend(dt),... 
        repeat(extend(...
        [wp0.extend(dt) wp0.extend(t_dmd_frame_orange) ...
        wp0.extend(t_dmd_frame_blue*(ii-2)) w_stim_frame_blue],...
        t_cam_frame),n_frames)
        ];
end
w_cam = extend([wp0 w_cam_frame.repeat(n_frames)],t_total);
w_gs = extend([wp0 w_gs_frame.repeat(n_gs_frames)],t_total);  % extra trigger to initiate recording
w_dmd = extend([wp0 wp0.extend(dt) w_dmd_frame.repeat(n_frames)],t_total);
w_ir = (wp0.extend(t_total)+1) * v_ir;
%% waveform definition section
w_orange = flexRampTrainWave(dt,t_total);
w_orange.tBefore = 0;
w_orange.tRamp = t_total;
w_orange.period = t_total;

% intens = [v_orange;v_blue;max(v_blue*0.2,0.2)];
intens = [v_orange;v_blue];


%----1 constant stim sequence------
w_stim1 = flexRampTrainWave(dt,t_total);
w_stim1.tBefore = 0.495; % duration off 
w_stim1.tRamp = 100e-3; % duration on
w_stim1.period = 1; % seconds

w_stim = [w_orange; w_stim1]*intens;
%------------end---------------------

%----2 constant stim sequence------
% w_stim1 = flexRampTrainWave(dt,t_total);
% w_stim1.tBefore = 1.5; % duration off 
% w_stim1.tRamp = .5; % duration on
% w_stim1.period = t_total; % seconds
% 
% w_stim2 = flexRampTrainWave(dt,t_total);
% w_stim2.tBefore = 2; % duration off 
% w_stim2.tRamp = 1.5; % duration on
% w_stim2.period = t_total; % seconds
% % w_stim2 = [baseWave(dt,0.35) w_stim2 baseWave(dt,0.25)];
% 
% w_stim = [w_orange; w_stim1; w_stim2]*intens;
%------------end---------------------

%----sequential step stim sequence------
% w_stim1 = flexRampTrainWave(dt,t_total);
% w_stim1.tBefore = 1.5; % duration off 
% w_stim1.tRamp = 3; % duration on
% w_stim1.period = 40; % seconds
% 
% % w_stim2 = flexRampTrainWave(dt,t_total);
% % w_stim2.tBefore = .5; % duration off 
% % w_stim2.tRamp = 2; % duration on
% % w_stim2.period = t_total; % seconds
% 
% 
% w_stim = [w_orange; w_stim1]*intens;
%------------end---------------------

%---- 1 ramp + 2 const sequence------
% w_const = flexRampTrainWave(dt,t_total);
% w_const.tBefore = 0.5; % duration off
% w_const.tRamp = 0.5; % duration on
% w_const.period = t_total; % seconds
% 
% w_ramp = rampWave(dt,t_total);
% w_ramp.tBefore = 0.6;
% w_ramp.tRamp = 0.6;
% 
% w_stim = [w_orange;w_const;w_ramp]*intens([1 2 2]);
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
% w_stim = [w_orange; w_ramp_up;w_ramp_down;w_const]*intens([1 2 2 2]);
%------------end---------------------

%----optopatch stim sequence---
% w_stim = baseWave(dt,.5);
% for amp_factor = [.1 .2 .5 1]
% w_stim = [w_stim extend(baseWave(dt,.5)+1*amp_factor,1)];
% end
% w_ramp = rampWave(dt,5); w_ramp.tRamp = 5;
% w_stim = [w_stim w_ramp];
% w_stim = w_stim.extend(t_total);
% w_stim = [w_orange;w_stim].*intens;
%---------- end --------------

%----n ramp sequences---------
% n = 1;
% w_ramp_all = baseWave(dt,t_total);
% t_ramp = 5;
% t_wait = 0.25;
% for tbefore = t_wait:t_ramp:((n-1)*t_ramp+t_wait)
% w_ramp_temp = rampWave(dt,t_total); 
% w_ramp_temp.tBefore = tbefore;
% w_ramp_temp.tRamp = t_ramp;
% 
% w_ramp_all = w_ramp_all+w_ramp_temp;
% end
% 
% w_stim = [w_orange; w_ramp_all]*intens;
%---------- end --------------

%----2 ramp sequences---
% w_ramp1 = rampWave(dt,t_total); 
% w_ramp2 = rampWave(dt,t_total); 
% w_env = flexRampTrainWave(dt,t_total);
% w_env.tBefore = .1; % duration off 
% w_env.tRamp = .25; % duration on
% w_env.period = 0.5; % seconds
% 
% w_ramp1.tBefore = 0;
% w_ramp1.tRamp = 5;
% 
% w_ramp2.tBefore = 5;
% w_ramp2.tRamp = 5;
% 
% 
% w_stim = [w_orange; w_ramp1; w_ramp2]*intens*w_env;
%---------- end --------------

% ----3 ramp sequences---
% w_ramp1 = rampWave(dt,t_total); 
% w_ramp2 = rampWave(dt,t_total); 
% w_ramp3 = rampWave(dt,t_total); 
% w_env = flexRampTrainWave(dt,t_total);
% w_env.tBefore = .1; % duration off 
% w_env.tRamp = .25; % duration on
% w_env.period = 0.5; % seconds
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
% w_stim = [w_orange; w_ramp1; w_ramp2; w_ramp3]*intens*w_env;
%---------- end --------------


%----2 constant stim sequence------
% w_stim1 = flexRampTrainWave(dt,t_total);
% w_stim1.tBefore = 0; % duration off 
% w_stim1.tRamp = t_total; % duration on
% w_stim1.period = t_total; % seconds
% 
% w_stim2 = baseWave(dt,.5);
% for amp_factor = [.1 .2 .5 1]
% w_stim2 = [w_stim2 extend(baseWave(dt,.5)+1*amp_factor,1)];
% end
% w_ramp = rampWave(dt,5); w_ramp.tRamp = 5;
% w_stim2 = [w_stim2 w_ramp];
% w_stim2 = w_stim2.extend(t_total);
% 
% w_stim = [w_orange; w_stim1; w_stim2]*intens;

%------------end---------------------


w_stim_all = extend(w_pat,t_total).*w_stim;
w_blue = baseWave(dt,t_total);
w_orange = baseWave(dt,t_total);

w_orange = w_orange+w_stim_all.amplitude(1,:);
for ii=2:n_pat
%     w_blue = w_blue+AOTF_LUT(w_stim_all.amplitude(ii,:));
    w_blue = w_blue+w_stim_all.amplitude(ii,:);
end
% w_blue_correction = zeros(size(w_blue.amplitude));
% w_blue_correction(w_blue.amplitude<.2 & w_blue.amplitude > 0) = ...
%     .2 - w_blue.amplitude(w_blue.amplitude<.2 & w_blue.amplitude > 0);
% w_blue = w_blue+w_blue_correction;
%% DAQ waveform definition

w_out{1,1} = 'cam';
w_out{1,2} = w_cam;

w_out{2,1} = 'blueIntensity';
w_out{2,2} = w_blue;

w_out{4,1} = 'orangeIntensity';
w_out{4,2} = w_orange;

w_out{3,1} = 'dmd';
w_out{3,2} = w_dmd;

w_out{5,1} = 'grashopper';
w_out{5,2} = w_gs;

w_out{6,1} = 'IR';
w_out{6,2} = w_ir;
% 
% w_out{4,1} = 'blueShutter';
% w_out{4,2} = wp0.extend(t_total)+1;
% 
% w_out{5,1} = 'pololu';
% w_out{5,2} = wp0.extend(t_total)+0;

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

% figure
% plot(allWaves+[0:size(allWaves.amplitude,1)-1]'*2)
% legend(daq_channels)
end