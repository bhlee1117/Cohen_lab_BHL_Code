function allWaves = dragonflyHadamardWaveform(nCells,recordFrames,maxBlueIntensity,centeredROI,exposureTime)
%
%   2018-2019 Vicente Parot
%   Cohen Lab - Harvard university
%
fname = fullfile('Waveforms','allRcamptopatchWaves.mat');
if false && exist(fname,'file')
    load(fname)
%     figure windowstyle docked
%     plot(allWaves - [.2 .3 .4 .5 .6]')
%     legend 'camera trig' 'DMD trig' 'blue intensity' 'lime shutter' 'red shutter'
else
    %%
    if ~exist('nCells','var')
        nCells = [];
    end
    if ~exist('recordFrames','var')
        recordFrames = [];
    end
    if ~exist('maxBlueIntensity','var')
        maxBlueIntensity = 400;
    end
    maxBlueAmp = blueLUT_mWcm2_2_V(maxBlueIntensity);
    if ~exist('centeredROI','var')
        centeredROI = [];
    end
    if ~exist('exposureTime','var')
        exposureTime = [];
    end
    
    dt                  =   1e-5; % s
    t_pulse_duration    =   4e-4; % s % use dt max precision
    t_cam_frame         =  24e-3; % s
    t_pat_exposure      =  12e-3; % s % use dt max precision
    n_periods_sequence  = 60;
    n_extra_dmd_frames  = 10;
    t_total = t_cam_frame*(n_periods_sequence + n_extra_dmd_frames + 3);

    assert(~mod(t_pulse_duration,dt),'pulse duration must be multiple of dt')
    assert(~mod(t_pat_exposure,dt),'pattern exposure must be multiple of dt')
    assert(~mod(t_cam_frame,dt),'camera frame must be multiple of dt')

    %
    wp0 = baseWave(dt,t_pulse_duration); % 0-valued constant wave of pulse duration
    wp1 = wp0+1; % 1-valued constant wave of pulse duration

    % 4 excess camera frames, 2 at beggining 2 at end
    w_cam_frame = wp1.extend(t_cam_frame);
    w_cam_sequence = [w_cam_frame.repeat(n_extra_dmd_frames-1)*0 w_cam_frame.repeat(n_periods_sequence+4)];

    % extra dmd frames at beggining, set n-2 to pseudorandom
    w_dmd_frame = extend([wp0.extend(t_cam_frame-t_pat_exposure-1e-3) wp1.extend(t_pat_exposure) wp1],t_cam_frame);
    w_dmd_sequence = [w_dmd_frame.repeat(1)*0 w_dmd_frame.repeat(n_extra_dmd_frames) w_dmd_frame.repeat(n_periods_sequence) w_dmd_frame.repeat(2)*0];

    w_blue_sequence = w_cam_sequence*0 + maxBlueAmp;
    w_red_sequence  = w_cam_sequence*0;
    w_lime_sequence = w_cam_sequence*0;
    
    allWaves = [
            w_cam_sequence  wp1 wp0          
            w_dmd_sequence  wp0 wp0           
            w_blue_sequence wp0 wp0
            w_red_sequence  wp0 wp0
            w_lime_sequence wp0 wp0
            ];
    
%     figure windowstyle docked
%     plot(allWaves - ((1:5)+3)'/10)
%     legend 'camera trig' 'DMD trig' 'blue intensity' 'red shutter' 'lime shutter'
%     f = gcf;
%     f.PaperPosition = f.PaperPosition.*[0 0 1 1];
%     f.PaperSize = f.PaperPosition(3:4);
%     saveas(gcf,fullfile('Waveforms',[datestr(now,'HHMMSS') 'waves.pdf']))

end

