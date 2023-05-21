function allWaves = dragonflyHadamardWaveform_yq(nCells,recordFrames,icolor,maxBlueIntensity,centeredROI,exposureTime)
%
%   2018-2020 Vicente Parot/Yitong Qi
%   Cohen Lab - Harvard university
%

    %%
    if ~exist('nCells','var')
        nCells = [];
    end
    if ~exist('recordFrames','var')
        recordFrames = 10;
    end
    if ~exist('maxBlueIntensity','var')
        maxBlueIntensity = .1;
    end
%     maxBlueAmp = blueLUT_mWcm2_2_V(maxBlueIntensity);
    if ~exist('centeredROI','var')
        centeredROI = [];
    end
    if ~exist('exposureTime','var')
        exposureTime = 30;
    end
    
    % rolling shutter, ~10*(n-1) us delay for the n-th horizontal line
    % wait till all lines are exposed, then turn on dmd
    % turn off dmd after set exposure time
    
    slines = 1024; %number of horizontal lines the camera needs to scan for each frame
    
    dt                  =   1e-5; % s
    t_pulse_duration    =   4e-4; % s % use dt max precision
    t_cam_frame         =  exposureTime*1e-3 + 10e-6*slines; % s
    t_pat_wait          = 10e-6*slines;
    t_pat_exposure      =  exposureTime*1e-3; % s % use dt max precision
    n_periods_sequence  = recordFrames;
%     n_extra_dmd_frames  = 10;
%     t_total = t_cam_frame*(n_periods_sequence + n_extra_dmd_frames + 1);

    assert(~mod(t_pulse_duration,dt),'pulse duration must be multiple of dt')
    assert(~mod(t_pat_exposure,dt),'pattern exposure must be multiple of dt')
    assert(~mod(t_cam_frame,dt),'camera frame must be multiple of dt')

    %
    wp0 = baseWave(dt,t_pulse_duration); % 0-valued constant wave of pulse duration
    wp1 = wp0+1; % 1-valued constant wave of pulse duration
    
    w_cam_frame = extend([wp0.extend(2e-3) wp1],t_cam_frame);
    w_cam_sequence = [w_cam_frame.repeat(n_periods_sequence+1) extend(wp0,t_cam_frame)];

    % extra dmd frames at beggining, set n-2 to pseudorandom
    w_dmd_frame = extend([wp0.extend(1e-3) wp1 wp0.extend(t_pat_wait+1e-3) wp1],t_cam_frame);
    w_dmd_sequence = w_dmd_frame.repeat(n_periods_sequence+2);

    w_blue_sequence = w_cam_sequence*0;
    w_red_sequence  = w_cam_sequence*0;
    w_lime_sequence = w_cam_sequence*0;
    switch icolor
        case 'Blue'
            w_blue_sequence = w_cam_sequence*0 + maxBlueIntensity;
        case 'Lime'
            w_lime_sequence = w_cam_sequence*0+1;
    end
    
    
    allWaves = [
            w_cam_sequence
            w_dmd_sequence
            w_blue_sequence
            w_red_sequence 
            w_lime_sequence 
            ];
    
    figure
    plot(allWaves - ((1:5)+3)'/10)
    legend 'camera trig' 'DMD trig' 'blue intensity' 'red shutter' 'lime shutter'
%     f = gcf;
%     f.PaperPosition = f.PaperPosition.*[0 0 1 1];
%     f.PaperSize = f.PaperPosition(3:4);
%     saveas(gcf,fullfile('Waveforms',[datestr(now,'HHMMSS') 'waves.pdf']))

end