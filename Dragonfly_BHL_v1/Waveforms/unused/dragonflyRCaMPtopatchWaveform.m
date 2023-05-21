function allWaves = dragonflyRCaMPtopatchWaveform(nCells,recordFrames,maxBlueIntensity,centeredROI,exposureTime)
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
    fprintf('using max blue intensity %g mW/cm2, amp %g V\n',maxBlueIntensity,maxBlueAmp)
    
    dt                  =   1e-5; % s
    t_pulse_duration    =   4e-4; % s % use dt max precision
    t_activity_period   = 180e-3; % s
    t_cam_frame         =  15e-3; % s
    t_pat_exposure      = 7.8e-3; % s % use dt max precision
    t_total             = 57;   % s
    n_periods_sequence = 12;
    n_frames_period = t_activity_period/t_cam_frame;
    n_frames_stim = 28;
    n_periods_wait = 30; 
    blue_ints = blueLUTmW(2.^(-2:3)/8,maxBlueIntensity);

    assert(~mod(t_pulse_duration,dt),'pulse duration must be multiple of dt')
    assert(~mod(t_pat_exposure,dt),'pattern exposure must be multiple of dt')
    assert(~mod(t_cam_frame,dt),'camera frame must be multiple of dt')
    assert(~mod(t_activity_period,t_cam_frame),'activity period must be multiple of camera frame')

    pre_offset = .28;

    %
    wp0 = baseWave(dt,t_pulse_duration); % 0-valued constant wave of pulse duration
    wp1 = wp0+1; % 1-valued constant wave of pulse duration

    w_cam_frame = wp1.extend(t_cam_frame);
    w_cam_sequence = w_cam_frame.repeat(n_frames_period*n_periods_sequence+4);
    w_cam_epoch =  [w_cam_sequence     w_cam_frame.repeat(n_frames_stim) w_cam_sequence];

    w_dmd_frame = extend([wp0.extend(4.9e-3) wp1.extend(t_pat_exposure) wp1],t_cam_frame);
    w_dmd_sequence = w_dmd_frame.repeat(n_frames_period*n_periods_sequence+4);
    w_dmd_epoch =  [w_dmd_sequence   extend([wp0.extend(t_cam_frame) wp1.extend(8e-4).repeat(501)],t_cam_frame*n_frames_stim) w_dmd_sequence  ];

    w_blue_epoch = [w_dmd_sequence*0 extend([wp0.extend(t_cam_frame) wp0.extend(.4)+1 wp0],   t_cam_frame*n_frames_stim) w_dmd_sequence*0];
    w_lime_epoch = [w_cam_sequence*0+1 w_cam_frame.repeat(n_frames_stim)*0 w_cam_sequence*0+1];

    allWaves = extend([wp0;wp0;wp0],pre_offset); % offset to compare with legacy waves
    w_lime =   extend( wp0,         pre_offset); % offset to compare with legacy waves
%     w_cam_frame0 = wp0.extend(t_cam_frame);
%     allWaves = [allWaves [w_cam_frame;w_cam_frame0;w_cam_frame0] [w_cam_frame0;w_cam_frame0;w_cam_frame0]]; %#ok to grow
    for it = 1:numel(blue_ints)
        allWaves = [allWaves [
            w_cam_epoch                 
            w_dmd_epoch                 
            w_blue_epoch.*blue_ints(it) 
            ]]; %#ok to grow
        w_lime = [w_lime [
            w_lime_epoch                
            ]]; %#ok to grow
        if it < numel(blue_ints)
            allWaves = allWaves.extend(allWaves.duration + t_activity_period*n_periods_wait + w_cam_frame.duration*4);
            w_lime   =   w_lime.extend(  w_lime.duration + t_activity_period*n_periods_wait + w_cam_frame.duration*4);
        else
            % add extra pulse to camera waveform, to close last frame
            allWaves = [allWaves [[wp1 wp0]; repmat([wp0 wp0],allWaves.numChannels-1)]]; %#ok<AGROW>
            w_lime = [w_lime [wp0 wp0]]; %#ok<AGROW>
        end
    end
    w_lime = w_lime.advanceRisingEdges(10e-3);
    allWaves = [allWaves; w_lime*0; w_lime]; % append red and lime lines
    allWaves =  allWaves.extend(t_total);
    
    %
    % f = baseWave(dt,nCells+6); % a constant wave, zero valued by default
    % do_redw = [z+1 f+1 z];

    % w_cam = w_cam.extend(told(numel(told)));
    % w_dmd = w_dmd.extend(w_cam.duration);
    % allWaves = [
    %     w_cam - .2 % - wcam'
    %     w_dmd - .3 % - wdmd' 
    %     ];
%     return
%     save(fname,'allWaves')
    % %% test
    srcdir = 'X:\Lab\Labmembers\Vicente Parot\Data\2018-03-21 CaP87 PFC with tracers\Hadamard runtime\wf had conf waves\10s wait 945 cycles backup\';
    w008 = dlmread(fullfile(srcdir,'8'));
    w010 = dlmread(fullfile(srcdir,'10'));
    wcam = dlmread(fullfile(srcdir,'21'));
    wdmd = dlmread(fullfile(srcdir,'22'));
    told = (1:numel(wcam))*1e-4;
    w008(told>t_total) = [];
    w010(told>t_total) = [];
    wcam(told>t_total) = [];
    wdmd(told>t_total) = [];
    told(told>t_total) = [];
    % figure windowstyle docked
    % plot(told,[wcam wdmd-.1])
    % xlim([.2 .4])
    
    figure windowstyle docked
%     clf
    plot(told,[wcam [wdmd w008 w010]-[1:3]/10])
    hold on
    plot(allWaves - ((1:5)+3)'/10)
%     % xlim([.12 .4])
%     % xlim([0 5])
% %     xlim([2.3 2.8])
% 
%     legend 'old camera trig' 'old DMD trig' old old 'camera trig' 'DMD trig' 'blue intensity' 'red shutter' 'lime shutter'
%     f = gcf;
%     f.PaperPosition = f.PaperPosition.*[0 0 1 1];
%     f.PaperSize = f.PaperPosition(3:4);
%     saveas(gcf,fullfile('Waveforms','allRcamptopatchWaves.pdf'))

end

