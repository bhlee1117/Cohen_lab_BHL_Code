function allWaves = dragonflySegmentedWideFieldWaveform_yq(recordFrames,icolor,maxBlueIntensity,exposureTime)
    
    global seg_wide_field_params
    
    opts.Interpreter = 'tex';
    cell_input = inputdlg('\fontsize{10} Please specify [seg_size active_area]','',1,{''},opts);
    seg_wide_field_params = str2num(cell2mat(cell_input));
    
    seg_size = seg_wide_field_params(1);
    FOV_size = seg_wide_field_params(2:3);
    id11 = [1:seg_size:FOV_size(2)];
    id12 = id11+seg_size-1; id12(id12>FOV_size(2))=FOV_size(2);
    id12 = unique(id12);id11 = id11(1:length(id12));
    id21 = [1:seg_size:FOV_size(1)];
    id22 = id21+seg_size-1; id22(id22>FOV_size(1))=FOV_size(1);
    id22 = unique(id22);id21 = id21(1:length(id22));
    
    n_seg = length(id11);
            
    dt                  =   1e-5; % s
    t_pulse_duration    =   4e-4; % s % use dt max precision
    t_cam_frame         =  exposureTime*1e-3; % s
    n_frames            = recordFrames;
    t_total             = n_frames*t_cam_frame*n_seg+.1;   % s

    wp0 = baseWave(dt,t_pulse_duration); % 0-valued constant wave of pulse duration
    wp1 = wp0+1; % 1-valued constant wave of pulse duration

    w_cam_frame = wp1.extend(t_cam_frame);
    w_cam = w_cam_frame.repeat(n_frames*n_seg);
    w_cam = extend([wp0 w_cam wp1],t_total);

    w_stim = flexRampTrainWave(dt,t_total);
    w_stim.tBefore = 0; % duration off
    w_stim.tRamp = 1; % duration on
    w_stim.period = 1; % seconds

    w_lime = flexRampTrainWave(dt,t_total);
    w_lime.tBefore = 0.1; % duration off
    w_lime.tRamp = 0.2; % duration on
    w_lime.period = 0.5; % seconds

    w_dmd = repeat(extend([wp0 wp1 wp0],t_cam_frame*n_frames),n_seg);
    w_dmd = extend(w_dmd,t_total);
    
    %% 
    
    switch icolor
        case 'Blue'
            w_blue = w_stim;
            w_lime = w_stim.*0;
            w_red = w_cam*0+0;
        case 'Lime'
            w_blue = w_stim.*0;
            w_lime = w_stim.*1;
            w_red = w_cam*0;
    end

    allWaves = [
        w_cam                 
        w_dmd
        w_blue
        w_red
        w_lime
        ];
    figure
    line(dt:dt:allWaves(1).duration,allWaves.amplitude+[1:5]'*0.3)
    legend cam dmd blue red lime
end