function allWaves = dragonflyL1PatternCalWaveform_yq(recordFrames,icolor,maxBlueIntensity,exposureTime)
    global dmd_frames_1 dmd_frames_2;
    

    fig = uifigure;
    wave_type = uiconfirm(fig,'Please select L1 stimulation pattern sequence','','option',...
                {'Step',...
                'Ramp'},...
                'closefcn',@(h,e) close(fig));
    opts.Interpreter = 'tex';
    cell_input = inputdlg('\fontsize{10} Please specify pattern time sequence (fmt: [duration1 duration2 ...])','',1,{''},opts);
    input = str2num(cell2mat(cell_input));
    
    opts.Interpreter = 'tex';
    cell_input = inputdlg('\fontsize{10} Blue voltage amplitude sequence (fmt: [...])','',1,{''},opts);
    
    maxBlueAmp = str2num(cell2mat(cell_input));

    fprintf('using max blue amp %g V\n',maxBlueAmp)
    
    opts.Interpreter = 'tex';
    cell_input = inputdlg('\fontsize{10} ["DMD on/off-ratio" "DMD on-time (ms)"]','',1,{''},opts);
    
    dmd_timings = str2num(cell2mat(cell_input));

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
    
    on_off_ratio = dmd_timings(1);
    on_time = dmd_timings(2)*1e-3;
    w_dmd_frame = [extend([wp0.extend(10e-3) wp1 wp0.extend(on_time) wp1],t_cam_frame) ...
                   repeat(wp0.extend(t_cam_frame),on_off_ratio)];
    
    switch wave_type
        case 'Step'
            t_i = floor(length(dt:dt:input(1))/length(maxBlueAmp))*dt;
            w_stim = [];
            w_dmd = [];
            for i=1:length(maxBlueAmp)
                w_stim_i = flexRampTrainWave(dt,t_i);
                w_stim_i.tBefore = 0; % duration off
                w_stim_i.tRamp = 1; % duration on
                w_stim_i.period = 1; % seconds
                w_stim_i = w_stim_i*maxBlueAmp(i);
                w_stim = [w_stim w_stim_i];
                
                
            end
            w_stim_i = flexRampTrainWave(dt,input(2));
            w_stim_i.tBefore = 0; % duration off
            w_stim_i.tRamp = 1; % duration on
            w_stim_i.period = 1; % seconds
            w_stim_i = w_stim_i*maxBlueIntensity;
            w_stim = [w_stim w_stim_i];
            w_stim = extend(w_stim,t_total);

            dmd_frames_1 = floor(input(1)/t_cam_frame/(1+on_off_ratio));
            dmd_frames_2 = floor(input(2)/t_cam_frame/(1+on_off_ratio));
            w_dmd = [w_dmd w_dmd_frame.repeat(dmd_frames_1) w_dmd_frame.repeat(dmd_frames_2)];
            w_dmd = w_dmd.extend(t_total);
        case 'Ramp'
            w_stim = [];
            for i=1:length(input)
                w = rampWave(dt,input(i));
                w.tBefore = 0; w.tRamp = input(i);
                w_stim = extend([w_stim w],sum(input(1:i)));
            end
            w_stim = w_stim.extend(t_total);
    end

    % w_lime = flexRampTrainWave(dt,t_total);
    % w_lime.tBefore = 0.1; % duration off
    % w_lime.tRamp = 0.2; % duration on
    % w_lime.period = 0.5; % seconds

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
    plot(allWaves+[1:5]'*0.3)
    legend cam dmd blue red lime
end