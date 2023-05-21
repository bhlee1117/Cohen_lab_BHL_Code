function allWaves = dragonflyL1PatternWaveform_yq(recordFrames,icolor,maxBlueIntensity,exposureTime)
    

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
    %% camera waveform
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
    
    %% blue waveform
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
        case 'Ramp'
            w_stim = [];
            for i=1:length(input)
                w = rampWave(dt,input(i));
                w.tBefore = 0; w.tRamp = input(i);
                w_stim = extend([w_stim w],sum(input(1:i)));
            end
            w_stim = w_stim.extend(t_total);
    end
    %% dmd
    w_dmd = wp0;
    for i=1:length(input)
        w_dmd = extend([w_dmd wp1],sum(input(1:i)));
    end
    w_dmd = w_dmd.extend(t_total);
    %% repeat waveform
    fig = uifigure;
    repeat_q = uiconfirm(fig,'Would you like to repeat?','','option',...
                            {'Yes',...
                            'No'},...
                            'closefcn',@(h,e) close(fig));
    if strcmp(repeat_q,'Yes')
        opts.Interpreter = 'tex';
        cell_input = inputdlg('\fontsize{10} Please specify repeat and wait time (fmt: repeat duration ...])','',1,{''},opts);
        global input_repeat_wait;
        input_repeat_wait = str2num(cell2mat(cell_input));
        
        w_cam = [w_cam extend(wp0,input_repeat_wait(2))];
        w_cam = w_cam.repeat(input_repeat_wait(1));
        
        w_stim = [w_stim extend(wp0,input_repeat_wait(2))];
        w_stim = w_stim.repeat(input_repeat_wait(1));
        
        w_dmd = [w_dmd extend(wp0,input_repeat_wait(2))];
        w_dmd = w_dmd.repeat(input_repeat_wait(1));
    end

    
    %% 
    
    switch icolor
        case 'Blue'
            w_blue = w_stim;
            w_lime = w_stim.*0;
            w_red = w_cam*0+0;
        case 'Lime'
            w_blue = w_stim;
            w_lime = w_stim./max(w_stim.amplitude);
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
    line(dt:dt:allWaves(1).duration,w_lime.amplitude)
%     line(dt:dt:allWaves(1).duration,allWaves.amplitude+[1:5]'*0.3)
%     legend cam dmd blue red lime
end