function allWaves = dragonflyPiezoTestWaveform_yq(nCells,recordFrames,icolor,maxBlueIntensity,centeredROI,exposureTime)
%
%   2018-2019 Vicente Parot
%   Cohen Lab - Harvard university
%
if ~exist('nCells','var')
    nCells = [];
end
if ~exist('recordFrames','var')
    recordFrames = 10000;
end
if ~exist('maxBlueIntensity','var')
    maxBlueIntensity =0.1;
end
maxBlueAmp = maxBlueIntensity;
if ~exist('centeredROI','var')
    centeredROI = [];
end
if ~exist('exposureTime','var')
    exposureTime = 1;
end

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

w_stim = flexRampTrainWave(dt,t_total);
w_stim.tBefore = 4; % duration off
w_stim.tRamp = 1; % duration on
w_stim.period = 5; % seconds

w_lime = flexRampTrainWave(dt,t_total);
w_lime.tBefore = 0.1; % duration off
w_lime.tRamp = 0.2; % duration on
w_lime.period = 0.5; % seconds

switch icolor
    case 'Blue'
        fig = uifigure;
        wave_type = uiconfirm(fig,'Please select illumination waveform sequence','','option',...
                    {'Constant',...
                    'Stair'},...
                    'closefcn',@(h,e) close(fig));
        fprintf('using blue amp %g V\n',maxBlueAmp)
        
        switch wave_type
            case 'Stair'
                opts.Interpreter = 'tex';
                cell_input = inputdlg('\fontsize{10} Blue voltage amplitude sequence (fmt: [...])','',1,{''},opts);

                maxBlueAmp = str2num(cell2mat(cell_input));

                t_i = floor(length(dt:dt:t_total)/length(maxBlueAmp))*dt;
                w_stim = [];
                for i=1:length(maxBlueAmp)
                    w_stim_i = flexRampTrainWave(dt,t_i);
                    w_stim_i.tBefore = 0; % duration off
                    w_stim_i.tRamp = 1; % duration on
                    w_stim_i.period = 1; % seconds
                    w_stim_i = w_stim_i*maxBlueAmp(i);
                    w_stim = [w_stim w_stim_i];
                end
                w_blue = w_stim;
            otherwise
                w_blue = w_stim.*maxBlueAmp;
        end
        
        w_lime = w_stim.*0;
        w_red = w_cam*0+0;
        w_orange = w_cam*0;
    case 'Lime'
%         w_blue = w_stim.*0;
%         w_lime = w_lime.*0+1;
    w_lime = w_lime.*0+1;
        w_blue = w_stim*maxBlueAmp;
        w_red = w_cam*0;
        w_orange = w_cam*0;
    case 'Optopatch'
         
        w_lime = w_lime.*0+1;
        w_blue = w_stim*maxBlueAmp;
        w_red = w_cam*0;
end
w_dmd = extend([wp0 wp1],t_total);

w_piezo = flexRampTrainWave(dt,t_total);
w_piezo.tBefore = 0.999; % duration off
w_piezo.tRamp = 0.001; % duration on
w_piezo.period = 1; % seconds

allWaves = [
    w_cam                 
    w_dmd
    w_blue
    w_red
    w_lime
    w_orange
    w_piezo
    ];
figure(17);clf
plot(allWaves+[1:7]'*0.3)
legend cam dmd blue red lime orange piezo
end