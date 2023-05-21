function allWaves = dragonflySimpleWaveform_yq(nCells,recordFrames,icolor,maxBlueIntensity,centeredROI,exposureTime)
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




%----cAmp long blue only sequence -----
% w_stim = flexRampTrainWave(dt,t_total);
% w_stim.tBefore = 0; % duration off
% w_stim.tRamp = 60; % duration on
% w_stim.period = 240; % seconds
% w_stim = extend(w_stim,t_total);
% 
% w_dmd = extend([wp0 wp1],60);
% w_dmd = extend([w_dmd wp1],t_total);

%------------end---------------------




%----cAmp switching sequence -----
% w_stim = flexRampTrainWave(dt,t_total);
% w_stim.tBefore = 0; % duration off
% w_stim.tRamp = 5; % duration on
% w_stim.period = 6; % seconds
% w_stim = extend(w_stim,t_total);
% % w_stim = extend([w_stim extend(wp0,180)],t_total);
% 
% w_dmd = extend([wp0 wp1],5);
% w_dmd = extend([w_dmd wp1],6);
% w_dmd = repeat(w_dmd,30);
% w_dmd = extend(w_dmd,t_total);
%------------end---------------------

%----constant stim sequence------
w_stim = flexRampTrainWave(dt,t_total);
w_stim.tBefore = 1.5; % duration off 
w_stim.tRamp = 2; % duration on
w_stim.period = 5; % seconds

w_dmd = extend([wp0 wp1],t_total);
%------------end---------------------

%----square train stim sequence------
% w_stim = flexRampTrainWave(dt,t_total);
% w_stim.tBefore = 0; % duration off
% w_stim.tRamp = 5; % duration on
% w_stim.period = 5; % seconds
% w_dmd = extend([wp0],t_total);
% w_dmd = flexRampTrainWave(dt,t_total);
% w_dmd.tBefore = 0;
% w_dmd.tRamp = .01;
% w_dmd.period = 2;
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
% ---------- end --------------

%----stair sequence---
% w_stim = [];
% w_wait = baseWave(dt,.5); 
% for amp_factor = logspace(-1.92,0,10) %0.1-5 V
%     w_stim = [w_stim w_wait baseWave(dt,1.5)+1*amp_factor];
% end
% w_stim = w_stim.extend(t_total);
% w_dmd = extend([wp0],t_total);
%---------- end --------------


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
        w_orange = w_cam*0;
    case 'Lime'
%         w_blue = w_stim.*0;
%         w_lime = w_lime.*0+1;
        w_lime = w_lime.*0;
        w_orange = w_lime.*0+1;
        w_pololu = w_lime.*0+1;
        w_blue = w_stim*maxBlueAmp;
%         w_red = w_cam*0;
%         w_orange = w_cam*0;
% 
% w_blue = w_stim*maxBlueAmp;
% w_lime = w_stim*0+1-w_stim;
% w_blue = extend(w_blue,t_total);
% w_lime = extend(w_lime,t_total);

    case 'Optopatch'
         
        w_lime = w_lime.*0+1;
        w_blue = w_stim*maxBlueAmp;
        
end

w_red = w_cam*0;
w_piezo = w_cam*0;
w_odwheel = w_cam*0;


allWaves = [
    w_cam                 
    w_dmd
    w_blue
    w_red
    w_lime
    w_orange
    w_piezo
    w_pololu
    w_odwheel
    ];
figure
plot(allWaves+[1:9]')
legend cam dmd blue red lime orange
end

