% NI DAQ channel definition
% outChannelSettings = {
%     'cameraTrigger',    {'Dev1','Port0/Line31', 'OutputOnly'}
%     'dmdTrigger',       {'Dev1','Port0/Line8',  'OutputOnly'}
%     'blueIntensity',    {'Dev1','ao2',          'Voltage'   }
%     'redShutter',       {'Dev1','Port0/Line22', 'OutputOnly'}
%     'limeShutter',      {'Dev1','Port0/Line23', 'OutputOnly'}
%     'orangeShutter',    {'Dev1','Port0/Line28', 'OutputOnly'}
%     'piezoTrigger',     {'Dev1','Port0/Line0',  'OutputOnly'}
%     };
%% ----------- magnet test ----------------------
s = daq.createSession('NI');
s.addDigitalChannel('Dev1','Port0/Line27', 'OutputOnly');
s.addDigitalChannel('Dev1','Port0/Line14', 'OutputOnly');
dt = 1e-3; % digitization rate
% s = daq.createSession('NI');
% s.addDigitalChannel('Dev1','Port0/Line28', 'OutputOnly');
s.Rate = 1/dt;
%%
s.outputSingleScan([1 1])
%%
while true
s.outputSingleScan([1 0])
pause(.1)
s.outputSingleScan([1 1])
pause(1)
end
%% ----------- pointing noise test ----------------------
s = daq.createSession('NI');
s.addDigitalChannel('Dev1','Port0/Line28', 'OutputOnly');
s.addDigitalChannel('Dev1','Port0/Line25', 'OutputOnly');
dt = 1e-3; % digitization rate
% s = daq.createSession('NI');
s.addAnalogInputChannel('Dev1','ai0', 'Voltage'); 
% s.addDigitalChannel('Dev1','Port0/Line28', 'OutputOnly');
s.Rate = 1/dt;
%%
t_total = 10; % duration
w_on = flexRampTrainWave(dt,t_total);
w_on.tBefore = 0; % duration off
w_on.tRamp = t_total; % duration on
w_on.period = t_total; % seconds
%%
s.queueOutputData([w_on.amplitude' w_on.amplitude']);

[data_in,t] = s.startForeground;

figure;
plot(t,data_in)
figure
plot_freq_spec(data_in,dt)
title('on')
%%
s.outputSingleScan([1 1])
%%----------- pointing noise test END ----------------------
%%
%%-----------594 shutter control start------------------

s = daq.createSession('NI');
s.addDigitalChannel('Dev1','Port0/Line28', 'OutputOnly');

%%
s.outputSingleScan(1)
%%
%%-----------594 shutter control END -------------------

%% --------- OD wheel test --------------------
dt = 1e-4; % digitization rate
s = daq.createSession('NI');
s.addDigitalChannel('Dev1','Port0/Line0', 'OutputOnly');
s.addAnalogOutputChannel('Dev1','ao2',          'Voltage')
s.Rate = 1/dt;
%%
f = thorlabsFilterWheel(s);
%%
moveTo(f,'0')
% out = [ones(round(.1/dt),1);zeros(round(50e-3/dt),1);ones(round(2/dt),1)];
% all_out = zeros(length(out),2);
% all_out(:,1) = out;
% s.queueOutputData(all_out)
% s.startForeground
%%
s.release
clear
%%--------- OD wheel test END --------------------

%% --------- Pololu test --------------------
dt = 1e-4; % digitization rate
s = daq.createSession('NI');
s.addDigitalChannel('Dev1','Port0/Line25', 'OutputOnly');
s.Rate = 1/dt;


%--------- Pololu test END--------------------

%% --------SLM cover glass voltage test------

dt = 1e-4; % digitization rate
s = daq.createSession('NI');
s.addAnalogInputChannel('Dev1','ai0', 'Voltage'); 
s.addDigitalChannel('Dev1','Port0/Line28', 'OutputOnly');
s.Rate = 1/dt;
% s.DurationInSeconds = t_total;

%% 594 shutter waveform
t_total = .5; % duration
w_on = flexRampTrainWave(dt,t_total);
w_on.tBefore = 0; % duration off
w_on.tRamp = t_total; % duration on
w_on.period = t_total; % seconds

%% set cover glass voltage
SDK_Filepath='C:\Program Files\Meadowlark Optics\Blink 1920 HDMI\SDK';
LUT_Filepath=fullfile('C:\Program Files\Meadowlark Optics\Blink 1920 HDMI\LUT Files\','slm4633_at594.lut');
loadlibrary(fullfile(SDK_Filepath,'Blink_C_wrapper.dll'), fullfile(SDK_Filepath,'Blink_C_wrapper.h'));
calllib('Blink_C_wrapper','Create_SDK');
%%


im = imread('C:\Program Files\Meadowlark Optics\Blink 1920 HDMI\Image Files\square_2x1_GS.bmp');

% disp('Blink SDK was successfully constructed')
calllib('Blink_C_wrapper', 'Load_lut', char(LUT_Filepath));

% disp('SLM Ready to go')
calllib('Blink_C_wrapper', 'Write_image', flip(im,1)', 1);

%%
volts_all = 6.01:.01:6.399;
for it_v = 1:length(volts_all)
    volts = volts_all(it_v)
if volts>6 && volts<6.4
    x=calllib('Blink_C_wrapper','Set_SLMVCom',volts);
    if x==false
        error('Coverglass Voltage Write Failed')
    end
else
    warning('Did not adjust coverglass voltage. Requested value out of range')
end

for it = 1:20
s.queueOutputData(w_on.amplitude');

[data_in,t] = s.startForeground;

L = t_total/dt;
data_in = data_in(L*5e-2+1:L*95e-2);
figure(98);clf;plot(data_in)
% calculate fft
 L = length(data_in);
X = fft(data_in); 
P = abs(X/L); P = P(1:L/2+1); P(2:end-1) = 2*P(2:end-1);
F = 1/dt*(0:L/2)/L;
[amp,loc] = findpeaks(P(F>400 & F<450),F(F>400 & F<450),'NPeaks',1,'sortstr','descend');
figure(99);clf;plot(F,P);ylim([0 max(P(F>100))*1.5])
hold on;plot(loc,amp,'v')
noise_data(it_v,it) = amp;
end
end

%% close SLM and daq
s.release
calllib('Blink_C_wrapper', 'Delete_SDK');
unloadlibrary('Blink_C_wrapper');


% --------SLM cover glass voltage test END------
%%
session = daq.createSession('NI');
session.addDigitalChannel('Dev1','Port0/Line28', 'OutputOnly'); 
%%
session = daq.createSession('NI');
session.addDigitalChannel('dev1','port0/line0','outputonly'); 
%%
dt = 1e-3;t_total = 600;
for i=1:length(0:dt:t_total)
if mod(i,100)==0
session.outputSingleScan(1)
elseif mod(i,100)==1
session.outputSingleScan(0)
end
pause(dt)
end

session.outputSingleScan(0)

%%

s = daq.createSession('ni');
s.addAnalogOutputChannel('Dev1','ao2',          'Voltage')
% s.addDigitalChannel('Dev1','Port0/Line23', 'OutputOnly')
s.Rate = 1e+5;
dt = 1e-5;
t_total = 60;
%%

t_wait=0;
% % w_stim = sin((dt:dt:t_total)/dt/2/pi/100/dt);
% w_stim = flexRampTrainWave(dt,t_total);
% w_stim.tBefore = 0.2e-3; % duration off
% w_stim.tRamp = 10e-3; % duration on
% w_stim.period = 50e-3; % seconds

w_wait = baseWave(dt,t_wait)*0;
w_stop = baseWave(dt,0.5)*0;

w_blue = baseWave(dt,t_total);
w_blue = [w_wait w_blue*0+1 w_stop];

% w_lime =  baseWave(dt,t_total);
% w_lime = [w_wait w_lime*0+1 w_stop];

all_waves = w_blue*10;
% all_waves = [w_blue*1.1;w_lime];
figure;plot(w_blue*1.1)


% figure;plot(w_stim)
%%
% s.queueOutputData(w_stim.amplitude'*0.1)
s.queueOutputData(all_waves.amplitude')
s.startForeground

%%
s.release