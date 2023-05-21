Dev1/port0/line22
d = daq.getDevices
d(1)
s = daq.createSession('NI');
addDigitalChannel(s,'Dev1','Port0/Line22','OutputOnly');
addAnalogInputChannel(s,'Dev1','ai0','Voltage');

% [ch, idx] = addDigitalChannel(s,'Dev1','Port0/Line22','OutputOnly');
% tic
% for it = 1:4
%     outputSingleScan(s,1)
%     pause(.01)
%     outputSingleScan(s,0)
%     pause(1)
% end
% toc
queueOutputData(s,[ones(1000,1);zeros(1,1)]);
tic
s.startForeground;
toc
