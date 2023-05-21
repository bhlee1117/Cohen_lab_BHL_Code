%% Test script as example and to try waveform definitions. 
% Execute line by line and watch output of commands. 
%
%   2018-2019 Vicente Parot
%   Cohen Lab - Harvard university
%
tic
nCells = 10; % 1 second pre, 1 second per cell, 5 seconds post
maxBlueIntensity = 25; % in mW/cm2
dt = 1e-5; % constructors can take a 2nd parameter with the duration in s
z = baseWave(dt,5e-3); % zero valued constant wave of duration 5 ms, to preopen shutter and buffer other waves
f = baseWave(dt,nCells+6); % a constant wave, zero valued by default
do_redw = [z+1 f+1 z];
f = rampWave(dt,1); f.tBefore = 0; f.tRamp = .5; % short repeated
g = rampWave(dt,5); g.tBefore = .5; g.tRamp = 4; % long segment
% ao_bluew = map( [z f*0 rep(f,nCells) g z], @(v)blueLUTmW(v,maxBlueIntensity) );
f = flexRampTrainWave(dt,nCells+6); f.period = 1e-3; f.dc = .1; f.phase = 0;
g = flexRampTrainWave(dt,1e-3);     g.period = 1e-3; g.dc = .1; g.phase = 0;
% do_camw = [z f g rep(g,4)*0]; % 1 extra pulse at the end
f = flexRampWave(dt,1); f.tBefore = .75; f.tRamp = .01; f.initialAmp = 1;
% do_dmdw = [z rep(f,nCells+1) 0*rep(f,5) z];

% clf
% hold on
figure
% plot(do_camw)
% plot(do_dmdw)
% plot(ao_bluew)
plot(do_redw)
axis tight
legend 'camera trig' 'DMD trig' 'blue intensity' 'red shutter'
toc