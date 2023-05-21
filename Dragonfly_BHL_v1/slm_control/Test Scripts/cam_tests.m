%% NI DAQ init
% s = daq.createSession('NI');
% addDigitalChannel(s,'Dev1','Port0/Line22','OutputOnly');
% addAnalogInputChannel(s,'Dev1','ai0','Voltage');
% % outputSingleScan(s,1)
% % outputSingleScan(s,0)
% % queueOutputData(s,[ones(200,1);zeros(1,1)]);
% % s.startForeground;
% % clear s

%% camera init
% imaqreset

vid = videoinput('hamamatsu', 1, 'MONO16_2048x2048_FastMode');
src = getselectedsource(vid);
%% set ROI
% square centered ROI
bin = 1;
roi_nr = 512;
roi_nc = 512;
roi_or = (2048/bin-roi_nr)/2;
roi_oc = (2048/bin-roi_nc)/2;
vid.ROIPosition = [roi_oc roi_or roi_nc roi_nr];
hImage = mypreview(vid);
% hImage.CDataMapping = 'scaled';
% hImage.Parent.CLim = [0 3];
% set(hImage,'cdatamapping','scaled')
% hold(hImage.Parent,'on')

%% PSF optimization with zernike coefficients
% refocus and avoid photobleaching
% preview(vid);
src.ExposureTime = .001;
return
%% release camera and daq
try stoppreview(vid); end
clear src
clear vid
clear s
imaqreset


