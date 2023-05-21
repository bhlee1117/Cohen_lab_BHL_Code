%% execute this code block to shutdown, compile, start, go live
try %#ok<TRYNC>
    cameraMex shutdown
    clear cameraMex
end
clc
mex cameraMex.c camControl.c graphicThreads.c dcimg2bin.c common.cpp ...
    -IC:\Users\Labmember\Desktop\Gui_data\libs\SDL2-vc\include  -LC:\Users\Labmember\Desktop\Gui_data\libs\SDL2-vc\lib\x64      -lSDL2 ...
    -IC:\Users\Labmember\Desktop\Gui_data\libs\dcamsdk4\inc\    -LC:\Users\Labmember\Desktop\Gui_data\libs\dcamsdk4\lib\win64   -ldcamapi ...
    -IC:\Users\Labmember\Desktop\Gui_data\libs\dcimgsdk\inc\    -LC:\Users\Labmember\Desktop\Gui_data\libs\dcimgsdk\lib\win64   -ldcimgapi
% cameraMex('startup') % default params
% cameraMex('startup',1,955,1174,2,24) % default params
% cameraMex('startup',0,400,600,300,100) % custom screen in monitor 0, with 400 wide 600 tall area and offset (300,100) (px)
cameraMex('startup',0,790,830,2,24) % custom screen in monitor 0, with 400 wide 600 tall area and offset (300,100) (px)
% cameraMex('startup',0) % wont start: must have 0 or 5 params
%
            exposureTime = .0001;
            inputROI = [2048 2048];
%%
            exposureTime = .1;
inputROI = [950 100 2048 94];
%
inputROI = [2048 2048];
binning = 2;
[a, b] = cameraMex('aq_live_restart',inputROI,binning,exposureTime), pause(.1)
% size(cameraMex('aq_snap'))
%%
% cameraMex aq_thread_snap
% cameraMex aq_snap;
%% live aq protocol
try, cameraMex startup, end
            exposureTime = .001;
            inputROI = [512 256];
[a, b] = cameraMex('aq_live_restart',inputROI,binning,exposureTime), pause(.1)
ims = imshow(cameraMex('aq_snap'),[]);
cameraMex verbose
tic
for it = 1:100
    ims.CData = cameraMex('aq_snap');
    drawnow
end
toc
cameraMex verbose
cameraMex aq_live_stop
cameraMex shutdown
clear cameraMex
%% sync aq protocol
cameraMex startup
%%
            exposureTime = 1e-3;
            binning = 2;
            inputROI = [2048 160]; % [2048 160] max dimensions at full throughput
% from here to repeat with updated exposure, roi, and length
try t.TimerFcn = ';'; end% disables the timer if running
cameraMex aq_sync_stop % stop to close the file if it was open
cameraMex('aq_sync_prepare',inputROI,binning,exposureTime)
                        % play waveforms synchronized with camera aq
                        clear nis
                        nis = daq.createSession('NI');
                        nis.addDigitalChannel('Dev1','Port0/Line31','OutputOnly'); % camera clock
                        nis.addDigitalChannel('Dev1','Port0/Line22','OutputOnly'); % shutter, for oscilloscope monitoring
                        nis.Rate = 1e5; % approx 102 kHz
                        nicl = nis.addClockConnection('External','Dev1/PFI0','ScanClock'); % camera programmable out clock 102 kHz
                        samplesPerFrame = nis.Rate*exposureTime;
% from here to repeat with updated length, but same exposure and roi
            recordFrames = 1000; % (N)
                        outputSignal1 = zeros(samplesPerFrame,1); % 100 samples for 1 khz rep
                        outputSignal1(1) = 1;
                        outputSignal1 = repmat([outputSignal1 outputSignal1*0],[recordFrames+1,1]);
% from here to repeat with same exposure, roi, and length
                        nis.queueOutputData(outputSignal1);
                        niDuration = nis.DurationInSeconds;
dcimg_fullpath = fullfile('R:\three.dcimg');
bin_fullpath = fullfile('R:\four.bin');
[fp, fn, fe] = fileparts(dcimg_fullpath);
dcimg_fpath = strrep(fullfile(fp,fn), filesep, [filesep filesep]);
dcimg_fext = strrep(fe, '.', '');
delete(dcimg_fullpath)
cameraMex('aq_sync_start',recordFrames,dcimg_fpath,dcimg_fext,dcimg_fpath,'bin',dcimg_fpath,'bin')
                        nis.startBackground
                        why
% for it = 1:1000
%     ims.CData = cameraMex('aq_snap');
%     drawnow
% end
% cmdstr = {
%     'cameraMex(''aq_live_restart'',inputROI,exposureTime),' % includes stop
%     'disp(''live mode running'')' % includes stop
%     };
% t = timer(...
%     'TimerFcn', sprintf('%s',cmdstr{:}),...
%     'StartDelay',niDuration,'ExecutionMode','singleShot');
% t.start % launch timer to close the file automatically

%
% cameraMex aq_sync_stop
cameraMex('aq_live_restart',inputROI,binning,exposureTime), pause(.1)
%%
cameraMex shutdown
%%
clc
mex dcimg2bin.c -ID:\Code\dcimgsdk\inc\ -LD:\Code\dcimgsdk\lib\win64 -ldcimgapi
dcimg2bin 'R:\testfi ble.dcimg'






