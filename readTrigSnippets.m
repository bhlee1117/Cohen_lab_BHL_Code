function [mov, info] = readTrigSnippets(matFile, which, events)
% READTRIGSNIPPETS  Spike-triggered voltage or glutamate movie snippets written by Section 6 of
% Analysis_20261002_SD_VoltronST_plus_Glu.m (<fp>/TriggerSnippets).
%   [mov, info] = readTrigSnippets(matFile, 'glu')          % all events of one bin
%   [mov, info] = readTrigSnippets(matFile, 'volt', events) % only these events (indices in the bin)
%
%   matFile : Trig_ROI<nn>_<kkk>.mat (the GluTrig_/VoltTrig_ .bin of the same bin is read)
%   mov     : H x W x nFrames x nEvents single (corrected fluorescence)
%   info    : bin information; info.events = the events read, info.(which).frameTimes_s = time of
%             each frame from the event
%
%   Spike-triggered average, e.g.:
%     [g, info] = readTrigSnippets('Trig_ROI02_001.mat', 'glu');
%     F0 = mean(g(:,:,info.glu.frameTimes_s < 0,:), 3);    % pre-event baseline of each event
%     staGlu = mean((g - F0) ./ F0, 4);                      % dF/F movie
%     [v, info] = readTrigSnippets('Trig_ROI02_001.mat', 'volt');
%     F0 = mean(v(:,:,info.volt.frameTimes_s < -0.01,:), 3);
%     staVolt = -mean((v - F0) ./ F0, 4);                    % -dF/F (Voltron dims with spikes)
if nargin < 2 || isempty(which), which = 'glu'; end
S = load(matFile, 'info'); info = S.info;
if nargin < 3 || isempty(events), events = 1:info.nEvents; end
switch lower(which)
    case 'glu',  m = info.glu;  binFile = info.gluFile;
    case 'volt', m = info.volt; binFile = info.voltFile;
    otherwise, error('which must be ''glu'' or ''volt''');
end
H = m.frameSize(1); W = m.frameSize(2); nFr = m.nFramesPerEvent;
evN = H * W * nFr;                                       % values per event
mov = zeros(H, W, nFr, numel(events), 'single');
fid = fopen(fullfile(fileparts(matFile), binFile), 'r');
if fid < 0, error('cannot open %s', binFile); end
c = onCleanup(@() fclose(fid));
for i = 1:numel(events)
    fseek(fid, (events(i) - 1) * evN * 4, 'bof');
    mov(:,:,:,i) = reshape(fread(fid, evN, 'single=>single'), H, W, nFr);
end
info.events = info.events(events);
end
