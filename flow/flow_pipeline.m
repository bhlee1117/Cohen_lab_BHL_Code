function R = flow_pipeline(data, cfg, R, upto)
% R = flow_pipeline(data, cfg)
% R = flow_pipeline(data, cfg, R)          % re-run only what the config change affects
% R = flow_pipeline(data, cfg, R, 'segments')   % ...and stop after that stage
%
% Runs the vascular-flow stages in order and caches them by config stamp. Pass the
% previous R back in and only the stages whose parameters (or whose upstream) changed are
% recomputed; everything else is reused untouched.
%
% WHY THIS EXISTS. As a six-cell script the stages communicated through ~100 loose
% workspace variables. Nothing recorded which cell owned which variable, so after changing
% a parameter the only safe move was to re-run the whole thing -- and the tempting move
% was to re-run one cell, which silently mixes new settings with stale results. Here the
% dependency is explicit: each stage's stamp folds in its own parameters AND its
% upstream's stamp, so changing the band-pass necessarily invalidates the tracking, the
% segments and everything after, with no one having to remember that.
%
% INPUT
%   data  struct describing the recording, built once by the caller:
%     .mov_res  [nrow x ncol x nT]  motion-residual movie (mean out, mcTrace regressed)
%     .fs       Hz          .pxUm  um per pixel
%     .badfr    [1 x nT] logical, peri-spike and blue-on frames
%     .spike    [1 x nT] logical soma spike train
%     .ftprnt   [nrow x ncol x nROI]  neuron footprints
%     .bvBW     [nrow x ncol] logical vessel mask
%     .neuronMap[nrow x ncol] max over ftprnt
%     .voltTr   [nV x nT] voltage traces to project out ([] = skip)
%     .tag      char, identifies the recording (guards against reusing R across datasets)
%   cfg   from flow_defaults()
%   R     previous result ([] or omitted to start fresh)
%   upto  stop after this stage: 'bandpass'|'track'|'segments'|'regressors'|'apply'|
%         'compare'  (default: run everything)
%
% OUTPUT  R.<stage> for each stage that ran, plus R.stamp.<stage> and R.tag.
%
% 2026.08.28, Byung Hun Lee

data = flow_checkData(data);    % fail here, with the field name, not four stages in
cfg  = flow_checkCfg(cfg);      % ...and the same for a cfg saved before a knob existed
if nargin < 3, R = struct; end
if nargin < 4 || isempty(upto), upto = 'compare'; end
if isempty(R), R = struct; end
stages = {'bandpass','track','segments','regressors','apply','compare'};
iStop  = find(strcmp(stages, upto), 1);
assert(~isempty(iStop), 'upto must be one of: %s', strjoin(stages, ', '));

% A cached result from a DIFFERENT recording is not a cache hit, however the config looks.
if ~isfield(data,'tag'), data.tag = ''; end
if isfield(R,'tag') && ~strcmp(R.tag, data.tag)
    fprintf('[flow] data tag changed (%s -> %s): discarding the cached result\n', ...
        R.tag, data.tag);
    R = struct;
end
R.tag = data.tag;
if ~isfield(R,'stamp'), R.stamp = struct; end

sz = [size(data.mov_res,1) size(data.mov_res,2)];
ref = struct('bvBW', data.bvBW, 'pxUm', data.pxUm, 'fs', data.fs, ...
             'neuronMap', data.neuronMap);

want = struct();
want.bandpass   = flow_stamp(cfg.bandpass, cfg.common);
want.track      = flow_stamp(cfg.track,       want.bandpass);
want.segments   = flow_stamp(cfg.segments,    want.track);
want.regressors = flow_stamp(cfg.regressors,  want.segments);
want.apply      = flow_stamp(cfg.apply,       want.regressors);
want.compare    = flow_stamp(cfg.compare,     want.apply);

for k = 1:iStop
    st = stages{k};
    if isfield(R, st) && isfield(R.stamp, st) && strcmp(R.stamp.(st), want.(st))
        fprintf('[flow] %-11s cached (%s)\n', st, want.(st)(1:8));
        continue
    end
    t0 = tic;
    switch st
        case 'bandpass'
            R.(st) = flow_bandpass(data, cfg);
        case 'track'
            R.(st) = flow_track(R.bandpass, cfg, ref);
        case 'segments'
            R.(st) = flow_segments(R.track, sz, cfg, ref);
        case 'regressors'
            R.(st) = flow_regressors(R.bandpass, R.track, R.segments, data, cfg);
        case 'apply'
            R.(st) = flow_apply(data, R.bandpass, R.segments, R.regressors, cfg);
        case 'compare'
            R.(st) = flow_compare(data, R.apply, R.regressors, cfg);
    end
    R.stamp.(st) = want.(st);
    fprintf('[flow] %-11s done in %.1f s (%s)\n', st, toc(t0), want.(st)(1:8));

    % a stage that re-ran invalidates everything after it
    for j = k+1:numel(stages)
        if isfield(R, stages{j}), R = rmfield(R, stages{j}); end
        if isfield(R.stamp, stages{j}), R.stamp = rmfield(R.stamp, stages{j}); end
    end
end
end
