function events = find_coactive_events(GluSpike, varargin)
% find_coactive_events  Detect coactive glutamate events across ROIs.
%
% Usage:
%   events = find_coactive_events(GluSpike)
%   events = find_coactive_events(GluSpike, 'Tolerance', 1, 'MinSize', 2)
%
% Inputs:
%   GluSpike  - N x T  binary matrix  (rows = ROIs, cols = time frames)
%
% Optional name-value pairs:
%   'Tolerance'  - time lag tolerance in frames (default: 1)
%                  ROIs within ±Tolerance frames of each other are
%                  considered coactive
%   'MinSize'    - minimum number of coactive ROIs to report (default: 2)
%
% Output:
%   events  - struct array, one entry per coactive event, with fields:
%     .time      - representative time frame of the event (earliest spike)
%     .rois      - vector of ROI indices that fired
%     .size      - number of coactive ROIs
%     .time_each - vector of actual spike times for each ROI in .rois

%% ---- parse inputs ------------------------------------------------------
p = inputParser;
addRequired(p,  'GluSpike');
addParameter(p, 'Tolerance', 1,  @(x) isnumeric(x) && x >= 0);
addParameter(p, 'MinSize',   2,  @(x) isnumeric(x) && x >= 2);
parse(p, GluSpike, varargin{:});

tol     = p.Results.Tolerance;
minSize = p.Results.MinSize;

[N, T] = size(GluSpike);

%% ---- expand each spike by ±tol to create tolerance windows -------------
% For each ROI, dilate its spike train so any frame within ±tol of a spike
% becomes 1. This way, a simple AND across ROIs detects coactivity.
GluExpanded = false(N, T);
for n = 1:N
    sp = find(GluSpike(n,:));
    for s = 1:numel(sp)
        t1 = max(1,   sp(s) - tol);
        t2 = min(T,   sp(s) + tol);
        GluExpanded(n, t1:t2) = true;
    end
end

%% ---- find frames where >= minSize ROIs are coactive -------------------
coact_count = sum(GluExpanded, 1);          % 1 x T: how many ROIs active
coact_frames = find(coact_count >= minSize);% frames with enough coactivity

if isempty(coact_frames)
    events = struct('time',{}, 'rois',{}, 'size',{}, 'time_each',{});
    fprintf('No coactive events found (MinSize = %d, Tolerance = ±%d).\n', ...
        minSize, tol);
    return
end

%% ---- cluster contiguous coactive frames into single events -------------
% Consecutive frames that are coactive belong to the same event.
breaks     = find(diff(coact_frames) > 1);
seg_starts = [1,            breaks + 1];
seg_ends   = [breaks,       numel(coact_frames)];

nEvents = numel(seg_starts);
events  = struct('time',    cell(nEvents,1), ...
                 'rois',    cell(nEvents,1), ...
                 'size',    cell(nEvents,1), ...
                 'time_each', cell(nEvents,1));

for e = 1:nEvents
    seg_frames = coact_frames(seg_starts(e) : seg_ends(e));

    % Which ROIs fire within this segment (check original spikes, not expanded)
    roi_list    = [];
    time_each   = [];
    for n = 1:N
        sp_in_seg = find(GluSpike(n, seg_frames(1)-tol : min(T, seg_frames(end)+tol)));
        if ~isempty(sp_in_seg)
            % Convert back to global frame index
            t_global = seg_frames(1) - tol - 1 + sp_in_seg;
            t_global = t_global(t_global >= 1 & t_global <= T);
            if ~isempty(t_global)
                roi_list  = [roi_list,  n];
                time_each = [time_each, t_global(1)];  % earliest spike
            end
        end
    end

    events(e).rois      = roi_list;
    events(e).size      = numel(roi_list);
    events(e).time      = min(time_each);       % earliest spike in the group
    events(e).time_each = time_each;
end

%% ---- filter by minSize (post clustering, some events may shrink) -------
keep   = [events.size] >= minSize;
events = events(keep);
nKept  = numel(events);

%% ---- sort by coactivity size (largest first) ---------------------------
[~, order] = sort([events.size], 'descend');
events = events(order);

%% ---- print summary table -----------------------------------------------
fprintf('\n=== Coactive Events (Tolerance = ±%d frames, MinSize = %d) ===\n', ...
    tol, minSize);
fprintf('%-6s  %-10s  %-8s  %s\n', 'Rank', 'Time (fr)', 'Size', 'ROIs');
fprintf('%s\n', repmat('-', 1, 55));
for e = 1:nKept
    roi_str = sprintf('%d ', events(e).rois);
    fprintf('%-6d  %-10d  %-8d  [%s]\n', e, events(e).time, events(e).size, strtrim(roi_str));
end
fprintf('\nTotal: %d events found.\n\n', nKept);

end
