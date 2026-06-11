function events = find_spatiotemporal_events(GluSpike, DistMat, coords, varargin)
% find_spatiotemporal_events  Scan frame-by-frame and report coactive events
%                              classified as clustered or isolated.
%
% Usage:
%   events = find_spatiotemporal_events(GluSpike, DistMat, coords)
%   events = find_spatiotemporal_events(GluSpike, DistMat, coords, ...
%               'Tolerance', 1, 'MaxDist', 50)
%
% Inputs:
%   GluSpike - N x T  binary spike matrix (rows=ROIs, cols=frames)
%   DistMat  - N x N  pairwise ROI distance matrix (µm)
%   coords   - N x 2  [x y] centroid coordinates of each ROI (pixels or µm)
%
% Optional name-value pairs:
%   'Tolerance' - time lag tolerance in frames (default: 1)
%                 ROI i (firing at t) and ROI j (firing at t') are coactive
%                 if |t - t'| <= Tolerance
%   'MaxDist'   - max spatial distance (µm) to be "clustered" (default: 50)
%
% Algorithm:
%   For each frame t where at least one ROI fires:
%     1. Collect all ROIs that fired within [t-tol, t+tol] that have NOT
%        already been assigned to a previous event (no double-counting).
%     2. Build a pairwise adjacency: ROIs i,j are linked if BOTH
%        |time_i - time_j| <= tol  AND  DistMat(i,j) <= MaxDist.
%     3. Union-Find gives spatial clusters among the coactive ROIs.
%     4. Each spatial cluster (size>=2) is one "clustered" event entry.
%        Each singleton is one "isolated" event entry.
%   Events are ordered by time (earliest spike in each group).
%
% Output:
%   events - struct array sorted by time, one row per spatial group, fields:
%     .time          - earliest spike frame in this group
%     .time_each     - spike frame of each ROI in this group
%     .rois          - ROI indices in this group
%     .size          - number of ROIs
%     .type          - 'clustered' or 'isolated'
%     .centroid      - [x y] mean centroid of active ROIs in this group

%% ---- parse inputs -------------------------------------------------------
p = inputParser;
addRequired(p,  'GluSpike');
addRequired(p,  'DistMat');
addRequired(p,  'coords');
addParameter(p, 'Tolerance', 1,  @(x) isnumeric(x) && x >= 0);
addParameter(p, 'MaxDist',   50, @(x) isnumeric(x) && x >  0);
parse(p, GluSpike, DistMat, coords, varargin{:});

tol     = p.Results.Tolerance;
maxDist = p.Results.MaxDist;

[N, T] = size(GluSpike);
assert(isequal(size(DistMat), [N N]), 'DistMat must be N x N.');
assert(size(coords,1) == N,           'coords must be N x 2.');

%% ---- Step 1: build spike time list per ROI ----------------------------
% spike_times{n} = vector of frame indices where ROI n fired
spike_times = cell(N, 1);
for n = 1:N
    spike_times{n} = find(GluSpike(n,:));
end

%% ---- Step 2: scan each spike and group coactive ROIs ------------------
% Track which (ROI, spike_index) pairs have already been assigned
assigned = false(N, T);   % assigned(n,t) = true once spike at t is used

events = empty_events();
nOut   = 0;

% Get all spike frames in chronological order
[all_rois, all_frames] = find(GluSpike);
[all_frames_sorted, sort_idx] = sort(all_frames);
all_rois_sorted = all_rois(sort_idx);


for si = 1:numel(all_frames_sorted)
    t_seed = all_frames_sorted(si);
    n_seed = all_rois_sorted(si);

    % Skip if this spike already assigned to an event
    if assigned(n_seed, t_seed), continue; end

    %-- Collect all ROIs with a spike within [t_seed-tol, t_seed+tol]
    %   that are not yet assigned
    t_lo = max(1, t_seed - tol);
    t_hi = min(T, t_seed + tol);

    group_rois  = [];
    group_times = [];

    for n = 1:N
        sp = spike_times{n};
        sp_in_win = sp(sp >= t_lo & sp <= t_hi & ~assigned(n, sp));
        % Take earliest unassigned spike within window
        sp_in_win_unassigned = sp_in_win(~assigned(n, sp_in_win));
        if ~isempty(sp_in_win_unassigned)
            t_hit = sp_in_win_unassigned(1);   % earliest spike in window
            group_rois  = [group_rois,  n];
            group_times = [group_times, t_hit];
        end
    end

    if isempty(group_rois), continue; end

    %-- Step 3: pairwise adjacency — both temporally AND spatially close
    nR  = numel(group_rois);
    adj = false(nR, nR);
    for a = 1:nR
        for b = a+1:nR
            dt = abs(group_times(a) - group_times(b));
            dd = DistMat(group_rois(a), group_rois(b));
            if dt <= tol && dd <= maxDist
                adj(a,b) = true;
                adj(b,a) = true;
            end
        end
    end

    %-- Step 4: union-find → spatial clusters within this coactive group
    labels        = union_find_clusters(adj, nR);
    unique_labels = unique(labels);

    for c = 1:numel(unique_labels)
        mem_local = find(labels == unique_labels(c));
        mem_rois  = group_rois(mem_local);
        mem_times = group_times(mem_local);

        % Mark these spikes as assigned
        for k = 1:numel(mem_rois)
            assigned(mem_rois(k), mem_times(k)) = true;
        end

        % Classify
        if numel(mem_rois) >= 2
            ev_type = 'clustered';
        else
            ev_type = 'isolated';
        end

        % Centroid of active ROIs
        ev_centroid = mean(coords(mem_rois, :), 1);

        nOut = nOut + 1;
        events(nOut).time      = min(mem_times);
        events(nOut).time_each = mem_times;
        events(nOut).rois      = mem_rois;
        events(nOut).size      = numel(mem_rois);
        events(nOut).type      = ev_type;
        events(nOut).centroid  = ev_centroid;
    end
end

%% ---- Step 5: sort by time ----------------------------------------------
if isempty(events), return; end
[~, order] = sort([events.time], 'ascend');
events = events(order);

%% ---- Step 6: print summary ---------------------------------------------
nKept = numel(events);
nClust = sum(strcmp({events.type}, 'clustered'));
nIso   = sum(strcmp({events.type}, 'isolated'));

% fprintf('\n=== Spatiotemporal Events  (Tol=±%d fr | MaxDist=%.0f µm) ===\n', ...
%     tol, maxDist);
% fprintf('Total: %d events  (%d clustered, %d isolated)\n\n', ...
%     nKept, nClust, nIso);
% fprintf('%-6s  %-12s  %-10s  %-6s  %-20s  %-15s  %s\n', ...
%     'Event','Type','Time (fr)','Size','ROIs','Time each','Centroid [x y]');
% fprintf('%s\n', repmat('-', 1, 95));

% for e = 1:nKept
%     ev      = events(e);
%     roi_str = sprintf('%d ', ev.rois);
%     t_str   = sprintf('%d ', ev.time_each);
%     cen_str = sprintf('%.1f %.1f', ev.centroid(1), ev.centroid(2));
%     fprintf('%-6d  %-12s  %-10d  %-6d  %-20s  %-15s  %s\n', ...
%         e, ev.type, ev.time, ev.size, ...
%         strtrim(roi_str), strtrim(t_str), cen_str);
% end
fprintf('\n');
end

%% =========================================================================
function labels = union_find_clusters(adj, n)
parent = 1:n;
    function r = find_root(x)
        while parent(x) ~= x
            parent(x) = parent(parent(x));
            x = parent(x);
        end
        r = x;
    end
    function union_nodes(x, y)
        rx = find_root(x);
        ry = find_root(y);
        if rx ~= ry, parent(ry) = rx; end
    end
for i = 1:n
    for j = i+1:n
        if adj(i,j), union_nodes(i,j); end
    end
end
roots = arrayfun(@find_root, 1:n);
[~, ~, labels] = unique(roots);
labels = labels';
end

%% =========================================================================
function ev = empty_events()
ev = struct('time',{}, 'time_each',{}, 'rois',{}, 'size',{}, ...
            'type',{}, 'centroid',{});
end
