function out = matchSglu2SWC(dendrites, spotsXY, varargin)
%matchSglu2SWC Assign iGluSnFR spots to neurons by proximity to dendrite nodes.
%
% out = matchSglu2SWC(dendrites, spotsXY)
% out = matchSglu2SWC(dendrites, spotsXY, 'MaxDist', 5)
%
% INPUTS
%   dendrites : 1x4 cell array, each cell is [m_i x 2] (XY dendrite node coords) for neuron i
%               Example: dendrites = {xy1, xy2, xy3, xy4};
%   spotsXY   : [N x 2] (XY coords of iGluSnFR spots)
%
% Name-Value options
%   'MaxDist' : (default = Inf) maximum distance allowed to assign a spot.
%               Spots farther than MaxDist are labeled 0 (unassigned).
%   'UseToolbox' : (default = true) use knnsearch if available. If false, uses fallback.
%
% OUTPUT (struct) fields
%   out.neuronID      : [N x 1] assigned neuron index (1..4), or 0 if unassigned by MaxDist
%   out.minDist       : [N x 1] distance to assigned neuron (Inf if unassigned)
%   out.closestNode   : [N x 1] index of closest dendrite node within assigned neuron
%   out.closestNodeXY : [N x 2] coords of closest node (NaN if unassigned)
%
% Notes:
%   - Proximity is computed to discrete dendrite nodes, not line segments.
%     If you want point-to-segment distance along dendrite edges, tell me and I'll adapt it.

% -------------------- parse inputs --------------------
p = inputParser;
p.addRequired('dendrites', @(c) iscell(c) );
p.addRequired('spotsXY', @(x) isnumeric(x) && size(x,2)==2);
p.addParameter('MaxDist', Inf, @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('UseToolbox', true, @(x) islogical(x) && isscalar(x));
p.parse(dendrites, spotsXY, varargin{:});

maxDist = p.Results.MaxDist;
useToolbox = p.Results.UseToolbox;

N = size(spotsXY,1);

% Validate dendrites
for i = 1:length(dendrites)
    if ~isnumeric(dendrites{i}) || size(dendrites{i},2) ~= 2
        error('dendrites{%d} must be numeric [m_i x 2].', i);
    end
end

% -------------------- compute nearest node per neuron --------------------
minDistPerNeuron = Inf(N,length(dendrites));
idxNodePerNeuron = NaN(N,length(dendrites));

hasKNN = false;
if useToolbox
    hasKNN = exist('knnsearch','file') == 2;
end

for i = 1:length(dendrites)
    Xi = dendrites{i};

    if isempty(Xi)
        % No nodes -> keep Inf distances
        continue;
    end

    if hasKNN
        % Nearest dendrite node for each spot
        [idx, d] = knnsearch(Xi, spotsXY, 'K', 1);
        idxNodePerNeuron(:,i) = idx;
        minDistPerNeuron(:,i) = d;
    else
        % Fallback: compute min distance to all nodes (memory-safe, loop spots)
        for k = 1:N
            diffs = Xi - spotsXY(k,:);
            dsq = sum(diffs.^2, 2);
            [minVal, minIdx] = min(dsq);
            minDistPerNeuron(k,i) = sqrt(minVal);
            idxNodePerNeuron(k,i) = minIdx;
        end
    end
end

% -------------------- choose winning neuron --------------------
[minDist, neuronID] = min(minDistPerNeuron, [], 2);

% Apply MaxDist threshold
unassigned = minDist > maxDist;
neuronID(unassigned) = 0;
minDist(unassigned)  = Inf;

% closest node index within assigned neuron
closestNode = NaN(N,1);
closestNodeXY = NaN(N,2);

assignedMask = neuronID > 0;
rows = find(assignedMask);
for r = rows'
    ni = neuronID(r);
    closestNode(r) = idxNodePerNeuron(r, ni);
    closestNodeXY(r,:) = dendrites{ni}(closestNode(r), :);
end

% -------------------- output --------------------
out = struct();
out.neuronID = neuronID;
out.minDist = minDist;
out.closestNode = closestNode;
out.closestNodeXY = closestNodeXY;
out.minDistPerNeuron = minDistPerNeuron;  % useful for debugging/QA
end