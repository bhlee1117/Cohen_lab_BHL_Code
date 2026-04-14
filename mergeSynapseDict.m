function Dict = mergeSynapseDict(tileResults, mergeRadius)
% mergeSynapseDict
%   Merge synapses from multiple aligned tiles (same H×W coordinate system)
%
% INPUT
%   tileResults{k}.S_glu        : H×W×Nsyn_k
%   tileResults{k}.synCentroid  : Nsyn_k×2  (global coordinates already)
%   mergeRadius                 : px radius for merging (e.g., 3)
%
% OUTPUT
%   Dict.S_glu      : H×W×Nsyn_merged
%   Dict.centroid   : Nsyn_merged×2
%   Dict.grp        : mapping from original synapses to merged groups

if nargin < 2
    mergeRadius = 3;
end

% Collect all synapses
S_cells = {};
C_all = [];

for k = 1:numel(tileResults)
    R = tileResults{k};
    if isempty(R) || isempty(R.S_glu), continue; end

    [H,W,Nsyn] = size(R.S_glu);

    for s = 1:Nsyn
        S_cells{end+1} = R.S_glu(:,:,s); %#ok<AGROW>
        C_all(end+1,:) = R.synCentroid(s,:); %#ok<AGROW>
    end
end

if isempty(S_cells)
    Dict = struct('S_glu', [], 'centroid', [], 'grp', []);
    return
end

S_stack = cat(3, S_cells{:});  % H×W×N0
N0 = size(S_stack,3);

%% --- Merge by centroid distance ---
r2 = mergeRadius^2;
A = false(N0,N0);

for i = 1:N0
    dx = C_all(i,1) - C_all(:,1);
    dy = C_all(i,2) - C_all(:,2);
    A(i,:) = (dx.^2 + dy.^2) <= r2;
end

G = graph(A);
grp = conncomp(G)';     % group labels
Ng = max(grp);

fprintf('Merged %d synapses into %d (radius=%g px)\n', N0, Ng, mergeRadius);

%% --- Combine footprints within each group ---
[H,W,~] = size(S_stack);
S_full = zeros(H,W,Ng,'single');
C_full = zeros(Ng,2);

for g = 1:Ng
    idx = find(grp == g);

    % average footprint
    S_full(:,:,g) = mean(S_stack(:,:,idx), 3);

    % centroid average
    C_full(g,:) = mean(C_all(idx,:), 1);
end

% Normalize each footprint (L2)
for g = 1:Ng
    v = S_full(:,:,g);
    n = sqrt(sum(v(:).^2)) + eps('single');
    S_full(:,:,g) = v / n;
end

Dict = struct();
Dict.S_glu = S_full;
Dict.centroid = C_full;
Dict.grp = grp;

end
