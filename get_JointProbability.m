function out = get_JointProbability(B, x, tau, varargin)
%get_JointProbability Joint probability of synchronous activity with ±tau tolerance.
%
% out = get_JointProbability(B, x, tau)
% out = get_JointProbability(..., 'Method','anybin')
%
% INPUTS
%   B   : [N x T] binary (logical or 0/1), synapse x time
%   x   : tuple size (2 for pairs, 3 for triples, ...)
%   tau : nonnegative integer, temporal tolerance in frames (±tau)
%
% Name-Value
%   'MaxComb' : cap number of combinations (default Inf). Useful when N is large.
%   'NeuronIdx' : vector of neuron indices to include (default 1:N)
%   'ReturnDilated' : true/false (default false) include dilated matrix in output
%
% OUTPUT (struct)
%   out.combos        : [C x x] indices of synapses in each tuple
%   out.P_joint       : [C x 1] joint probability with tolerance
%   out.P_ind         : [C x 1] independent estimate (product of marginals with tolerance)
%   out.P_marg_tau    : [M x 1] marginal probabilities after dilation (for included neurons)
%   out.tau           : tau used
%   out.x             : tuple size
%
% Definition:
%   Let D be B dilated in time by ±tau (i.e., D(i,t)=1 if B(i, t±tau window) has any 1).
%   Then P_joint(combo) = mean( all(D(combo,:), 1) )
%        P_ind(combo)   = prod( mean(D(i,:)) for i in combo )
%
% Notes:
%   - This counts synchrony across ANY time bin; it's a "co-occupancy" probability.
%   - If you want an event-triggered definition (e.g., around events of one synapse),
%     tell me and I'll adapt it.

% ------------------ parse ------------------
p = inputParser;
p.addParameter('MaxComb', Inf, @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.addParameter('NeuronIdx', [], @(v)isnumeric(v)&&isvector(v));
p.addParameter('ReturnDilated', false, @(v)islogical(v)&&isscalar(v));
p.parse(varargin{:});

maxComb = p.Results.MaxComb;
idx = p.Results.NeuronIdx;
retDil = p.Results.ReturnDilated;

B = logical(B);
[N,T] = size(B);

if isempty(idx), idx = 1:N; end
idx = idx(:)';
M = numel(idx);

if x < 2 || x > M
    error('x must be between 2 and number of selected neurons (%d).', M);
end
tau = max(0, round(tau));

% ------------------ dilate in time by ±tau ------------------
if tau == 0
    D = B(idx,:);
else
    % dilation by convolution along time dimension
    kernel = ones(1, 2*tau + 1);
    D = conv2(double(B(idx,:)), kernel, 'same') > 0;  % [M x T] logical
end

% marginal probability under same tolerance
P_marg_tau = mean(D, 2);  % [M x 1]

% ------------------ enumerate combos ------------------
combos = nchoosek(1:M, x);       % indices into idx
C = size(combos,1);
if isfinite(maxComb) && C > maxComb
    % random subsample combos if too many
    rp = randperm(C, maxComb);
    combos = combos(rp,:);
    C = size(combos,1);
end

% map back to original neuron indices
fprintf('Number of combination: %2.0f\n',C);
combosNeuron = idx(combos);

% ------------------ compute joint probability ------------------
P_joint = zeros(C,1);
P_ind   = zeros(C,1);

for c = 1:C
    rows = combos(c,:);
    jointVec = all(D(rows,:), 1);      % 1 x T
    P_joint(c) = mean(jointVec);

    P_ind(c) = prod(P_marg_tau(rows));
end

% ------------------ output ------------------
out = struct();
out.combos      = combosNeuron;     % in original neuron IDs
out.P_joint     = P_joint;
out.P_ind       = P_ind;
out.P_marg_tau  = P_marg_tau;
out.idxIncluded = idx;
out.tau         = tau;
out.x           = x;

if retDil
    out.D = D;
end
end