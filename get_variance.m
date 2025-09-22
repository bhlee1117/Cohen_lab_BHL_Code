function [fraction] = get_variance(X, m)
% X: N x T matrix (N rows, T columns)
% m: 1 x T vector (length T)
% Returns the fraction of variance in X along T explained by m

% Get dimensions
[N, T] = size(X);
if length(m) ~= T
    error('Vector m must have length T matching columns of X');
end
validind= sum(isnan(X),1)==0 & ~isnan(m);
X_valid=X(:,validind);
m_valid=m(validind);
% Compute per-row sample variances (ddof=1 equivalent)
vars = var(X_valid, 0, 2); % Shape (N,1), variance along T (dimension 2)

% Compute per-row correlations
cors = zeros(N, 1); % Preallocate
for i = 1:N
    % Correlation between row i and m
    R = corrcoef(X_valid(i,:), m_valid);
    cors(i) = R(1,2); % Extract corr(X_i, m)
end

% Compute per-row R^2
rsqs = cors .^ 2; % Shape (N,1)

% Compute weighted average of R^2, weighted by variances
if sum(vars) == 0
    fraction = 0; % Avoid division by zero
else
    fraction = sum(rsqs .* vars) / sum(vars);
end

end