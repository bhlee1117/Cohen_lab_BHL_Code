function [R2_total] = get_variance(M, L)
% X: N x T matrix (N rows, T columns)
% m: 1 x T vector (length T)
% Returns the fraction of variance in X along T explained by m

% % Get dimensions
% [N, T] = size(X);
% if length(m) ~= T
%     error('Vector m must have length T matching columns of X');
% end
% validind= sum(isnan(X),1)==0 & ~isnan(m);
% X_valid=X(:,validind);
% m_valid=m(validind);
% % Compute per-row sample variances (ddof=1 equivalent)
% vars = var(X_valid, 0, 2); % Shape (N,1), variance along T (dimension 2)
% 
% % Compute per-row correlations
% cors = zeros(N, 1); % Preallocate
% for i = 1:N
%     % Correlation between row i and m
%     R = corrcoef(X_valid(i,:), m_valid);
%     cors(i) = R(1,2); % Extract corr(X_i, m)
% end
% 
% % Compute per-row R^2
% rsqs = cors .^ 2; % Shape (N,1)
% 
% % Compute weighted average of R^2, weighted by variances
% if sum(vars) == 0
%     fraction = 0; % Avoid division by zero
% else
%     prod=rsqs .* vars;
%     fraction = sum(prod(~isnan(prod))) / sum(vars(~isnan(prod)));
% end
% 
% end


% M: N x T matrix (e.g., dendritic voltage traces)
% L: 1 x T vector (predictor trace)

[N, T] = size(M);
validind= sum(isnan(M),1)==0 & ~isnan(L);
M=M(:,validind);
L=L(validind);

% Ensure L is a row vector
L = L(:)';  

% Compute beta_n for each row: beta_n = (M_n * L') / (L * L')
beta = (M * L') / (L * L');

% Predicted matrix M_hat = beta * L  (each row = beta_n * L)
M_hat = beta * L;

% Residuals
Res = M - M_hat;

% Compute variance explained per row
var_total = sum((M - mean(M,2)).^2, 2);      % N x 1
var_res   = sum(Res.^2, 2);                  % N x 1

R2_rows = 1 - (var_res ./ var_total);        % N x 1 variance explained

% Compute total variance explained across whole matrix
R2_total = 1 - (sum(var_res) / sum(var_total));
end

