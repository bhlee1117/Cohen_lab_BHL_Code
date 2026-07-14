function [y_fit, t_consts, coeff] = expfit_bnd(x, y, xs, t_guess, tau_lb, tau_ub)
% expfit_bnd  Exponential fit with NO offset and bounded time constants.
%
%   [y_fit, t_consts, coeff] = expfit_bnd(x, y, xs, t_guess, tau_lb, tau_ub)
%
% Fits   y ≈ sum_i coeff_i * exp(-x / tau_i)     (NO constant term c0)
% so the model decays to 0 as x -> Inf. The number of exponentials equals
% numel(t_guess): scalar -> single, 2-vector -> double exponential.
%
% Drop-in for expfitDM_2 (same first 4 args and 3 outputs) except:
%   * no offset (c0 = 0),
%   * the time constants are constrained to [tau_lb, tau_ub] (bounded fmincon).
%
% Inputs:
%   x, y      column data to fit
%   xs        x values at which to evaluate the fitted curve (y_fit)
%   t_guess   row vector of 1/e guesses (length = #exponentials)
%   tau_lb    lower bound(s) on tau  (scalar or per-exponential; default 0)
%   tau_ub    upper bound(s) on tau  (scalar or per-exponential; default Inf)
%
% Outputs:
%   y_fit     fitted curve evaluated at xs
%   t_consts  fitted 1/e time constants [1 x nExp]
%   coeff     exponential amplitudes    [nExp x 1]
%
% Notes:
%   * For a DECAYING fit keep tau_lb >= 0 (positive tau). Setting tau_ub <= 0
%     forces negative tau, i.e. a GROWING exponential (for rising data) — in
%     that case also set tau_lb < 0, otherwise the bounds collapse onto 0.
%
% See also: expfitDM_2, fmincon

x = x(:);  y = y(:);  xs = xs(:);
t_guess = t_guess(:)';
nExp = numel(t_guess);

if nargin < 5 || isempty(tau_lb), tau_lb = 0;   end
if nargin < 6 || isempty(tau_ub), tau_ub = Inf; end
if isscalar(tau_lb), tau_lb = repmat(tau_lb, 1, nExp); end
if isscalar(tau_ub), tau_ub = repmat(tau_ub, 1, nExp); end
t_guess = min(max(t_guess, tau_lb), tau_ub);   % keep the start point feasible

% Given tau, amplitudes are the linear least-squares solution (no offset column)
    function [sse, coef] = resid(tau)
        tau(abs(tau) < eps) = eps;              % avoid 1/0 -> NaN at x = 0
        M    = exp(-x * (1 ./ tau));            % [nSamp x nExp]
        coef = M \ y;
        sse  = sum((y - M * coef).^2);
    end

opts    = optimoptions('fmincon', 'Display', 'off');
tau_hat = fmincon(@(t) resid(t), t_guess, [], [], [], [], tau_lb, tau_ub, [], opts);

[~, coeff] = resid(tau_hat);
t_consts   = tau_hat;
t_eval     = tau_hat;  t_eval(abs(t_eval) < eps) = eps;
y_fit      = exp(-xs * (1 ./ t_eval)) * coeff;
end
