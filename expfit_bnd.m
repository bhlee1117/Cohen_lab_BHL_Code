function [y_fit, t_consts, coeff, c0] = expfit_bnd(x, y, xs, t_guess, tau_lb, tau_ub, use_offset)
% expfit_bnd  Exponential fit with optional offset and bounded time constants.
%
%   [y_fit, t_consts, coeff, c0] = expfit_bnd(x, y, xs, t_guess, tau_lb, tau_ub, use_offset)
%
% Fits   y ≈ sum_i coeff_i * exp(-x / tau_i) + c0
% The number of exponentials equals numel(t_guess): scalar -> single,
% 2-vector -> double exponential.
%
% use_offset (default FALSE) controls c0:
%   false - c0 is forced to 0, so the model decays to 0 as x -> Inf (the
%           original behaviour; omit the argument and nothing changes).
%   true  - c0 is a free parameter, i.e. y = a*exp(-x/b) + c. Use this when the
%           data plateaus at a non-zero level instead of decaying to baseline.
% c0 is recovered by separable least squares exactly like the amplitudes: for a
% given tau the model is linear in [coeff; c0], so adding the offset is one extra
% column of ones in the design matrix and costs no extra nonlinear parameter.
%
% Drop-in for expfitDM_2 (same first 4 args and 3 outputs) except:
%   * no offset by default (c0 = 0 unless use_offset is set),
%   * the time constants are constrained to [tau_lb, tau_ub] (bounded fmincon).
%
% Inputs:
%   x, y      column data to fit
%   xs        x values at which to evaluate the fitted curve (y_fit)
%   t_guess   row vector of 1/e guesses (length = #exponentials)
%   tau_lb    lower bound(s) on tau  (scalar or per-exponential; default 0)
%   tau_ub    upper bound(s) on tau  (scalar or per-exponential; default Inf)
%   use_offset  fit a free constant term c0 (default false)
%
% Outputs:
%   y_fit     fitted curve evaluated at xs (offset included)
%   t_consts  fitted 1/e time constants [1 x nExp]
%   coeff     exponential amplitudes    [nExp x 1]
%   c0        fitted constant term (0 when use_offset is false)
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
if nargin < 7 || isempty(use_offset), use_offset = false; end
if isscalar(tau_lb), tau_lb = repmat(tau_lb, 1, nExp); end
if isscalar(tau_ub), tau_ub = repmat(tau_ub, 1, nExp); end
t_guess = min(max(t_guess, tau_lb), tau_ub);   % keep the start point feasible

% Given tau, amplitudes (and c0) are the linear least-squares solution
    function [sse, coef] = resid(tau)
        tau(abs(tau) < eps) = eps;              % avoid 1/0 -> NaN at x = 0
        M    = exp(-x * (1 ./ tau));            % [nSamp x nExp]
        if use_offset, M = [M, ones(numel(x),1)]; end   % extra column -> free c0
        coef = M \ y;
        sse  = sum((y - M * coef).^2);
    end

opts    = optimoptions('fmincon', 'Display', 'off');
tau_hat = fmincon(@(t) resid(t), t_guess, [], [], [], [], tau_lb, tau_ub, [], opts);

[~, coef_all] = resid(tau_hat);
if use_offset
    coeff = coef_all(1:nExp);  c0 = coef_all(end);
else
    coeff = coef_all;          c0 = 0;
end
t_consts   = tau_hat;
t_eval     = tau_hat;  t_eval(abs(t_eval) < eps) = eps;
y_fit      = exp(-xs * (1 ./ t_eval)) * coeff + c0;
end
