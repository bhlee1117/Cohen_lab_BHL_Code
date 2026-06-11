%% validate_nanXCorr.m
% Self-contained validation script for nanXCorr.
% Each test prints PASS or FAIL with the reason.
% Run this script after any changes to nanXCorr.m.
%
% Tests covered:
%   1. Known-lag sinusoid: xcorr peak should be at the true lag
%   2. Identical signals: zero-lag correlation should be 1.0
%   3. Anti-correlated signals: zero-lag correlation should be -1.0
%   4. NaN handling: NaN frames should be ignored, not corrupt the result
%   5. Unnormalized output: compare against manual sum
%   6. Zero-variance signal: should return NaN, not error
%   7. Single valid pair: should return NaN, not error
%   8. Symmetry: xcorr(x,y,lag) == xcorr(y,x,-lag)
%   9. Staggered NaNs: x and y have NaNs at different time points; only
%      overlapping valid pairs should be used at each lag

fprintf('\n===== nanXCorr validation =====\n\n');
tol  = 1e-6;   % floating point tolerance for equality checks
nPts = 500;    % signal length
rng(42);       % reproducible random numbers

pass = 0;
fail = 0;

%% ---- Test 1: Known-lag sinusoid ----
% x is a sine wave; y is x shifted by true_lag samples.
% The xcorr peak should land exactly at true_lag.
true_lag = 20;
t = (0:nPts-1)';
x1 = sin(2*pi*t/100);
y1 = sin(2*pi*(t - true_lag)/100);

[cc1, lags1] = nanXCorr(x1, y1, 60, true);
[~, peak_idx] = max(cc1);
detected_lag  = lags1(peak_idx);

if detected_lag == true_lag
    fprintf('[PASS] Test 1: Known-lag sinusoid      (detected lag = %d)\n', detected_lag);
    pass = pass + 1;
else
    fprintf('[FAIL] Test 1: Known-lag sinusoid      (expected %d, got %d)\n', true_lag, detected_lag);
    fail = fail + 1;
end

%% ---- Test 2: Identical signals -> zero-lag corr == 1 ----
x2 = randn(nPts, 1);
y2 = x2;

[cc2, lags2] = nanXCorr(x2, y2, 10, true);
zero_lag_corr = cc2(lags2 == 0);

if abs(zero_lag_corr - 1.0) < tol
    fprintf('[PASS] Test 2: Identical signals       (corr at lag-0 = %.6f)\n', zero_lag_corr);
    pass = pass + 1;
else
    fprintf('[FAIL] Test 2: Identical signals       (expected 1.0, got %.6f)\n', zero_lag_corr);
    fail = fail + 1;
end

%% ---- Test 3: Anti-correlated signals -> zero-lag corr == -1 ----
x3 = randn(nPts, 1);
y3 = -x3;

[cc3, lags3] = nanXCorr(x3, y3, 10, true);
zero_lag_corr3 = cc3(lags3 == 0);

if abs(zero_lag_corr3 - (-1.0)) < tol
    fprintf('[PASS] Test 3: Anti-correlated         (corr at lag-0 = %.6f)\n', zero_lag_corr3);
    pass = pass + 1;
else
    fprintf('[FAIL] Test 3: Anti-correlated         (expected -1.0, got %.6f)\n', zero_lag_corr3);
    fail = fail + 1;
end

%% ---- Test 4: NaN handling (same positions in x and y) ----
% Insert NaNs at the same frames in both signals.
% The detected peak should still be at true_lag.
nan_lag  = 15;
x4 = sin(2*pi*t/100);
y4 = sin(2*pi*(t - nan_lag)/100);

nan_mask = false(nPts, 1);
nan_mask(randperm(nPts, round(nPts*0.2))) = true;   % 20% NaNs, same positions
x4(nan_mask) = NaN;
y4(nan_mask) = NaN;

[cc4, lags4] = nanXCorr(x4, y4, 40, true);
[~, peak_idx4] = max(cc4);
detected_lag4  = lags4(peak_idx4);

if detected_lag4 == nan_lag
    fprintf('[PASS] Test 4: NaN same positions      (detected lag = %d with 20%% NaNs)\n', detected_lag4);
    pass = pass + 1;
else
    fprintf('[FAIL] Test 4: NaN same positions      (expected %d, got %d)\n', nan_lag, detected_lag4);
    fail = fail + 1;
end

%% ---- Test 5: Unnormalized output vs manual sum ----
% For lag=0 and no NaNs, unnormalized xcorr should equal sum(x.*y)/N
x5 = randn(nPts, 1);
y5 = randn(nPts, 1);

[cc5, lags5] = nanXCorr(x5, y5, 5, false);
manual_lag0   = sum(x5 .* y5) / nPts;
func_lag0     = cc5(lags5 == 0);

if abs(func_lag0 - manual_lag0) < tol
    fprintf('[PASS] Test 5: Unnormalized lag-0      (%.6f vs manual %.6f)\n', func_lag0, manual_lag0);
    pass = pass + 1;
else
    fprintf('[FAIL] Test 5: Unnormalized lag-0      (got %.6f, expected %.6f)\n', func_lag0, manual_lag0);
    fail = fail + 1;
end

%% ---- Test 6: Zero-variance signal -> NaN, no error ----
x6 = ones(nPts, 1);   % flat signal, zero variance
y6 = randn(nPts, 1);

try
    [cc6, lags6] = nanXCorr(x6, y6, 5, true);
    zero_lag6 = cc6(lags6 == 0);
    if isnan(zero_lag6)
        fprintf('[PASS] Test 6: Zero-variance           (correctly returned NaN)\n');
        pass = pass + 1;
    else
        fprintf('[FAIL] Test 6: Zero-variance           (expected NaN, got %.6f)\n', zero_lag6);
        fail = fail + 1;
    end
catch ME
    fprintf('[FAIL] Test 6: Zero-variance           (threw error: %s)\n', ME.message);
    fail = fail + 1;
end

%% ---- Test 7: Single valid pair after NaN removal -> NaN, no error ----
x7        = randn(nPts, 1);
y7        = randn(nPts, 1);
x7(2:end) = NaN;   % only index 1 is valid
y7(2:end) = NaN;

try
    [cc7, lags7] = nanXCorr(x7, y7, 5, true);
    zero_lag7 = cc7(lags7 == 0);
    if isnan(zero_lag7)
        fprintf('[PASS] Test 7: Single valid pair       (correctly returned NaN)\n');
        pass = pass + 1;
    else
        fprintf('[FAIL] Test 7: Single valid pair       (expected NaN, got %.6f)\n', zero_lag7);
        fail = fail + 1;
    end
catch ME
    fprintf('[FAIL] Test 7: Single valid pair       (threw error: %s)\n', ME.message);
    fail = fail + 1;
end

%% ---- Test 8: Symmetry xcorr(x,y,lag) == xcorr(y,x,-lag) ----
x8 = randn(nPts, 1);
y8 = randn(nPts, 1);
maxLag8 = 30;

[cc8_xy, lags8] = nanXCorr(x8, y8, maxLag8, true);
[cc8_yx, ~]     = nanXCorr(y8, x8, maxLag8, true);

% xcorr(x,y) at lag L == xcorr(y,x) at lag -L
cc8_yx_flipped = fliplr(cc8_yx);
sym_err = max(abs(cc8_xy - cc8_yx_flipped));

if sym_err < tol
    fprintf('[PASS] Test 8: Symmetry                (max asymmetry = %.2e)\n', sym_err);
    pass = pass + 1;
else
    fprintf('[FAIL] Test 8: Symmetry                (max asymmetry = %.2e, tol = %.2e)\n', sym_err, tol);
    fail = fail + 1;
end

%% ---- Test 9: Staggered NaNs (x and y NaN at different time points) ----

% x and y have NaN masks that do NOT overlap — every frame is valid in at
% least one signal but may be missing in the other.
%
% At each lag, the function must use ONLY frames where BOTH the lagged x
% and the lagged y are simultaneously valid. The ground truth is computed
% by manually applying the same pairwise-valid mask for lag=0.
%
% Sub-test A: lag=0, unnormalized — result must match manual pair-wise sum.
% Sub-test B: peak detection — staggered NaNs should not shift the peak.

nan_lag9   = 10;
x9_clean   = sin(2*pi*t/80);
y9_clean   = sin(2*pi*(t - nan_lag9)/80);

% Build two NON-OVERLAPPING NaN masks (each covers ~20% of frames)
all_idx    = randperm(nPts);
nan_x_idx  = all_idx(1             : round(nPts*0.2));       % first 20%
nan_y_idx  = all_idx(round(nPts*0.2)+1 : round(nPts*0.4));  % next 20%

x9 = x9_clean;
y9 = y9_clean;
x9(nan_x_idx) = NaN;
y9(nan_y_idx) = NaN;

% Sub-test A: unnormalized lag-0, manual vs function
valid9      = ~isnan(x9) & ~isnan(y9);
n_valid9    = sum(valid9);
manual9_lag0 = sum(x9(valid9) .* y9(valid9)) / n_valid9;

[cc9, lags9] = nanXCorr(x9, y9, 30, false);
func9_lag0   = cc9(lags9 == 0);

subtestA_pass = abs(func9_lag0 - manual9_lag0) < tol;

% Sub-test B: peak should still be near nan_lag9
[~, peak_idx9] = max(nanXCorr(x9, y9, 30, true));
[~, lags9_tmp] = nanXCorr(x9, y9, 30, true);   % get lags vector
detected_lag9  = lags9_tmp(peak_idx9);
subtestB_pass  = detected_lag9 == nan_lag9;

if subtestA_pass && subtestB_pass
    fprintf('[PASS] Test 9: Staggered NaNs          (lag-0 unnorm=%.6f=manual, peak at %d)\n', ...
        func9_lag0, detected_lag9);
    pass = pass + 1;
else
    if ~subtestA_pass
        fprintf('[FAIL] Test 9A: Staggered NaNs unnorm  (got %.6f, expected %.6f)\n', ...
            func9_lag0, manual9_lag0);
    end
    if ~subtestB_pass
        fprintf('[FAIL] Test 9B: Staggered NaNs peak    (expected %d, got %d)\n', ...
            nan_lag9, detected_lag9);
    end
    fail = fail + 1;
end

%% ---- Summary ----

fprintf('\n===== Results: %d / %d passed =====\n\n', pass, pass+fail);



%% ---- Plot: Test 1, Test 4, Test 9 for visual inspection ----
figure('Name', 'nanXCorr validation', 'Position', [100 100 1100 350]);

subplot(1,3,1);
plot(lags1, cc1, 'b-', 'LineWidth', 1.5); hold on;
xline(true_lag, 'r--', 'LineWidth', 1.5);
xlabel('Lag (samples)'); ylabel('Correlation');
title(sprintf('Test 1: Known lag=%d (peak at %d)', true_lag, detected_lag));
legend('xcorr', 'true lag'); grid on;

subplot(1,3,2);
plot(lags4, cc4, 'b-', 'LineWidth', 1.5); hold on;
xline(nan_lag, 'r--', 'LineWidth', 1.5);
xlabel('Lag (samples)'); ylabel('Correlation');
title(sprintf('Test 4: Same NaNs 20%%, lag=%d (peak at %d)', nan_lag, detected_lag4));
legend('xcorr', 'true lag'); grid on;

subplot(1,3,3);
[cc9_norm, lags9_norm] = nanXCorr(x9, y9, 30, true);
plot(lags9_norm, cc9_norm, 'b-', 'LineWidth', 1.5); hold on;
xline(nan_lag9, 'r--', 'LineWidth', 1.5);
xlabel('Lag (samples)'); ylabel('Correlation');
title(sprintf('Test 9: Staggered NaNs, lag=%d (peak at %d)', nan_lag9, detected_lag9));
legend('xcorr', 'true lag'); grid on;
