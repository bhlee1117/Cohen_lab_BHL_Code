function results = interactive_expfit(X, Y, opt)
% interactive_expfit - Interactive exponential fit GUI (single or double exp).
%
% Fits Y = sum_i c_i*exp(-X/tau_i) (NO offset, c0 = 0) with expfit_bnd (same
% fitter as the batch analysis, so interactive refinement stays consistent).
% The number of exponentials is set by opt.t_guess: a scalar -> single
% exponential, a 2-element vector -> double exponential. Time constants are
% constrained to [opt.tau_lb, opt.tau_ub].
%
% Shows a scatter of X vs Y. Use Automatic (all points) or Manual (lasso-
% select a subset) to choose which points are used for the fit. The fit curve
% and parameters update live. Click Finish (or close the window) to return.
%
% Usage:
%   results = interactive_expfit(X, Y)
%   results = interactive_expfit(X, Y, opt)
%
% X, Y      vectors of equal length (data to fit)
%
% opt fields:
%   .t_guess        1/e guess(es); scalar -> single exp, [g1 g2] -> double exp
%                                                     (default: range(X)/2)
%   .tau_lb         lower bound(s) on tau              (default: 0)
%   .tau_ub         upper bound(s) on tau              (default: Inf)
%   .x_fit          x values at which to sample the fit curve (default: dense)
%   .predefined_idx indices into X/Y to use for the fit; skips interactive GUI
%
% results fields:
%   .a          amplitude of the slowest exponential (back-compat)
%   .b          slowest 1/e decay constant           (back-compat)
%   .tau        all fitted 1/e constants  [1 x nExp]
%   .amp        all exponential amplitudes [1 x nExp]
%   .offset     constant term c0
%   .Rsq        R^2 of the fit
%   .Stds       std of residuals
%   .p          [a b]
%   .idxSel     indices (into X/Y) used for the fit
%   .xFit,.yFit fit curve (dense, for plotting)
%   .done, .aborted
%
% See also: expfitDM_2

% -------------------------------------------------------------------------
% Input checks / default options
% -------------------------------------------------------------------------
X = X(:);
Y = Y(:);
if numel(X) ~= numel(Y)
    error('interactive_expfit:sizeMismatch', 'X and Y must have the same number of elements.');
end

if nargin < 3 || isempty(opt), opt = struct(); end
if ~isfield(opt, 't_guess') || isempty(opt.t_guess), opt.t_guess = max(range(X)/2, eps); end
if ~isfield(opt, 'tau_lb') || isempty(opt.tau_lb), opt.tau_lb = 0;   end
if ~isfield(opt, 'tau_ub') || isempty(opt.tau_ub), opt.tau_ub = Inf; end
if ~isfield(opt, 'x_fit'),          opt.x_fit          = [];    end
if ~isfield(opt, 'predefined_idx'), opt.predefined_idx = [];    end

nPts  = numel(X);
nExp  = numel(opt.t_guess);

% -------------------------------------------------------------------------
% Pre-allocate result
% -------------------------------------------------------------------------
results = struct('a', NaN, 'b', NaN, 'tau', nan(1,nExp), 'amp', nan(1,nExp), ...
                 'offset', NaN, 'Rsq', NaN, 'Stds', NaN, 'p', [NaN NaN], ...
                 'idxSel', [], 'xFit', [], 'yFit', [], 'done', false, 'aborted', false);

% -------------------------------------------------------------------------
% Predefined-index mode: skip GUI entirely
% -------------------------------------------------------------------------
if ~isempty(opt.predefined_idx)
    idxSel = opt.predefined_idx(:);
    [p, Rsq_val, Stds_val, xFit, yFit, tcon, amp, offset] = fitSelection(idxSel);
    storeResult(p, Rsq_val, Stds_val, idxSel, xFit, yFit, tcon, amp, offset);
    results.done = true;
    return;
end

% -------------------------------------------------------------------------
% Build figure
% -------------------------------------------------------------------------
hFig = figure('Name', 'interactive_expfit', ...
    'NumberTitle', 'off', ...
    'Position',    [80 120 760 640], ...
    'CloseRequestFcn', @onClose);

tl = tiledlayout(hFig, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
tl.OuterPosition = [0 0.13 1 0.87];

ax1 = nexttile(tl, 1);

% -------------------------------------------------------------------------
% Buttons
% -------------------------------------------------------------------------
BW = 108;  BH = 32;  G = 8;  BY = 14;
x  = G;

hAutoBtn = uicontrol(hFig, 'Style', 'togglebutton', 'String', 'Automatic', ...
    'Value', 1, 'Position', [x BY BW BH], ...
    'BackgroundColor', [0.55 0.80 0.55], ...
    'Callback', @cbAuto);
x = x + BW + G;

hManBtn = uicontrol(hFig, 'Style', 'togglebutton', 'String', 'Manual', ...
    'Value', 0, 'Position', [x BY BW BH], ...
    'Callback', @cbManual);
x = x + BW + G*4;

hFinBtn = uicontrol(hFig, 'Style', 'pushbutton', 'String', 'Finish', ...    %#ok<NASGU>
    'Position', [x BY BW BH], ...
    'BackgroundColor', [0.85 0.65 0.45], ...
    'Callback', @(~,~) cbFinish());
x = x + BW + G;

uicontrol(hFig, 'Style', 'pushbutton', 'String', 'Abort', ...
    'Position', [x BY BW BH], ...
    'BackgroundColor', [0.90 0.40 0.40], ...
    'ForegroundColor', [1 1 1], ...
    'FontWeight', 'bold', ...
    'Callback', @(~,~) cbAbort());
x = x + BW + G*2;

hLbl = uicontrol(hFig, 'Style', 'text', ...
    'Position', [x BY+2 260 BH-4], ...
    'String', sprintf('%d points | %d-exp', nPts, nExp), ...
    'FontSize', 12, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left');                                          %#ok<NASGU>

% -------------------------------------------------------------------------
% Mutable state (shared across nested functions via closure)
% -------------------------------------------------------------------------
curMode = 'auto';   % 'auto' | 'manual'
isDone  = false;    % set true only by Finish button or window close
idxSel  = (1:nPts)';

% -------------------------------------------------------------------------
% Initial draw — auto mode, all points selected
% -------------------------------------------------------------------------
doAutoSelect();

while ~isDone && isvalid(hFig)
    pause(0.05);
end
if isvalid(hFig), delete(hFig); end

% =========================================================================
%  NESTED FUNCTIONS
% =========================================================================

    % ---------------------------------------------------------------------
    function [p, Rsq_val, Stds_val, xFit, yFit, tcon, amp, offset] = fitSelection(idx)
    % Fit Y = sum c_i*exp(-X/tau_i) (no offset) on the selected indices via expfit_bnd.
        if numel(idx) < 3
            p = [NaN NaN]; Rsq_val = NaN; Stds_val = NaN; xFit = []; yFit = [];
            tcon = nan(1,nExp); amp = nan(1,nExp); offset = 0;
            return;
        end

        xs = X(idx);
        ys = Y(idx);
        tg = opt.t_guess(:)';                 % row of guesses (1 or 2 -> single/double)

        % Fit on the selected points, evaluate at the same points for R^2
        [yhat_sel, tcon, coeff] = expfit_bnd(xs, ys, xs, tg, opt.tau_lb, opt.tau_ub);
        offset = 0;                           % no offset in this model (c0 = 0)
        amp    = coeff(:)';                   % 1 x nExp exponential amplitudes

        SS_res = sum((ys - yhat_sel).^2);
        SS_tot = sum((ys - mean(ys)).^2);
        if SS_tot > 0, Rsq_val = 1 - SS_res / SS_tot; else, Rsq_val = NaN; end
        Stds_val = std(ys - yhat_sel);

        % Dense fit curve (re-uses the same fitter)
        if ~isempty(opt.x_fit)
            xFit = opt.x_fit(:);
        else
            xFit = linspace(min(xs), max(xs), 200)';
        end
        yFit = expfit_bnd(xs, ys, xFit, tg, opt.tau_lb, opt.tau_ub);

        % Back-compat scalars: slowest exponential component
        [~, islow] = max(tcon);
        p = [amp(islow), tcon(islow)];
    end

    % ---------------------------------------------------------------------
    function storeResult(p, Rsq_val, Stds_val, idx, xFit, yFit, tcon, amp, offset)
        results.p      = p;
        results.a      = p(1);
        results.b      = p(2);
        results.tau    = tcon;
        results.amp    = amp;
        results.offset = offset;
        results.Rsq    = Rsq_val;
        results.Stds   = Stds_val;
        results.idxSel = idx;
        results.xFit   = xFit;
        results.yFit   = yFit;
    end

    % ---------------------------------------------------------------------
    function refreshAx(idx, p, Rsq_val, xFit, yFit, tcon)
        cla(ax1);
        hold(ax1, 'on');

        % All points (grey)
        scatter(ax1, X, Y, 16, [0.76 0.76 0.76], 'filled');

        % Selected points (blue)
        if ~isempty(idx)
            scatter(ax1, X(idx), Y(idx), 16, [0.20 0.50 0.90], 'filled');
        end

        % Fit curve (red)
        if ~any(isnan(p)) && ~isempty(xFit)
            plot(ax1, xFit, yFit, 'r-', 'LineWidth', 2);
        end

        hold(ax1, 'off');

        % Add 15% margin around data extents
        xRange = range(X); yRange = range(Y);
        if xRange == 0, xRange = abs(X(1)); end
        if yRange == 0, yRange = abs(Y(1)); end
        xlim(ax1, [min(X) - 0.15*xRange, max(X) + 0.15*xRange]);
        ylim(ax1, [min(Y) - 0.15*yRange, max(Y) + 0.15*yRange]);

        xlabel(ax1, 'X');
        ylabel(ax1, 'Y');
        if ~any(isnan(p))
            tauStr = strjoin(compose('%.3g', tcon(:)'), ', ');
            title(ax1, sprintf('\\tau = [%s]  |  R^2 = %.3f  |  n = %d', ...
                tauStr, Rsq_val, numel(idx)), 'FontSize', 10);
        else
            title(ax1, 'Select at least 3 points to fit', 'FontSize', 10);
        end
    end

    % ---------------------------------------------------------------------
    function doAutoSelect()
        idxSel = (1:nPts)';
        [p, Rsq_val, Stds_val, xFit, yFit, tcon, amp, offset] = fitSelection(idxSel);
        storeResult(p, Rsq_val, Stds_val, idxSel, xFit, yFit, tcon, amp, offset);
        refreshAx(idxSel, p, Rsq_val, xFit, yFit, tcon);
    end

    % ---------------------------------------------------------------------
    function cbAuto(~, ~)
        curMode = 'auto';                                                    %#ok<NASGU>
        hAutoBtn.Value           = 1;
        hManBtn.Value            = 0;
        hAutoBtn.BackgroundColor = [0.55 0.80 0.55];
        hManBtn.BackgroundColor  = get(0, 'defaultUicontrolBackgroundColor');
        doAutoSelect();
    end

    % ---------------------------------------------------------------------
    function cbManual(~, ~)
        curMode = 'manual';                                                  %#ok<NASGU>
        hAutoBtn.Value           = 0;
        hManBtn.Value            = 1;
        hAutoBtn.BackgroundColor = get(0, 'defaultUicontrolBackgroundColor');
        hManBtn.BackgroundColor  = [0.55 0.80 0.55];
        doManualSelect();
    end

    % ---------------------------------------------------------------------
    function doManualSelect()
        title(ax1, 'Draw a lasso around points to select. Double-click or right-click to confirm.', ...
            'FontSize', 9);
        drawnow;

        % --- try drawpolygon (Image Processing Toolbox R2019b+) ----------
        xv = [];  yv = [];
        try
            roi = drawpolygon(ax1);
            wait(roi);
            if isvalid(roi) && ~isempty(roi.Position)
                xv = roi.Position(:,1);
                yv = roi.Position(:,2);
                delete(roi);
            end
        catch
            % Fallback: click-to-build polygon, right-click to close
            [xv, yv] = ginputPolygon(ax1);
        end

        if numel(xv) < 3
            % Cancelled — redisplay current stored state without changes
            refreshAx(results.idxSel, results.p, results.Rsq, results.xFit, results.yFit, results.tau);
            return;
        end

        in  = inpolygon(X, Y, xv(:), yv(:));
        idx = find(in);

        if numel(idx) < 3
            title(ax1, 'Too few points inside polygon — try again.', 'FontSize', 9);
            return;
        end

        idxSel = idx;
        [p, Rsq_val, Stds_val, xFit, yFit, tcon, amp, offset] = fitSelection(idxSel);
        storeResult(p, Rsq_val, Stds_val, idxSel, xFit, yFit, tcon, amp, offset);
        refreshAx(idxSel, p, Rsq_val, xFit, yFit, tcon);
    end

    % ---------------------------------------------------------------------
    function cbFinish()
        results.done = true;
        isDone = true;
    end

    % ---------------------------------------------------------------------
    function cbAbort()
        results = struct('a', NaN, 'b', NaN, 'tau', nan(1,nExp), 'amp', nan(1,nExp), ...
                         'offset', NaN, 'Rsq', NaN, 'Stds', NaN, 'p', [NaN NaN], ...
                         'idxSel', [], 'xFit', [], 'yFit', [], 'done', false, 'aborted', true);
        isDone = true;
    end

    % ---------------------------------------------------------------------
    function onClose(~, ~)
        results.done = true;
        isDone = true;
        delete(hFig);
    end

    % ---------------------------------------------------------------------
    function [xv, yv] = ginputPolygon(ax)
    % Click to add vertices; right-click (btn==3) or Enter to close polygon.
        xv = [];  yv = [];  hLines = gobjects(0);
        axes(ax);                                                            %#ok<LAXES>
        title(ax, 'Click to add vertices. Right-click to finish.', 'FontSize', 9);

        while true
            try
                [xi, yi, btn] = ginput(1);
            catch
                break;
            end
            if isempty(xi) || btn == 3, break; end
            xv(end+1) = xi;   yv(end+1) = yi;                              %#ok<AGROW>
            if numel(xv) > 1
                hLines(end+1) = line(ax, xv(end-1:end), yv(end-1:end), ...  %#ok<AGROW>
                    'Color', [0.9 0.2 0.2], 'LineWidth', 1.5);
            end
        end

        if numel(xv) > 2
            hLines(end+1) = line(ax, [xv(end) xv(1)], [yv(end) yv(1)], ... %#ok<AGROW,NASGU>
                'Color', [0.9 0.2 0.2], 'LineWidth', 1.5);
        end
        drawnow;
    end

end
