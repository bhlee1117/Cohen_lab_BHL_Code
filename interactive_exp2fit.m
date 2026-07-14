function results = interactive_exp2fit(X, Y, opt)
% interactive_exp2fit - Interactive double-exponential fit GUI.
%
%   Fits Y = a1*exp(-x/b1) + a2*exp(-x/b2) to selected data points.
%   Shows a scatter of X vs Y. Use Automatic (all points) or Manual
%   (lasso-select a subset) to choose which points are used for the fit.
%   The fit curve and parameters update live. Click Finish (or close the
%   window) to return results, or Abort to cancel.
%
% Usage:
%   results = interactive_exp2fit(X, Y)
%   results = interactive_exp2fit(X, Y, opt)
%
% opt fields:
%   .a10      initial guess for a1               (default: max(Y)*0.75)
%   .b10      initial guess for b1               (default: range(X)/4)
%   .a20      initial guess for a2               (default: max(Y)*0.25)
%   .b20      initial guess for b2               (default: range(X))
%   .lb       lower bounds [a1 b1 a2 b2]         (default: [-Inf 0 -Inf 0])
%   .ub       upper bounds [a1 b1 a2 b2]         (default: [Inf Inf Inf Inf])
%   .x_fit    custom x vector for fit curve      (default: linspace over selected range)
%   .weights  weight vector (same length as X)   (default: uniform)
%
% results fields:
%   .a1, .b1, .a2, .b2   fitted parameters
%   .Rsq                  R^2 of the fit
%   .Stds                 std of residuals
%   .p                    [a1 b1 a2 b2]
%   .idxSel               indices used for the fit
%   .xFit, .yFit          fit curve (for plotting)
%   .done                 true once accepted
%   .aborted              true if Abort was clicked

% -------------------------------------------------------------------------
% Input checks / default options
% -------------------------------------------------------------------------
X = X(:);  Y = Y(:);
if numel(X) ~= numel(Y)
    error('interactive_exp2fit:sizeMismatch', 'X and Y must have the same number of elements.');
end

if nargin < 3 || isempty(opt), opt = struct(); end
if ~isfield(opt, 'a10'),     opt.a10     = max(Y) * 0.75;            end
if ~isfield(opt, 'b10'),     opt.b10     = max(range(X)/4, eps);     end
if ~isfield(opt, 'a20'),     opt.a20     = max(Y) * 0.25;            end
if ~isfield(opt, 'b20'),     opt.b20     = max(range(X), eps);       end
if ~isfield(opt, 'lb'),      opt.lb      = [-Inf, 0, -Inf, 0];       end
if ~isfield(opt, 'ub'),      opt.ub      = [ Inf, Inf,  Inf, Inf];   end
if ~isfield(opt, 'x_fit'),   opt.x_fit   = [];                       end
if ~isfield(opt, 'weights'), opt.weights = [];                        end

nPts = numel(X);

% -------------------------------------------------------------------------
% Pre-allocate result
% -------------------------------------------------------------------------
results = struct('a1', NaN, 'b1', NaN, 'a2', NaN, 'b2', NaN, ...
                 'Rsq', NaN, 'Stds', NaN, 'p', [NaN NaN NaN NaN], ...
                 'idxSel', [], 'xFit', [], 'yFit', [], ...
                 'done', false, 'aborted', false);

% -------------------------------------------------------------------------
% Build figure
% -------------------------------------------------------------------------
hFig = figure('Name', 'interactive_exp2fit', ...
    'NumberTitle', 'off', ...
    'Position',    [80 120 780 640], ...
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

uicontrol(hFig, 'Style', 'pushbutton', 'String', 'Finish', ...
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

uicontrol(hFig, 'Style', 'text', ...
    'Position', [x BY+2 280 BH-4], ...
    'String', sprintf('%d points', nPts), ...
    'FontSize', 12, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left');

% -------------------------------------------------------------------------
% Mutable state
% -------------------------------------------------------------------------
curMode = 'auto';
isDone  = false;
idxSel  = (1:nPts)';

% -------------------------------------------------------------------------
% Initial draw
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
    function [p, Rsq_val, Stds_val, xFit, yFit] = fitSelection(idx)
        if numel(idx) < 5
            p = [NaN NaN NaN NaN]; Rsq_val = NaN; Stds_val = NaN;
            xFit = []; yFit = [];
            return;
        end

        xs = X(idx);
        ys = Y(idx);

        % Weights
        if ~isempty(opt.weights)
            w    = opt.weights(idx);
            w    = w(:) / sum(w(~isnan(w)));
            sqrtW = sqrt(w);
        else
            sqrtW = ones(size(xs));
        end

        modelFun   = @(c, xx) c(1).*exp(-xx./c(2)) + c(3).*exp(-xx./c(4));
        modelFun_w = @(c, xx) sqrtW .* modelFun(c, xx);
        c0 = [opt.a10, opt.b10, opt.a20, opt.b20];

        try
            fitOpts = optimoptions('lsqcurvefit', 'Display', 'off');
            c = lsqcurvefit(modelFun_w, c0, xs, sqrtW.*ys, opt.lb, opt.ub, fitOpts);
        catch
            sseFun = @(c) sum(sqrtW.^2 .* (modelFun(c, xs) - ys).^2);
            c = fminsearch(sseFun, c0);
        end

        p = [c(1), c(2), c(3), c(4)];

        y_fit  = modelFun(c, xs);
        SS_res = sum((ys - y_fit).^2);
        SS_tot = sum((ys - mean(ys)).^2);
        Rsq_val  = 1 - SS_res / max(SS_tot, eps);
        Stds_val = std(ys - y_fit);

        if ~isempty(opt.x_fit)
            xFit = opt.x_fit(:);
        else
            xFit = linspace(min(xs), max(xs), 300)';
        end
        yFit = modelFun(c, xFit);
    end

    % ---------------------------------------------------------------------
    function storeResult(p, Rsq_val, Stds_val, idx, xFit, yFit)
        results.p      = p;
        results.a1     = p(1);  results.b1 = p(2);
        results.a2     = p(3);  results.b2 = p(4);
        results.Rsq    = Rsq_val;
        results.Stds   = Stds_val;
        results.idxSel = idx;
        results.xFit   = xFit;
        results.yFit   = yFit;
    end

    % ---------------------------------------------------------------------
    function refreshAx(idx, p, Rsq_val, xFit, yFit)
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

        % 15% margin
        xRange = range(X); yRange = range(Y);
        if xRange == 0, xRange = abs(X(1)); end
        if yRange == 0, yRange = abs(Y(1)); end
        xlim(ax1, [min(X) - 0.15*xRange, max(X) + 0.15*xRange]);
        ylim(ax1, [min(Y) - 0.15*yRange, max(Y) + 0.15*yRange]);

        xlabel(ax1, 'X');
        ylabel(ax1, 'Y');

        if ~any(isnan(p))
            title(ax1, sprintf( ...
                'Y = %.3g·exp(-X/%.3g) + %.3g·exp(-X/%.3g)  |  R^2 = %.3f  |  n = %d', ...
                p(1), p(2), p(3), p(4), Rsq_val, numel(idx)), 'FontSize', 10);
        else
            title(ax1, 'Select at least 5 points to fit', 'FontSize', 10);
        end
    end

    % ---------------------------------------------------------------------
    function doAutoSelect()
        idxSel = (1:nPts)';
        [p, Rsq_val, Stds_val, xFit, yFit] = fitSelection(idxSel);
        storeResult(p, Rsq_val, Stds_val, idxSel, xFit, yFit);
        refreshAx(idxSel, p, Rsq_val, xFit, yFit);
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
            [xv, yv] = ginputPolygon(ax1);
        end

        if numel(xv) < 3
            refreshAx(results.idxSel, results.p, results.Rsq, results.xFit, results.yFit);
            return;
        end

        in  = inpolygon(X, Y, xv(:), yv(:));
        idx = find(in);

        if numel(idx) < 5
            title(ax1, 'Too few points inside polygon (need at least 5) — try again.', 'FontSize', 9);
            return;
        end

        idxSel = idx;
        [p, Rsq_val, Stds_val, xFit, yFit] = fitSelection(idxSel);
        storeResult(p, Rsq_val, Stds_val, idxSel, xFit, yFit);
        refreshAx(idxSel, p, Rsq_val, xFit, yFit);
    end

    % ---------------------------------------------------------------------
    function cbFinish()
        results.done = true;
        isDone = true;
    end

    % ---------------------------------------------------------------------
    function cbAbort()
        results = struct('a1', NaN, 'b1', NaN, 'a2', NaN, 'b2', NaN, ...
                         'Rsq', NaN, 'Stds', NaN, 'p', [NaN NaN NaN NaN], ...
                         'idxSel', [], 'xFit', [], 'yFit', [], ...
                         'done', false, 'aborted', true);
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
                hLines(end+1) = line(ax, xv(end-1:end), yv(end-1:end), ... %#ok<AGROW>
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
