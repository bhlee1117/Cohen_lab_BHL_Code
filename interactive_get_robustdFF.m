function results = interactive_get_robustdFF(STAmov, ftprnt, F0image, opt)
% interactive_get_robustdFF - Cell-by-cell interactive dF/F fitting GUI.
%
% At startup ALL cells are auto-fitted and the full slope colormap is shown
% in panel 3.  Click any ROI in that panel to jump to it — panels 1 and 2
% update immediately.  Then use Auto / Manual to refine the selection for
% that ROI, and use Prev / Next or just keep clicking other ROIs.
%
% For each cell shows 3 panels:
%   Panel 1 (left)    F0 vs dF scatter with fit; toggle Automatic / Manual
%                     pixel selection. Title updates slope, intercept, R².
%   Panel 2 (center)  Selected pixels (px2) overlaid on F0image, zoomed to cell.
%   Panel 3 (right)   Slope map — ALL cells colored by Fslope from startup.
%                     Current cell is outlined in white.
%                     Click any cell to jump to it.
%
% Navigation buttons let you accept a cell and advance, or go back and redo.
% Clicking directly on a cell in the slope map also navigates to that cell.
%
% Usage:
%   results = interactive_get_robustdFF(STAmov, ftprnt, F0image)
%   results = interactive_get_robustdFF(STAmov, ftprnt, F0image, opt)
%
% opt fields (same as get_robustdFF):
%   .pct_low    lower percentile for auto pixel selection  (default: 2)
%   .pct_high   upper percentile for auto pixel selection  (default: 80)
%   .slope_lb   [slope_lb  intercept_lb]  for lsqlin      (default: [-Inf -Inf])
%   .slope_ub   [slope_ub  intercept_ub]  for lsqlin      (default: [ Inf  Inf])
%
% results(n) fields:
%   .Fslope     slope of the dF ~ F0 fit
%   .Intercept  intercept of the fit
%   .Rsq        R² of the fit
%   .Stds       std of residuals
%   .p          [slope intercept]
%   .px2        pixel indices used for the fit
%   .done       true once the user has accepted and moved past this cell

% -------------------------------------------------------------------------
% Default options
% -------------------------------------------------------------------------
if nargin < 4 || isempty(opt), opt = struct(); end
if ~isfield(opt, 'pct_low'),  opt.pct_low  = 2;             end
if ~isfield(opt, 'pct_high'), opt.pct_high = 80;            end
if ~isfield(opt, 'slope_lb'), opt.slope_lb = [-Inf -Inf];   end
if ~isfield(opt, 'slope_ub'), opt.slope_ub = [ Inf  Inf];   end

% -------------------------------------------------------------------------
% Vectorize inputs
% -------------------------------------------------------------------------
[nY, nX] = size(F0image);
nCells   = size(ftprnt, 3);

F0    = tovec(F0image);   % [nPix x 1]
dFall = tovec(STAmov);    % [nPix x 1]  (single STA frame assumed)

% -------------------------------------------------------------------------
% Pre-allocate results
% -------------------------------------------------------------------------
emptyResult = struct('Fslope', NaN, 'Intercept', NaN, 'Rsq', NaN, ...
                     'Stds', NaN, 'p', [NaN NaN], 'px2', [], 'done', false);
results = repmat(emptyResult, 1, nCells);

% -------------------------------------------------------------------------
% Pre-compute auto fit for ALL cells so the slope map is full at startup
% -------------------------------------------------------------------------
fprintf('Pre-computing auto fit for %d cells...\n', nCells);
for k = 1:nCells
    px_k  = find(tovec(ftprnt(:,:,k)) > 0);
    px2_k = autoPx2(px_k);
    [p_k, Rsq_k, Stds_k] = fitSelection(px2_k);
    storeResult(k, p_k, Rsq_k, Stds_k, px2_k);
end
fprintf('Done pre-computing.\n');

% -------------------------------------------------------------------------
% Build pixel-to-cell lookup for click navigation in ax3
% -------------------------------------------------------------------------
cellMap = zeros(nY, nX);
for k = 1:nCells
    cellMap(ftprnt(:,:,k) > 0) = k;
end

% -------------------------------------------------------------------------
% Build figure
% -------------------------------------------------------------------------
hFig = figure('Name', 'interactive_get_robustdFF', ...
    'NumberTitle', 'off', ...
    'Position',    [40 60 1560 680], ...
    'CloseRequestFcn', @onClose);

% Tile layout — leave bottom strip for buttons
tl = tiledlayout(hFig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
tl.OuterPosition = [0 0.13 1 0.87];

ax1 = nexttile(tl, 1);
ax2 = nexttile(tl, 2);
ax3 = nexttile(tl, 3);

% -------------------------------------------------------------------------
% Buttons  (x grows left-to-right in two logical groups)
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

hPrevBtn = uicontrol(hFig, 'Style', 'pushbutton', 'String', '< Prev', ...   %#ok<NASGU>
    'Position', [x BY BW BH], 'Callback', @(~,~) navigate(-1));
x = x + BW + G;

hNextBtn = uicontrol(hFig, 'Style', 'pushbutton', 'String', 'Next >', ...   %#ok<NASGU>
    'Position', [x BY BW BH], 'Callback', @(~,~) navigate(+1));
x = x + BW + G*4;

hFinBtn = uicontrol(hFig, 'Style', 'pushbutton', 'String', 'Finish', ...    %#ok<NASGU>
    'Position', [x BY BW BH], ...
    'BackgroundColor', [0.85 0.65 0.45], ...
    'Callback', @(~,~) cbFinish());
x = x + BW + G*3;

hLbl = uicontrol(hFig, 'Style', 'text', ...
    'Position', [x BY+2 200 BH-4], ...
    'String', sprintf('Cell 1 / %d', nCells), ...
    'FontSize', 12, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left');

% -------------------------------------------------------------------------
% Mutable state (shared across nested functions via closure)
% -------------------------------------------------------------------------
curN    = 1;
curMode = 'auto';   % 'auto' | 'manual'
isDone  = false;    % set true only by Finish button or window close

% -------------------------------------------------------------------------
% Initial draw — slope map in ax3 is already fully populated
% -------------------------------------------------------------------------
drawCell(curN);

while ~isDone && isvalid(hFig)
    pause(0.05);
end
if isvalid(hFig), delete(hFig); end

% =========================================================================
%  NESTED FUNCTIONS
% =========================================================================

    % ---------------------------------------------------------------------
    function drawCell(n)
    % Use stored result (pre-computed or last manual), refresh all panels.
        curN = n;
        hLbl.String = sprintf('Cell %d / %d', n, nCells);

        px  = find(tovec(ftprnt(:,:,n)) > 0);
        px2 = results(n).px2;
        if isempty(px2)
            px2 = autoPx2(px);
            [p_new, Rsq_new, Stds_new] = fitSelection(px2);
            storeResult(n, p_new, Rsq_new, Stds_new, px2);
        end

        refreshAx1(n, px, px2, results(n).p, results(n).Rsq);
        refreshAx2(n, px2);
        refreshAx3();
    end

    % ---------------------------------------------------------------------
    function px2 = autoPx2(px)
        px2 = px(F0(px) > prctile(F0(px), opt.pct_low) & ...
                 F0(px) < prctile(F0(px), opt.pct_high));
    end

    % ---------------------------------------------------------------------
    function [p, Rsq_val, Stds_val] = fitSelection(px2)
        if numel(px2) < 3
            p = [NaN NaN]; Rsq_val = NaN; Stds_val = NaN; return;
        end

        w  = (F0(px2)).^0.3;
        C  = [F0(px2), ones(numel(px2), 1)];
        d  = dFall(px2);
        sw = sqrt(w(:));
        Cw = sw .* C;
        dw = sw .* d;

        lb = opt.slope_lb;
        ub = opt.slope_ub;

        if all(isinf(lb) & lb < 0) && all(isinf(ub) & ub > 0)
            p_vec = Cw \ dw;
        else
            p_vec = lsqlin(Cw, dw, [], [], [], [], lb(:), ub(:), [], ...
                           optimoptions('lsqlin', 'Display', 'off'));
        end
        p = [p_vec(1), p_vec(2)];

        y_fit    = polyval(p, F0(px2));
        SS_res   = sum((dFall(px2) - y_fit).^2);
        SS_tot   = sum((dFall(px2) - mean(dFall(px2))).^2);
        Rsq_val  = 1 - SS_res / SS_tot;
        Stds_val = std(dFall(px2) - y_fit);
    end

    % ---------------------------------------------------------------------
    function storeResult(n, p, Rsq_val, Stds_val, px2)
        results(n).p         = p;
        results(n).Fslope    = p(1);
        results(n).Intercept = p(2);
        results(n).Rsq       = Rsq_val;
        results(n).Stds      = Stds_val;
        results(n).px2       = px2;
    end

    % ---------------------------------------------------------------------
    function refreshAx1(n, px, px2, p, Rsq_val)
        cla(ax1);
        hold(ax1, 'on');

        % All pixels in cell (grey)
        scatter(ax1, F0(px), dFall(px), 12, [0.76 0.76 0.76], 'filled');

        % Selected pixels (blue)
        if ~isempty(px2)
            scatter(ax1, F0(px2), dFall(px2), 12, [0.20 0.50 0.90], 'filled');
        end

        % Fit line (red)
        if ~any(isnan(p)) && ~isempty(px2)
            xrng = linspace(min(F0(px2)), max(F0(px2)), 200)';
            plot(ax1, xrng, polyval(p, xrng), 'r-', 'LineWidth', 2);
        end

        hold(ax1, 'off');
        xlabel(ax1, 'F_0');
        ylabel(ax1, '\DeltaF');
        title(ax1, sprintf('ROI #%d  |  \\DeltaF = %.3g \\cdot F_0 + %.3g  |  R^2 = %.3f', ...
            n, p(1), p(2), Rsq_val), 'FontSize', 9);
    end

    % ---------------------------------------------------------------------
    function refreshAx2(n, px2)
        cla(ax2);

        % Bounding box with padding
        [rows, cols] = find(ftprnt(:,:,n) > 0);
        pad = max(8, round(0.25 * max(range(rows), range(cols))));
        r1 = max(1, min(rows)-pad);  r2 = min(nY, max(rows)+pad);
        c1 = max(1, min(cols)-pad);  c2 = min(nX, max(cols)+pad);

        % F0 background (grey)
        imagesc(ax2, F0image(r1:r2, c1:c2));
        colormap(ax2, gray);
        axis(ax2, 'image'); axis(ax2, 'off');
        hold(ax2, 'on');

        % All cell pixels (yellow, dim)
        allMask = false(nY, nX);
        allMask(tovec(ftprnt(:,:,n)) > 0) = true;
        [ay, ax_] = find(allMask(r1:r2, c1:c2));
        if ~isempty(ay)
            scatter(ax2, ax_, ay, 10, [1 0.9 0.1], 'filled', 'MarkerFaceAlpha', 0.3);
        end

        % Selected pixels (green)
        if ~isempty(px2)
            selMask = false(nY, nX);
            selMask(px2) = true;
            [sy, sx] = find(selMask(r1:r2, c1:c2));
            if ~isempty(sy)
                scatter(ax2, sx, sy, 18, [0.15 0.82 0.25], 'filled', ...
                    'MarkerFaceAlpha', 0.70);
            end
        end

        hold(ax2, 'off');
        title(ax2, sprintf('px2 overlay  (ROI #%d)', n), 'FontSize', 9);
    end

    % ---------------------------------------------------------------------
    function refreshAx3()
        cla(ax3);

        % Build slope image from ALL cells (pre-computed or updated)
        slopeImg = nan(nY, nX);
        for k = 1:nCells
            if ~isnan(results(k).Fslope)
                linIdx = find(tovec(ftprnt(:,:,k)) > 0);
                slopeImg(linIdx) = results(k).Fslope;
            end
        end

        validIdx  = find(~isnan(slopeImg));
        bgNorm    = mat2gray(F0image);
        rgbOut    = repmat(bgNorm, [1 1 3]);

        if ~isempty(validIdx)
            allSlopes = slopeImg(validIdx);
            clim_     = [min(allSlopes) max(allSlopes)];
            if diff(clim_) < eps, clim_ = clim_ + [-1 1] * 1e-9; end

            if exist('vec2cmap', 'file')
                rgb_pts = vec2cmap(allSlopes, jet(256), clim_);
            else
                rgb_pts = localVec2cmap(allSlopes, jet(256), clim_);
            end

            for ch = 1:3
                layer           = rgbOut(:,:,ch);
                layer(validIdx) = rgb_pts(:, ch);
                rgbOut(:,:,ch)  = layer;
            end

            hImg = imagesc(ax3, rgbOut);
            colormap(ax3, jet(256));
            clim(ax3, clim_);
            cb = colorbar(ax3);
            ylabel(cb, 'slope');
        else
            hImg = imagesc(ax3, rgbOut);
            colormap(ax3, gray);
        end

        % Make the image transparent to clicks so they reach the axes
        hImg.PickableParts = 'none';

        axis(ax3, 'image'); axis(ax3, 'off');

        % Highlight the current cell with a white convex-hull outline
        hold(ax3, 'on');
        [ry, rx] = find(ftprnt(:,:,curN) > 0);
        if ~isempty(rx)
            try
                k_hull = convhull(double(rx), double(ry));
                hHull  = plot(ax3, rx(k_hull), ry(k_hull), 'w-', 'LineWidth', 2);
                hHull.PickableParts = 'none';
            catch
                hDots = scatter(ax3, rx, ry, 6, 'w', 'filled');
                hDots.PickableParts = 'none';
            end
        end
        hold(ax3, 'off');

        nDone = sum([results.done]);
        title(ax3, sprintf('Slope map  —  %d / %d accepted  |  Click ROI to jump', ...
            nDone, nCells), 'FontSize', 9);

        % Re-attach click handler (axes ButtonDownFcn survives cla,
        % but reassigning here keeps it explicit after any rebuild)
        ax3.ButtonDownFcn = @clickAx3;
    end

    % ---------------------------------------------------------------------
    function clickAx3(~, ~)
    % Map a click in ax3 to the ROI under the cursor and jump to it.
        pt  = ax3.CurrentPoint;       % [2 x 3] front/back points
        col = round(pt(1, 1));        % x = column
        row = round(pt(1, 2));        % y = row
        if row < 1 || row > nY || col < 1 || col > nX, return; end
        k = cellMap(row, col);
        if k == 0, return; end        % clicked on background

        results(curN).done = true;    % mark current cell before leaving
        navigateTo(k);
    end

    % ---------------------------------------------------------------------
    function navigateTo(k)
    % Jump to cell k and reset button state to Auto.
        if k == curN
            % Same cell — just refresh in case slope updated
            drawCell(k);
            return;
        end
        curMode = 'auto';
        hAutoBtn.Value           = 1;
        hManBtn.Value            = 0;
        hAutoBtn.BackgroundColor = [0.55 0.80 0.55];
        hManBtn.BackgroundColor  = get(0, 'defaultUicontrolBackgroundColor');
        drawCell(k);
    end

    % ---------------------------------------------------------------------
    function cbAuto(~, ~)
        curMode = 'auto';
        hAutoBtn.Value           = 1;
        hManBtn.Value            = 0;
        hAutoBtn.BackgroundColor = [0.55 0.80 0.55];
        hManBtn.BackgroundColor  = get(0, 'defaultUicontrolBackgroundColor');
        % Force recompute auto px2 for the current cell
        px  = find(tovec(ftprnt(:,:,curN)) > 0);
        px2 = autoPx2(px);
        [p, Rsq_val, Stds_val] = fitSelection(px2);
        storeResult(curN, p, Rsq_val, Stds_val, px2);
        refreshAx1(curN, px, px2, p, Rsq_val);
        refreshAx2(curN, px2);
        refreshAx3();
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
        n  = curN;
        px = find(tovec(ftprnt(:,:,n)) > 0);

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
            refreshAx1(n, px, results(n).px2, results(n).p, results(n).Rsq);
            return;
        end

        in  = inpolygon(F0(px), dFall(px), xv(:), yv(:));
        px2 = px(in);

        if numel(px2) < 3
            title(ax1, 'Too few points inside polygon — try again.', 'FontSize', 9);
            return;
        end

        [p, Rsq_val, Stds_val] = fitSelection(px2);
        storeResult(n, p, Rsq_val, Stds_val, px2);
        refreshAx1(n, px, px2, p, Rsq_val);
        refreshAx2(n, px2);
        refreshAx3();
    end

    % ---------------------------------------------------------------------
    function navigate(dir)
        results(curN).done = true;
        newN = max(1, min(nCells, curN + dir));
        navigateTo(newN);
    end

    % ---------------------------------------------------------------------
    function cbFinish()
        results(curN).done = true;
        isDone = true;
    end

    % ---------------------------------------------------------------------
    function onClose(~, ~)
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

    % ---------------------------------------------------------------------
    function rgb = localVec2cmap(v, cmap, clim_)
    % Fallback colormap lookup used when vec2cmap is not on the path.
        v      = double(v(:));
        v_norm = (v - clim_(1)) / (clim_(2) - clim_(1));
        v_norm = min(max(v_norm, 0), 1);
        idx    = round(v_norm * (size(cmap,1) - 1)) + 1;
        rgb    = cmap(idx, :);
    end

end
