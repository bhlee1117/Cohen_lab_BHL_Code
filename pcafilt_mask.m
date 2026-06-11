function [moviefilt, eigvecs, eigvals] = pcafilt_mask(movie, mask)
% function [moviefilt, eigvecs, eigvals] = pcafilt_mask(movie, mask)
% PCA denoising restricted to masked pixels.
% movie : X x Y x T
% mask  : X x Y logical (or numeric, thresholded at >0)
%
% Shows all PCs (10 per page) with eigentraces + eigenimages in one window.
% Use Prev / Next buttons to navigate; click Done to proceed to PC selection.
% Unmasked pixels are returned unchanged.

[ysize, xsize, nframes] = size(movie);
mask = logical(mask);

% --- extract masked pixels: [npix x T] ---
mov2d      = reshape(movie, [ysize*xsize, nframes]);
mov2d_mask = mov2d(mask(:), :);          % [npix x T]

avgvec = mean(mov2d_mask, 2);            % [npix x 1]
dmov   = mov2d_mask - avgvec;

% --- PCA via T x T covariance ---
[V, D] = eig(dmov' * dmov);
eigvals = flipud(diag(D));
V = V(:, end:-1:1);
V = V .* sign(max(V) - max(-V));
eigvecs = V;                             % [T x T]

nShow = nframes;                         % all PCs

% --- compute all eigenimages and eigentraces ---
spatialMaps = dmov * V;                  % [npix x nShow]
eigenImgs   = zeros(ysize, xsize, nShow);
for k = 1:nShow
    tmp          = zeros(ysize*xsize, 1);
    tmp(mask(:)) = spatialMaps(:, k);
    eigenImgs(:,:,k) = reshape(tmp, [ysize, xsize]);
end
eigenTraces = V;                         % [T x nShow]  (temporal coefficients)

var_total = sum(eigvals);
fprintf('\nTotal PCs: %d\n', nShow);

% --- interactive figure ---
hFig = figure(997);
set(hFig, 'Name', 'PCA Inspector  —  Prev / Next to browse, Done to continue', ...
    'NumberTitle', 'off', 'Color', [0.15 0.15 0.15]);

d.page        = 1;
d.nShow       = nShow;
d.eigenImgs   = eigenImgs;
d.eigenTraces = eigenTraces;
d.eigvals     = eigvals;
d.var_total   = var_total;
d.nframes     = nframes;
d.ysize       = ysize;
d.xsize       = xsize;
setappdata(hFig, 'd', d);

draw_pca_page(hFig);
uiwait(hFig);

% collect PC selection after figure is done
n = input('PCs to keep (e.g. [1 2 3]): \n');
if ishandle(hFig), close(hFig); end

% --- variance report ---
var_kept = sum(eigvals(n));
fprintf('Variance explained by selected PCs: %.2f%% (of total)\n', ...
    100 * var_kept / var_total);
fprintf('Variance fractions of selected PCs: %s\n', ...
    num2str(eigvals(n)' ./ var_total, '%.4f '));

% --- reconstruct filtered signal in masked region ---
coeffs     = dmov * V(:, n);             % [npix x nKept]
dmov_filt  = coeffs * V(:, n)';         % [npix x T]
mov2d_out  = mov2d;
mov2d_out(mask(:), :) = dmov_filt + avgvec;
moviefilt  = reshape(mov2d_out, [ysize, xsize, nframes]);
end

% -------------------------------------------------------------------------
function draw_pca_page(hFig)
d        = getappdata(hFig, 'd');
nPerPage = 10;
nPages   = ceil(d.nShow / nPerPage);
page     = d.page;
idx      = (page-1)*nPerPage + (1:nPerPage);
idx      = idx(idx <= d.nShow);
nPlot    = numel(idx);

clf(hFig);

% ---- layout constants ----
btnH      = 0.055;
btnMargin = 0.01;
plotArea  = 1 - btnH - 2*btnMargin;   % vertical space for rows
rowH      = plotArea / nPlot;
traceL    = 0.05;  traceW = 0.60;
imgL      = 0.68;  imgW   = 0.28;

% aspect ratio scaling for eigenimage box height
imgAR = d.xsize / d.ysize;            % width/height of image

for k = 1:nPlot
    pcIdx     = idx(k);
    rowBottom = btnH + 2*btnMargin + (nPlot - k) * rowH;
    pct       = d.eigvals(pcIdx) / d.var_total * 100;

    % --- trace ---
    ax1 = axes('Parent', hFig, ...
        'Position', [traceL, rowBottom + rowH*0.1, traceW, rowH*0.78]);
    plot(ax1, d.eigenTraces(:, pcIdx), 'Color', [0.4 0.8 1], 'LineWidth', 0.6);
    xlim(ax1, [1 d.nframes]);
    set(ax1, 'XTick', [], 'YTick', [], 'Color', 'k', ...
        'XColor', [0.6 0.6 0.6], 'YColor', [0.6 0.6 0.6]);
    ylabel(ax1, sprintf('PC%d\n%.2f%%', pcIdx, pct), ...
        'Color', [0.9 0.9 0.9], 'FontSize', 7, 'Rotation', 0, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');

    % --- eigenimage ---
    % scale box so image is not distorted; fit within imgW x rowH*0.78
    boxH  = rowH * 0.78;
    boxW  = min(imgW, boxH * imgAR);
    ax2 = axes('Parent', hFig, ...
        'Position', [imgL, rowBottom + rowH*0.1, boxW, boxH]);
    img = d.eigenImgs(:,:,pcIdx);
    clim = max(abs(img(:)));
    if clim == 0, clim = 1; end
    imagesc(ax2, img, [-clim clim]);
    axis(ax2, 'image', 'off');
    bwr = [linspace(0,1,128)' linspace(0,1,128)' ones(128,1); ...
           ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];
    colormap(ax2, bwr);
end

% ---- navigation buttons ----
btnY = btnMargin;
prevEn = 'off'; if page > 1,      prevEn = 'on'; end
nextEn = 'off'; if page < nPages, nextEn = 'on'; end

uicontrol(hFig, 'Style', 'pushbutton', 'String', '< Prev', ...
    'Units', 'normalized', 'Position', [0.03 btnY 0.15 btnH], ...
    'Enable', prevEn, 'FontSize', 10, ...
    'Callback', @(~,~) set_page(hFig, page - 1));

uicontrol(hFig, 'Style', 'text', ...
    'String', sprintf('Page  %d / %d   (PC %d – %d)', page, nPages, idx(1), idx(end)), ...
    'Units', 'normalized', 'Position', [0.22 btnY 0.38 btnH], ...
    'FontSize', 10, 'BackgroundColor', [0.15 0.15 0.15], 'ForegroundColor', 'w');

uicontrol(hFig, 'Style', 'pushbutton', 'String', 'Next >', ...
    'Units', 'normalized', 'Position', [0.63 btnY 0.15 btnH], ...
    'Enable', nextEn, 'FontSize', 10, ...
    'Callback', @(~,~) set_page(hFig, page + 1));

uicontrol(hFig, 'Style', 'pushbutton', 'String', 'Done — select PCs', ...
    'Units', 'normalized', 'Position', [0.81 btnY 0.17 btnH], ...
    'FontSize', 10, 'BackgroundColor', [0.2 0.6 0.2], 'ForegroundColor', 'w', ...
    'Callback', @(~,~) uiresume(hFig));
end

% -------------------------------------------------------------------------
function set_page(hFig, newPage)
d = getappdata(hFig, 'd');
d.page = newPage;
setappdata(hFig, 'd', d);
draw_pca_page(hFig);
end
