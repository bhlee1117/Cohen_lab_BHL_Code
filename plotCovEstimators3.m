function plotCovEstimators3(E, PT, opts)
% plotCovEstimators3  Figure for covEstimators3: correlation matrix, PC maps,
% cumulative variance -- one row per estimator.
%
%   plotCovEstimators3(E, PT, opts)
%
% Inputs
%   E    : output of covEstimators3.
%   PT   : painting info for the PC maps, struct with
%            .aMax  H x W, index of the ROI owning each pixel (from max(ft,[],3))
%            .supp  H x W logical, pixels owned by any ROI
%            .siz   [H W]
%   opts : struct, all optional
%            .figName   figure name                       (default 'covEstimators3')
%            .titleStr  title for the whole layout         (default '')
%            .nPC       PC maps per estimator              (default 5)
%            .ord       ROI order for the matrices         (default 1:nROI)
%            .corrClim  correlation colour limits          (default [-0.25 1])
%            .paintPrc  percentile setting the PC map range(default 99)
%            .cmap      diverging colormap                 (default blue-white-red)
%
% See also: covEstimators3

if nargin < 3, opts = struct(); end
d = @(f,v) ternaryField(opts, f, v);
figName  = d('figName',  'covEstimators3');
titleStr = d('titleStr', '');
nPC      = d('nPC',      5);
ord      = d('ord',      1:E.nROI);
corrClim = d('corrClim', [-0.25 1]);
paintPrc = d('paintPrc', 99);
if isfield(opts,'cmap'), cmap = opts.cmap;
else
    try
        cmap = gen_colormap([0 0.5 1; 1 1 1; 1 0 0], 256);
    catch
        cmap = [linspace(0,1,128)' linspace(.5,1,128)' ones(128,1); ...
                ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];
    end
end
if numel(ord) ~= E.nROI, ord = 1:E.nROI; end

nCol = 1 + nPC;
try
    get_figure(figName);
catch
    figure('Name', figName, 'NumberTitle','off'); clf;
end
tl = tiledlayout(4, nCol, 'TileSpacing','compact','Padding','compact');

for e = 1:3
    ax = nexttile((e-1)*nCol + 1);
    imagesc(ax, E.Corr{e}(ord,ord), corrClim); axis(ax,'image');
    colormap(ax, cmap); set(ax,'XTick',[],'YTick',[]);
    ylabel(ax, E.names{e}, 'FontSize', 7);
    if e == 1, title(ax, 'correlation matrix', 'FontSize', 9); end
    if e == 3, colorbar(ax,'southoutside'); end

    tot = sum(E.Eval{e}(E.Eval{e} > 0));
    for k = 1:nPC
        ax = nexttile((e-1)*nCol + 1 + k);
        v = E.Evec{e}(:,k);
        [~, im] = max(abs(v));  v = v * sign(v(im));        % SVD sign is arbitrary
        IM = nan(PT.siz);
        IM(PT.supp) = v(PT.aMax(PT.supp));
        cl = prctile(abs(v), paintPrc);  if cl <= 0, cl = eps; end
        imagesc(ax, IM, 'AlphaData', PT.supp); set(ax,'Color',[.12 .12 .12]);
        clim(ax, [-cl cl]); axis(ax,'image','off'); colormap(ax, cmap);
        title(ax, sprintf('PC %d  (%.1f%%)', k, 100*E.Eval{e}(k)/tot), 'FontSize', 8);
    end
end

ax = nexttile(3*nCol + 1, [1 nCol]); hold(ax,'on');
cm = [0 0 0; 0.85 0.33 0.10; 0.00 0.45 0.74];
for e = 1:3
    D = E.Eval{e};  pos = D > 0;
    plot(ax, 1:nnz(pos), 100*cumsum(D(pos))/sum(D(pos)), '-o', 'Color', cm(e,:), ...
         'LineWidth',1.5, 'MarkerSize',3, 'MarkerFaceColor', cm(e,:), 'DisplayName', E.names{e});
    if nnz(pos) < numel(D)
        xline(ax, nnz(pos), ':', 'Color', cm(e,:), 'HandleVisibility','off');
    end
end
grid(ax,'on'); xlim(ax,[1 E.nROI]); ylim(ax,[0 100]);
xlabel(ax,'PC #'); ylabel(ax,'cumulative variance (%)');
legend(ax,'Location','southeast','Box','off','FontSize',7);
title(ax, ['cumulative variance over POSITIVE eigenvalues; dotted line = where an ' ...
           'unbiased estimator runs out of them'], 'FontSize', 9);
if ~isempty(titleStr), title(tl, titleStr, 'Interpreter','none'); end
end

function v = ternaryField(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
