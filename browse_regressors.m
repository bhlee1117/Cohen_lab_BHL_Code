function browse_regressors(traces, maps, sz, names, opts)
% browse_regressors(traces, maps, sz, names, opts)
%
% Step through a set of regressors ONE AT A TIME: its time course on top, its spatial
% map below, and -- when opts.stats is supplied -- a scatter of EVERY regressor in a
% chosen pair of statistics with the current one marked.
%
% Showing many regressors as a tall stack makes each one too small to judge; this shows
% one at full size instead. The scatter puts that one back in context: a regressor is
% suspicious not in isolation but relative to the others (an outlier in frequency, an
% amplitude far above the rest, a spatial map that alone correlates with the neuron).
%
% INPUT
%   traces  [nReg x T]              one regressor per row
%   maps    [nPix x nReg]           its spatial map, nPix = sz(1)*sz(2)
%                                   (pass [] for no map panel)
%   sz      [nrow ncol]             image size the maps reshape to
%   names   {1 x nReg} cellstr      label for each regressor
%   opts    struct, all optional:
%     .preTrace  [nReg x T]  the trace each regressor was DERIVED FROM, before any
%                            filtering. Drawn in its own panel above the regressor, so
%                            you can see what the bandpass + Hilbert actually extracted.
%     .score     [nReg x 1]  a quality score shown in the title
%     .scoreMin  scalar      scores below this turn the title red
%     .f0        [nReg x 1]  frequency shown in the title (Hz)
%     .seg       [nReg x 1]  which source each regressor came from (shown in the title)
%     .segMask   {1 x nSrc}  pixel indices per source, outlined in green on the map
%     .refMask   [nrow x ncol] logical, outlined in black (e.g. a vessel mask)
%     .cmap      colormap for the map          (default: blue-white-red)
%     .clim      [lo hi] for the map           (default: symmetric 99.5th percentile)
%     .tWin      frames plotted                (default: first 1000)
%     .figName   figure name                   (default: 'regressor browser')
%     .stats     struct of [nReg x 1] numeric vectors, one field per statistic
%                (e.g. f0, oscScore, rVolt, ampCov, bandLo, bandHi). Adds the scatter
%                panel and a per-regressor readout of every field. Fields that are not
%                numeric or not nReg long are ignored.
%     .statXY    {'xfield','yfield'} initially plotted   (default: first two fields)
%     .statC     field name used to colour the points    (default: none)
%     .statLog   cellstr of fields to draw on a log axis (default: none)
%     .eigImg    [nPix x nReg] the SPATIAL MODE each regressor came from, painted on the
%                frame (NaN off the segment). May be COMPLEX, in which case it is drawn
%                as phase on a cyclic map with |value| as opacity -- for a travelling
%                wave the propagation is the half of that mode a signed image throws
%                away. Adds a panel beside the covariance map.
%     .mask      [nrow x ncol] logical. Pixels outside it are not drawn on the covariance
%                map, and the colour limits are taken from inside it only -- otherwise a
%                bright background pixel sets the scale and the neuron washes out.
%
% The scatter is clickable: clicking near a point jumps to that regressor.
%
% 2026.08.27, Byung Hun Lee

if nargin < 5, opts = struct; end
def = struct('preTrace',[], 'score',[], 'scoreMin',[], 'f0',[], 'seg',[], 'segMask',{{}}, ...
             'refMask',[], 'cmap',[], 'clim',[], 'tWin',[], 'figName','regressor browser', ...
             'stats',[], 'statXY',{{}}, 'statC','', 'statLog',{{}}, 'eigImg',[], 'mask',[]);
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(opts, fn{k}) || isempty(opts.(fn{k})), opts.(fn{k}) = def.(fn{k}); end
end
nReg = size(traces,1);
assert(nReg > 0, 'No regressors to browse.');
if isempty(opts.tWin), opts.tWin = 1:min(1000, size(traces,2)); end
if isempty(opts.cmap)
    opts.cmap = [linspace(0,1,128)' linspace(0.5,1,128)' ones(128,1); ...
                 ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];
end
if isempty(opts.clim) && ~isempty(maps)
    mv = maps;
    if ~isempty(opts.mask) && numel(opts.mask) == size(maps,1)
        mv = maps(opts.mask(:), :);          % scale from the pixels that matter
    end
    c = prctile(abs(mv(:)), 99.5);
    if ~isfinite(c) || c <= 0, c = max([abs(mv(:)); eps]); end
    opts.clim = c*[-1 1];
end

% ---- keep only the stat fields that are usable ----
statNm = {};
if ~isempty(opts.stats) && isstruct(opts.stats)
    f = fieldnames(opts.stats);
    for k = 1:numel(f)
        v = opts.stats.(f{k});
        if isnumeric(v) && numel(v) == nReg
            opts.stats.(f{k}) = double(v(:));
            statNm{end+1} = f{k};   %#ok<AGROW>
        else
            opts.stats = rmfield(opts.stats, f{k});
        end
    end
end
hasStat = numel(statNm) >= 2;
if hasStat
    ix = 1;  iy = min(2, numel(statNm));
    if numel(opts.statXY) == 2
        jx = find(strcmp(statNm, opts.statXY{1}), 1);
        jy = find(strcmp(statNm, opts.statXY{2}), 1);
        if ~isempty(jx), ix = jx; end
        if ~isempty(jy), iy = jy; end
    end
end

S.traces = traces;  S.maps = maps;  S.sz = sz;  S.names = names;
S.opts = opts;      S.i = 1;        S.n = nReg;
S.statNm = statNm;  S.hasStat = hasStat;
S.clickFcn = @clickScatter;   % set on the scatter axes by draw()
if hasStat, S.ix = ix;  S.iy = iy; end

fh = findobj(groot, 'Type','figure', 'Name', opts.figName);
if isempty(fh), fh = figure('Name', opts.figName, 'NumberTitle','off');
else,           fh = fh(1); figure(fh); end
clf(fh);
set(fh, 'KeyPressFcn', @onKey);

% ---- layout: traces and map on the left, stats on the right when present ----
if hasStat, wL = 0.60; else, wL = 0.88; end
hasPre = ~isempty(opts.preTrace);
if hasPre
    S.axP = axes(fh, 'Position', [0.06 0.760 wL 0.155]);
    S.axT = axes(fh, 'Position', [0.06 0.575 wL 0.155]);
else
    S.axP = [];
    S.axT = axes(fh, 'Position', [0.06 0.62 wL 0.28]);
end
S.axE = [];
if ~isempty(maps)
    if ~isempty(opts.eigImg)
        S.axE = axes(fh, 'Position', [0.06 0.075 wL*0.46 0.44]);
        S.axM = axes(fh, 'Position', [0.06+wL*0.54 0.075 wL*0.46 0.44]);
    else
        S.axM = axes(fh, 'Position', [0.06 0.075 wL 0.44]);
    end
else
    S.axM = [];
end
S.axS = [];  S.txt = [];
if hasStat
    S.axS = axes(fh, 'Position', [0.73 0.575 0.24 0.32]);
    S.txt = annotation(fh, 'textbox', [0.70 0.075 0.29 0.42], 'String','', ...
        'Interpreter','none', 'FontName','Courier', 'FontSize', 8, ...
        'EdgeColor',[0.8 0.8 0.8], 'VerticalAlignment','top');
    uicontrol(fh,'Style','text','String','x','Units','normalized', ...
        'Position',[0.70 0.522 0.02 0.03],'BackgroundColor',get(fh,'Color'));
    S.popX = uicontrol(fh,'Style','popupmenu','String',statNm,'Value',ix, ...
        'Units','normalized','Position',[0.725 0.522 0.11 0.032], ...
        'Callback',@(src,~) setAx('x', src));
    uicontrol(fh,'Style','text','String','y','Units','normalized', ...
        'Position',[0.845 0.522 0.02 0.03],'BackgroundColor',get(fh,'Color'));
    S.popY = uicontrol(fh,'Style','popupmenu','String',statNm,'Value',iy, ...
        'Units','normalized','Position',[0.868 0.522 0.11 0.032], ...
        'Callback',@(src,~) setAx('y', src));
end
uicontrol(fh,'Style','pushbutton','String','< Prev','Units','normalized', ...
    'Position',[0.06 0.945 0.10 0.04],'Callback',@(~,~) step(-1));
uicontrol(fh,'Style','pushbutton','String','Next >','Units','normalized', ...
    'Position',[0.17 0.945 0.10 0.04],'Callback',@(~,~) step(+1));
S.box = uicontrol(fh,'Style','edit','String','1','Units','normalized', ...
    'Position',[0.29 0.945 0.06 0.04],'Callback',@(src,~) jumpTo(src));
uicontrol(fh,'Style','text','String',sprintf('of %d   (arrow keys work too)', nReg), ...
    'Units','normalized','Position',[0.36 0.943 0.35 0.038], ...
    'HorizontalAlignment','left','BackgroundColor',get(fh,'Color'));
guidata(fh, S);
draw(fh);

    function step(d)
        S = guidata(fh);
        S.i = max(1, min(S.n, S.i + d));
        guidata(fh, S);  draw(fh);
    end
    function jumpTo(src)
        S = guidata(fh);
        v = round(str2double(get(src,'String')));
        if isfinite(v), S.i = max(1, min(S.n, v)); end
        guidata(fh, S);  draw(fh);
    end
    function setAx(which, src)
        S = guidata(fh);
        if strcmp(which,'x'), S.ix = get(src,'Value'); else, S.iy = get(src,'Value'); end
        guidata(fh, S);  draw(fh);
    end
    function onKey(~, e)
        switch e.Key
            case {'rightarrow','downarrow'}, step(+1);
            case {'leftarrow','uparrow'},    step(-1);
            case 'home', S = guidata(fh); S.i = 1;    guidata(fh,S); draw(fh);
            case 'end',  S = guidata(fh); S.i = S.n;  guidata(fh,S); draw(fh);
        end
    end
    function clickScatter(~, ~)
        S = guidata(fh);
        p = get(S.axS, 'CurrentPoint');
        xv = S.opts.stats.(S.statNm{S.ix});  yv = S.opts.stats.(S.statNm{S.iy});
        xl = xlim(S.axS);  yl = ylim(S.axS);          % normalise so both axes count
        dx = (xv - p(1,1)) / max(diff(xl), eps);
        dy = (yv - p(1,2)) / max(diff(yl), eps);
        [~, j] = min(dx.^2 + dy.^2);
        S.i = j;  guidata(fh, S);  draw(fh);
    end
end

function draw(fh)
S = guidata(fh);  o = S.opts;  i = S.i;
set(S.box, 'String', num2str(i));

% ---- title text, assembled from whatever opts provided ----
txt = sprintf('%d/%d   %s', i, S.n, S.names{i});
if ~isempty(o.score), txt = sprintf('%s   score %.3f', txt, o.score(i)); end
if ~isempty(o.f0),    txt = sprintf('%s   f0 %.0f Hz', txt, o.f0(i));    end
if ~isempty(o.seg),   txt = sprintf('%s   from source %d', txt, o.seg(i)); end
tc = [0 0 0];
if ~isempty(o.score) && ~isempty(o.scoreMin) && o.score(i) < o.scoreMin
    tc = [0.75 0 0];  txt = [txt '   (below threshold)'];
end

if ~isempty(S.axP)
    cla(S.axP);
    plot(S.axP, o.tWin, o.preTrace(i, o.tWin), '-', 'Color', [0.45 0.45 0.45], ...
        'LineWidth', 0.7);
    xlim(S.axP, [o.tWin(1) o.tWin(end)]);  box(S.axP,'off');  set(S.axP,'XTick',[]);
    ylabel(S.axP, 'source');
    title(S.axP, txt, 'Color', tc, 'Interpreter','none');
end

cla(S.axT);
plot(S.axT, o.tWin, S.traces(i, o.tWin), 'k-', 'LineWidth', 0.7);
xlim(S.axT, [o.tWin(1) o.tWin(end)]);  box(S.axT,'off');
xlabel(S.axT, 'frame');  ylabel(S.axT, 'regressor (z)');
if isempty(S.axP)
    title(S.axT, txt, 'Color', tc, 'Interpreter','none');
else
    title(S.axT, 'after bandpass + Hilbert');
end

% ---- scatter of every regressor, current one marked ----
if S.hasStat
    xNm = S.statNm{S.ix};  yNm = S.statNm{S.iy};
    xv = o.stats.(xNm);    yv = o.stats.(yNm);
    cla(S.axS);  hold(S.axS, 'on');
    cv = [];
    if ~isempty(o.statC) && isfield(o.stats, o.statC), cv = o.stats.(o.statC); end
    if isempty(cv)
        scatter(S.axS, xv, yv, 26, [0.55 0.6 0.7], 'filled', 'MarkerFaceAlpha', 0.65);
    else
        scatter(S.axS, xv, yv, 26, cv, 'filled', 'MarkerFaceAlpha', 0.85);
        cb = colorbar(S.axS);  cb.Label.String = o.statC;
    end
    % the current regressor: a ring plus crosshairs, so it is findable in a dense cloud
    plot(S.axS, xv(i), yv(i), 'o', 'MarkerSize', 13, 'LineWidth', 1.8, ...
        'MarkerEdgeColor', [0.85 0.1 0.1]);
    xline(S.axS, xv(i), ':', 'Color', [0.85 0.1 0.1]);
    yline(S.axS, yv(i), ':', 'Color', [0.85 0.1 0.1]);
    if ~isempty(o.scoreMin) && any(strcmp({xNm, yNm}, 'score'))
        if strcmp(xNm,'score'), xline(S.axS, o.scoreMin, 'k--');
        else,                   yline(S.axS, o.scoreMin, 'k--'); end
    end
    if any(strcmp(o.statLog, xNm)) && all(xv > 0), set(S.axS,'XScale','log'); end
    if any(strcmp(o.statLog, yNm)) && all(yv > 0), set(S.axS,'YScale','log'); end
    box(S.axS,'off');  xlabel(S.axS, xNm, 'Interpreter','none');
    ylabel(S.axS, yNm, 'Interpreter','none');
    title(S.axS, 'all regressors (click to jump)');
    % clicks must reach the AXES, not the markers drawn on top of it
    set(S.axS, 'ButtonDownFcn', S.clickFcn, 'HitTest', 'on', 'PickableParts', 'all');
    set(allchild(S.axS), 'HitTest', 'off');

    % ---- the current regressor's numbers, every field ----
    L = {sprintf('#%d of %d', i, S.n), S.names{i}, ''};
    for q = 1:numel(S.statNm)
        v = o.stats.(S.statNm{q})(i);
        if abs(v) >= 1e4 || (v ~= 0 && abs(v) < 1e-3)
            L{end+1} = sprintf('%-12s %10.3e', S.statNm{q}, v);   %#ok<AGROW>
        else
            L{end+1} = sprintf('%-12s %10.4f', S.statNm{q}, v);   %#ok<AGROW>
        end
    end
    set(S.txt, 'String', L);
end

% ---- the spatial mode this regressor came from ----
if ~isempty(S.axE)
    cla(S.axE);
    E = reshape(o.eigImg(:,i), S.sz(1), S.sz(2));
    if ~isreal(E)
        hE = imagesc(S.axE, angle(E), [-pi pi]);
        set(hE, 'AlphaData', abs(E) .* ~isnan(abs(E)));
        colormap(S.axE, hsv);
        cbE = colorbar(S.axE, 'Ticks', [-pi 0 pi], 'TickLabels', {'-\pi','0','\pi'});
        cbE.Label.String = 'phase';
        ttlE = 'mode: phase (colour) x |U| (opacity)';
    else
        hE = imagesc(S.axE, E, [-1 1]);
        set(hE, 'AlphaData', ~isnan(E));
        colormap(S.axE, o.cmap);  colorbar(S.axE);
        ttlE = 'spatial mode (eigenvector)';
    end
    axis(S.axE, 'image', 'off');  hold(S.axE, 'on');
    if ~isempty(o.refMask) && any(o.refMask(:))
        contour(S.axE, o.refMask, [0.5 0.5], 'k', 'LineWidth', 0.7);
    end
    title(S.axE, ttlE);
end

if isempty(S.axM), return; end
cla(S.axM);
Mi = reshape(S.maps(:,i), S.sz(1), S.sz(2));
if ~isempty(o.mask) && isequal(size(o.mask), S.sz)
    Mi(~o.mask) = NaN;            % outside the analysis mask the covariance is not asked
end
hM = imagesc(S.axM, Mi, o.clim);  set(hM, 'AlphaData', ~isnan(Mi));
axis(S.axM, 'image', 'off');  colormap(S.axM, o.cmap);  colorbar(S.axM);
hold(S.axM, 'on');
if ~isempty(o.refMask)
    contour(S.axM, o.refMask, [0.5 0.5], 'k', 'LineWidth', 0.7);
end
if ~isempty(o.segMask) && ~isempty(o.seg)
    m = false(S.sz);  m(o.segMask{o.seg(i)}) = true;
    contour(S.axM, m, [0.5 0.5], 'g', 'LineWidth', 1.4);
end
title(S.axM, 'covariance with mov\_res   (green = source segment)');
end
