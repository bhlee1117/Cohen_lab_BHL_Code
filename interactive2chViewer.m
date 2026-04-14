function interactive2chViewer(GluTrace, GluFootprint, GluTime, VoltTrace, VoltFootprint, VoltTime)
% interactive2chViewer
%
% Inputs:
%   GluTrace      : N x T
%   GluFootprint  : X x Y x N  (positive pixels define ROI)
%   GluTime       : 1 x T
%   VoltTrace     : M x T2
%   VoltFootprint : X x Y x M
%   VoltTime      : 1 x T2
%
% Interaction:
%   Left-click  on footprint panel : toggle trace on/off (add if hidden, remove if shown)
%   Right-click on footprint panel : remove trace
%   Footprint axes are linked (zoom/pan together)
%   Trace x-axes are linked

    %% ---------------- Input checks ----------------
    validateattributes(GluTrace,      {'numeric'}, {'2d',     'nonempty'});
    validateattributes(GluFootprint,  {'numeric'}, {'3d',     'nonempty'});
    validateattributes(GluTime,       {'numeric'}, {'vector', 'nonempty'});
    validateattributes(VoltTrace,     {'numeric'}, {'2d',     'nonempty'});
    validateattributes(VoltFootprint, {'numeric'}, {'3d',     'nonempty'});
    validateattributes(VoltTime,      {'numeric'}, {'vector', 'nonempty'});

    [N,  T ] = size(GluTrace);
    [X1, Y1, Nf] = size(GluFootprint);
    [M,  T2] = size(VoltTrace);
    [X2, Y2, Mf] = size(VoltFootprint);

    if N  ~= Nf,          error('size(GluTrace,1) must equal size(GluFootprint,3).'); end
    if M  ~= Mf,          error('size(VoltTrace,1) must equal size(VoltFootprint,3).'); end
    if numel(GluTime)~=T, error('numel(GluTime) must equal size(GluTrace,2).'); end
    if numel(VoltTime)~=T2, error('numel(VoltTime) must equal size(VoltTrace,2).'); end
    if X1~=X2 || Y1~=Y2, error('GluFootprint and VoltFootprint must share X x Y dimensions.'); end

    GluTime  = GluTime(:).';
    VoltTime = VoltTime(:).';

    %% ---------------- Colours ----------------
    gluColors  = turbo(N);
    voltColors = turbo(M);

    %% ---------------- Figure ----------------
    fig = figure( ...
        'Name',        'Interactive Glu / Voltage Viewer', ...
        'NumberTitle', 'off', ...
        'Color',       'k');

    axGluFP  = subplot(2,2,1, 'Parent', fig);
    axGluTr  = subplot(2,2,2, 'Parent', fig);
    axVoltFP = subplot(2,2,3, 'Parent', fig);
    axVoltTr = subplot(2,2,4, 'Parent', fig);

    for ax = [axGluFP axGluTr axVoltFP axVoltTr]
        ax.Color  = 'k';
        ax.XColor = 'w';
        ax.YColor = 'w';
        hold(ax, 'on');
    end

    %% ---------------- Reference images ----------------
    gluRefIm  = mean(max(GluFootprint,  0), 3);
    voltRefIm = mean(max(VoltFootprint, 0), 3);

    imshow2(gluRefIm,  [], 'Parent', axGluFP);   hold(axGluFP,  'on');
    imshow2(voltRefIm, [], 'Parent', axVoltFP);  hold(axVoltFP, 'on');

    axis(axGluFP,  'image');
    axis(axVoltFP, 'image');
    set([axGluFP axVoltFP], 'YDir','reverse');

    %% ---------------- Draw contours (bwboundaries style) ----------------
    gluContourHandles  = cell(N, 1);   % each cell holds array of line handles
    voltContourHandles = cell(M, 1);
    gluCoord  = footprint_centroids(GluFootprint);
    voltCoord = footprint_centroids(VoltFootprint);

    for k = 1:N
        mask = GluFootprint(:,:,k) > 0;
        bnds = bwboundaries(mask);
        gluContourHandles{k} = gobjects(numel(bnds), 1);
        for j = 1:numel(bnds)
            gluContourHandles{k}(j) = plot(axGluFP, ...
                bnds{j}(:,2), bnds{j}(:,1), ...
                'Color',     gluColors(k,:), ...
                'LineWidth', 1.5, ...
                'HitTest',       'off', ...
                'PickableParts', 'none', ...
                'Tag', sprintf('glu_roi_%d', k));
        end
        text(axGluFP, gluCoord(k,1), gluCoord(k,2), num2str(k), ...
            'Color','w', 'FontSize',10, 'FontWeight','bold', ...
            'HorizontalAlignment','center', ...
            'HitTest','off', 'PickableParts','none');
    end

    for k = 1:M
        mask = VoltFootprint(:,:,k) > 0;
        bnds = bwboundaries(mask);
        voltContourHandles{k} = gobjects(numel(bnds), 1);
        for j = 1:numel(bnds)
            voltContourHandles{k}(j) = plot(axVoltFP, ...
                bnds{j}(:,2), bnds{j}(:,1), ...
                'Color',     voltColors(k,:), ...
                'LineWidth', 1.5, ...
                'HitTest',       'off', ...
                'PickableParts', 'none', ...
                'Tag', sprintf('volt_roi_%d', k));
        end
        text(axVoltFP, voltCoord(k,1), voltCoord(k,2), num2str(k), ...
            'Color','w', 'FontSize',10, 'FontWeight','bold', ...
            'HorizontalAlignment','center', ...
            'HitTest','off', 'PickableParts','none');
    end

    %% ---------------- Invisible click-capture layers ----------------
    % Must come AFTER drawing so they sit on top and catch clicks
    hClick1 = imagesc(axGluFP,  zeros(X1,Y1));
    set(hClick1, 'AlphaData',0, 'HitTest','on', 'PickableParts','all', ...
        'ButtonDownFcn', @(s,e) footprintClickCb(s,e,fig,'glu'));

    hClick2 = imagesc(axVoltFP, zeros(X1,Y1));
    set(hClick2, 'AlphaData',0, 'HitTest','on', 'PickableParts','all', ...
        'ButtonDownFcn', @(s,e) footprintClickCb(s,e,fig,'volt'));

    %% ---------------- Trace panel styling ----------------
    title(axGluFP,  'Glu Footprints',  'Color','w');
    title(axGluTr,  'Glu Traces',      'Color','w');
    title(axVoltFP, 'Volt Footprints', 'Color','w');
    title(axVoltTr, 'Volt Traces',     'Color','w');

    xlabel(axGluTr,  'Time', 'Color','w');  ylabel(axGluTr,  'Signal', 'Color','w');
    xlabel(axVoltTr, 'Time', 'Color','w');  ylabel(axVoltTr, 'Signal', 'Color','w');
    set(axGluTr,  'XLim', [GluTime(1)  GluTime(end)]);
    set(axVoltTr, 'XLim', [VoltTime(1) VoltTime(end)]);

    %% ---------------- Link axes ----------------
    % Link footprint panels (zoom/pan together)
    linkaxes([axGluFP, axVoltFP], 'xy');
    % Link trace panels x-axis
    linkaxes([axGluTr, axVoltTr], 'x');

    %% ---------------- State ----------------
    state.GluTrace      = GluTrace;
    state.GluFootprint  = GluFootprint;
    state.GluTime       = GluTime;
    state.VoltTrace     = VoltTrace;
    state.VoltFootprint = VoltFootprint;
    state.VoltTime      = VoltTime;

    state.axGluFP  = axGluFP;
    state.axGluTr  = axGluTr;
    state.axVoltFP = axVoltFP;
    state.axVoltTr = axVoltTr;

    state.gluColors  = gluColors;
    state.voltColors = voltColors;

    state.gluContourHandles  = gluContourHandles;
    state.voltContourHandles = voltContourHandles;

    state.gluCoord  = gluCoord;
    state.voltCoord = voltCoord;

    state.gluTraceLines  = gobjects(N, 1);
    state.voltTraceLines = gobjects(M, 1);

    state.gluSelected  = false(N, 1);
    state.voltSelected = false(M, 1);

    guidata(fig, state);
end

%% =========================================================================
function footprintClickCb(~, ~, fig, mode)
    state = guidata(fig);

    if strcmpi(mode, 'glu')
        ax = state.axGluFP;
        F  = state.GluFootprint;
        coord = state.gluCoord;
    else
        ax = state.axVoltFP;
        F  = state.VoltFootprint;
        coord = state.voltCoord;
    end

    cp     = get(ax, 'CurrentPoint');
    xClick = cp(1,2);   % row
    yClick = cp(1,1);   % col

    [X, Y, K] = size(F);
    xIdx = round(xClick);
    yIdx = round(yClick);
    if xIdx < 1 || xIdx > X || yIdx < 1 || yIdx > Y, return; end

    idx = pickROI(F, coord, xIdx, yIdx);
    if isempty(idx), return; end

    clickType = get(fig, 'SelectionType');
    switch clickType
        case 'normal'   % left click -> toggle
            toggleComponent(fig, mode, idx);
        case 'alt'      % right click -> always remove
            removeComponent(fig, mode, idx);
    end
end

%% =========================================================================
function idx = pickROI(F, coord, xIdx, yIdx)
% 1) Any ROI mask that contains the clicked pixel -> highest value wins
% 2) Fallback: nearest centroid within 30 px

    K = size(F, 3);
    score = -inf(K, 1);

    for k = 1:K
        A = double(F(:,:,k));
        if A(xIdx, yIdx) > 0
            mx = max(A(:));
            score(k) = A(xIdx, yIdx) / max(mx, eps);
        end
    end

    if any(isfinite(score) & score > -inf)
        [~, idx] = max(score);
        return;
    end

    % Fallback: nearest centroid
    dists = sqrt((coord(:,1) - yIdx).^2 + (coord(:,2) - xIdx).^2);
    [~, idx] = min(dists);
    if dists(idx) > 30, idx = []; end
end

%% =========================================================================
function toggleComponent(fig, mode, idx)
% Left-click: if already selected -> remove; if not -> add

    state = guidata(fig);

    switch lower(mode)
        case 'glu',  already = state.gluSelected(idx);
        case 'volt', already = state.voltSelected(idx);
    end

    if already
        removeComponent(fig, mode, idx);
    else
        addComponent(fig, mode, idx);
    end
end

%% =========================================================================
function addComponent(fig, mode, idx)
    state = guidata(fig);

    switch lower(mode)
        case 'glu'
            if state.gluSelected(idx), return; end
            c = state.gluColors(idx,:);
            state.gluTraceLines(idx) = plot(state.axGluTr, ...
                state.GluTime, state.GluTrace(idx,:), ...
                'Color', c, 'LineWidth', 1.5, ...
                'DisplayName', sprintf('Glu %d', idx));
            state.gluSelected(idx) = true;
            setContourWidth(state.gluContourHandles{idx}, 2.8);
            refreshLegend(state.axGluTr);

        case 'volt'
            if state.voltSelected(idx), return; end
            c = state.voltColors(idx,:);
            state.voltTraceLines(idx) = plot(state.axVoltTr, ...
                state.VoltTime, state.VoltTrace(idx,:), ...
                'Color', c, 'LineWidth', 1.5, ...
                'DisplayName', sprintf('Volt %d', idx));
            state.voltSelected(idx) = true;
            setContourWidth(state.voltContourHandles{idx}, 2.8);
            refreshLegend(state.axVoltTr);
    end

    guidata(fig, state);
end

%% =========================================================================
function removeComponent(fig, mode, idx)
    state = guidata(fig);

    switch lower(mode)
        case 'glu'
            if ~state.gluSelected(idx), return; end
            if isgraphics(state.gluTraceLines(idx))
                delete(state.gluTraceLines(idx));
            end
            state.gluTraceLines(idx) = gobjects(1);
            state.gluSelected(idx)   = false;
            setContourWidth(state.gluContourHandles{idx}, 1.5);
            refreshLegend(state.axGluTr);

        case 'volt'
            if ~state.voltSelected(idx), return; end
            if isgraphics(state.voltTraceLines(idx))
                delete(state.voltTraceLines(idx));
            end
            state.voltTraceLines(idx) = gobjects(1);
            state.voltSelected(idx)   = false;
            setContourWidth(state.voltContourHandles{idx}, 1.5);
            refreshLegend(state.axVoltTr);
    end

    guidata(fig, state);
end

%% =========================================================================
function setContourWidth(handles, lw)
% Set LineWidth on all boundary line handles for one ROI

    for j = 1:numel(handles)
        if isgraphics(handles(j))
            set(handles(j), 'LineWidth', lw);
        end
    end
end

%% =========================================================================
function refreshLegend(ax)
    h = findobj(ax, 'Type','line');
    h = h(isgraphics(h));
    if isempty(h)
        legend(ax, 'off'); return;
    end
    % get() returns char when h is scalar, cell when h is array -- normalise
    names = get(h, 'DisplayName');
    if ischar(names), names = {names}; end
    keep = ~cellfun(@isempty, names);
    h    = h(keep);
    if isempty(h)
        legend(ax, 'off');
    else
        legend(ax, flipud(h(:)), 'TextColor','w', 'Color','none', ...
            'EdgeColor','w', 'Location','best');
    end
end

%% =========================================================================
function coord = footprint_centroids(fp)
    [Y, X] = size(fp(:,:,1));
    [xg, yg] = meshgrid(1:X, 1:Y);
    N = size(fp, 3);
    coord = zeros(N, 2);
    for i = 1:N
        w = max(fp(:,:,i), 0);
        s = sum(w(:));
        if s > 0
            coord(i,1) = sum(sum(w .* xg)) / s;   % x (col)
            coord(i,2) = sum(sum(w .* yg)) / s;   % y (row)
        else
            coord(i,:) = [X/2, Y/2];
        end
    end
end
