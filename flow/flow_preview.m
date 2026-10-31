function flow_preview(parent, stage, R, cfg, data)
% flow_preview(parent, stage, R, cfg, data)
%
% Draw a stage's result into the GUI's right-hand panel. Each stage lays out its own
% subplots inside `parent`, so a stage can show as many views as it needs -- and a movie,
% via flow_movieAxes, when the thing being judged is temporal rather than spatial.
%
% Adding a stage means adding a case. A stage with no case lists what its result
% contains, which is enough to confirm it ran.
%
% 2026.08.28, Byung Hun Lee

delete(allchild(parent));
if ~isfield(R, stage)
    ax = axes(parent, 'Position',[0 0 1 1]);
    text(ax, 0.5, 0.5, sprintf('%s has not been run', stage), ...
        'HorizontalAlignment','center', 'Color',[0.5 0.5 0.5], 'FontSize', 11);
    axis(ax,'off');  return
end
res = R.(stage);

switch stage
% =====================================================================================
case 'bandpass'
    % what the filter IS, what it ACHIEVED, and what it did to the movie
    a1 = axes(parent,'Position',[0.07 0.70 0.36 0.26]);
    if strcmpi(res.method,'dog')
        imagesc(a1, res.kernel); axis(a1,'image');
        clim(a1, max(abs(res.kernel(:)))*[-1 1]);
        colormap(a1, gen_colormap([0 0.5 1; 1 1 1; 1 0 0], 256));
        title(a1, 'DoG kernel', 'FontSize', 9);
    else
        imagesc(a1, res.kernel); axis(a1,'image'); colormap(a1, turbo);
        title(a1, 'Butterworth band (2D)', 'FontSize', 9);
    end
    set(a1,'XTick',[],'YTick',[]);  colorbar(a1);

    a2 = axes(parent,'Position',[0.58 0.70 0.36 0.26]);
    semilogx(a2, res.radial.perUm, res.radial.gain/max(res.radial.gain), ...
        '-', 'LineWidth', 1.4);  hold(a2,'on');
    xline(a2, cfg.bandpass.loUm, 'r--');  xline(a2, cfg.bandpass.hiUm, 'r--');
    yline(a2, 0.5, 'k:');
    set(a2,'XDir','reverse');  xlim(a2,[1 200]);  grid(a2,'on');  box(a2,'off');
    xlabel(a2,'spatial period (\mum)','FontSize',8);  ylabel(a2,'gain','FontSize',8);
    title(a2, sprintf('achieved %.1f-%.1f \\mum (asked %g-%g)', res.bandUm(1), ...
        res.bandUm(2), cfg.bandpass.loUm, cfg.bandpass.hiUm), 'FontSize', 9);

    a3 = axes(parent,'Position',[0.07 0.40 0.36 0.24]);
    showMap(a3, res.varRes, data, 'var(residual)  [input]');
    a4 = axes(parent,'Position',[0.58 0.40 0.36 0.24]);
    % When the ring reference was used, paint WHERE it was taken from over the output
    % variance. Which pixels supplied the "voltage" is the one judgement call this stage
    % makes that a number cannot settle -- a ring sitting on a vessel is obvious in a
    % picture and invisible in a variance ratio.
    if isfield(res,'volt') && strcmp(res.volt.mode,'ring') && ~isempty(res.volt.W)
        showMap(a4, res.varS, data, sprintf('var(S)  [output]   cyan = the %d ring ROIs removed', ...
            res.volt.n));
        hold(a4,'on');
        for q = 1:size(res.volt.W, 2)
            contour(a4, reshape(res.volt.W(:,q) > 0, size(res.varS)), [0.5 0.5], ...
                'LineWidth', 0.6, 'LineColor', [0 1 1]);
        end
    else
        showMap(a4, res.varS, data, 'var(S)  [output]');
    end

    % the movie: input or output, one fixed intensity scale either way
    pop = uicontrol(parent,'Style','popupmenu','String',{'S (band-passed)','mov\_nv (input)'}, ...
        'Units','normalized','Position',[0.07 0.335 0.26 0.035], ...
        'Callback',@(s,~) swapMovie(s));
    drawMovie(parent, res.S, 'S (band-passed)');
    set(pop,'UserData',res);

% =====================================================================================
case 'track'
    % PANEL 1 -- every detection, coloured by polarity. This is the raw material: if the
    % points do not trace out vessels here, nothing downstream can rescue it, and the
    % threshold (cfg.track.thrSD) is what to change.
    a1 = axes(parent,'Position',[0.07 0.56 0.86 0.39]);
    bg = res.speedMap;  bg(:) = 0;
    imagesc(a1, zeros(size(res.speedMap)));  axis(a1,'image','off');
    colormap(a1, gray);  hold(a1,'on');
    if isfield(data,'bvBW') && any(data.bvBW(:))
        contour(a1, data.bvBW, [0.5 0.5], 'Color', [0.2 0.8 0.2], 'LineWidth', 0.8);
    end
    xlim(a1, [0.5 size(res.speedMap,2)+0.5]);      % the frame, not the outliers
    ylim(a1, [0.5 size(res.speedMap,1)+0.5]);
    if isfield(res,'det') && ~isempty(res.det)
        dk = res.det(res.det(:,4) < 0, :);      % dark   -> red
        br = res.det(res.det(:,4) > 0, :);      % bright -> blue
        if ~isempty(dk), plot(a1, dk(:,1), dk(:,2), '.', 'Color',[0.85 0.15 0.15], 'MarkerSize', 3); end
        if ~isempty(br), plot(a1, br(:,1), br(:,2), '.', 'Color',[0.15 0.35 0.90], 'MarkerSize', 3); end
        title(a1, sprintf('detections: %d dark (red), %d bright (blue), %d frames', ...
            size(dk,1), size(br,1), numel(res.frames)), 'FontSize', 9);
    else
        title(a1, 'detections (re-run flow\_track to populate)', 'FontSize', 9);
    end

    % PANEL 2 -- S itself, zoomed on where the tracking actually found something, with
    % the trajectories drawn ON the frame they belong to. Scrolling through it is the
    % only way to tell a real cell train from a threshold picking up noise: real cells
    % move smoothly from frame to frame, noise blinks in and out at random places.
    Sm = bpMovie(R);
    if isempty(Sm)
        axes(parent,'Position',[0.07 0.05 0.86 0.42]); axis off
        text(0.5,0.5,'run step 1 to see the movie','HorizontalAlignment','center');
        return
    end
    fr = res.frames;
    Sm = Sm(:,:,fr);
    kept = find(res.keep);
    if isempty(kept), kept = 1:numel(res.all); end
    cen = res.cen(kept,:);
    % zoom to the bounding box of the kept tracklets, with a margin, so the vessels fill
    % the panel instead of being a few pixels in a wide field
    pad = 0.25;
    x1 = max(1, floor(min(cen(:,1)) - pad*range(cen(:,1)) - 8));
    x2 = min(size(Sm,2), ceil(max(cen(:,1)) + pad*range(cen(:,1)) + 8));
    y1 = max(1, floor(min(cen(:,2)) - pad*range(cen(:,2)) - 8));
    y2 = min(size(Sm,1), ceil(max(cen(:,2)) + pad*range(cen(:,2)) + 8));
    if x2 - x1 < 16, x2 = min(size(Sm,2), x1+16); end
    if y2 - y1 < 16, y2 = min(size(Sm,1), y1+16); end
    Zm = -Sm(y1:y2, x1:x2, :);       % show it spike-positive, as everything else is

    % flatten the kept tracklets once: [x y frame trackId], shifted into the zoom box
    P = [];
    for q = kept(:)'
        r = res.all{q};
        P = [P; r(:,1)-x1+1, r(:,2)-y1+1, r(:,5), repmat(q, size(r,1), 1)]; %#ok<AGROW>
    end
    ov = @(ax, k) drawTracks(ax, P, k, 10);
    ax2 = flow_movieAxes(parent, [0.07 0.045 0.86 0.44], Zm, struct( ...
        'title', sprintf('S zoomed to x %d-%d, y %d-%d  (points bold, 10-frame trail light)', ...
                         x1, x2, y1, y2), ...
        'cmap', gray(256), 'overlay', ov));
    set(ax2, 'Tag', 'flowMovie');

% =====================================================================================
case 'segments'
    st = res.stat;  cu = res.cut;

    % ---- LEFT TOP: what the sort did, in the two quantities that decide it -----------
    % Duration against net travel, each tracklet a point, coloured by the FIRST cut it
    % failed. Colouring by the first failure rather than a single grey "rejected" is what
    % makes it possible to see WHICH threshold is doing the work.
    a1 = axes(parent,'Position',[0.09 0.58 0.36 0.35]);  hold(a1,'on');
    grp = {res.keep,                                [0.85 0.15 0.15], 'kept'; ...
           ~cu.net,                                 [0.55 0.55 0.55], 'went nowhere'; ...
            cu.net & ~cu.len,                       [0.45 0.30 0.65], 'too brief'; ...
            cu.net &  cu.len & ~cu.amp,             [0.95 0.65 0.15], 'too dim'; ...
            cu.net &  cu.len &  cu.amp & ~cu.per,   [0.30 0.45 0.85], 'wanders'};
    lg = {};
    for q = 1:size(grp,1)
        if ~any(grp{q,1}), continue, end
        scatter(a1, st.len(grp{q,1}), st.net(grp{q,1}), 16, grp{q,2}, 'filled', ...
            'MarkerFaceAlpha', 0.6);
        lg{end+1} = sprintf('%s (%d)', grp{q,3}, nnz(grp{q,1})); %#ok<AGROW>
    end
    xline(a1, cfg.segments.minLenKeep, 'k--');
    yline(a1, cfg.segments.minNetUm,  'k--');
    box(a1,'off'); xlabel(a1,'track length (frames)','FontSize',8);
    ylabel(a1,'displacement, first to last (\mum)','FontSize',8);
    legend(a1, lg, 'Location','southeast','FontSize',7);
    title(a1, sprintf('SORT: %d of %d kept', nnz(res.keep), numel(res.keep)), 'FontSize', 9);

    % ---- RIGHT TOP: the surviving tracklets, on the dendrite ------------------------
    % Colour = which segment the tracklet was assigned to, so this panel and the segment
    % map below can be read against each other. Grey = kept by the sort but landing in no
    % accepted segment, which is how you see a vessel that the density mask has cut up.
    a2 = axes(parent,'Position',[0.57 0.58 0.36 0.35]);
    bgm = [];
    if isfield(data,'neuronMap') && any(data.neuronMap(:)), bgm = data.neuronMap; end
    if isempty(bgm), bgm = res.densMap; end
    imagesc(a2, bgm, [0 prctile(bgm(:), 99.5)]); axis(a2,'image','off');
    colormap(a2, gray); hold(a2,'on');
    if isfield(data,'bvBW') && any(data.bvBW(:))
        contour(a2, data.bvBW, [0.5 0.5], 'Color',[0.2 0.9 0.9], 'LineWidth', 0.8);
    end
    cmapS = lines(max(1, numel(res.px)));
    % the paths live in the track result, not the segment one -- read them from there
    % rather than duplicating every tracklet into seg
    nOrph = 0;
    if isfield(R,'track') && isfield(R.track,'all')
        for q = res.trkUse
            r = R.track.all{q};
            segOfQ = find(res.useIdx == res.trkSeg(q), 1);
            if isempty(segOfQ)
                col = [0.6 0.6 0.6];  nOrph = nOrph + 1;
            else
                col = cmapS(segOfQ,:);
            end
            plot(a2, r(:,1), r(:,2), '-', 'Color', col, 'LineWidth', 1.1);
        end
        title(a2, sprintf('%d valid tracklets on the dendrite (%d in no segment)', ...
            numel(res.trkUse), nOrph), 'FontSize', 9);
    else
        title(a2, 'tracklet paths need step 2 in the same session', 'FontSize', 9);
    end

    % ---- LEFT BOTTOM: is a segment's speed backed by enough tracklets? --------------
    % A segment's median speed becomes its frequency band via v/d, so a speed resting on
    % three tracklets is a band resting on three tracklets. Points left of the minTrk
    % line were rejected for exactly that reason.
    a3 = axes(parent,'Position',[0.09 0.09 0.36 0.35]);  hold(a3,'on');
    scatter(a3, res.nTrkSeg, res.segSpeed, 46, cmapS(1:numel(res.px),:), 'filled');
    for k = 1:numel(res.px)
        text(a3, res.nTrkSeg(k), res.segSpeed(k), sprintf('  s%d', res.useIdx(k)), ...
            'FontSize', 8, 'VerticalAlignment','middle');
    end
    xline(a3, cfg.segments.minTrk, 'r--', 'minTrk');
    box(a3,'off'); xlabel(a3,'tracklets in the segment','FontSize',8);
    ylabel(a3,'median speed (\mum/s)','FontSize',8);
    xl = xlim(a3); xlim(a3, [0 max(xl(2), cfg.segments.minTrk*1.3)]);
    title(a3, 'segment speed vs how much evidence it has', 'FontSize', 9);

    % ---- RIGHT BOTTOM: the segments themselves -------------------------------------
    a4 = axes(parent,'Position',[0.57 0.09 0.36 0.35]);
    L = zeros(size(res.densBW));
    for k = 1:numel(res.px), L(res.px{k}) = k; end
    h = imagesc(a4, L, [0 max(1,numel(res.px))]);
    set(h,'AlphaData', 0.25 + 0.75*(L>0));  axis(a4,'image','off');
    colormap(a4, [0.85 0.85 0.85; cmapS]);  hold(a4,'on');
    if isfield(data,'bvBW') && any(data.bvBW(:))
        contour(a4, data.bvBW, [0.5 0.5], 'k', 'LineWidth', 1);
    end
    for k = 1:numel(res.px)
        [rr,cc] = ind2sub(size(L), res.px{k});
        text(a4, mean(cc), mean(rr), sprintf('%d', res.useIdx(k)), ...
            'Color','k','FontWeight','bold','HorizontalAlignment','center');
    end
    title(a4, sprintf('SEGMENTS: %d used of %d patches', numel(res.px), res.nSegAll), ...
        'FontSize', 9);

% =====================================================================================
case 'regressors'
    % TOP: every regressor as a trace, strongest first, ID on the y axis. The scatters
    % that used to live here are in the browser, where a point can be clicked.
    nR = size(res.design,1);
    a1 = axes(parent,'Position',[0.10 0.58 0.85 0.37]);  hold(a1,'on');
    tAx = (1:size(res.design,2)) / cfg.common.fs;
    tW  = 1:min(numel(tAx), round(2*cfg.common.fs));   % 2 s: enough to judge an oscillator
    cmR = turbo(max(nR,2));
    for q = 1:nR
        y = res.design(q, tW);  y = 0.42 * y / max(max(abs(y)), eps);
        plot(a1, tAx(tW), y - (q-1), '-', 'Color', cmR(q,:), 'LineWidth', 0.7);
    end
    ylim(a1, [-(nR-1)-0.6, 0.6]);  xlim(a1, tAx(tW([1 end])));
    % YTick must be ASCENDING, and the traces run downwards, so both the positions and
    % the labels are flipped together -- flipping only one silently mislabels every row.
    step = max(1, ceil(nR/25));
    ids  = 1:step:nR;
    set(a1, 'YTick', -fliplr(ids-1), 'YTickLabel', fliplr(ids), 'FontSize', 7.5, ...
            'YTickMode','manual', 'YTickLabelMode','manual');
    box(a1,'off');  xlabel(a1,'time (s)','FontSize',8);
    ylabel(a1,'regressor ID  (1 = most variance explained)','FontSize',8);
    title(a1, sprintf('all %d regressors from %d segments, ordered by variance explained', ...
        nR, numel(res.useSeg)), 'FontSize', 9);

    % BOTTOM: what the design as a whole takes out of the movie
    % BOTTOM LEFT: the spatial modes, on the dendrite. Segments are disjoint, so every
    % segment's mode of a given rank shares one image; the popup steps through the ranks.
    nMode = max(res.modeOf);
    mStr = arrayfun(@(k) sprintf('mode %d', k), 1:nMode, 'uni', 0);
    popE = uicontrol(parent,'Style','popupmenu','String',mStr, ...
        'Units','normalized','Position',[0.05 0.475 0.13 0.035], ...
        'Callback',@(s2,~) drawEig(parent, res, data, get(s2,'Value')));
    set(popE,'UserData',res);
    drawEig(parent, res, data, 1);

    % BOTTOM RIGHT: what the design as a whole takes out
    pop = uicontrol(parent,'Style','popupmenu', ...
        'String',{'variance removed','fraction of variance removed'}, ...
        'Units','normalized','Position',[0.55 0.475 0.24 0.035], ...
        'Callback',@(s2,~) drawVar(parent, res, data, get(s2,'Value')));
    set(pop, 'UserData', res);
    drawVar(parent, res, data, 1);

    uicontrol(parent,'Style','pushbutton','String','Browse regressors one at a time', ...
        'Units','normalized','Position',[0.20 0.475 0.30 0.038],'FontSize',9, ...
        'Callback',@(~,~) openBrowser(res, cfg, data));

% =====================================================================================
case 'apply'
    % LEFT: the movie, three panels STACKED, one intensity scale. Stacked rather than
    % side by side because the frame is wide: laid out in a row the three panels shrink
    % to a strip each, while stacked they keep the full width.
    applyMovie(parent, res, data);

    % RIGHT: var(fitted) for the two fits, on the SAME upper limit so they can be
    % compared by eye. Both are measured on the same held-out frames -- the only footing
    % on which they are comparable, since a windowed fit removes more in-sample whatever
    % the data says.
    a1 = axes(parent,'Position',[0.55 0.55 0.42 0.38]);
    a2 = axes(parent,'Position',[0.55 0.07 0.42 0.38]);
    if ~isempty(res.split) && isfield(res.split,'varFitG')
        sp = res.split;
        hi = prctile([sp.varFitG(:); sp.varFitS(:)], 99);
        if ~isfinite(hi) || hi <= 0, hi = 1; end
        showFix(a1, sp.varFitG, hi, data, 'var(fitted) GLOBAL, held out');
        showFix(a2, sp.varFitS, hi, data, sprintf('var(fitted) SPLIT, held out\n%s', ...
            firstLine(sp.verdict)));
    else
        showFix(a1, res.varFit, [], data, 'var(fitted), global fit');
        axis(a2,'off');
        text(a2, 0.5, 0.5, sprintf(['split test not run\n(cfg.apply.split = false)\n\n' ...
            'turn it on to compare a time-varying\ncoefficient against this one']), ...
            'HorizontalAlignment','center', 'Color',[0.5 0.5 0.5]);
    end

% =====================================================================================
case 'compare'
    % The voltage traces themselves, before and after, as kymographs: one row per ROI of
    % data.ftprnt, time across. Both share ONE colour scale taken from the RAW panel --
    % scaling them separately would renormalise the corrected one and hide exactly the
    % change being looked for.
    nR  = size(res.dRaw,1);
    tAx = (1:size(res.dRaw,2)) / res.fs;
    sc  = res.scale(:);                       % per ROI, so depth does not dominate
    XR  = res.dRaw .* sc;
    XC  = res.dCor .* sc;
    cl  = prctile(abs(XR(:)), 99.5);
    if ~isfinite(cl) || cl <= 0, cl = 1; end
    tick = unique(round(linspace(1, nR, min(10, nR))));
    cmT  = gen_colormap([0 0.5 1; 1 1 1; 1 0 0], 256);

    XD  = res.dBld .* sc;                     % what was REMOVED: raw minus corrected

    axA = axes(parent,'Position',[0.08 0.755 0.86 0.185]);
    imagesc(axA, tAx, 1:nR, XR, [-cl cl]);  colormap(axA, cmT);
    set(axA,'YTick',tick,'XTickLabel',[]);
    ylabel(axA, {'mov\_res','ROI'}, 'FontSize',8);
    title(axA, sprintf('voltage traces . ftprnt, %d ROIs   [%s]', nR, res.normLbl), 'FontSize', 9);

    axB = axes(parent,'Position',[0.08 0.560 0.86 0.185]);
    imagesc(axB, tAx, 1:nR, XC, [-cl cl]);  colormap(axB, cmT);
    set(axB,'YTick',tick,'XTickLabel',[]);
    ylabel(axB, {'corrected','ROI'}, 'FontSize',8);
    title(axB, 'after the correction', 'FontSize', 9);

    % The residual on the SAME colour scale as the two above, which is the whole point:
    % scaled to its own range it would look dramatic however little was taken out.
    axD = axes(parent,'Position',[0.08 0.365 0.86 0.185]);
    imagesc(axD, tAx, 1:nR, XD, [-cl cl]);  colormap(axD, cmT);
    set(axD,'YTick',tick);
    ylabel(axD, {'removed','ROI'}, 'FontSize',8);
    xlabel(axD, 'time (s)', 'FontSize', 8);
    title(axD, 'residual = mov\_res - corrected   (same colour scale as above)', 'FontSize', 9);
    linkaxes([axA axB axD], 'x');

    % What was REMOVED from each ROI, painted back onto the dendrite. Read against the
    % segment outline: an ROI far from every segment that still loses variance is the
    % correction reaching somewhere it was not derived from.
    axC = axes(parent,'Position',[0.22 0.035 0.56 0.30]);
    vB  = res.varBld(:) ./ max(res.varRaw(:), eps) * 100;
    L   = roiLabel(data);
    M   = nan(res.sz);
    ok  = L > 0 & L <= numel(vB);
    M(ok) = vB(L(ok));
    hi = prctile(M(:), 99);  if ~isfinite(hi) || hi <= 0, hi = max([M(isfinite(M)); eps]); end
    hM = imagesc(axC, M, [0 hi]);  set(hM,'AlphaData', ~isnan(M));
    axis(axC,'image','off');  colormap(axC, turbo);
    cb = colorbar(axC);  cb.Label.String = '% of var(mov\_res) removed';
    hold(axC,'on');
    if isfield(res,'segMask') && any(res.segMask(:)) && ~all(res.segMask(:))
        contour(axC, res.segMask, [0.5 0.5], 'w', 'LineWidth', 1.3);
    end
    if isfield(data,'bvBW') && any(data.bvBW(:)) && ~all(data.bvBW(:))
        contour(axC, data.bvBW, [0.5 0.5], 'Color',[0.2 0.9 0.9], 'LineWidth', 0.8);
    end
    onS = res.onSeg(:);
    if any(onS) && any(~onS)
        ttl = sprintf(['variance removed per ROI, on the dendrite   |   %.2f%% median ' ...
            'on a segment, %.2f%% off'], median(vB(onS)), median(vB(~onS)));
    else
        ttl = 'variance removed per ROI, on the dendrite';
    end
    title(axC, sprintf('%s\n(white = regressor segments, cyan = bvMask)', ttl), 'FontSize', 9);

% =====================================================================================
otherwise
    ax = axes(parent,'Position',[0.05 0.05 0.9 0.9]);
    f = fieldnames(res);
    Ls = {sprintf('%s produced:', stage), ''};
    for k = 1:min(numel(f),16)
        v = res.(f{k});
        if isnumeric(v) || islogical(v)
            Ls{end+1} = sprintf('  %-14s %s %s', f{k}, class(v), mat2str(size(v))); %#ok<AGROW>
        elseif iscell(v)
            Ls{end+1} = sprintf('  %-14s cell %s', f{k}, mat2str(size(v)));         %#ok<AGROW>
        else
            Ls{end+1} = sprintf('  %-14s %s', f{k}, class(v));                      %#ok<AGROW>
        end
    end
    text(ax, 0.02, 0.98, Ls, 'VerticalAlignment','top', 'FontName','Courier', ...
        'FontSize', 9, 'Interpreter','none');
    axis(ax,'off');
end

% -------------------------------------------------------------------------------------
    function swapMovie(src)
        r = get(src,'UserData');
        if get(src,'Value') == 1, drawMovie(parent, r.S, 'S (band-passed)');
        else,                     drawMovie(parent, r.mov_nv, 'mov\_nv (input)'); end
    end
end

function applyMovie(parent, res, data)
% [ input | fitted | corrected ] concatenated into one frame, one intensity scale.
%
% Only a window is built: the correction is a linear projection, so any frame follows
% from C, and materialising three panels over a whole recording would be tens of GB for
% something nobody scrolls end to end.
delete(findobj(parent, 'Tag', 'flowMovie'));
nT  = size(data.mov_res, 3);
nQ  = min(600, nT);                       % three panels wide, so keep the window short
q   = 1:nQ;
sz  = res.sz;
I    = [ones(1,nQ); res.design(:,q) - mean(res.design(:,res.good), 2)];
fitQ = reshape(single(res.C) * single(I), sz(1), sz(2), nQ);
raw  = -single(data.mov_res(:,:,q));      % spike-positive, as everywhere else
fitD = -fitQ;
corD = raw - fitD;

cl = double(prctile(raw(:), [1 99]));
if ~(cl(2) > cl(1)), cl = double([min(raw(:)) max(raw(:))]); end
gap = 3;
sep = single(cl(2)) * ones(gap, sz(2), nQ, 'single');   % bright divider
M   = cat(1, raw, sep, fitD, sep, corD);                % STACKED

ax = flow_movieAxes(parent, [0.04 0.06 0.44 0.88], M, struct( ...
    'title', 'input / fitted / corrected   (one shared scale)', ...
    'cmap', gray(256), 'clim', cl));
% label each panel on the image itself, so no legend is needed
h = sz(1);
lbl = {'mov\_res', 'fitted', 'corrected'};
for k = 1:3
    text(ax, 6, (k-1)*(h+gap) + 10, lbl{k}, 'Color',[1 1 0.3], ...
        'HorizontalAlignment','left', 'FontWeight','bold', 'FontSize', 9);
end
set(ax, 'Tag', 'flowMovie');
u = ax.UserData;  set([u.sl u.ed], 'Tag', 'flowMovie');
end

function showFix(ax, M, hi, data, ttl)
% a variance map on a FIXED upper limit, so two of them can be compared by eye
if isfield(data,'anaMask') && isequal(size(data.anaMask), size(M)), M(~data.anaMask) = NaN; end
if isempty(hi)
    hi = prctile(M(:), 99);
    if ~isfinite(hi) || hi <= 0, hi = max([M(isfinite(M)); eps]); end
end
h = imagesc(ax, M, [0 hi]);  set(h, 'AlphaData', ~isnan(M));
axis(ax,'image','off');  colormap(ax, turbo);  colorbar(ax);  hold(ax,'on');
if isfield(data,'bvBW') && any(data.bvBW(:)) && ~all(data.bvBW(:)) ...
        && isequal(size(data.bvBW), size(M))
    contour(ax, data.bvBW, [0.5 0.5], 'w', 'LineWidth', 0.8);
end
title(ax, ttl, 'FontSize', 9);
end

function drawEig(parent, res, data, rank)
% Every segment's rank-th spatial mode, painted over the dendrite.
%
% Built as TRUECOLOR rather than an image plus an overlay, because one axes has one
% colormap: a grey backdrop and a cyclic phase map cannot both use it. The dendrite is
% converted to grey RGB and the mode blended in on top, weighted by |U| -- so the colour
% is the propagation phase and the opacity is how strongly the mode lives there.
delete(findobj(parent, 'Tag', 'flowEigAx'));
ax = axes(parent, 'Position', [0.05 0.055 0.42 0.39]);
sz = res.sz;

bg = zeros(sz);
if isfield(data,'neuronMap') && isequal(size(data.neuronMap), sz), bg = double(data.neuronMap); end
if ~any(bg(:)) && isfield(res,'anaMask'), bg = double(res.anaMask); end
hi = prctile(bg(:), 99.5);  if ~isfinite(hi) || hi <= 0, hi = max([bg(:); eps]); end
bg = min(max(bg/hi, 0), 1) * 0.55;                 % dim, so the modes read on top
RGB = repmat(bg, 1, 1, 3);

% one row per regressor of this rank; segments are disjoint so they cannot collide
rows = find(res.modeOf == rank);
seen = false(1, numel(res.useSeg));
cmapP = hsv(256);
nPaint = 0;
for r = rows(:)'
    k = res.segOf(r);
    if seen(k), continue, end                       % cos and sin share a mode
    seen(k) = true;
    E = reshape(res.eigImg(:,r), sz);
    m = ~isnan(E);
    if ~any(m(:)), continue, end
    nPaint = nPaint + 1;
    if isreal(E)
        val = (real(E)+1)/2;                        % signed map -> 0..1
        aW  = abs(real(E));
    else
        val = (angle(E)+pi)/(2*pi);                 % phase -> 0..1 on a cyclic map
        aW  = abs(E);
    end
    aW(~m) = 0;  aW = aW / max(aW(:));
    idx = max(1, min(256, round(val*255)+1));
    for ch = 1:3
        C = RGB(:,:,ch);
        Cq = cmapP(idx(m), ch);
        C(m) = (1-aW(m)).*C(m) + aW(m).*Cq;
        RGB(:,:,ch) = C;
    end
end
image(ax, min(max(RGB,0),1));  axis(ax,'image','off');  hold(ax,'on');
% imagesc/image run newplot, which CLEARS the axes Tag -- so it is set here, after the
% plotting call. Tagging at creation looks right and silently stops the next redraw from
% finding this axes to delete, leaving a new one stacked on top every time.
set(ax, 'Tag', 'flowEigAx');
if isfield(data,'bvBW') && ~isempty(data.bvBW) && any(data.bvBW(:)) && ~all(data.bvBW(:)) ...
        && isequal(size(data.bvBW), sz)
    contour(ax, data.bvBW, [0.5 0.5], 'w', 'LineWidth', 0.8);
end
if nPaint == 0
    title(ax, sprintf('no segment has a mode %d', rank), 'FontSize', 9);
else
    lbl = 'signed eigenvector';
    if ~isreal(res.eigImg), lbl = 'phase (hue) x |U| (opacity)'; end
    title(ax, sprintf('mode %d of %d segment(s) on the dendrite   [%s]', ...
        rank, nPaint, lbl), 'FontSize', 9);
end
end

function drawVar(parent, res, data, which)
% the removed-variance image, drawn large. Outside the footprint mask nothing is drawn:
% a variance figure there is not small, it is undefined.
delete(findobj(parent, 'Tag', 'flowVarAx'));
ax = axes(parent, 'Position', [0.53 0.055 0.44 0.39]);
if which == 1
    M = res.varFit;   ttl = 'variance removed by all regressors';
else
    M = res.dropFrac; ttl = 'fraction of the variance removed';
end
if isfield(res,'anaMask') && isequal(size(res.anaMask), size(M))
    M(~res.anaMask) = NaN;
end
hi = prctile(M(:), 99);
if ~isfinite(hi) || hi <= 0, hi = max([M(isfinite(M)); eps]); end
h = imagesc(ax, M, [0 hi]);  set(h, 'AlphaData', ~isnan(M));
axis(ax, 'image', 'off');  colormap(ax, turbo);  colorbar(ax);  hold(ax, 'on');
set(ax, 'Tag', 'flowVarAx');       % after the plot: newplot clears Tag
% contour needs a mask with an inside AND an outside; one that covers everything (a
% footprint that never reaches zero, say) has no boundary to draw
drawable = @(m) ~isempty(m) && any(m(:)) && ~all(m(:)) && isequal(size(m), size(M));
if isfield(data,'bvBW') && drawable(data.bvBW)
    contour(ax, data.bvBW, [0.5 0.5], 'w', 'LineWidth', 0.9);
end
onNeu = false(size(M));
if isfield(data,'neuronMap') && isequal(size(data.neuronMap), size(M))
    onNeu = data.neuronMap > 0;
    if drawable(onNeu)
        contour(ax, onNeu, [0.5 0.5], 'Color',[0.2 0.9 0.2], 'LineWidth', 0.7);
    end
    if isfield(data,'bvBW') && isequal(size(data.bvBW), size(M))
        onNeu = onNeu & ~data.bvBW;
    end
end
sub = '';
if any(onNeu(:)) && isfield(data,'bvBW') && any(data.bvBW(:))
    sub = sprintf('   |   %.1f%% on bvMask, %.1f%% on the dendrite', ...
        100*mean(res.dropFrac(data.bvBW & ~isnan(res.dropFrac))), ...
        100*mean(res.dropFrac(onNeu & ~isnan(res.dropFrac))));
end
title(ax, sprintf('%s%s   (white = bvMask, green = neuron)', ttl, sub), 'FontSize', 9);
end

function openBrowser(res, cfg, data)
% the per-regressor view: source trace, filtered trace, spatial mode, covariance map
segMask = cellfun(@(v) v, res.segPx, 'uni', 0);
stats = struct('varExpl', res.varExpl, 'rSTA', res.rVolt, 'rAlt', res.rVoltAlt, ...
    'oscScore', res.oscScore, 'score', res.score, ...
    'f0Hz', res.f0, 'segment', res.useSeg(res.segOf)', 'mode', res.modeOf);
browse_regressors(res.design, res.maps, res.sz, res.names, struct( ...
    'preTrace', res.srcTrace(res.segOf,:), ...
    'eigImg',   res.eigImg, ...
    'score',    res.score, 'scoreMin', cfg.regressors.scoreMin, ...
    'f0',       res.f0,    'seg', res.segOf, 'segMask', {segMask}, ...
    'refMask',  data.bvBW, 'mask', res.anaMask, 'stats', stats, ...
    'statXY',   {{'rSTA','varExpl'}}, 'statC', 'f0Hz', ...
    'tWin',     1:min(1000, size(res.design,2)), ...
    'figName',  'flow regressors'));
flow_theme(findobj(groot,'Type','figure','Name','flow regressors'));
end

function L = roiLabel(data)
% pixel -> ROI index. Footprints overlap, so a pixel is assigned to the ROI that weights
% it most; taking any() instead would let one ROI's tail overwrite another's centre.
L = zeros(0);
if ~isfield(data,'ftprnt') || isempty(data.ftprnt), return, end
F = double(data.ftprnt);
[mx, L] = max(F, [], 3);
L(mx <= 0) = 0;
end

function t = firstLine(s)
c = strsplit(char(s), {':','.'});
t = strtrim(c{1});
end

function M = bpMovie(R)
% the band-passed movie, if step 1 has been run
M = [];
if isfield(R,'bandpass') && isfield(R.bandpass,'S'), M = R.bandpass.S; end
end

function drawTracks(ax, P, k, trail)
% current detections bold, the preceding `trail` frames as a light line per tracklet
if isempty(P), return, end
old = P(P(:,3) >= k-trail & P(:,3) < k, :);
for q = unique(old(:,4))'
    r = old(old(:,4)==q, :);
    if size(r,1) < 2, continue, end
    plot(ax, r(:,1), r(:,2), '-', 'Color',[0.55 0.75 1.0], 'LineWidth', 1.0, ...
        'Tag','flowOverlay');
end
now = P(P(:,3) == k, :);
if ~isempty(now)
    plot(ax, now(:,1), now(:,2), 'o', 'MarkerSize', 6, 'LineWidth', 1.6, ...
        'MarkerEdgeColor',[0.05 0.25 0.85], 'Tag','flowOverlay');
end
end

function drawMovie(parent, mov, ttl)
% replace whatever movie panel is there with a new one
old = findobj(parent, 'Tag', 'flowMovie');
delete(old);
ax = flow_movieAxes(parent, [0.07 0.045 0.86 0.28], mov, ...
    struct('title', ttl, 'cmap', gray(256)));
set(ax, 'Tag', 'flowMovie');
u = ax.UserData;
set([u.sl u.ed], 'Tag', 'flowMovie');
end

function showMap(ax, M, data, ttl)
% outside the footprints a variance figure has no meaning, so it is not drawn
if isfield(data,'anaMask') && isequal(size(data.anaMask), size(M))
    M(~data.anaMask) = NaN;
end
hi = prctile(M(:), 99);
if ~isfinite(hi) || hi <= 0, hi = max([M(isfinite(M)); eps]); end
h = imagesc(ax, M, [0 hi]);  set(h, 'AlphaData', ~isnan(M));
axis(ax,'image','off');  colormap(ax, turbo);  colorbar(ax);
hold(ax,'on');
if isfield(data,'bvBW') && any(data.bvBW(:))
    contour(ax, data.bvBW, [0.5 0.5], 'w', 'LineWidth', 0.8);
end
title(ax, ttl, 'FontSize', 9);
end
