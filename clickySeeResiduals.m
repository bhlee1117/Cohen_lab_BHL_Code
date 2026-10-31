function out = clickySeeResiduals(movFit, movSrc, varargin)
% out = clickySeeResiduals(mov)
% out = clickySeeResiduals(movFit, movSrc)
% out = clickySeeResiduals(movFit, movSrc, 'fs', 1000, 'band', [40 200], ...)
%
% Pick regressors by hand, from ROIs drawn on one movie, and watch what removing them
% does to another.
%
% TWO MOVIES, TWO JOBS. movSrc is where the REGRESSOR comes from and movFit is what it is
% regressed OUT of. They are separated because the movie that shows a component most
% clearly is usually not the movie you want to correct: a spatially band-passed movie (S)
% makes a vessel obvious and has already had the voltage suppressed, while the movie you
% actually care about is the raw motion residual with everything still in it. Drawing on
% one and subtracting from the other is the whole point. Pass one movie and both jobs fall
% to it.
%
% WHAT IS ON SCREEN
%   top left      movFit, playing
%   top right     movFit after SeeResiduals against the regressor list -- same frame, same
%                 colour scale, so before and after can be read off against each other
%   middle left   a still of movSrc -- var, mean or median over time. DRAW ROIs HERE.
%   middle right  the trace of the ROI selected in the left-hand list, with its band-passed
%                 version and that version shifted 90 degrees overlaid
%   bottom left   left list: every trace available.  right list: the ones being regressed
%                 out.  Move them across with Add / Remove, then Apply.
%   bottom right  var(movFit - residual): WHERE the selected regressors actually took
%                 variance from. Bright away from the ROIs is the design reaching
%                 somewhere it should not, and no other panel shows that.
%
% IT RUNS UNTIL YOU PRESS DONE. Not until the first ROI is closed: drawpolygon does its
% own blocking with uiresume on this figure, so a plain uiwait here would return the
% moment a polygon was finished. Done (or closing the window) is the only thing that ends
% it, and the output is read at that point.
%
% WHY THE 90-DEGREE COPY IS OFFERED. A band-passed trace and the same trace shifted 90
% degrees are a quadrature pair, and a pair spans EVERY phase: a cos b sin reaches any
% A cos(phi - theta). So regressing the pair out removes that oscillation from a pixel
% whatever its phase offset from the ROI, while regressing the band-passed trace alone
% only removes the part that happens to be in phase with it. If you are chasing an
% oscillation, take the pair.
%
% THE TRACES ARE Z-SCORED before the fit. It does not change the residual -- least squares
% is invariant to the scale of a regressor -- but it conditions the design, which matters
% once several ROIs overlap and their traces are correlated.
%
% INPUT
%   movFit  [nrow x ncol x nT]  the movie that gets corrected
%   movSrc  [nrow x ncol x nT]  the movie the traces are read from (default: movFit)
%   options
%     'fs'      Hz, for the band-pass and the time axis          (default 1000)
%     'band'    'auto' (default) detects a band PER ROI from its own trace: the PSD peak
%               inside 'search', opened +/- 'bandFrac' around it. One band for the whole
%               window is only right if every ROI carries the same rhythm, and two vessels
%               at different flow speeds do not -- a band centred on one turns the other's
%               row into filtered noise. Pass [lo hi] Hz to fix one band for all ROIs;
%               typing in the band boxes does the same.
%     'search'  [lo hi] Hz the peak is looked for in    (default [20, min(400, fs/2-1)])
%     'bandFrac' half-width as a fraction of the peak            (default 0.5)
%               Deliberately generous: a passage rate wanders, and a narrow band centred
%               on its mean misses most of it.
%     'stat'    'var' | 'mean' | 'median', the still on the right (default 'var')
%     'clim'    display limits for the movies       (default: 1st-99th percentile of movFit)
%     'mask'    logical [nrow x ncol] drawn as a contour on the still, for reference
%     'auto'    recompute the residual on every change           (default true when the
%               movie is small enough to be quick, false otherwise)
%
% OUTPUT (on closing the window)
%   out.traces   [nTrace x nT]  every trace built, z-scored
%   out.names    {1 x nTrace}   what each one is
%   out.regIdx   indices into .traces of the ones selected as regressors
%   out.reg      [nSel x nT]    just those, so SeeResiduals(movFit, out.reg) reproduces
%                               the bottom-right movie exactly
%   out.roi      {1 x nROI}     logical masks
%   out.roiXY    {1 x nROI}     the polygon vertices
%   out.fs .band .stat
%
% 2026.08.30, Byung Hun Lee

% ---------------------------------------------------------------------------- input ---
if nargin < 2 || isempty(movSrc), movSrc = movFit; end
if ~isnumeric(movSrc) && ~islogical(movSrc)      % called as (mov,'opt',val,...)
    varargin = [{movSrc} varargin];  movSrc = movFit;
end
p = struct('fs',1000, 'band','auto', 'search',[], 'bandfrac',0.5, 'stat','var', ...
           'clim',[], 'mask',[], 'auto',[]);
for k = 1:2:numel(varargin)
    f = lower(varargin{k});
    assert(isfield(p,f), 'Unknown option ''%s''.', varargin{k});
    p.(f) = varargin{k+1};
end
assert(ndims(movFit) == 3, 'movFit must be [nrow x ncol x nT].');
assert(isequal(size(movFit), size(movSrc)), ...
    'The two movies must be the same size: movFit is %s, movSrc is %s.', ...
    mat2str(size(movFit)), mat2str(size(movSrc)));
[nrow, ncol, nT] = size(movFit);
if isempty(p.auto), p.auto = numel(movFit) <= 2e8; end
if isempty(p.search), p.search = [20, min(400, p.fs/2 - 1)]; end
assert(p.search(1) > 0 && p.search(2) > p.search(1) && p.search(2) < p.fs/2, ...
    'search must be [lo hi] with 0 < lo < hi < fs/2 (fs = %g).', p.fs);
if ischar(p.band) || isstring(p.band)
    assert(strcmpi(p.band, 'auto'), 'band must be [lo hi] or ''auto''.');
    autoBand0 = true;   band0 = p.search;
else
    assert(numel(p.band) == 2 && p.band(1) < p.band(2) && p.band(2) < p.fs/2, ...
        'band must be [lo hi] with hi < fs/2 (fs = %g), or ''auto''.', p.fs);
    autoBand0 = false;  band0 = p.band(:)';
end

S = struct();
S.movFit = movFit;  S.movSrc = movSrc;
S.nrow = nrow;  S.ncol = ncol;  S.nT = nT;
S.fs = p.fs;  S.band = band0;  S.mask = p.mask;  S.auto = p.auto;
S.autoBand = autoBand0;  S.search = p.search(:)';  S.bandFrac = p.bandfrac;
S.bandOf = zeros(0,2);  S.fPeak = [];    % per ROI: its own band, and the peak it came from
S.traces = zeros(0, nT);  S.names = {};  S.roi = {};  S.roiXY = {};  S.ofROI = [];
S.regIdx = [];  S.movRes = [];  S.frame = 1;  S.playing = false;
S.varFit0 = [];                 % cached: mean var of movFit, computed once on first Apply
S.popStatStr = lower(p.stat);   % set here too, so the output is defined however we exit
if isempty(p.clim)
    smp = double(movFit(:, :, unique(round(linspace(1, nT, min(nT,40))))));
    S.clim = prctile(smp(:), [1 99]);
    if ~(S.clim(2) > S.clim(1)), S.clim = [min(smp(:)) max(smp(:)) + eps]; end
else
    S.clim = p.clim;
end

% --------------------------------------------------------------------------- figure ---
% Opens at 4:3 and is then left alone -- a normalized Position is a fraction of the
% SCREEN, so without this the window's shape would follow whatever display it opened on
% (wide and short on a 21:9 monitor). Only the opening size is set; resizing afterwards is
% the user's business.
scr = get(groot, 'ScreenSize');
figW = min(1400, 0.92*scr(3));
figH = figW * 3/4;
if figH > 0.92*scr(4)
    figH = 0.92*scr(4);  figW = figH * 4/3;
end
fh = figure('Name','clickySeeResiduals','NumberTitle','off','Color',[0.94 0.94 0.94], ...
    'Units','pixels', ...
    'Position',[scr(1) + (scr(3)-figW)/2, scr(2) + (scr(4)-figH)/2, figW, figH], ...
    'CloseRequestFcn',@(~,~) onClose());
% R2025a+ figures follow the DESKTOP theme, so on a dark desktop this comes up dark and
% the listboxes, edit boxes and axes go with it -- and the theme is not something the
% plotting code below can override, because MATLAB refuses a groot default for it. Set it
% per figure. flow_theme does exactly this and is the one place in the repo that knows
% how; inlined as a fallback so this file still works on its own if lifted out.
if ~isempty(which('flow_theme'))
    flow_theme(fh);
elseif isprop(fh, 'Theme')
    try, set(fh, 'Theme', 'light'); catch, end   %#ok<NOCOM>
end

% Three bands: movies on top, traces and residual in the middle, the two lists across the
% The two MOVIES sit side by side on the top row so before and after can be compared
% frame by frame without the eye travelling. Below them, the still you draw on and the
% trace that comes out of it. Bottom left the lists, bottom right the map of what was
% actually taken -- which is the answer to "did that do what I wanted", and belongs
% next to the controls that decide it.
S.axFit = axes(fh,'Units','normalized','Position',[0.050 0.680 0.40 0.25]);
S.axRes = axes(fh,'Units','normalized','Position',[0.550 0.680 0.40 0.25]);
S.axSt  = axes(fh,'Units','normalized','Position',[0.050 0.380 0.40 0.25]);
S.axTr  = axes(fh,'Units','normalized','Position',[0.560 0.400 0.37 0.21]);
S.axRem = axes(fh,'Units','normalized','Position',[0.550 0.055 0.40 0.26]);

S.hFit = imagesc(S.axFit, movFit(:,:,1), S.clim);
axis(S.axFit,'image','off'); colormap(S.axFit, gray(256));
title(S.axFit, 'movFit  (input)', 'FontSize', 10);
S.hRes = imagesc(S.axRes, movFit(:,:,1), S.clim);
axis(S.axRes,'image','off'); colormap(S.axRes, gray(256));
title(S.axRes, 'residual  (nothing selected yet)', 'FontSize', 10);
S.hRem = imagesc(S.axRem, zeros(nrow, ncol));
axis(S.axRem,'image','off'); colormap(S.axRem, turbo);
S.cbRem = colorbar(S.axRem);
title(S.axRem, 'var(removed)  -- nothing selected', 'FontSize', 10);

% ---- controls ----
uicontrol(fh,'Style','text','String','still:','Units','normalized', ...
    'Position',[0.655 0.955 0.04 0.025],'HorizontalAlignment','right', ...
    'BackgroundColor',get(fh,'Color'));
S.popStat = uicontrol(fh,'Style','popupmenu','String',{'var','mean','median'}, ...
    'Units','normalized','Position',[0.70 0.955 0.09 0.028], ...
    'Callback',@(~,~) drawStill());
if ~strcmpi(p.stat,'var')
    S.popStat.Value = find(strcmpi({'var','mean','median'}, p.stat), 1);
end
uicontrol(fh,'Style','pushbutton','String','Draw ROI','Units','normalized', ...
    'Position',[0.80 0.955 0.085 0.028],'FontWeight','bold','Callback',@(~,~) addROI());
uicontrol(fh,'Style','pushbutton','String','Clear ROIs','Units','normalized', ...
    'Position',[0.893 0.955 0.09 0.028],'Callback',@(~,~) clearROIs());

S.btnPlay = uicontrol(fh,'Style','togglebutton','String','Play','Units','normalized', ...
    'Position',[0.04 0.955 0.06 0.028],'Callback',@(~,~) togglePlay());
S.sld = uicontrol(fh,'Style','slider','Min',1,'Max',max(nT,2),'Value',1, ...
    'SliderStep',[1/max(nT-1,1) 20/max(nT-1,1)],'Units','normalized', ...
    'Position',[0.105 0.955 0.24 0.028],'Callback',@(src,~) setFrame(round(src.Value)));
S.txtFrm = uicontrol(fh,'Style','text','String','1','Units','normalized', ...
    'Position',[0.35 0.955 0.09 0.025],'HorizontalAlignment','left', ...
    'BackgroundColor',get(fh,'Color'));

S.chkAutoB = uicontrol(fh,'Style','checkbox','String','auto band','Value',S.autoBand, ...
    'Units','normalized','Position',[0.445 0.955 0.072 0.028], ...
    'BackgroundColor',get(fh,'Color'),'Callback',@(src,~) setAutoBand(src.Value));
S.edLo = uicontrol(fh,'Style','edit','String',num2str(band0(1)),'Units','normalized', ...
    'Position',[0.520 0.955 0.032 0.028],'Callback',@(~,~) setBand());
S.edHi = uicontrol(fh,'Style','edit','String',num2str(band0(2)),'Units','normalized', ...
    'Position',[0.554 0.955 0.032 0.028],'Callback',@(~,~) setBand());
uicontrol(fh,'Style','text','String','Hz','Units','normalized', ...
    'Position',[0.589 0.953 0.02 0.025],'HorizontalAlignment','left', ...
    'BackgroundColor',get(fh,'Color'));

% ---- the lists and their buttons, bottom left ----
uicontrol(fh,'Style','text','String','all traces','Units','normalized', ...
    'Position',[0.035 0.300 0.185 0.022],'BackgroundColor',get(fh,'Color'), ...
    'HorizontalAlignment','left','FontWeight','bold');
S.lstAll = uicontrol(fh,'Style','listbox','String',{},'Units','normalized', ...
    'Position',[0.035 0.135 0.185 0.165],'Max',2,'Callback',@(~,~) showTrace());
uicontrol(fh,'Style','pushbutton','String','Add  >>','Units','normalized', ...
    'Position',[0.228 0.230 0.078 0.038],'Callback',@(~,~) addReg());
uicontrol(fh,'Style','pushbutton','String','<<  Remove','Units','normalized', ...
    'Position',[0.228 0.185 0.078 0.038],'Callback',@(~,~) rmReg());
uicontrol(fh,'Style','text','String','regressors','Units','normalized', ...
    'Position',[0.315 0.300 0.185 0.022],'BackgroundColor',get(fh,'Color'), ...
    'HorizontalAlignment','left','FontWeight','bold');
S.lstReg = uicontrol(fh,'Style','listbox','String',{},'Units','normalized', ...
    'Position',[0.315 0.135 0.185 0.165],'Max',2,'Callback',@(~,~) showTrace());
uicontrol(fh,'Style','pushbutton','String','Apply / update residual','Units','normalized', ...
    'Position',[0.035 0.082 0.185 0.042],'FontWeight','bold','Callback',@(~,~) applyReg());
% Under the button it governs, not up in the toolbar beside "auto band" -- two controls
% both labelled "auto", each a screen away from what it does, is how the toolbar became
% unreadable. This one decides whether the residual recomputes on every Add / Remove or
% waits for the button; it starts off on a movie big enough for that to hurt.
S.chkAuto = uicontrol(fh,'Style','checkbox','String','auto apply','Value',S.auto, ...
    'Units','normalized','Position',[0.035 0.048 0.185 0.026], ...
    'BackgroundColor',get(fh,'Color'),'Callback',@(src,~) setAuto(src.Value));
S.status = uicontrol(fh,'Style','text','String','draw an ROI on the still, middle left', ...
    'Units','normalized','Position',[0.228 0.045 0.130 0.075], ...
    'HorizontalAlignment','left','BackgroundColor',get(fh,'Color'),'FontSize',8);
uicontrol(fh,'Style','pushbutton','String','Done','Units','normalized', ...
    'Position',[0.368 0.082 0.132 0.042],'FontWeight','bold','Callback',@(~,~) finish());

S.done = false;
guidata(fh, S);
drawStill();
setFrame(1);
showTrace();

% WHY THIS IS A LOOP AND NOT A BARE uiwait. drawpolygon implements its own blocking
% interaction with uiresume on this same figure, so finishing an ROI released the uiwait
% below and the function returned the moment the first polygon was closed. Anything else
% that blocks interactively (getline, waitfor on a child) does the same. A plain uiwait
% cannot tell "someone else resumed me" from "the user is finished", so the intent is
% recorded explicitly: nothing exits until S.done is set, by the Done button or by
% closing the window.
while isvalid(fh)
    uiwait(fh);
    if ~isvalid(fh), break, end
    S = guidata(fh);
    if S.done, break, end
end

if isvalid(fh), S = guidata(fh); delete(fh); end
out = struct('traces', S.traces, 'names', {S.names}, 'regIdx', S.regIdx, ...
    'reg', S.traces(S.regIdx,:), 'roi', {S.roi}, 'roiXY', {S.roiXY}, ...
    'ofROI', S.ofROI, 'fs', S.fs, 'band', S.band, 'bandOf', S.bandOf, ...
    'fPeak', S.fPeak, 'autoBand', S.autoBand, 'search', S.search, ...
    'bandFrac', S.bandFrac, 'stat', S.popStatStr);

% =====================================================================================
    function drawStill()
        S = guidata(fh);
        nm = {'var','mean','median'};
        st = nm{S.popStat.Value};
        S.popStatStr = st;
        switch st
            case 'var',    I = var(double(S.movSrc), 0, 3);
            case 'mean',   I = mean(double(S.movSrc), 3);
            case 'median', I = median(double(S.movSrc), 3);
        end
        S.still = I;
        cla(S.axSt);
        lim = prctile(I(:), [1 99]);
        if ~(lim(2) > lim(1)), lim = [min(I(:)) max(I(:)) + eps]; end
        imagesc(S.axSt, I, lim); axis(S.axSt,'image','off');
        colormap(S.axSt, gray(256)); colorbar(S.axSt);
        hold(S.axSt,'on');
        if ~isempty(S.mask)
            contour(S.axSt, double(S.mask), [0.5 0.5], 'w', 'LineWidth', 0.8);
        end
        guidata(fh, S);
        redrawROIs();
        title(S.axSt, sprintf('movSrc: %s over time   -- draw ROIs here', st), 'FontSize', 10);
    end

    function redrawROIs()
        S = guidata(fh);
        delete(findobj(S.axSt, 'Tag', 'roiline'));
        delete(findobj(S.axSt, 'Tag', 'roitext'));
        co = lines(max(numel(S.roi), 7));
        for k = 1:numel(S.roi)
            xy = S.roiXY{k};
            plot(S.axSt, [xy(:,1); xy(1,1)], [xy(:,2); xy(1,2)], '-', ...
                'Color', co(k,:), 'LineWidth', 1.5, 'Tag', 'roiline');
            text(S.axSt, mean(xy(:,1)), mean(xy(:,2)), sprintf('%d', k), ...
                'Color', co(k,:), 'FontWeight','bold', 'Tag','roitext', ...
                'HorizontalAlignment','center');
        end
    end

    function addROI()
        S = guidata(fh);
        set(S.status, 'String', 'click vertices, double-click / Enter to close');
        drawnow;
        xy = [];
        try
            h = drawpolygon(S.axSt, 'Color', [1 1 1], 'LineWidth', 1.2);
            if isvalid(h), xy = h.Position; delete(h); end
        catch
            % no Image Processing Toolbox: fall back to ginput
            [x, y] = getline(S.axSt); %#ok<GTLIN>
            xy = [x(:) y(:)];
        end
        if size(xy,1) < 3
            set(S.status, 'String', 'need at least 3 vertices');
            return
        end
        m = poly2mask(xy(:,1), xy(:,2), S.nrow, S.ncol);
        if ~any(m(:))
            set(S.status, 'String', 'that polygon contains no pixel');
            return
        end
        k = numel(S.roi) + 1;
        S.roi{k} = m;  S.roiXY{k} = xy;
        guidata(fh, S);
        addTracesFor(k);
        redrawROIs();
        S = guidata(fh);
        set(S.status, 'String', sprintf('ROI %d: %d px', k, nnz(m)));
    end

    function addTracesFor(k)
        % raw, band-passed, and the band-passed trace shifted 90 degrees. The pair spans
        % every phase, so together they take an oscillation out of a pixel whatever its
        % phase offset from this ROI -- the band-passed row alone only takes the part
        % already in phase with it.
        S = guidata(fh);
        idx = find(S.roi{k});
        % Chunked over frames, and only this ROI's pixels are ever pulled out. Reshaping
        % the whole movie to double first is the obvious way to write this and costs
        % 8 bytes x nPix x nT -- several GB for a real recording, for a trace that needs
        % a few hundred pixels.
        raw = zeros(1, S.nT);
        for c0 = 1:500:S.nT
            ci = c0:min(c0+499, S.nT);
            B  = reshape(S.movSrc(:,:,ci), S.nrow*S.ncol, numel(ci));
            raw(ci) = mean(double(B(idx, :)), 1);
        end
        % EACH ROI GETS ITS OWN BAND, measured from its own trace, when auto is on. One
        % band for the whole window is only right if every ROI carries the same rhythm;
        % two vessels at different flow speeds do not, and a band centred on one of them
        % turns the other's row into filtered noise. The peak is taken over `search`, and
        % the band opened +/- bandFrac around it -- generously, because a passage rate
        % wanders and a narrow band centred on its mean misses most of it.
        if S.autoBand
            % pwelch, not nanPSD: nanPSD parallelises over rows and prints progress, which
            % for the ONE row here means starting a parallel pool and a page of output to
            % spectrum-estimate a single trace. raw is a mean over pixels of a finite
            % movie, so there are no NaNs for nanPSD to handle anyway.
            nfft = min(1024, 2^floor(log2(max(S.nT,4)/2)));
            [Pw, Fw] = pwelch(raw(:) - mean(raw), hann(nfft), round(nfft/2), nfft, S.fs);
            Pw = Pw(:);  Fw = Fw(:);
            in = Fw >= S.search(1) & Fw <= S.search(2) & isfinite(Pw);
            if ~any(in)
                bnd = S.search;  fpk = NaN;
            else
                Pm = Pw;  Pm(~in) = -inf;
                [~, ip] = max(Pm);  fpk = Fw(ip);
                bnd = fpk * [1-S.bandFrac, 1+S.bandFrac];
                bnd = [max(bnd(1), 0.5), min(bnd(2), S.fs/2 - 1)];
                if ~all(isfinite(bnd)) || bnd(2) <= bnd(1)*1.05
                    bnd = S.search;  fpk = NaN;
                end
            end
        else
            bnd = S.band;  fpk = NaN;
        end
        bp = get_bandpass(raw, S.fs, bnd);
        bp(~isfinite(bp)) = 0;
        z  = hilbert(bp);
        new = [raw; real(z); imag(z)];
        new = (new - mean(new,2)) ./ max(std(new,0,2), eps);
        nm  = {sprintf('ROI%d raw', k), ...
               sprintf('ROI%d bp %.0f-%.0f Hz', k, bnd(1), bnd(2)), ...
               sprintf('ROI%d quad (90 deg)', k)};
        S.traces = [S.traces; new];
        S.names  = [S.names, nm];
        S.ofROI  = [S.ofROI, k k k];
        S.bandOf(k,:) = bnd;
        S.fPeak(k)    = fpk;
        guidata(fh, S);
        refreshLists();
    end

    function refreshLists()
        S = guidata(fh);
        set(S.lstAll, 'String', S.names, ...
            'Value', min(max(S.lstAll.Value,1), max(numel(S.names),1)));
        set(S.lstReg, 'String', S.names(S.regIdx), ...
            'Value', min(max(S.lstReg.Value,1), max(numel(S.regIdx),1)));
    end

    function addReg()
        S = guidata(fh);
        sel = S.lstAll.Value;
        sel = sel(sel >= 1 & sel <= numel(S.names));
        S.regIdx = unique([S.regIdx(:); sel(:)], 'stable')';
        guidata(fh, S);
        refreshLists();
        if S.auto, applyReg(); else, markStale(); end
    end

    function rmReg()
        S = guidata(fh);
        sel = S.lstReg.Value;
        sel = sel(sel >= 1 & sel <= numel(S.regIdx));
        S.regIdx(sel) = [];
        guidata(fh, S);
        refreshLists();
        if S.auto, applyReg(); else, markStale(); end
    end

    function clearROIs()
        S = guidata(fh);
        S.roi = {};  S.roiXY = {};  S.traces = zeros(0, S.nT);
        S.names = {};  S.ofROI = [];  S.regIdx = [];  S.movRes = [];
        S.bandOf = zeros(0,2);  S.fPeak = [];
        guidata(fh, S);
        refreshLists();  redrawROIs();  showTrace();
        set(S.hRes, 'CData', S.movFit(:,:,S.frame));
        title(S.axRes, 'residual  (nothing selected)', 'FontSize', 10);
        set(S.hRem, 'CData', zeros(S.nrow, S.ncol));
        title(S.axRem, 'var(removed)  -- nothing selected', 'FontSize', 10);
    end

    function setBand()
        S = guidata(fh);
        lo = str2double(S.edLo.String);  hi = str2double(S.edHi.String);
        if ~isfinite(lo) || ~isfinite(hi) || lo <= 0 || hi <= lo || hi >= S.fs/2
            set(S.status,'String', sprintf('band must be 0 < lo < hi < %g', S.fs/2));
            S.edLo.String = num2str(S.band(1));  S.edHi.String = num2str(S.band(2));
            return
        end
        S.band = [lo hi];
        if S.autoBand           % typing a band means you want that band, not a detected one
            S.autoBand = false;  S.chkAutoB.Value = 0;
            set(S.status,'String','auto band off: typing a band switches to manual');
        end
        guidata(fh, S);
        rebuildTraces();
    end

    function setAutoBand(v)
        S = guidata(fh);
        S.autoBand = logical(v);
        guidata(fh, S);
        rebuildTraces();
    end

    function rebuildTraces()
        % Every trace comes from its ROI and the band in force, so changing the band
        % rebuilds them all rather than trying to patch the list in place. The regressor
        % selection is dropped: the rows it pointed at no longer exist.
        S = guidata(fh);
        S.traces = zeros(0, S.nT);  S.names = {};  S.ofROI = [];  S.regIdx = [];
        S.bandOf = zeros(0,2);      S.fPeak = [];
        guidata(fh, S);
        for k = 1:numel(S.roi), addTracesFor(k); end
        showTrace();
        S = guidata(fh);
        if ~isempty(S.roi), markStale(); end
    end

    function setAuto(v)
        S = guidata(fh);  S.auto = logical(v);  guidata(fh, S);
    end

    function showTrace()
        S = guidata(fh);
        cla(S.axTr);
        if isempty(S.traces)
            title(S.axTr, 'no ROI yet', 'FontSize', 10);  return
        end
        i = S.lstAll.Value(1);
        k = S.ofROI(i);
        ii = find(S.ofROI == k);              % raw / bp / quad of the same ROI
        bnd = S.bandOf(k,:);
        tt = (0:S.nT-1) / S.fs;
        hold(S.axTr,'on');
        plot(S.axTr, tt, S.traces(ii(1),:), '-', 'Color',[0.25 0.25 0.25], ...
            'LineWidth', 0.7, 'DisplayName','raw');
        plot(S.axTr, tt, S.traces(ii(2),:), 'r-', 'LineWidth', 1.1, ...
            'DisplayName', sprintf('band-pass %.0f-%.0f Hz', bnd(1), bnd(2)));
        plot(S.axTr, tt, S.traces(ii(3),:), 'b-', 'LineWidth', 1.0, ...
            'DisplayName','90 deg shifted');
        box(S.axTr,'off');  xlabel(S.axTr,'time (s)');  ylabel(S.axTr,'z-score');
        xlim(S.axTr, [0 min(tt(end), 2)]);
        legend(S.axTr, 'Location','northeast', 'FontSize', 8);
        lbl = {'raw','bp','quad'};
        inReg = lbl(ismember(ii, S.regIdx));
        if isempty(inReg), inRegStr = 'none'; else, inRegStr = strjoin(inReg, ', '); end
        % The edit boxes always show the band actually in force for the ROI on screen,
        % detected or typed -- a box showing something the trace was not filtered with
        % would be worse than no box at all.
        S.edLo.String = sprintf('%.0f', bnd(1));
        S.edHi.String = sprintf('%.0f', bnd(2));
        if S.autoBand && isfinite(S.fPeak(k))
            pkStr = sprintf(', peak %.0f Hz (auto)', S.fPeak(k));
        elseif S.autoBand
            pkStr = ', no peak found -> search band';
        else
            pkStr = ' (manual)';
        end
        title(S.axTr, sprintf('ROI %d (%d px)%s   in regressors: %s', k, nnz(S.roi{k}), ...
            pkStr, inRegStr), 'FontSize', 10);
    end

    function markStale()
        S = guidata(fh);
        title(S.axRes, sprintf('residual STALE -- %d regressors selected, press Apply', ...
            numel(S.regIdx)), 'FontSize', 10, 'Color', [0.75 0.35 0]);
    end

    function applyReg()
        S = guidata(fh);
        if isempty(S.regIdx)
            S.movRes = [];  guidata(fh, S);
            set(S.hRes, 'CData', S.movFit(:,:,S.frame));
            title(S.axRes, 'residual  (nothing selected)', 'FontSize', 10, 'Color','k');
            set(S.hRem, 'CData', zeros(S.nrow, S.ncol));
            title(S.axRem, 'var(removed)  -- nothing selected', 'FontSize', 10, 'Color','k');
            return
        end
        set(S.status,'String','regressing ...');  drawnow;
        t0 = tic;
        S.movRes = SeeResiduals(S.movFit, S.traces(S.regIdx,:));
        if isempty(S.varFit0)                     % the input's variance never changes
            S.varFit0 = mean(var(double(S.movFit), 0, 3), 'all');
        end

        % var of what was TAKEN, per pixel: movFit - movRes. Accumulated in frame chunks
        % rather than differencing the two movies whole -- that would be two more full
        % double copies, and this map is the one thing here that has to be trustworthy on
        % a real recording. Where it is bright is where the regressors actually acted;
        % brightness away from the ROIs is the design reaching somewhere it should not.
        npx = S.nrow * S.ncol;
        s1 = zeros(npx,1);  s2 = zeros(npx,1);
        for c0 = 1:500:S.nT
            ci = c0:min(c0+499, S.nT);
            D  = double(reshape(S.movFit(:,:,ci), npx, numel(ci))) ...
               - double(reshape(S.movRes(:,:,ci), npx, numel(ci)));
            s1 = s1 + sum(D, 2);  s2 = s2 + sum(D.^2, 2);
        end
        varRem = reshape((s2 - s1.^2/S.nT) / max(S.nT-1,1), S.nrow, S.ncol);
        varRem = max(varRem, 0);
        guidata(fh, S);

        v1 = mean(var(double(S.movRes), 0, 3), 'all');
        set(S.hRes, 'CData', S.movRes(:,:,S.frame));
        title(S.axRes, sprintf('residual: %d regressors, %.1f%% of the variance removed', ...
            numel(S.regIdx), 100*(1 - v1/max(S.varFit0, eps))), 'FontSize', 10, 'Color','k');
        set(S.hRem, 'CData', varRem);
        hi = prctile(varRem(:), 99);
        if ~(hi > 0), hi = max(varRem(:)) + eps; end
        clim(S.axRem, [0 hi]);
        title(S.axRem, sprintf('var(removed)   peak %.3g, 99th pct %.3g', ...
            max(varRem(:)), hi), 'FontSize', 10, 'Color','k');
        set(S.status,'String', sprintf('done in %.1f s', toc(t0)));
    end

    function setFrame(fr)
        S = guidata(fh);
        fr = max(1, min(S.nT, round(fr)));
        S.frame = fr;  guidata(fh, S);
        set(S.hFit, 'CData', S.movFit(:,:,fr));
        if isempty(S.movRes)
            set(S.hRes, 'CData', S.movFit(:,:,fr));
        else
            set(S.hRes, 'CData', S.movRes(:,:,fr));
        end
        set(S.sld, 'Value', fr);
        set(S.txtFrm, 'String', sprintf('frame %d / %d   (%.3f s)', fr, S.nT, (fr-1)/S.fs));
    end

    function togglePlay()
        S = guidata(fh);
        S.playing = logical(S.btnPlay.Value);
        S.btnPlay.String = merge_(S.playing, 'Stop', 'Play');
        guidata(fh, S);
        while true
            if ~isvalid(fh), return, end
            S = guidata(fh);
            if ~S.playing, break, end
            setFrame(mod(S.frame, S.nT) + 1);
            drawnow limitrate;
        end
    end

    function finish()
        S = guidata(fh);  S.done = true;  S.playing = false;  guidata(fh, S);
        uiresume(fh);
    end

    function onClose()
        % Closing the window means the same thing as pressing Done -- otherwise the wait
        % loop above would spin on an invalid figure.
        S = guidata(fh);  S.done = true;  S.playing = false;  guidata(fh, S);
        if isvalid(fh), uiresume(fh); end
    end
end

% -------------------------------------------------------------------------------------
function v = merge_(c, a, b)
if c, v = a; else, v = b; end
end
