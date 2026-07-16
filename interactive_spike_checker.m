function [sp, CS_tr] = interactive_spike_checker(trace, sp, CS_tr, viewWidth)
% INTERACTIVE_SPIKE_CHECKER  Manually curate a spike train and complex-spike
% regions on a single voltage trace.  Both sp and CS_tr are edited together
% and returned, so they stay consistent.
%
%   [sp, CS_tr] = interactive_spike_checker(trace, sp, CS_tr, viewWidth)
%
% INPUTS
%   trace     : 1 x T voltage trace (spike-positive)
%   sp        : 1 x T spike train (0/1).            Default: zeros(1,T)
%   CS_tr     : 1 x T complex-spike region (0/1).   Default: zeros(1,T)
%   viewWidth : # of frames shown at once.          Default: 2000
%
% OUTPUTS
%   sp    : edited spike train  (0/1 row vector)
%   CS_tr : edited complex-spike region (0/1 row vector)
%
% CONTROLS  (pick the click-action with the buttons)
%   Add spike     : click near a peak  -> adds a spike (snaps to local max)
%   Remove spike  : click near a spike -> removes the nearest spike
%   Mark CS       : click TWO points   -> sets that span as a complex spike
%   Unmark CS     : click TWO points   -> clears complex spike over that span
% Keys:  <-/-> jump to prev/next event (spike or CS start), shift+<-/-> pan,
%        up/down zoom out/in,  a/r/c/u pick mode,
%        esc cancels a pending CS click,  d or Done finishes.
% Black markers = simple spikes, red markers = spikes inside a CS region;
% CS spans are shaded pink.

    trace = trace(:)';  T = numel(trace);
    if nargin < 2 || isempty(sp),    sp    = zeros(1,T); end
    if nargin < 3 || isempty(CS_tr), CS_tr = zeros(1,T); end
    if nargin < 4 || isempty(viewWidth), viewWidth = 2000; end
    sp    = double(sp(:)')    > 0;
    CS_tr = double(CS_tr(:)') > 0;

    snapWin   = 5;      % add: search +/- this many frames for the local max
    removeWin = 10;     % remove: max distance (frames) to snap to a spike
    mode      = 'add';
    csAnchor  = NaN;    % pending first click for a CS span
    viewStart = 1;

    yl = [min(trace) max(trace)];
    yl = yl + [-1 1]*0.05*range(yl) + [-eps eps];

    % ---- figure / axes ----
    fig = figure('Name','Interactive Spike Checker','NumberTitle','off', ...
                 'Position',[150 150 1100 500],'KeyPressFcn',@onKey);
    ax  = axes('Parent',fig,'Position',[0.06 0.22 0.92 0.72]); hold(ax,'on'); box(ax,'on');
    hTrace  = plot(ax, 1:T, trace, 'k'); set(hTrace,'HitTest','off');
    hSS     = plot(ax, nan, nan, 'ko','MarkerFaceColor','k','MarkerSize',5); set(hSS,'HitTest','off');
    hCS     = plot(ax, nan, nan, 'ro','MarkerFaceColor','r','MarkerSize',6); set(hCS,'HitTest','off');
    hAnchor = plot(ax, [nan nan],yl,'b--','LineWidth',1.5); set(hAnchor,'HitTest','off');
    hCurEvt = plot(ax, [nan nan],yl,'-','Color',[1 0 1 0.6],'LineWidth',1.2); set(hCurEvt,'HitTest','off');
    csPatches = gobjects(0);
    ylim(ax, yl); xlabel(ax,'Frame'); ylabel(ax,'Amplitude');
    set(ax,'ButtonDownFcn',@onClick);

    % ---- mode buttons ----
    mkBtn(0.06,'Add spike',   @() setMode('add'));
    mkBtn(0.20,'Remove spike',@() setMode('remove'));
    mkBtn(0.34,'Mark CS',     @() setMode('cs'));
    mkBtn(0.48,'Unmark CS',   @() setMode('uncs'));
    uicontrol(fig,'Style','pushbutton','String','Done','Units','normalized', ...
        'Position',[0.88 0.03 0.10 0.08],'Callback','uiresume(gcbf)');

    updateView(); updateMarkers(); refreshTitle();
    uiwait(fig);

    if isvalid(fig), close(fig); end
    sp = double(sp); CS_tr = double(CS_tr);

    % ================= nested helpers =================
    function mkBtn(x,str,cb)
        uicontrol(fig,'Style','pushbutton','String',str,'Units','normalized', ...
            'Position',[x 0.03 0.12 0.08],'Callback',@(~,~) cb());
    end

    function setMode(m)
        mode = m; csAnchor = NaN; set(hAnchor,'XData',[nan nan]); refreshTitle();
    end

    function refreshTitle()
        pend = ''; if ~isnan(csAnchor), pend = '  [CS: click 2nd point]'; end
        title(ax, sprintf('Mode: %s   |   %d spikes, %d CS frames%s', ...
            mode, sum(sp), sum(CS_tr), pend));
    end

    function updateMarkers()
        ssIdx = find(sp & ~CS_tr);
        csIdx = find(sp &  CS_tr);
        set(hSS,'XData',ssIdx,'YData',trace(ssIdx));
        set(hCS,'XData',csIdx,'YData',trace(csIdx));
        delete(csPatches(isgraphics(csPatches))); csPatches = gobjects(0);
        b = bwlabel(CS_tr);
        for r = 1:max(b)
            fr = find(b==r); x0 = fr(1)-0.5; x1 = fr(end)+0.5;
            p = patch(ax,[x0 x1 x1 x0],[yl(1) yl(1) yl(2) yl(2)], ...
                      [1 0.82 0.82],'EdgeColor','none','FaceAlpha',0.5,'HitTest','off');
            uistack(p,'bottom'); csPatches(end+1) = p; %#ok<AGROW>
        end
        refreshTitle();
    end

    function updateView()
        viewStart = max(1, min(viewStart, max(1,T-viewWidth+1)));
        xlim(ax, [viewStart, min(T, viewStart+viewWidth-1)]);
    end

    function ev = getEvents()
        csStart = find(diff([0 CS_tr])==1);       % start frame of each CS span
        ev = unique([find(sp), csStart]);
    end

    function jumpEvent(dir)
        ev = getEvents();
        if isempty(ev), return; end
        cur = mean(xlim(ax));                      % current view center
        if dir > 0
            nxt = ev(find(ev > cur+0.5, 1, 'first'));
            if isempty(nxt), nxt = ev(end); end
        else
            nxt = ev(find(ev < cur-0.5, 1, 'last'));
            if isempty(nxt), nxt = ev(1); end
        end
        viewStart = round(nxt - viewWidth/2); updateView();
        set(hCurEvt,'XData',[nxt nxt]);           % highlight current event
    end

    function onClick(~,~)
        cp = get(ax,'CurrentPoint'); xc = round(cp(1,1)); xc = max(1,min(T,xc));
        switch mode
            case 'add'
                w = max(1,xc-snapWin):min(T,xc+snapWin);
                [~,mi] = max(trace(w)); pk = w(mi);
                sp(pk) = true; updateMarkers();
            case 'remove'
                spi = find(sp);
                if ~isempty(spi)
                    [d,mi] = min(abs(spi-xc));
                    if d <= removeWin, sp(spi(mi)) = false; end
                end
                updateMarkers();
            case {'cs','uncs'}
                if isnan(csAnchor)
                    csAnchor = xc; set(hAnchor,'XData',[xc xc]); refreshTitle();
                else
                    a = min(csAnchor,xc); b2 = max(csAnchor,xc);
                    CS_tr(a:b2) = strcmp(mode,'cs');
                    csAnchor = NaN; set(hAnchor,'XData',[nan nan]);
                    updateMarkers();
                end
        end
    end

    function onKey(~,evt)
        panMode = any(strcmp(evt.Modifier,'shift'));
        switch evt.Key
            case 'rightarrow'
                if panMode, viewStart = viewStart + round(viewWidth/2); updateView();
                else, jumpEvent(+1); end
            case 'leftarrow'
                if panMode, viewStart = viewStart - round(viewWidth/2); updateView();
                else, jumpEvent(-1); end
            case 'uparrow',    viewWidth = min(T, viewWidth*2);      updateView();
            case 'downarrow',  viewWidth = max(50, round(viewWidth/2)); updateView();
            case 'a', setMode('add');
            case 'r', setMode('remove');
            case 'c', setMode('cs');
            case 'u', setMode('uncs');
            case 'escape', csAnchor = NaN; set(hAnchor,'XData',[nan nan]); refreshTitle();
            case {'d','q','return'}, uiresume(fig);
        end
    end
end
