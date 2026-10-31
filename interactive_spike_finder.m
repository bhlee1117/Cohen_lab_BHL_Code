function [spikeIdx, threshold1, threshold2, zFloor, detWin] = interactive_spike_finder(trace, thr0, zTrace, zFloor0, detFn, detWin0)
% interactive_spike_finder Visual spike detector with adjustable thresholds.
%   [spikeIdx, threshold1, threshold2] = interactive_spike_finder(trace)
%   [spikeIdx, threshold1, threshold2] = interactive_spike_finder(trace, [thr1 thr2])  % initial thresholds
%   [spikeIdx, threshold1, threshold2, zFloor] = interactive_spike_finder(trace, thr0, zTrace, zFloor0)
%   [..., detWin] = interactive_spike_finder(trace, thr0, zTrace, zFloor0, detFn, detWin0)
%
%   INPUT:
%       trace   : 1D numeric array (spikes positive-going)
%       thr0    : initial [threshold1 threshold2] (default [5 2])
%       zTrace  : optional noise-normalised trace before any gain correction (e.g. before
%                 dividing by a bleaching envelope). A spike must also exceed the noise floor
%                 on it: zTrace(spike) > zFloor. -trace detections use -zTrace > zFloor.
%       zFloor0 : initial noise floor (default 4)
%       detFn   : optional @(win) -> detection trace, rebuilt from the 'Detection HP' box
%                 (e.g. trace minus its win-frame running median, then a matched filter);
%                 when given, the first input is ignored and detFn(detWin0) is shown
%       detWin0 : initial window for detFn (frames)
%   OUTPUT:
%       spikeIdx   : indices of detected spikes
%       threshold1 : final threshold1 value
%       threshold2 : final threshold2 value
%       zFloor     : final noise floor (NaN without zTrace)
%       detWin     : final detection window (NaN without detFn)
%
%   Right panels: find_spike_bh over a grid of threshold1 (2:0.5:6) x threshold2 (0.5:0.25:2.5)
%   on the trace and on the sign-flipped trace (with the noise floor, when given). Spikes are
%   positive-going, so detections on -trace estimate the false positives at each setting:
%       top    : # detections on +trace
%       bottom : false-positive estimate, # on -trace / # on +trace (%)
%   Click a cell to use that threshold pair. Blue x on the trace = -trace detections at the
%   current thresholds.

    % Default values
    threshold1 = 5;
    threshold2 = 2;
    if nargin > 1 && numel(thr0) == 2, threshold1 = thr0(1); threshold2 = thr0(2); end
    useFloor = nargin > 2 && ~isempty(zTrace);
    zFloor = NaN;
    if useFloor
        zTrace = double(zTrace(:)'); zTrace(~isfinite(zTrace)) = 0;
        zFloor = 4; if nargin > 3 && ~isempty(zFloor0), zFloor = zFloor0; end
    end
    th1Grid = 2:0.5:6;
    th2Grid = 0.5:0.25:2.5;
    useDet = nargin > 4 && ~isempty(detFn);
    detWin = NaN;
    if useDet
        detWin = 50; if nargin > 5 && ~isempty(detWin0), detWin = detWin0; end
        trace = detFn(detWin);
    end

    trace = double(trace(:)');
    trace(~isfinite(trace)) = 0;

    % --- candidates for the threshold scan (exactly find_spike_bh: local maxima above
    %     threshold1 whose 7-pt high-passed value exceeds threshold2), found once per trace
    [pkP, prP, locP, pkN, prN, locN] = deal([]);
    findCandidates();
    scanFloor = -Inf;

    % Create figure
    f = figure('Name','Interactive Spike Finder','NumberTitle','off','Position',[100 150 1400 600]);
    try, theme(f, 'light'); catch, end

    % Trace axes
    ax = axes('Parent',f,'Position',[0.05 0.25 0.55 0.68]);
    hTrace = plot(ax, trace, 'k'); hold(ax,'on'); grid(ax,'on');
    plotSpikes = plot(ax, nan, nan, 'ro', 'MarkerSize', 6);
    plotNeg    = plot(ax, nan, nan, 'bx', 'MarkerSize', 6);
    ylabel(ax,'Amplitude'); xlabel(ax,'Time');
    legend(ax, {'Trace','Spikes','Detections on -trace'}, 'Location','northeast');

    % Scan heatmaps
    axN = axes('Parent',f,'Position',[0.68 0.58 0.27 0.34]);
    axF = axes('Parent',f,'Position',[0.68 0.13 0.27 0.34]);
    selN = gobjects(1); selF = gobjects(1);

    % UI: Threshold 1
    uicontrol('Style','text','String','Threshold1:',...
              'Units','normalized','Position',[0.03 0.1 0.06 0.05],'HorizontalAlignment','right');
    th1Edit = uicontrol('Style','edit','String',num2str(threshold1),...
              'Units','normalized','Position',[0.10 0.1 0.05 0.05],'Callback',@updateThresholds);

    % UI: Threshold 2
    uicontrol('Style','text','String','Threshold2:',...
              'Units','normalized','Position',[0.16 0.1 0.06 0.05],'HorizontalAlignment','right');
    th2Edit = uicontrol('Style','edit','String',num2str(threshold2),...
              'Units','normalized','Position',[0.23 0.1 0.05 0.05],'Callback',@updateThresholds);

    % UI: noise floor (only with zTrace)
    uicontrol('Style','text','String','Noise floor:',...
              'Units','normalized','Position',[0.29 0.1 0.065 0.05],'HorizontalAlignment','right');
    flEdit = uicontrol('Style','edit','String',num2str(zFloor),...
              'Units','normalized','Position',[0.36 0.1 0.05 0.05],'Callback',@updateThresholds);
    if ~useFloor, set(flEdit, 'Enable','off', 'String','-'); end

    % UI: detection high-pass window (only with detFn)
    uicontrol('Style','text','String','Detection HP (frames):',...
              'Units','normalized','Position',[0.03 0.165 0.1 0.05],'HorizontalAlignment','right');
    dwEdit = uicontrol('Style','edit','String',num2str(detWin),...
              'Units','normalized','Position',[0.14 0.165 0.05 0.05],'Callback',@updateThresholds);
    if ~useDet, set(dwEdit, 'Enable','off', 'String','-'); end

    % Update button
    uicontrol('Style','pushbutton','String','Update',...
              'Units','normalized','Position',[0.43 0.1 0.07 0.05],...
              'Callback',@updateThresholds);

    % Done button
    uicontrol('Style','pushbutton','String','Done',...
              'Units','normalized','Position',[0.52 0.1 0.07 0.05],...
              'Callback','uiresume(gcbf);');

    % Count readout
    infoTxt = uicontrol('Style','text','String','', 'Units','normalized', ...
              'Position',[0.05 0.03 0.55 0.05],'HorizontalAlignment','left');

    % Initial detection
    updateThresholds();

    % Wait for user
    uiwait(f);

    % Return output
    threshold1 = str2double(th1Edit.String);
    threshold2 = str2double(th2Edit.String);
    if useFloor, zFloor = str2double(flEdit.String); end
    if useDet, detWin = str2double(dwEdit.String); end
    spikeIdx = get(plotSpikes, 'XData');

    % Close figure
    if isvalid(f)
        close(f);
    end

    % --- Nested update function
    function updateThresholds(~, ~)
        threshold1 = str2double(th1Edit.String);
        threshold2 = str2double(th2Edit.String);
        if useFloor, zFloor = str2double(flEdit.String); end
        if useDet                                          % rebuild the detection trace on a new window
            w = round(str2double(dwEdit.String));
            if isfinite(w) && w >= 0 && w ~= detWin
                detWin = w; set(dwEdit, 'String', num2str(w));
                trace = double(detFn(detWin)); trace = trace(:)'; trace(~isfinite(trace)) = 0;
                set(hTrace, 'YData', trace);
                findCandidates(); scanFloor = -Inf;
            end
        end
        if ~isequaln(scanFloor, zFloor), drawScans(); end
        spikes_frame = find(find_spike_bh(trace, threshold1, threshold2));
        neg_frame = find(find_spike_bh(-trace, threshold1, threshold2));
        if useFloor
            spikes_frame = spikes_frame(zTrace(spikes_frame) > zFloor);
            neg_frame = neg_frame(-zTrace(neg_frame) > zFloor);
        end
        set(plotSpikes, 'XData', spikes_frame, 'YData', trace(spikes_frame));
        set(plotNeg, 'XData', neg_frame, 'YData', trace(neg_frame));
        dstr = ''; if useDet, dstr = sprintf(', detection HP %d', detWin); end
        set(infoTxt, 'String', sprintf('threshold [%.2f %.2f]%s%s: %d spikes, %d on -trace (false positives ~%.1f%%)', ...
            threshold1, threshold2, floorStr(), dstr, numel(spikes_frame), numel(neg_frame), ...
            100*numel(neg_frame)/max(numel(spikes_frame),1)));
        % mark the nearest grid cell
        [~, i] = min(abs(th1Grid - threshold1)); [~, j] = min(abs(th2Grid - threshold2));
        set([selN selF], 'Position', [j-0.5 i-0.5 1 1]);
        drawnow;
    end

    function findCandidates()
        [pkP, prP, locP] = candidatePeaks(trace,  th1Grid(1));
        [pkN, prN, locN] = candidatePeaks(-trace, th1Grid(1));
    end

    function s = floorStr()
        if useFloor, s = sprintf(', noise floor %.2f', zFloor); else, s = ''; end
    end

    function drawScans()
        gP = true(size(pkP)); gN = true(size(pkN));
        if useFloor, gP = zTrace(locP) > zFloor; gN = -zTrace(locN) > zFloor; end
        [nPos, nNeg] = deal(zeros(numel(th1Grid), numel(th2Grid)));
        for i1 = 1:numel(th1Grid)
            for j1 = 1:numel(th2Grid)
                nPos(i1,j1) = sum(gP & pkP > th1Grid(i1) & prP > th2Grid(j1));
                nNeg(i1,j1) = sum(gN & pkN > th1Grid(i1) & prN > th2Grid(j1));
            end
        end
        fpr = 100 * nNeg ./ max(nPos, 1);
        drawScan(axN, nPos, '%d', ['# detections on +trace' floorStr()]);
        drawScan(axF, fpr, '%.0f', 'false positives: # on -trace / # on +trace (%)');
        clim(axF, [0 min(100, max(fpr(:))) + eps]);
        selN = rectangle(axN, 'Position',[0 0 1 1], 'EdgeColor','r', 'LineWidth',2);
        selF = rectangle(axF, 'Position',[0 0 1 1], 'EdgeColor','r', 'LineWidth',2);
        scanFloor = zFloor;
    end

    function drawScan(a, M, fmt, ttl)
        cla(a);
        imagesc(a, M, 'ButtonDownFcn', @onCellClick); axis(a, 'xy'); axis(a, 'tight'); hold(a, 'on');
        colormap(a, parula); colorbar(a);
        set(a, 'XTick', 1:numel(th2Grid), 'XTickLabel', compose('%g', th2Grid), ...
               'YTick', 1:numel(th1Grid), 'YTickLabel', compose('%g', th1Grid), 'FontSize', 8);
        xlabel(a, 'threshold2'); ylabel(a, 'threshold1'); title(a, ttl, 'FontSize', 9);
        for ii = 1:size(M,1)
            for jj = 1:size(M,2)
                text(a, jj, ii, sprintf(fmt, M(ii,jj)), 'HorizontalAlignment','center', ...
                    'FontSize', 7, 'Color', [0.5 0.5 0.5], 'HitTest','off');
            end
        end
    end

    function onCellClick(src, ~)
        cp = get(src.Parent, 'CurrentPoint');
        j = min(max(round(cp(1,1)), 1), numel(th2Grid));
        i = min(max(round(cp(1,2)), 1), numel(th1Grid));
        set(th1Edit, 'String', num2str(th1Grid(i)));
        set(th2Edit, 'String', num2str(th2Grid(j)));
        updateThresholds();
    end
end

function [pk, pr, loc] = candidatePeaks(x, thr)
% Local maxima above thr and their 7-pt high-passed value, as in find_spike_bh
[pk, loc] = findpeaks(x, 'MinPeakHeight', thr);
x2 = x - movmedian(x, 7, 2);
pr = x2(loc);
end
