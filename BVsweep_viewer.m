function BVsweep_viewer(src)
% BVsweep_viewer  Interactive viewer for a bleach-dr blood-correction sweep.
%
%   BVsweep_viewer(sweepFile)   % path to R3_bleachdr_sweep_f<f>.mat
%   BVsweep_viewer(D)           % struct holding the sweep results
%
% One window with everything needed to compare parameter conditions:
%
%   left   : condition list (multi-select) + what to display + playback
%   centre : the movie, panels stacked in the selected order, one fixed
%            intensity scale taken from the ORIGINAL so that "less flicker"
%            reads as real removal rather than rescaling
%   right  : for the PRIMARY (first-selected) condition, its unsmoothed
%            correlation matrix and the difference from the original
%            sub_win-smoothed matrix; underneath, bAP survival along the
%            dendrite for EVERY selected condition overlaid
%   bottom : the metrics of the selected conditions
%
% Selecting different conditions updates the movie and every plot together,
% so a setting that removes more blood can be checked against the branch-
% signal cost in the same glance.
%
% Required fields (all produced by section (5) of the analysis script):
%   SW, okRun, Cmat_in, Cmat_raw_in, Cmat_raw_sw, Cmat_sw, dendaxisS, scr
% Optional:
%   prev_in, prev_corr, prev_blood   (without them the movie panel is disabled)
%   labelvecS, classNm, Dpair
%
% See also: get_corrTensor, moviefixsc

% ---------------------------------------------------------------- load ----
if ischar(src) || isstring(src)
    fprintf('[BVviewer] loading %s ...\n', src);
    D = load(char(src));
else
    D = src;
end
req = {'SW','okRun','Cmat_in','Cmat_raw_in','Cmat_raw_sw','scr'};
missing = req(~isfield(D, req));
assert(isempty(missing), 'BVsweep_viewer:missingFields', ...
       'missing field(s): %s', strjoin(missing, ', '));

S.SW      = D.SW;
S.runIdx  = find(D.okRun);
assert(~isempty(S.runIdx), 'BVsweep_viewer:noRuns', 'no combination completed');
S.Cmat_in     = D.Cmat_in;
S.Cmat_raw_in = D.Cmat_raw_in;
S.Cmat_raw_sw = D.Cmat_raw_sw;
S.Cmat_sw     = getfielddef(D, 'Cmat_sw', cell(1, numel(D.SW)));
S.subwin      = D.scr.sub_win;
S.dendaxis    = getfielddef(D, 'dendaxisS', []);
S.labelvec    = getfielddef(D, 'labelvecS', []);
S.classNm     = getfielddef(D, 'classNm', {});
S.nROI        = size(S.Cmat_in, 1);

S.hasMovie = isfield(D,'prev_in') && isfield(D,'prev_corr') && ~isempty(D.prev_in);
if S.hasMovie
    S.prev_in    = D.prev_in;
    S.prev_corr  = D.prev_corr;
    S.prev_blood = getfielddef(D, 'prev_blood', cell(1, numel(D.SW)));
    S.nFrame     = size(S.prev_in, 3);
    S.clim       = double(prctile(S.prev_in(:), [0.5 99.5]));
    if S.clim(2) <= S.clim(1), S.clim = double([min(S.prev_in(:)) max(S.prev_in(:))+eps]); end
else
    S.nFrame = 0;
end

% leaf order from the original smoothed matrix, shared by every matrix panel
D0 = 1 - S.Cmat_in;  D0 = (D0+D0')/2;  D0(1:S.nROI+1:end) = 0;  D0(isnan(D0)) = 1;
dv = squareform(D0);  Z = linkage(dv, 'average');
try
    S.ord = optimalleaforder(Z, dv);
catch
    fT = figure('Visible','off'); [~,~,S.ord] = dendrogram(Z, 0); close(fT);
end
S.triu = triu(true(S.nROI), 1);

try
    S.cmapCorr = gen_colormap([0 0.5 1; 1 1 1; 1 0 0], 256);
catch
    S.cmapCorr = [linspace(0,1,128)' linspace(0.5,1,128)' ones(128,1); ...
                  ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];
end

% ---------------------------------------------------------------- figure --
% The viewer is always LIGHT, whatever the session's theme is set to: the movie
% panels are greyscale and the correlation maps use a white-centred colormap,
% both of which only read correctly on a light ground.
fh = figure('Name','bleach-dr sweep viewer', 'NumberTitle','off', ...
            'Units','normalized', 'Position',[0.05 0.08 0.9 0.82], ...
            'WindowStyle','normal', 'Toolbar','figure', ...
            'CloseRequestFcn', @onClose);
try                                       % R2025a+; harmless (and unneeded) before
    theme(fh, 'light');
catch
end
CBG = [1 1 1];                           % figure / static-control background
CFG = [0.13 0.13 0.13];                  % text
CIN = [1 1 1];                           % editable + list controls
set(fh, 'Color', CBG);                   % after theme(): it resets Color itself
S.fh = fh;
S.col = struct('bg',CBG,'fg',CFG,'in',CIN);

condNames = arrayfun(@(k) sprintf('%d  %s', k, S.SW(k).tag), S.runIdx, 'UniformOutput', false);

% ---- left column: selection + controls ----
uicontrol(fh,'Style','text','Units','normalized','Position',[.008 .955 .17 .03], ...
    'String','CONDITIONS  (ctrl/shift = multi)','BackgroundColor',CBG, ...
    'ForegroundColor',CFG,'HorizontalAlignment','left','FontWeight','bold','FontSize',9);
S.hList = uicontrol(fh,'Style','listbox','Units','normalized', ...
    'Position',[.008 .70 .17 .25], 'String',condNames, 'Min',0,'Max',2, ...
    'Value',1, 'FontName','Menlo','FontSize',9, 'Callback',@onChange, ...
    'BackgroundColor',CIN,'ForegroundColor',CFG);

uicontrol(fh,'Style','text','Units','normalized','Position',[.008 .663 .08 .028], ...
    'String','Show','BackgroundColor',CBG,'ForegroundColor',CFG, ...
    'HorizontalAlignment','left','FontSize',9);
S.hShow = uicontrol(fh,'Style','popupmenu','Units','normalized', ...
    'Position',[.055 .665 .123 .03], 'FontSize',9, ...
    'String',{'Corrected','Blood','Corrected - Original','Original only'}, ...
    'Callback',@onChange,'BackgroundColor',CIN,'ForegroundColor',CFG);

S.hBlood = uicontrol(fh,'Style','checkbox','Units','normalized', ...
    'Position',[.008 .628 .17 .028], 'String','also show blood panel', ...
    'BackgroundColor',CBG,'ForegroundColor',CFG,'FontSize',9,'Callback',@onChange);
S.hOrig = uicontrol(fh,'Style','checkbox','Units','normalized', ...
    'Position',[.008 .598 .17 .028], 'String','ORIGINAL on top','Value',1, ...
    'BackgroundColor',CBG,'ForegroundColor',CFG,'FontSize',9,'Callback',@onChange);
S.hHoriz = uicontrol(fh,'Style','checkbox','Units','normalized', ...
    'Position',[.008 .568 .17 .028], 'String','tile horizontally', ...
    'BackgroundColor',CBG,'ForegroundColor',CFG,'FontSize',9,'Callback',@onChange);

uicontrol(fh,'Style','text','Units','normalized','Position',[.008 .528 .17 .028], ...
    'String','display scale','BackgroundColor',CBG,'ForegroundColor',CFG, ...
    'HorizontalAlignment','left','FontWeight','bold','FontSize',9);
S.hClimSld = uicontrol(fh,'Style','slider','Units','normalized', ...
    'Position',[.008 .50 .17 .026], 'Min',0.1,'Max',4,'Value',1, ...
    'Callback',@onChange);
S.hClimTxt = uicontrol(fh,'Style','text','Units','normalized', ...
    'Position',[.008 .472 .17 .026],'String','','BackgroundColor',CBG, ...
    'ForegroundColor',CFG,'HorizontalAlignment','left','FontName','Menlo','FontSize',8);

S.hPlay = uicontrol(fh,'Style','togglebutton','Units','normalized', ...
    'Position',[.008 .425 .052 .038], 'String','Play','FontSize',10, ...
    'Callback',@onPlay);
uicontrol(fh,'Style','text','Units','normalized','Position',[.066 .425 .04 .03], ...
    'String','fps','BackgroundColor',CBG,'ForegroundColor',CFG, ...
    'HorizontalAlignment','left','FontSize',9);
S.hFps = uicontrol(fh,'Style','edit','Units','normalized', ...
    'Position',[.094 .427 .04 .034],'String','15','FontSize',9, ...
    'BackgroundColor',CIN,'ForegroundColor',CFG);
uicontrol(fh,'Style','pushbutton','Units','normalized', ...
    'Position',[.14 .425 .038 .038],'String','MP4','FontSize',9, ...
    'TooltipString','write the current stack to an .mp4 next to the sweep file', ...
    'Callback',@onExport);

S.hInfo = uicontrol(fh,'Style','edit','Units','normalized', ...
    'Position',[.008 .01 .17 .40], 'Max',2,'Min',0,'Enable','inactive', ...
    'HorizontalAlignment','left','BackgroundColor',[.97 .97 .97], ...
    'ForegroundColor',CFG,'FontName','Menlo','FontSize',8,'String','');

% ---- centre: movie ----
S.axMov = axes(fh,'Units','normalized','Position',[.20 .13 .44 .82], ...
                 'Color',CBG,'XColor',CFG,'YColor',CFG);
axis(S.axMov,'off');
S.hSld = uicontrol(fh,'Style','slider','Units','normalized', ...
    'Position',[.20 .065 .38 .028],'Min',1,'Max',max(2,S.nFrame),'Value',1, ...
    'SliderStep',[1/max(1,S.nFrame-1) 10/max(1,S.nFrame-1)],'Callback',@onFrame);
S.hFrameTxt = uicontrol(fh,'Style','text','Units','normalized', ...
    'Position',[.585 .063 .055 .028],'String','','BackgroundColor',CBG, ...
    'ForegroundColor',CFG,'FontName','Menlo','FontSize',9);
if ~S.hasMovie
    set([S.hSld S.hPlay S.hShow S.hBlood S.hOrig S.hHoriz S.hClimSld],'Enable','off');
    text(S.axMov,0.5,0.5,{'no preview movies in this sweep file', ...
         'rerun section (5) with bdr.savePreview = true'}, ...
         'HorizontalAlignment','center','Color',[.5 .5 .5]);
end

% ---- right: plots ----
axArgs = {'Units','normalized','Color',CIN,'XColor',CFG,'YColor',CFG, ...
          'GridColor',[.75 .75 .75],'Box','on'};
S.axCorr = axes(fh, axArgs{:}, 'Position',[.675 .70 .135 .25]);
S.axDiff = axes(fh, axArgs{:}, 'Position',[.845 .70 .135 .25]);
S.axBAP  = axes(fh, axArgs{:}, 'Position',[.675 .40 .305 .22]);
S.axFast = axes(fh, axArgs{:}, 'Position',[.675 .09 .305 .22]);
S.axArgs = axArgs;

S.timer = [];
guidata(fh, S);
redraw(fh);

% ======================================================================
    function onChange(~,~)
        redraw(S.fh);
    end

    function onFrame(~,~)
        Sx = guidata(S.fh);
        showFrame(Sx, round(get(Sx.hSld,'Value')));
    end

    function onPlay(src_,~)
        Sx = guidata(S.fh);
        if get(src_,'Value') && Sx.hasMovie
            set(src_,'String','Stop');
            fps = max(1, min(60, str2double(get(Sx.hFps,'String'))));
            if isnan(fps), fps = 15; end
            stopTimer(Sx);
            Sx.timer = timer('ExecutionMode','fixedRate','Period',round(1/fps,3), ...
                             'TimerFcn',@(~,~) tick(), 'BusyMode','drop');
            guidata(S.fh, Sx);
            start(Sx.timer);
        else
            set(src_,'String','Play');
            stopTimer(guidata(S.fh));
        end
    end

    function tick()
        if ~isvalid(S.fh), return; end
        Sx = guidata(S.fh);
        fr = round(get(Sx.hSld,'Value')) + 1;
        if fr > Sx.nFrame, fr = 1; end
        set(Sx.hSld,'Value',fr);
        showFrame(Sx, fr);
    end

    function onExport(~,~)
        Sx = guidata(S.fh);
        if ~Sx.hasMovie, return; end
        [panels, labels] = buildPanels(Sx);                       %#ok<ASGLU>
        [fn, pth] = uiputfile('*.mp4','Write the current stack to', 'BVsweep_stack.mp4');
        if isequal(fn,0), return; end
        mont = assembleAll(Sx, panels);
        writeMov(fullfile(pth,fn), mont, [], ...
                 max(1,str2double(get(Sx.hFps,'String'))), 1, climNow(Sx));
        fprintf('[BVviewer] wrote %s\n', fullfile(pth,fn));
    end

    function onClose(~,~)
        stopTimer(guidata(S.fh));
        delete(S.fh);
    end
end

% ========================================================================
function redraw(fh)
S = guidata(fh);
sel = get(S.hList,'Value');
if isempty(sel), sel = 1; set(S.hList,'Value',1); end
S.sel = S.runIdx(sel);            % combination indices, in list order
guidata(fh, S);

if S.hasMovie
    fr = round(get(S.hSld,'Value'));
    fr = min(max(1,fr), S.nFrame);
    showFrame(S, fr);
end
drawPlots(S);
drawInfo(S);
end

% ------------------------------------------------------------------------
function [panels, labels] = buildPanels(S)
% Return the movies (cell) and labels for the current selection, in order.
panels = {};  labels = {};
mode   = get(S.hShow,'Value');            % 1 corrected, 2 blood, 3 diff, 4 original
if get(S.hOrig,'Value') || mode == 4
    panels{end+1} = S.prev_in;  labels{end+1} = 'ORIGINAL';
end
if mode == 4, return; end
for k = S.sel(:)'
    switch mode
        case 1, m = S.prev_corr{k};                  lb = 'corrected';
        case 2, m = S.prev_blood{k};                 lb = 'blood';
        case 3, m = S.prev_corr{k} - S.prev_in;      lb = 'corrected - original';
    end
    if isempty(m), continue; end
    panels{end+1} = m;                    labels{end+1} = sprintf('%s: %s', lb, S.SW(k).tag); %#ok<AGROW>
    if get(S.hBlood,'Value') && mode ~= 2 && ~isempty(S.prev_blood{k})
        panels{end+1} = S.prev_blood{k};  labels{end+1} = sprintf('blood: %s', S.SW(k).tag);  %#ok<AGROW>
    end
end
end

% ------------------------------------------------------------------------
function cl = climNow(S)
g  = get(S.hClimSld,'Value');
c0 = S.clim;
mid = mean(c0);  half = diff(c0)/2 / g;
cl = [mid-half mid+half];
if get(S.hShow,'Value') == 3       % difference: centre on zero
    cl = [-half half];
end
end

% ------------------------------------------------------------------------
function mont = assembleFrame(S, panels, fr)
sep = 2;
if get(S.hHoriz,'Value'), dcat = 2; else, dcat = 1; end
cl  = climNow(S);
n   = numel(panels);
mont = [];
for i = 1:n
    p = panels{i}(:,:,fr);
    if i > 1
        if dcat == 1, sepBlk = repmat(single(cl(2)), sep, size(p,2));
        else,         sepBlk = repmat(single(cl(2)), size(p,1), sep);
        end
        mont = cat(dcat, mont, sepBlk, p);
    else
        mont = p;
    end
end
end

% ------------------------------------------------------------------------
function mont = assembleAll(S, panels)
mont = [];
for fr = 1:S.nFrame
    mont = cat(3, mont, assembleFrame(S, panels, fr));
end
end

% ------------------------------------------------------------------------
function showFrame(S, fr)
[panels, labels] = buildPanels(S);
if isempty(panels), cla(S.axMov); axis(S.axMov,'off'); return; end
img = assembleFrame(S, panels, fr);
cl  = climNow(S);

ax = S.axMov;
h  = getappdata(ax,'imgHandle');
if isempty(h) || ~isgraphics(h) || ~isequal(size(get(h,'CData')), size(img))
    cla(ax);
    h = imagesc(ax, img, cl); axis(ax,'image'); axis(ax,'off');
    colormap(ax, gray(256));
    setappdata(ax,'imgHandle',h);
    % panel labels
    if get(S.hHoriz,'Value')
        w = size(panels{1},2);
        for i = 1:numel(labels)
            text(ax, (i-1)*(w+2)+6, 12, labels{i}, 'Color','y','FontSize',8, ...
                 'Interpreter','none','VerticalAlignment','top', ...
                 'BackgroundColor',[0 0 0 .45], 'Margin',1);
        end
    else
        hgt = size(panels{1},1);
        for i = 1:numel(labels)
            text(ax, 6, (i-1)*(hgt+2)+10, labels{i}, 'Color','y','FontSize',8, ...
                 'Interpreter','none','VerticalAlignment','top', ...
                 'BackgroundColor',[0 0 0 .45], 'Margin',1);
        end
    end
else
    set(h,'CData',img);
end
set(ax,'CLim',cl);
set(S.hFrameTxt,'String',sprintf('%d/%d', fr, S.nFrame));
set(S.hClimTxt,'String',sprintf('[%.1f %.1f]', cl(1), cl(2)));
drawnow limitrate                      % keeps playback smooth, drops frames if the
end                                    % renderer cannot keep up with the timer

% ------------------------------------------------------------------------
function drawPlots(S)
kp = S.sel(1);                               % primary condition

% --- unsmoothed correlation matrix of the primary condition -------------
cla(S.axCorr); set(S.axCorr, S.axArgs{3:end});
C = S.Cmat_raw_sw{kp};
if isempty(C), C = nan(S.nROI); end
imagesc(S.axCorr, C(S.ord,S.ord), [-0.25 1]); axis(S.axCorr,'image');
colormap(S.axCorr, S.cmapCorr); set(S.axCorr,'XTick',[],'YTick',[]);
title(S.axCorr, sprintf('%s\nraw corr,  r vs %dms = %.3f', S.SW(kp).tag, ...
      S.subwin, S.SW(kp).corrRaw_vs_in20), 'FontSize',8,'Interpreter','none', ...
      'Color', S.col.fg);

% --- difference from the original smoothed matrix -----------------------
cla(S.axDiff); set(S.axDiff, S.axArgs{3:end});
Dd = C - S.Cmat_in;
imagesc(S.axDiff, Dd(S.ord,S.ord), [-0.3 0.3]); axis(S.axDiff,'image');
colormap(S.axDiff, S.cmapCorr); set(S.axDiff,'XTick',[],'YTick',[]);
title(S.axDiff, sprintf('\\Delta vs ORIGINAL %d ms\nmean|\\Delta| = %.3f', ...
      S.subwin, mean(abs(Dd(S.triu)),'omitnan')), 'FontSize',8, 'Color', S.col.fg);

% --- bAP survival along the dendrite, all selected conditions -----------
cla(S.axBAP); set(S.axBAP, S.axArgs{3:end}); hold(S.axBAP,'on');
cm = lines(max(4, numel(S.sel)));
xax = S.dendaxis;  if isempty(xax), xax = (1:S.nROI)'; end
for i = 1:numel(S.sel)
    k = S.sel(i);
    if isempty(S.SW(k).bAP_ratio), continue; end
    plot(S.axBAP, xax, S.SW(k).bAP_ratio, '.-', 'Color',cm(i,:), 'MarkerSize',9, ...
         'LineWidth',1, 'DisplayName', S.SW(k).tag);
end
yline(S.axBAP,1,'k--','HandleVisibility','off');
yline(S.axBAP,0.9,'r:','HandleVisibility','off');
ylim(S.axBAP,[0 1.3]); grid(S.axBAP,'on');
xlabel(S.axBAP,'dendritic position (\mum)'); ylabel(S.axBAP,'bAP kept');
title(S.axBAP,'branch-local fast signal survival (want flat at 1)','FontSize',9,'Color',S.col.fg);
lg = legend(S.axBAP,'show','Location','southwest','Box','off','FontSize',7);
set(lg,'Interpreter','none','TextColor',S.col.fg,'Color',S.col.in);

% --- per-ROI trade-off: fast power kept vs bAP kept ---------------------
cla(S.axFast); set(S.axFast, S.axArgs{3:end}); hold(S.axFast,'on');
for i = 1:numel(S.sel)
    k = S.sel(i);
    if isempty(S.SW(k).fast_keep), continue; end
    scatter(S.axFast, S.SW(k).fast_keep, S.SW(k).bAP_ratio, 20, cm(i,:), 'filled', ...
            'MarkerFaceAlpha',0.6, 'DisplayName', S.SW(k).tag);
end
yline(S.axFast,1,'k--','HandleVisibility','off');
ylim(S.axFast,[0 1.3]); grid(S.axFast,'on');
xlabel(S.axFast,'fast power kept  (blood removal lowers this)');
ylabel(S.axFast,'bAP kept');
title(S.axFast,'per ROI: blood removed vs signal lost','FontSize',9,'Color',S.col.fg);
end

% ------------------------------------------------------------------------
function drawInfo(S)
% Kept narrow on purpose: the panel is ~34 characters wide, so rows are keyed
% by the condition number shown in the listbox rather than by the full tag.
L = {};
L{end+1} = sprintf('%3s %6s %6s %7s %6s', 'k','bAPmin','bAPmed','bloodRm','r|20ms');
L{end+1} = repmat('-',1,33);
for k = S.sel(:)'
    L{end+1} = sprintf('%3d %6.3f %6.3f %7.3f %6.3f', ...
        k, S.SW(k).bAP_min, S.SW(k).bAP_median, ...
        S.SW(k).blood_removed, S.SW(k).corrRaw_vs_in20);                   %#ok<AGROW>
end
L{end+1} = '';
kp = S.sel(1);
L{end+1} = sprintf('PRIMARY  k=%d', kp);
L{end+1} = sprintf('  %s', S.SW(kp).tag);
L{end+1} = sprintf('  bandpass     %g-%g Hz', S.SW(kp).bandpass(1), S.SW(kp).bandpass(2));
L{end+1} = sprintf('  preset       %d', S.SW(kp).preset);
L{end+1} = sprintf('  blood blocks %s', mat2str(S.SW(kp).blockSize*(1:S.SW(kp).preset+1)));
L{end+1} = sprintf('  blockSize    %d', S.SW(kp).blockSize);
L{end+1} = sprintf('  maxComp      %d', S.SW(kp).maxComp);
L{end+1} = sprintf('  rank1Method  %s', S.SW(kp).rank1);
L{end+1} = '';
L{end+1} = sprintf('  worst class  %s', S.SW(kp).bAP_worstclass);
L{end+1} = sprintf('   bAP kept    %.3f', S.SW(kp).bAP_worstclass_val);
L{end+1} = sprintf('  spike slope  %.3f', S.SW(kp).spike_slope);
L{end+1} = sprintf('  HF removed   %.3f', S.SW(kp).hf_removed);
L{end+1} = sprintf('  blood in/out %.2f', S.SW(kp).bloodratio);
L{end+1} = sprintf('  r(HF,blood)  %.2f', S.SW(kp).hf_vs_blood);
L{end+1} = '';
L{end+1} = '  unsmoothed corrected vs:';
L{end+1} = sprintf('   ORIGINAL 20ms  %.3f', S.SW(kp).corrRaw_vs_in20);
L{end+1} = sprintf('   ORIGINAL raw   %.3f', S.SW(kp).corrRaw_vs_inRaw);
L{end+1} = '';
L{end+1} = sprintf('  runtime      %.0f s%s', S.SW(kp).seconds, ...
                   repmat(' (cached)',1,S.SW(kp).fromCache));
set(S.hInfo,'String',L);
end

% ------------------------------------------------------------------------
function stopTimer(S)
if isfield(S,'timer') && ~isempty(S.timer) && isvalid(S.timer)
    stop(S.timer); delete(S.timer);
end
end

function v = getfielddef(D, name, dflt)
if isfield(D, name), v = D.(name); else, v = dflt; end
end

function s = trunc(s, n)
if numel(s) > n, s = s(1:n); end
end
