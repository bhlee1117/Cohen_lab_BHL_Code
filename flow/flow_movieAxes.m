function ax = flow_movieAxes(parent, pos, mov, opts)
% ax = flow_movieAxes(parent, pos, mov, opts)
%
% A scrollable movie panel: one axes plus a frame slider, embedded inside a GUI rather
% than opening its own window. The moviefixsc idea -- ONE intensity scale held fixed for
% every frame, so brightness changes you see are real and not the display rescaling
% itself -- but as a component.
%
% The movie is referenced, not copied: MATLAB's copy-on-write means storing a 4 GB S in
% the axes UserData costs nothing as long as nobody writes to it, and nothing here does.
%
% INPUT
%   parent  uipanel or figure to draw into
%   pos     [x y w h] normalized within parent
%   mov     [nrow x ncol x nT]
%   opts    optional struct:
%     .clim    [lo hi]      default: 1st-99th percentile of a frame sample, held fixed
%     .cmap    colormap     default gray
%     .title   char         drawn above, with the frame number appended
%     .frame   initial frame (default 1)
%     .step    slider arrow step in frames (default 1)
%     .overlay function handle @(ax, k) called after each frame is drawn, so the caller
%              can paint per-frame annotation on top. Objects it creates must be tagged
%              'flowOverlay'; everything with that tag is deleted before the next call,
%              which is what keeps the overlay from accumulating as you scroll.
%
% OUTPUT  ax  the image axes (ax.UserData holds the movie and the handles)
%
% 2026.08.28, Byung Hun Lee

if nargin < 4, opts = struct; end
d = struct('clim',[], 'cmap',gray(256), 'title','', 'frame',1, 'step',1, 'overlay',[]);
f = fieldnames(d);
for k = 1:numel(f)
    if ~isfield(opts,f{k}) || isempty(opts.(f{k})), opts.(f{k}) = d.(f{k}); end
end
nT = size(mov,3);

if isempty(opts.clim)
    % sample rather than scan every frame: the scale only has to be representative, and
    % percentiles over the whole of a big movie are slow enough to be noticeable
    smp = mov(:,:,unique(round(linspace(1, nT, min(nT, 40)))));
    opts.clim = double(prctile(single(smp(:)), [1 99]));
    if ~(opts.clim(2) > opts.clim(1)), opts.clim = double([min(smp(:)) max(smp(:))]); end
    if ~(opts.clim(2) > opts.clim(1)), opts.clim = opts.clim(1) + [0 1]; end
end

% axes on top, slider strip underneath
hSl = 0.10;
ax  = axes(parent, 'Units','normalized', ...
    'Position', [pos(1) pos(2)+hSl*pos(4) pos(3) pos(4)*(1-hSl)]);
im  = imagesc(ax, mov(:,:,opts.frame), opts.clim);
axis(ax, 'image', 'off');  colormap(ax, opts.cmap);

sl = uicontrol(parent, 'Style','slider', 'Units','normalized', ...
    'Position', [pos(1) pos(2) pos(3)*0.78 hSl*pos(4)*0.85], ...
    'Min', 1, 'Max', max(nT, 1.0001), 'Value', opts.frame, ...
    'SliderStep', [opts.step/max(nT-1,1), max(10*opts.step,1)/max(nT-1,1)]);
ed = uicontrol(parent, 'Style','edit', 'Units','normalized', ...
    'Position', [pos(1)+pos(3)*0.80 pos(2) pos(3)*0.20 hSl*pos(4)*0.85], ...
    'String', num2str(opts.frame), 'HorizontalAlignment','center');

hold(ax, 'on');
ax.UserData = struct('mov',mov, 'im',im, 'sl',sl, 'ed',ed, 'nT',nT, ...
    'title',opts.title, 'overlay',opts.overlay);
set(sl, 'Callback', @(s,~) go(ax, round(get(s,'Value'))));
set(ed, 'Callback', @(s,~) go(ax, str2double(get(s,'String'))));
if nT < 2, set(sl, 'Enable','off'); end
go(ax, opts.frame);
end

function go(ax, k)
u = ax.UserData;
if ~isfinite(k), k = 1; end
k = max(1, min(u.nT, round(k)));
set(u.im, 'CData', u.mov(:,:,k));
delete(findobj(ax, 'Tag', 'flowOverlay'));      % previous frame's annotation
if ~isempty(u.overlay)
    try
        u.overlay(ax, k);
    catch ME
        warning('flow_movieAxes:overlay', 'overlay failed at frame %d: %s', k, ME.message);
    end
end
set(u.sl, 'Value', k);
set(u.ed, 'String', num2str(k));
if isempty(u.title)
    title(ax, sprintf('frame %d / %d', k, u.nT), 'FontSize', 9);
else
    title(ax, sprintf('%s   (frame %d / %d)', u.title, k, u.nT), 'FontSize', 9);
end
end
