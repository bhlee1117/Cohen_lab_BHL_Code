function interactive_glu_volt_viewer(t_volt, volt, t_glu, glu, glu_pts, ref_im, varargin)
% interactive_glu_volt_viewer  Interactive viewer: glutamate spots + voltage trace with frame slider.
%
%   viz_volt_glu(t_volt, volt, t_glu, glu, glu_pts, ref_im)
%   viz_volt_glu(t_volt, volt, t_glu, glu, glu_pts, ref_im, 'clim', [lo hi])
%
%   INPUTS
%     t_volt   : 1 x T1   voltage time axis (seconds)
%     volt     : N x T1   voltage traces (N ROIs)
%     t_glu    : 1 x T2   glutamate time axis (seconds)
%     glu      : M x T2   glutamate signal (M synapses)
%     glu_pts  : M x 2    glutamate spot positions [x, y] in ref_im pixels
%     ref_im   : H x W    reference image (grayscale)
%
%   OPTIONAL NAME-VALUE
%     'clim'   : [lo hi]  color limits for glutamate colormap (default: auto)
%
%   CONTROLS
%     Slider   : scrub through voltage frames
%     Click on voltage axes : jump to that time point
%
%   EXAMPLE
%     viz_volt_glu(t_vol, VoltTraces, t_glu, GluTraces, glu_xy, ref_img, ...
%                  'clim', [-0.05 0.2]);

% --- parse optional args ---
p = inputParser;
addParameter(p, 'clim', []);
parse(p, varargin{:});
clim_range = p.Results.clim;

% --- interpolate glu onto voltage time axis ---
glu_interp = interp1(t_glu, glu', t_volt, 'linear', 'extrap')';  % M x T1

% --- auto clim ---
if isempty(clim_range)
    clim_range = [prctile(glu_interp(:), 1), prctile(glu_interp(:), 99)];
end

T1 = numel(t_volt);
N  = size(volt, 1);
M  = size(glu_pts, 1);

% --- build figure ---
fig = figure('Name', 'Voltage + Glutamate Viewer', ...
             'NumberTitle', 'off', ...
             'Color', [0.12 0.12 0.12], ...
             'Units', 'normalized', ...
             'OuterPosition', [0.05 0.05 0.9 0.9]);

% Layout: image panel top 62%, volt panel middle 28%, slider bottom 5%
ax_img  = axes('Parent', fig, 'Position', [0.05 0.33 0.90 0.62]);
ax_volt = axes('Parent', fig, 'Position', [0.05 0.10 0.90 0.22]);

% --- draw reference image in grayscale ---
frame = 1;
axes(ax_img);
imshow(mat2gray(ref_im), 'Parent', ax_img);
hold(ax_img, 'on');
ax_img.YDir = 'reverse';

% --- draw glutamate spots using vec2cmap for colors ---
cvals  = glu_interp(:, frame);
colors = vec2cmap(cvals, 'hot', clim_range);   % M x 3 RGB
sc = scatter(ax_img, glu_pts(:,1), glu_pts(:,2), 60, colors, 'filled', ...
             'MarkerEdgeColor', 'none');

% --- colorbar: separate invisible axes so it doesn't affect the image ---
ax_cb = axes('Parent', fig, 'Position', ax_img.Position, 'Visible', 'off');
colormap(ax_cb, hot);
clim(ax_cb, clim_range);
cb = colorbar(ax_cb, 'Color', [0.85 0.85 0.85]);
cb.Label.String = 'dF/F';
cb.Label.Color  = [0.85 0.85 0.85];

ax_img.Visible = 'off';
title(ax_img, sprintf('Frame %d  |  t = %.3f s', frame, t_volt(frame)), ...
      'Color', [0.9 0.9 0.9], 'FontSize', 11);

% --- draw voltage traces ---
axes(ax_volt);
scale = max(prctile(volt, 99, 2)) * 3;
t_plot = t_volt;
plot(ax_volt, t_plot, volt' + (1:N)*scale, 'Color', [0.5 0.8 1.0], 'LineWidth', 0.5);
ax_volt.Color            = [0.12 0.12 0.12];
ax_volt.XColor           = [0.8 0.8 0.8];
ax_volt.YColor           = [0.8 0.8 0.8];
ax_volt.GridColor        = [0.3 0.3 0.3];
ax_volt.YTick            = (1:N) * scale;
ax_volt.YTickLabel       = arrayfun(@num2str, 1:N, 'UniformOutput', false);
ax_volt.FontSize         = 8;
xlabel(ax_volt, 'Time (s)', 'Color', [0.8 0.8 0.8]);
ylabel(ax_volt, 'ROI',      'Color', [0.8 0.8 0.8]);
hold(ax_volt, 'on');

% red frame indicator line
yl = ylim(ax_volt);
h_line = plot(ax_volt, [t_volt(frame) t_volt(frame)], yl, ...
              'r-', 'LineWidth', 1.5);

% --- slider ---
sl = uicontrol('Style', 'slider', ...
               'Min', 1, 'Max', T1, 'Value', frame, ...
               'SliderStep', [1/(T1-1), 50/(T1-1)], ...
               'Units', 'normalized', ...
               'Position', [0.05 0.04 0.90 0.03], ...
               'BackgroundColor', [0.25 0.25 0.25], ...
               'Callback', @onSlider);

% click on voltage axes to jump
ax_volt.ButtonDownFcn = @onAxClick;

% --- update function ---
    function update(fr)
        frame = max(1, min(T1, round(fr)));

        % update spot colors via vec2cmap
        cvals_new      = glu_interp(:, frame);
        sc.CData       = vec2cmap(cvals_new, 'hot', clim_range);

        % update red line
        h_line.XData = [t_volt(frame) t_volt(frame)];

        % update title
        title(ax_img, sprintf('Frame %d  |  t = %.3f s', frame, t_volt(frame)), ...
              'Color', [0.9 0.9 0.9], 'FontSize', 11);

        % sync slider
        sl.Value = frame;
        drawnow limitrate;
    end

    function onSlider(src, ~)
        update(round(src.Value));
    end

    function onAxClick(~, evt)
        t_clicked = evt.IntersectionPoint(1);
        [~, fr]   = min(abs(t_volt - t_clicked));
        update(fr);
    end

end