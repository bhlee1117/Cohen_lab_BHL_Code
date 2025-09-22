function ScrollXaxis(ax)
% ScrollXaxis Enable mouse wheel x-axis scrolling for a given axis
%
%   ScrollXaxis(gca) enables horizontal panning via mouse scroll.

    % Input validation
    if nargin < 1 || ~isgraphics(ax, 'axes')
        error('Input must be a valid axes handle.');
    end
    fig = ancestor(ax, 'figure');

    % Store initial scroll data
    winSize = diff(xlim(ax));
    setappdata(fig, 'ScrollData', struct('WindowSize', winSize, 'Axis', ax));

    % Set scroll callback
    fig.WindowScrollWheelFcn = @(src, event) scrollCallback(src, event);
end

function scrollCallback(fig, event)
    data = getappdata(fig, 'ScrollData');
    ax = data.Axis;

    xl = xlim(ax);
    step = data.WindowSize * 0.1;

    if event.VerticalScrollCount > 0
        new_xlim = xl + step;
    else
        new_xlim = xl - step;
    end

    % Optional: clamp to some bounds (e.g., based on data range) here

    % Apply new limits
    xlim(ax, new_xlim);

    % Update stored window size in case user zooms manually
    data.WindowSize = diff(new_xlim);
    setappdata(fig, 'ScrollData', data);
end
