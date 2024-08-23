function errorbarxy(x, y, xerr, yerr, cmap,varargin)
    % Plot the data points
    h = plot(x, y, varargin{:}); hold on;

    % Add horizontal error bars
    for k = 1:length(x)
        plot([x(k) - xerr(k), x(k) + xerr(k)], [y(k), y(k)], 'Color',cmap);
    end

    % Add vertical error bars
    for k = 1:length(y)
        plot([x(k), x(k)], [y(k) - yerr(k), y(k) + yerr(k)], 'Color',cmap);
    end

    hold off;
end