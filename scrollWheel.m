function scrollWheel(src, event)
    ax = src.CurrentAxes;
    xlims = ax.XLim;
    range = diff(xlims);
    step = range * 0.1; % Adjust the step size as needed

    if event.VerticalScrollCount > 0
        new_xlims = xlims + step;
    else
        new_xlims = xlims - step;
    end

    ax.XLim = new_xlims;
end