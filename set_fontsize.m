function set_fontsize(ftsz)

axes = findall(gcf, 'Type', 'axes');
for ax = axes'
    ax.FontSize = ftsz; % Tick labels
    ax.Title.FontSize = ftsz; % Title
    ax.XLabel.FontSize = ftsz; % X-axis label
    ax.YLabel.FontSize = ftsz; % Y-axis label
    % Handle legend if present
    leg = findobj(ax, 'Type', 'legend');
    if ~isempty(leg)
        leg.FontSize = ftsz;
    end
end

end