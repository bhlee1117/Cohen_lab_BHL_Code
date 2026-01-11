function set_font(Fontname)

axes = findall(gcf, 'Type', 'axes');
for ax = axes'
    ax.FontName = Fontname; % Tick labels
    ax.Title.FontName = Fontname; % Title
    ax.XLabel.FontName = Fontname; % X-axis label
    ax.YLabel.FontName = Fontname; % Y-axis label
    % Handle legend if present
    leg = findobj(ax, 'Type', 'legend');
    if ~isempty(leg)
        leg.FontName = Fontname;
    end
end

end