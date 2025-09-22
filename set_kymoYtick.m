function [yticks ylabels]=set_kymoYtick(dendaxis)
% set_kymoYtick Set Y-ticks and labels based on dendritic distances.
%
%   set_kymoYtick(dendaxis)
%   - dendaxis: vector of distances along the dendritic axis (1D array)

    idx = find(dendaxis == 0, 1);  % use first match only
    N = length(dendaxis);

    if isempty(idx) || idx == 1
        yticks = [1, N];
        yvals = [0, dendaxis(end)];
    else
        yticks = [1, idx, N];
        yvals = [dendaxis(1), 0, dendaxis(end)];
    end

    % Ensure uniqueness in tick positions/labels
    [yticks, ia] = unique(yticks, 'stable');
    ylabels = num2str(yvals(ia)', 3);

    set(gca, 'YTick', yticks, 'YTickLabel', ylabels);
end
